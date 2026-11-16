#!/usr/bin/env python3

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy import sparse
from scipy.io import mmwrite
from tdc.single_pred import ADME


ROOT = Path.cwd()
CANDIDATE = ROOT / "data" / "candidate" / "tdc_adme"
PROCESSED = ROOT / "data" / "processed"
CANDIDATE.mkdir(parents=True, exist_ok=True)
PROCESSED.mkdir(parents=True, exist_ok=True)

SPECS = [
    ("caco2", "Caco2_Wang", "Half_Life_Obach"),
    ("vdss", "VDss_Lombardo", "Caco2_Wang"),
    ("ppbr", "PPBR_AZ", "Solubility_AqSolDB"),
    ("half_life", "Half_Life_Obach", "PPBR_AZ"),
    ("clearance_hepatocyte", "Clearance_Hepatocyte_AZ", "VDss_Lombardo"),
    ("clearance_microsome", "Clearance_Microsome_AZ", "Lipophilicity_AstraZeneca"),
    ("solubility", "Solubility_AqSolDB", "Clearance_Hepatocyte_AZ"),
    ("lipophilicity", "Lipophilicity_AstraZeneca", "Clearance_Microsome_AZ"),
    ("hydration", "HydrationFreeEnergy_FreeSolv", "Half_Life_Obach"),
]
MAX_SOURCE_N = 800
FP_SIZE = 1024
fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FP_SIZE)


def load_dataset(name: str) -> pd.DataFrame:
    frame = ADME(name=name, path=str(CANDIDATE)).get_data(format="df")
    frame = frame[["Drug_ID", "Drug", "Y"]].dropna().copy()
    frame["Y"] = pd.to_numeric(frame["Y"], errors="coerce")
    frame = frame[np.isfinite(frame["Y"])].copy()
    frame["canonical_smiles"] = frame["Drug"].map(canonical_smiles)
    frame = frame.dropna(subset=["canonical_smiles"]).reset_index(drop=True)
    frame = (
        frame.groupby("canonical_smiles", as_index=False, sort=True)
        .agg(
            Drug_ID=("Drug_ID", "first"),
            Y=("Y", "mean"),
            replicate_count=("Y", "size"),
        )
    )
    frame["Drug"] = frame["canonical_smiles"]
    frame["compound_id"] = frame["canonical_smiles"].map(compound_identifier)
    if frame["compound_id"].duplicated().any():
        raise RuntimeError(f"Compound identifier collision in {name}")
    return frame


def canonical_smiles(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else Chem.MolToSmiles(mol, canonical=True)


def compound_identifier(smiles: str) -> str:
    return hashlib.sha256(smiles.encode("utf-8")).hexdigest()[:20]


def compound_bucket(smiles: str) -> int:
    return int(hashlib.md5(smiles.encode("utf-8")).hexdigest()[:8], 16) % 5


def sample_source(frame: pd.DataFrame, seed: int) -> pd.DataFrame:
    if len(frame) <= MAX_SOURCE_N:
        return frame.copy()
    return frame.sample(MAX_SOURCE_N, random_state=seed).reset_index(drop=True)


all_rows = []
responses = []
tasks = []
diagnostics = []

for family_index, (family, target_name, heterogeneous_name) in enumerate(SPECS):
    target_data = load_dataset(target_name)
    target_data["bucket"] = target_data["canonical_smiles"].map(compound_bucket)
    target = target_data[target_data["bucket"] == 0].copy()
    related = target_data[target_data["bucket"] != 0].copy()
    heterogeneous = load_dataset(heterogeneous_name)
    heterogeneous_before_overlap = len(heterogeneous)
    heterogeneous = heterogeneous[
        ~heterogeneous["compound_id"].isin(target["compound_id"])
    ].reset_index(drop=True)
    heterogeneous_overlap_removed = (
        heterogeneous_before_overlap - len(heterogeneous)
    )
    related = sample_source(related, 20261114 + family_index * 10 + 1)
    heterogeneous = sample_source(
        heterogeneous, 20261114 + family_index * 10 + 2
    )
    if len(target) < 80 or len(related) < 200 or len(heterogeneous) < 200:
        raise RuntimeError(
            f"Insufficient samples for {family}: "
            f"{len(target)}, {len(related)}, {len(heterogeneous)}"
        )

    domains = [
        (f"{family}_target", target, target_name, 0),
        (f"{family}_related", related, target_name, 0),
        (
            f"{family}_heterogeneous",
            heterogeneous,
            heterogeneous_name,
            heterogeneous_overlap_removed,
        ),
    ]
    for domain, frame, endpoint, overlap_removed in domains:
        for _, row in frame.reset_index(drop=True).iterrows():
            expression_id = f"{domain}_{row['compound_id']}"
            all_rows.append((expression_id, row["canonical_smiles"]))
            responses.append({
                "drug": family,
                "domain": domain,
                "expression_id": expression_id,
                "compound_id": row["compound_id"],
                "response_value": float(row["Y"]),
                "response_measure": endpoint,
                "replicate_count": int(row["replicate_count"]),
            })
        diagnostics.append({
            "family": family,
            "domain": domain,
            "endpoint": endpoint,
            "n": len(frame),
            "raw_rows": int(frame["replicate_count"].sum()),
            "duplicates_aggregated": int(
                frame["replicate_count"].sum() - len(frame)
            ),
            "target_overlap_removed": overlap_removed,
            "response_sd": float(frame["Y"].astype(float).std()),
        })

    rule = (
        "responses aggregated by canonical compound; target canonical-compound hash "
        "bucket 0/5; related source is remaining compounds from the same assay; "
        "heterogeneous source is the locked different assay after removing target "
        "compounds; sources capped at 800 without method-performance filtering"
    )
    tasks.extend([
        {
            "task": f"adme_{family}_related",
            "drug": family,
            "target": f"{family}_target",
            "sources": f"{family}_related",
            "source_regime": "same assay, disjoint compounds",
            "selection_rule": rule,
        },
        {
            "task": f"adme_{family}_heterogeneous",
            "drug": family,
            "target": f"{family}_target",
            "sources": f"{family}_heterogeneous",
            "source_regime": "different ADME assay",
            "selection_rule": rule,
        },
    ])

response_frame = pd.DataFrame(responses)
if response_frame.duplicated(["domain", "compound_id"]).any():
    raise RuntimeError("Duplicate compound remains within an ADME domain")
if response_frame["expression_id"].duplicated().any():
    raise RuntimeError("Duplicate ADME expression identifier")
for task in tasks:
    target_compounds = set(response_frame.loc[
        response_frame["domain"] == task["target"], "compound_id"
    ])
    source_compounds = set(response_frame.loc[
        response_frame["domain"] == task["sources"], "compound_id"
    ])
    if target_compounds & source_compounds:
        raise RuntimeError(f"Target/source compound overlap in {task['task']}")

fingerprints = np.zeros((len(all_rows), FP_SIZE), dtype=np.uint8)
valid = np.ones(len(all_rows), dtype=bool)
for i, (_, smiles) in enumerate(all_rows):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        valid[i] = False
        continue
    fingerprints[i, :] = fp_generator.GetFingerprintAsNumPy(mol)
if not valid.all():
    raise RuntimeError(f"Fingerprint generation failed for {(~valid).sum()} rows")

mmwrite(PROCESSED / "adme_fingerprints.mtx", sparse.csr_matrix(fingerprints))
pd.DataFrame({"expression_id": [row[0] for row in all_rows]}).to_csv(
    PROCESSED / "adme_expression_ids.csv", index=False
)
response_frame.to_csv(PROCESSED / "adme_response.csv", index=False)
pd.DataFrame(tasks).to_csv(PROCESSED / "adme_tasks.csv", index=False)
pd.DataFrame(diagnostics).to_csv(
    PROCESSED / "adme_task_diagnostics.csv", index=False
)
print(
    f"Prepared {len(all_rows)} assay rows, {FP_SIZE} fingerprint features, "
    f"and {len(tasks)} locked tasks."
)
