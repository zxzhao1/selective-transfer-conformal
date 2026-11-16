# Generated ADME data

Raw and processed ADME files are not stored in the repository. Run the
following commands from `real_data/` to download the public TDC datasets and
reconstruct the analysis inputs:

```bash
python3 22_prepare_tdc_adme.py
Rscript 23_finalize_tdc_adme.R
```

The generated `candidate/` and `processed/` directories are ignored by Git.

