#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

root <- normalizePath(getwd(), mustWork = TRUE)
processed <- file.path(root, "data", "processed")
x <- as.matrix(readMM(file.path(processed, "adme_fingerprints.mtx")))
ids <- fread(file.path(processed, "adme_expression_ids.csv"))$expression_id
stopifnot(nrow(x) == length(ids), ncol(x) == 1024L, !anyNA(x))
storage.mode(x) <- "double"
rownames(x) <- ids
colnames(x) <- sprintf("morgan_%04d", seq_len(ncol(x)))
saveRDS(x, file.path(processed, "adme_expression.rds"), compress = FALSE)
cat(sprintf("Saved %d x %d ADME expression matrix.\n", nrow(x), ncol(x)))
