"""
02_sample_level_qc.py
For each sample, load the 10x .h5 file and compute basic QC metrics:
  - n_cells
  - mean_counts_per_cell
  - mean_detected_genes_per_cell
Flags samples with fewer than 1000 cells as low_cell_count outliers.
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

os.environ["PANDAS_USE_NUMEXPR"] = "0"

MANIFEST = "meta/sample_manifest_paired_with_files.fixed.tsv"
OUT_QC = "results/sample_level_qc.tsv"
LOW_CELL_THRESHOLD = 1000


def compute_sample_qc(h5_path: str) -> dict:
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    X = adata.X
    if sparse.issparse(X):
        cell_counts = np.asarray(X.sum(axis=1)).ravel()
        detected = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        cell_counts = X.sum(axis=1)
        detected = (X > 0).sum(axis=1)
    return {
        "n_cells": int(adata.n_obs),
        "mean_counts_per_cell": float(np.mean(cell_counts)),
        "mean_detected_genes_per_cell": float(np.mean(detected)),
    }


def main():
    manifest = pd.read_csv(MANIFEST, sep="\t")
    rows = []
    for _, r in manifest.iterrows():
        print(f"  QC: {r['sample_title']} ({r['raw_file']})")
        qc = compute_sample_qc(r["raw_file"])
        rows.append({
            "gsm": r["gsm"],
            "sample_title": r["sample_title"],
            "subject_id": r.get("subject_id", ""),
            "label": r.get("label", ""),
            "day": r.get("timepoint_from_title", r.get("day", "")),
            "sex": r.get("sex", ""),
            "age_bin": r.get("age_bin", ""),
            "chromium_batch": r.get("chromium_batch", ""),
            **qc,
        })

    qc_df = pd.DataFrame(rows)
    qc_df["low_cell_count"] = qc_df["n_cells"] < LOW_CELL_THRESHOLD
    qc_df.to_csv(OUT_QC, sep="\t", index=False)

    flagged = qc_df[qc_df["low_cell_count"]]
    flagged.to_csv("results/sample_level_qc_flagged.tsv", sep="\t", index=False)

    print(f"\nWrote {len(qc_df)} rows to {OUT_QC}")
    print(f"Flagged {len(flagged)} low-cell-count samples (<{LOW_CELL_THRESHOLD} cells)")


if __name__ == "__main__":
    main()
