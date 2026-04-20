"""
03_build_pseudobulk.py
Aggregate all cells within each sample into a pseudobulk expression vector
(gene-level sum across all cells), producing a genes × samples count matrix.
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

os.environ["PANDAS_USE_NUMEXPR"] = "0"

MANIFEST = "meta/sample_manifest_paired_with_files.fixed.tsv"
QC_FILE  = "results/sample_level_qc.tsv"
OUT_COUNTS = "results/pseudobulk_counts.tsv"
OUT_META   = "results/pseudobulk_sample_metadata.tsv"


def main():
    manifest = pd.read_csv(MANIFEST, sep="\t")
    qc = pd.read_csv(QC_FILE, sep="\t").set_index("gsm")

    counts_dict = {}
    meta_rows = []
    gene_index = None

    for _, r in manifest.iterrows():
        print(f"  Pseudobulk: {r['sample_title']}")
        adata = sc.read_10x_h5(r["raw_file"])
        adata.var_names_make_unique()
        X = adata.X

        gene_sum = (
            np.asarray(X.sum(axis=0)).ravel()
            if sparse.issparse(X)
            else X.sum(axis=0)
        )
        genes = adata.var_names.tolist()

        if gene_index is None:
            gene_index = genes
        elif genes != gene_index:
            raise ValueError(f"Gene order mismatch for sample {r['sample_title']}")

        counts_dict[r["sample_title"]] = gene_sum

        q = qc.loc[r["gsm"]]
        meta_rows.append({
            "sample_name": r["sample_title"],
            "gsm": r["gsm"],
            "subject_id": r.get("subject_id", ""),
            "label": r.get("label", ""),
            "day": r.get("timepoint_from_title", r.get("day", "")),
            "sex": r.get("sex", ""),
            "age_bin": r.get("age_bin", ""),
            "chromium_batch": r.get("chromium_batch", ""),
            "n_cells": int(q["n_cells"]),
            "mean_counts_per_cell": float(q["mean_counts_per_cell"]),
            "mean_detected_genes_per_cell": float(q["mean_detected_genes_per_cell"]),
            "low_cell_count": bool(q["low_cell_count"]),
        })

    counts_df = pd.DataFrame(counts_dict, index=gene_index)
    counts_df.index.name = "gene"
    meta_df = pd.DataFrame(meta_rows)

    counts_df.to_csv(OUT_COUNTS, sep="\t")
    meta_df.to_csv(OUT_META, sep="\t", index=False)

    print(f"\nPseudobulk matrix: {counts_df.shape} (genes × samples)")
    print(f"Wrote {OUT_COUNTS}")
    print(f"Wrote {OUT_META}")


if __name__ == "__main__":
    main()
