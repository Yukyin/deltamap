"""
04_build_delta.py
Compute within-subject temporal delta: D2 − D1 pseudobulk expression.
Only subjects with both D1 and D2 samples are retained.
Propagates QC flags (has_low_cell_count_member) for downstream filtering.
"""
import os
import pandas as pd

os.environ["PANDAS_USE_NUMEXPR"] = "0"

COUNTS_FILE = "results/pseudobulk_counts.tsv"
META_FILE   = "results/pseudobulk_sample_metadata.tsv"
OUT_DELTA   = "results/pseudobulk_delta.tsv"
OUT_META    = "results/pseudobulk_delta_metadata.tsv"


def main():
    counts = pd.read_csv(COUNTS_FILE, sep="\t").set_index("gene")
    meta   = pd.read_csv(META_FILE, sep="\t")

    delta_cols = {}
    meta_rows  = []
    skipped    = []

    for sid in sorted(meta["subject_id"].unique()):
        sub  = meta[meta["subject_id"] == sid].copy()
        days = set(sub["day"].tolist())

        if not {"D1", "D2"}.issubset(days):
            skipped.append(sid)
            continue

        d1 = sub[sub["day"] == "D1"].iloc[0]
        d2 = sub[sub["day"] == "D2"].iloc[0]
        s1, s2 = d1["sample_name"], d2["sample_name"]

        if s1 not in counts.columns or s2 not in counts.columns:
            skipped.append(sid)
            continue

        name = f"{sid}_DELTA"
        delta_cols[name] = counts[s2] - counts[s1]

        meta_rows.append({
            "delta_name": name,
            "subject_id": sid,
            "label": d1["label"],
            "sex": d1["sex"],
            "age_bin": d1["age_bin"],
            "D1_sample": s1,
            "D2_sample": s2,
            "D1_n_cells": d1["n_cells"],
            "D2_n_cells": d2["n_cells"],
            "D1_low_cell_count": d1["low_cell_count"],
            "D2_low_cell_count": d2["low_cell_count"],
            "has_low_cell_count_member": bool(
                d1["low_cell_count"] or d2["low_cell_count"]
            ),
        })

    delta_df = pd.DataFrame(delta_cols)
    delta_df.index.name = "gene"
    meta_df = pd.DataFrame(meta_rows)

    delta_df.to_csv(OUT_DELTA, sep="\t")
    meta_df.to_csv(OUT_META, sep="\t", index=False)

    print(f"Delta matrix: {delta_df.shape}  (genes × subjects)")
    print(f"Paired subjects: {len(meta_df)}")
    if skipped:
        print(f"Skipped (missing D1 or D2): {skipped}")
    print(f"Subjects with low-cell-count member: "
          f"{meta_df['has_low_cell_count_member'].sum()}")
    print(f"Wrote {OUT_DELTA}")
    print(f"Wrote {OUT_META}")


if __name__ == "__main__":
    main()
