"""
09_neighbor_purity_outlier_audit.py
For each subject, compute:
  - Distance-based outlier score (z-score of mean distance to all others)
  - Nearest-neighbor label purity at k = 3, 5, 10
    (proportion of k-NN with the same case/control label)

A value near 0.50 at all k indicates the evidence space is intermixed —
no natural case/control clustering at that scale.
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

os.environ["PANDAS_USE_NUMEXPR"] = "0"

ZSCORE_FILE = "results/evidence_modules_v3.no_lowcell_matrix.zscore.tsv"
META_COLS   = ["subject_id", "label", "sex", "age_bin", "has_low_cell_count_member"]
K_VALUES    = [3, 5, 10]


def main():
    df       = pd.read_csv(ZSCORE_FILE, sep="\t")
    features = [c for c in df.columns if c not in META_COLS]
    X        = df[features].values.astype(float)
    D        = pairwise_distances(X, metric="euclidean")
    n        = len(df)

    # Outlier audit
    out_rows = []
    for i in range(n):
        d = np.array([D[i, j] for j in range(n) if j != i])
        out_rows.append({
            "subject_id": df.loc[i, "subject_id"],
            "label": df.loc[i, "label"],
            "mean_distance_to_all":   float(np.mean(d)),
            "median_distance_to_all": float(np.median(d)),
            "min_distance_to_any":    float(np.min(d)),
            "max_distance_to_any":    float(np.max(d)),
        })
    out_df = pd.DataFrame(out_rows)
    for c in ["mean_distance_to_all", "median_distance_to_all",
              "min_distance_to_any",  "max_distance_to_any"]:
        mu, sd = out_df[c].mean(), out_df[c].std(ddof=0)
        out_df[f"{c}_z"] = 0.0 if (sd == 0 or pd.isna(sd)) else (out_df[c] - mu) / sd
    out_df.sort_values(["mean_distance_to_all_z", "median_distance_to_all_z"],
                       ascending=False).to_csv(
        "results/evidence_modules_v3_outlier_audit.tsv", sep="\t", index=False
    )

    # Neighbor purity
    purity_rows = []
    for i in range(n):
        order = [j for j in np.argsort(D[i]) if j != i]
        for k in K_VALUES:
            take = order[:k]
            same = sum(df.loc[j, "label"] == df.loc[i, "label"] for j in take)
            purity_rows.append({
                "subject_id":            df.loc[i, "subject_id"],
                "label":                 df.loc[i, "label"],
                "k":                     k,
                "same_label_neighbors":  same,
                "diff_label_neighbors":  k - same,
                "same_label_fraction":   same / k,
            })
    purity_df = pd.DataFrame(purity_rows)
    purity_df.to_csv("results/evidence_modules_v3_neighbor_purity.tsv", sep="\t", index=False)

    # Summary
    sum_rows = []
    for (lab, k), sub in purity_df.groupby(["label", "k"]):
        vals = sub["same_label_fraction"].values
        sum_rows.append({
            "label": lab, "k": int(k),
            "n_subjects": int(len(sub)),
            "mean_same_label_fraction":   float(np.mean(vals)),
            "median_same_label_fraction": float(np.median(vals)),
            "n_subjects_gte_0.8": int(np.sum(vals >= 0.8)),
            "n_subjects_lte_0.2": int(np.sum(vals <= 0.2)),
        })
    summary = pd.DataFrame(sum_rows)
    summary.to_csv("results/evidence_modules_v3_neighbor_labelmix_summary.tsv",
                   sep="\t", index=False)

    print("Neighbor purity summary (k=5):")
    print(summary[summary.k == 5][
        ["label", "n_subjects", "mean_same_label_fraction"]
    ].to_string(index=False))
    print("\n(~0.50 = random mixing; evidence space is intermixed)")


if __name__ == "__main__":
    main()
