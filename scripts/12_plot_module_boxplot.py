"""
12_plot_module_boxplot.py
Side-by-side boxplots of module scores for case vs. control subjects.
One panel per evidence module.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.environ["PANDAS_USE_NUMEXPR"] = "0"

MATRIX_FILE = "results/evidence_modules_v3.no_lowcell_matrix.tsv"
META_COLS   = ["subject_id", "label", "sex", "age_bin", "has_low_cell_count_member"]


def main():
    df       = pd.read_csv(MATRIX_FILE, sep="\t")
    features = [c for c in df.columns if c not in META_COLS]

    # Save summary stats
    rows = []
    for feat in features:
        for lab, sub in df.groupby("label"):
            vals = sub[feat].dropna().values
            rows.append({
                "module": feat, "label": lab,
                "n": len(vals),
                "mean":   float(np.mean(vals)),
                "median": float(np.median(vals)),
                "std":    float(np.std(vals, ddof=1)) if len(vals) > 1 else float("nan"),
                "q25":    float(np.percentile(vals, 25)),
                "q75":    float(np.percentile(vals, 75)),
                "min":    float(np.min(vals)),
                "max":    float(np.max(vals)),
            })
    pd.DataFrame(rows).to_csv(
        "results/evidence_modules_v3_group_box_summary.tsv", sep="\t", index=False
    )

    # Plot
    fig, axes = plt.subplots(1, len(features), figsize=(4 * len(features), 5), squeeze=False)
    axes = axes[0]

    for i, feat in enumerate(features):
        ax   = axes[i]
        ctrl = df.loc[df.label == "control", feat].dropna().values
        case = df.loc[df.label == "case",    feat].dropna().values
        ax.boxplot([ctrl, case], tick_labels=["control", "case"])
        ax.set_title(feat.replace("_v3", "").replace("_", "\n"), fontsize=11)
        ax.set_ylabel("module score (D2 − D1)")
        ax.axhline(0, linestyle="--", linewidth=1, color="gray")
        ax.tick_params(axis="x", labelsize=10)

    plt.suptitle("Evidence module scores: case vs. control", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig("results/evidence_modules_v3_case_control_boxplot.png",
                dpi=200, bbox_inches="tight")
    plt.close()
    print("Wrote results/evidence_modules_v3_case_control_boxplot.png")


if __name__ == "__main__":
    main()
