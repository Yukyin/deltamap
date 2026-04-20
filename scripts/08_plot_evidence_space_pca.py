"""
08_plot_evidence_space_pca.py
PCA of the 4D evidence module score matrix (z-scored).
Saves subject coordinates, pairwise distances, nearest-neighbor table,
and a case/control scatter plot.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances

os.environ["PANDAS_USE_NUMEXPR"] = "0"

MATRIX_FILE = "results/evidence_modules_v3.no_lowcell_matrix.tsv"
META_COLS   = ["subject_id", "label", "sex", "age_bin", "has_low_cell_count_member"]

COLORS = {"control": "#1f77b4", "case": "#ff7f0e"}
CENTROID_COLORS = {"control": "#2ca02c", "case": "#d62728"}


def main():
    df       = pd.read_csv(MATRIX_FILE, sep="\t")
    features = [c for c in df.columns if c not in META_COLS]

    # Z-score
    z = df.copy()
    for c in features:
        mu, sd = z[c].mean(), z[c].std(ddof=0)
        z[c]   = 0.0 if (sd == 0 or pd.isna(sd)) else (z[c] - mu) / sd
    z.to_csv("results/evidence_modules_v3.no_lowcell_matrix.zscore.tsv", sep="\t", index=False)

    # PCA
    X   = z[features].values.astype(float)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)

    pd.DataFrame({
        "module": features,
        "PC1_loading": pca.components_[0],
        "PC2_loading": pca.components_[1],
    }).to_csv("results/evidence_modules_v3_pca_loadings.tsv", sep="\t", index=False)

    pca_df = z[META_COLS].copy()
    pca_df["PC1"] = pcs[:, 0]
    pca_df["PC2"] = pcs[:, 1]
    pca_df["distance_from_origin"] = (pca_df["PC1"] ** 2 + pca_df["PC2"] ** 2) ** 0.5
    pca_df.to_csv("results/evidence_modules_v3_subject_pca.tsv", sep="\t", index=False)

    # Pairwise distances
    D        = pairwise_distances(pcs, metric="euclidean")
    dist_rows, nn_rows = [], []
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            dist_rows.append({
                "subject1": df.loc[i, "subject_id"], "label1": df.loc[i, "label"],
                "subject2": df.loc[j, "subject_id"], "label2": df.loc[j, "label"],
                "euclidean_distance": float(D[i, j]),
            })
        neighbors = [j for j in np.argsort(D[i]) if j != i][:3]
        for rank, j in enumerate(neighbors, start=1):
            nn_rows.append({
                "subject_id": df.loc[i, "subject_id"],
                "label": df.loc[i, "label"],
                "neighbor_rank": rank,
                "neighbor_subject_id": df.loc[j, "subject_id"],
                "neighbor_label": df.loc[j, "label"],
                "distance": float(D[i, j]),
            })

    pd.DataFrame(dist_rows).sort_values("euclidean_distance", ascending=False).to_csv(
        "results/evidence_modules_v3_subject_pairwise_distance.tsv", sep="\t", index=False
    )
    pd.DataFrame(nn_rows).to_csv(
        "results/evidence_modules_v3_subject_nearest_neighbors.tsv", sep="\t", index=False
    )

    # Summary statistics
    summary = pd.DataFrame([
        {"item": "explained_variance_ratio_PC1", "value": float(pca.explained_variance_ratio_[0])},
        {"item": "explained_variance_ratio_PC2", "value": float(pca.explained_variance_ratio_[1])},
        {"item": "case_PC1_mean",   "value": float(pca_df.loc[pca_df.label == "case",    "PC1"].mean())},
        {"item": "case_PC2_mean",   "value": float(pca_df.loc[pca_df.label == "case",    "PC2"].mean())},
        {"item": "control_PC1_mean","value": float(pca_df.loc[pca_df.label == "control", "PC1"].mean())},
        {"item": "control_PC2_mean","value": float(pca_df.loc[pca_df.label == "control", "PC2"].mean())},
    ])
    summary.to_csv("results/evidence_modules_v3_subject_pca_summary.tsv", sep="\t", index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    for lab in ["control", "case"]:
        sub = pca_df[pca_df.label == lab]
        ax.scatter(sub["PC1"], sub["PC2"], s=70, alpha=0.9, label=lab, color=COLORS[lab])
    for lab, color in CENTROID_COLORS.items():
        sub = pca_df[pca_df.label == lab]
        ax.scatter(
            [sub["PC1"].mean()], [sub["PC2"].mean()],
            s=180, marker="X", color=color, label=f"{lab} centroid",
        )
    ax.axhline(0, linestyle="--", linewidth=1)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_title("Evidence space — module score PCA")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig("results/evidence_modules_v3_subject_pca.png", dpi=250, bbox_inches="tight")
    plt.close()

    print("PCA explained variance:", pca.explained_variance_ratio_)
    print("Wrote results/evidence_modules_v3_subject_pca.png")


if __name__ == "__main__":
    main()
