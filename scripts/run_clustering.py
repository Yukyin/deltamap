"""
run_clustering.py
KMeans and hierarchical clustering on the z-scored evidence module matrix.
Searches k=2..6, selects best k by silhouette score.
Outputs cluster assignments, centroids, and Fisher enrichment for case/control.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

os.environ["PANDAS_USE_NUMEXPR"] = "0"

MATRIX_FILE = "results/evidence_modules_v3.no_lowcell_matrix.tsv"
META_COLS   = ["subject_id", "label", "sex", "age_bin", "has_low_cell_count_member"]


def summarize_clusters(df_assign: pd.DataFrame, cluster_col: str, feature_cols: list) -> tuple:
    summary = (
        df_assign.groupby(cluster_col)
        .agg(
            n_subjects=("subject_id", "count"),
            n_case=(   "label", lambda x: int((x == "case").sum())),
            n_control=("label", lambda x: int((x == "control").sum())),
        )
        .reset_index()
    )
    summary["case_fraction"]    = summary["n_case"]    / summary["n_subjects"]
    summary["control_fraction"] = summary["n_control"] / summary["n_subjects"]

    total_case    = int((df_assign["label"] == "case").sum())
    total_control = int((df_assign["label"] == "control").sum())

    enrich_rows = []
    for cl in sorted(df_assign[cluster_col].unique()):
        sub = df_assign[df_assign[cluster_col] == cl]
        a   = int((sub["label"] == "case").sum())
        b   = int((sub["label"] == "control").sum())
        c   = total_case    - a
        d   = total_control - b
        odds_ratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        enrich_rows.append({
            cluster_col: cl,
            "cluster_case":    a, "cluster_control": b,
            "outside_case":    c, "outside_control": d,
            "odds_ratio":      float(odds_ratio),
            "fisher_pvalue":   float(pvalue),
        })

    centroids = df_assign.groupby(cluster_col)[feature_cols].mean().reset_index()
    enrich    = pd.DataFrame(enrich_rows)
    merged    = summary.merge(enrich, on=cluster_col, how="left")
    return merged, enrich, centroids


def main():
    df           = pd.read_csv(MATRIX_FILE, sep="\t")
    feature_cols = [c for c in df.columns if c not in META_COLS]
    Xz           = StandardScaler().fit_transform(df[feature_cols].values)

    # KMeans grid search
    grid_rows, models = [], {}
    for k in [2, 3, 4, 5, 6]:
        km     = KMeans(n_clusters=k, random_state=42, n_init=50)
        labels = km.fit_predict(Xz)
        sil    = silhouette_score(Xz, labels)
        grid_rows.append({"k": k, "inertia": float(km.inertia_), "silhouette": float(sil)})
        models[k] = (km, labels)

    kgrid  = pd.DataFrame(grid_rows).sort_values("k")
    best_k = int(kgrid.sort_values(["silhouette", "k"], ascending=[False, True]).iloc[0]["k"])
    kgrid.to_csv("results/evidence_modules_v3_kmeans_k_grid.tsv", sep="\t", index=False)
    print(f"KMeans grid:\n{kgrid.to_string(index=False)}\nSelected k={best_k}")

    # KMeans with best k
    km, km_labels = models[best_k]
    df_km         = df.copy()
    df_km["cluster"] = km_labels.astype(int)
    df_km = df_km.sort_values(["cluster", "label", "subject_id"]).reset_index(drop=True)
    km_summary, km_enrich, km_centroids = summarize_clusters(df_km, "cluster", feature_cols)

    df_km.to_csv(        "results/evidence_modules_v3_kmeans_cluster_assignment.tsv", sep="\t", index=False)
    km_summary.to_csv(   "results/evidence_modules_v3_kmeans_cluster_summary.tsv",    sep="\t", index=False)
    km_enrich.to_csv(    "results/evidence_modules_v3_kmeans_cluster_enrichment.tsv", sep="\t", index=False)
    km_centroids.to_csv( "results/evidence_modules_v3_kmeans_cluster_module_centroids.tsv", sep="\t", index=False)

    # Hierarchical (same k)
    hc        = AgglomerativeClustering(n_clusters=best_k, linkage="ward")
    hc_labels = hc.fit_predict(Xz)
    df_hc     = df.copy()
    df_hc["cluster"] = hc_labels.astype(int)
    df_hc = df_hc.sort_values(["cluster", "label", "subject_id"]).reset_index(drop=True)
    hc_summary, hc_enrich, hc_centroids = summarize_clusters(df_hc, "cluster", feature_cols)

    df_hc.to_csv(        "results/evidence_modules_v3_hclust_cluster_assignment.tsv", sep="\t", index=False)
    hc_summary.to_csv(   "results/evidence_modules_v3_hclust_cluster_summary.tsv",    sep="\t", index=False)
    hc_enrich.to_csv(    "results/evidence_modules_v3_hclust_cluster_enrichment.tsv", sep="\t", index=False)
    hc_centroids.to_csv( "results/evidence_modules_v3_hclust_cluster_module_centroids.tsv", sep="\t", index=False)

    print("\nKMeans cluster summary:")
    print(km_summary.to_string(index=False))
    print("\nHierarchical cluster summary:")
    print(hc_summary.to_string(index=False))


if __name__ == "__main__":
    main()
