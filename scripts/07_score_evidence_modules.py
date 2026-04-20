"""
07_score_evidence_modules.py
Score each subject on the four champion evidence modules using mean delta
expression of curated gene lists. Produces a subject × module score matrix
(the core evidence representation), plus summary statistics and inter-module
correlations. Outputs both full-cohort and low-cell-count-excluded versions.

Module definitions (v3 champion):
  - myeloid_inflammation: S100A8/9/12, VCAN, LYZ, FCN1, CD14, TYROBP, CTSD, ANXA2, PSAP, MNDA
  - apc_mhcII: CD74, HLA-DQA1/DQB1/DPA1/DRA/DRB1/DPB1, CST3, CD37
  - stress_response: NFKBIA/Z, PPP1R15A, KLF6, ZFP36, BTG2, FOSB, JUND, DUSP1, GADD45B
  - humoral_bcell: IGKC, JCHAIN, IGHA1, IGHM, IGLC1/2/3
"""
import os
import pandas as pd

os.environ["PANDAS_USE_NUMEXPR"] = "0"

DELTA_FILE = "results/pseudobulk_delta.tsv"
META_FILE  = "results/pseudobulk_delta_metadata.tsv"

EVIDENCE_MODULES = {
    "myeloid_inflammation_v3": [
        "S100A8", "S100A9", "S100A12", "VCAN", "LYZ",
        "FCN1", "CD14", "TYROBP", "CTSD", "ANXA2", "PSAP", "MNDA",
    ],
    "apc_mhcII_v3": [
        "CD74", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1",
        "HLA-DRA", "HLA-DRB1", "HLA-DPB1", "CST3", "CD37",
    ],
    "stress_response_v3": [
        "NFKBIA", "NFKBIZ", "PPP1R15A", "KLF6", "ZFP36",
        "BTG2", "FOSB", "JUND", "DUSP1", "GADD45B",
    ],
    "humoral_bcell_v3": [
        "IGKC", "JCHAIN", "IGHA1", "IGHM",
        "IGLC1", "IGLC2", "IGLC3",
    ],
}

META_COLS = ["subject_id", "label", "sex", "age_bin", "has_low_cell_count_member"]


def score_modules(delta: pd.DataFrame, modules: dict) -> tuple[list, list]:
    """Return (definition_rows, score_rows)."""
    def_rows, score_rows = [], []
    for module, genes in modules.items():
        present = [g for g in genes if g in delta.index]
        missing = set(genes) - set(present)
        if missing:
            print(f"  [{module}] genes not found: {missing}")
        for rank, g in enumerate(present, start=1):
            def_rows.append({"module": module, "rank": rank, "gene": g})
        score = delta.loc[present].mean(axis=0)
        for sample_name, val in score.items():
            score_rows.append({
                "module": module,
                "delta_name": sample_name,
                "n_genes_used": len(present),
                "module_score": float(val),
            })
    return def_rows, score_rows


def write_outputs(score_df: pd.DataFrame, prefix: str) -> None:
    features = [c for c in score_df.columns if c not in META_COLS + ["delta_name", "module",
                                                                       "n_genes_used", "module_score"]]
    matrix = (
        score_df
        .pivot_table(
            index=META_COLS,
            columns="module",
            values="module_score",
        )
        .reset_index()
    )
    matrix.columns.name = None

    summary = (
        score_df
        .groupby(["module", "label"])["module_score"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )

    p = summary.pivot(index="module", columns="label", values="mean").reset_index()
    p["case_minus_control"] = p["case"] - p["control"]
    m = summary.pivot(index="module", columns="label", values="median").reset_index()
    eff = p.merge(m, on="module", suffixes=("_mean", "_median"))
    n1 = summary[summary.label == "case"][["module", "n"]].rename(columns={"n": "n_case"})
    n0 = summary[summary.label == "control"][["module", "n"]].rename(columns={"n": "n_control"})
    eff = eff.merge(n1, on="module").merge(n0, on="module")
    effect_cols = ["module", "case_mean", "control_mean", "case_minus_control",
                   "case_median", "control_median", "n_case", "n_control"]

    feat_cols = [c for c in matrix.columns if c not in META_COLS]
    corr = matrix[feat_cols].corr()
    corr_rows = [
        {"module1": r, "module2": c, "pearson_r": corr.loc[r, c]}
        for r in corr.index for c in corr.columns
    ]

    score_df.to_csv(f"results/{prefix}_scores.tsv", sep="\t", index=False)
    matrix.to_csv(f"results/{prefix}_matrix.tsv", sep="\t", index=False)
    summary.to_csv(f"results/{prefix}_summary.tsv", sep="\t", index=False)
    eff[effect_cols].sort_values("case_minus_control", ascending=False).to_csv(
        f"results/{prefix}_effects.tsv", sep="\t", index=False
    )
    pd.DataFrame(corr_rows).to_csv(f"results/{prefix}_correlation.tsv", sep="\t", index=False)


def main():
    delta = pd.read_csv(DELTA_FILE, sep="\t").set_index("gene")
    meta  = pd.read_csv(META_FILE, sep="\t")

    # Write module definitions
    def_rows, score_rows = score_modules(delta, EVIDENCE_MODULES)
    pd.DataFrame(def_rows).to_csv(
        "results/evidence_modules_v3_definition.tsv", sep="\t", index=False
    )

    score_df = (
        pd.DataFrame(score_rows)
        .merge(meta, on="delta_name", how="left")
    )
    score_df.to_csv("results/evidence_modules_v3_scores.tsv", sep="\t", index=False)

    # Full cohort outputs
    write_outputs(score_df, "evidence_modules_v3")

    # Exclude low-cell-count samples
    score_clean = score_df[score_df["has_low_cell_count_member"] != True].copy()
    score_clean.to_csv("results/evidence_modules_v3_scores.no_lowcell.tsv", sep="\t", index=False)
    write_outputs(score_clean, "evidence_modules_v3.no_lowcell")

    # Summary for quick inspection
    matrix = pd.read_csv("results/evidence_modules_v3.no_lowcell_matrix.tsv", sep="\t")
    print(f"\nEvidence matrix (no_lowcell): {matrix.shape}")
    effects = pd.read_csv("results/evidence_modules_v3.no_lowcell_effects.tsv", sep="\t")
    print("\nCase − control effect per module:")
    print(effects[["module", "case_minus_control", "n_case", "n_control"]].to_string(index=False))


if __name__ == "__main__":
    main()
