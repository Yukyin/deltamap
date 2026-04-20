"""
06_build_candidate_modules.py
Expand biologically motivated seed genes into candidate evidence modules
using cross-subject correlation in the paired-delta space.

Approach:
  1. Filter noise genes (MT-, RPL, RPS, globin, lncRNA)
  2. Z-score the top 300 high-variance genes across subjects
  3. For each seed gene, compute Pearson correlation with all other genes
  4. Take the top-15 positively and negatively correlated genes per seed

Seed genes are chosen to span four hypothesized immune programs:
  - Myeloid inflammation: S100A8, S100A9, LYZ
  - APC / MHC-II: CD74, HLA-DRA, HLA-DPA1
  - Stress / immediate-early: NFKBIA, TXNIP, FOS, JUNB
  - Humoral B-cell: IGKC, IGHM, JCHAIN
"""
import os
import re
import numpy as np
import pandas as pd

os.environ["PANDAS_USE_NUMEXPR"] = "0"

DELTA_FILE = "results/pseudobulk_delta.tsv"
OUT_VAR    = "results/delta_gene_variability.tsv"
OUT_CORR   = "results/delta_gene_seed_correlations.tsv"
OUT_MODS   = "results/delta_candidate_modules.tsv"

TOP_VARIABLE_GENES = 300
TOP_CORRELATED     = 15

NOISE_PATTERNS = [
    r"^MT-", r"^RPL", r"^RPS", r"^HBA", r"^HBB$", r"^HBD$",
    r"^HBG", r"^HBM$", r"^MALAT1$", r"^NEAT1$", r"^TMSB10$", r"^TMSB4X$",
]

SEED_GENES = [
    "S100A8", "S100A9", "LYZ",
    "CD74", "HLA-DRA", "HLA-DPA1",
    "NFKBIA", "TXNIP", "FOS", "JUNB",
    "IGKC", "IGHM", "JCHAIN",
]


def is_noise(gene: str) -> bool:
    return any(re.match(p, str(gene)) for p in NOISE_PATTERNS)


def main():
    delta = pd.read_csv(DELTA_FILE, sep="\t").set_index("gene")

    # Remove noise genes
    keep_mask = ~pd.Series(delta.index, index=delta.index).apply(is_noise)
    delta_f   = delta.loc[keep_mask].copy()

    # Variability ranking
    var_df = pd.DataFrame({
        "gene": delta_f.index,
        "mean_delta": delta_f.mean(axis=1).values,
        "std_delta": delta_f.std(axis=1).values,
        "mean_abs_delta": delta_f.abs().mean(axis=1).values,
    }).sort_values(["std_delta", "mean_abs_delta"], ascending=False)
    var_df.to_csv(OUT_VAR, sep="\t", index=False)

    # Z-score top variable genes
    top_genes = var_df.head(TOP_VARIABLE_GENES)["gene"].tolist()
    X  = delta_f.loc[top_genes].copy()
    Xz = (
        X.sub(X.mean(axis=1), axis=0)
         .div(X.std(axis=1).replace(0, np.nan), axis=0)
         .fillna(0.0)
    )

    # Seed gene correlation
    seeds_present = [g for g in SEED_GENES if g in Xz.index]
    if len(seeds_present) < len(SEED_GENES):
        missing = set(SEED_GENES) - set(seeds_present)
        print(f"Warning: seed genes not in top-variable set: {missing}")

    corr_rows = []
    for seed in seeds_present:
        s = Xz.loc[seed].values
        for gene in Xz.index:
            r = np.corrcoef(s, Xz.loc[gene].values)[0, 1]
            corr_rows.append({"seed_gene": seed, "gene": gene, "pearson_r": float(r)})

    corr_df = (
        pd.DataFrame(corr_rows)
        .sort_values(["seed_gene", "pearson_r"], ascending=[True, False])
    )
    corr_df.to_csv(OUT_CORR, sep="\t", index=False)

    # Build candidate module gene lists
    module_rows = []
    for seed in seeds_present:
        sub = corr_df[corr_df["seed_gene"] == seed].copy()
        sub = sub[sub["gene"] != seed]

        # Top positive correlates
        for rank, (_, row) in enumerate(sub.head(TOP_CORRELATED).iterrows(), start=1):
            module_rows.append({
                "seed_gene": seed, "direction": "positive",
                "rank": rank, "gene": row["gene"], "pearson_r": row["pearson_r"],
            })
        # Top negative correlates
        for rank, (_, row) in enumerate(
            sub.sort_values("pearson_r", ascending=True).head(TOP_CORRELATED).iterrows(),
            start=1
        ):
            module_rows.append({
                "seed_gene": seed, "direction": "negative",
                "rank": rank, "gene": row["gene"], "pearson_r": row["pearson_r"],
            })

    pd.DataFrame(module_rows).to_csv(OUT_MODS, sep="\t", index=False)

    print(f"Variability file: {OUT_VAR}")
    print(f"Seed correlations: {OUT_CORR}  ({len(corr_df)} rows)")
    print(f"Candidate modules: {OUT_MODS}  ({len(module_rows)} entries)")


if __name__ == "__main__":
    main()
