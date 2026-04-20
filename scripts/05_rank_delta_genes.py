"""
05_rank_delta_genes.py
Rank genes by the magnitude of the case-minus-control difference in paired delta.
Also produces a clean version with mitochondrial, ribosomal, globin, and high-noise
genes removed — these dominate variance but are uninformative for immune programs.
"""
import os
import re
import numpy as np
import pandas as pd

os.environ["PANDAS_USE_NUMEXPR"] = "0"

DELTA_FILE = "results/pseudobulk_delta.tsv"
META_FILE  = "results/pseudobulk_delta_metadata.tsv"
OUT_RANK   = "results/delta_gene_ranking.no_lowcell.tsv"
OUT_CLEAN  = "results/delta_gene_ranking.no_lowcell.clean.tsv"

# Genes excluded from evidence module construction
NOISE_PATTERNS = [
    r"^MT-",       # mitochondrial
    r"^RPL",       # ribosomal large
    r"^RPS",       # ribosomal small
    r"^HBA",       # hemoglobin alpha
    r"^HBB$",      # hemoglobin beta
    r"^HBD$",
    r"^HBG",
    r"^HBM$",
    r"^MALAT1$",   # lncRNA, high abundance
    r"^NEAT1$",    # lncRNA
    r"^TMSB10$",   # cytoskeletal, high abundance
    r"^TMSB4X$",
]


def is_noise_gene(gene: str) -> bool:
    return any(re.match(p, str(gene)) for p in NOISE_PATTERNS)


def main():
    delta = pd.read_csv(DELTA_FILE, sep="\t").set_index("gene")
    meta  = pd.read_csv(META_FILE, sep="\t")

    # Exclude low-cell-count samples
    clean_meta = meta[meta["has_low_cell_count_member"] != True].copy()
    case_cols  = [c for c in clean_meta.loc[clean_meta.label == "case", "delta_name"]
                  if c in delta.columns]
    ctrl_cols  = [c for c in clean_meta.loc[clean_meta.label == "control", "delta_name"]
                  if c in delta.columns]

    print(f"Case subjects: {len(case_cols)}, Control subjects: {len(ctrl_cols)}")

    case_mean = delta[case_cols].mean(axis=1)
    ctrl_mean = delta[ctrl_cols].mean(axis=1)
    case_med  = delta[case_cols].median(axis=1)
    ctrl_med  = delta[ctrl_cols].median(axis=1)
    effect    = case_mean - ctrl_mean

    ranking = pd.DataFrame({
        "gene": delta.index,
        "case_mean_delta": case_mean.values,
        "control_mean_delta": ctrl_mean.values,
        "case_median_delta": case_med.values,
        "control_median_delta": ctrl_med.values,
        "delta_diff_case_minus_control": effect.values,
        "abs_delta_diff": np.abs(effect.values),
    }).sort_values("abs_delta_diff", ascending=False)

    ranking.to_csv(OUT_RANK, sep="\t", index=False)

    # Clean version: remove noise genes
    ranking["filtered_out"] = ranking["gene"].apply(is_noise_gene)
    ranking[ranking["filtered_out"] == False].to_csv(OUT_CLEAN, sep="\t", index=False)

    print(f"Total ranked genes: {len(ranking)}")
    print(f"After noise filtering: {(~ranking['filtered_out']).sum()}")
    print(f"Wrote {OUT_RANK}")
    print(f"Wrote {OUT_CLEAN}")
    print("\nTop 20 case-minus-control genes:")
    print(ranking[~ranking["filtered_out"]].head(20)[
        ["gene", "delta_diff_case_minus_control", "abs_delta_diff"]
    ].to_string(index=False))


if __name__ == "__main__":
    main()
