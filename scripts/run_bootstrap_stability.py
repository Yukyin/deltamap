"""
run_bootstrap_stability.py
Bootstrap stability analysis for evidence module rankings.
Repeatedly resamples subjects (with replacement) and re-ranks by module score.
Outputs per-subject rank mean, std, and quantiles to assess stability.
"""
import os
import numpy as np
import pandas as pd

os.environ["PANDAS_USE_NUMEXPR"] = "0"

MATRIX_FILE  = "results/evidence_modules_v3.no_lowcell_matrix.tsv"
META_COLS    = ["subject_id", "label", "sex", "age_bin", "has_low_cell_count_member"]
N_BOOTSTRAP  = 500
RANDOM_SEED  = 42


def bootstrap_ranks(scores: np.ndarray, n_boot: int, rng: np.random.Generator) -> np.ndarray:
    """Return (n_subjects × n_boot) rank matrix."""
    n        = len(scores)
    rank_mat = np.zeros((n, n_boot), dtype=float)
    for b in range(n_boot):
        idx          = rng.choice(n, size=n, replace=True)
        boot_scores  = scores[idx]
        # Rank in full (non-bootstrap) space using bootstrap-estimated mean order
        order        = np.argsort(boot_scores)
        temp         = np.empty_like(order)
        temp[order]  = np.arange(n)
        rank_mat[idx, b] = temp
    return rank_mat


def main():
    df           = pd.read_csv(MATRIX_FILE, sep="\t")
    feature_cols = [c for c in df.columns if c not in META_COLS]
    rng          = np.random.default_rng(RANDOM_SEED)

    all_rows = []
    for module in feature_cols:
        scores   = df[module].values.astype(float)
        n        = len(scores)
        rank_mat = np.zeros((n, N_BOOTSTRAP), dtype=float)

        for b in range(N_BOOTSTRAP):
            idx         = rng.choice(n, size=n, replace=True)
            boot_scores = scores[idx]
            order       = np.argsort(boot_scores)
            ranks       = np.empty_like(order, dtype=float)
            ranks[order] = np.arange(1, n + 1, dtype=float)
            rank_mat[:, b] = ranks

        for i in range(n):
            r = rank_mat[i]
            all_rows.append({
                "feature":      module,
                "delta_name":   df.loc[i, "subject_id"] + "_DELTA",
                "subject_id":   df.loc[i, "subject_id"],
                "label":        df.loc[i, "label"],
                "rank_mean":    float(np.mean(r)),
                "rank_std":     float(np.std(r, ddof=1)),
                "rank_q05":     float(np.percentile(r, 5)),
                "rank_q25":     float(np.percentile(r, 25)),
                "rank_median":  float(np.median(r)),
                "rank_q75":     float(np.percentile(r, 75)),
                "rank_q95":     float(np.percentile(r, 95)),
            })

    result = pd.DataFrame(all_rows)
    result.to_csv(
        "results/evidence_modules_champion_v1_bootstrap_rank_summary.tsv",
        sep="\t", index=False,
    )
    print(f"Bootstrap stability complete ({N_BOOTSTRAP} iterations, {len(feature_cols)} modules)")
    print("Wrote results/evidence_modules_champion_v1_bootstrap_rank_summary.tsv")

    # Quick summary: mean rank_std per module
    print("\nMean rank std per module (lower = more stable):")
    print(
        result.groupby("feature")["rank_std"]
        .mean()
        .sort_values()
        .to_string()
    )


if __name__ == "__main__":
    main()
