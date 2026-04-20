"""
11_plot_neighbor_purity_bar.py
Bar chart of mean nearest-neighbor label purity at k=5.
Dashed line marks the expected value under random mixing (~0.50).
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

os.environ["PANDAS_USE_NUMEXPR"] = "0"


def main():
    df      = pd.read_csv("results/evidence_modules_v3_neighbor_labelmix_summary.tsv", sep="\t")
    plot_df = df[df["k"] == 5].copy()
    plot_df["label_cap"] = plot_df["label"].str.capitalize()
    plot_df = plot_df.sort_values("label_cap")

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    bars = ax.bar(
        plot_df["label_cap"],
        plot_df["mean_same_label_fraction"],
        color=["#7c746a", "#b0a48f"],
        width=0.55,
    )
    ax.axhline(0.5, linestyle="--", linewidth=1.6, color="#9b4343")
    ax.text(0.82, 0.505, "random mixing ~0.50",
            color="#9b4343", fontsize=12, transform=ax.get_yaxis_transform())
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("mean same-label neighbor proportion\n(k = 5)", fontsize=14)
    ax.set_title("Nearest-neighbor purity", fontsize=18, pad=12)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=12)
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar, val in zip(bars, plot_df["mean_same_label_fraction"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val + 0.02, f"{val:.2f}",
            ha="center", va="bottom", fontsize=16, fontweight="bold",
        )
    plt.tight_layout()
    plt.savefig("results/evidence_modules_v3_neighbor_purity_k5.png",
                dpi=300, bbox_inches="tight")
    plt.close()
    print("Wrote results/evidence_modules_v3_neighbor_purity_k5.png")


if __name__ == "__main__":
    main()
