"""
10_plot_module_effect_bar.py
Horizontal bar chart of case − control difference per evidence module.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

os.environ["PANDAS_USE_NUMEXPR"] = "0"

EFFECTS_FILE = "results/evidence_modules_v3.no_lowcell_effects.tsv"

MODULE_ORDER = [
    "humoral_bcell_v3",
    "myeloid_inflammation_v3",
    "apc_mhcII_v3",
    "stress_response_v3",
]
LABEL_MAP = {
    "humoral_bcell_v3":       "Humoral\nB-cell",
    "myeloid_inflammation_v3":"Myeloid\ninflammation",
    "apc_mhcII_v3":           "APC /\nMHC-II",
    "stress_response_v3":     "Stress\nresponse",
}
COLORS = {
    "humoral_bcell_v3":       "#96413f",
    "myeloid_inflammation_v3":"#627f69",
    "apc_mhcII_v3":           "#607b92",
    "stress_response_v3":     "#7d7d7d",
}


def main():
    df      = pd.read_csv(EFFECTS_FILE, sep="\t").set_index("module")
    plot_df = df.loc[MODULE_ORDER].reset_index()
    plot_df["display_label"] = plot_df["module"].map(LABEL_MAP)

    fig, ax = plt.subplots(figsize=(10, 6.5))
    bars = ax.barh(
        plot_df["display_label"],
        plot_df["case_minus_control"],
        color=[COLORS[m] for m in plot_df["module"]],
        height=0.58,
    )
    ax.axvline(0, color="black", linewidth=1.2)
    ax.set_title("Case − control difference in evidence modules", fontsize=20, pad=14)
    ax.set_xlabel("case minus control (paired delta module score)", fontsize=16)
    ax.tick_params(axis="y", labelsize=16, pad=10)
    ax.tick_params(axis="x", labelsize=13)
    ax.grid(axis="x", linestyle="-", alpha=0.22)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    x_min = plot_df["case_minus_control"].min()
    x_max = max(0, plot_df["case_minus_control"].max())
    span  = x_max - x_min
    ax.set_xlim(x_min - 0.12 * span, x_max + 0.03 * span)

    for bar, val in zip(bars, plot_df["case_minus_control"]):
        y  = bar.get_y() + bar.get_height() / 2
        x  = val - 0.10 * span if val < 0 else val + 0.02 * span
        ha = "right" if val < 0 else "left"
        ax.text(x, y, f"{val:.1f}", va="center", ha=ha,
                fontsize=14, color="black", fontweight="bold")

    plt.tight_layout()
    plt.savefig("results/evidence_modules_v3_case_minus_control_barplot.png",
                dpi=300, bbox_inches="tight")
    plt.close()
    print("Wrote results/evidence_modules_v3_case_minus_control_barplot.png")


if __name__ == "__main__":
    main()
