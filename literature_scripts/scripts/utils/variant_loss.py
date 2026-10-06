import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Patch


# ---- VARIANT LOSS FUNCTIONS ----


# Function to check for the reason of variant loss during harmonization
def check_variant_loss(summary_df):

    loss_columns = [
        "VARIANTS_LOSS_BADALLELES",
        "VARIANTS_LOSS_BADSTATS",
        "VARIANTS_LIFTOVER_UNMAPPED",
    ]

    # Test whether total variant loss is the sum of variants remove due to:
    # - bad alleles
    # - bad statistics
    # - unmapped variants during liftover (for build 37 -> 38)
    summary_df["VARIANTS_LOSS_CALCULATED"] = summary_df[loss_columns].fillna(0).sum(axis=1)
    mismatch_mask = summary_df["VARIANTS_LOSS_TOTAL"]  != summary_df["VARIANTS_LOSS_CALCULATED"]

    if mismatch_mask.any():
        mismatch_df = summary_df.loc[
            mismatch_mask,
            ["COHORT", "VARIANTS_LOSS_TOTAL", "VARIANTS_LOSS_CALCULATED", *loss_columns]
        ]

        logging.warning("Total variant loss does not match the sum of its components for %d cohort(s):\n%s",
                        len(mismatch_df), mismatch_df.to_string(index=False))
    else:
        logging.info("Total variant loss matches bad alleles, bad statistics and unmapped variants removal.")


# Function to plot variant loss during harmonization
def plot_variant_loss(summary_df, output_path):

    # Calculate percentages for each removal reason
    loss_columns = {
        "VARIANTS_LOSS_BADALLELES": "Bad alleles",
        "VARIANTS_LOSS_BADSTATS": "Bad statistics",
        "VARIANTS_LIFTOVER_UNMAPPED": "Liftover unmapped",
    }
    plot_df = summary_df.copy()
    for source_column, label in loss_columns.items():
        plot_df[label] = (
            plot_df[source_column]
            .fillna(0)
            .div(plot_df["VARIANTS_RAW"])
            .mul(100)
        )
    plot_columns = list(loss_columns.values())

    # Prepare stacked bar
    plot_df["VARIANT_LOSS_PCT"] = plot_df[plot_columns].sum(axis=1)
    plot_df["VARIANT_LOSS_NR"] = (
        plot_df[list(loss_columns)]
        .fillna(0)
        .sum(axis=1)
        .astype(int)
    )
    plot_df = (
        plot_df
        .sort_values("VARIANT_LOSS_PCT", ascending=False)
        .reset_index(drop=True)
    )

    # Plot settings
    acadia = [
        "#A4BED5",
        "#FED789",
        "#72874E",
        "#023743",
        "#476F84",
        "#453947",
    ]
    loss_colors = dict(zip(["Bad alleles", "Bad statistics", "Liftover unmapped"], acadia[:3]))
    fig, ax = plt.subplots(figsize=(14, 6))
    x_positions = np.arange(len(plot_df))
    bar_bottom = np.zeros(len(plot_df))

    # Plot stacked bars
    for label in plot_columns:
        values = plot_df[label].to_numpy()
        ax.bar(
            x_positions,
            values,
            bottom=bar_bottom,
            color=loss_colors[label],
            edgecolor="black",
            linewidth=0.5,
            label=label,
        )
        bar_bottom += values

    # Absolute total variant loss above bars
    for x, percentage, count in zip(x_positions, plot_df["VARIANT_LOSS_PCT"], plot_df["VARIANT_LOSS_NR"]):
        ax.annotate(
            f"{count:,}",
            xy=(x, percentage),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    # Labels and formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(plot_df["COHORT"], rotation=45, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("Variant loss (%)")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)

    maximum = plot_df["VARIANT_LOSS_PCT"].max()
    if pd.notna(maximum):
        ax.set_ylim(0, max(maximum * 1.15, 1))

    # Legend
    legend_handles = [
        Patch(facecolor=color, edgecolor="black", label=label)
        for label, color in loss_colors.items()
    ]
    ax.legend(handles=legend_handles, title="Variants Removed by")

    # Save plot
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)