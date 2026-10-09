import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ---- STRAND ALIGNMENT FUNCTIONS ----


# Function to plot strang alignment matrics
def plot_strand_alignment_metrics(summary_df, output_path):

    plot_df = summary_df.copy()

    # Format to plot
    plot_df["VARIANTS_RAW"] = pd.to_numeric(plot_df["VARIANTS_RAW"], errors="coerce")
    plot_df["Reference Match"] = pd.to_numeric(plot_df["REF_MATCH"], errors="coerce")
    plot_df = (
        plot_df
        .sort_values("Reference Match", ascending=False)
        .reset_index(drop=True)
    )

    # Calculate percentage for absolute value metrics
    # (Reference Match is alreaady a percentage)
    count_columns = {
        "REF_FLIP_VARIANT_NR": "Flipped Variants",
        "REF_AMBIGUOUS_INDEL_NR": "Ambiguous Indels",
    }
    for source_column, label in count_columns.items():
        plot_df[source_column] = pd.to_numeric(plot_df[source_column], errors="coerce")
        plot_df[label] = np.where(plot_df["VARIANTS_RAW"].gt(0),(
            plot_df[source_column]
            .fillna(0)
            .div(plot_df["VARIANTS_RAW"])
            .mul(100)),
            np.nan
        )

    # Plot settings
    plot_columns = [
        "Reference Match",
        "Flipped Variants",
        "Ambiguous Indels",
    ]
    acadia = [
        "#A4BED5",
        "#FED789",
        "#72874E",
        "#023743",
        "#476F84",
        "#453947",
    ]
    metric_colors = dict(zip(plot_columns, acadia))
    fig, ax = plt.subplots(figsize=(14, 6))
    x_positions = np.arange(len(plot_df))
    bar_width = 0.25
    offsets = (np.arange(len(plot_columns)) - (len(plot_columns) - 1) / 2) * bar_width

    # Plot offset bars
    for offset, label in zip(offsets, plot_columns):
        values = plot_df[label].to_numpy(dtype=float)
        bars = ax.bar(
            x_positions + offset,
            values,
            width=bar_width,
            color=metric_colors[label],
            edgecolor="black",
            linewidth=0.5,
            label=label,
        )

        # Do not label Reference Match
        # (Reference Match is already a percentage)
        if label == "Reference Match":
            continue

        # Absolute counts above the other bars
        source_column = next(
            column
            for column, display_label in count_columns.items()
            if display_label == label
        )
        absolute_counts = plot_df[source_column]
        bar_labels = [
            f"{int(count):,}" if pd.notna(count) else ""
            for count in absolute_counts
        ]
        ax.bar_label(
            bars,
            labels=bar_labels,
            padding=3,
            fontsize=8,
            rotation=90,
        )

    # Labels and formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(plot_df["COHORT"], rotation=45, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("Strand alignment metric (%)")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)

    maximum = plot_df[plot_columns].max().max()
    if pd.notna(maximum):
        ax.set_ylim(0, max(maximum * 1.15, 1))

    # Save plot
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
