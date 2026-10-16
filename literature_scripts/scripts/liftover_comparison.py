import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils.ncbi_validation import add_ncbi_validation

from matplotlib.patches import Patch
from pathlib import Path
from paths import PathManager


KEY_COLUMNS = ["SEQID", "UNIPROT", "rsID", "CHR", "POS37"]
#KEY_COLUMNS = ["SEQID", "OLINKID", "rsID", "CHR", "POS37"] #TODO Use OlinkID after clean update
COMPARISON_COLUMNS = KEY_COLUMNS + ["POS"]
METHODS = (
    ("bcftools", "literature_table_harmonized_bcftools", "#A4BED5"),
    ("GWASLab (bridge)", "literature_table_harmonized_gwaslab_bridge", "#FED789"),
    ("GWASLab (standard)", "literature_table_harmonized_gwaslab_standard", "#72874E"),
)


# ---- HELPER FUNCTIONS ----

# Function to load and format harmonization summary
def load_summary(directory: Path, method: str) -> pd.DataFrame:
    path = directory / "harmonization_summary.tsv"
    df = pd.read_csv(path, sep="\t", usecols=["COHORT", "VARIANTS_RAW", "VARIANTS_LIFTOVER_UNMAPPED"])
    df["VARIANTS_RAW"] = pd.to_numeric(df["VARIANTS_RAW"], errors="raise")
    df["VARIANTS_LIFTOVER_UNMAPPED"] = pd.to_numeric(df["VARIANTS_LIFTOVER_UNMAPPED"], errors="raise")
    df["METHOD"] = method
    df["UNMAPPED_PCT"] = df["VARIANTS_LIFTOVER_UNMAPPED"] / df["VARIANTS_RAW"] * 100
    return df


# Function to load harmonized and mapped summary statistics
def load_sumstats(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", 
                     usecols=COMPARISON_COLUMNS,
                     dtype={"SEQID": "string", "rsID": "string", "CHR": "string"}
    )
    for col in ("POS37", "POS"):
        df[col] = pd.to_numeric(df[col], errors="raise").astype("Int64")
    #if df[KEY_COLUMNS].duplicated().any():
    #    examples = df.loc[df[KEY_COLUMNS].duplicated(keep=False), KEY_COLUMNS].head(5)
    #    raise ValueError(
    #        f"{path}: key {KEY_COLUMNS} is not unique; inspect:\n"
    #        f"{examples.to_string(index=False)}"
    #    )
    return df


# Function to plot the number of unmapped variants per liftover method
def plot_unmapped(summaries: pd.DataFrame, output: Path) -> None:
    
    # Format to plot
    cohort_order = summaries.groupby("COHORT")["UNMAPPED_PCT"].max().sort_values(ascending=False).index.tolist()
    x_positions = np.arange(len(cohort_order))
    width = 0.26
    maximum = pd.to_numeric(summaries["UNMAPPED_PCT"], errors="coerce").max()
    y_maximum = (max(maximum * 1.30, 1) if pd.notna(maximum) else 1)
    minimum_label_height = y_maximum * 0.02
    fig, ax = plt.subplots(figsize=(max(14, len(cohort_order) * 0.65), 6))

    for index, (method, _, color) in enumerate(METHODS):
        method_df = summaries.loc[summaries["METHOD"].eq(method)].set_index("COHORT")
        percentages = method_df["UNMAPPED_PCT"].reindex(cohort_order)
        absolute_counts = method_df["VARIANTS_LIFTOVER_UNMAPPED"].reindex(cohort_order)
        bar_positions = x_positions + (index - (len(METHODS) - 1) / 2) * width

        # Plot bars
        ax.bar(
            bar_positions,
            percentages,
            width=width,
            label=method,
            color=color,
            edgecolor="black",
            linewidth=0.5,
        )

        # Absolute values above bars
        for (bar_x, percentage, count) in zip(bar_positions, percentages, absolute_counts):
            count = 0 if pd.isna(count) else count
            percentage = (float(percentage) if pd.notna(percentage) else 0)
            label_height = max(percentage, minimum_label_height)
            ax.annotate(
                f"{int(count):,}",
                xy=(bar_x, label_height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
                color=color,
                rotation=90,
            )

    # Labels and formatting
    ax.set_xticks(x_positions, cohort_order, rotation=45, ha="right")
    ax.set_ylabel("Unmapped variants (%)")
    ax.set_xlabel("")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)
    ax.set_ylim(0, y_maximum)
    ax.legend(title="", loc="upper right")
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


# Function to plot mismatched mapped variants (bcftools VS GWASLab bridge)
# plus NCBI validation
def plot_mismatches(summary: pd.DataFrame, output: Path) -> None:

    count_columns = [
        "POSITION_MISMATCHES",
        "NCBI_BCFTOOLS_MATCH",
        "NCBI_GWASLAB_BRIDGE_MATCH",
    ]
    match_columns = {
        "NCBI_NOT_FOUND": "Not Found",
        "NCBI_GWASLAB_BRIDGE_MATCH": "GWASLab bridge",
        "NCBI_BCFTOOLS_MATCH": "bcftools",
    }

    # Format to calculate percentage of total mismatches and NCBI matches
    plot_df = summary.copy()
    plot_df[count_columns] = plot_df[count_columns].fillna(0)
    plot_df["NCBI_NOT_FOUND"] = plot_df["POSITION_MISMATCHES"] - plot_df["NCBI_BCFTOOLS_MATCH"] - plot_df["NCBI_GWASLAB_BRIDGE_MATCH"]
    for source_column, label in match_columns.items():
        plot_df[label] = plot_df[source_column].div(plot_df["TOTAL_VARIANTS"]).mul(100)
    plot_df["POSITION_MISMATCH_PCT"] = plot_df["POSITION_MISMATCHES"].div(plot_df["TOTAL_VARIANTS"]).mul(100)
    plot_df = plot_df.sort_values("POSITION_MISMATCH_PCT", ascending=False, kind="stable").reset_index(drop=True)

    # Plot settings
    match_colors = {
        "Not Found": "#A4BED5",
        "GWASLab bridge": "#FED789",
        "bcftools": "#72874E",
    }
    fig, ax = plt.subplots(figsize=(14, 6))
    x_positions = np.arange(len(plot_df))
    bar_bottom = np.zeros(len(plot_df))

    # Plot stacked bars
    for label in match_columns.values():
        values = plot_df[label].to_numpy()
        ax.bar(
            x_positions,
            values,
            bottom=bar_bottom,
            color=match_colors[label],
            edgecolor="black",
            linewidth=0.5,
            label=label,
        )
        bar_bottom += values

    # Absolute mismatch value above bars
    for x, percentage, count in zip(x_positions, plot_df["POSITION_MISMATCH_PCT"], plot_df["POSITION_MISMATCHES"]):
        ax.annotate(
            f"{int(count):,}",
            xy=(x, percentage),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    # Labels and formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(plot_df["COHORT"], rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("POS38 mismatches (bcftools VS GWASLab bridge, %)")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)
    ax.set_ylim(0, max(plot_df["POSITION_MISMATCH_PCT"].max() * 1.2, 1))

    # Legend
    legend_handles = [
        Patch(facecolor=color, edgecolor="black", label=label)
        for label, color in match_colors.items()
    ]
    ax.legend(handles=legend_handles, title="NCBI validation")

    # Save plot
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


# Function to compare the mapped positions between bcftools and GWASLab bridge
def compare_positions(bcftools_dir: Path, bridge_dir: Path, bcftools_summary: pd.DataFrame, bridge_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    mismatch_frames = []
    bcftools_files = sorted(bcftools_dir.glob("*.gwaslab.liftover.tsv"))
    unmapped_bcftools = bcftools_summary.set_index("COHORT")["VARIANTS_LIFTOVER_UNMAPPED"]
    unmapped_bridge = bridge_summary.set_index("COHORT")["VARIANTS_LIFTOVER_UNMAPPED"]
    suffix = ".gwaslab.liftover.tsv"
    for bcftools_path in bcftools_files:

        # Get bcftools and GWASLab bridge sumstat pair
        cohort = bcftools_path.name.removesuffix(suffix)
        bridge_path = bridge_dir / f"{cohort}.gwaslab.tsv"
        bcftools = load_sumstats(bcftools_path)
        bridge = load_sumstats(bridge_path)
        compared = bcftools.merge(
            bridge,
            on=KEY_COLUMNS, 
            how="outer", 
            #validate="one_to_one", #TODO Uncomment after clean update
            suffixes=("_BCFTOOLS", "_BRIDGE"), indicator=True,
        )

        # Compare bcftools and GWASLab bridge
        shared = compared.loc[compared["_merge"].eq("both")]
        bcftools_position = shared["POS_BCFTOOLS"]
        bridge_position = shared["POS_BRIDGE"]

        # Both methods produced a mapped position, but positions differ
        position_mismatch = bcftools_position.notna() & bridge_position.notna() & bcftools_position.ne(bridge_position)

        # Key is shared, but only bcftools has a mapped position
        mapped_only_bcftools = bcftools_position.notna() & bridge_position.isna()

        # Key is shared, but only GWASLab has a mapped position
        mapped_only_bridge = bcftools_position.isna() & bridge_position.notna()

        # Both methods retained the key but neither has a position
        both_positions_missing = bcftools_position.isna() & bridge_position.isna()

        # Out mismatch summary report
        summary_rows.append({
            "COHORT": cohort,
            "TOTAL_VARIANTS": int(len(bcftools)),
            "UNMAPPED_VARIANTS_BCFTOOLS": int(unmapped_bcftools.loc[cohort]),
            "UNMAPPED_GWASLAB_BRIDGE": int(unmapped_bridge.loc[cohort]),
            "POSITION_MISMATCHES": int(position_mismatch.sum()),
            "SHARED_MAPPED_ONLY_BCFTOOLS": int(mapped_only_bcftools.sum()),
            "SHARED_MAPPED_ONLY_BRIDGE": int(mapped_only_bridge.sum()),
            "SHARED_BOTH_POSITIONS_MISSING": int(both_positions_missing.sum()),
            "ONLY_BCFTOOLS": int(compared["_merge"].eq("left_only").sum()),
            "ONLY_GWASLAB_BRIDGE": int(compared["_merge"].eq("right_only").sum())
        })
        summary_df = pd.DataFrame(summary_rows)

        # Out mismatched rowsß
        cohort_mismatches = shared.loc[position_mismatch].copy()

        if not cohort_mismatches.empty:
            cohort_mismatches.insert(0, "COHORT", cohort)
            cohort_mismatches.insert(1, "COMPARISON_STATUS", "POSITION_MISMATCH")
            cohort_mismatches = cohort_mismatches.rename(columns={
                "POS_BCFTOOLS": "POS38_BCFTOOLS",
                "POS_BRIDGE": "POS38_GWASLAB",
            })
            cohort_mismatches = cohort_mismatches.drop(columns="_merge")
            mismatch_frames.append(cohort_mismatches)

    if mismatch_frames:
        mismatch_rows_df = pd.concat(mismatch_frames, ignore_index=True)
    else:
        mismatch_rows_df = pd.DataFrame(columns=["COHORT", "COMPARISON_STATUS"] + KEY_COLUMNS + ["POS38_BCFTOOLS", "POS38_GWASLAB"])

    return (summary_df, mismatch_rows_df)



# ---- MAIN FUNCTION FOR LIFTOVER REPORT ----

# Report of the comparison between liftover methods by:
# the number of unmapped variants
# the number of variants whose mapped coordinates differ between bcftools and GWASLab bridge
def main() -> None:

    # ---- PATHS ----
    pm = PathManager()
    paths = pm.get_outputs()
    directories = {method: paths[key].parent for method, key, _ in METHODS}
    bcftools_dir = directories["bcftools"]
    bridge_dir = directories["GWASLab (bridge)"]
    liftover_dir = bcftools_dir.parent


    # ---- LOAD HARMONZIATION SUMMARIES ----
    summaries = pd.concat([load_summary(directories[method], method) for method, _, _ in METHODS], ignore_index=True)
    raw_totals = summaries.groupby("COHORT")["VARIANTS_RAW"].nunique(dropna=False)
    if raw_totals.gt(1).any():
        raise ValueError(f"VARIANTS_RAW differs between methods:\n{raw_totals[raw_totals.gt(1)]}")


    # ---- REPORT LIFTOVER COMPARISON ----
    plot_unmapped_path = liftover_dir / "liftover_unmapped.png"
    plot_mismatches_path = liftover_dir / "liftover_mismatches.png"
    mismatches_path = liftover_dir / "liftover_mismatches.tsv"
    mismatches_rows_path = liftover_dir / "liftover_mismatches_rows.tsv"

    mismatch_summary, mismatch_rows = compare_positions(bcftools_dir, bridge_dir, 
                                                        summaries.loc[summaries["METHOD"].eq("bcftools")],
                                                        summaries.loc[summaries["METHOD"].eq("GWASLab (bridge)")])
    mismatch_rows, mismatch_summary = add_ncbi_validation(mismatch_rows, mismatch_summary) # NCBI validation
    mismatch_summary.to_csv(mismatches_path, sep="\t", index=False)
    mismatch_rows.to_csv(mismatches_rows_path, sep="\t", index=False)

    plot_unmapped(summaries, plot_unmapped_path)
    plot_mismatches(mismatch_summary, plot_mismatches_path)

    print(f"Saved {mismatches_path}\nSaved {mismatches_rows_path}\nSaved {plot_unmapped_path}\nSaved {plot_mismatches_path}")



if __name__ == "__main__":
    main()
