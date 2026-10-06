import pandas as pd
import logging
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Patch


# ---- MISSING SEQIDs & UNIPROTs FUNCTIONS ----


# Check missing SeqIDs and UniProts against BELIEVE and Literature Protein Panels
def missing_seqid_uniprot(df, cohort, panels_map, missing_ids_df):

    # Check for both raw and harmonized datasets
    column_aliases = {}
    if "SeqID" not in df.columns and "SEQID" in df.columns:
        column_aliases["SEQID"] = "SeqID"
    if "UniProt" not in df.columns and "UNIPROT" in df.columns:
        column_aliases["UNIPROT"] = "UniProt"
    df = df.rename(columns=column_aliases)


    # SomaScan Panel
    if not df["SeqID"].dropna().empty:


        # ---- CHECK MISSING SEQID ----

        # Get nr. SeqIDs/variants unmatching BELIEVE reference
        seqid_df = df.loc[
            ~df["SeqID"].isin(panels_map["SeqID"]) & df["SeqID"].notna()
        ].copy()
        missing_seqids_var_nr = len(seqid_df)
        missing_seqids = seqid_df["SeqID"].unique()
        missing_seqids_nr = len(missing_seqids)

         # If any missing SeqIDs...
        if missing_seqids_nr > 0:
            logging.info(
                f"> {cohort}: Missing SeqIDs: {missing_seqids_nr} out of {len(set(df.SeqID))}."
                f" Variants affected: {missing_seqids_var_nr} out of {len(df)}."
            )

            # ...Get nr. missing variants per SeqID
            seqid_summary = (
                seqid_df[["COHORT", "SeqID", "UniProt"]]
                .groupby(["COHORT", "SeqID"])
                .agg(
                    VARIANT_NR=("SeqID", "size"),
                    UNIPROT=("UniProt", lambda x: ", ".join(x.dropna().astype(str).unique()))
                )
                .reset_index()
            )

            seqid_summary["COHORT"] = cohort
            seqid_summary["PANEL"] = "SomaScan"
            seqid_summary["ID_MISSING"] = "SeqID"
            seqid_summary = seqid_summary[
                [
                    "COHORT",
                    "PANEL",
                    "ID_MISSING",
                    "SeqID",
                    "UNIPROT",
                    "VARIANT_NR",
                ]
            ].rename(columns={"SeqID": "SEQID"})
            seqid_summary["VARIANT_TOT_NR"] = len(df)

            missing_ids_df.append(seqid_summary)

        else:
            logging.info(f"> {cohort}: No Missing SeqID.")


    # Olink Panel
    else:


        # ---- CHECK MISSING UNIPROT ----

        # Get nr. proteins/variants unmatching BELIEVE reference
        uniprots_df = df.loc[
            ~df["UniProt"].isin(panels_map["UniProt_BELIEVE"]) & df["UniProt"].notna()
        ].copy()
        missing_uniprots_var_nr = len(uniprots_df)
        missing_uniprots = uniprots_df["UniProt"].unique()
        missing_uniprots_nr = len(missing_uniprots)

        # If any missing UniProts...
        if missing_uniprots_nr > 0:
            logging.info(
                f"> {cohort}: Missing UniProts: {missing_uniprots_nr} out of {len(set(df.UniProt))}."
                f" Variants affected: {missing_uniprots_var_nr} out of {len(df)}."
            )

            # ...Get nr. missing variants per UniProt
            uniprot_summary = (
                uniprots_df[["COHORT", "SeqID", "UniProt"]]
                .groupby(["COHORT", "UniProt"])
                .agg(
                    VARIANT_NR=('UniProt', 'size'),
                    SEQID=('SeqID', lambda x: ', '.join(x.dropna().astype(str).unique()))
                )
                .reset_index()
            )

            # Format and Store results
            uniprot_summary["COHORT"] = cohort
            uniprot_summary["PANEL"] = "Olink"
            uniprot_summary["ID_MISSING"] = "UniProt"
            uniprot_summary = uniprot_summary[
                [
                    "COHORT",
                    "PANEL",
                    "ID_MISSING",
                    "SEQID",
                    "UniProt",
                    "VARIANT_NR",
                ]
            ].rename(columns={"UniProt": "UNIPROT"})
            uniprot_summary["VARIANT_TOT_NR"] = len(df)

            missing_ids_df.append(uniprot_summary)

        else:
            logging.info(f"> {cohort}: No Missing UniProt.")


# Report missing SeqIDs and UniProts
def report_missing_ids(missing_ids_df):

    panel_id_type = {
        "SomaScan": "SeqID",
        "Olink": "UniProt",
    }

    missing_ids_df["ID_TYPE"] = missing_ids_df["PANEL"].map(panel_id_type)
    missing_ids_df["MISSING_ID_VALUE"] = missing_ids_df["SEQID"].where(
        missing_ids_df["ID_TYPE"].eq("SeqID"),
        missing_ids_df["UNIPROT"],
    )

    missing_summary_df = (
        missing_ids_df.loc[missing_ids_df["ID_MISSING"].eq(missing_ids_df["ID_TYPE"])]
        .groupby(["COHORT", "PANEL", "ID_TYPE"], as_index=False)
        .agg(
            ID_MISSING=("MISSING_ID_VALUE", "nunique"),
            ID_MISSING_VARS=("VARIANT_NR", "sum"),
            TOTAL_VARS=("VARIANT_TOT_NR", "first"),
        )
        .sort_values(["COHORT", "PANEL"])
        .reset_index(drop=True)
    )

    count_columns = ["ID_MISSING", "ID_MISSING_VARS", "TOTAL_VARS"]
    missing_summary_df[count_columns] = missing_summary_df[count_columns].fillna(0).astype(int)

    return missing_summary_df


# Plot  missing SeqIDs and UniProts
def plot_missing_ids(missing_summary_df, output_path):

    # Calculate percentage of variants affected by missing identifiers
    missing_summary_df["MISSING_ID_VARS_PCT"] = np.where(
        missing_summary_df["TOTAL_VARS"].gt(0),
        (missing_summary_df["ID_MISSING_VARS"] / missing_summary_df["TOTAL_VARS"]) * 100,
        np.nan,
    )

    # Order cohorts by percentage
    plot_df = missing_summary_df.sort_values("MISSING_ID_VARS_PCT", ascending=False).reset_index(drop=True)

    # Plot settings
    acadia = [
        "#A4BED5",
        "#FED789",
        "#72874E",
        "#023743",
        "#476F84",
        "#453947",
    ]
    panel_colors = dict(zip(["SomaScan", "Olink"], acadia[:2]))
    bar_colors = plot_df["PANEL"].map(panel_colors).fillna("#A0A0A0")
    fig, ax = plt.subplots(figsize=(14, 6))

    # Bars
    bars = ax.bar(
        range(len(plot_df)),
        plot_df["MISSING_ID_VARS_PCT"],
        color=bar_colors,
        edgecolor="black",
        linewidth=0.5,
    )

    # Absolute total missingness above bars
    ax.bar_label(
        bars,
        labels=[f"{value:,}" for value in plot_df["ID_MISSING_VARS"]],
        padding=3,
        fontsize=9,
    )
    ax.set_xticks(range(len(plot_df)))

    # Labels
    ax.set_xticklabels(plot_df["COHORT"], rotation=45, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("Variants with missing IDs (%)")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    maximum = plot_df["MISSING_ID_VARS_PCT"].max()
    if pd.notna(maximum): 
        ax.set_ylim(0, max(maximum * 1.15, 1))

    # Legend
    legend_handles = [Patch(facecolor=color, edgecolor="black", label=panel)
        for panel, color in panel_colors.items()
        if panel in plot_df["PANEL"].unique()
    ]
    ax.legend(handles=legend_handles, title="Panel")

    # Save plot
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
