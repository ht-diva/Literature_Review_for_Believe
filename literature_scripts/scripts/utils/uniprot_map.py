import pandas as pd
import logging

from .helper import best_id_match, is_uniprot_match, swap_uniprots


# ---- MAP PROTEIN PANELS FUNCTION ----


# Function to map BELIEVE and Literature Protein Panels
def make_panels_mapping(believe_metadata_path, literature_panel_path, panels_map_path):

    logging.info("=== Map BELIEVE and Literature Protein Panels ===")

    # Read BELIEVE Panel derived from SomaScan Annotated Panel 7k
    believe_metadata = pd.read_csv(believe_metadata_path, sep="\t", usecols=["trait_desc", "trait_seqid", "trait_gene_ids", "trait_protein_ids"])
    believe_metadata = believe_metadata.rename(columns={
        "trait_desc": "Target_Name",
        "trait_seqid": "SeqID",
        "trait_gene_ids": "Ensembl_Gene_ID",
        "trait_protein_ids": "UniProt",
        })
    believe_metadata["UniProt"] = (believe_metadata["UniProt"].str.strip().str.upper())

    # Read Literature Protein Panel
    literature_panel = pd.read_csv(literature_panel_path, sep="\t", usecols=["Target_Name", "SeqID", "Ensembl_Gene_ID", "UniProt"])
    literature_panel.loc[:, "SeqID"] = literature_panel.loc[:, "SeqID"].replace("_", "-", regex=True)

    # Get best ID to match BELIEVE and Literature Protein Panels 
    best_match = best_id_match(believe_metadata, literature_panel)

    merged = pd.merge(
        believe_metadata[["Target_Name", best_match, "UniProt"]],
        literature_panel[[best_match, "UniProt"]],
        on=best_match,
        how="inner"
    ).drop_duplicates().reset_index(drop=True)

    # Align swapped multi-Prots
    # Example: P29460|Q9NPF7 <-> Q9NPF7|P29460
    mask = swap_uniprots(merged, "UniProt_x", "UniProt_y")
    merged.loc[mask, "UniProt_y"] = merged.loc[mask, "UniProt_x"]
    merged = merged.drop_duplicates().reset_index(drop=True)

    # Check multiple UniProts per SeqID
    multi_uniprots_per_seqid = len(merged[merged.duplicated(subset=best_match, keep="first")])

    # Group by SeqID and aggregate multiple UniProt (literature) values
    if multi_uniprots_per_seqid > 0:
        logging.info(f"> SeqIDs with multiple UniProts: {multi_uniprots_per_seqid}.")
        merged = (
            merged.groupby(best_match)
            .agg({
                "Target_Name": "first",
                "UniProt_x": "first",
                "UniProt_y": lambda x: ", ".join(sorted(set(filter(None, x))))
            })
            .reset_index()
        )

    # Check matching UniProts
    merged["UniProt_Match"] = merged.apply(
        lambda row: is_uniprot_match(row, "UniProt_x", "UniProt_y"),
        axis=1
    )
    uniprot_match_df = merged[merged["UniProt_Match"]].reset_index(drop=True)
    uniprot_mismatch_df = merged[~merged["UniProt_Match"]].reset_index(drop=True)
    merged = pd.concat([uniprot_match_df, uniprot_mismatch_df]).reset_index(drop=True)
    uniprot_mismatch_var_nr = len(uniprot_mismatch_df)
    uniprot_mismatch_nr = len(set(uniprot_mismatch_df["UniProt_y"]))
    if uniprot_mismatch_nr > 0:
        logging.info(f"> Mismatched UniProts: {uniprot_mismatch_nr}. Variants affected: {uniprot_mismatch_var_nr}.")
    else:
        logging.info("> No Mismatched UniProt.")

    # Format and Save
    panels_map = merged.rename(columns={
        "UniProt_x": "UniProt_BELIEVE",
        "UniProt_y": "UniProt_Literature"
    })
    panels_map["SeqID"] = "seq." + panels_map["SeqID"].str.replace("-", ".", regex=False)

    panels_map.to_csv(panels_map_path, sep="\t", index=False)
    logging.info(f"> Written mapping file of BELIEVE and Literature Protein Panels to: {panels_map_path}")

    return panels_map


# Check (& Harmonize) UniProt against BELIEVE and Literature Protein Panels
def uniprot_check(df, cohort, panels_map, harmonized_uniprots_df):


    # ---- CLEAN UNIPROT ----
    df["UniProt"] = (
        df["UniProt"]
        .fillna("")  # Fill NaN values with empty string
        .astype(str)
        .str.strip()  # Remove leading and trailing spaces
        .str.replace(r"[,;|\s]+", "|", regex=True)  # Replace commas, semicolons, pipes and whitespace with |
        .str.replace(r"\s*\|\s*", "|", regex=True)  # Remove spaces around |
        .str.replace(r"\s+", "|", regex=True)  # Replace any remaining spaces with |
        # Remove duplicate UniProt identifiers (e.g. P47929;P47929 -> P47929)
        .apply(
            lambda value: "|".join(
                dict.fromkeys(value.split("|"))
            ) if value else ""
        )
    )


    # ---- CHECK UNIPROT MATCH cf. BELIEVE ----

    # Only SomaScan
    if not df["SeqID"].dropna().empty:

        # Merge with BELIEVE Metadata by SEQID
        merged = df.merge(
            panels_map[["SeqID", "UniProt_BELIEVE", "UniProt_Literature"]], 
            on="SeqID", 
            how="left"
        )
        
        # Align swapped multi-Prots
        mask = swap_uniprots(merged, "UniProt", "UniProt_BELIEVE")
        merged.loc[mask, "UniProt"] = merged.loc[mask, "UniProt_BELIEVE"]
        merged = merged.drop_duplicates().reset_index(drop=True)

        # Update aligned swapped multi-Prots
        df.loc[mask, "UniProt"] = merged.loc[mask, "UniProt"]

        # Fill empty UniProts
        mask = merged["UniProt"] == ""
        merged.loc[mask, "UniProt"] = merged.loc[mask, "UniProt_BELIEVE"]
        df.loc[mask, "UniProt"] = merged.loc[mask, "UniProt_BELIEVE"]

        # Mismatched UniProts (exclude NaN)
        merged["UniProt_Match"] = (
            (merged["UniProt"] == merged["UniProt_BELIEVE"]) |
            (merged["UniProt"] == "") |
            (merged["UniProt_BELIEVE"] == "")
        )
        uniprot_mismatch_df = merged[~merged["UniProt_Match"]].reset_index(drop=True)
        uniprot_mismatch_var_nr = len(uniprot_mismatch_df)
        uniprot_mismatch_nr = len(set(uniprot_mismatch_df["UniProt"]))

        # If any mismatched UniProts...
        if uniprot_mismatch_nr > 0:
            logging.info(f"> Mismatched UniProts: {uniprot_mismatch_nr}. Variants affected: {uniprot_mismatch_var_nr}.")

            # ...Store mismatched UniProts
            uniprot_mismatch_df = merged.loc[~merged["UniProt_Match"]].reset_index(drop=True)
            uniprot_mismatch_df.loc[:, "UniProt_Raw"] = uniprot_mismatch_df.loc[:, "UniProt"]

            # ...Store mismatched UniProts to BELIEVE UniProts
            uniprot_mismatch_df.loc[:, "UniProt"] = uniprot_mismatch_df.loc[:, "UniProt_BELIEVE"]

            # ...Update mismatched UniProts to BELIEVE UniProts
            df = df.merge(
                uniprot_mismatch_df[["SeqID", "UniProt"]],
                on="SeqID",
                how="left",
                suffixes=("", "_updated")
            )
            df.loc[df["UniProt_updated"].notna(), "UniProt"] = df.loc[df["UniProt_updated"].notna(), "UniProt_updated"]
            logging.info(f"  |-> {uniprot_mismatch_nr} Mismatched UniProts ({len(df.UniProt_updated.notna())} variants) updated to BELIEVE UniProts.")
            df = df.drop(columns=["UniProt_updated"])

            # Store results
            uniprot_mismatch_df = uniprot_mismatch_df[
                ["COHORT", "SeqID", "UniProt_Raw", "UniProt"]
            ].drop_duplicates().reset_index(drop=True)

            uniprot_mismatch_df["COHORT"] = cohort
            uniprot_mismatch_df = uniprot_mismatch_df.rename(columns={
                "SeqID" : "SEQID",
                "UniProt_Raw": "UNIPROT_RAW",
                "UniProt" : "UNIPROT"
            })
            harmonized_uniprots_df.append(uniprot_mismatch_df)

        else:
            logging.info("> No Mismatched UniProts.")

    return df
