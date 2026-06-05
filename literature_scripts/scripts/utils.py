import pandas as pd
import numpy as np
import subprocess
import logging
import re

from scipy.stats import norm



# ---- HELPER FUNCTIONS ----

# Helper function to make unique identifier (variant key)
def make_variant_key(df, id_col, chr_col, pos_col, a1_col, a2_col):
    alleles = (
        df[[a1_col, a2_col]]
        .astype(str)
        .apply(lambda x: "_".join(sorted(x)), axis=1)
    )

    return (
        df[id_col].astype(str)
        + "_"
        + df[chr_col].astype(str)
        + "_"
        + df[pos_col].astype(str)
        + "_"
        + alleles
    )


# Helper function to align swapped multi-Prots
# Example: P29460|Q9NPF7 <-> Q9NPF7|P29460
def swap_uniprots(df, uniprot1, uniprot2):
    df[uniprot1] = df[uniprot1].fillna("")
    df[uniprot2] = df[uniprot2].fillna("")
    uniprot1_set = df[uniprot1].str.split("|").apply(lambda x: set(sorted(x)))
    uniprot2_set = df[uniprot2].str.split("|").apply(lambda x: set(sorted(x)))
    mask = (
        (uniprot1_set == uniprot2_set) &
        (df[uniprot1] != df[uniprot2])
    )
    return mask


# Helper function to check (multi-)UniProt match
# Example P0DMV8-P16519 matches P0DMV8 and P16519
def is_uniprot_match(row, uniprot1, uniprot2):
    uniprot1set = set(str(row[uniprot1]).split(", ") if pd.notna(row[uniprot1]) else [])
    uniprot2set = set(str(row[uniprot2]).split(", ") if pd.notna(row[uniprot2]) else [])
    return uniprot1set.issubset(uniprot2set) or uniprot2set.issubset(uniprot1set)


# Set column format and data types
def format_and_dtype(df, dtype_map, numeric_cols):

    expected_cols = list(dtype_map.keys())

    # Add missing columns
    missing_cols = [c for c in expected_cols if c not in df.columns]
    df = df.assign(**{c: pd.NA for c in missing_cols})

    # Reorder columns: expected first, extras last
    ordered_cols = expected_cols + [c for c in df.columns if c not in expected_cols]
    df = df[ordered_cols]

    # Format chromosome
    df["chr"] = (
        df["chr"]
        .astype(str)
        .replace({"X": "23", "Y": "24"})
    )

    # Numeric coercion
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(
                df[col].astype(str).str.replace(",", ".", regex=False),
                errors="coerce"
            )

    # Apply dtypes
    for col, dtype in dtype_map.items():
        if col in df.columns:
            try:
                df[col] = df[col].astype(dtype)
            except Exception as e:
                logging.warning(f"Could not cast {col} to {dtype}: {e}")

    return df


# Helper function to find the best ID to match BELIEVE and Literature Protein Panels
# IDs are:
#   "Target_Name"
#   "SeqID"
#   "Ensembl_Gene_ID"
#   "UniProt"
def best_id_match(believe_metadata, literature_panel):

    # Count matches
    matches = {
        "Target_Name": 0,
        "SeqID": 0,
        "Ensembl_Gene_ID": 0,
        "UniProt": 0
    }

    for col in matches.keys():
        merged = pd.merge(
            believe_metadata[[col]],
            literature_panel[[col]],
            on=col,
            how="inner"
        ).drop_duplicates()
        matches[col] = len(merged)

    for col, count in matches.items():
        logging.info(f"> Number of matches for {col}: {count} on {len(believe_metadata)}.")

    # Best match
    best_match = max(matches, key=matches.get)
    best_count = matches[best_match]
    logging.info(f"> ID with the maximum matches: {best_match} ({best_count} matches).")

    # Missing match
    all_matches = pd.merge(
        believe_metadata[[best_match]],
        literature_panel[[best_match]],
        on=best_match,
        how="inner"
    ).drop_duplicates()
    missing_matches = believe_metadata[~believe_metadata[best_match].isin(all_matches[best_match])]

    if len(missing_matches) > 0:
        logging.info(f"> Missing {best_match}: {missing_matches}.")
    else:
        logging.info(f"> All matched by {best_match}.")

    return best_match



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



# ---- MAIN SANITY CHECK FUNCTION ----


# Allele check: "S" or "!" are excluded
# SeqID format check
# UniProt check against BELIEVE and Literature Protein Panels
# Find missing SeqID and UniProt against BELIEVE and Literature Protein Panels
def sanity_check(df, cohort, panels_map, uniprot_check_df, missing_seqid_df, missing_uniprot_df):

    logging.info(f"Sanity check for {cohort}")


    # ---- ALLELE CHECK ----
    n_alleles_orig = len(df)
    mask_bad_alleles = (
        df["EFFECT_ALLELE"].astype(str).str.contains(
            r"!|\.|NAN|\*",
            case=False,
            regex=True
        ) |
        df["OTHER_ALLELE"].astype(str).str.contains(
            r"!|\.|NAN|\*",
            case=False,
            regex=True
        )
    )
    n_removed = mask_bad_alleles.sum()
    df = df.loc[~mask_bad_alleles].copy()

    if n_removed > 0:
        logging.warning(
            f"{cohort}: Removed {n_removed} rows (out of {n_alleles_orig}) "
            "with '!' or '.' in EA or NEA"
        )

    n_dropped = (df["EFFECT_ALLELE"] == "S").sum()
    df = df.loc[df["EFFECT_ALLELE"] != "S"].copy()
    df.reset_index(drop=True, inplace=True)
    if n_dropped > 0:
        logging.info(f"Dropped {n_dropped} rows with EFFECT_ALLELE == 'S'")


    # ---- CIS-TRANS FORMATTING ----
    df["cis_trans"] = df["cis_trans"].str.strip().str.lower()


    # ---- SEQID CHECK ----
    if not df["SeqID"].dropna().empty:

        # Eliminate null SeqIDs
        seqna_mask = df["SeqID"] == "seq.NA"
        n_seqna = seqna_mask.sum()
        if n_seqna > 0:
            df = df.loc[~seqna_mask].reset_index(drop=True)
            logging.info(f"> Eliminated {n_seqna} seq.NA entries")

        # Detect malformed SeqIDs
        malformed_mask = df["SeqID"].str.match(r"^seq\.\d+\.\d+\.\d+$", na=False)
        n_fixed = malformed_mask.sum()
        if n_fixed > 0:
            
            # Fix SeqID
            df.loc[malformed_mask, "SeqID"] = (df.loc[malformed_mask, "SeqID"].str.replace(r"^(seq\.\d+\.\d+)\.\d+$", r"\1", regex=True))

            # Update pqtlID
            df.loc[malformed_mask, "pqtlID"] = (
                df.loc[malformed_mask, "rsID"]
                + "_" + df.loc[malformed_mask, "SeqID"]
                + "_" + df.loc[malformed_mask, "PMID"].astype(str)
                + "_" + df.loc[malformed_mask, "COHORT"]
            )
            logging.info(f"> Updated {n_fixed} malformed SEQIDs and pqtlIDs")


    # ---- UNIPROT CHECK / FORMATTING ----
    df = uniprot_check(df, cohort, panels_map, uniprot_check_df)


    # ---- MISSING SEQID AND UNIPROT ----
    missing_seqid_uniprot(df, cohort, panels_map, missing_seqid_df, missing_uniprot_df)

    return df



# ---- UNIPROT FORMAT & MATCH FUNCTION ----


# Check UniProt against BELIEVE and Literature Protein Panels
def uniprot_check(df, cohort, panels_map, uniprot_check_df):


    # ---- CLEAN UNIPROT ----
    df["UniProt"] = (
        df["UniProt"]
        .fillna("")  # Fill NaN values with empty string
        .str.strip()  # Remove leading and trailing spaces
        .str.replace(",", "|", regex=False)  # Replace commas with |
        .str.replace(r'\s*\|\s*', '|', regex=True)  # Remove spaces around |
        .str.replace(r'\s+', '|', regex=True)  # Replace any remaining spaces with |
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
            uniprot_check_df.append(uniprot_mismatch_df)

        else:
            logging.info("> No Mismatched UniProts.")

    return df



# ---- MISSING SEQIDs & UNIPROTs FUNCTION ----


# Check  missing SeqIDs and UniProts against BELIEVE and Literature Protein Panels
def missing_seqid_uniprot(df, cohort, panels_map, missing_seqid_df, missing_uniprot_df):


    # ---- CHECK MISSING SEQID ----

    # Only SomaScan
    if not df["SeqID"].dropna().empty:

        # Get nr. SeqIDs/variants unmatching BELIEVE reference
        seqid_df = df.loc[
            ~df["SeqID"].isin(panels_map["SeqID"]) & df["SeqID"].notna()
        ]
        missing_seqids_var_nr = len(seqid_df)
        missing_seqids = seqid_df["SeqID"].unique()
        missing_seqids_nr = len(missing_seqids)

         # If any missing SeqIDs...
        if missing_seqids_nr > 0:
            logging.info(
                f"> Missing SeqIDs: {missing_seqids_var_nr} out of {len(set(df.SeqID))}."
                f" Variants affected: {missing_seqids_var_nr} out of {len(df)}."
            )

            # ...Get nr. missing variants per SeqId
            seqid_df = (
                seqid_df[["COHORT", "SeqID", "UniProt"]]
                .groupby(["COHORT", "SeqID"])
                .agg(
                    VARIANTS_NR=('SeqID', 'size'),
                    UNIPROT=('UniProt', lambda x: ', '.join(x.unique()))
                )
                .reset_index()
            )

            # Format and Store results
            seqid_df["COHORT"] = cohort
            seqid_df = seqid_df.rename(columns={
                "SeqID" : "SEQID_MISSING"
            })
            missing_seqid_df.append(seqid_df)

        else:
            logging.info("> No Missing SeqID.")


    # ---- CHECK MISSING UNIPROT ----

    # Get nr. proteins/variants unmatching BELIEVE reference
    uniprots_df = df.loc[
        ~df["UniProt"].isin(panels_map["UniProt_BELIEVE"]) & df["UniProt"].notna()
    ]
    missing_uniprots_var_nr = len(uniprots_df)
    missing_uniprots = uniprots_df["UniProt"].unique()
    missing_uniprots_nr = len(missing_uniprots)

    # If any missing UniProts...
    if missing_uniprots_nr > 0:
        logging.info(
            f"> Missing UniProts: {missing_uniprots_nr} out of {len(set(df.UniProt))}."
            f" Variants affected: {missing_uniprots_var_nr} out of {len(df)}."
        )

        # ...Get nr. missing variants per UniProt
        uniprots_df = (
            uniprots_df[["COHORT", "SeqID", "UniProt"]]
            .groupby(["COHORT", "UniProt"])
            .agg(
                VARIANTS_NR=('UniProt', 'size'),
                SEQID=('SeqID', lambda x: ', '.join(x.dropna().astype(str).unique()))
            )
            .reset_index()
        )

        # Format and Store results
        uniprots_df["COHORT"] = cohort
        uniprots_df = uniprots_df.rename(columns={
            "UniProt" : "UNIPROT_MISSING"
        })
        missing_uniprot_df.append(uniprots_df)

    else:
        logging.info("> No Missing UniProt.")



# ---- GIT COMMIT FUNCTIONS ----


def get_last_commit_id():
    # Run git log -1 --format=%H and capture its output
    last_commit_id = subprocess.check_output(
        ['git', 'log', '-1', '--format=%H']).decode('utf-8').strip()

    return last_commit_id

def save_last_commit_id_to_file(file_name):
    # Get the last commit ID
    last_commit_id = get_last_commit_id()
    msg = (f"This folder contains data produced by this commit id {last_commit_id} of the code.\n"
           f"Check out: https://github.com/ht-diva/Literature_Review_for_Believe/main")

    # Save it to a file
    with open(file_name, 'w') as f:
        f.write(msg)



# ---- CONVERT TO VCF FUNCTION ----

def write_vcf(df, output_filename, build="GRCh37"):
    with open(output_filename, "w") as vcf_file:

        # Ensure POS is valid integer
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
        df = df.dropna(subset=["POS"])
        df["POS"] = df["POS"].astype(int)

        # Convert CHR for bcftools processing
        df["CHR"] = df["CHR"].astype(str).replace({
            "23": "X",
            "24": "Y",
            "25": "MT"
        })
        if build == "GRCh38":
            df["CHR"] = df["CHR"].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)

        # Write the VCF header
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write("##source=PythonScript\n")
        vcf_file.write(f"##reference={build}\n")
        vcf_file.write('##INFO=<ID=BETA,Number=1,Type=Float,Description=Effect Size Estimate>\n')
        vcf_file.write('##INFO=<ID=SE,Number=1,Type=Float,Description=Standard Error>\n')
        vcf_file.write('##INFO=<ID=N,Number=1,Type=Integer,Description=Sample Size>\n')
        vcf_file.write('##INFO=<ID=MLOG10P,Number=1,Type=Float,Description=Negative Log10 P-value>\n')
        chroms = df["CHR"].unique()
        for c in chroms:
            vcf_file.write(f"##contig=<ID={c}>\n")
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Iterate through rows
        for _, row in df.iterrows():
            chrom = row["CHR"]
            pos = row["POS"]
            vid = row["rsID"]
            ref = row["EA"]
            alt = row["NEA"]
            qual = "."
            filt = "."

            info = (
                f"BETA={row['BETA']};"
                f"SE={row['SE']};"
                f"N={row['N']};"
                f"MLOG10P={row['MLOG10P']}"
            )

            line = f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\n"
            vcf_file.write(line)
