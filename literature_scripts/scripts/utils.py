import pandas as pd
import numpy as np
import subprocess
import logging

from scipy.stats import norm



# ---- SANITY CHECK FUNCTION ----
def sanity_check(df, cohort, believe_metadata, uniprot_check_df):


    # ---- ALLELE CHECK ----
    n_alleles_orig = len(df)
    mask_bad_alleles = (
        df["EFFECT_ALLELE"].astype(str).str.contains("!", regex=False) |
        df["OTHER_ALLELE"].astype(str).str.contains("!", regex=False)
    )
    n_removed = mask_bad_alleles.sum()
    df = df.loc[~mask_bad_alleles].copy()

    if n_removed > 0:
        logging.warning(
            f"{cohort}: Removed {n_removed} rows (out of {n_alleles_orig}) "
            "with '!' in EA or NEA"
        )

    n_dropped = (df["EFFECT_ALLELE"] == "S").sum()
    df = df.loc[df["EFFECT_ALLELE"] != "S"].copy()
    df.reset_index(drop=True, inplace=True)
    if n_dropped > 0:
        logging.info(f"Dropped {n_dropped} rows with EFFECT_ALLELE == 'S'")


    # ---- SEQID CHECK ----
    if not df["SeqID"].dropna().empty:

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
            logging.info(f"Updated {n_fixed} malformed SEQIDs and pqtlIDs")


        # ---- UNIPROT CHECK ----
        uniprot_check(df, believe_metadata, uniprot_check_df)

    return df


# ---- SANITY CHECK HELPER FUNCTIONS ----

# Check UniProt against BELIEVE annotated metadata
def uniprot_check(df, believe_metadata, uniprot_check_df):
    cohort = df["COHORT"].unique()
    
    df["UniProt_orig"] = df["UniProt"].copy()
    df["UniProt"] = df["UniProt"].str.strip().str.upper()
    merged = df.merge(believe_metadata[["SeqID", "trait_protein_ids"]], on="SeqID", how="left")
    uniprot_mismatch_mask = (merged["trait_protein_ids"].notna() & (merged["UniProt"] != merged["trait_protein_ids"]))
    n_fixed_uniprot = uniprot_mismatch_mask.sum()
    mask = uniprot_mismatch_mask.to_numpy()
    df.loc[mask, "UniProt"] = (merged.loc[mask, "trait_protein_ids"].to_numpy())
    
    if n_fixed_uniprot > 0:
        uniprot_check_df.append(df.loc[df["UniProt"] != df["UniProt_orig"]].copy())
        logging.info(f"{cohort}: Fixed {n_fixed_uniprot} UniProt mismatches using BELIEVE metadata")
    df.drop(columns=["UniProt_orig"], inplace=True)


# Set column format and data types
def format_and_dtype(df, dtype_map, numeric_cols):

    expected_cols = list(dtype_map.keys())

    # Add missing columns
    missing_cols = [c for c in expected_cols if c not in df.columns]
    df = df.assign(**{c: pd.NA for c in missing_cols})

    # Reorder columns: expected first, extras last
    ordered_cols = expected_cols + [c for c in df.columns if c not in expected_cols]
    df = df[ordered_cols]

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
           f"Check https://github.com/ht-diva/believe_dm/commits/main/")

    # Save it to a file
    with open(file_name, 'w') as f:
        f.write(msg)

