import logging

from .missing_ids import missing_seqid_uniprot
from .uniprot_map import uniprot_check


# ---- MAIN SANITY CHECK FUNCTION ----


# Allele check: "S" or "!" are excluded
# SeqID format check
# UniProt check against BELIEVE and Literature Protein Panels
# Find missing SeqID and UniProt against BELIEVE and Literature Protein Panels
def sanity_check(df, cohort, panels_map, harmonized_uniprots_df, missing_ids_df):

    logging.info(f"Sanity check for {cohort}")


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
    df = uniprot_check(df, cohort, panels_map, harmonized_uniprots_df)


    # ---- MISSING SEQID AND UNIPROT ----
    missing_seqid_uniprot(df, cohort, panels_map, missing_ids_df)

    return df

