import pandas as pd
import logging

from paths import PathManager
from utils import sanity_check, format_and_dtype, make_panels_mapping


# ---- PATHS ----
pm = PathManager()
LITERATURE_INPUT = pm.get_inputs()["literature_table_liftover"]
LITERATURE_INPUT_DIR = LITERATURE_INPUT.parent
BELIEVE_METADATA = pm.get_config()["believe_metadata"]
LITERATURE_PANEL = pm.get_config()["literature_panel"]
PANELS_MAP = pm.get_config()["panels_map"]
OUTPUT = pm.get_inputs()["literature_table_cleaned"]


# ---- FORMATS ----
NUMERIC_COLS = ["BETA", "SE", "minuslog10pval", "chr", "pos37", "pos38"]
DTYPE_MAP = {
    "pqtlID": "object",
    "rsID": "object",
    "chr": "int64",
    "pos37": "int64",
    "pos38": "int64",
    "SeqID": "object",
    "OlinkID": "object",
    "UniProt": "object",
    "OTHER_ALLELE": "category",
    "EFFECT_ALLELE": "category",
    "cis_trans": "category",
    "PMID": "int64",
    "BETA": "float64",
    "SE": "float64",
    "minuslog10pval": "float64",
    "SAMPLE_SIZE": "int64",
    "COHORT": "object",
    "TECHNOLOGY": "object",
    "Unit": "object",
}


# ---- LOGGING ----
log_file = LITERATURE_INPUT_DIR / "literature_table_all_somalogic_cleaned.log"
logging.basicConfig(
    filename=log_file,
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)


# ---- MAP BELIEVE & LITERATURE PANELS ----
if PANELS_MAP.exists():
    panels_map = pd.read_csv(PANELS_MAP, sep="\t")
else:
    panels_map = make_panels_mapping(BELIEVE_METADATA, LITERATURE_PANEL, PANELS_MAP)


# ---- READ LITERATURE TABLES ----
skip_sheets = {"credits", "variant", "protein", "olink", "cohort", "study"}
xls = pd.ExcelFile(LITERATURE_INPUT)
uniprot_check_df = []
uniprot_check_out = LITERATURE_INPUT_DIR / "uniprots_change.tsv"
missing_seqid_df = []
missing_seqid_out = LITERATURE_INPUT_DIR / "seqids_missing.tsv"
missing_uniprot_df = []
missing_uniprot_out = LITERATURE_INPUT_DIR / "uniprots_missing.tsv"

with pd.ExcelWriter(OUTPUT) as writer:
    for sheet in xls.sheet_names:
        new_sheet = sheet
        df = pd.read_excel(xls, sheet_name=sheet)
        if sheet.lower() in skip_sheets:
            print(f"Writing sheet: {sheet}")
            df.to_excel(writer, sheet_name=sheet, index=False)
            continue
        cohort = sheet
        logging.info(f"=== Processing {[str(cohort)]} ===")
        logging.info(f"Extracting: {cohort}")


        # ---- FORMAT ----
        df = format_and_dtype(df, DTYPE_MAP, NUMERIC_COLS)


        # ---- SANITY CHECK ----
        # Allele check: "S" or "!" are excluded
        # SeqID format check
        # UniProt check against BELIEVE and Literature Protein Panels
        # Find missing SeqID and UniProt against BELIEVE and Literature Protein Panels
        df = sanity_check(df, cohort, panels_map, uniprot_check_df, missing_seqid_df, missing_uniprot_df)


        # ---- SAVE ----
        print(f"Writing sheet: {sheet}")
        df.to_excel(writer, sheet_name=sheet, index=False)



# ---- SAVE OUTPUTS ----
logging.info(f"=== Saving outputs ===")
logging.info(f"> Written cleaned literature table to: {OUTPUT}")


# ---- SAVE UNIPROT CHECK ----
if uniprot_check_df:
    uniprot_check_df = pd.concat(uniprot_check_df, ignore_index=True)
    uniprot_check_df.to_csv(uniprot_check_out, sep="\t", index=False)
    logging.info(f"> Written UniProt correction variant table to: {uniprot_check_out}")


# ---- SAVE MISSING SEQIDs ----
if missing_seqid_df:
    print("\n=== MISSING SEQIDs ===")

    missing_seqid_df = pd.concat(missing_seqid_df, ignore_index=True)
    missing_seqid_df = (
        missing_seqid_df[["COHORT", "SEQID_MISSING", "UNIPROT", "VARIANTS_NR"]]
        .sort_values(by=["COHORT", "VARIANTS_NR"], ascending=[True, False])
    )
    missing_seqid_df.to_csv(missing_seqid_out, sep="\t", index=False)
    logging.info(f"> Written Missing SEQIDs table to: {missing_seqid_out}")

    summary_df = missing_seqid_df.groupby("COHORT").agg(
        SEQID_MISSING=("SEQID_MISSING", 'nunique'),
        VARIANTS_NR=("VARIANTS_NR", 'sum')
    ).reset_index()
    missing_seqid_summary_out = str(missing_seqid_out).replace(".tsv","_summary.tsv")
    summary_df.to_csv(missing_seqid_summary_out, sep="\t", index=False)
    logging.info(f"> Written Missing SEQIDs summary to: {missing_seqid_summary_out}")

    print(summary_df)


# ---- SAVE MISSING UNIPROTs ----
if missing_uniprot_df:
    print("\n=== MISSING UNIPROTs ===")

    missing_uniprot_df = pd.concat(missing_uniprot_df, ignore_index=True)
    missing_uniprot_df = (
        missing_uniprot_df[["COHORT", "UNIPROT_MISSING", "SEQID", "VARIANTS_NR"]]
        .sort_values(by=["COHORT", "VARIANTS_NR"], ascending=[True, False])
    )
    missing_uniprot_df.to_csv(missing_uniprot_out, sep="\t", index=False)
    logging.info(f"> Written Missing UniProts table to: {missing_uniprot_out}")

    summary_df = missing_uniprot_df.groupby("COHORT").agg(
        UNIPROT_MISSING=("UNIPROT_MISSING", 'nunique'),
        VARIANTS_NR=("VARIANTS_NR", 'sum')
    ).reset_index()
    missing_uniprot_summary_out = str(missing_uniprot_out).replace(".tsv","_summary.tsv")
    summary_df.to_csv(missing_uniprot_summary_out, sep="\t", index=False)
    logging.info(f"> Written Missing UniProts summary to: {missing_uniprot_summary_out}")

    print(summary_df)