import pandas as pd
import logging

from paths import PathManager
from utils.sanity_check import sanity_check
from utils.helper import format_and_dtype
from utils.uniprot_map import make_panels_mapping
from utils.missing_ids import plot_missing_ids, report_missing_ids


# ---- PATHS ----
pm = PathManager()
LITERATURE_INPUT = pm.get_inputs()["literature_table"]
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
harmonized_uniprots_df = []
harmonized_uniprots_out = LITERATURE_INPUT_DIR / "harmonized_uniprots.tsv"
harmonized_uniprots_summary_out = LITERATURE_INPUT_DIR / "harmonized_uniprots_summary.tsv"
missing_ids_df = []
missing_ids_out = LITERATURE_INPUT_DIR / "missing_ids.tsv"
missing_ids_summary_out = LITERATURE_INPUT_DIR / "missing_ids_summary.tsv"
missing_ids_plot_out = LITERATURE_INPUT_DIR / "missing_ids.png"

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
        df = sanity_check(df, cohort, panels_map, harmonized_uniprots_df, missing_ids_df)


        # ---- SAVE ----
        print(f"Writing sheet: {sheet}")
        df = df.drop_duplicates().reset_index(drop=True)
        df.to_excel(writer, sheet_name=sheet, index=False)



# ---- REPORT HARMONIZED UNIPROTs ----
if harmonized_uniprots_df:

    # Save harmonized Uniprots
    harmonized_uniprots_df = pd.concat(harmonized_uniprots_df, ignore_index=True)
    harmonized_uniprots_df.to_csv(harmonized_uniprots_out, sep="\t", index=False)
    logging.info(f"> Written UniProt harmonization table to: {harmonized_uniprots_out}")

    # Summary of harmonized Uniprots
    harmonized_uniprots_summary = (harmonized_uniprots_df[["SEQID", "UNIPROT_RAW", "UNIPROT"]]
        .drop_duplicates()
        .sort_values(["SEQID", "UNIPROT_RAW", "UNIPROT"])
        .reset_index(drop=True)
    )
    harmonized_uniprots_summary.to_csv(harmonized_uniprots_summary_out, sep="\t", index=False)


# ---- REPORT MISSING IDs ----
if missing_ids_df:

    # Save missing IDs
    missing_ids_df = pd.concat(missing_ids_df, ignore_index=True)
    missing_ids_df = missing_ids_df.sort_values(
        by=["COHORT", "PANEL", "VARIANT_NR"],
        ascending=[True, True, False],
    )
    missing_ids_df.to_csv(missing_ids_out, sep="\t", index=False)
    logging.info(f"> Written Missing SeqIDs & UniProts table to: {missing_ids_out}")


    # Reports for missing IDs
    missing_summary_df = report_missing_ids(missing_ids_df)
    missing_summary_df.to_csv(missing_ids_summary_out, sep="\t", index=False)
    logging.info(f"> Written missing identifier summary to: {missing_ids_summary_out}")

    print("\n=== MISSING IDENTIFIERS ===")
    print(missing_summary_df)

    # Plot missing IDs
    plot_missing_ids(missing_summary_df, missing_ids_plot_out)
