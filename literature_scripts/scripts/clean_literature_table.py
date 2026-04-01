import pandas as pd
import logging

from paths import PathManager
from utils import sanity_check, format_and_dtype, make_panels_mapping


# ---- PATHS ----
pm = PathManager()
LITERATURE_INPUT = pm.get_inputs()["literature_table_raw"]
LITERATURE_INPUT_DIR = LITERATURE_INPUT.parent
BELIEVE_METADATA = pm.get_config()["believe_metadata"]
LITERATURE_PANEL = pm.get_config()["literature_panel"]
PANELS_MAP = pm.get_config()["panels_map"]
LITERATURE_FILES = pm.get_files()
SUN_UKB_NONEU = LITERATURE_FILES["pqtl_sun_ukb_csa"]
INTERVAL_CHRIS_META = LITERATURE_FILES["pqtl_interval_chris_meta"]
DECODE_2023 = LITERATURE_FILES["pqtl_decode_2023"]
OUTPUT = pm.get_inputs()["literature_table"]


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


        # ---- SANITY CHECK ----
        # Allele check: "S" or "!" are excluded
        # SeqID format check
        # UniProt check against BELIEVE and Literature Protein Panels
        # Find missing SeqID and UniProt against BELIEVE and Literature Protein Panels
        df = sanity_check(df, cohort, panels_map, uniprot_check_df, missing_seqid_df, missing_uniprot_df)


        # ---- SAVE ----
        print(f"Writing sheet: {sheet}")
        df.to_excel(writer, sheet_name=sheet, index=False)


        # ---- ADD MISSING STUDIES ----

        # deCODE 2023 to literature review format
        if sheet == "pqtl_decode":
            new_sheet = "pqtl_decode_2023"
            logging.info(f"=== Processing {[str(new_sheet)]} ===")
            logging.info(f"Extracting: {new_sheet}")

            df = pd.read_csv(DECODE_2023, sep="\t", 
                             usecols=["rsID","SeqID","UniProt","cis_trans",
                                      "Chrom","Pos","effectAllele","otherAllele",
                                      "Beta","SE","minus_log10_pval"])
            df = df.rename(columns={
                "Chrom":"chr",
                "Pos":"pos38",
                "effectAllele":"EFFECT_ALLELE",
                "otherAllele":"OTHER_ALLELE",
                "Beta": "BETA",
                "minus_log10_pval":"minuslog10pval"
                })
            df["chr"] = df["chr"].str.replace("chr", "", regex=False)
            df["chr"] = df["chr"].astype(str).replace({"X": "23", "Y": "24"})
            df["chr"] = df["chr"].astype(int)
            df["SeqID"] = "seq." + df["SeqID"].str.replace("_", ".", regex=False)
            df["PMID"] = 37794188
            df["SAMPLE_SIZE"] = 35892
            df["COHORT"] = "deCODE_2023"
            df["TECHNOLOGY"] = "SOMAscan"
            df["pqtlID"] = (
                df["rsID"].astype(str) + "_" + 
                df["SeqID"].astype(str) + "_" + 
                df["PMID"].astype(str) + 
                df["COHORT"].astype(str)
            )


        # UKB-PPP to literature review format
        if sheet == "pqtl_sun_ukb":
            new_sheet = sheet + "_csa"
            logging.info(f"=== Processing {[str(new_sheet)]} ===")
            logging.info(f"Extracting: {new_sheet}")

            df = pd.read_csv(SUN_UKB_NONEU, sep=";")
            df = df[df["Ancestry"] == "CSA"].copy() # filter for only Central South Asian
            df.reset_index(drop=True, inplace=True)
            split_varID = df["variantID"].str.split(":", expand=True)
            df["pos37"] = split_varID[1]
            df["OTHER_ALLELE"] = split_varID[2]
            df["EFFECT_ALLELE"] = split_varID[3]
            df["PMID"] = 37794186
            df["SAMPLE_SIZE"] = 920
            df["COHORT"] = "UKB_CSA"
            df["TECHNOLOGY"] = "Olink"
            df["Unit"] = "INVRN"
            df["pqtlID"] = df["rsID"].astype(str) + "__" + df["PMID"].astype(str) + "_" + df["COHORT"].astype(str)


        # Meta-analysis to literature review format
        if sheet == "pqtl_sun":
            new_sheet = "pqtl_interval_chris_meta"
            logging.info(f"=== Processing {[str(new_sheet)]} ===")
            logging.info(f"Extracting: {new_sheet}")

            df = pd.read_csv(INTERVAL_CHRIS_META, sep=";",
                             usecols=["SeqID","UniProt","cis_trans","chr","pos37","pos38",
                                      "EFFECT_ALLELE","OTHER_ALLELE","BETA","SE","minuslog10pval"])
            df["PMID"] = 0
            df["SAMPLE_SIZE"] = 13445
            df["COHORT"] = "INTERVAL_CHRIS_META"
            df["TECHNOLOGY"] = "SOMAscan"
            df["pqtlID"] = "__" + df["SeqID"].astype(str) + "__" + df["COHORT"].astype(str)


        # Format, sanity check and save new sheets
        if new_sheet != sheet:
            cohort = new_sheet
            df = format_and_dtype(df, DTYPE_MAP, NUMERIC_COLS)
            df = sanity_check(df, cohort, panels_map, uniprot_check_df, missing_seqid_df, missing_uniprot_df)
            print(f"Writing sheet: {new_sheet}")
            df.to_excel(writer, sheet_name=new_sheet, index=False)


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