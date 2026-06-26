import pandas as pd
import numpy as np
import logging

from paths import PathManager
from utils import format_and_dtype


# ---- PATHS ----
pm = PathManager()
LITERATURE_INPUT = pm.get_inputs()["literature_table_raw"]
LITERATURE_INPUT_DIR = LITERATURE_INPUT.parent
LITERATURE_FILES = pm.get_files()
SUN_UKB_NONEU = LITERATURE_FILES["pqtl_sun_ukb_csa"]
INTERVAL_CHRIS_META = LITERATURE_FILES["pqtl_interval_chris_meta"]
DECODE_2023 = LITERATURE_FILES["pqtl_decode_2023"]
CKB_SOMASCAN = LITERATURE_FILES["pqtl_CKB_SomaScan"]
CKB_OLINK = LITERATURE_FILES["pqtl_CKB_Olink"]
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


# ---- READ LITERATURE TABLES ----
skip_sheets = {"credits", "variant", "protein", "olink", "cohort", "study"}
xls = pd.ExcelFile(LITERATURE_INPUT)

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


        # CKB SomaScan to literature review format
        if sheet == "pqtl_Brain":
            new_sheet = "pqtl_CKB_SomaScan"
            logging.info(f"=== Processing {[str(new_sheet)]} ===")
            logging.info(f"Extracting: {new_sheet}")

            df = pd.read_csv(CKB_SOMASCAN, sep=";",
                             usecols=["SeqID","cis_trans","chr","pos38",
                                      "EFFECT_ALLELE","OTHER_ALLELE","BETA","SE","P"])
            df["minuslog10pval"] = np.where(
                df["P"] > 0,
                -np.log10(df["P"]),
                pd.NA
            )
            df.drop(columns=["P"], inplace=True)
            df["pos37"] = pd.NA
            df["PMID"] = 39984443
            df["SAMPLE_SIZE"] = 3976
            df["COHORT"] = "CKB_SOMASCAN"
            df["TECHNOLOGY"] = "SOMAscan"
            df["pqtlID"] = (
                "__" + 
                df["SeqID"].astype(str) + "_" + 
                df["PMID"].astype(str) + 
                df["COHORT"].astype(str)
            )


        # CKB Olink to literature review format
        if sheet == "pqtl_CSF":
            new_sheet = "pqtl_CKB_Olink"
            logging.info(f"=== Processing {[str(new_sheet)]} ===")
            logging.info(f"Extracting: {new_sheet}")

            # Olink panel
            olink_panel = pd.read_excel(
                xls,
                sheet_name="OLINK",
                usecols=["Target_Name", "OlinkID", "UniProt"]
            ).drop_duplicates(subset="Target_Name")

            # Load original table
            df = pd.read_csv(CKB_OLINK, sep=";",
                             usecols=["Target_Name","cis_trans","chr","pos38",
                                      "EFFECT_ALLELE","OTHER_ALLELE","BETA","SE","P"])

            # Add UniProt from panel
            df = df.merge(olink_panel, on="Target_Name", how="left", validate="many_to_one")
            missing = df[df["OlinkID"].isna()]
            print(f"{len(missing)} genes not found in the OLINK panel.")
            if not missing.empty:
                print("Missing Target_Name:")
                print(missing["Target_Name"].drop_duplicates().sort_values().to_list())
            df.drop(columns=["Target_Name"], inplace=True)

            # Add other missing statistics
            df["minuslog10pval"] = np.where(
                df["P"] > 0,
                -np.log10(df["P"]),
                pd.NA
            )
            df.drop(columns=["P"], inplace=True)
            df["pos37"] = pd.NA
            df["PMID"] = 39984443
            df["SAMPLE_SIZE"] = 3976
            df["COHORT"] = "CKB_OLINK"
            df["TECHNOLOGY"] = "Olink"
            df["pqtlID"] = (
                "__" + 
                df["OlinkID"].astype(str) + "_" + 
                df["PMID"].astype(str) + 
                df["COHORT"].astype(str)
            )


        # Format, sanity check and save new sheets
        if new_sheet != sheet:
            cohort = new_sheet
            df = format_and_dtype(df, DTYPE_MAP, NUMERIC_COLS)
            print(f"Writing sheet: {new_sheet}")
            df.to_excel(writer, sheet_name=new_sheet, index=False)
