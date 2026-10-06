import gzip
import shutil
import subprocess
import pandas as pd
import numpy as np
import re
import logging
import os

from pathlib import Path
from ruamel.yaml import YAML
from paths import PathManager
from utils.git import save_last_commit_id_to_file
from utils.helper import extract_removed_number
from utils.missing_ids import missing_seqid_uniprot, plot_missing_ids, report_missing_ids
from utils.variant_loss import check_variant_loss, plot_variant_loss


# ---- PATHS & CONFIG ----
pm = PathManager()
LITERATURE_INPUT = pm.get_inputs()["literature_table_cleaned"]
LITERATURE_INPUT_DIR = LITERATURE_INPUT.parent
CONFIGS = pm.get_config()
CONFIG_HARMONIZE_BUILD38 = CONFIGS["config_harmonize_build38"]
CONFIG_HARMONIZE_BUILD37 = CONFIGS["config_harmonize_build37"]
METADATA = CONFIGS["believe_metadata"]
PANELS_MAP = pm.get_config()["panels_map"]
OUTDIR = pm.get_output("literature_harmonized", exists=False)
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTPUT = pm.get_inputs()["literature_table_harmonized"]

FORMAT = "literature_rev"
SEP = "\t"


# ---- SELECT LIFTOVER CONFIG and OUTPUT (liftover_test) ----

liftover_key = os.environ.get("LIFTOVER_KEY", "default")
allowed_keys = {
    "default",
    "bcftools",
    "gwaslab_standard",
    "gwaslab_bridge",
}
if liftover_key not in allowed_keys:
    raise ValueError(f"Unknown liftover config: {liftover_key}")

if liftover_key != "default":
    CONFIG_HARMONIZE_BUILD37 = CONFIGS[f"config_harmonize_build37_{liftover_key}"]
    OUTPUT = pm.get_outputs()[f"literature_table_harmonized_{liftover_key}"]
    OUTDIR = OUTPUT.parent
    OUTDIR.mkdir(parents=True, exist_ok=True)


# ---- LOGGING ----
log_file = OUTDIR / "literature_harmonized.log"
logging.basicConfig(
    filename=log_file,
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)


# ---- MAP BELIEVE & LITERATURE PANELS ----
panels_map = pd.read_csv(PANELS_MAP, sep="\t")


# ---- HELPER FUNCTIONS ----
def decompress_gz(gz_file: Path, out_file: Path):
    """Decompress .gz to .tsv"""
    with gzip.open(gz_file, "rt") as f_in, open(out_file, "wt") as f_out:
        shutil.copyfileobj(f_in, f_out)

def build_tmp_config(base_config: Path, tmp_name: str) -> Path:
    yaml = YAML()
    yaml.preserve_quotes = True
    with open(base_config) as f:
        config_data = yaml.load(f)
    formatbook_dir = base_config.parent
    config_data["formatbook_path"] = str(formatbook_dir / "formatbook.json")
    tmp_config = Path.cwd() / tmp_name
    with open(tmp_config, "w") as f:
        yaml.dump(config_data, f)
    return tmp_config


# ---- FORMATBOOKS ----
tmp_harmonize_build38_config = build_tmp_config(
    CONFIG_HARMONIZE_BUILD38,
    "config_harmonize_build38_tmp.yml",
)
tmp_harmonize_build37_config = build_tmp_config(
    CONFIG_HARMONIZE_BUILD37,
    "config_harmonize_build37_tmp.yml",
)


# ---- HARMONIZE LITERATURE TABLES ----
skip_sheets = {"credits", "variant", "protein", "olink", "cohort", "study"}
xls = pd.ExcelFile(LITERATURE_INPUT)
print(f"INPUT: {LITERATURE_INPUT}")
studies = pd.read_excel(xls, sheet_name="STUDY")
pqtl_studies = "pqtl_" + studies["StudyNAME"].str.lower()
summary_rows = []
missing_ids_df = []
missing_ids_plot_out = OUTDIR / "missing_ids.png"
variant_loss_plot_out = OUTDIR / "variant_loss.png"

with pd.ExcelWriter(OUTPUT) as writer:
    for sheet in xls.sheet_names:
        df = pd.read_excel(xls, sheet_name=sheet)
        if sheet.lower() in skip_sheets:
            df.to_excel(writer, sheet_name=sheet, index=False)
            continue


        # ---- EXTRACT LITERATURE TABLES ----
        cohort = sheet
        print(f"\n=== Processing {[str(cohort)]} ===")
        print(f"Extracting: {cohort}")


        # ---- REFERENCE GENOME ----
        refgenome = studies.loc[pqtl_studies == sheet.lower(), "ReferenceGenome"].item()
        print(f"{sheet} Reference Genome: {refgenome}")

        # In case of liftover_test, we process only GRCh37
        if liftover_key != "default" and refgenome != "GRCh37":
            print(f"Skip liftover test for {sheet}")
            continue

        # For GRCh37, swap positions (37 <-> 38):
        # For strand alignment with ref. GRCh37, POS is the target (mapped from pos38 in harmonization)
        # For liftover in the next step, pos38 is needed to merge with liftovered output
        if refgenome == "GRCh37":
            pos37 = df["pos37"]
            pos38 = df["pos38"]
            df["pos38"] = pos37
            df["pos37"] = pos38
            print("Swap POS for strand alignment...")

        fname = LITERATURE_INPUT_DIR / f"{sheet}.tsv"
        df.to_csv(fname, sep="\t", index=False)


        # ---- BACK-UP SE (pqtl_QMDiab) ----
        if cohort == "pqtl_QMDiab":
            df_raw = pd.read_csv(fname, sep="\t")
            df_raw["SE_orig"] = df_raw["SE"]
            df_raw.loc[df_raw["SE"] < -1e-07, "SE"] = -0.99e-07 #-1e-07 < SE < inf
            df_raw.to_csv(fname, sep="\t", index=False)


        # ---- BACK-UP MLOG10P (pqtl_interval_chris_meta) ----
        if cohort == "pqtl_interval_chris_meta":
            df_raw = pd.read_csv(fname, sep="\t")
            df_raw["MLOG10P_orig"] = df_raw["minuslog10pval"]
            df_raw.loc[df_raw["minuslog10pval"] > 999.0, "minuslog10pval"] = 999.0 #-1e-07 < MLOG10P < 9999.0000001
            df_raw.to_csv(fname, sep="\t", index=False)


        # ---- Re-calculate MLOG10P from GWASLab (pqtl_CKB) ----
        if cohort in ("pqtl_CKB_SomaScan", "pqtl_CKB_Olink"):
            df_raw = pd.read_csv(fname, sep="\t")
            if "minuslog10pval" in df_raw.columns:
                df_raw.drop(columns="minuslog10pval", inplace=True)
            df_raw.to_csv(fname, sep="\t", index=False)


        # ---- RUN GWASPIPE HARMONIZATION ----
        if refgenome == "GRCh37":
            cmd = [
                "gwaspipe",
                "-f", FORMAT,
                "-s", SEP,
                "-c", str(tmp_harmonize_build37_config),
                "-i", str(fname),
                "-o", str(OUTDIR)
            ]
            print("Running Harmonization for GRCh37:", " ".join(cmd))
        else:
            cmd = [
                "gwaspipe",
                "-f", FORMAT,
                "-s", SEP,
                "-c", str(tmp_harmonize_build38_config),
                "-i", str(fname),
                "-o", str(OUTDIR)
            ]
            print("Running Harmonization for GRCh38:", " ".join(cmd))

        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        #subprocess.run(cmd, check=True)

        gz_out = OUTDIR / f"{cohort}.gwaslab.tsv.gz"
        tsv_out = OUTDIR / f"{cohort}.gwaslab.tsv"
        print("Decompressing into:", tsv_out)
        decompress_gz(gz_out, tsv_out)
        if gz_out.exists():
            gz_out.unlink()
        df_raw = pd.read_csv(fname, sep="\t")
        df_harm = pd.read_csv(tsv_out, sep="\t")


        # ---- COUNT VARIANT LOSS (TOTAL) ----
        n_raw = len(df_raw)
        n_harm = len(df_harm)
        loss_tot = n_raw - n_harm


        # ---- BACK-UP SE (pqtl_QMDiab) ----
        if cohort == "pqtl_QMDiab":
            df_harm = df_harm.merge(
                df_raw[["pqtlID", "chr", "pos38", "SE_orig"]],
                left_on=["PQTLID", "CHR", "POS"],
                right_on=["pqtlID", "chr", "pos38"],
                how="left",
            ).drop_duplicates().reset_index(drop=True)
            df_harm["SE"] = df_harm["SE_orig"]
            df_harm.drop(columns=["pqtlID", "chr", "pos38", "SE_orig"], inplace=True)


        # ---- BACK-UP MLOG10P (pqtl_interval_chris_meta) ----
        if cohort == "pqtl_interval_chris_meta":
            df_harm = df_harm.merge(
                df_raw[["pqtlID", "chr", "pos38", "MLOG10P_orig"]],
                left_on=["PQTLID", "CHR", "POS"],
                right_on=["pqtlID", "chr", "pos38"],
                how="left",
            ).drop_duplicates().reset_index(drop=True)
            df_harm["MLOG10P"] = df_harm["MLOG10P_orig"]
            df_harm.drop(columns=["pqtlID", "chr", "pos38", "MLOG10P_orig"], inplace=True)


        # ---- SAVE TSV ----
        df_harm = df_harm.drop_duplicates().reset_index(drop=True)
        df_harm.to_csv(tsv_out, sep="\t", index=False)
        print(f"Saving: {tsv_out}")


        # ---- SAVE SHEET ----
        print(f"Writing sheet: {sheet}")
        df_harm.to_excel(writer, sheet_name=sheet, index=False)


        # ---- CHECK: EMPTY POS37 ----
        pos_nans = df_harm["POS"].isna().sum() if "POS" in df_harm.columns else np.nan
        if pos_nans > 0:
            pos37_vals = df_harm.loc[df_harm["POS"].isna(), "POS37"].dropna().unique()
            print(
                f"WARNING {cohort}: POS is NA for {pos_nans} rows. "
                f"Unique POS37 values: {pos37_vals}"
            )


        # ---- EXTRACT HARMONIZATION INFO ----

        # Extract log information
        loss_badalleles = None
        loss_badstats = None
        liftover_unmapped = None
        ref_match_log = None
        ref_strand_flip_log = None
        palindromic_snps_log = None

        log_out = tsv_out.with_suffix(".log")
        timestamp_pat = re.compile(r'^[\d:/\s-]+-\s*')

        with log_out.open("r") as fp:
            for line in fp:
                line = timestamp_pat.sub("", line).strip()

                value = extract_removed_number(line, "removed variants with na alleles")
                if value is not None and loss_badalleles is None:
                    loss_badalleles = value
                value = extract_removed_number(line, "variants with bad statistics in total")
                if value is not None and loss_badstats is None:
                    loss_badstats = value
                value = extract_removed_number(line, "unmapped variants")
                if value is not None:
                    liftover_unmapped = value
                value = extract_removed_number(line, "raw matching rate")
                if value is not None:
                    ref_match_log = value
                value = extract_removed_number(line, "variants flipped")
                if value is not None:
                    ref_strand_flip_log = value
                value = extract_removed_number(line, "both allele on genome + unable to distinguish")
                if value is not None:
                    palindromic_snps_log = value

        # Count multi-allelic variants
        multiallelic_snps_mask = df_harm.groupby(["CHR", "POS"])["SNPID"].transform("nunique").gt(1)
        nr_multiallelic_snps = multiallelic_snps_mask.sum()

        # Report missing IDs after harmonization
        missing_seqid_uniprot(df_harm, cohort, panels_map, missing_ids_df)

        # ---- SUMMARY ----
        summary_rows.append([
            cohort,
            n_raw,
            n_harm,
            loss_tot,
            loss_badalleles,
            loss_badstats,
            liftover_unmapped,
            ref_match_log,
            ref_strand_flip_log,
            palindromic_snps_log,
            nr_multiallelic_snps
        ])


        # ---- CLEAN ----
        print("Removing:", fname)
        fname.unlink()
        


# ---- CLEAN ----
tmp_harmonize_build38_config.unlink(missing_ok=True)
tmp_harmonize_build37_config.unlink(missing_ok=True)


# ---- REPORT MISSING IDs ----
if missing_ids_df:
    missing_ids_df = pd.concat(missing_ids_df, ignore_index=True)
    missing_ids_df = missing_ids_df.sort_values(
        by=["COHORT", "PANEL", "VARIANT_NR"],
        ascending=[True, True, False],
    )
    missing_summary_df = report_missing_ids(missing_ids_df)
    plot_missing_ids(missing_summary_df, missing_ids_plot_out)


# ---- SAVE HARMONIZATION SUMMARY ----
summary_df = pd.DataFrame(
    summary_rows,
    columns=[
        "COHORT",
        "VARIANTS_RAW",
        "VARIANTS_HARM",
        "VARIANTS_LOSS_TOTAL",
        "VARIANTS_LOSS_BADALLELES",
        "VARIANTS_LOSS_BADSTATS",
        "VARIANTS_LIFTOVER_UNMAPPED",
        "REF_MATCH",
        "REF_FLIP_VARIANT_NR",
        "REF_PALINDROMIC_NR",
        "VARIANTS_MULTIALLELIC"
    ]
)
summary_df.to_csv(OUTDIR / "harmonization_summary.tsv", sep="\t", index=False)

# ---- REPORT VARIANT LOSS ----
check_variant_loss(summary_df)
plot_variant_loss(summary_df, variant_loss_plot_out)

save_last_commit_id_to_file(OUTDIR / "release.txt")
print("\n=== DONE ===")
