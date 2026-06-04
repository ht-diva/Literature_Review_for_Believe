import gzip
import shutil
import subprocess
import pandas as pd
import numpy as np

from pathlib import Path
from ruamel.yaml import YAML
from paths import PathManager
from utils import save_last_commit_id_to_file, make_variant_key


# ---- PATHS & CONFIG ----
pm = PathManager()
LITERATURE_INPUT = pm.get_inputs()["literature_table_cleaned"]
LITERATURE_INPUT_DIR = LITERATURE_INPUT.parent
CONFIGS = pm.get_config()
CONFIG_HARMONIZE = CONFIGS["config_harmonize"]
CONFIG_LIFTOVER_HARMONIZE = CONFIGS["config_liftover_harmonize"]
METADATA = CONFIGS["believe_metadata"]
FORMAT = "literature_rev"
SEP = "\t"
OUTDIR = pm.get_output("literature_harmonized", exists=False)
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTPUT = pm.get_inputs()["literature_table_harmonized"]


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
tmp_harmonize_config = build_tmp_config(
    CONFIG_HARMONIZE,
    "config_harmonize_tmp.yml",
)
tmp_liftover_harmonize_config = build_tmp_config(
    CONFIG_LIFTOVER_HARMONIZE,
    "config_liftover_harmonize_tmp.yml",
)


# ---- HARMONIZE LITERATURE TABLES ----
skip_sheets = {"credits", "variant", "protein", "olink", "cohort", "study"}
xls = pd.ExcelFile(LITERATURE_INPUT)
print(f"INPUT: {LITERATURE_INPUT}")
studies = pd.read_excel(xls, sheet_name="STUDY")
pqtl_studies = "pqtl_" + studies["StudyNAME"].str.lower()
summary_rows = []

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


        # ---- RUN GWASPIPE HARMONIZATION ----
        if refgenome == "GRCh37":
            cmd = [
                "gwaspipe",
                "-f", FORMAT,
                "-s", SEP,
                "-c", str(tmp_liftover_harmonize_config),
                "-i", str(fname),
                "-o", str(OUTDIR)
            ]
            print("Running + Strand Alignment Harmonization:", " ".join(cmd))
        else:
            cmd = [
                "gwaspipe",
                "-f", FORMAT,
                "-s", SEP,
                "-c", str(tmp_harmonize_config),
                "-i", str(fname),
                "-o", str(OUTDIR)
            ]
            print("Running Harmonization:", " ".join(cmd))

        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        gz_out = OUTDIR / f"{cohort}.gwaslab.tsv.gz"
        tsv_out = OUTDIR / f"{cohort}.gwaslab.tsv"
        print("Decompressing into:", tsv_out)
        decompress_gz(gz_out, tsv_out)
        if gz_out.exists():
            gz_out.unlink()
        df_raw = pd.read_csv(fname, sep="\t")
        df_harm = pd.read_csv(tsv_out, sep="\t")


        # ---- COUNT VARIANT LOSS FROM BAD STATISTICS ----
        n_raw = len(df_raw)
        n_harm = len(df_harm)
        loss = n_raw - n_harm


        # ---- REMOVE D-I ALLELES ----
        di_mask = (df_harm["EA"].isin(["D", "I"])) | (df_harm["NEA"].isin(["D", "I"]))
        print(f"D-I allele removal: {di_mask.sum()}")
        df_harm = df_harm.loc[~di_mask]
        df_harm.reset_index(drop=True, inplace=True)


        # ---- COUNT VARIANT LOSS FROM D-I REMOVAL ----
        n_harm_di = n_harm
        n_harm = len(df_harm)
        loss_di = n_harm_di - n_harm


        # ---- ADD UNIQUE VARIANT KEY ----
        if cohort in ["pqtl_QMDiab", "pqtl_interval_chris_meta"]:
            df_raw["merge_key"] = make_variant_key(
                df_raw,
                "pqtlID",
                "chr",
                "pos38",
                "EFFECT_ALLELE",
                "OTHER_ALLELE",
            )

            df_harm["merge_key"] = make_variant_key(
                df_harm,
                "PQTLID",
                "CHR",
                "POS",
                "EA",
                "NEA",
            )


        # ---- BACK-UP SE (pqtl_QMDiab) ----
        if cohort == "pqtl_QMDiab":
            df_harm = df_harm.merge(
                df_raw[["merge_key", "SE_orig"]],
                on="merge_key",
                how="left",
            ).drop_duplicates().reset_index(drop=True)
            df_harm["SE"] = df_harm["SE_orig"]
            df_harm.drop(columns=["merge_key", "SE_orig"], inplace=True)


        # ---- BACK-UP MLOG10P (pqtl_interval_chris_meta) ----
        if cohort == "pqtl_interval_chris_meta":
            df_harm = df_harm.merge(
                df_raw[["merge_key", "MLOG10P_orig"]],
                on="merge_key",
                how="left",
            ).drop_duplicates().reset_index(drop=True)
            df_harm["MLOG10P"] = df_harm["MLOG10P_orig"]
            df_harm.drop(columns=["merge_key", "MLOG10P_orig"], inplace=True)


        # ---- CHECK SEQID & UNIPROT ----
        believe_metadata = pd.read_csv(METADATA, sep="\t")
        believe_metadata = believe_metadata.rename(columns={"notes_source_id": "SEQID"})

        harm_seqids = df_harm["SEQID"].dropna()
        if not harm_seqids.empty:
            harm_seqids = set(harm_seqids)

            # Find not-matching SEQIDs
            believe_seqids = set(believe_metadata["SEQID"].dropna())
            non_matching = [s for s in harm_seqids if s not in believe_seqids]
            if non_matching:
                print(
                    f"WARNING {cohort}: {len(non_matching)} SEQIDs not found. "
                    f"SEQIDs not found: {non_matching[:10]}"
                )

            # Check UniProt consistency for matched SEQIDs
            believe_metadata["trait_protein_ids"] = believe_metadata["trait_protein_ids"].str.strip().str.upper()
            df_harm["UNIPROT"] = df_harm["UNIPROT"].str.strip().str.upper()
            merged = believe_metadata.merge(df_harm, on="SEQID", how="inner")
            uniprot_mismatch = merged[merged["trait_protein_ids"] != merged["UNIPROT"]][["SEQID", "trait_protein_ids", "UNIPROT"]].drop_duplicates()
            if not uniprot_mismatch.empty:
                print(
                    f"WARNING {cohort}: {len(uniprot_mismatch)} SEQIDs with UNIPROT mismatches. "
                    f"SEQIDs with UNIPROT mismatches: {uniprot_mismatch}"
                )


        # ---- SAVE ----
        df_harm = df_harm.drop_duplicates().reset_index(drop=True)
        df_harm.to_csv(tsv_out, sep="\t", index=False)
        print(f"Saving: {tsv_out}")


        # ---- SAVE SHEET ----
        print(f"Writing sheet: {sheet}")
        df_harm.to_excel(writer, sheet_name=sheet, index=False)


        # ---- COUNT VARIANT LOSS FROM LIFTOVER (TOTAL) ----
        loss_liftover = df_harm["POS"].isna().sum() if "POS" in df_harm.columns else np.nan
        if loss_liftover > 0:
            pos37_vals = df_harm.loc[df_harm["POS"].isna(), "POS37"].dropna().unique()
            print(
                f"WARNING {cohort}: POS is NA for {loss_liftover} rows. "
                f"Unique POS37 values: {pos37_vals}"
            )


        # ---- COUNT MULTI-ALLELIC SNPS/LOCI ----
        multiallelic_snps_mask = df_harm.groupby(["CHR", "POS"])["SNPID"].transform("nunique").gt(1)
        nr_multiallelic_snps = multiallelic_snps_mask.sum()
        nr_multiallelic_loci = df_harm.groupby(["CHR", "POS"])["SNPID"].nunique().gt(1).sum()
        if nr_multiallelic_snps > 0:
            multiallelic_snps_df = df_harm[multiallelic_snps_mask][["PQTLID", "SEQID", "UNIPROT", "SNPID"]]
            multiallelic_snps_df = multiallelic_snps_df.drop_duplicates().reset_index(drop=True)
            nr_multiallelic_snps = len(multiallelic_snps_df)


        # ---- SUMMARY ----
        summary_rows.append([
            cohort,
            n_raw,
            n_harm,
            loss,
            loss_di,
            loss_liftover,
            nr_multiallelic_snps,
            nr_multiallelic_loci
        ])


        # ---- CLEAN ----
        print("Removing:", fname)
        fname.unlink()
        


# ---- CLEAN ----
tmp_harmonize_config.unlink(missing_ok=True)
tmp_liftover_harmonize_config.unlink(missing_ok=True)


# ---- SAVE SUMMARY ----
summary_df = pd.DataFrame(
    summary_rows,
    columns=[
        "COHORT",
        "VARIANT_NR_RAW",
        "VARIANT_NR_HARM",
        "VARIANT_LOSS_STATS",
        "VARIANT_LOSS_DI",
        "VARIANT_LOSS_LIFTOVER",
        "MULTI-ALLELIC_SNPS",
        "MULTI-ALLELIC_LOCI"
    ]
)
summary_df.to_csv(OUTDIR / "harmonization_summary.tsv", sep="\t", index=False)
print("\n=== DONE ===")
print(summary_df)

save_last_commit_id_to_file(OUTDIR / "release.txt")
