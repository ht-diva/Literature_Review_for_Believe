import pandas as pd
import numpy as np
import logging
import subprocess
import pysam
import gzip

from paths import PathManager
from utils import write_vcf


# ---- PATHS ----
pm = PathManager()
LITERATURE_INPUT = pm.get_inputs()["literature_table_harmonized"]
OUTPUT = pm.get_inputs()["literature_table_liftover"]
OUTDIR = OUTPUT.parent
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTDIR_HARM = pm.get_output("literature_harmonized", exists=False)
singularity_image = "/ssu/gassu/singularity/bcftools_latest.sif"


# ---- HELPER FUNCTIONS ----
def run(cmd):
    subprocess.run(cmd, check=True)

def bcftools(image, args):
    return ["singularity", "run", image, "bcftools"] + [str(a) for a in args]

def run(cmd, stdout=None):
    subprocess.run(cmd, check=True, stdout=stdout)

def bcftools(img, args):
    return ["singularity", "exec", img, "bcftools"] + args


# ---- REFERENCES ----
hg19_fa = "/group/diangelantonio/public_data/liftOver/human_g1k_v37.fasta"
hg19_fai = "/group/diangelantonio/public_data/liftOver/human_g1k_v37.fasta.fai"
hg38_fa = "/group/diangelantonio/public_data/liftOver/hg38.fa"
hg38_fai = "/group/diangelantonio/public_data/liftOver/hg38.fa.fai"
chain = "/group/diangelantonio/public_data/liftOver/hg19ToHg38.over.chain.gz"


# ---- LIFTOVER PIPELINE ----
skip_sheets = {"credits", "variant", "protein", "olink", "cohort", "study"}
xls = pd.ExcelFile(LITERATURE_INPUT)
studies = pd.read_excel(xls, sheet_name="STUDY")
pqtl_studies = "pqtl_" + studies["StudyNAME"].str.lower()

with pd.ExcelWriter(OUTPUT) as writer:
    for sheet in xls.sheet_names:
        df = pd.read_excel(xls, sheet_name=sheet)
        if sheet.lower() in skip_sheets:
            df.to_excel(writer, sheet_name=sheet, index=False)
            continue

        cohort = sheet
        logging.info(f"=== Processing {[str(cohort)]} ===")
        logging.info(f"Extracting: {cohort}")

        # Reference Genome
        refgenome = studies.loc[pqtl_studies == sheet.lower(), "ReferenceGenome"].item()
        print(f"{sheet} Reference Genome: {refgenome}")

        # Define file names
        vcf_file = OUTDIR / f"{sheet}.vcf"
        stem = vcf_file.stem
        vcf_gz = OUTDIR / f"{stem}.vcf.gz"
        reheader_vcf = OUTDIR / f"{stem}.reheader.vcf.gz"
        sorted_vcf = OUTDIR / f"{stem}.sorted.vcf.gz"
        standard_vcf = OUTDIR / f"{stem}.standard.vcf.gz"
        liftover_vcf = OUTDIR / f"{stem}.liftover.vcf"
        liftover_sorted = OUTDIR / f"{stem}.liftover.vcf.gz"
        fixref_vcf = OUTDIR / f"{stem}.liftover.fixref.vcf.gz"
        log_file = OUTDIR_HARM / f"{stem}.liftover.log"


        # ---- GRCh38: CHECK REF-CONSISTENCY ONLY ----
        if refgenome in ["GRCh38", "GRCh37/GRCh38"]:
            df.to_excel(writer, sheet_name=sheet, index=False)
            print(" ...Only Ref-Consistency Check.")

            # 0. Convert to VCF
            write_vcf(df, vcf_file, build="GRCh38")

            # 1. BGZIP
            run(bcftools(singularity_image, [
                "view", "-Oz",
                "-o", vcf_gz,
                vcf_file
            ]))

            # 2. Re-header
            run(bcftools(singularity_image, [
                "reheader",
                "-f", hg38_fai,
                "-o", reheader_vcf,
                vcf_gz
            ]))

            # 3. Sort
            run(bcftools(singularity_image, [
                "sort",
                "-Oz",
                "-o", sorted_vcf,
                reheader_vcf
            ]))
            pysam.tabix_index(str(sorted_vcf), preset="vcf", force=True)

            # 4. Norm hg38
            run(bcftools(singularity_image, [
                "norm",
                "-f", hg38_fa,
                "-c", "s",
                "-Oz",
                "-o", standard_vcf,
                sorted_vcf
            ]))

            # 5. Check Ref-consistency
            with open(log_file, "w") as log:
                subprocess.run(
                    bcftools(singularity_image, [
                        "+fixref",
                        standard_vcf,
                        "-Oz",
                        "-o", fixref_vcf,
                        "--",
                        "-f", hg38_fa
                    ]),
                    stdout=log,
                    stderr=log,
                    check=True
                )
            continue


        # ---- CONVERT TO VCF FOR LIFTOVER ----
        write_vcf(df, vcf_file)


        # ---- GRCh37: LIFTOVER STEPS ----

        # 1. BGZIP
        run(bcftools(singularity_image, [
            "view", "-Oz",
            "-o", vcf_gz,
            vcf_file
        ]))

        # 2. Re-header
        run(bcftools(singularity_image, [
            "reheader",
            "-f", hg19_fai,
            "-o", reheader_vcf,
            vcf_gz
        ]))

        # 3. Sort
        run(bcftools(singularity_image, [
            "sort",
            "-Oz",
            "-o", sorted_vcf,
            reheader_vcf
        ]))
        pysam.tabix_index(str(sorted_vcf), preset="vcf", force=True)

        # 4. Norm hg19
        run(bcftools(singularity_image, [
            "norm",
            "-f", hg19_fa,
            "-c", "s",
            "-Oz",
            "-o", standard_vcf,
            sorted_vcf
        ]))

        # 5. Liftover Pipe → Sort
        with open(log_file, "w") as log:

            p1 = subprocess.Popen(
                bcftools(singularity_image, [
                    "+liftover",
                    "--no-version",
                    "-Ou",
                    standard_vcf,
                    "--",
                    "-s", hg19_fa,
                    "-f", hg38_fa,
                    "-c", chain
                ]),
                stdout=subprocess.PIPE,
                stderr=log
            )

            p2 = subprocess.Popen(
                bcftools(singularity_image, [
                    "sort",
                    "-Oz",
                    "-o", liftover_sorted
                ]),
                stdin=p1.stdout
            )

            p1.stdout.close()
            p2.communicate()

        pysam.tabix_index(str(liftover_sorted), preset="vcf", force=True)

        # 6. Check Ref-consistency
        with open(log_file, "a") as log:
            subprocess.run(
                bcftools(singularity_image, [
                    "+fixref",
                    liftover_sorted,
                    "-Oz",
                    "-o", fixref_vcf,
                    "--",
                    "-f", hg38_fa
                ]),
                stdout=log,
                stderr=log,
                check=True
            )


        # ---- LOAD AND FORMAT LIFTOVERED FILES ----
        
        # Read the liftovered VCF file
        vcf_rows = []
        with gzip.open(liftover_sorted, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    vcf_rows.append(line.strip().split("\t")[:5])

        vcf_df = pd.DataFrame(vcf_rows, columns=["CHROM", "POS", "ID", "REF", "ALT"]).drop_duplicates()
        vcf_df["CHROM"] = (
            vcf_df["CHROM"]
            .str.replace("^chr", "", regex=True)
            .replace({"X": "23", "Y": "24"})
        )
        df["CHR"] = df["CHR"].astype(str)
        df["POS"] = df["POS"].astype(str)
        df["POS37"] = df["POS37"].astype(str)
        df["rsID"] = df["rsID"].astype(str)
        vcf_df["POS"] = vcf_df["POS"].astype(str)
        vcf_df["ID"] = vcf_df["ID"].astype(str)

        vcf_df = vcf_df.rename(columns={
            "CHROM": "CHR_lift",
            "POS": "POS_lift",
            "ID": "rsID",
        })

        # Merge
        merged_df = df.merge(
            vcf_df,
            left_on=["CHR", "POS37", "rsID"],
            right_on=["CHR_lift", "POS_lift", "rsID"],
            how="left"
        )

        # Allele formatting
        merged_df["REF"] = merged_df["REF"].replace("nan", np.nan)
        merged_df["ALT"] = merged_df["ALT"].replace("nan", np.nan)
        for col in ["NEA", "EA", "REF", "ALT"]:
            merged_df[col] = merged_df[col].astype(str).str.upper().str.strip()

        # Update SNPs
        merged_df = (
            merged_df
            .assign(
                CHR=merged_df["CHR_lift"].combine_first(merged_df["CHR"]),
                POS37=merged_df["POS"],
                POS=merged_df["POS_lift"].combine_first(merged_df["POS"]),
                SNPID=merged_df["CHR"] + ":" + merged_df["POS"] + merged_df["EA"] + merged_df["NEA"],
            )
            .drop(columns=["CHR_lift", "POS_lift", "REF", "ALT"])
        )

        # Save the liftovered formatted file
        merged_df.to_excel(writer, sheet_name=sheet, index=False)

        # Update harmonized files for GWASStudio file built
        tsv_out = OUTDIR_HARM / f"{cohort}.gwaslab.tsv"
        merged_df.to_csv(tsv_out, sep="\t", index=False)


        # ---- CLEAN ----

        # Remove all intermediate files
        paths_to_clean = [
            vcf_file,
            vcf_gz,
            reheader_vcf,
            sorted_vcf,
            standard_vcf,
            liftover_vcf,
            liftover_sorted,
        ]

        for p in paths_to_clean:
            if p.exists():
                print("Removing:", p)
                p.unlink()

        # Remove all .tbi files
        for p in OUTDIR.glob("*.tbi"):
            print("Removing:", p)
            p.unlink()
