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
LITERATURE_INPUT = pm.get_inputs()["literature_table"]
OUTPUT = pm.get_inputs()["literature_table_liftover"]
OUTDIR = OUTPUT.parent
OUTDIR.mkdir(parents=True, exist_ok=True)
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

        # Skip GRCh38
        refgenome = studies.loc[pqtl_studies == sheet.lower(), "ReferenceGenome"].item()
        print(f"{sheet} Reference Genome: {refgenome}")
        if refgenome == "GRCh38":
            df.to_excel(writer, sheet_name=sheet, index=False)
            print(" ...skip liftover.")
            continue

        cohort = sheet
        logging.info(f"=== Processing {[str(cohort)]} ===")
        logging.info(f"Extracting: {cohort}")


        # ---- CONVERT TO VCF FOR LIFTOVER ----
        vcf_file = OUTDIR / f"{sheet}.vcf"
        write_vcf(df, vcf_file)


        # ---- LIFTOVER STEPS ----

        # Define file names
        stem = vcf_file.stem
        vcf_gz = OUTDIR / f"{stem}_to_convert_from37_to38.vcf.gz"
        reheader_vcf = OUTDIR / f"{stem}_to_convert_from37_to38.reheader.vcf.gz"
        sorted_vcf = OUTDIR / f"{stem}_to_convert_from37_to38.sorted.vcf.gz"
        standard_vcf = OUTDIR / f"{stem}_to_convert_from37_to38_standard.vcf.gz"
        liftover_vcf = OUTDIR / f"{stem}_to_convert_from37_to38_liftover.vcf"
        liftover_txt = OUTDIR / f"{stem}_liftover.vcf.txt"
        liftover_sorted = OUTDIR / f"{stem}_to_convert_from37_to38_liftover.vcf.gz"
        norm_split = OUTDIR / f"{stem}_to_convert_from37_to38_liftover.hg38.norm.split.vcf.gz"
        atom_vcf = OUTDIR / f"{stem}_to_convert_from37_to38_liftover.hg38.norm.split.atom.vcf.gz"

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

        # 5. Liftover
        with open(liftover_vcf, "wb") as out:
            run(bcftools(singularity_image, [
                "+liftover",
                "--no-version",
                "-Ou",
                standard_vcf,
                "--",
                "-s", hg19_fa,
                "-f", hg38_fa,
                "-c", chain
            ]), stdout=out)

        # 6. View → txt
        with open(liftover_txt, "wb") as out:
            run(bcftools(singularity_image, [
                "view",
                liftover_vcf
            ]), stdout=out)

        # 7. Liftover Pipe → Sort
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
            stdout=subprocess.PIPE
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

        # 8. Norm hg38 + Split
        run(bcftools(singularity_image, [
            "norm",
            "-f", hg38_fa,
            "-c", "s",
            "-m", "-both",
            "-Oz",
            "-o", norm_split,
            liftover_sorted
        ]))
        pysam.tabix_index(str(norm_split), preset="vcf", force=True)

        # 9. Atomize
        run(bcftools(singularity_image, [
            "norm",
            "-f", hg38_fa,
            "-a",
            "-Oz",
            "-o", atom_vcf,
            norm_split
        ]))
        pysam.tabix_index(str(atom_vcf), preset="vcf", force=True)


        # ---- LOAD AND FORMAT LIFTOVERED FILES ----
        
        # Read the liftovered VCF file
        vcf_rows = []
        with gzip.open(atom_vcf, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    vcf_rows.append(line.strip().split("\t")[:5])

        vcf_df = pd.DataFrame(vcf_rows, columns=["CHROM", "POS", "ID", "REF", "ALT"]).drop_duplicates()
        
        vcf_df["CHROM"] = vcf_df["CHROM"].str.replace("^chr", "", regex=True)
        df["chr"] = df["chr"].astype(str)
        df["pos38"] = df["pos38"].astype(str)
        df["rsID"] = df["rsID"].astype(str)
        vcf_df["POS"] = vcf_df["POS"].astype(str)
        vcf_df["ID"] = vcf_df["ID"].astype(str)

        # Merge
        merged_df = df.merge(
            vcf_df,
            left_on=["chr", "pos38", "rsID"],
            right_on=["CHROM", "POS", "ID"],
            how="left"
        )

        # Allele formatting
        merged_df["REF"] = merged_df["REF"].replace("nan", np.nan)
        merged_df["ALT"] = merged_df["ALT"].replace("nan", np.nan)
        for col in ["OTHER_ALLELE", "EFFECT_ALLELE", "REF", "ALT"]:
            merged_df[col] = merged_df[col].astype(str).str.upper().str.strip()

        # Allele status
        merged_df["allele_status"] = np.select(
            [
                merged_df["REF"].isna() | merged_df["ALT"].isna(),
                (merged_df["OTHER_ALLELE"] == merged_df["REF"]) &
                (merged_df["EFFECT_ALLELE"] == merged_df["ALT"]),
                (merged_df["OTHER_ALLELE"] == merged_df["ALT"]) &
                (merged_df["EFFECT_ALLELE"] == merged_df["REF"]),
            ],
            ["no_match_in_vcf", "match", "swapped"],
            default="different"
        )

        # Update SNPs
        merged_df = (
            merged_df
            .assign(
                CHROM_orig=merged_df["chr"],
                POS_orig=merged_df["pos37"],
                OTHER_ALLELE_orig=merged_df["OTHER_ALLELE"],
                EFFECT_ALLELE_org=merged_df["EFFECT_ALLELE"],
                chr=merged_df["CHROM"],
                pos37=merged_df["POS"],
                OTHER_ALLELE=merged_df["REF"],
                EFFECT_ALLELE=merged_df["ALT"],
            )
            .drop(columns=["CHROM", "POS", "REF", "ALT"])
        )

        # Flip the sign of BETA for swapped alleles
        merged_df.loc[merged_df["allele_status"] == "swapped", "BETA"] *= -1
        
        # Save the liftovered formatted file
        merged_df.to_excel(writer, sheet_name=sheet, index=False)


        # ---- CLEAN ----

        # Remove all intermediate files
        paths_to_clean = [
            vcf_file,
            vcf_gz,
            reheader_vcf,
            sorted_vcf,
            standard_vcf,
            liftover_vcf,
            liftover_txt,
            liftover_sorted,
            norm_split,
            atom_vcf,
        ]

        for p in paths_to_clean:
            if p.exists():
                print("Removing:", p)
                p.unlink()

        # Remove all .tbi files
        for p in OUTDIR.glob("*.tbi"):
            print("Removing:", p)
            p.unlink()
