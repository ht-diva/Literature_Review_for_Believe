from pathlib import Path
import subprocess
import pandas as pd
import click

from paths import PathManager


# ---- PATHS ----
pm = PathManager()
SEARCH_DIR = pm.get_output("literature_harmonized")
OUTDIR = pm.get_output("literature_gwasstudio_files", exists=False)
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTPUT_ROOT = pm.get_output("gwasstudio_output", exists=False)

SEARCH_PROJECT = "hdsc"
SEARCH_STUDY = "believe"


def main():

    # ---- BUILD GWASSTUDIO FILES ----
    for tsv in sorted(SEARCH_DIR.glob("pqtl*.gwaslab.tsv")):
        cohort = tsv.stem.replace(".gwaslab", "")
        search_file = OUTDIR / f"search_file_{cohort}.yml"

        # Choose METADATA field to use for lead SNP search
        df = pd.read_csv(tsv, sep="\t", usecols=["SEQID", "UNIPROT"])

        seqid_all = df["SEQID"].notna().all()
        uniprot_all = df["UNIPROT"].notna().all()

        if seqid_all:
            search_meta = "SEQID"
        elif uniprot_all:
            search_meta = "UNIPROT"
        else:
            raise ValueError(
                f"{cohort}: Missing metadata — "
                f"SEQID non-NA={df['SEQID'].notna().sum()}, "
                f"UNIPROT non-NA={df['UNIPROT'].notna().sum()}, "
                f"rows={len(df)}"
            )
        print(f"\nMetadata field for {cohort}: {search_meta}")

        cmd = [
            "python",
            "scripts/build_gwasstudio_files.py",
            "--search-table", str(tsv),
            "--search-project", SEARCH_PROJECT,
            "--search-study", SEARCH_STUDY,
            "--search-meta-value", search_meta,
            "--search-file-prefix", str(search_file),
            "--output-root", str(OUTPUT_ROOT),
        ]

        print(f"Running ({cohort}):", " ".join(cmd))
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
