from pathlib import Path
import subprocess
import pandas as pd

from paths import PathManager


# ---- SEARCH SETTINGS ----
SEARCH_FIELDS = ["SEQID", "UNIPROT"]
SEARCH_PROJECT = "hdsc"
SEARCH_STUDY = "believe"


# ---- PATHS ----
pm = PathManager()

BELIEVE_METADATA = pm.get_config()["believe_metadata"]

LITERATURE_INPUT = pm.get_inputs()["literature_table_raw"]
LITERATURE_INPUT_DIR = LITERATURE_INPUT.parent
missing_seqid_out = LITERATURE_INPUT_DIR / "seqids_missing.tsv"
missing_uniprot_out = LITERATURE_INPUT_DIR / "uniprots_missing.tsv"

SEARCH_DIR = pm.get_output("literature_harmonized")

OUTDIR = pm.get_output("literature_gwasstudio_files", exists=False)
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTPUT_ROOT = pm.get_output("literature_gwasstudio_output", exists=False)


def run_gwasstudio_file_builder(
    tsv: str,
    cohort: str,
    search_meta: str,
    output_prefix: Path,
    search_file: Path,
    missing_seqids: Path,
    missing_uniprots: Path,
    believe_metadata: Path,
    logger_path: Path,
):
    """Run the script to build GWASStudio files based on the search field."""
    cmd = [
        "python",
        "scripts/build_gwasstudio_files.py",
        "--search-table", str(tsv),
        "--search-project", SEARCH_PROJECT,
        "--search-study", SEARCH_STUDY,
        "--cohort", cohort,
        "--search-meta-value", search_meta,
        "--search-file-prefix", str(search_file),
        "--output-root", str(output_prefix),
        "--missing-seqids", str(missing_seqids),
        "--missing-uniprots", str(missing_uniprots),
        "--believe-metadata", str(believe_metadata),
        "--logger-path", str(logger_path),
    ]

    print(f"\nMetadata field for {cohort}: {search_meta}")
    print(f"Running ({cohort}):", " ".join(cmd))

    output_prefix.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)


def main():


    # ---- LOGGING ----
    log_out = OUTDIR / "literature_summary.tsv"
    if Path(log_out).exists():
        log_out.unlink()
    with open(log_out, 'w') as f:
        f.write("COHORT\tSEARCH_FIELD\tTOT_VAR\tMISSING_VAR\tFOUND_VAR\tFLAT_PROT\n")


    # ---- BUILD GWASSTUDIO FILES ----
    for tsv in sorted(SEARCH_DIR.glob("pqtl*.gwaslab.tsv")):
        cohort = tsv.stem.replace(".gwaslab", "")

        # Choose technology & METADATA field to use for lead SNP search
        df = pd.read_csv(tsv, sep="\t", usecols=["SEQID", "UNIPROT"])
        seqid_all = df["SEQID"].notna().all()
        uniprot_all = df["UNIPROT"].notna().all()

        # No technology found
        if not seqid_all and not uniprot_all:
            raise ValueError(
                f"{cohort}: Missing metadata — "
                f"SEQID non-NA={df['SEQID'].notna().sum()}, "
                f"UNIPROT non-NA={df['UNIPROT'].notna().sum()}, "
                f"rows={len(df)}"
            )

        # SOMAscan technologies
        if seqid_all:
            SEARCH_META = SEARCH_FIELDS[0]
            uniprot_all = False

        # Olink technologies
        else:
            SEARCH_META = SEARCH_FIELDS[1]

        # Make search file
        search_file = OUTDIR / f"search_file_{cohort}.yml"

        # Build GWASStudio files based on the search field
        run_gwasstudio_file_builder(tsv, cohort, SEARCH_META, OUTPUT_ROOT, search_file, 
                                    missing_seqid_out, missing_uniprot_out, BELIEVE_METADATA,
                                    log_out)



if __name__ == "__main__":
    main()
