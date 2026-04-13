import pandas as pd
import numpy as np
import click
import cloup
import math
import textwrap

from pathlib import Path
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.scalarstring import DoubleQuotedScalarString


help_doc = """
Build the search file (and formatted search table) used for querying metadata
and extract summary statistics with the method extract_regions_leadsnps
"""


@cloup.command(
    "build-search-file-leadsnps",
    no_args_is_help=True,
    help=help_doc,
)
@cloup.option(
    "--search-table",
    required=True,
    help="The search table used to build the search file",
)
@cloup.option(
    "--search-project",
    required=True,
    help="The project to search into (ref. BELIEVE)",
)
@cloup.option(
    "--search-study",
    required=True,
    help="The study to search into (ref. BELIEVE)",
)
@cloup.option(
    "--cohort",
    required=True,
    help="The cohort from the literature review used for the search",
)
@cloup.option(
    "--search-meta-value",
    required=True,
    help="The metadata value to search",
)
@cloup.option(
    "--search-file-prefix",
    default="search_file.yml",
    help="Prefix to be used for naming the search file",
)
@cloup.option(
    "--output-fields",
    default=["build", "population", "trait_protein_ids", "trait_desc"],
    multiple=True,
    help="List of metadata fields to include in the output",
)
@cloup.option(
    "--output-root",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help="Optional root directory where gwasstudio_output will be written. "
         "Defaults to the parent directory of the search table.",
)
@cloup.option(
    "--missing-seqids",
    type=click.Path(path_type=Path, file_okay=True),
    required=True,
    help="Path to SEQIDs missing from BELIEVE.",
)
@cloup.option(
    "--missing-uniprots",
    type=click.Path(path_type=Path, file_okay=True),
    required=True,
    help="Path to UNIPROTs missing from BELIEVE.",
)
@cloup.option(
    "--believe-metadata",
    type=click.Path(path_type=Path, file_okay=True),
    required=True,
    help="Path to BELIEVE annotated Protein Panel.",
)
@cloup.option(
    "--logger-path",
    type=click.Path(path_type=Path, file_okay=True),
    required=True,
    help="Path to logger.",
)
def build_search_file(
    search_table,
    search_project,
    search_study,
    cohort,
    search_meta_value,
    search_file_prefix,
    output_fields,
    output_root,
    missing_seqids,
    missing_uniprots,
    believe_metadata,
    logger_path,
):
    """
    Build the files used to search for trait-specific lead SNPs with the method
    --extract-regions-leadsnps:

        - the search file used for querying metadata
        - the SNP lists specific to the trait/metadata field (search_meta_value)
        - the formatted search table

    Necessary columns in the input table:
        "CHR", "POS", "EA", "NEA"
        plus one of: "UNIPROT", "OLINKID", "SEQID"
    """

    # ---- INPUT CHECKS ----
    if not Path(search_table).is_file():
        raise FileNotFoundError(f"Search table not found: {search_table}")

    if not Path(missing_seqids).is_file():
        raise FileNotFoundError(f"BELIEVE missing SEQIDs table not found: {missing_seqids}")

    if not Path(missing_uniprots).is_file():
        raise FileNotFoundError(f"BELIEVE missing UNIPROTs table not found: {missing_uniprots}")

    if not Path(believe_metadata).is_file():
        raise FileNotFoundError(f"BELIEVE annotated Protein Panel not found: {believe_metadata}")
    

    # ---- LOAD INPUTS ----
    search_table_df = pd.read_csv(search_table, sep="\t")
    tot_var = len(search_table_df)
    missing_var = found_var = flat_prot = 0
    missing_seqids_df = pd.read_csv(missing_seqids, sep="\t")
    missing_seqids_df = missing_seqids_df.loc[missing_seqids_df["COHORT"] == cohort]
    missing_uniprots_df = pd.read_csv(missing_uniprots, sep="\t")
    missing_uniprots_df = missing_uniprots_df.loc[missing_uniprots_df["COHORT"] == cohort]
    believe_metadata_df = pd.read_csv(believe_metadata, sep="\t", usecols=["notes_source_id", "trait_protein_ids"])
    believe_metadata_df = believe_metadata_df.rename(columns={
        "notes_source_id": "SEQID",
        "trait_protein_ids": "UNIPROT",
    })


    # ---- OUTPUT PATHS ----
    search_path = Path(search_file_prefix)
    search_dir = search_path.parent if search_path.parent != Path("") else Path(".")
    search_dir.mkdir(parents=True, exist_ok=True)

    formatted_name = Path(search_table).stem + "_formatted.csv"
    formatted_path = search_dir / formatted_name


    # ---- REQUIRED COLUMNS ----
    required_cols = ["CHR", "POS", "EA", "NEA"]
    if not all(c in search_table_df.columns for c in required_cols):
        raise ValueError("CHR, POS, EA, NEA are missing from the search table.")


    # ---- VALIDATE METADATA ----
    VALID_META = {
        "UNIPROT": "trait_protein_ids",
        "OLINKID": "trait_olink_id",
        "SEQID": "notes_source_id",
    }

    if search_meta_value not in VALID_META:
        raise ValueError(
            f"Invalid search_meta_value '{search_meta_value}'. "
            f"Must be one of: {', '.join(VALID_META.keys())}"
        )

    if search_meta_value not in search_table_df.columns:
        raise ValueError(
            f"Column '{search_meta_value}' is missing from the search table."
        )


    # ---- ADD PROJECT & STUDY ----
    search_table_df["project"] = search_project
    search_table_df["study"] = search_study


    # ---- FORMAT for GWASSTUDIO ----
    search_table_df["CHR"] = (
        search_table_df["CHR"]
        .astype(str)
        .str.replace("chr", "", case=False)
        .replace({"X": "23", "Y": "24"})
    )
    search_table_df = (
        search_table_df[
            search_table_df["CHR"].str.isnumeric()
        ]
        .assign(
            CHR=lambda d: d["CHR"].astype(int),
            POS=lambda d: pd.to_numeric(d["POS"], errors="coerce"),
        )
        .dropna(subset=["POS"])
        .assign(
            POS=lambda d: d["POS"].astype(int),
            EA=lambda d: d["EA"].astype(str),
            NEA=lambda d: d["NEA"].astype(str),
        )
        .reset_index(drop=True)
    )

    # Back-up original SEQIDs and UNIPROTs
    search_table_df["ORIG_SEQID"] = search_table_df["SEQID"].astype(str)
    search_table_df["ORIG_UNIPROT"] = search_table_df["UNIPROT"].astype(str)


    # ---- MISSING SEQIDs VIA UNIPROT ----
    if search_meta_value == "SEQID" and len(missing_seqids_df) > 0:

        # Get missing SEQIDs
        missing_mask = search_table_df["SEQID"].isin(missing_seqids_df["SEQID_MISSING"])
        missing_df = search_table_df.loc[missing_mask]
        search_table_df = search_table_df.loc[~missing_mask]
        missing_var = len(missing_df)
        print(f"   > {cohort}: {missing_var} variants from missing SEQIDs.")

        # Search missing SEQIDs in BELIEVE via UNIPROT
        found_missing_df = missing_df.loc[missing_df["UNIPROT"].isin(believe_metadata_df["UNIPROT"])]
        found_var = len(found_missing_df)
        print(f"   > {cohort}: {found_var} variants from missing SEQIDs found via UNIPROT.")

        # Flatten multi-UNIPROTs if present
        # P0C0L5|P0C0L4 -> P0C0L5 and P0C0L4 in two separate rows
        mask_multi = found_missing_df["UNIPROT"].str.contains(r"\|", na=False)
        multi_df = found_missing_df.loc[mask_multi].copy()
        if mask_multi.sum() > 0:
            multi_df["UNIPROT"] = multi_df["UNIPROT"].str.split(r"\|")
            multi_df = multi_df.explode("UNIPROT").reset_index(drop=True)
            found_missing_df = pd.concat([found_missing_df, multi_df]).reset_index(drop=True)
            flat_prot = len(multi_df)
            print(f"   > {cohort}: {flat_prot} flattened multi-UNIPROTs.")

        # If found, update SEQIDs
        # Note: to one UNIPROT, multiple SEQIDs might be associated
        if len(found_missing_df) > 0:
            merged = found_missing_df.merge(
                believe_metadata_df,
                on="UNIPROT",
                how="left",
                suffixes=("", "_believe")
            )
            mask = (merged["SEQID"] != merged["SEQID_believe"]) & (merged["SEQID_believe"].notna())
            merged["SEQID"] = np.where(mask, merged["SEQID_believe"], merged["SEQID"])
            merged.drop(columns="SEQID_believe", inplace=True)
            found_missing_df = merged
            print(f"   > {cohort}: {len(found_missing_df)} variants from missing SEQIDs found via UNIPROT after adjusting SEQIDs.")

        # Add variants recovered by UNIPROT
        search_table_df = pd.concat([search_table_df, found_missing_df]).drop_duplicates().reset_index(drop=True)


    # ---- HANDLE UNIPROTs ----
    if search_meta_value == "UNIPROT":


        # ---- HANDLE MULTI-UNIPROTs ----

        # Flatten multi-UniProts
        # P0C0L5|P0C0L4 -> P0C0L5 and P0C0L4 in two separate rows
        mask_multi = search_table_df["UNIPROT"].str.contains(r"\|", na=False)
        multi_df = search_table_df.loc[mask_multi].copy()
        multi_df["UNIPROT"] = multi_df["UNIPROT"].str.split(r"\|")
        multi_df = multi_df.explode("UNIPROT").reset_index(drop=True)
        search_table_df = pd.concat([search_table_df, multi_df]).drop_duplicates().reset_index(drop=True)
        flat_prot = len(multi_df)
        print(f"   > {cohort}: {flat_prot} flattened multi-UNIPROTs.")


        # ---- MISSING UNIPROTs ----

        # Drop missing UNIPROTs
        missing_mask = search_table_df["UNIPROT"].isin(missing_uniprots_df["UNIPROT_MISSING"])
        search_table_df = search_table_df.loc[~missing_mask]
        missing_var = missing_mask.sum()
        print(f"   > {cohort}: {missing_var} variants from missing UNIPROTs.")


    # Flags for (un)matching UNIPROTs and SEQIDs
    # SEQID_MATCH == False -> missing SEQID found via UNIPROT
    # UNIPROT_MATCH == False -> flattened multi-UNIPROT
    search_table_df["SEQID_MATCH"] = (
        search_table_df["ORIG_SEQID"] == search_table_df["SEQID"]
    ) | (
        search_table_df["SEQID"].isna()
    )
    search_table_df["UNIPROT_MATCH"] = (
        search_table_df["ORIG_UNIPROT"] == search_table_df["UNIPROT"]
    )


    # ---- BUILD SOURCEID_SNP ----
    search_table_df["SOURCE_ID"] = search_table_df[search_meta_value].astype(str)
    search_table_df["SOURCEID_SNP"] = (
        search_table_df["SOURCE_ID"]
        + ":" + search_table_df["CHR"].astype(str)
        + ":" + search_table_df["POS"].astype(str)
        + ":" + search_table_df["EA"].astype(str)
        + ":" + search_table_df["NEA"].astype(str)
    )


    # ---- WRITE FORMATTED TABLE ----
    search_table_df.to_csv(formatted_path, index=False)


    # ---- BUILD YAML ----
    yaml_data = CommentedMap()
    yaml_data["project"] = search_table_df["project"].iloc[0]
    yaml_data["study"] = search_table_df["study"].iloc[0]

    prefix, subkey = VALID_META[search_meta_value].split("_", 1)
    yaml_data[prefix] = CommentedSeq()

    for meta_value in sorted(search_table_df[search_meta_value].astype(str).unique()):
        yaml_data[prefix].append(CommentedMap({subkey: meta_value}))

    yaml_data["output"] = CommentedSeq(output_fields)

    yaml_data.yaml_set_comment_before_after_key(prefix, before="\n")
    yaml_data.yaml_set_comment_before_after_key("output", before="\n")

    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)

    with open(search_file_prefix, "w") as f:
        yaml.dump(yaml_data, f)


    # ---- BUILD SLURM SCRIPT ----
    cohort_name = Path(search_table).stem.replace(".gwaslab", "")
    run_script = search_dir / f"run_gwasstudio_{cohort_name}.sbatch"
    output_prefix = output_root / cohort_name
    output_prefix.mkdir(parents=True, exist_ok=True)
    output_prefix.chmod(0o2775)

    MIN_PER_GWAS = 1.5
    SAFETY_FACTOR = 1.15
    MAX_TIME_MIN = 1440
    DEFAULT_WORKERS = 4
    MAX_WORKERS = 16

    workers = DEFAULT_WORKERS

    if search_meta_value == "SEQID":
        gwas_nr = search_table_df["SEQID"].nunique()
    elif search_meta_value == "UNIPROT":
        gwas_nr = search_table_df["UNIPROT"].nunique()
    else:
        gwas_nr = len(search_table_df)

    def estimate_time_min(proteins, workers):
        return (proteins / workers) * MIN_PER_GWAS * SAFETY_FACTOR

    est_time = estimate_time_min(gwas_nr, workers)
    while est_time > MAX_TIME_MIN and workers < MAX_WORKERS:
        workers *= 2
        est_time = estimate_time_min(gwas_nr, workers)

    if est_time > MAX_TIME_MIN:
        raise ValueError(
            f"Estimated runtime too long even with {workers} workers "
            f"({est_time:.1f} min for {gwas_nr} proteins)"
        )

    hours = math.ceil(est_time / 60)
    time_str = f"{hours:02d}:00:00"

    print(
        f"Estimated wall-time for {cohort_name} "
        f"({gwas_nr} GWAS, {workers} workers): {time_str}"
    )

    slurm_script = textwrap.dedent(
        f"""\
        #!/bin/bash
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=${{USER}}@fht.org
        #SBATCH --job-name=gwasstudio_leadsnps_{cohort_name}
        #SBATCH --output=logs/%j_gwasstudio_leadsnps_{cohort_name}.log
        #SBATCH --partition=cpuq
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=2G
        #SBATCH --time={time_str}

        source /exchange/healthds/singularity_functions

        gwasstudio --dask-deployment slurm \\
            --workers {workers} \\
            --cores-per-worker 1 \\
            --memory-per-worker 2GiB \\
            export \\
            --search-file {search_file_prefix} \\
            --get-regions-leadsnps {formatted_path} \\
            --output-prefix {output_prefix}/{cohort_name}
        """
    )

    with open(run_script, "w") as f:
        f.write(slurm_script)

    run_script.chmod(0o755)


    # ---- LOGGER ----
    with open(logger_path, 'a') as f:
        f.write(f"{cohort}\t{search_meta_value}\t{tot_var}\t{missing_var}\t{found_var}\t{flat_prot}\n")



if __name__ == "__main__":
    build_search_file()
