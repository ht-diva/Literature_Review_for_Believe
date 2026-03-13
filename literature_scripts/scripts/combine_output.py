from pathlib import Path
import subprocess
import pandas as pd

from paths import PathManager


# ---- PATHS ----
pm = PathManager()
INPUTS_DIR = pm.get_output("literature_gwasstudio_files", exists=False)
OUTPUT_DIR = Path(INPUTS_DIR.parent / "gwasstudio_output")


for cohort_dir in OUTPUT_DIR.iterdir():

    # ---- READ COHORT OUTPUTs & COMBINE ----
    cohort_name = cohort_dir.name
    files = sorted(cohort_dir.glob(f"{cohort_name}_hdsc_believe_*.csv.gz"))
    if not files:
        continue
    print(f"\nCombining outputs for {cohort_name} ({len(files)} files)")
    combined_output = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)


    # ---- MERGE THE FORMATTED INPUT SEARCH TABLE & OUTPUT ----
    formatted_path = INPUTS_DIR / f"{cohort_name}.gwaslab_formatted.csv"
    formatted_input = pd.read_csv(formatted_path)
    merged = formatted_input.merge(
        combined_output,
        on="SOURCEID_SNP",
        how="outer",
        indicator=True
    )


    # ---- SANITY CHECKS ----

    # Check for mismatches in SOURCEID_SNP
    missing_in_combined = merged.loc[merged["_merge"] == "left_only", "SOURCEID_SNP"]
    missing_in_formatted = merged.loc[merged["_merge"] == "right_only", "SOURCEID_SNP"]
    if not missing_in_combined.empty:
        print(
            f"WARNING [{cohort_name}]: "
            f"{len(missing_in_combined)} SOURCEID_SNP missing in output"
        )
    if not missing_in_formatted.empty:
        raise ValueError(
            f"{cohort_name}: {len(missing_in_formatted)} extra SOURCEID_SNP in output"
        )
    missing = merged[merged["_merge"] == "left_only"].drop(columns="_merge")
    merged = merged[merged["_merge"] != "left_only"]
    merged = merged.drop(columns="_merge")

    # Check for mismatches in meta_link_id
    tech = merged["TECHNOLOGY"].dropna().unique()
    if tech[0] != "Olink":
        bad_link = merged[
            merged["meta_link_id"] != merged["SOURCEID_SNP"].str.split(":", n=1).str[0]
        ]
        if not bad_link.empty:
            bad_values = bad_link[["meta_link_id","meta_trait_protein_ids","SOURCEID_SNP","UNIPROT","ORIG_UNIPROT","IS_FLATPROT"]]
            raise ValueError(
                f"{cohort_name}: meta_link_id mismatch for {len(bad_link)} rows\n"
                f"bad_values"
            )

    # Check for mismatches in UniProt
    # Normalize UNIPROTs (order- and space-insensitive)
    # Note: "P0C0L4 | P0C0L5" is the same UniProt as "P0C0L5|P0C0L4"
    u1 = (merged["UNIPROT"].fillna("").str.replace(" ", "", regex=False).str.split("|").apply(lambda x: set(filter(None, x))))
    u2 = (merged["meta_trait_protein_ids"].fillna("").str.replace(" ", "", regex=False).str.split("|").apply(lambda x: set(filter(None, x))))
    exact_match = u1 == u2
    contained = u1.combine(u2, lambda a, b: a.issubset(b))
    error_mask = ~contained
    warn_mask = contained & ~exact_match
    bad_uniprot = merged[error_mask]
    if not bad_uniprot.empty:
        raise ValueError(
            f"{cohort_name}: UNIPROT not contained in meta_trait_protein_ids "
            f"for {len(bad_uniprot)} rows\n"
            f"{bad_uniprot[['UNIPROT', 'meta_trait_protein_ids']]}"
        )
    warn_uniprot = merged[warn_mask]
    if not warn_uniprot.empty:
        print(
            f"WARNING [{cohort_name}]: UNIPROT differs but is contained in "
            f"meta_trait_protein_ids for {len(warn_uniprot)} rows"
        )
    merged["TILEDB_UNIPROT"] = merged["meta_trait_protein_ids"]
    merged["UNIPROT_MATCH"] = merged["meta_trait_protein_ids"] == merged["UNIPROT"]


    # ---- FORMATTING & SAVE ----

    # Append missing IDs
    merged = pd.concat([merged, missing], ignore_index=True)

    # Drop duplicate information
    cols_to_drop = [
        "meta_project",
        "meta_study",
        "meta_category",
        "meta_notes_source_id",
        "meta_build",
        "meta_population",
        "meta_trait_protein_ids",
        "meta_link_id",
    ]
    final = (
        merged
        .drop(columns=cols_to_drop, errors="ignore")
        .rename(columns={"meta_trait_desc": "TRAIT_DESC"})
    )
    final = final.drop_duplicates(keep="first")

    # Save
    output_file = cohort_dir / f"{cohort_name}_hdsc_believe.csv"
    final.to_csv(output_file, index=False)
    print(f"Written to {output_file}")
