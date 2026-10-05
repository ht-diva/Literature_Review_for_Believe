import pandas as pd
import logging
import re


# ---- HELPER FUNCTIONS ----


# Helper function to make unique identifier (variant key)
def make_variant_key(df, id_col, chr_col, pos_col, a1_col, a2_col):
    alleles = (
        df[[a1_col, a2_col]]
        .astype(str)
        .apply(lambda x: "_".join(sorted(x)), axis=1)
    )

    return (
        df[id_col].astype(str)
        + "_"
        + df[chr_col].astype(str)
        + "_"
        + df[pos_col].astype(str)
        + "_"
        + alleles
    )


# Helper function to align swapped multi-Prots
# Example: P29460|Q9NPF7 <-> Q9NPF7|P29460
def swap_uniprots(df, uniprot1, uniprot2):
    df[uniprot1] = df[uniprot1].fillna("")
    df[uniprot2] = df[uniprot2].fillna("")
    uniprot1_set = df[uniprot1].str.split("|").apply(lambda x: set(sorted(x)))
    uniprot2_set = df[uniprot2].str.split("|").apply(lambda x: set(sorted(x)))
    mask = (
        (uniprot1_set == uniprot2_set) &
        (df[uniprot1] != df[uniprot2])
    )
    return mask


# Helper function to check (multi-)UniProt match
# Example P0DMV8-P16519 matches P0DMV8 and P16519
def is_uniprot_match(row, uniprot1, uniprot2):
    uniprot1set = set(str(row[uniprot1]).split(", ") if pd.notna(row[uniprot1]) else [])
    uniprot2set = set(str(row[uniprot2]).split(", ") if pd.notna(row[uniprot2]) else [])
    return uniprot1set.issubset(uniprot2set) or uniprot2set.issubset(uniprot1set)


# Set column format and data types
def format_and_dtype(df, dtype_map, numeric_cols):

    expected_cols = list(dtype_map.keys())

    # Add missing columns
    missing_cols = [c for c in expected_cols if c not in df.columns]
    df = df.assign(**{c: pd.NA for c in missing_cols})

    # Reorder columns: expected first, extras last
    ordered_cols = expected_cols + [c for c in df.columns if c not in expected_cols]
    df = df[ordered_cols]

    # Format chromosome
    df["chr"] = (
        df["chr"]
        .astype(str)
        .replace({"X": "23", "Y": "24"})
    )

    # Numeric coercion
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(
                df[col].astype(str).str.replace(",", ".", regex=False),
                errors="coerce"
            )

    # Apply dtypes
    for col, dtype in dtype_map.items():
        if col in df.columns:
            try:
                df[col] = df[col].astype(dtype)
            except Exception as e:
                logging.warning(f"Could not cast {col} to {dtype}: {e}")

    return df


# Helper function to find the best ID to match BELIEVE and Literature Protein Panels
# IDs are:
#   "Target_Name"
#   "SeqID"
#   "Ensembl_Gene_ID"
#   "UniProt"
def best_id_match(believe_metadata, literature_panel):

    # Count matches
    matches = {
        "Target_Name": 0,
        "SeqID": 0,
        "Ensembl_Gene_ID": 0,
        "UniProt": 0
    }

    for col in matches.keys():
        merged = pd.merge(
            believe_metadata[[col]],
            literature_panel[[col]],
            on=col,
            how="inner"
        ).drop_duplicates()
        matches[col] = len(merged)

    for col, count in matches.items():
        logging.info(f"> Number of matches for {col}: {count} on {len(believe_metadata)}.")

    # Best match
    best_match = max(matches, key=matches.get)
    best_count = matches[best_match]
    logging.info(f"> ID with the maximum matches: {best_match} ({best_count} matches).")

    # Missing match
    all_matches = pd.merge(
        believe_metadata[[best_match]],
        literature_panel[[best_match]],
        on=best_match,
        how="inner"
    ).drop_duplicates()
    missing_matches = believe_metadata[~believe_metadata[best_match].isin(all_matches[best_match])]

    if len(missing_matches) > 0:
        logging.info(f"> Missing {best_match}: {missing_matches}.")
    else:
        logging.info(f"> All matched by {best_match}.")

    return best_match


# Helper function to extract number for log lines
def extract_removed_number(line, description):

    number_pattern = r"(?P<value>\d[\d,]*(?:\.\d+)?)"

    description = re.escape(description)
    patterns = [
        # Removed 10 variants with bad statistics...
        rf"removed\s+{number_pattern}\s+{description}",

        # Raw matching rate: 99.5%
        # Removed variants with NA alleles or ...: 11
        rf"{description}[^:\n]*:\s*{number_pattern}\s*%?",
    ]

    for pattern in patterns:
        match = re.search(pattern, line, flags=re.IGNORECASE)
        if match:
            return float(match.group("value").replace(",", ""))

    return None
