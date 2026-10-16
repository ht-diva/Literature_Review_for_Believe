import pandas as pd
import re
import time
import requests

from urllib.parse import quote
from html.parser import HTMLParser


class _SearchPageText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []
        self.ignored = False

    def handle_starttag(self, tag, attrs):
        if tag in {"script", "style"}:
            self.ignored = True

    def handle_endtag(self, tag):
        if tag in {"script", "style"}:
            self.ignored = False

    def handle_data(self, data):
        if not self.ignored:
            self.parts.append(data)



# ---- NCBI LIFTOVER VALIDATION FUNCTIONS ----


# Helper function to convert a human RefSeq chromosome accession into a chromosome label
def chromosome_from_accession(sequence_id: str) -> int | None:

    # Mitochondrial reference sequence
    if sequence_id.startswith("NC_012920"):
        return "MT"

    # Extract the accession number, ignoring leading zeros and the version
    # Example: NC_000009.12 -> "9"
    match = re.match(r"NC_0*(\d+)\.", sequence_id)
    if not match:
        return None

    # Canonical CHR
    chromosome_number = int(match.group(1))
    if 1 <= chromosome_number <= 24:
        return chromosome_number

    return None


# Helper function to extract chromosome-position pairs for an assembly from RefSNP JSON
def extract_ncbi_placements(response_data: dict, assembly_prefix: str) -> set[tuple[str, int]]:
    placements = response_data.get("primary_snapshot_data", {}).get("placements_with_allele", [])

    # Deduplicate coordinates shared by multiple alleles or placements
    # Note: Keep placements belonging to the requested assembly and marked as top-level
    coordinates = set()
    for placement in placements:
        annotation = placement.get("placement_annot", {})
        assembly_traits = annotation.get("seq_id_traits_by_assembly", [])
        assembly_matches = any(
            str(trait.get("assembly_name", "")).startswith(assembly_prefix)
            and trait.get("is_top_level", False)
            for trait in assembly_traits
        )
        if not assembly_matches:
            continue

        # Convert the sequence accession into a canonical chromosome label
        sequence_id = placement.get("seq_id", "")
        chromosome = chromosome_from_accession(sequence_id)
        if chromosome is None:
            continue

        # Read the SPDI position for each allele in this placement
        for allele_entry in placement.get("alleles", []):
            spdi = allele_entry.get("allele", {}).get("spdi", {})
            position = spdi.get("position")
            if position is None:
                continue

            # Convert zero-based SPDI coordinate to one-based
            # Note: GWAS coordinates are one-based, so 1 is added
            coordinates.add((chromosome, int(position) + 1))

    return coordinates


# Helper function to fetch NCBI GRCh37 and GRCh38 coordinates (plus fetch status) for one rsID
def fetch_ncbi_coordinates(rsid, session: requests.Session) -> tuple[set[tuple[str, int]], set[tuple[str, int]], str]:
    if pd.isna(rsid):
        return set(), set(), "MISSING_RSID"

    # Get valid rsID
    rsid = str(rsid).strip()
    match = re.fullmatch( r"rs(\d+)", rsid, flags=re.IGNORECASE)
    if not match:
        return (set(), set(), "INVALID_RSID")

    # Search for rsID url
    numeric_rsid = match.group(1)
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{numeric_rsid}"
    try:
        response = session.get(url, timeout=30)
        response.raise_for_status()
        response_data = response.json()

    # Capture errors in query
    except requests.RequestException as error:
        return (set(), set(), f"REQUEST_ERROR: {error}")
    except ValueError:
        return (set(), set(), "INVALID_JSON_RESPONSE")

    # Extract coordinates from rsID url
    grch37_coordinates = extract_ncbi_placements(response_data, "GRCh37")
    grch38_coordinates = extract_ncbi_placements(response_data, "GRCh38")

    return (grch37_coordinates, grch38_coordinates, "OK")


# Helper function to fetch NCBI GRCh37 and GRCh38 coordinates from CHR:GRCh37
def fetch_ncbi_coordinates_by_position(chromosome: str, pos37: int, session: requests.Session) -> tuple[
    set[tuple[str, int]], set[tuple[str, int]], str, str | None
]:
    # Search directly using CHR:POS37
    try:
        response = session.get("https://www.ncbi.nlm.nih.gov/snp/",
            params={"term": f"{chromosome}:{pos37}"},
            timeout=30,
        )
        response.raise_for_status()

    # Capture errors in query
    except requests.RequestException as error:
        return set(), set(), f"SEARCH_ERROR: {error}", None

    # Extract readable text from the HTML
    parser = _SearchPageText()
    parser.feed(response.text)
    text = " ".join(parser.parts)

    # Identify result titles, including merged rsIDs
    titles = list(re.finditer(
        r"rs\d+(?:\s+has\s+merged\s+into\s+rs\d+)?"
        r"\s*\[\s*Homo\s+sapiens\s*\]",
        text,
        flags=re.IGNORECASE,
    ))

    if not titles:
        return set(), set(), "NO_READABLE_SEARCH_RESULT", None

    # Restrict extraction to the first result
    first = titles[0]
    end = titles[1].start() if len(titles) > 1 else len(text)
    first_result = text[first.start():end]
    resolved_rsid = re.findall(r"rs\d+", first.group(), flags=re.IGNORECASE)[-1].lower()

    # Read assembly-labelled coordinates
    coordinates37 = set()
    coordinates38 = set()
    for match in re.finditer(
        r"\b(\d{1,2}|X|Y|MT)\s*:\s*([\d,]+)"
        r"\s*\(\s*GRCh(37|38)\s*\)",
        first_result,
        flags=re.IGNORECASE,
    ):
        chr_value, position, assembly = match.groups()
        coordinate = chr_value.upper(), int(position.replace(",", ""))

        # Displayed coordinates are already one-based: no +1.
        if assembly == "37":
            coordinates37.add(coordinate)
        else:
            coordinates38.add(coordinate)

    # Only accept this result if its GRCh37 coordinate matches the input
    if (chromosome, pos37) not in coordinates37:
        return (coordinates37, set(), "GRCH37_POSITION_NOT_CONFIRMED", resolved_rsid)

    return (coordinates37, coordinates38, "OK_POSITION_LOOKUP", resolved_rsid)


# Main function to validate bcftools and GWASLab bridge mismatched mapped variants cf. NCBI
def add_ncbi_validation(mismatch_df: pd.DataFrame, mismatch_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:

    output_df = mismatch_df.copy()
    session = requests.Session()
    coordinate_cache = {} # avoid querying the same rsID
    validation_rows = []

    for _, row in output_df.iterrows():
        rsid = row["rsID"]
        chromosome = int(row["CHR"]) if pd.notna(row["CHR"]) else None
        pos37 = pd.to_numeric(row["POS37"], errors="coerce")
        missing_rsid = pd.isna(rsid) or not str(rsid).strip()

        # NCBI: Return GRCh37 and GRCh38 coordinates from rsID search url
        if not missing_rsid:
            normalized_rsid = str(rsid).strip().lower()
            cache_key = ("rsid", normalized_rsid)
            if cache_key not in coordinate_cache:
                coordinates37, coordinates38, status = fetch_ncbi_coordinates(normalized_rsid, session)
                coordinate_cache[cache_key] = (coordinates37, coordinates38, status, normalized_rsid)
                time.sleep(0.35)

        # NCBI: Return GRCh37 and GRCh38 coordinates from CHR:POS37 search url
        elif chromosome is not None and pd.notna(pos37):
            cache_key = ("grch37", chromosome, int(pos37))
            if cache_key not in coordinate_cache:
                coordinate_cache[cache_key] = fetch_ncbi_coordinates_by_position(chromosome, int(pos37), session)
                time.sleep(0.35)

        else:
            cache_key = ("missing_coordinates",)
            coordinate_cache[cache_key] = (set(), set(), "MISSING_RSID_AND_COORDINATES", None)

        # Get NCBI coordinates and search status
        coordinates37, coordinates38, request_status, queried_rsid = coordinate_cache[cache_key]

        # Get coordinates from bcftools and GWASLab liftover
        pos37 = pd.to_numeric(row["POS37"], errors="coerce")
        pos38_bcftools = pd.to_numeric(row["POS38_BCFTOOLS"], errors="coerce")
        pos38_gwaslab = pd.to_numeric(row["POS38_GWASLAB"], errors="coerce")

        # Check whether GRCh37 position is found in rsID url
        if pd.notna(pos37):
            expected37 = (chromosome, int(pos37))
            grch37_confirmed = (expected37 in coordinates37)
        else:
            grch37_confirmed = False

        # NCBI GRCh38 placements on the same chromosome
        same_chromosome38 = sorted(position for ncbi_chromosome, position in coordinates38 if ncbi_chromosome == chromosome)

        # Retain all NCBI positions in case dbSNP reports more than one top-level placement
        ncbi_pos38 = ";".join(str(position) for position in same_chromosome38) if same_chromosome38  else pd.NA

        # Check NCBI-liftover coordinates' match
        matches_bcftools = pd.notna(pos38_bcftools) and int(pos38_bcftools) in same_chromosome38
        matches_gwaslab = pd.notna(pos38_gwaslab) and int(pos38_gwaslab) in same_chromosome38
        
        # Report matches
        if not grch37_confirmed:
            position_match = "GRCH37_NOT_CONFIRMED"
        elif matches_bcftools:
            position_match = "BCFTOOLS"
        elif matches_gwaslab:
            position_match = "GWASLAB_BRIDGE"
        elif not same_chromosome38:
            position_match = "NO_NCBI_GRCH38_POSITION"
        else:
            position_match = "NEITHER"

        # Search url
        if chromosome is not None and pd.notna(pos37):
            search_term = f"{int(pos37)}[POSITION_GRCH37] AND {chromosome}[CHR]"
            ncbi_search_url = "https://www.ncbi.nlm.nih.gov/snp/?term=" + quote(search_term)
        else:
            ncbi_search_url = pd.NA

        validation_rows.append({
            "NCBI_REQUEST_STATUS": request_status,
            "NCBI_GRCH37_CONFIRMED": grch37_confirmed,
            "POS38_NCBI": ncbi_pos38,
            "NCBI_POSITION_MATCH": position_match,
            "NCBI_SEARCH_URL": ncbi_search_url,
            "NCBI_QUERIED_RSID": queried_rsid,
            "NCBI_LOOKUP_SOURCE": "GRCH37_POSITION" if missing_rsid else "RSID",
        })

    # Combine output and add NCBI validation fields to individual mismatch rows
    validation_columns = [
        "NCBI_REQUEST_STATUS",
        "NCBI_GRCH37_CONFIRMED",
        "POS38_NCBI",
        "NCBI_POSITION_MATCH",
        "NCBI_SEARCH_URL",
        "NCBI_QUERIED_RSID",
        "NCBI_LOOKUP_SOURCE",
    ]
    validation_df = pd.DataFrame(validation_rows, index=output_df.index, columns=validation_columns)
    validated_rows = pd.concat([output_df, validation_df], axis=1)

    # Count each match category within its cohort
    count_columns = ["NCBI_BCFTOOLS_MATCH", "NCBI_GWASLAB_BRIDGE_MATCH"]
    cohort_counts = (
        validated_rows.assign(
            NCBI_BCFTOOLS_MATCH=validated_rows["NCBI_POSITION_MATCH"].eq("BCFTOOLS"),
            NCBI_GWASLAB_BRIDGE_MATCH=validated_rows["NCBI_POSITION_MATCH"].eq("GWASLAB_BRIDGE"),
        )
        .groupby("COHORT")[count_columns]
        .sum()
    )

    # Update liftover comparison summary
    updated_summary = mismatch_summary.copy()
    for column in count_columns:
        updated_summary[column] = (
            updated_summary["COHORT"]
            .map(cohort_counts[column])
            .fillna(0)
            .astype("int64")
        )

    return validated_rows, updated_summary
