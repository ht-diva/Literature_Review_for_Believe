import pandas as pd
from pathlib import Path
import glob
import re


# ---- LOAD INPUTS & FORMAT ----
path_file = "decode_2023_st19.csv"
base_path = Path(
   "/exchange/healthds/public_data/sumstats/decode/largescaleplasma-2023/final_somascan_smp"
)
df = pd.read_csv(path_file)


# ---- EXTRACT SUMSTATS FILE PATHS ----
# Build file paths (pQTL_ID_prot + gene_prot) and 
# extract supporting information from the deCODE et al. (2023) ST19

path_files = []
total = len(df)
for i, (_, row) in enumerate(df.iterrows(), start=1):
    percent = (i / total) * 100
    print(f"\rProcessed: {percent:.2f}%", end="")

    pqtl_id = row["pQTL ID_prot"]
    gene = row["gene_prot"]
    gene = re.sub(r"^SEP(\d+)", r"SEPT\1", gene)
    gene = gene.replace("-", "_")
    chrom = row["chr"]
    pos = row["pos38"]
    amin = row["Amin"]
    amaj = row["Amaj"]
    beta = row["beta_unadj"]
    mlog10p = row["mLog10pval_unadj"]
    uniprot = row["UniProt"]
    seqid = row["SeqID"]
    cis_trans = row["cis_trans"]
    rsID = row["rsID"]

    parts = pqtl_id.split("_")
    if len(parts) < 2:
        path_files.append((pqtl_id, gene, "invalid_pQTL_ID_format"))
        continue

    prefix = f"{parts[0]}_{parts[1]}_{gene}"
    pattern = base_path / f"Proteomics_SMP_PC0_{prefix}_*.txt.gz"
    matches = glob.glob(str(pattern))

    if not matches:
        path_files.append((pqtl_id, gene, pd.NA, False, uniprot, seqid, cis_trans, rsID, chrom, pos, amin, amaj, beta, mlog10p))
    else:
        full_path = matches[0]
        path_files.append((pqtl_id, gene, full_path, True, uniprot, seqid, cis_trans, rsID, chrom, pos, amin, amaj, beta, mlog10p))


# ---- FORMAT & SAVE ----
path_files_df = pd.DataFrame(
    path_files,
    columns=["pQTL_ID_prot", "gene_prot", "file_name",
             "match_name", "UniProt", "SeqID", "cis_trans", 
             "rsID", "chr", "pos38", "Amin", "Amaj", 
             "beta_unadj", "mLog10pval_unadj"]
)

print(f"\nPath files nr.: {len(path_files_df)}")
path_files_df = path_files_df.drop_duplicates()
path_files_df = path_files_df.sort_values(by="match_name")
path_files_df.to_csv("decode_2023_filenames.csv", sep="\t", index=False)
print(f"Path files nr. (uniq): {len(path_files_df)}")
