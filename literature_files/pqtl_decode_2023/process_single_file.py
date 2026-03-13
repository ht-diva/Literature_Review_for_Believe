import sys
import pandas as pd
from pathlib import Path


# ---- GWAS FILE PATH ----
fp = sys.argv[1]


# ---- LOAD ST19 WITH MATCHING GWAS FILE PATHS & FORMAT ----
decode2023_orig = pd.read_csv("decode_2023_filenames.csv", sep="\t")
decode2023_orig_filt = decode2023_orig[decode2023_orig["match_name"]].copy()

decode2023_orig_filt["pos38"] = pd.to_numeric(decode2023_orig_filt["pos38"], errors="coerce")
decode2023_orig_filt["chr_num"] = decode2023_orig_filt["chr"].str.replace("chr", "", regex=False)
chr_map = {"X": 23, "Y": 24, "MT": 25}
decode2023_orig_filt["chr_num"] = decode2023_orig_filt["chr_num"].replace(chr_map)
decode2023_orig_filt["chr_num"] = pd.to_numeric(decode2023_orig_filt["chr_num"], errors="coerce")
decode2023_orig_filt = decode2023_orig_filt.sort_values(
    by=["file_name","chr_num","pos38"]
).reset_index(drop=True)


# ---- SUBSET FOR THE CURRENT GWAS ----
decode2023_fp = decode2023_orig_filt.loc[
    decode2023_orig_filt["file_name"] == fp, 
    ["pQTL_ID_prot","gene_prot","UniProt","cis_trans","rsID","chr","pos38",
     "Amin","Amaj","beta_unadj","mLog10pval_unadj"]
].copy()
decode2023_fp["chr"] = decode2023_fp["chr"].astype(str)
decode2023_fp["pos38"] = decode2023_fp["pos38"].astype(int)


# ---- LOAD GWAS ----
gwas_fp = pd.read_csv(
    fp,
    sep="\t",
    usecols=["Chrom","Pos","effectAllele","otherAllele","Beta",
             "minus_log10_pval","SE","ImpMAF"],
    engine="pyarrow"
)

# Merge
merged = gwas_fp.merge(
    decode2023_fp,
    left_on=["Chrom","Pos"],
    right_on=["chr","pos38"],
    how="inner"
)


# ---- SAVE OUTPUT ----
fp_path = Path(fp)
base_name = fp_path.stem.replace(".txt","")
out_fp = Path(f"decode_2023_leadsnp_sumstats/{base_name}_leadsnps.csv")
merged.to_csv(out_fp, index=False)
