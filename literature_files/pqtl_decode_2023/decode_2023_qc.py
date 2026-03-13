import pandas as pd
import numpy as np


# ---- QC STATISTICS THRESHOLD ----
# Maximum allowed total difference of BETA and MLOG10P between:
# lead SNPs published in ST19 and
# original summary statistics
MISMATCH_THR = 1.0


# ---- READ AND FORMAT ----
df = pd.read_csv("decode_2023_leadsnp_sumstats/pqtl_decode_2023_all_chr_pos.csv")
df_fn = pd.read_csv("decode_2023_filenames.csv", sep="\t")
num_cols = ["Beta", "beta_unadj", "minus_log10_pval", "mLog10pval_unadj"]
for c in num_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# Add seqID
df = df_fn[["UniProt","SeqID", "gene_prot", "pQTL_ID_prot"]].drop_duplicates().merge(
    df,
    on=["UniProt", "gene_prot", "pQTL_ID_prot"],
    how="right",
    indicator=True
).drop(columns=["_merge"])


# ---- STATISTICS DIFFERENCE ----
df["beta_diff"] = (df["Beta"].abs() - df["beta_unadj"].abs()).abs()
raw_pval_diff = (df["minus_log10_pval"] - df["mLog10pval_unadj"]).abs()
inf_mask = (np.isinf(df["minus_log10_pval"]) | np.isinf(df["mLog10pval_unadj"]))
df["pval_diff"] = np.where(inf_mask, 0.0, raw_pval_diff) # if MLOG10P is inf, keep 0.0 difference
df["total_diff"] = df["beta_diff"] + df["pval_diff"]

# Keep closest statistics per Chrom-Pos
df_sorted = df.sort_values(
    by=["gene_prot", "pQTL_ID_prot", "Chrom", "Pos", "total_diff", "rsID", "Amin", "Amaj"],
    ascending=[True, True, True, True, True, True, True, True]
)
df_final = (
    df_sorted
    .drop_duplicates(subset=["gene_prot", "pQTL_ID_prot", "Chrom", "Pos", "rsID", "Amin", "Amaj"], keep="first")
)


# ---- QC MISMATCHES ----

# Statistical mismatches after filtering?
mismatch_filt_mask = df_final["total_diff"] > MISMATCH_THR
df_mismatch_filt = df_final.loc[mismatch_filt_mask].copy()
df_mismatch_filt = df_mismatch_filt[
    ["gene_prot", "pQTL_ID_prot", "Chrom", "Pos", "rsID", 
     "effectAllele", "otherAllele", "Amin", "Amaj", "Beta", "minus_log10_pval", "total_diff"]].drop_duplicates()
nr_stat_mismatch = len(df_mismatch_filt)

if nr_stat_mismatch > 0:
    df_mismatch_filt = df_fn.merge(
        df_mismatch_filt,
        left_on=["pQTL_ID_prot", "gene_prot", "chr", "pos38", "rsID", "Amin", "Amaj"],
        right_on=["pQTL_ID_prot", "gene_prot", "Chrom", "Pos", "rsID", "Amin", "Amaj"],
        how="right",
        indicator=True
    ).query('_merge == "both"').drop(columns=["_merge", "Chrom", "Pos"])
    df_mismatch_filt.to_csv("decode_2023_stats_mismatch.csv", sep="\t", index=False)

# Drop statistical mismatches
df_final = df_final.loc[~mismatch_filt_mask]
print(f"\nDropped variants:")
print(f"    |-> Statistical mismatch: {nr_stat_mismatch}")

# General allele check
def assign_status(row):
    effect_allele = row["effectAllele"]
    other_allele = row["otherAllele"]
    amin = row["Amin"]
    amaj = row["Amaj"]

    if ((effect_allele == amin and other_allele == amaj) or
        (effect_allele == amaj and other_allele == amin)):
        return "OK"
    elif ((effect_allele == amin and amaj in ["!", "*"]) or
          (effect_allele == amaj and amin in ["!", "*"]) or
          (other_allele == amin and amaj in ["!", "*"]) or
          (other_allele == amaj and amin in ["!", "*"])):
        return "MISSING_1"
    elif amin in ["!", "*"] and amaj in ["!", "*"]:
        return "MISSING_ALL"
    else:
        return "MISMATCH"

# Check for allele mismatches
df_final["STATUS"] = df_final.apply(assign_status, axis=1)
df_mismatch_alleles = df_final.loc[df_final["STATUS"] == "MISMATCH"].copy()
df_mismatch_alleles = df_mismatch_alleles[
    ["gene_prot", "pQTL_ID_prot", "Chrom", "Pos", "rsID", 
     "effectAllele", "otherAllele", "Amin", "Amaj", "Beta", "minus_log10_pval", "total_diff"]].drop_duplicates()

if len(df_mismatch_alleles) > 0:
    df_mismatch_alleles = df_fn.merge(
        df_mismatch_alleles,
        left_on=["pQTL_ID_prot", "gene_prot", "chr", "pos38", "rsID", "Amin", "Amaj"],
        right_on=["pQTL_ID_prot", "gene_prot", "Chrom", "Pos", "rsID", "Amin", "Amaj"],
        how="right",
        indicator=True
    ).query('_merge == "both"').drop(columns=["_merge", "Chrom", "Pos"])
    df_mismatch_alleles.to_csv("decode_2023_alleles_mismatch.csv", sep="\t", index=False)

# Drop allele mismatches
df_final = df_final.loc[df_final["STATUS"] != "MISMATCH"].reset_index(drop=True)
print(f"    |-> Alleles mismatch: {len(df_mismatch_alleles)}")

# Drop monomorphic alleles
nr_monomorphic_alleles = len(df_final.loc[df_final["effectAllele"] == df_final["otherAllele"]])
df_final = df_final.loc[df_final["effectAllele"] != df_final["otherAllele"]]
print(f"    |-> Monomorphic alleles: {nr_monomorphic_alleles}")

# Report allele status
print(f"\nAllele status:")
nr_ok_alleles = len(df_final.loc[df_final["STATUS"] == "OK"])
nr_miss1_alleles = len(df_final.loc[df_final["STATUS"] == "MISSING_1"])
nr_missall_alleles = len(df_final.loc[df_final["STATUS"] == "MISSING_ALL"])
print(f"    |-> OK alleles: {nr_ok_alleles}")
print(f"    |-> One allele with */!: {nr_miss1_alleles}")
print(f"    |-> Both alleles with */!: {nr_missall_alleles}")


# ---- MISSING VARIANTS ----
leads_keys = df_final[["gene_prot", "pQTL_ID_prot", "Chrom", "Pos", "rsID", "Amin", "Amaj"]].drop_duplicates()
df_missing = df_fn.merge(
    leads_keys,
    left_on=["pQTL_ID_prot", "gene_prot", "chr", "pos38", "rsID", "Amin", "Amaj"],
    right_on=["pQTL_ID_prot", "gene_prot", "Chrom", "Pos", "rsID", "Amin", "Amaj"],
    how="left",
    indicator=True
).query('_merge == "left_only"').drop(columns="_merge")

df_missing.to_csv("decode_2023_missing.csv", index=False)


# ---- MLOG10P Inf ----
mask_mlog10p_inf = np.isinf(df_final["minus_log10_pval"])
df_final.loc[mask_mlog10p_inf, "minus_log10_pval"] = df_final.loc[mask_mlog10p_inf, "mLog10pval_unadj"]
print(f"\nMLOG10P former Inf: {len(df_final.loc[mask_mlog10p_inf])}")


# ---- REPORT ----
print(f"\nTotal original variants: {len(df_fn)}")
print(f"Total found variants: {len(df_final)}")
print(f"Total missing variants: {len(df_missing)}")
print(f"    |-> Missing file paths: {len(df_missing[df_missing["match_name"] == False])}")
print(f"    |-> Missing variants in GWAS: {len(df_missing[df_missing["match_name"] != False])}")
print(f"        |-> Missing due to absence in GWAS: {len(df_missing[df_missing["match_name"] != False]) - nr_stat_mismatch - len(df_mismatch_alleles) - nr_monomorphic_alleles}")
print(f"        |-> Missing due to statistical mismatch: {nr_stat_mismatch}")
print(f"        |-> Missing due to allele mismatch: {len(df_mismatch_alleles)}")
print(f"        |-> Missing due to monomorphic status: {nr_monomorphic_alleles}")


# ---- SAVE CLEANED LITERATURE TABLE ----
df_final.to_csv("pqtl_decode_2023_leadsnps.csv", sep="\t", index=False)
