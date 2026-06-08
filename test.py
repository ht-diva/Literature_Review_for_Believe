import pandas as pd

harm_ref_check = pd.read_csv("harmonization_summary.tsv", sep="\t")
bcftools_ref_check = pd.read_csv("summary.txt", sep="\t")

col_filt = ["COHORT", "REF_MATCH", "REF_FLIP_VARIANT_NR", "REF_PALINDROMIC_NR"]
harm_ref_check = harm_ref_check[col_filt].copy()
harm_ref_check["COHORT"] = harm_ref_check["COHORT"].str.replace("pqtl_", "", regex=True)
harm_ref_check = harm_ref_check.rename(columns={
    "COHORT" : "study",
    "REF_MATCH" : "gwaslab_ref_match",
    "REF_FLIP_VARIANT_NR" : "gwaslab_variants_flip",
    "REF_PALINDROMIC_NR": "gwaslab_palindromic",
    })
print(harm_ref_check)

bcftools_ref_check = bcftools_ref_check.rename(columns={
    "ref_match": "bcftools_ref_match",
    "ref_mismatch": "bcftools_ref_mismatch",
    })

merged = harm_ref_check.merge(
    bcftools_ref_check,
    on = "study",
    how = "inner"
    )
print(merged)

merged.to_csv("ref_match_summary.tsv", sep="\t", index=False)
