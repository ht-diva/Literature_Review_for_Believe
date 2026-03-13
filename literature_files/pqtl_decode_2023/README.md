# Make Literature Table for deCODE 2023 

We extracted the variants listed in the Supplementary Table 19 from Eldjarn et al. (2023) (ST19) to:
- correctly classify alleles as effect and non-effect alleles (ST19 defines alleles as major and minor without specifying the effect allele);
- obtain full statistics with higher decimal accuracy (SE and ImpMAF are missing in ST19).

Additionally, we performed a quality check to exclude monomorphic alleles and verify the correct statistical and allele match after SNP extraction.

## Step 1: Extract matching file paths of summary statistics

**Step 1.A.** Build file paths (pQTL_ID_prot + gene_prot) and extract supporting information from ST19
- Input: decode_2023_st19.csv
- Output: decode_2023_filenames.csv

```
python decode_2023_extract_file_paths.py 
```

**Step 1.B.** Extract only file paths
- Output: decode_2023_filenames_list.txt 

```
./decode_2023_get_file_paths.sh
```

## Step 2: Extract lead SNPs from summary statistics

**Step 2.A.** Extract matching ID:CHR:POS of lead SNPs listed in ST19 from summary statistics
- Output: All filtered summary statistics in decode_2023_leadsnp_sumstats

```
sbatch decode_2023_extract_alleles_stats.sbatch
```

**Step 2.B.** Combine all outputs in one table
- Output: decode_2023_leadsnp_sumstats/pqtl_decode_2023_all_chr_pos.csv 

```
./decode_2023_combo.sh
```

## Step 3: Quality check & Missing variants

Final list of extracted variants (pQTIL_ID-Chr-Pos) matched for the closest statistics
- Output: pqtl_decode_2023_leadsnps.csv

Lists of statstical and allele mismatches:
- Output: decode_2023_stats_mismatch.csv
- Output: decode_2023_alleles_mismatch.csv

List of missing variants with ST19 information:
- Output: decode_2023_missing.csv

```
python decode_2023_qc.py
```

## Final QC Report

Starting from the 39,520 variants of ST19, we found and retained 34,829 for the literature review.

```
Dropped variants:
    |-> Statistical mismatch: 15
    |-> Alleles mismatch: 153
    |-> Monomorphic alleles: 4333

Allele status:
    |-> OK alleles: 29998
    |-> One allele with */!: 4831
    |-> Both alleles with */!: 0

Total original variants: 39520
Total found variants: 34829
Total missing variants: 4691
    |-> Missing file paths: 110
    |-> Missing variants in GWAS: 4581
        |-> Missing due to absence in GWAS: 80
        |-> Missing due to statistical mismatch: 15
        |-> Missing due to allele mismatch: 153
        |-> Missing due to monomorphic status: 4333
```

## References

decode_2023_st19.csv corresponds to the Supplemetary Table 19 from:

```
Eldjarn, G.H., Ferkingstad, E., Lund, S.H. et al. Large-scale plasma proteomics comparisons through genetics and disease associations. Nature 622, 348–358 (2023). https://doi.org/10.1038/s41586-023-06563-x
```