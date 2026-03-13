2025-02-26 Literature Review for Believe

**Authors:**
- Giulia Pontali, Lucia Piubeni, Solène Cadiou, Claudia Giambartolomei

## Literature Review Aim ##
The objective of this work is to develop a standardized literature table to evaluate whether the cis-protein quantitative trait loci (pQTL) signals identified in the INTERVAL-CHRIS meta-analysis represent novel discoveries. To achieve this, we focus on four of the largest European pQTL studies in the literature, all of which utilized the Somalogic proteomic platform -the same platform used in our analysis-.

- Sun et al., 2018 (DOI: 10.1038/s41586-018-0175-2)
- Pietzner et al., 2021 (DOI: 10.1126/science.abj1541)
- Ferkingstad et al., 2021 (DOI: 10.1038/s41588-021-00978-w)
- Zhang et al., 2022 (DOI: 10.1038/s41588-022-01051-w)

Additionally, we included the Sun et al. 2023 UK Biobank study (DOI: 10.1038/s41586-023-06592-6), which employed the Olink proteomic platform and represents the largest pQTL study conducted to date in a European population (54,219 individuals and 2,923 protein targets). Although this study used a different proteomic technology, it presents an important opportunity to compare and cross-validate pQTL signals across distinct platforms. </br>
(/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/literature_table/literature_file.xlsx)

We also aim to apply the same approach for the Believe study by including a non-European population to ensure a broader spectrum for determining whether a signal is new or not. 

## Starting point ##
We still need to discuss with Adam on the 19th of March to finalize which other studies to include in our analysis. However, in the meantime, the work can begin with the ARIC study as a non-European dataset to start familiarizing with the process (you can download the ARIC table here: https://github.com/ht-diva/Literature_Review_for_Believe/tree/main/literature_files).
</br>
</br>
All literature files must be standardized using the following column names:
- **pqtlID** : This column contains the following information rsID_SeqID_PMID_COHORT separated by a dash (e.g. *rs3766509_seq.5742.14_34857953_deCODE*). If any information is missing, such as rsID, leave the space blank (e.g. *_seq.3484.60_34857953_deCODE*).
- **rsID** : If the information is present in the study, include it; otherwise, leave it out. In the meta-analysis, the format chr:pos:a1:a2 is used, so the absence of an rsID is not an issue
- **chr** 
- **pos37** : This information must be present. If pos38 is available, perform a liftover. To do this, use bcftools liftover as it provides 100% coverage in both directions
```
step1: Convert the file into a VCF file (The function needs to be adapted)
write_vcf <- function(df, output_filename) {
  # Create a connection to write to a file
  vcf_file <- file(output_filename, "w")

  # Use tryCatch to ensure the file is closed properly
  tryCatch({
    # Write the VCF header
    writeLines("##fileformat=VCFv4.2", vcf_file)
    writeLines("##source=RScript", vcf_file)
    writeLines("##reference=GRCh38", vcf_file)
    writeLines("##INFO=<ID=EAF,Number=1,Type=Float,Description=Effect Allele Frequency>", vcf_file)
    writeLines("##INFO=<ID=BETA,Number=1,Type=Float,Description=Effect Size Estimate>", vcf_file)
    writeLines("##INFO=<ID=SE,Number=1,Type=Float,Description=Standard Error>", vcf_file)
    writeLines("##INFO=<ID=N,Number=1,Type=Integer,Description=Sample Size>", vcf_file)
    writeLines("##INFO=<ID=MLOG10P,Number=1,Type=Float,Description=Negative Log10 P-value>", vcf_file)
    writeLines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", vcf_file)

    # Write each row as a VCF entry
    for (i in 1:nrow(df)) {
      # Extract chromosome, position, reference allele (NEA), and alternate allele (EA)
      chrom <- df$CHR[i]
      pos <- df$POS[i]
      id <- df$SNPID[i]
      ref <- df$EA[i]
      alt <- df$NEA[i]
      qual <- "."
      filter <- "."

      # Create the INFO field (customize as needed)
      info <- paste0("EAF=", df$EAF[i], ";BETA=", df$BETA[i], ";SE=", df$SE[i],
                     ";N=", df$N[i], ";MLOG10P=", df$MLOG10P[i])

      # Write the line to the VCF file
      line <- paste(chrom, pos, id, ref, alt, qual, filter, info, sep = "\t")
      writeLines(line, vcf_file)
    }
  }, finally = {
    # Ensure the file connection is closed
    close(vcf_file)
  })
}

write_vcf(df, vcf_path)

step2: 
export BCFTOOLS_PLUGINS=/group/diangelantonio/software/liftOver_plugins/score_1.20-20240505 && \
        bgzip -c {input_vcf} > {input_vcf}.gz && \
        tabix -p vcf {input_vcf}.gz && \
        bcftools norm -f {input.hg37} -c s -Oz -o {output.output_norm} {input_vcf}.gz && \
        bcftools +liftover --no-version -Ou {output.output_norm} -- -s {input.hg37} -f {input.hg38} -c {input.chain_file} > {output.output_vcf} && \
        bcftools view {output.output_vcf} > {output.output_txt}
```
- **pos38** : This information must be present. If pos38 is available, perform a liftover. To do this, use bcftools liftover as it provides 100% coverage in both directions
```
step1: convert the file into a vcf-file
step2: use BCFTOOLS
```
- **SeqID** : The following format should be used: seq.17333.20. If the format is 17333-20, start by adding 'seq.' to the string, then replace the dash (-) with a dot (.). If the operation is reversed, the final 0 will be lost
- **OlinkID** : If this information is not available, leave it blank
- **UniProt** 
- **OTHER_ALLELE**
- **EFFECT_ALLELE**
- **cis_trans**
- **PMID**
- **BETA**
- **SE**
- **minuslog10pval** : Please use the following code to retrieve it
```
library(rtracklayer)
pval <- 2 * pnorm(mpfr(-abs(data$BETA / data$SE), 120))
data$minuslog10pval <- as.numeric(-log10(pval))
```
- **SAMPLE_SIZE**
- **COHORT**
- **TECHNOLOGY**
- **Unit**


## Literature Harmonization & Lead SNPs Export ##

### Aim

Literature harmonization is designed to format summary statistics based on GWASLab through the pipeline [gwaspipe](https://github.com/ht-diva/gwaspipe/tree/main). This allows to compare literature tables with the harmonized BELIEVE summary statistics hosted on tileDB (GWASStudio).

### Harmonization steps

**Harmonization** includes the following **steps**:

- Check SNP identifiers (SNPID/rsID).
- Fix chromosome notation (CHR), basepair positions (POS) and alleles (EA and NEA).
- Sanity check on statistics.
- Infer genome reference build version.
- Align alleles sorting them alphabetically. Precisely EA is the effect allele for which BETA is estimated, and EA is alphabetically lower than NEA (in case, flip the alleles to match the alphabetical order).
- Flip allele-specific statistics for mismatches: BETA = - BETA; Z = - Z; EAF = 1 - EAF.
- Build SNPID column (CHR:POS:EA:NEA).
- Re-name and re-order columns based on GWASLab format.

### Lead SNPs Export

[GWASStudio](https://github.com/ht-diva/gwasstudio) was used to conduct locus-level replication analyses for each published lead SNP reported in the harmonized literature tables. The GWASStudio command `--get-regions-leadsnps` defined a genomic windows of given width (cis = gene ±500 kb; trans, or if not specified in the source study = lead SNP ±1 Mb) around each published lead SNP and extracted from this window the BELIEVE GWAS summary statistics (MLOG10P, BETA and SE) of:

* the exact SNP, i.e. the exact CHR:POS:EA:NEA of the input
* the lead SNP, i.e. the SNPID with the most significant P-value

### How to run

To run all the harmonization and GWASStudio analyses, please see https://github.com/ht-diva/Literature_Review_for_Believe/tree/main/literature_scripts/README.md
