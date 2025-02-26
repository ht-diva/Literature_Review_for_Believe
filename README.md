2025-02-26 Literature Review for Believe

**Authors:**
- Giulia Pontali, Solène Cadiou, Lucia Piubeni, Luisa Zuccolo

### Aim ###
The objective of this work is to develop a standardized literature table to evaluate whether the cis-protein quantitative trait loci (pQTL) signals identified in the INTERVAL-CHRIS meta-analysis represent novel discoveries. To achieve this, we focus on four of the largest European pQTL studies in the literature, all of which utilized the Somalogic proteomic platform -the same platform used in our analysis-.

- Sun et al., 2018 (DOI: 10.1038/s41586-018-0175-2)
- Pietzner et al., 2021 (DOI: 10.1126/science.abj1541)
- Ferkingstad et al., 2021 (DOI: 10.1038/s41588-021-00978-w)
- Zhang et al., 2022 (DOI: 10.1038/s41588-022-01051-w)

Additionally, we included the Sun et al. 2023 UK Biobank study (DOI: 10.1038/s41586-023-06592-6), which employed the Olink proteomic platform and represents the largest pQTL study conducted to date in a European population (54,219 individuals and 2,923 protein targets). Although this study used a different proteomic technology, it presents an important opportunity to compare and cross-validate pQTL signals across distinct platforms.

We also aim to apply the same approach for the Believe study by including a non-European population to ensure a broader spectrum for determining whether a signal is true or not. 

## Starting point ##
As a starting point, we have decided to include the ARIC study, which focuses on the African-American population (download file here: https://github.com/ht-diva/Literature_Review_for_Believe/tree/main/files)

**step to need to followed:**
All literature files must be standardized using the following column names:
- **pqtlID** : This column contains the following information, separated by a dash (_): rsID_SeqID_PMID_COHORT (e.g. *rs3766509_seq.5742.14_34857953_deCODE*). If some information is missing, such as rsID, leave the space blank (e.g. *_seq.3484.60_34857953_deCODE*).
- **rsID** : If it is present in the study report, include it; otherwise, leave it out. In the meta-analysis, chr:pos:a1:a2 is used, so not having an rsID is not a problem
- **chr**
_ **pos37** : This information must be present. If pos38 is available, perform a liftover. To do this, use bcftools liftover as it provides 100% coverage in both directions
- **pos38** : This information must be present. If pos38 is available, perform a liftover. To do this, use bcftools liftover as it provides 100% coverage in both directions
- **SeqID** : we need to have the following format seq.17333.20. Pay attention if we have this format 17333-20 start with adding a seq. so now we have a string and then replace - (dash) with . (dot). If we do in the otherway around we will lose the final 0 
- **OlinkID** : If this information is not available, leave it blank
- **UniProt** 
- ****
