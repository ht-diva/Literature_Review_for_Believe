# BELIEVE Harmonization & Lead SNPs Export

## 0. Install dependencies

To install all required dependencies (first time), run the following command:

```
make dependencies
```

## 1. Add Studies

To add the studies:

* `pqtl_interval_chris_meta` 
* `pqtl_sun_ukb_csa`
* `pqtl_decode_2023` (fore more details on table generation, see [here](https://github.com/ht-diva/Literature_Review_for_Believe/tree/main/literature_files/pqtl_decode_2023))

run the following command:

```
make add_studies
```

The table `literature_table_all_somalogic_allstudies.xlsx` will be generated.

## 2. Clean Literature Table

Before harmonization, the literature table `literature_table_all_somalogic_allstudies.xlsx` is cleaned to `literature_table_all_somalogic_cleaned.xlsx` via:

```
make clean_table
```

Cleaning includes: 

1. removal of variants with bad alleles ("!" or "." or "*" or "NAN")
2. preliminary sanity checks for SeqIDs and UniProt IDs
3. change mismatched UniProts to UniProt format in BELIEVE Metadata
4. report literature's SEQIDs and UniProts missing in BELIEVE Metadata

## 3. Literature Table Harmonization

To harmonize the literature table, run the following command:

```
make harmonization
```

### Notes

* The harmonization is done using the [gwaspipe pipeline](https://github.com/ht-diva/gwaspipe/tree/main). 

* All harmonized tables will be generated in the `literature_harmonized` folder, along with cohort-specific log files and a summmary table `harmonization_summary.tsv`.

* Along with the table files, the command creates a release file **release.txt** with the code commit ID that generated them.

## 4. Literature Liftover

To perform bcftool liftover (pos37 to pos38), run the following command:

```
make liftover
```

The table `literature_table_all_somalogic_liftover.xlsx` will be generated.

## 5. GWASStudio Files Generation

To generate all required GWASStudio files, run the following command:

```
make gwasstudio_files
```

### Notes

The generated `literature_gwasstudio_files` folder will contain all required GWASStudio files:

* cohort-specific YAML files with search items (`search_file_*.yml`)
* cohort-specific search tables formatted as GWASStudio inputs (`*..gwaslab_formatted.csv`)
* cohort-specific sbatch scripts to execute GWASStudio exports (`run_gwasstudio_*.sbatch`)

## 6. GWASStudio Lead SNPs Export

To execute the GWASStudio export of cohort-specific lead SNPs, run the following command:

```
sbatch literature_gwasstudio_files/run_gwasstudio_*.sbatch
```

### Notes

The GWASStudio command `--get-regions-leadsnps` creates a window of given width `--region-width` (default 500,000) around trait-specific SNPs and extracts from this region the statistics MLOG10P, BETA and SE of:

* the lead SNP, i.e. the SNPID with the most significant P-value
* the exact SNP, i.e. the exact CHR:POS:EA:NEA of the input

## 7. GWASStudio Output Format

To combine and format the final GWASStudio output, run the following commmand:

```
make combine_output
```

### Notes

The final formatted outputs can be found in the `gwasstudio_output\pqtl_*` folder as `pqtl_*_hdsc_believe.csv`.
