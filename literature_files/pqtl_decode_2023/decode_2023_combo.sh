head -n 1 decode_2023_leadsnp_sumstats/Proteomics_SMP_PC0_9999_1_IRF6_IRF6_10032022_leadsnps.csv > pqtl_decode_2023_all_chr_pos.csv
for f in decode_2023_leadsnp_sumstats/Proteomics_SMP_PC0_*_leadsnps.csv; do
    tail -n +2 "$f" >> pqtl_decode_2023_all_chr_pos.csv
done

