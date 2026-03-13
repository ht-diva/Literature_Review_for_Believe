awk -F'\t' 'NR==1 {next} $4 != "False" && $3 != "" {print $3}' decode_2023_filenames.csv | sort -u > decode_2023_filenames_list.txt
