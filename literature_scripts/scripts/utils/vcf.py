import pandas as pd


# ---- CONVERT TO VCF FUNCTION ----


def write_vcf(df, output_filename, build="GRCh37"):
    with open(output_filename, "w") as vcf_file:

        # Ensure POS is valid integer
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
        df = df.dropna(subset=["POS"])
        df["POS"] = df["POS"].astype(int)

        # Convert CHR for bcftools processing
        df["CHR"] = df["CHR"].astype(str).replace({
            "23": "X",
            "24": "Y",
            "25": "MT"
        })
        if build == "GRCh38":
            df["CHR"] = df["CHR"].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)

        # Write the VCF header
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write("##source=PythonScript\n")
        vcf_file.write(f"##reference={build}\n")
        vcf_file.write('##INFO=<ID=BETA,Number=1,Type=Float,Description=Effect Size Estimate>\n')
        vcf_file.write('##INFO=<ID=SE,Number=1,Type=Float,Description=Standard Error>\n')
        vcf_file.write('##INFO=<ID=N,Number=1,Type=Integer,Description=Sample Size>\n')
        vcf_file.write('##INFO=<ID=MLOG10P,Number=1,Type=Float,Description=Negative Log10 P-value>\n')
        chroms = df["CHR"].unique()
        for c in chroms:
            vcf_file.write(f"##contig=<ID={c}>\n")
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Iterate through rows
        for _, row in df.iterrows():
            chrom = row["CHR"]
            pos = row["POS"]
            vid = row["rsID"]
            ref = row["EA"]
            alt = row["NEA"]
            qual = "."
            filt = "."

            info = (
                f"BETA={row['BETA']};"
                f"SE={row['SE']};"
                f"N={row['N']};"
                f"MLOG10P={row['MLOG10P']}"
            )

            line = f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\n"
            vcf_file.write(line)
