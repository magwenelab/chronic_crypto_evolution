import pandas as pd
import WeavePop_scripts.query_database_commands as qd
import os

# Define the metadata file path and database
metadata_file = "~/WeavePop_Chronic_Crypto/config/chronic_metadata.csv"

# Database file produced during the WeavePop analysis step on all chronic samples
database = "~/WeavePop_Chronic_Crypto/results_chronic_1/02.Dataset/database.db"

# Import chronic sample metadata
chronic_metadata_df = pd.read_csv(metadata_file)

# Define output directory
output_dir = "database_calling/"
os.makedirs(output_dir, exist_ok=True)

# Process for both HIGH and MODERATE impact values
for impact in ["HIGH", "MODERATE"]:
    print(f"Processing variants with impact: {impact}")

    # Loop through each sample and extract the variants
    for sample in chronic_metadata_df['sample']:
        print(f"Processing sample: {sample}")
        df = qd.effects_chrom_by_sample_location(db=database, sample=sample, impact=impact)

        # Save the extracted variants to a CSV file
        output_file = os.path.join(output_dir, f"{sample}_{impact}_variants.csv")
        df.to_csv(output_file, index=False)
        print(f"Saved {impact} variants for sample {sample} to {output_file}")

    # Concatenate all variant files for the current impact level
    concatenated_file = os.path.join(output_dir, f"concatenated_{impact}_variants.csv")
    header_added = False

    with open(concatenated_file, "w") as outfile:
        for sample in chronic_metadata_df['sample']:
            input_file = os.path.join(output_dir, f"{sample}_{impact}_variants.csv")
            if not os.path.exists(input_file):
                continue

            with open(input_file, "r") as infile:
                lines = infile.readlines()
                if not header_added:
                    # Write the header from the first file
                    outfile.write(lines[0].strip() + ",sample\n")
                    header_added = True
                # Write the data rows with the sample name appended
                for line in lines[1:]:
                    outfile.write(line.strip() + f",{sample}\n")

    print(f"Concatenated file for {impact} variants saved to {concatenated_file}")

    # Remove individual variant files
    for sample in chronic_metadata_df['sample']:
        input_file = os.path.join(output_dir, f"{sample}_{impact}_variants.csv")
        if os.path.exists(input_file):
            os.remove(input_file)

    # Read the concatenated file and merge with metadata
    concat_df = pd.read_csv(concatenated_file)
    merged_df = concat_df.merge(chronic_metadata_df, on="sample", how="left")
    merged_df = merged_df.rename(columns={"plate": "strain_id"})

    # Fill missing values with 'n/a'
    merged_df = merged_df.fillna("n/a")

    # Drop the sample column
    merged_df = merged_df.drop(columns=["sample"])

    # Collapse identical rows within patients and keep the first DOI occurrence
    collapsed_df = merged_df.groupby(
        merged_df.columns.difference(["strain_id", "doi"]).tolist()
    ).agg({"strain_id": ", ".join, "doi": "min"}).reset_index()


    # Remove VNI and VNBII from accession numbers for hybrids
   # collapsed_df["accession"] = collapsed_df["accession"].str.replace("VNI_", "")
    #collapsed_df["accession"] = collapsed_df["accession"].str.replace("VNBII_", "")

    # Reorder columns and sort
    collapsed_df = collapsed_df.reindex(
        columns=[
            "lineage",
            "strain",
            "doi",
            "strain_id",
            "accession",
            "gene_id",
            "gene_coding",
            "effect_type",
            "effect",
            "transcript_biotype",
            "pos",
            "ref",
            "alt",
            "amino_acid_change",
            "var_id",
        ]
    )
    sorted_df = collapsed_df.sort_values(by=["lineage", "strain", "accession", "gene_id"])

    # Save the final output
    final_output_file = os.path.join(output_dir, f"all_cnags_{impact.lower()}_variants.csv")
    sorted_df.to_csv(final_output_file, index=False)
    print(f"Final output for {impact} variants saved to {final_output_file}")