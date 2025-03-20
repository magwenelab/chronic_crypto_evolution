import pandas as pd
import os

# Define the metadata file path
metadata_file = "chronic_metadata.csv"

# Read the metadata file
metadata = pd.read_csv(metadata_file)

# Get unique patient names
patients = metadata['strain'].unique()

# Loop over each patient
for patient in patients:
    print(f"Processing patient: {patient}")

    # Get all samples for the current patient
    patient_metadata = metadata[metadata['strain'] == patient]

    # Define the input and output file paths
    input_file = f"{patient}_combined_output.tsv"
    output_file = f"{patient}_SV_genes.csv"

    # Check if the input file exists
    if not os.path.exists(input_file):
        print(f"Input file {input_file} not found. Skipping patient {patient}.")
        continue

    # Read the combined TSV file for the patient
    df = pd.read_csv(input_file, sep='\t')

    # Rename columns to match expected format
    df.columns.values[1] = 'sample'
    df.columns.values[0] = 'doi'

    # Merge with metadata to get the plate information
    df = df.merge(patient_metadata[['sample', 'plate']], left_on='sample', right_on='sample', how='left')

    # Fill missing values with 'n/a'
    df = df.fillna('n/a')

    # If a Gene_ID column value is 'n/a', replace it with the ID value
    df['Gene_ID'] = df.apply(lambda x: x['ID'] if x['Gene_ID'] == 'n/a' else x['Gene_ID'], axis=1)

    # Drop the sample column (already merged with plate)
    df = df.drop(columns=['sample'])

    # Collapse identical rows within patients and keep the first DOI occurrence
    collapsed_df = df.groupby(df.columns.difference(['plate', 'doi']).tolist()).agg({'plate': ', '.join, 'doi': 'first'}).reset_index()

    # Remove VNI, VNII, and VNBII from accession numbers for hybrids
    collapsed_df['Chromosome'] = collapsed_df['Chromosome'].str.replace('VNI_', '', regex=False)
    collapsed_df['Chromosome'] = collapsed_df['Chromosome'].str.replace('VNII_', '', regex=False)
    collapsed_df['Chromosome'] = collapsed_df['Chromosome'].str.replace('VNBII_', '', regex=False)

    # Rename columns to be lowercase
    collapsed_df.columns = map(str.lower, collapsed_df.columns)

    # Rename column 'chromosome' to 'accession' and 'gene_id' to 'locus'
    collapsed_df = collapsed_df.rename(columns={'chromosome': 'accession', 'gene_id': 'locus'})

    # Reorder columns and sort
    collapsed_df = collapsed_df.reindex(columns=['accession', 'start', 'end', 'type', 'plate', 'doi', 'feature', 'locus', 'id', 'predicted_function'])
    sorted_df = collapsed_df.sort_values(by=['accession', 'start', 'end'])

    # Save the output to a CSV file
    os.makedirs(patient, exist_ok=True)
    sorted_df.to_csv(output_file, index=False)

    print(f"Output saved to {output_file}")