#!/bin/bash

# Define the metadata file path
metadata_file="chronic_metadata.csv"

# Read the metadata file and process each unique patient and lineage
awk -F',' 'NR > 1 { print $2, $3 }' "$metadata_file" | sort -u | while read -r lineage patient; do
    echo "Processing patient: $patient with lineage: $lineage"

    # Define the GFF file for the lineage
    gff_file="/data/references/${lineage}.gff"
    protein_coding_genes_file_suffix="_protein_coding_genes.tsv"

    # Create a combined TSV file for the patient
    combined_file="${patient}_combined_output.tsv"

    # Check if the combined file already exists; if so, remove it
    if [ -f "$combined_file" ]; then
        rm "$combined_file"
    fi

    # Get all samples associated with the patient
    awk -F',' -v pat="$patient" '$3 == pat { print $3 }' "$metadata_file" | while read -r name; do
        input_file="${name}_cnv.tsv"
        tsv_file="filtered_${name}_cnv.tsv"

        # Read in the cnv_calls.tsv file from FungalPop output
        cp ../${name}/cnv_calls.tsv $input_file

        # Use awk to filter the rows where the value in the 10th column (repeat fraction) is less than 0.5
        # This removes duplications and deletions in highly repetitive regions
        awk -F'\t' '$10 < 0.5' "$input_file" > "$tsv_file"

        echo "Filtered data has been written to $tsv_file"

        bed_file="${name}_cnags_features.bed"
        intersected_file="${name}_cnags_intersected.tsv"
        output_file="${name}_cnags_output.tsv"
        filtered_output_file="${name}_cnags_filtered_output.tsv"
        gene_ids_file="${name}_cnags_unique_gene_ids.txt"
        protein_coding_genes_file="${name}${protein_coding_genes_file_suffix}"

        # Step 1: Extract relevant columns from TSV and create a BED file
        awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$7}' "$tsv_file" > "$bed_file"

        # Step 2: Intersect BED with GFF using bedtools
        bedtools intersect -a "$bed_file" -b "$gff_file" -wa -wb > "$intersected_file"

        # Step 3: Extract relevant information from intersected results
        awk -F'\t' '
        BEGIN { OFS="\t"; print "Chromosome", "Start", "End", "Type", "Feature", "ID", "Gene_ID", "Predicted_Function" }
        {
            split($13, gff_info, ";");
            gene_name="";
            predicted_function="";
            gene_id="";
            for(i=1; i<=length(gff_info); i++) {
                if (gff_info[i] ~ /^ID=/) {
                    gene_name = substr(gff_info[i], 4);
                }
                if (gff_info[i] ~ /^description=/) {
                    predicted_function = substr(gff_info[i], 13);
                }
                if (gff_info[i] ~ /^(gene_id=|locus=|Parent=)/) {
                    gene_id = substr(gff_info[i], index(gff_info[i], "=") + 1);
                }
            }
            print $1, $2, $3, $4, $7, gene_name, gene_id, predicted_function
        }' "$intersected_file" > "$output_file"

        echo "Output written to $output_file"

        # Step 4: Filter out exons if their Gene_ID matches an ID of a protein_coding_gene
            # rename liftOvered genes back to the CNAG nomenclature
        awk -F'\t' '
        BEGIN { OFS="\t" }
        {
            if ($5 == "protein_coding_gene" || $5 == "gene" ||$5 == "mRNA"||$5 == "ncRNA_gene") {
                ids[$6] = 1;
            } else if (($5 == "exon" || $5 == "CDS" || $5 == "five_prime_UTR" || $5 == "three_prime_UTR") && ids[$7] == 1) {
                next;  # Skip exon or CDS if its Gene_ID matches a protein_coding_gene ID
            } 
            print $0;  # Print all other lines
        }
        ' "$output_file" | awk -F'\t' '!($5 == "mRNA")' | awk -F'\t' '{ gsub(/(CBAG|CTAG|CDAG)/ , "CNAG", $0); print }' > "$filtered_output_file"

        echo "Filtered output written to $filtered_output_file"

        cp $filtered_output_file /results/patient_plots_by_doi/${patient}/

        # Step 4: Extract unique Gene_IDs from the filtered output
        awk -F'\t' '$7 != "" && !seen[$7]++ { print $7 }' "$filtered_output_file" > "$gene_ids_file"

        echo "Unique Gene IDs extracted to $gene_ids_file"

        # Add DOI and sample columns, then append to the patient's combined file
        DOI=$(awk -F',' -v pat="$patient" -v sample="$name" '$3 == sample { print $9 }' "$metadata_file")
        awk -v doi="$DOI" -v sample="$name" 'BEGIN {FS="\t"; OFS="\t"} {print doi, sample, $0}' "$filtered_output_file" >> "$combined_file"
    done

    echo "Combined file for patient $patient written to $combined_file"
done