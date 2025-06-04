# chronic_crypto_evolution
Code to accompany Montoya et al. study of fungal genome evolution during chronic cryptococcus infection.

This study uses WeavePop (https://github.com/magwenelab/WeavePop) to analyze the genomes of the sequenced isolates.

## WeavePop_parameters
This folder contains the config parameters and metadata files used to run the WeavePop pipeline.

### Additional files
The WeavePop pipeline was run with the following repeat database file (not included in this repository):
- RepBase library fasta file - version 29.01

### Reference genomes
The WeavePop pipeline was run with the following reference genomes (not included in this repository):
* VNI - CP003820.1, version=2015-04-01 
    - Genome assembly = FungiDB release 65 Cryptococcus neoformans var. grubii H99 
        - GenBank = GCF_000149245.1 
    - Genome annotation = FungiDB release 65 Cryptococcus neoformans var. grubii H99 
        - GenBank = GCA_000149245.3

All other reference genomes were aligned to the VNI reference to normalize the chromosome numbering and orientation between references. The mapping from the original accessions to the VNI-aligned chromosome numbers and whether the chromosomes were reverse complemented can be found in the genome_accession_mapping.csv file in WeavePop_parameters. The renamed and reoriented reference genomes were annotated via liftover from the VNI genome annotation through WeavePop.

* VNII - CP091247.1, version=2022-04-07
    - Genome assembly = FungiDB release 65 Cryptococcus neoformans strain:VNII
        - GenBank = GCA_022832995.1

* VNBII - CP097910.1, version=2022-04-07
    - Genome assembly = Cryptococcus neoformans Bt65
        - GenBank = GCA_023650535.1

* VNIV - AE017341.1 , version=2016-06-16 
    - Genome assembly = FungiDB release 65 Cryptococcus neoformans var. neoformans JEC21
        - GenBank = GCF_000091045.1

## WeavePop_scripts
This folder contains the scripts used to process the files created by the WeavePop pipeline. 

### Extracting genes within structural variants
* patient_sv_pipeline.py - Extracts the WeavePop output file containing gene ID's that fall within duplicated or deleted regions, and collapses them across samples for each patient.
    * -> {patient_ID}_collapsed_sv_genes.tsv

### Extracting variant calls and predicted effects
* query_database_script.py - queries the DuckDB database of called SNPs and predicted effects produced by WeavePop to obtain the high-impact and moderate-impact variants for each patient.
* query_database_commands.py - contains the commands used to query the DuckDB database.
    * -> chronic_high_impact_variants.csv
    * -> chronic_moderate_impact_variants.csv

## WeavePop_output
This folder contains subfolders for each patient with read depth plots produced by the WeavePop pipeline. Each plot file is named with the DOI when the associated sample was collected. 

This folder also contains output files from post-processing the WeavePop output using the scripts in the WeavePop_scripts folder.

