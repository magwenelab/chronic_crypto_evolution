# chronic_crypto_evolution
Code to accompany Montoya et al. study of fungal genome evolution during chronic cryptococcus infection.


## FungalPop_parameters
This folder contains the config parameters and metadata files used to run the FungalPop pipeline.

### Additional files
The FungalPop pipeline was run with the following repeat database file (not included in this repository):
- RepBase library fasta - version 29.01

### Reference genomes
The FungalPop pipeline was run with the following reference genomes (not included in this repository):
* VNI - CP003820.1, version=2015-04-01 
- Genome assembly = FungiDB release 65 Cryptococcus neoformans var. grubii H99 
    - GenBank = GCF_000149245.1 
- Genome annotation = FungiDB release 65 Cryptococcus neoformans var. grubii H99 
    - GenBank = GCA_000149245.3

All other genomes were annotated via liftover from the VNI genome annotation.

* VNII - CP091247.1, version=2022-04-07
- Genome assembly = FungiDB release 65 Cryptococcus neoformans strain:VNII
    - GenBank = GCA_022832995.1

* VNBII - CP097910.1, version=2022-04-07
- Genome assembly = Cryptococcus neoformans Bt65
    - GenBank = GCA_023650535.1

* VNIV - AE017341.1 , version=2016-06-16 
- Genome assembly = FungiDB release 65 Cryptococcus neoformans var. neoformans JEC21
    - GenBank = GCF_000091045.1


## FungalPop_output
This folder contains subfolders for each patient with read depth plots produced by the FungalPop pipeline. Each plot file is named with the DOI when the associated sample was collected. 

This folder also contains output files from post-processing the FungalPop output using the scripts in the FungalPop_scripts folder.


## FungalPop_scripts
This folder contains the scripts used to process the files created by the FungalPop pipeline. 

### Extracting genes within structural variants
* sv_genes_intersect.sh - obtains gene IDs located within structural variants by intersecting the called deletion and duplication regions with the gene annotation file.
* collapse_intersect_files.py - collapses the files created by sv_genes_intersect.sh into a readable single file for each patient.
    * -> {patient_ID}_SV_genes.csv

### Extracting variant calls and predicted effects
* query_database_variants.py - queries the DuckDB database of called SNPs and predicted effects produced by FungalPop to obtain the high-impact and moderate-impact variants for each patient.
* query_database_commands.py - contains the commands used to query the DuckDB database.
    * -> chronic_high_impact_variants.csv
    * -> chronic_moderate_impact_variants.csv

## manuscript_figures
This folder contains the R scripts and other code used to generate the figures in the manuscript.
