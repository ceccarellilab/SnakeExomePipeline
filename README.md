# ExomePipeline

### Input files:
- options.yaml (yaml file): 
	- **batch** : name of batch a tab-delimited csv file

	- **fastq_1_suffix** : fastq file suffix
	- **fastq_2_suffix** : fastq file suffix (If using Illumina paired data)
	- **platform** : name of the sequencing platform used to sequence the DNA. Possible options:
        - ILLUMINA when the fastq files have been produced in an Illumina machine;
        - IONTORRENT when the fastq files have been produced in a Life Technologies Ion-Torrent machine.

	- **fasta** : path fasta file 
	- **dbsnp** : path dbsnp file  
	- **known_Mills_indels** : path known_Mills_indels file  
	- **known_1000G_indels** : path known_1000G_indels file  

	- **panel_of_normal_TNsnv** : path panel of normal (recommended to create the normal file panel with the corresponding algorithm that you plan to use for somatic mutation calling)
	- **panel_of_normal_TNhaplotyper** : path panel of normal
	- **cosmic_db** : path CosmicDB vcf files 

	- **ML_MODEL_N** : Path DNAscope model 
	- **ML_MODEL_T** : Path TNscope model

	- **SENTIEON_INSTALL_DIR** : Path Sentieon software 
	- **SENTIEON_LICENSE** : Sentieon license file/server

	- **bcfdir** : Path bcftools

	- **bgzipdir** : Path bgzipdir

	- **samplesParallel** : Numbers of samples in parallel (set 1000 for all in parallel)
	- **nt** : Number of threads to use for each command
 
- BATCH.tsv (tab-delimited csv file): 
	- **COLUMN1** : name of sample
	- **COLUMN2** : name of tumor
	- **COLUMN3** : folder path to fastq of tumor
	- **COLUMN4** : name of normal
	- **COLUMN5** : folder path to fastq of normal

### Output file:
- **html_pipeline.zip**: Contains a dynamic HTML report executable locally