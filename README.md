# ExomePipeline

### Input files:
- ***options.yaml*** (yaml file): 
	- **batch** : name of batch a tab-delimited csv file

	- **fastq_1_suffix** : FASTQ file suffix
	- **fastq_2_suffix** : FASTQ file suffix (If using Illumina paired data)
	- **platform** : name of the sequencing platform used to sequence the DNA. Possible options:
        - ILLUMINA when the fastq files have been produced in an Illumina machine;
        - IONTORRENT when the fastq files have been produced in a Life Technologies Ion-Torrent machine.

	- **fasta** : Path to FASTA file containing the nucleotide sequence of the reference genome corresponding to the sample you will analyze  
	- **dbsnp** : Path to dbSNP (Single Nucleotide Polymorphism database) that will be used to label known variants (Only one file is supported)  
	- **known_Mills_indels** : Path to Mills indels (INsertion-DELetions) file.
	- **known_1000G_indels** : Path to 1000G indels (INsertion-DELetions) file.

	- **panel_of_normal_TNsnv** : Path to panel of normal (recommended to create the normal file panel with the corresponding algorithm that you plan to use for somatic mutation calling)
	- **panel_of_normal_TNhaplotyper** : path panel of normal
	- **cosmic_db** : path to Catalogue of Somatic Mutations in Cancer (COSMIC) VCF file used to create the panel of normal file (Only one file is supported)  

	- **ML_MODEL_N** : Path to DNAscope machine learning model (second step of variant calling)
	- **ML_MODEL_T** : Path to TNscope model (variant filtration)

	- **SENTIEON_INSTALL_DIR** : Path Sentieon software 
	- **SENTIEON_LICENSE** : Sentieon license file/server

	- **bcfdir** : Path to bcftools a set of utilities that manipulate variant calls in the Variant Call Format (VCF) (https://samtools.github.io/bcftools/bcftools)

	- **bgzipdir** : Path to bgzip (VCF files compressed with bgzip)

	- **samplesParallel** : Numbers of samples in parallel (set 1000 for all in parallel)
	- **nt** : Number of threads to use for each command
 
- ***BATCH.tsv*** (tab-delimited csv file): 
	- **SAMPLE** : Name of Sample
	- **TUMOR** : Name of Tumor
	- **FASTQ_TUMOR_PATH** : Folder path to FASTQ of Tumor
	- **NORMAL** : Name of Normal
	- **FASTQ_NORMAL_PATH** : Folder path to FASTQ of Normal

### Output file:
- **html_pipeline.zip**: Contains a dynamic HTML report executable locally