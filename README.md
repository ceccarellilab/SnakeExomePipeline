# ExomePipeline

### Input files:
- ***options.yaml*** (yaml file): 

    Set the batch
	- **batch** : Name of batch a tab-delimited csv file

    Samples
	- **fastq_1_suffix** : FASTQ file suffix
	- **fastq_2_suffix** : FASTQ file suffix (If using Illumina paired data)
	- **platform** : Name of the sequencing platform used to sequence the DNA. Possible options:
        - ILLUMINA when the fastq files have been produced in an Illumina machine;
        - IONTORRENT when the fastq files have been produced in a Life Technologies Ion-Torrent machine.

    Reference data files
	- **fasta** : Path to FASTA file containing the nucleotide sequence of the reference genome corresponding to the sample you will analyze  
	- **dbsnp** : Path to dbSNP (Single Nucleotide Polymorphism database) that will be used to label known variants (Only one file is supported)  
	- **known_Mills_indels** : Path to Mills indels (INsertion-DELetions) file, a set of known sites used to help identify likely sites where the realignment is necessary; only indel variants in the file will be used.  (can include multiple collections of known sites)
	- **known_1000G_indels** : Path to 1000G indels file.

    Panel of normal and CosmicDB vcf files
	- **panel_of_normal_TNsnv** : Path to the Panel of Normal (PON) for the TNsnv algorithm. File containing the variants detected in the Panel of Normal analysis that will be used to remove false positives (Only one file is supported and recommended to create the normal file panel with the corresponding algorithm that you plan to use for somatic mutation calling). TNsnv algorithm performs the somatic variant calling on the tumor-normal matched pair or the tumor and panel of normal data, using a Genotyper algorithm.
	- **panel_of_normal_TNhaplotyper** : Path to the Panel of Normal (PON) for the TNhaplotyper algorithm. The TNhaplotyper algorithm performs the somatic variant calling on the tumor-normal matched pair or the tumor and panel of normal data, using a Haplotyper algorithm.
    
	- **cosmic_db** : Path to Catalogue of Somatic Mutations in Cancer (COSMIC) VCF file used to create the panel of normal file (Only one file is supported)  

    DNAscope and TNscope models
	- **ML_MODEL_N** : Path to DNAscope machine learning model (second step of variant calling)
	- **ML_MODEL_T** : Path to TNscope model (variant filtration)

    Software
	- **SENTIEON_INSTALL_DIR** : Path Sentieon software 
	- **SENTIEON_LICENSE** : Sentieon license file/server
	- **bcfdir** : Path to bcftools a set of utilities that manipulate variant calls in the Variant Call Format (VCF) (https://samtools.github.io/bcftools/bcftools)
	- **bgzipdir** : Path to bgzip (VCF files compressed with bgzip)

    Parallelism settings
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