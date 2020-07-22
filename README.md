# ExomePipeline

- 1a Mapping reads with BWA-MEM, sorting for normal and tumor samples
- 2a Metrics for normal and tumor samples
- 3a Remove Duplicate Reads for normal and tumor samples
- 4a Indel realigner for normal and tumor samples
- 5a Base recalibration for normal and tumor samples
- 7a HC Variant caller normal and tumor samples
- 8a Variant calling DNAscope normal and tumor samples
- 9a Variant Annotation normal and tumor samples
<!-- -->
- 7b Somatic Variant Calling TNseq
- 8b Somatic Variant calling TNscope
- 9a Somatic Variant Annotation

### Input files:
- ***options.yaml*** (yaml file): 

    (Set the batch)
	- **BATCH** : Name of batch a tab-delimited csv file (without extension)

    (Samples)
	- **FASTQ_1_SUFFIX** : FASTQ file suffix
	- **FASTQ_2_SUFFIX** : FASTQ file suffix (If using Illumina paired data)
	- **PLATFORM** : Name of the sequencing platform used to sequence the DNA. Possible options:
        - ILLUMINA when the fastq files have been produced in an Illumina machine;
        - IONTORRENT when the fastq files have been produced in a Life Technologies Ion-Torrent machine.

    (Reference data files)
	- **FASTA** : Path to FASTA file containing the nucleotide sequence of the reference genome corresponding to the sample you will analyze  
	- **DBSNP** : Path to dbSNP (Single Nucleotide Polymorphism database) that will be used to label known variants (Only one file is supported)  
	- **KNOWN_MILLS_INDELS** : Path to Mills indels (INsertion-DELetions) file, a set of known sites used to help identify likely sites where the realignment is necessary; only indel variants in the file will be used.  (can include multiple collections of known sites)
	- **KNOWN_1000G_INDELS** : Path to 1000G indels file.

    (Panel of normal and CosmicDB vcf files)
	- **PON_TNsnv** : Path to the Panel of Normal (PON) for the TNsnv algorithm. File containing the variants detected in the Panel of Normal analysis that will be used to remove false positives (Only one file is supported and recommended to create the normal file panel with the corresponding algorithm that you plan to use for somatic mutation calling). TNsnv algorithm performs the somatic variant calling on the tumor-normal matched pair or the tumor and panel of normal data, using a Genotyper algorithm.
	- **PON_TNhaplotyper** : Path to the Panel of Normal (PON) for the TNhaplotyper algorithm. The TNhaplotyper algorithm performs the somatic variant calling on the tumor-normal matched pair or the tumor and panel of normal data, using a Haplotyper algorithm.
    
	- **COSMIC_DB** : Path to Catalogue of Somatic Mutations in Cancer (COSMIC) VCF file used to create the panel of normal file (Only one file is supported)  

    (DNAscope and TNscope models)
	- **ML_MODEL_N** : Path to DNAscope ML model to apply on the results of DNAscope
	- **ML_MODEL_T** : Path to ML model to apply on the results of TNscope 

    (Software)
	- **SENTIEON_DIR** : Path Sentieon software 
	- **SENTIEON_LICENSE** : Sentieon license file/server
	- **BCF_DIR** : Path to bcftools a set of utilities that manipulate variant calls in the Variant Call Format (VCF) (https://samtools.github.io/bcftools/bcftools)
	- **BGZIP_DIR** : Path to bgzip (VCF files compressed with bgzip)

    (Parallelism settings)
	- **SAMPLES_PARALLEL** : Numbers of samples in parallel (set 1000 for all in parallel)
	- **NUM_THREADS** : Number of threads to use for each command
 
- ***BATCH.tsv*** (tab-delimited csv file): 
	- **SAMPLE** : Name of Sample
	- **TUMOR** : Name of Tumor
	- **FASTQ_TUMOR_PATH** : Folder path to FASTQ of Tumor
	- **NORMAL** : Name of Normal
	- **FASTQ_NORMAL_PATH** : Folder path to FASTQ of Normal

### Output file:
- **html_pipeline.zip**: Contains a dynamic HTML report executable locally