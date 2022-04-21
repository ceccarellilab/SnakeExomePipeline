# conda install pyyaml numpy
report: "report/workflow.rst"

configfile: "options.yaml"

my_basedir = workflow.current_basedir

# --- functions to create lists of snakemake target files ---

import glob
import re
import os
import shutil

#Set working directory

config['home'] = os.getcwd()+"/"
#config['workdir'] = "/home/adefalco/"+ "/" + config['BATCH'] + "/"
config['workdir'] = config['home'] + config['BATCH'] + "/"

if not os.path.isdir(config['workdir']): 
	os.mkdir(config['workdir']) 

shutil.copy("options.yaml", config['workdir']+"options.yaml")

os.environ["SENTIEON_DIR"] = config['SENTIEON_DIR']
os.environ["SENTIEON_LICENSE"]= config['SENTIEON_LICENSE']
os.environ["BCFTOOLS_PLUGINS"]= config['BCF_DIR']

#from IPython.display import HTML

#Read sample list

#config['samplelist'] = config['home']+config['BATCH']+".tsv"

import csv
import numpy as np

batch = []
with open (config['samplelist'], 'r') as f:
	for row in csv.reader(f,delimiter='\t'):
			batch.append(row)
			
print(type(batch))

batch = np.array(batch)
batch = np.delete(batch, 0, 0)
batchSize = int(np.size(batch)/np.size(batch,1)) 

fastqNorm = np.column_stack((batch[:,0],batch[:,3],batch[:,4]))
fastqTum = batch[:,2]

ERROR = "error"
commands = []


tumors = [row[1] for row in batch]
fastqNorm = np.unique(fastqNorm, axis=0)
print(fastqNorm[:,0])

uniqNorm = fastqNorm[:,0]#np.unique(normals)

normSize = np.size(uniqNorm)

if config['SAMPLES_PARALLEL']==1000:
    config['SAMPLES_PARALLEL'] = batchSize+normSize
   
config['KNOWN_SITES_STR'] = " "
if config['KNOWN_SITES'] is not None:
    for key in config['KNOWN_SITES']:
        config['KNOWN_SITES_STR'] = config['KNOWN_SITES_STR'] +"-k " + config['KNOWN_SITES'][key] +" "

config['KNOWN_SITES_STR']

""" i = 0
for file in fastqNorm[:,2]:
    if isinstance(file,np.ndarray):
        file = file[0]
        #norm = fastqNorm[i,1]
    #else:
        #norm = fastqNorm[0,1]      
    pathFile = file+fastqNorm[i,1]+'_'
    
    for key in config['FASTQ_SUFFIXES']:
        print(pathFile+config['FASTQ_SUFFIXES'][key])
        if not os.path.isfile(pathFile+config['FASTQ_SUFFIXES'][key]):
          raise Exception("Sorry, fastq not found :" + pathFile+config['FASTQ_SUFFIXES'][key])
    i = i+1

i = 0
for file in fastqTum:
    #if isinstance(file,np.ndarray):
    #    file = file[0]
    pathFile = file+tumors[i]+'_'
    
    for key in config['FASTQ_SUFFIXES']:
        print(pathFile+config['FASTQ_SUFFIXES'][key])
        if not os.path.isfile(pathFile+config['FASTQ_SUFFIXES'][key]):
            raise Exception("Sorry, fastq not found :" + pathFile+config['FASTQ_SUFFIXES'][key])
    i = i+1 """


os.chdir(config['workdir'])


def get_fastqN(namefile):
	ls = []
	for i in range(normSize):
		SAMPLE = uniqNorm[i]
		#display(Markdown("### Normal: " + SAMPLE))
		fastqNor = fastqNorm[i,2]
		norm = fastqNorm[i,1]
		
		#fastq_str = " "
		#for key in config['FASTQ_SUFFIXES']:
		#        fastq_str = fastq_str + fastqNor + norm + "_"+ config['FASTQ_SUFFIXES'][key] +" "  
		#ls.append(fastqNor + norm + "_"+ config['FASTQ_SUFFIXES']["1_SUFFIX"]) 
		ls.append(norm + namefile) 
	return ls    

def get_fastqT(namefile):
	ls = []
	for i in range(batchSize):
		TUMOR = batch[i][1]
		#display(Markdown("### Tumor: "+TUMOR))
		
		fastq_str = " "
		#for key in config['FASTQ_SUFFIXES']:
		#    fastq_str = fastq_str + fastqTum[i] + TUMOR + "_"+ config['FASTQ_SUFFIXES'][key] +" " 
		ls.append(TUMOR + namefile) 

	return ls


def get_sortbam(wildcards):
	ls = get_fastqN("N.bam") + get_fastqT("T.bam")
	return ls

def get_metrics(wildcards):
	#ls = get_fastqN("N_gc_metrics.txt") + get_fastqT("T_gc_metrics.txt")
	ls = get_fastqN("N_gc-report.pdf") + get_fastqT("T_gc-report.pdf")
	return ls

def get_duplicate(wildcards):
	ls = get_fastqN("N_deduped.bam") + get_fastqT("T_deduped.bam")
	return ls

def get_indelRealigner(wildcards):
	ls = get_fastqN("N_realigned.bam") + get_fastqT("T_realigned.bam")
	return ls

def get_baseRecalibration(wildcards):
	ls = get_fastqN("N_recal_plots.pdf") + get_fastqT("T_recal_plots.pdf")
	return ls

def get_HC_VariantCaller(wildcards):
	ls = get_fastqN("N-output-hc.vcf.gz") + get_fastqT("T-output-hc.vcf.gz")
	return ls

def get_DNAscope_variantCaller(wildcards):
	ls = get_fastqN("N-filtDNAscope.vcf.gz") + get_fastqT("T-filtDNAscope.vcf.gz")
	return ls

def get_variantAnnotation(wildcards):
	ls = get_fastqN("N-filtDNAscope.snpEff.vcf.gz") + get_fastqT("T-filtDNAscope.snpEff.vcf.gz")
	return ls

def get_somaticVariantCallingTNseq(wildcards):
	ls = get_fastqT("-TNhaplotyper.vcf.gz")
	return ls

def get_somaticVariantCallingTNscope(wildcards):
	ls = get_fastqT("-filtTNscope.vcf.gz")
	return ls

def get_somaticVariantAnnotation(wildcards):
	ls = get_fastqT("-filtTNscope.snpEff.vcf.gz")
	return ls

def get_fastqFile(wildcards):
	ls = []
	sample = unpack({wildcards.sample})[0]
	NorT = unpack({wildcards.NorT})[0]
	if	NorT=='N':
		ind = np.where(batch[:,3] == sample)[0][0]
		path = batch[ind][4]
	else:
		ind = np.where(batch[:,1] == sample)[0][0]
		path = batch[ind][2]
	ls.append(path + sample + "_"+ config['FASTQ_SUFFIXES']["1_SUFFIX"]) 
	ls.append(path + sample + "_"+ config['FASTQ_SUFFIXES']["2_SUFFIX"]) 
	return ls

def getControlSample(wildcards):
	ls = []
	sample = unpack({wildcards.sample})[0]
	ind = np.where(batch[:,1] == sample)[0][0]
	path = batch[ind][3]
	ls.append(path) 
	return ls

rule all:
	input:
		get_somaticVariantAnnotation

rule metrics:
	input:
		get_sortbam
	output:
		report("{sample}{NorT}_gc-report.pdf", category="metrics")
		#"{sample}{NorT}_gc-report.pdf"
		#"{sample}{NorT}_gc_metrics.txt"
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		nthreads = config['NUM_THREADS']
	shell:
		"""
	    {params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}.bam --algo MeanQualityByCycle {wildcards.sample}{wildcards.NorT}_mq_metrics.txt --algo QualDistribution {wildcards.sample}{wildcards.NorT}_qd_metrics.txt --algo GCBias --summary {wildcards.sample}{wildcards.NorT}_gc_summary.txt {wildcards.sample}{wildcards.NorT}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' {wildcards.sample}{wildcards.NorT}_aln_metrics.txt --algo InsertSizeMetricAlgo {wildcards.sample}{wildcards.NorT}_is_metrics.txt
		{params.sentionpath}/bin/sentieon plot GCBias -o {wildcards.sample}{wildcards.NorT}_gc-report.pdf {wildcards.sample}{wildcards.NorT}_gc_metrics.txt
		"""

rule removeDuplicate:
	input:
		get_metrics
	output:
		"{sample}{NorT}_deduped.bam"
	params:
		sentionpath =config['SENTIEON_DIR'],
		nthreads = config['NUM_THREADS']
	shell:
		"""
		{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}.bam --algo LocusCollector --fun score_info {wildcards.sample}{wildcards.NorT}_score.txt
		{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}.bam --algo Dedup --rmdup --score_info {wildcards.sample}{wildcards.NorT}_score.txt --metrics {wildcards.sample}{wildcards.NorT}_dedup_metrics.txt {wildcards.sample}{wildcards.NorT}_deduped.bam
		"""


rule indelRealigner:
	input:
		get_duplicate
	output:
		"{sample}{NorT}_realigned.bam"
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		knowsitestr = config['KNOWN_SITES_STR'],
		nthreads = config['NUM_THREADS']
	shell:
		"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_deduped.bam --algo Realigner{params.knowsitestr}{wildcards.sample}{wildcards.NorT}_realigned.bam"

rule baseRecalibration:
	input:
		get_indelRealigner
	output:
		#"{sample}{NorT}_recal.bam"
		report("{sample}{NorT}_recal_plots.pdf", category="Recal plot")
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		knowsitestr = config['KNOWN_SITES_STR'],
		dbsnppath = config['DBSNP'],
		nthreads = config['NUM_THREADS']
	shell:
		"""
		{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_realigned.bam --algo QualCal -k {params.dbsnppath}{params.knowsitestr}{wildcards.sample}{wildcards.NorT}_recal_data.table
		{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_realigned.bam -q {wildcards.sample}{wildcards.NorT}_recal_data.table --algo QualCal -k {params.dbsnppath}{params.knowsitestr}{wildcards.sample}{wildcards.NorT}_recal_data.table.post
		{params.sentionpath}/bin/sentieon driver -t {params.nthreads} --algo QualCal --plot --before {wildcards.sample}{wildcards.NorT}_recal_data.table --after {wildcards.sample}{wildcards.NorT}_recal_data.table.post {wildcards.sample}{wildcards.NorT}_recal.csv
		{params.sentionpath}/bin/sentieon plot QualCal -o {wildcards.sample}{wildcards.NorT}_recal_plots.pdf {wildcards.sample}{wildcards.NorT}_recal.csv
		{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_realigned.bam -q {wildcards.sample}{wildcards.NorT}_recal_data.table --algo ReadWriter {wildcards.sample}{wildcards.NorT}_recal.bam
		"""


rule HC_VariantCaller:
	input:
		get_baseRecalibration
	output:
		"{sample}{NorT}-output-hc.vcf.gz"
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		dbsnppath = config['DBSNP'],
		nthreads = config['NUM_THREADS']
	shell:
		"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_recal.bam --algo Haplotyper -d {params.dbsnppath} --emit_conf=30 --call_conf=30 {wildcards.sample}{wildcards.NorT}-output-hc.vcf.gz"


rule DNAscope_variantCaller:
	input:
		get_HC_VariantCaller
	output:
		"{sample}{NorT}-filtDNAscope.vcf.gz"
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		dbsnppath = config['DBSNP'],
		bcfpath = config['BCF_DIR'],
		MLmodel = config['ML_MODEL_N'],
		nthreads = config['NUM_THREADS']
	shell:
		"""
		{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -r {params.fastapath} -i {wildcards.sample}{wildcards.NorT}_recal.bam --algo DNAscope -d {params.dbsnppath} --model {params.MLmodel} {wildcards.sample}{wildcards.NorT}-tmpDNAscope.vcf.gz
		{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -r {params.fastapath} --algo DNAModelApply --model {params.MLmodel} -v {wildcards.sample}{wildcards.NorT}-tmpDNAscope.vcf.gz {wildcards.sample}{wildcards.NorT}-DNAscope.vcf.gz
		{params.bcfpath} filter -s ML_FAIL -i INFO/ML_PROB > 0.81 {wildcards.sample}{wildcards.NorT}-DNAscope.vcf.gz -O z -m x -o {wildcards.sample}{wildcards.NorT}-filtDNAscope.vcf.gz
		"""

rule variantAnnotation:
	input:
		get_DNAscope_variantCaller
	output:
		"{sample}{NorT}-filtDNAscope.snpEff.vcf.gz"
	params:
		sentionpath =config['SENTIEON_DIR'],
		bgzippath = config['BGZIP_DIR'],
		nthreads = config['NUM_THREADS'],
		snpEffpath = config['SNPEFF_DIR']
	conda:
		"envs/environment.yml"
	shell:
		"""
		/storage/qnap_vol1/bcbio/anaconda/pkgs/java-jdk-8.0.92-1/bin/java -Xmx8g -jar /home/adabbo/gdc/SnakeFromBAM/SnakeExomePipeline/snpEff/snpEff.jar eff -noStats -t -noLog -dataDir /storage/qnap_vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 {wildcards.sample}{wildcards.NorT}-output-hc.vcf.gz | {params.bgzippath} --threads {params.nthreads} -c > {wildcards.sample}{wildcards.NorT}-output-hc.snpEff.vcf.gz || true
		/storage/qnap_vol1/bcbio/anaconda/pkgs/java-jdk-8.0.92-1/bin/java -Xmx8g -jar /home/adabbo/gdc/SnakeFromBAM/SnakeExomePipeline/snpEff/snpEff.jar eff -noStats -t -noLog -dataDir /storage/qnap_vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 {wildcards.sample}{wildcards.NorT}-filtDNAscope.vcf.gz | {params.bgzippath} --threads {params.nthreads} -c > {wildcards.sample}{wildcards.NorT}-filtDNAscope.snpEff.vcf.gz || true
		"""

rule somaticVariantCallingTNseq:
	input:
		get_variantAnnotation
	output:
		"{sample}-TNhaplotyper.vcf.gz"
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		nthreads = config['NUM_THREADS'],
		controlSample = lambda wildcards: getControlSample(wildcards),
		PON_TNsvn = config['PON_TNsnv'],
		dbsnppath = config['DBSNP'],
		cosmicdb = config['COSMIC_DB'],
		PON_TNhaplotyper = config['PON_TNhaplotyper']
	shell:
		"""
		{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}T_recal.bam -i {params.controlSample}N_recal.bam --algo TNsnv --tumor_sample {wildcards.sample} --normal_sample {params.controlSample} --pon {params.PON_TNsvn} --cosmic {params.cosmicdb} --dbsnp {params.dbsnppath} --call_stats_out {wildcards.sample}-call.stats {wildcards.sample}-TNsnv.vcf.gz
		{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}T_recal.bam -i {params.controlSample}N_recal.bam --algo TNhaplotyper --tumor_sample {wildcards.sample} --normal_sample {params.controlSample} --pon {params.PON_TNhaplotyper} --cosmic {params.cosmicdb} --dbsnp {params.dbsnppath} {wildcards.sample}-TNhaplotyper.vcf.gz
		"""

rule somaticVariantCallingTNscope:
	input:
		get_somaticVariantCallingTNseq
	output:
		"{sample}-filtTNscope.vcf.gz"
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		nthreads = config['NUM_THREADS'],
		controlSample = lambda wildcards: getControlSample(wildcards),
		dbsnppath = config['DBSNP'],
		mlmodelT = config['ML_MODEL_T'],
		bcfpath = config['BCF_DIR']
	shell:
		"""
		{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}T_recal.bam -i {params.controlSample}N_recal.bam --algo TNscope --tumor_sample {wildcards.sample} --normal_sample {params.controlSample} --dbsnp {params.dbsnppath} --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 --min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.00005 {wildcards.sample}-tmpTNscope.vcf.gz
		{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -r {params.fastapath} --algo TNModelApply --model {params.mlmodelT} -v {wildcards.sample}-tmpTNscope.vcf.gz {wildcards.sample}-TNscope.vcf.gz
		{params.bcfpath} filter -s ML_FAIL -i \INFO/ML_PROB > 0.81 {wildcards.sample}-TNscope.vcf.gz -O z -m x -o {wildcards.sample}-filtTNscope.vcf.gz
		"""


rule somaticVariantAnnotation:
	input:
		get_somaticVariantCallingTNscope
	output:
		"{sample}-filtTNscope.snpEff.vcf.gz"
	params:
		bgzippath = config['BGZIP_DIR'],
		nthreads = config['NUM_THREADS'],
		snpEffpath = config['SNPEFF_DIR']
	conda:
		"envs/environment.yml"
	shell:
		"""
        /storage/qnap_vol1/bcbio/anaconda/pkgs/java-jdk-8.0.92-1/bin/java -Xmx8g -jar /home/adabbo/gdc/SnakeFromBAM/SnakeExomePipeline/snpEff/snpEff.jar eff -noStats -t -noLog -dataDir /storage/gluster/vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 {wildcards.sample}-TNsnv.vcf.gz | {params.bgzippath} --threads {params.nthreads} -c > {wildcards.sample}-TNsnv.snpEff.vcf.gz || true
		/storage/qnap_vol1/bcbio/anaconda/pkgs/java-jdk-8.0.92-1/bin/java -Xmx8g -jar /home/adabbo/gdc/SnakeFromBAM/SnakeExomePipeline/snpEff/snpEff.jar eff -noStats -t -noLog -dataDir /storage/gluster/vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 {wildcards.sample}-TNhaplotyper.vcf.gz | {params.bgzippath} --threads {params.nthreads} -c > {wildcards.sample}-TNhaplotyper.snpEff.vcf.gz || true
		/storage/qnap_vol1/bcbio/anaconda/pkgs/java-jdk-8.0.92-1/bin/java -Xmx8g -jar /home/adabbo/gdc/SnakeFromBAM/SnakeExomePipeline/snpEff/snpEff.jar eff -noStats -t -noLog -dataDir /storage/gluster/vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 {wildcards.sample}-filtTNscope.vcf.gz | {params.bgzippath} --threads {params.nthreads} -c > {wildcards.sample}-filtTNscope.snpEff.vcf.gz || true
        """
