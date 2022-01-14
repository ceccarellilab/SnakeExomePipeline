# conda install pyyaml numpy

configfile: "options.yaml"

my_basedir = workflow.current_basedir

# --- functions to create lists of snakemake target files ---

import glob
import re
import os


#Set working directory

config['home'] = os.getcwd()+"/"
#config['workdir'] = "/home/adefalco/"+ "/" + config['BATCH'] + "/"
config['workdir'] = config['home'] + config['BATCH'] + "/"

if not os.path.isdir(config['workdir']): 
	os.mkdir(config['workdir']) 

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
   
i = 0
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
    i = i+1


os.chdir(config['workdir'])


def get_fastqN():
	ls = []
	for i in range(normSize):
		SAMPLE = uniqNorm[i]
		#display(Markdown("### Normal: " + SAMPLE))
		fastqNor = fastqNorm[i,2]
		norm = fastqNorm[i,1]
		
		#fastq_str = " "
		#for key in config['FASTQ_SUFFIXES']:
		#        fastq_str = fastq_str + fastqNor + norm + "_"+ config['FASTQ_SUFFIXES'][key] +" "  
		ls.append(fastqNor + norm + "_"+ config['FASTQ_SUFFIXES']["1_SUFFIX"]) 
	return ls    

def get_fastqT():
	ls = []
	for i in range(batchSize):
		TUMOR = batch[i][1]
		#display(Markdown("### Tumor: "+TUMOR))
		
		fastq_str = " "
		#for key in config['FASTQ_SUFFIXES']:
		#    fastq_str = fastq_str + fastqTum[i] + TUMOR + "_"+ config['FASTQ_SUFFIXES'][key] +" " 
		ls.append(TUMOR + "T_sorted.bam") 

	return ls


def get_sortbam(wildcards):
	ls = get_fastqN() + get_fastqT()
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

rule all:
	input:
		get_sortbam

# 1 Mapping reads with BWA-MEM, sorting for normal sample
rule mapping:
	input:
		get_fastqFile
	output:
		"{sample}{NorT}_sorted.bam"
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		platform = config['PLATFORM'],
		nthreads = config['NUM_THREADS']
	shell:
		"({params.sentionpath}/bin/sentieon bwa mem -M -R '@RG\\tID:{wildcards.sample}{wildcards.NorT}\\tSM:{wildcards.sample}{wildcards.NorT}\\tPL:{params.platform}' -t {params.nthreads} -K 10000000 {params.fastapath} {input} || echo -n error ) | {params.sentionpath}/bin/sentieon util sort -o {wildcards.sample}{wildcards.NorT}_sorted.bam -t {params.nthreads} --sam2bam -i -" 

rule metrics:
input:
		get_fastqFile
	output:
		"{sample}{NorT}_sorted.bam"
	params:
		sentionpath =config['SENTIEON_DIR'],
		fastapath = config['FASTA'],
		platform = config['PLATFORM'],
		nthreads = config['NUM_THREADS']
	shell:
		""
	    "{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_sorted.bam --algo MeanQualityByCycle {wildcards.sample}{wildcards.NorT}_mq_metrics.txt --algo QualDistribution {wildcards.sample}{wildcards.NorT}_qd_metrics.txt --algo GCBias --summary {wildcards.sample}{wildcards.NorT}_gc_summary.txt {wildcards.sample}{wildcards.NorT}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' {wildcards.sample}{wildcards.NorT}_aln_metrics.txt --algo InsertSizeMetricAlgo {wildcards.sample}{wildcards.NorT}_is_metrics.txt"
		
		#ONLYTUMOR
		"{params.sentionpath}/bin/sentieon plot GCBias -o {wildcards.sample}{wildcards.NorT}_gc-report.pdf {wildcards.sample}{wildcards.NorT}_gc_metrics.txt"

rule removeDuplicate:

"{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_sorted.bam --algo LocusCollector --fun score_info {wildcards.sample}{wildcards.NorT}_score.txt"
"{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_sorted.bam --algo Dedup --rmdup --score_info {wildcards.sample}{wildcards.NorT}_score.txt --metrics {wildcards.sample}{wildcards.NorT}_dedup_metrics.txt {wildcards.sample}{wildcards.NorT}_deduped.bam"

rule indelRealigner:

"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_deduped.bam --algo Realigner"+know_site_str+SAMPLE+"N_realigned.bam"

rule baseRecalibration:

"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_realigned.bam --algo QualCal -k "+data['DBSNP']+know_site_str+SAMPLE+"N_recal_data.table"
"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_realigned.bam -q {wildcards.sample}{wildcards.NorT}_recal_data.table --algo QualCal -k "+data['DBSNP']+know_site_str+SAMPLE+"N_recal_data.table.post"
"{params.sentionpath}/bin/sentieon driver -t {params.nthreads} --algo QualCal --plot --before {wildcards.sample}{wildcards.NorT}_recal_data.table --after {wildcards.sample}{wildcards.NorT}_recal_data.table.post {wildcards.sample}{wildcards.NorT}_recal.csv"
"{params.sentionpath}/bin/sentieon plot QualCal -o {wildcards.sample}{wildcards.NorT}_recal_plots.pdf {wildcards.sample}{wildcards.NorT}_recal.csv"
"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_realigned.bam -q {wildcards.sample}{wildcards.NorT}_recal_data.table --algo ReadWriter {wildcards.sample}{wildcards.NorT}_recal.bam"

rule HC_VariantCaller:

"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}{wildcards.NorT}_recal.bam --algo Haplotyper -d "+data['DBSNP']+" --emit_conf=30 --call_conf=30 {wildcards.sample}N-output-hc.vcf.gz"

rule DNAscope_variantCaller:

"{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -r {params.fastapath} -i {wildcards.sample}{wildcards.NorT}_recal.bam --algo DNAscope -d "+data['DBSNP']+" --model "+data['ML_MODEL_N']+" {wildcards.sample}N-tmpDNAscope.vcf.gz"
"{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -r {params.fastapath} --algo DNAModelApply --model "+data['ML_MODEL_N']+" -v {wildcards.sample}N-tmpDNAscope.vcf.gz {wildcards.sample}N-DNAscope.vcf.gz"
"{params.bcfpath} filter -s ML_FAIL -i INFO/ML_PROB > 0.81 {wildcards.sample}N-DNAscope.vcf.gz -O z -m x -o {wildcards.sample}N-filtDNAscope.vcf.gz"

rule variantAnnotation:

"/storage/gluster/vol1/bcbio/anaconda/bin/snpEff -Xms1000m -Xmx36400m -Djava.io.tmpdir="+data['SENTIEON_TMPDIR']+" eff -noStats -t -noLog -dataDir /storage/gluster/vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 {wildcards.sample}N-output-hc.vcf.gz | "+data['BGZIP_DIR']+" --threads {params.nthreads} -c > {wildcards.sample}N-output-hc.snpEff.vcf.gz"
"/storage/gluster/vol1/bcbio/anaconda/bin/snpEff -Xms1000m -Xmx36400m -Djava.io.tmpdir="+data['SENTIEON_TMPDIR']+" eff -noStats -t -noLog -dataDir /storage/gluster/vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 {wildcards.sample}N-filtDNAscope.vcf.gz | "+data['BGZIP_DIR']+" --threads {params.nthreads} -c > {wildcards.sample}N-filtDNAscope.snpEff.vcf.gz"

rule somaticVariantCallingTNseq:

"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i "+TUMOR+"T_recal.bam -i {wildcards.sample}{wildcards.NorT}_recal.bam --algo TNsnv --tumor_sample "+TUMOR+"T --normal_sample {wildcards.sample}N --pon "+data['PON_TNsnv'] +" --cosmic "+data['COSMIC_DB']+" --dbsnp "+data['DBSNP']+" --call_stats_out "+TUMOR+"-call.stats "+TUMOR+"-TNsnv.vcf.gz"
"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i "+TUMOR+"T_recal.bam -i {wildcards.sample}{wildcards.NorT}_recal.bam --algo TNhaplotyper --tumor_sample "+TUMOR+"T --normal_sample {wildcards.sample}N --pon "+data['PON_TNhaplotyper']+" --cosmic "+data['COSMIC_DB']+" --dbsnp "+data['DBSNP']+" "+TUMOR+"-TNhaplotyper.vcf.gz"

rule somaticVariantCallingTNscope:

"{params.sentionpath}/bin/sentieon driver -r {params.fastapath} -t {params.nthreads} -i "+TUMOR+"T_recal.bam -i {wildcards.sample}{wildcards.NorT}_recal.bam --algo TNscope --tumor_sample "+TUMOR+"T --normal_sample {wildcards.sample}N --dbsnp "+data['DBSNP']+" --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 --min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.00005 "+TUMOR+"-tmpTNscope.vcf.gz"
"{params.sentionpath}/bin/sentieon driver -t {params.nthreads} -r {params.fastapath} --algo TNModelApply --model "+data['ML_MODEL_T'] +" -v "+TUMOR+"-tmpTNscope.vcf.gz "+TUMOR+"-TNscope.vcf.gz"
"{params.bcfpath} filter -s ML_FAIL -i \INFO/ML_PROB > 0.81 "+TUMOR+"-TNscope.vcf.gz -O z -m x -o "+TUMOR+ "-filtTNscope.vcf.gz"

rule somaticVariantAnnotation:

"/storage/gluster/vol1/bcbio/anaconda/bin/snpEff -Xms1000m -Xmx36400m -Djava.io.tmpdir="+data['SENTIEON_TMPDIR']+" eff -noStats -t -noLog -dataDir /storage/gluster/vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 "+TUMOR+"-TNsnv.vcf.gz | "+data['BGZIP_DIR']+" --threads {params.nthreads} -c > "+TUMOR+"-TNsnv.snpEff.vcf.gz"
"/storage/gluster/vol1/bcbio/anaconda/bin/snpEff -Xms1000m -Xmx36400m -Djava.io.tmpdir="+data['SENTIEON_TMPDIR']+" eff -noStats -t -noLog -dataDir /storage/gluster/vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 "+TUMOR+"-TNhaplotyper.vcf.gz | "+data['BGZIP_DIR']+" --threads {params.nthreads} -c > "+TUMOR+"-TNhaplotyper.snpEff.vcf.gz"
"/storage/gluster/vol1/bcbio/anaconda/bin/snpEff -Xms1000m -Xmx36400m -Djava.io.tmpdir="+data['SENTIEON_TMPDIR']+" eff -noStats -t -noLog -dataDir /storage/gluster/vol1/bcbio/genomes/Hsapiens/hg19/snpeff -hgvs -noLof -i vcf -o vcf -noInteraction -noMotif -noNextProt -strict GRCh37.75 "+TUMOR+"-filtTNscope.vcf.gz | "+data['BGZIP_DIR']+" --threads {params.nthreads} -c > "+TUMOR+"-filtTNscope.snpEff.vcf.gz"
