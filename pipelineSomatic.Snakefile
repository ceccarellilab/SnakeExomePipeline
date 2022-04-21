# conda install pyyaml numpy

configfile: "options.yaml"
configfile: "optionsSomaticMutation.yaml"

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

os.environ["BEDTOOLS_DIR"] = config['BEDTOOLS_DIR']

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

"""
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
"""

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

def get_locFromR(wildcards):
	ls = get_fastqT(".loc")
	ls = [str(my_basedir) + "/fpfilter/" + s for s in ls]
	return ls[0]

def get_readcount(wildcards):
	ls = get_fastqT(".fpfilter.pass")
	ls = [str(my_basedir) + "/fpfilter/" + s for s in ls]
	return ls

def get_annfile(wildcards):
	ls = str(my_basedir) + "/Analysis/ann.out.hg19_multianno.txt" #"/Analysis/ann.out.log"
	return ls

def get_finalannfile(wildcards):
	ls = str(my_basedir) + "/Analysis/Final.bed"
	return ls

def get_merscore(wildcards):
	ls = str(my_basedir) + "/Analysis/" + config["MER_SCORE"]
	return ls

def get_SomaticVariants(wildcards):
	ls = str(my_basedir) + "/Analysis/SomaticVariants.bed"
	return ls

def get_FinalSomaticVariants(wildcards):
	ls = str(my_basedir) + "/Analysis/SomaticVariants.txt"
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
		get_FinalSomaticVariants

rule locFromR:
	input:
		
	output:
		"{dir}/{sample}.loc"
	params:
		workdirr = config["WORK_DIR"],
		vcfdir =  config["VCF_DIR"], 
		vcffile = config["VCF_FILE"],
		threads= config["NUM_FORKS"]
	conda:
		"envs/environmentR.yml"
	shell:
		"""
		Rscript {workflow.basedir}/R/script1.R {workflow.basedir} {params.workdirr} {params.vcfdir} {params.vcffile} {params.threads}
		"""

rule readCount:
	input:
		get_locFromR
	output:
		"{dir}/{sample}.fpfilter.pass"
	params:
		fastaPath = config["FASTA_FILE"],
		workPath =  config["WORK_DIR"], 
		bamreadcount = config["BAMREADCOUNT_FILE"],
		fpfilterfile = config["FPFILTER_FILE"]
	shell:
		" sh {workflow.basedir}/scripts/bamreadcount.sh {workflow.basedir}/fpfilter {params.fastaPath} {params.workPath} {params.bamreadcount} {wildcards.sample}.loc {wildcards.sample}.var {params.fpfilterfile}"

rule ANNOVARAnnotation:
	input:
		get_readcount
	output:
		"{basedir}/Analysis/ann.out.hg19_multianno.txt" #.log
	shell:
		"""
		Rscript {workflow.basedir}/R/script2.R {workflow.basedir}
		cd {workflow.basedir}/Analysis/
		/storage/qnap_vol1/data/PUBLIC/Tools/annovar2017Jul16/table_annovar.pl ann.in /storage/qnap_vol1/data/PUBLIC/Tools/annovar2017Jul16/humandb/ -buildver hg19 -out ann.out -remove -protocol refGene,avsnp150,snp138NonFlagged,exac03,1000g2015aug_all,esp6500siv2_all,kaviar_20150923,hrcr1,cosmic80,clinvar_20170905,dbnsfp33a,dbscsnv11,dann,eigen,gerp++gt2,cadd  -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . --thread 80 || true
		"""

rule FinalAnnotation:
	input:
		get_annfile
	params:
		ponfilter=config["PONFILTER_FILE"],
		pondir=config["PON_DIR"]
	output:
		"{basedir}/Analysis/Final.bed"
	conda:
		"envs/environmentR.yml"
	shell:
		"""
		cd {workflow.basedir}/Analysis/
		Rscript {workflow.basedir}/R/script3.R {workflow.basedir} {params.ponfilter} {params.pondir}
		"""

rule merScore:
	input:
		get_finalannfile
	output:
		"{basedir}/Analysis/SomaticVariants.bed"
	params:
		merPath= config["MER_FILE"] ,
		MER_SCORE= config["MER_SCORE"], 
		bedtools= config["BEDTOOLS_PATH"], 
		FLAGS_FILE= config["FLAGS_FILE"], 
		cbioPath = config["CBIO_FILE"],
		CBIO_SCORE = config["CBIO_SCORE"],
		maskpath = config["MASK_FILE"],
		MASK_SCORE = config["MASK_SCORE"]
	shell:
		"""
		cd {workflow.basedir}/Analysis/
		{params.bedtools} intersect -wo -a Final.bed -b {params.merPath} > {params.MER_SCORE}
		Rscript {workflow.basedir}/R/script4.R {params.MER_SCORE}
		{params.bedtools} intersect -wao -a Finalmask.bed -b {params.maskpath} > {params.MASK_SCORE}
		Rscript {workflow.basedir}/R/script5.R {params.MER_SCORE} {params.FLAGS_FILE} {params.MASK_SCORE}
		{params.bedtools} intersect -wao -f 1 -a SomaticVariants.bed -b {params.cbioPath} > {params.CBIO_SCORE}
		"""

rule somaticVariant:
	input:
		get_SomaticVariants
	output:
		"{basedir}/Analysis/SomaticVariants.txt"
	params:
		CBIO_DB= config["CBIO_DB"] ,
		CBIO_SCORE = config["CBIO_SCORE"],
		NUM_FORKS= config["NUM_FORKS"], 
		CENSUS_FILE= config["CENSUS_FILE"], 
		NATURE_FILE= config["NATURE_FILE"], 
		ONCOKB_FILE = config["ONCOKB_FILE"],
		SMG_FILE = config["SMG_FILE"]
	conda:
		"envs/environmentR.yml"
	shell:
		"""
		cd {workflow.basedir}/Analysis/
		Rscript {workflow.basedir}/R/script6.R {params.CBIO_DB} {params.CBIO_SCORE} {params.NUM_FORKS} {params.CENSUS_FILE} {params.NATURE_FILE} {params.ONCOKB_FILE} {params.SMG_FILE}
		"""