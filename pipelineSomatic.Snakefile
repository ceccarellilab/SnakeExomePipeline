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
	return ls

def get_readcount(wildcards):
	ls = get_fastqT(".log")
	ls = [str(my_basedir) + "/fpfilter/" + s for s in ls]
	return ls

print(get_locFromR(""))

def getControlSample(wildcards):
	ls = []
	sample = unpack({wildcards.sample})[0]
	ind = np.where(batch[:,1] == sample)[0][0]
	path = batch[ind][3]
	ls.append(path) 
	return ls

rule all:
	input:
		get_readcount
		

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
		"{dir}/{sample}.log"
	params:
		fastaPath = config["FASTA_FILE"],
		workPath =  config["WORK_DIR"], 
		bamreadcount = config["BAMREADCOUNT_FILE"]
	shell:
		" sh {workflow.basedir}/scripts/bamreadcount.sh {workflow.basedir}/fpfilter {params.fastaPath} {params.workPath} {params.bamreadcount} {wildcards.sample}.loc {wildcards.sample}.var"
