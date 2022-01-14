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
		nthreads = config['NUM_THREADS'],
		fastqsuffix1 = config['FASTQ_SUFFIXES']["1_SUFFIX"],
		fastqsuffix2 = config['FASTQ_SUFFIXES']["2_SUFFIX"],
	shell:
		"({params.sentionpath}/bin/sentieon bwa mem -M -R '@RG\\tID:{wildcards.sample}{wildcards.NorT}\\tSM:{wildcards.sample}{wildcards.NorT}\\tPL:{params.platform}' -t {params.nthreads} -K 10000000 {params.fastapath} {input} || echo -n error ) | {params.sentionpath}/bin/sentieon util sort -o {wildcards.sample}{wildcards.NorT}_sorted.bam -t {params.nthreads} --sam2bam -i -"   	