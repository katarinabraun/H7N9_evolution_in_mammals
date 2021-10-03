#!/usr/bin/env python3
#Sniffles2
#Authors: Joseph Lalli and Kelsey Florek
#Performs SNP analysis of influenza genomes from raw reads.
import yaml
import argparse
import os, sys
import time, datetime
from shutil import copyfile
import multiprocessing as mp
from time import strftime

sys.path.append(os.path.join(os.getcwd(), 'supporting_code'))

import sc
from trim import trim
from mapping import mapping,indexing,average_depth
from consensus import consensus
import readcleaning as rc
from snpcaller import snpcaller
from RePlow import RePlow
from VCF_Results_Compiler import compareVCFs
from VCFannotater import VCFannotator
from SNPGenier import SNPgenier
import fileparser as fp
from send_notification_slack import sendSlack


#print main display title
sc.mainTitle()

os.environ["DOCKER_CLIENT_TIMEOUT"] = "240"
os.environ["COMPOSE_HTTP_TIMEOUT"]="240"

#determine command line arguments and get path
parser = argparse.ArgumentParser(description='Pipeline to examine SNPs from raw illumina reads')
parser.add_argument('-c',metavar='config', type=str,help="Location of configuration file",default='config.yml')
parser.add_argument('-i',metavar='input', type=str, help="Input fasta files directory - defaults to working directory")
parser.add_argument('-o',metavar='output', type=str,help="Output directory - defaults to sniffles_files")
parser.add_argument('-t',metavar='threads', type=int,help="number of cpus to use for pipeline",default=mp.cpu_count())
parser.add_argument('--slack', metavar='slack user', type=str, help='Username of slack account to update during long runs. If left blank, won\'t send any slack messages.', default = False)

if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
args = parser.parse_args()
numThreads = args.t
configFile = args.c

#get start time
start = time.time()

#open config file and store configuation
with open(configFile,'r') as ymlFile:
	cfg = yaml.safe_load(ymlFile)

cfg['slackUser'] = args.slack

#get input path
try:
	inDir = os.path.abspath(args.i)
except (AttributeError, TypeError) as err:
	inDir = os.getcwd()
	print(f"Raw reads directory {args.i} cannot be found. Sniffles will look for fasta files in the current working directory {inDir}.\n")

#create outdir
try:
	outDir = os.path.abspath(args.o)
except (AttributeError, TypeError) as err:
	outDir = os.getcwd()
	print(f"Output directory {args.o} cannot be found. Output will be placed in a separate folder in the current working directory {outDir}.\n")

sc.checkexists(outDir)

cfg['exec']['outdir'] = os.path.join(outDir,cfg['exec']['outdir'])

try:
	os.mkdir(cfg['exec']['outdir'])
except FileExistsError:
	cfg['exec']['outdir'] = cfg['exec']['outdir']+'_'+str(int(time.time()))
	os.mkdir(cfg['exec']['outdir'])
outDir = cfg['exec']['outdir']

logfile=os.path.join(outDir,cfg['exec']['logfile'])
cfg['exec']['logfile'] = logfile

startRunMessage = f"Beginning run at {strftime('%a, %d %b %Y %I:%M:%S %p', time.localtime())}"
sc.procTitle(startRunMessage, cfg)

with open(logfile,'a') as outlog:
	outlog.write(startRunMessage + "\n")

for reference, gtf in zip(list(cfg['exec']['referenceSequences']), cfg['postprocessing']['gtfFileNames']):
	sc.procTitle (f"Processing samples for reference sequence {reference.split('.')[0]}", cfg)
	
	#assign reference sequence value and gtf value to current reference and gtf
	cfg['exec']['referenceSequence'] = os.path.join(os.getcwd(), reference)
	cfg['postprocessing']['gtfFileName'] = os.path.join(os.getcwd(), gtf)

	try:
		inDir = os.path.abspath(args.i)
	except (AttributeError, TypeError) as err:
		inDir = os.getcwd()

	#check that fasta files are arranged appropriately
	if type(cfg['exec']['referenceSequences']) == list	:
		inDir = os.path.join(inDir, reference.split(".")[0])

		assert os.path.exists(inDir), "Error: Fasta files must be in subfolders named with the sequence to be referenced against"

		#create reference-specific outdir
		cfg['exec']['outdir'] = os.path.join(outDir, reference.split(".")[0])
		try:
			os.mkdir(cfg['exec']['outdir'])
		except FileExistsError:
			cfg['exec']['outdir'] = cfg['exec']['outdir']+'_'+str(int(time.time()))
			os.mkdir(cfg['exec']['outdir'])

	#copy reference to outdir
	copyfile(os.path.join(os.getcwd(), cfg['exec']['referenceSequence']),os.path.join(cfg['exec']['outdir'],reference))

	#parse and store read information from input directory
	readData = fp.RunFiles(inDir)
	
	#trim the reads
	trim(readData,cfg,numThreads)
	assert len(readData.runtime['trimmed']) > 0, "Cannot map reads: read trimming failed"
	#setup initial mapping jobs
	mapping_list = []
	
	
	for id in readData.runtime['trimmed']:
		mapping_list.append((id,readData.runtime['trimmed'][id][0],readData.runtime['trimmed'][id][1],os.path.abspath(cfg['exec']['referenceSequence'])))
	sc.checkexists(os.path.join(cfg['exec']['outdir']+'/initial_mapping'))


	#index reference sequence
	indexing(cfg,os.path.abspath(cfg['exec']['referenceSequence']))

	#run initial mapping jobs
	sc.procTitle("Beginning initial mapping", cfg)
	bam_list = mapping(cfg,mapping_list,cfg['exec']['outdir']+'/initial_mapping',numThreads)

	#determine average depth
	assert len(bam_list) >0, "Cannot normalize coverage: 'mapping' failed to produce a list of bam files."

	#Filter by average depth coverage
	if cfg['exec']['coverageFilter']:
		bam_list = average_depth(cfg,bam_list,cfg['exec']['outdir']+'/initial_mapping',cfg['exec']['outdir']+'/coverage')

	#normalize coverage
	if cfg['exec']['normalizeCoverage']:
		assert len(bam_list) >0, sc.procTitle("Cannot normalize coverage: 'average_depth' failed to produce a list of bam files.",cfg)
		fastq_list = rc.normCoverage(cfg,bam_list,numThreads)
		mapping_list = []
		for fastq in fastq_list:
			print (fastq_list)
			read1 = fastq[0]
			read2 = fastq[1]
			if cfg['exec']['unpaired']:
				unpairedReads = fastq[2]
			else:
				unpairedReads = ""
			id = "_".join(os.path.basename(read1).split('_')[0:-1])
			mapping_list.append((id,read1,read2,os.path.abspath(cfg['exec']['referenceSequence']),unpairedReads))
		sc.checkexists(os.path.join(cfg['exec']['outdir']+'/norm_mapping'))
		
		bam_list = mapping(cfg,mapping_list,cfg['exec']['outdir']+'/norm_mapping',numThreads)

	#generate consensus
	if cfg['exec']['generateConsensus']:
		fasta_list = consensus(cfg,bam_list,numThreads)
		# map reads to consensus
		if cfg['exec']['mapToConsensus']:
			mapping_list = []
			indexing(cfg,*fasta_list)
			for fastq in fastq_list:
				read1 = fastq[0]
				read2 = fastq[1]
				if cfg['exec']['unpaired']:
					unpairedReads = fastq[3]
				else:
					unpairedReads = ""
				id = "_".join(os.path.basename(read1).split('_')[0:-1])
				mapping_list.append((id,read1,read2,os.path.abspath(cfg['exec']['referenceSequence'])))
			sc.checkexists(os.path.join(cfg['exec']['outdir']+'/norm_mapping'))

			mapping_list.append((id,readData.runtime['trimmed'][id][0],readData.runtime['trimmed'][id][1],os.path.abspath(fasta)))
			bam_list = mapping(cfg,mapping_list,cfg['exec']['outdir']+'/map_to_consensus',numThreads)

	#call snps
	assert len(bam_list) > 0, "Cannot call SNPs: Mapping reads to consensus failed."
	
	if cfg['exec']['replicates']:
		if cfg['exec']['callSNPs'] == "Varscan":
			snpcaller(cfg,bam_list,numThreads)
		elif cfg['exec']['callSNPs'] == "RePlow":
			RePlow(cfg,bam_list,numThreads)
		elif cfg['exec']['callSNPs'] == "Compare":
			replow = RePlow(cfg,bam_list,numThreads)
			varscan = snpcaller(cfg,bam_list,numThreads)
			compareVCFs(varscan, replow) #creates merged vcf that compares the output of the two SNP callers
		else:
			raise ValueError('callSNPs in config.yml must be either \'Varscan\' or \'RePlow\'')
	else:
		snpcaller(cfg,bam_list,numThreads)

	vcfpath = os.path.join(cfg['exec']['outdir'],"snp_calls")

	if cfg['exec']['annotateSNPs']:
		VCFannotator(cfg, vcfpath)

	if cfg['exec']['SNPgenier']:
		SNPgenier(cfg, vcfpath, numThreads)
	

end = time.time()
runtime = round(end - start,2)
runtime = str(datetime.timedelta(seconds=runtime))
sc.procTitle(f'Sniffles: Finished with a total runtime of {runtime}.', cfg)
