# Sniffles2

Sniffles2 is a pipeline for the analysis of influenza genomes built on top of Kelsey Florek's Sniffles (https://github.com/k-florek/sniffles). It's a configurable pipeline that performs multiple functions to generate variant and consensus level information. It requires minimal dependencies as the pipeline relies on a docker container to host the software.

## Table of Contents
* [Requirements](#requirements)
* [Installing](#installing)
* [Usage](#usage)
* [FAQ](#FAQ)

### Requirements
* Linux or MacOS
* Python 3.6 or later
* Docker

### Installing
* Make sure you are using the correct version of python
* Use git to get the latest sniffles build `git clone https://github.com/JosephLalli/Sniffles2.git`
* Install the [Docker CE engine](https://docs.docker.com/install/)
* Install the required python libraries by pointing the pip installer to the requirements document in the sniffles project `pip3 install -r requirements.txt`

### Usage
Sniffles uses a configuration file to provide parameters to the pipeline. There is a default configuration file included in the sniffles project. To run the pipeline simple point the `sniffles.py` program at the correct config.yml and directory containing the raw Illumina reads.

Please see config.yml for more information.

```
 #####
#     #  #    #  #  ######  ######  #       ######   ####
#        ##   #  #  #       #       #       #       #
 #####   # #  #  #  #####   #####   #       #####    ####
      #  #  # #  #  #       #       #       #            #
#     #  #   ##  #  #       #       #       #       #    #
 #####   #    #  #  #       #       ######  ######   ####



usage: sniffles.py [-h] [-c config] [-i input] [-o output] [-t threads]

Pipeline to examine SNPs from raw illumina reads

optional arguments:
 -h, --help         show this help message and exit
 -c config          config file
 -i input           raw reads directory - defaults to working directory
 -o output          output directory - defaults to working directory
 -t threads         number of cpus to use for pipeline
 --slack slackuser  AVRL user who will be notified as Sniffles2 runs
```

Example usage:
`./sniffles.py -c config.yml -i ~/my_reads/ -t 8 -o my_results --slack joelalli`

###FAQ

What kind of input files do I need?
    - Sniffles2 uses raw FASTQ files from Illumina runs. You can name your input folder whatever you'd like.
        - If analyzing samples that need different references, those files must be in subfolders named for the reference they will be compared against.

Where should I put my reference files? Is there anything special I need to know about my reference files?
    - All reference fastas should be 2-line format (i.e., no newline characters in the sequence).
    - All reference fastas and gtfs should be in the same folder as sniffles.py

What should I do if I encounter an error?
    - Go talk to Joe. Sniffles2 is still in beta, and the quickest way to handle an issue is face-to-face.
