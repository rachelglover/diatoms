# diatoms
QIIME-free diatom pipeline

Diatom Pipeline
===============

These instructions assume you have a new installation of linux (or a new VM) and that all dependencies will require installation.

Installation
------------
1. Copy this folder ("diatoms") to a convenient location
2. Install the following packages with the command: sudo apt-get install python-pip zlib1g-dev parallel ncbi-blast+
3. Download PEAR from https://www.h-its.org/en/research/sco/software/#NextGenerationSequencingSequenceAnalysis and install the software in /usr/local/bin
4. Download usearch v5.2.236 from https://www.drive5.com/usearch/download.html. It is very important that you download version 5.2.236. When it is downloaded, do not change the name of the file, just move it to /usr/local/bin. Change the permissions of the file by typing: sudo chmod 755 /usr/local/bin/usearch5.2.236_i86linux32
5. Download sickle v1.33 from https://github.com/najoshi/sickle/archive/v1.33.tar.gz. Extract the folder, move into the folder on the command line and type: make. Move the sickle executable to /usr/local/bin.
6. Move into the directory above "diatoms" and install the pipeline with the command: pip install -e diatoms.

If you have had no errors at this point, the pipeline is ready to go! Please report any errors to rachel.glover@taxagenomics.com

Note: if you get an error on running python about cutadapt being missing, please try to install cutadapt with the command: sudo apt-get install python-cutadapt. It is supposed to be installed during step 9 but it appears to not be on some debian systems. sudo pip install --upgrade cutadapt also works in this scenario.

Usage
-----
Usage example 1 (easy):
sh runpipeline.sh mydata

Usage example 2 (more typing):
python pipeline.py --data mydata --dbseq database/diatoms.sequences.FINAL2017.fasta

All the forward and reverse raw MiSeq sequence files should be in the same folder. Just replace 'mydata' in the example above with the name of the folder you have put your data in. It should be within this folder or provide the full path to the folder. NOTE: The name of the folder should not contain full-stops or dashes.