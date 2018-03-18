from __future__ import division
import sys
import os
import pygal
import argparse
import re
import glob
import subprocess
import pandas
from Bio import SeqIO
from collections import Counter
from collections import defaultdict


def main():
    forwardprimer = "ATGCGTTGGAGAGARCGTTTC"
    reverseprimer = "GATCACCTTCTAATTTACCWACAACTG"

    # Installation checks. If any fail the software will not proceed
    checkResult = dependencyChecks()

    # Parse the arguments from the command line
    options = processArguments()
    dataDirectory = os.path.realpath(options.data) + "/"

    # Run the amplicon QC
    clusteringFile = ampliconQC(dataDirectory, options.threads, forwardprimer, reverseprimer)
    # run the OTU picking (with usearch rather than QIIME uclust. uclust is ONLY licensed
    # for use with QIIME so we can't use it.)
    clustered = clustering(clusteringFile, dataDirectory)

    # Create the repset sequences and counts data
    counts = repsetCounts(clusteringFile, dataDirectory)

    # Assign taxonomy to the repset sequences
    taxonomy = assignTaxonomy(options.threads, options.dbseq, dataDirectory)

    # Calculate abundances
    abundances = calculateAbundance(counts, dataDirectory)

    # Output pass/fail excel files.
    output = outputReport(abundances, dataDirectory)


def outputReport(abundances, directory):
    """
    This function outputs the pass.csv and fail.csv abundance files.
    """
    repeats = []
    passes = []
    for sample in abundances:
        total = abundances[sample].sum()
        if total <= 3000:
            repeats.append(sample)
        if total > 3000:
            passes.append(sample)

    passed = abundances[passes]
    repeat = abundances[repeats]

    passed.to_csv(directory + "analysis/Abundances.pass.csv")
    repeat.to_csv(directory + "analysis/Abundances.fail.csv")


def calculateAbundance(countsDF, directory):
    """
    This function calculates the abundance of each OTU in each sample and the
    associated taxonomy for output later.
    """
    abundancesDict = defaultdict(lambda: defaultdict(float))
    for line in open(directory + "analysis/repset.taxonomy.txt", "rU"):
        linesplit = line.split('\t')
        otu = linesplit[0].rstrip()
        assigned_taxonomy = linesplit[1].strip(';')
        otudata = countsDF.loc[[otu]]
        for sample in otudata:
            try:
                currentCount = abundancesDict[assigned_taxonomy][sample]
            except:
                currentCount = 0
            newCount = currentCount + otudata[sample][0]
            abundancesDict[assigned_taxonomy][sample] = newCount
    abundancesDF = pandas.DataFrame().from_dict(abundancesDict)
    abundancesDF = abundancesDF.transpose()

    return abundancesDF


def assignTaxonomy(threads, dbseq, directory):
    """
    This function assigns taxonomy using blastn
    """
    dbcmd = "makeblastdb -in " + dbseq + " -out " + directory + "analysis/diatoms -dbtype nucl"
    print dbcmd
    dbchild = subprocess.Popen(str(dbcmd),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True,
                               shell=(sys.platform != "win32"))
    dboutput, dberror = dbchild.communicate()
    print dboutput, dberror

    blastcmd = "blastn -db " + directory + "analysis/diatoms -query " + directory + "analysis/repset.fasta -out " + directory + "analysis/repset.diatoms.blastn -task blastn -max_target_seqs 1 -num_threads " + threads + " -outfmt 6 -evalue 0.01"
    print blastcmd
    blastchild = subprocess.Popen(str(blastcmd),
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  universal_newlines=True,
                                  shell=(sys.platform != "win32"))
    blastoutput, blasterror = blastchild.communicate()
    print blastoutput, blasterror

    outputfile = open(directory + "analysis/repset.taxonomy.txt", "w")
    for line in open(directory + "analysis/repset.diatoms.blastn", "rU"):
        linelist = re.split('\t', line)
        otu = linelist[0]
        percid = float(linelist[2])
        evalue = linelist[10]
        hit = linelist[1]
        if (percid >= 95.0):
            taxonomy = hit
            outputline = str(otu) + "\t" + str(taxonomy) + "\t" + str(evalue) + "\t" + str(percid) + "\n"
            outputfile.write(outputline)
        else:
            outputline = str(otu) + "\tNo_blast_hit;\tNone\tNone\n"
            outputfile.write(outputline)

    return True


def repsetCounts(clusteringFile, directory):
    """
    This function goes through the usearch results, determines the representative sequence
    and calculates the number of sequences from each sample per OTU.
    """
    ucfile = open(directory + "analysis/picked_otus_97/otu_clusters.uc", "rU")
    seeds = defaultdict(list)
    for line in ucfile:
        if line.startswith("S"):
            linesplit = line.split('\t')
            seed = linesplit[8].rstrip()
            seeds[seed].append(seed)
        if line.startswith("H"):
            linesplit = line.split('\t')
            seed = linesplit[9].rstrip()
            newseq = linesplit[8].rstrip()
            seeds[seed].append(newseq)
    otus = defaultdict(list)
    i = 0
    repset_names = {}
    testout = open(directory + "analysis/picked_otus_97/testseeds.txt","w")
    for seed in seeds:
        otuname = "denovo" + str(i)
        repset_names[otuname] = seed
        otus[otuname] = seeds[seed]
        testout.write(otuname + "," + seed + "\n")
        i = i + 1


    otu_names_file = open(directory + "analysis/picked_otus_97/otu_clusters.otus.txt", "w")
    countData = defaultdict(lambda: defaultdict(int))
    for otu in otus:
        otu_line = otu + "\t" + "\t".join(otus[otu]) + "\n"
        otu_names_file.write(otu_line)
        for seq in otus[otu]:
            seqsplit = re.split('_', seq)
            sample = seqsplit[0]
            sequence = seqsplit[1].rstrip()
            try:
                currentCount = countData[otu][sample]
            except:
                currentCount = 0
            newCount = currentCount + 1
            countData[otu][sample] = int(newCount)
    dataframe = pandas.DataFrame.from_dict(countData)
    dataframe = dataframe.transpose()
    dataframe = dataframe.fillna(0)
    dataframe.to_csv(directory + "analysis/picked_otus_97/otu_clusters.counts.txt", sep="\t")

    print "Retrieving the representative sequences..."
    repset_file = open(directory + "analysis/repset.fasta", "w")
    allsequences = SeqIO.index(directory + clusteringFile,"fasta")
    for otu, seedseq in repset_names.iteritems():
        if seedseq in allsequences.keys():
            seq = allsequences[seedseq]
            seq.id = otu
            seq.description = ""
            SeqIO.write(seq, repset_file, "fasta")
    return dataframe


def clustering(fastafile, directory):
    fastafile = fastafile.replace("analysis/", "")
    sortcmd = "uclustq1.2.22_i86linux32 --sort " + directory + "analysis/" + str(
        fastafile) + " --output " + directory + "analysis/" + str(fastafile) + ".sorted"
    sortchild = subprocess.Popen(str(sortcmd),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform != "win32"))
    sortoutput, sorterror = sortchild.communicate()
    print sortoutput, sorterror

    try:
        makeDir = os.mkdir(directory + "analysis/picked_otus_97")
    except OSError:
        pass
    usearchcmd = "uclustq1.2.22_i86linux32 --input " + directory + "analysis/" + str(
        fastafile) + ".sorted --id 0.97 --w 8 --tmpdir /tmp --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc " + directory + "analysis/picked_otus_97/otu_clusters.uc"
    usearchchild = subprocess.Popen(str(usearchcmd),
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True,
                                    shell=(sys.platform != "win32"))
    usearchoutput, usearcherror = usearchchild.communicate()
    print usearchoutput, usearcherror

    return True


def ampliconQC(dataDirectory, threads, forwardprimer, reverseprimer):
    """
    This function coordinates the amplicon quality control part of the pipeline
    """
    # Rename the files to teh format sampleName.readDirection.fastq.gz
    renamed = renameFiles(dataDirectory)

    # Move the 'undetermined' files to the raw_data folder
    moved = moveUndetermined(dataDirectory)

    # Trim off the primers from the sequences
    trimmed = runCutadapt(threads, dataDirectory, forwardprimer, reverseprimer)

    # Run sickle paired-end to trim off bad quality 3' bases
    peQC = runSicklePE(dataDirectory)

    # Run pear to merge the R1 and R2 reads
    joined = runPear(threads, dataDirectory)

    # Run sickle single-end to remove post-merging bad quality sequences
    seQC = runSickleSE(threads, dataDirectory)

    # Prepare the sequences for clustering
    clustering = clusteringPrep(threads, dataDirectory)

    # Draw sequence length histograms
    sequenceFiles = glob.glob(dataDirectory + "*.passedQC.fastq")
    for sequenceFile in sequenceFiles:
        graph = drawHistogram(sequenceFile)

    return "analysis/readyForClustering.allsamples.fasta"


def drawHistogram(sequenceFile):
    """
    Draws a histogram of the sequence length distribution of the sequences in the file.
    """
    filename, extension = os.path.splitext(sequenceFile)
    extension = re.sub('\.', '', extension)
    outfile = str(sequenceFile) + ".histogram.svg"
    sequences = [len(rec) for rec in SeqIO.parse(sequenceFile, extension)]
    counts = Counter(sequences)
    plot = pygal.XY(show_x_guides=True, show_legend=False, title="Sequence Length Histogram: " + sequenceFile,
                    x_title="Sequence length (nt)", y_title="Number of sequences")
    plot.add('Sequence Lengths', counts.items(), show_dots=False)
    plot.render_to_file(outfile)
    print "Histogram " + outfile + " finished"
    return True


def clusteringPrep(threads, directory):
    """
    Preparation of sequence files for clustering with usearch.
    """
    sequenceFiles = glob.glob(directory + "*.passedQC.fastq")

    sampleNames = []

    for fastqFile in sequenceFiles:
        filenameList = re.split('\.', fastqFile)
        originalSampleName = filenameList[0].replace(directory, '')
        newSampleName = originalSampleName.replace('-', '.')
        sampleNames.append(newSampleName)
        tempOutfile = open(directory + "temp." + newSampleName + ".fastq", "w")
        for seq in SeqIO.parse(fastqFile, "fastq"):
            name = seq.id
            seq.description = ""
            seq.id = newSampleName
            SeqIO.write(seq, tempOutfile, "fastq")
        tempOutfile.close()

        rename = os.rename(fastqFile, directory + newSampleName + ".passedQC.fastq")

    amendedSequenceFiles = glob.glob(directory + "temp.*")
    allseqs = open(directory + "allpassedQC.unnumbered.fastq", "w")
    for seqFile in amendedSequenceFiles:
        for line in open(seqFile, "rU"):
            allseqs.write(line)
    allseqs.close()

    for sample in sampleNames:
        tempdel = os.remove(directory + "temp." + sample + ".fastq")

    os.mkdir(directory + "analysis")
    outFasta = open(directory + "analysis/readyForClustering.allsamples.fasta", "w")
    i = 1
    for seq in SeqIO.parse(directory + "allpassedQC.unnumbered.fastq", "fastq"):
        name = seq.id
        seq.description = ""
        seq.id = name + "_" + str(i)
        SeqIO.write(seq, outFasta, "fasta")
        i = i + 1
    outFasta.close()
    print i, " sequences ready for clustering\n"

    deleted = os.remove(directory + "allpassedQC.unnumbered.fastq")

    graph = drawHistogram(directory + "analysis/readyForClustering.allsamples.fasta")

    return True


def runSickleSE(threads, directory):
    minLength = calculateMeanLength(directory) - 30  # forward primer that was trimmed off
    sickleSEcommand = "ls " + directory + "*.assembled.fastq | parallel -j " + str(
        threads) + " '/usr/local/bin/sickle se -f {} -t sanger -o {}.passedQC.fastq -n -l " + str(minLength) + " -q 30'"
    sickleSEchild = subprocess.Popen(str(sickleSEcommand),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True,
                                     shell=(sys.platform != "win32"))
    sickleSEoutput, sickleSEerror = sickleSEchild.communicate()
    print sickleSEoutput, sickleSEerror

    sequenceFiles = glob.glob(directory + "*.fastq")
    samples = []
    for fastqFile in sequenceFiles:
        filenameList = re.split('\.', fastqFile)
        samples.append(filenameList[0])
    sampleNames = list(set(samples))
    for sampleName in sampleNames:
        rm1 = os.remove(sampleName + ".assembled.fastq")
        rename = os.rename(sampleName + ".assembled.fastq.passedQC.fastq", sampleName + ".passedQC.fastq")
    return True


def calculateMeanLength(directory):
    mean = None
    means = []
    for fastqfile in glob.glob(directory + "*.assembled.fastq"):
        awk = "awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' " + fastqfile
        awkChild = subprocess.Popen(str(awk),
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True,
                                    shell=(sys.platform != "win32"))
        awkOutput, awkError = awkChild.communicate()
        output = awkOutput.rstrip()
        thismean = float(output)
        means.append(thismean)
    mean = sum(means) / float(len(means))
    return mean


def runPear(threads, directory):
    sequenceFiles = glob.glob(directory + "*.fastq.gz")
    samples = []
    for fastqFile in sequenceFiles:
        filenameList = re.split('\.', fastqFile)
        samples.append(filenameList[0].replace(directory, ''))
    sampleNames = list(set(samples))

    for sampleName in sampleNames:
        sampleName = directory + sampleName
        pearCommand = "/usr/local/bin/pear -e 1 -f " + sampleName + ".sickle.trimmed.R1.fastq.gz -r " + sampleName + ".sickle.trimmed.R2.fastq.gz -o " + sampleName + " -j " + str(
            threads)
        pearChild = subprocess.Popen(str(pearCommand),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True,
                                     shell=(sys.platform != "win32"))
        pearOutput, pearError = pearChild.communicate()
        print pearOutput, pearError
        rm1 = os.remove(sampleName + ".sickle.trimmed.R1.fastq.gz")
        rm2 = os.remove(sampleName + ".sickle.trimmed.R2.fastq.gz")
        rm3 = os.remove(sampleName + ".discarded.fastq")
        rm4 = os.remove(sampleName + ".unassembled.forward.fastq")
        rm5 = os.remove(sampleName + ".unassembled.reverse.fastq")
    return True


def runSicklePE(directory):
    sequenceFiles = glob.glob(directory + "*.fastq.gz")
    samples = []
    for fastqFile in sequenceFiles:
        filenameList = re.split('\.', fastqFile)
        samples.append(filenameList[0].replace(directory, ''))
    sampleNames = list(set(samples))

    for sampleName in sampleNames:
        sicklePEcommand = "/usr/local/bin/sickle pe -f " + directory + sampleName + ".R1.fastq.gz.trimmed.fastq.gz -r " + directory + sampleName + ".R2.fastq.gz.trimmed.fastq.gz -t sanger -o " + directory + sampleName + ".sickle.trimmed.R1.fastq.gz -p " + directory + sampleName + ".sickle.trimmed.R2.fastq.gz -s " + directory + sampleName + ".single.trimmed.fastq.gz -x"
        sicklePEchild = subprocess.Popen(str(sicklePEcommand),
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         universal_newlines=True,
                                         shell=(sys.platform != "win32"))
        sicklePEOutput, sicklePEError = sicklePEchild.communicate()
        print sicklePEOutput, sicklePEError

        delr1 = os.remove(directory + sampleName + ".R1.fastq.gz.trimmed.fastq.gz")
        delr2 = os.remove(directory + sampleName + ".R2.fastq.gz.trimmed.fastq.gz")
        delsingle = os.remove(directory + sampleName + ".single.trimmed.fastq.gz")
        move2 = os.rename(directory + sampleName + ".R1.fastq.gz",
                          directory + "raw_data/" + sampleName + ".R1.fastq.gz")
        move3 = os.rename(directory + sampleName + ".R2.fastq.gz",
                          directory + "raw_data/" + sampleName + ".R2.fastq.gz")
    return True


def runCutadapt(threads, directory, forwardprimer, reverseprimer):
    errorRate = 0.3
    cutadaptCommand = "ls " + directory + "*.fastq.gz | parallel -j " + str(threads) + " 'cutadapt -e " + str(
        errorRate) + " -b " + str(forwardprimer) + " -b " + str(reverseprimer) + " -o {}.trimmed.fastq.gz {}'"
    cutadaptChild = subprocess.Popen(str(cutadaptCommand),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True,
                                     shell=(sys.platform != "win32"))
    cutadaptOutput, cutadaptError = cutadaptChild.communicate()
    print cutadaptOutput, cutadaptError
    return True


def moveUndetermined(directory):
    try:
        makeDir = os.mkdir(directory + "raw_data")
    except OSError:
        pass
    try:
        moveR1 = os.rename(directory + "Undetermined.R1.fastq.gz", directory + "raw_data/Undetermined.R1.fastq.gz")
        moveR2 = os.rename(directory + "Undetermined.R2.fastq.gz", directory + "raw_data/Undetermined.R2.fastq.gz")
    except OSError:
        return True
    return True


def renameFiles(directory):
    files = glob.glob(directory + "*.fastq.gz")
    for filename in files:
        path, file = os.path.split(filename)
        filenameSplit = re.split('_', file)
        if (len(filenameSplit) == 1):
            return None
        else:
            sampleName = filenameSplit[0]
            readDirection = filenameSplit[3]
            newFilename = path + "/" + sampleName + "." + readDirection + ".fastq.gz"
            rename = os.rename(filename, newFilename)
            print "renamed:\t", filename, "\tto\t", newFilename
    return True


def processArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', help="Relative location of the data", required=True)
    parser.add_argument('--threads', help="Number of threads", default=4, required=False)
    parser.add_argument('--dbseq', help="Sequence database file", required=True)

    return parser.parse_args()


def dependencyChecks():
    if sys.version_info < (2, 7, 0):
        sys.stderr.write("!!ERROR!!\nYou need python 2.7 or later to run this script\n")
    if sys.version_info > (3, 0, 0):
        sys.stderr.write("!!ERROR!!\nYou need to run this script with python 2.7 rather than 3.0\n")

    sys.stdout.write("Checking the presence of dependencies...\n")

    if os.path.isfile("/usr/local/bin/uclustq1.2.22_i86linux32") is not True:
        sys.stderr.write("!!ERROR!!\nThe pipeline can't find usearch in /usr/local/bin.\n")
        sys.exit(1)
    else:
        sys.stdout.write("UCLUST has been found ---- PASS\n")

    if os.path.isfile("/usr/local/bin/sickle") is not True:
        sys.stderr.write("!!ERROR!!\nThe pipeline can't find sickle in /usr/local/bin.\n")
        sys.exit(1)
    else:
        sys.stdout.write("SICKLE has been found ---- PASS\n")

    if os.path.isfile("/usr/local/bin/pear") is not True:
        sys.stderr.write("!!ERROR!!\nThe pipeline can't find pear in /usr/local/bin.\n")
    else:
        sys.stdout.write("PEAR has been found ---- PASS\n")

    return True


if __name__ == "__main__":
    main()
