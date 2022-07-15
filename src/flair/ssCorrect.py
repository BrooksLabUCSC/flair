#!/usr/bin/env python3
from __future__ import print_function

########################################################################
# File: ssCorrect.py
#  executable: ssCorrect.py
# Purpose:
#
#
# Author: Cameron M. Soulette
# History:      cms 05/01/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
from multiprocessing import Pool
from tqdm import *
import subprocess
import shutil
import uuid
import pybedtools
import re

scriptPath = os.path.realpath(__file__)
path = "/".join(scriptPath.split("/")[:-1])
helperScript = path + "/" + "ssPrep.py"

########################################################################
# CommandLine
########################################################################

class CommandLine(object) :
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    and a standard usage and help,
    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']

    methods:

    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = ' ssCorrect.py - a tool to leverage annotation and short read data to correct misaligned splice junctions in short read data.',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s -i reads.bed -g annotations.gtf -j other_junctions.bed -o out_file.bed')
        # Add args
        self.parser.add_argument('-i', '--input_bed', action = 'store', required=True, help='Input reads in bed12 format.')
        self.parser.add_argument('-g', '--gtf', action = 'store', required=False, help='Gencode annotation file.')
        self.parser.add_argument('-j', '--junctionsBed', default=None, action = 'store', required=False, help='Short-read supported junctions in bed6 format (Optiona) [BED entries must be UNIQUE and have strand information].')
        self.parser.add_argument('-w', '--wiggleWindow', action = 'store', type=int, required=False, default = 15, help='Splice site correction window flank size.')
        self.parser.add_argument('-o', '--output_fname', action = 'store', required=True, help='Output file name.')
        self.parser.add_argument('-f', '--genome_fasta', action = 'store', required=True, help='Bedtools indexed genome fasta.')

        self.parser.add_argument('--correctStrand', action = 'store_true', required=False, default = False, help='Try to resolve read strand by using annotated splice site strand.')
        self.parser.add_argument('-p', '--threads', action = 'store', required=False, type=int, default = 2, help='Number of threads.')
        self.parser.add_argument('--progress', action = 'store_true', required=False, default = False, help='Display progress')
        self.parser.add_argument('--tempDir', action = 'store', required=False, default = None,   help='Output directory for temporary files.')
        self.parser.add_argument('--keepTemp', action = 'store_true', required=False, default = False, help='Keep temporary/intermediate files.')
        self.parser.add_argument('--print_check', action = 'store_true', required=False, default = False, help='Print workflow checking')

        #self.parser.add_argument('--keepZero', action = 'store_true', required=False, default = False, help='Keep alignments with no spliced junctions (single exon txns).')
        #self.parser.add_argument("--quiet", action = 'store_false', required=False, default = True, help='Do not display progress')

        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))


########################################################################
# Functions
########################################################################


def addOtherJuncs(juncs, bedJuncs, chromosomes, fa, printErrFname, known):

    lineNum = 0
    if verbose: sys.stderr.write("Step 2/5: Processing additional junction file  %s ..." % (bedJuncs))
    cols = None

    with open(bedJuncs) as l:
        for num,ll in enumerate(l,0):
            cols = ll.rstrip().split()
            if num >10:
                break


    # guess what kind of bedFile
    if cols is None:
        raise Exception("Empty junctions BED file, not supported")

    if cols[-1] == "+" or cols[-1] == "-":
        # normal bed
        reverseSS = "-"
        strandCol = -1
        starOffset = 0

    elif len(cols) == 12:
        # bed12
        bedType   = "bed12"
        raise Exception("Bed12 not currently supported for other_juncs.bed. Please convert to bed6.")

    elif cols[3] == "0" or cols[3] == "1" or cols[3] == "2":
        # star junc.tab
        reverseSS = "2"
        strandCol = 3
        starOffset = 1

    else:
        raise Exception("Cannot find strand info for %s. Is this bed6 or STAR_juncs.tab file?" % bedJuncs)


    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Adding other juncs, assuming file is %s" % "bed6" if strandCol == -1 else "STAR", file=fo)

    tempJuncs = list()
    with open(bedJuncs,'r') as bedLines:
        for line in bedLines:
            cols = line.rstrip().split()
            chrom, c1, c2, strand = cols[0], int(cols[1])-starOffset, int(cols[2]), cols[strandCol]

            if chrom not in juncs:
                juncs[chrom] = dict()

            if c2-c1 < 5:
                continue

            if starOffset:
                if strand == "1": strand = "+"
                elif strand == "2": strand = "-"
                else: continue

            chromosomes.add(chrom)
            key = (c1, c2, strand)
            if key in juncs[chrom]:
                juncs[chrom][key] = "both"
                continue
            tempJuncs.append((chrom,c1,c2,"%s,%s,%s,%s" % (chrom,c1,c2,strand),0,strand))

    try:
        btJuncs = pybedtools.BedTool(tempJuncs)
        dinucSeq = btJuncs.sequence(fi=fa, s=True, tab=True, name=True)
        with open(dinucSeq.seqfn) as fileObj:
            for i in fileObj:
                header,seq = i.rstrip().split()
                chrom,c1,c2,strand = header.split(",")
                c1,c2 = int(c1),int(c2)
                if "+" in strand:
                    strand = strand[strand.rfind('+')]
                elif "-" in strand:
                    strand = strand[strand.rfind('-')]
                key = (c1,c2, strand)
                known1,known2 = known.get((chrom,c1),None),known.get((chrom,c2),None)

                if known1 != None:
                    if known1 != strand:
                        continue
                    else:
                        pass
                else:
                    pass

                if known2 != None:
                    if known2 != strand:
                        continue
                    else:
                        pass
                else:
                    pass

                fivePrime = seq[:2]
                if key not in juncs[chrom]:
                    juncs[chrom][key] = "sr"
                elif fivePrime == "GT":
                    juncs[chrom][key] = "sr"

    except Exception as e:
        print(e,"Splice site motif filtering failed. Check pybedtools and bedtools is properly install and in $PATH",file=sys.stderr)
        sys.exit(1)

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** GTF Juncs + other juncs now total %s juncs from %s chromosomes." % (sum([len(x)for x in juncs.values()]), len(list(juncs.keys()))), file=fo)

    return juncs, chromosomes

def gtfToSSBed(file, knownSS, printErr, printErrFname, verbose):
    ''' Convenience function, reformats GTF to bed'''

    # First: get all exons per transcript.
    exons = dict()
    chromosomes = set()
    with open(file,'r') as lines:
        for l in lines:
            if l[0] == "#": # skip header lines
                continue

            cols = l.split("\t")

            if "exon" == cols[2]:

                # -1 for 1 to 0 based conversion
                chrom, c1, c2, strand =  cols[0], int(cols[3])-1, int(cols[4]), cols[6]
                chromosomes.add(chrom)
                #txn info is in the SECOND position of the shoutout column
                txn = re.search('transcript_id "([^\"]+)"', l).group(1)#
                #cols[-1].split(";")[1].split()[-1].replace('"','')

                key = (chrom, txn, strand)

                if key not in exons:
                    exons[key] = list()
                exons[key].append(c1)
                exons[key].append(c2)

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Read GTF. Got %s transcripts" % len(list(exons.keys())), file=fo)
            print("** Getting introns...Read GTF", file=fo)

    # Second: get junction and splice sites from transcript exons.
    txnList = list(exons.keys())
    juncs = dict()

    for exonInfo in tqdm(txnList, total=len(txnList), desc="Step 1/5: Splitting junctions from GTF by chromosome", dynamic_ncols=True, position=1) if verbose else txnList:
        chrom, txn, strand = exonInfo

        if chrom not in juncs:
            juncs[chrom] = dict()

        coords = list(exons[exonInfo])

        # assume lowest and highest as TSS and TES, and remove them
        coords.sort()
        coords = coords[1:-1]

        # Coords is list of exons, so a list less than 2 is a single exon gene.
        if len(coords)<2: continue

        for pos in range(0,len(coords)-1,2):
            c1 = coords[pos]
            c2 = coords[pos+1]

            if abs(c2 - c1) <= 5:
                continue

            juncs[chrom][(c1,c2,strand)] = "gtf"
            knownSS[(chrom, c1)] = strand
            knownSS[(chrom, c2)] = strand

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Created %s juncs from %s chromosomes." % (sum([len(x)for x in juncs.values()]), len(list(juncs.keys()))), file=fo)
    return juncs, chromosomes, knownSS


def runCMD(x):

    tDir, prefix,juncs,reads, rs, f, err, errFname = x
    cmd = "%s %s -i %s -j %s -o %s --workingDir %s -f %s " % (sys.executable, helperScript, reads,juncs,prefix, tDir, f)
    if rs:
        cmd += "--correctStrand "
    if err:
        cmd += "--check_file %s " % errFname
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std, perr = p.communicate()
    return perr

def main():
    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()
    bed           = myCommandLine.args['input_bed']
    gtf           = myCommandLine.args['gtf']
    otherJuncs    = myCommandLine.args['junctionsBed']
    wiggle        = myCommandLine.args['wiggleWindow']
    threads       = myCommandLine.args['threads']
    outFile       = myCommandLine.args['output_fname']
    keepTemp      = myCommandLine.args['keepTemp']
    resolveStrand = myCommandLine.args['correctStrand']
    tempDirName   = myCommandLine.args['tempDir']
    genomeFasta   = myCommandLine.args['genome_fasta']

    if os.path.isfile(genomeFasta+".fai"):
        pass
    else:
        testString =  """
            chrX 1    100   feature1  0 +
        """
        test = pybedtools.BedTool(testString, from_string=True)
        a = test.sequence(fi=genomeFasta)


    # make temp dir for dumping
    if tempDirName == None:
        tempDirName = "tmp_%s" % str(uuid.uuid4())
    try:
        current_directory = os.getcwd()
        tempDir = os.path.join(current_directory, tempDirName)
        os.mkdir(tempDir)
    except OSError:
        print ("Creation of the directory %s failed" % tempDirName)
        sys.exit(1)

    # There are a few functions that evaluate what verbose is defined as.
    # Instead of passing it around, just global it.
    global verbose
    global printErrFname
    global printErr
    verbose  = myCommandLine.args['progress']
    printErr = myCommandLine.args['print_check']
    printErrFname = "err_%s.txt" % tempDirName

    # Convert gtf to bed and split by cromosome.
    juncs, chromosomes, knownSS  = dict(), set(), dict() # initialize juncs for adding to db
    if gtf != None: juncs, chromosomes, knownSS = gtfToSSBed(gtf, knownSS, printErr, printErrFname, verbose)

    # Do the same for the other juncs file.
    if otherJuncs != None: juncs, chromosomes = addOtherJuncs(juncs, otherJuncs, chromosomes, genomeFasta, printErrFname, knownSS)
    knownSS = dict()

    # added to allow annotations not to be used.
    if len(list(juncs.keys()))<1:
        print("No junctions from GTF or junctionsBed to correct with. Exiting...", file=sys.stderr)
        sys.exit(1)

    annotations = dict()
    for chrom, data in tqdm(juncs.items(), desc="Step 3/5: Preparing annotated junctions to use for correction", total=len(list(juncs.keys())), dynamic_ncols=True, position=1) if verbose else juncs.items():
        annotations[chrom] = os.path.join(tempDir,"%s_known_juncs.bed" % chrom)
        with open(os.path.join(tempDir,"%s_known_juncs.bed" % chrom),"w") as out:
            sortedData = sorted(list(data.keys()), key=lambda item: item[0])
            for k in sortedData:
                annotation = data[k]
                c1, c2, strand = k
                print(chrom,c1,c2,annotation,".",strand, sep="\t", file=out)

    sortedData = None
    skippedChroms = set()
    readDict = dict()
    with open(bed) as lines:
        outDict = dict()
        for line in tqdm(lines, desc="Step 4/5: Preparing reads for correction", dynamic_ncols=True, position=1) if verbose else lines:
            cols  = line.rstrip().split()
            chrom = cols[0]
            if chrom not in chromosomes:
                if chrom not in skippedChroms:
                    skippedChroms.add(chrom)
                    #if verbose: tqdm.write("Reference sequence not found in annotations, skipping: %s" % (chrom), file=sys.stderr)
                    continue
            else:
                if chrom not in outDict:
                    readDict[chrom] = os.path.join(tempDir,"%s_temp_reads.bed" % chrom)
                    outDict[chrom] = open(os.path.join(tempDir,"%s_temp_reads.bed" % chrom),'w')
                print(line.rstrip(),file=outDict[chrom])

    cmds = list()
    for chrom in readDict:
        juncs = annotations[chrom]
        reads = readDict[chrom]

        outDict[chrom].close()

        cmds.append((tempDir, chrom, juncs,reads, resolveStrand, genomeFasta, printErr, printErrFname))

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Prepared correct commands for %s read files" % len(cmds), file=fo)


    juncs = None
    annotations = None
    p = Pool(threads)
    childErrs = set()
    for i in tqdm(p.imap(runCMD, cmds), total=len(cmds), desc="Step 5/5: Correcting Splice Sites", dynamic_ncols=True,position=1) if verbose else p.imap(runCMD,cmds):
        childErrs.add(i)
    if len(childErrs)>1:
        print(childErrs,file=sys.stderr)
        sys.exit(1)


    with open("%s_all_inconsistent.bed" % outFile,'wb') as inconsistent:
        for chrom in readDict:
            with open(os.path.join(tempDir, "%s_inconsistent.bed" % chrom),'rb') as fd:
                shutil.copyfileobj(fd, inconsistent, 1024*1024*10)


    with open("%s_all_corrected.bed" % outFile,'wb') as corrected:
        for chrom in readDict:
            with open(os.path.join(tempDir, "%s_corrected.bed" % chrom),'rb') as fd:
                shutil.copyfileobj(fd, corrected, 1024*1024*10)
    if keepTemp:
        pass
    else:
        try:
            shutil.rmtree(tempDir)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror), file=sys.stderr)

    print("\n")
if __name__ == "__main__":
    main()
