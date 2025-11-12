#!/usr/bin/env python3

import argparse
import sys
import multiprocessing as mp
import os
import shutil
import uuid
import pipettor
from flair.ssUtils import addOtherJuncs, gtfToSSBed
from flair.ssPrep import ssPrep
from flair import FlairInputDataError

import logging

def parseargs(aligned_reads=''):
    parser = argparse.ArgumentParser(description='take bed file of long RNA-seq reads and filter out those with anomalous splice junctions' \
                                     ' correct remaining to nearest orthogonally supported splice site')
    required = parser.add_argument_group('required named arguments')
    atleastone = parser.add_argument_group('at least one of the following arguments is required')
    if not aligned_reads:
        required.add_argument('-q', '--query', type=str, required=True,
                                                  help='uncorrected bed12 file')
    atleastone.add_argument('-f', '--gtf', default='',
                            help='GTF annotation file')
    mutexc = atleastone.add_mutually_exclusive_group(required=False)
    mutexc.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
                                               'Use this option if you aligned your short-reads with STAR, '
                                               'STAR will automatically output this file')
    mutexc.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from short-read alignment with junctions_from_sam)')
    parser.add_argument('--junction_support', type = int, default=1,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support. Default=1')
    parser.add_argument('-o', '--output', default='flair',
                                            help='output name base (default: flair)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                                            help='number of threads (4)')
    parser.add_argument('--nvrna', action='store_true', default=False,
                                            help='''specify this flag to make the strand of a read consistent with the annotation during correction''')
    parser.add_argument('-w', '--ss_window', type=int, default=15,
                                            help='window size for correcting splice sites (15)')
    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        parser.error('No arguments passed')

    args = parser.parse_args()

    if not (args.junction_tab or args.junction_bed):
        logging.warn('No short-read junctions provided. NO NOVEL SPLICE SITES WILL BE DETECTED.')
    return args

def correct(aligned_reads='', args=None):
    if not args: args = parseargs(aligned_reads)

    if aligned_reads:
        query = aligned_reads
    else:
        query = args.query

# TODO:This seems opposite the intended use, see what happens
    resolveStrand = False
    if not args.nvrna:
        resolveStrand = True

    # make temp dir for dumping
    tempDirName = str(uuid.uuid4())
    try:
        current_directory = os.getcwd()
        tempDir = os.path.join(current_directory, tempDirName)
        os.mkdir(tempDir)
    except OSError:
        raise OSError("Creation of the directory %s failed" % tempDirName)

    # There are a few functions that evaluate what verbose is defined as.
    # Instead of passing it around, just global it.
    global verbose
    global printErrFname
    global printErr
    verbose  = False # TODO
    printErr = False
    printErrFname = False
    if printErr:
        printErrFname = os.path.join(tempDirName, 'ssCorrect.err')

    # Convert gtf to bed and split by chromosome.
    juncs, chromosomes, knownSS  = dict(), set(), dict() # initialize juncs for adding to db
    if args.gtf:
        juncs, chromosomes, knownSS = gtfToSSBed(args.gtf, knownSS, printErr, printErrFname, verbose)

    # Do the same for the other juncs file.
    if args.junction_tab or args.junction_bed:
        if args.junction_tab:
            shortread, type = args.junction_tab, 'tab'
        else:
            shortread, type = args.junction_bed, 'bed'

        juncs, chromosomes, addFlag = addOtherJuncs(juncs, type, shortread, args.junction_support, chromosomes,
                printErrFname, knownSS, verbose, printErr)
        if addFlag == False:
            logging.info(f'ERROR Added no extra junctions from {shortread}\n')
    knownSS = dict()

    # added to allow annotations not to be used.
    if len(list(juncs.keys())) < 1:
        raise FlairInputDataError("No junctions from GTF or junctionsBed to correct with. Exiting...")

    annotations = dict()
    for chrom, data in juncs.items():
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
    notfound = False
    prevchrom = False
    tempoutfile = None
    with open(query) as lines, open("%s_cannot_verify.bed" % args.output,'w') as nochrom:
        outDict = dict()
        for line in lines:
            cols  = line.rstrip().split()
            chrom = cols[0]
            if chrom not in chromosomes:
                notfound = True
                nochrom.write(line)
                if chrom not in skippedChroms:
                    skippedChroms.add(chrom)
                    continue
            else:
                if chrom not in outDict:
                    readDict[chrom] = os.path.join(tempDir,"%s_temp_reads.bed" % chrom)
                    outDict[chrom] = os.path.join(tempDir,"%s_temp_reads.bed" % chrom)
                    #outDict[chrom] = open(os.path.join(tempDir,"%s_temp_reads.bed" % chrom),'w')
                if chrom != prevchrom:
                    if prevchrom is not False:
                        tempoutfile.close()
                    tempoutfile = open(outDict[chrom], 'a')
                    prevchrom = chrom
                print(line.rstrip(),file=tempoutfile)
                #with open(outDict[chrom], 'a') as tempoutfile:
                #print(line.rstrip(),file=outDict[chrom])
    nochrom.close()
    if tempoutfile:
        tempoutfile.close()

    if notfound is False:
        os.remove(f'{args.output}_cannot_verify.bed')

    cmds = list()
    for chrom in readDict:
        juncs = annotations[chrom]
        reads = readDict[chrom]

#               outDict[chrom].close()

        cmds.append([reads, juncs, args.ss_window, chrom, resolveStrand,
      tempDir, printErrFname])

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Prepared correct commands for %s read files" % len(cmds), file=fo)

    juncs = None
    annotations = None
    mp.set_start_method('fork')
    p = mp.Pool(args.threads)
    childErrs = set()
    for i in p.imap(ssPrep,cmds):
        childErrs.add(i)
    p.close()
    p.join()
    if len(childErrs) > 1:
        raise Exception(childErrs)

    with open("%s_all_inconsistent.bed" % args.output,'wb') as inconsistent:
        for chrom in readDict:
            with open(os.path.join(tempDir, "%s_inconsistent.bed" % chrom),'rb') as fd:
                shutil.copyfileobj(fd, inconsistent, 1024*1024*10)

    correct_bed = args.output + '_all_corrected.bed'
    with open(correct_bed,'wb') as corrected:
        for chrom in readDict:
            with open(os.path.join(tempDir, "%s_corrected.bed" % chrom),'rb') as fd:
                shutil.copyfileobj(fd, corrected, 1024*1024*10)
    pipettor.run([('bedtools', 'sort', '-i', correct_bed)], stdout=args.output + '_all_corrected.sorted.bed')
    shutil.move(args.output + '_all_corrected.sorted.bed', correct_bed)
    pipettor.run([('bedtools', 'sort', '-i', args.output + '_all_inconsistent.bed')], stdout=args.output + '_all_inconsistent.sorted.bed')
    shutil.move(args.output + '_all_inconsistent.sorted.bed', args.output + '_all_inconsistent.bed')
    if printErr:
        shutil.move(printErrFname, f'{args.output}.err')
    shutil.rmtree(tempDir)

    return correct_bed

if __name__ == '__main__':
    correct()
