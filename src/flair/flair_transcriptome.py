#! /usr/bin/env python3

import sys
import argparse
import os
import pipettor
import glob
import uuid
import shutil
import pysam
import logging
from flair.flair_align import inferMM2JuncStrand, intronChainToestarts
from flair.ssUtils import addOtherJuncs, gtfToSSBed
from flair.ssPrep import buildIntervalTree, ssCorrect, juncsToBed12
import multiprocessing as mp
import time
from collections import Counter
from flair import FlairInputDataError


# export PATH="/private/groups/brookslab/cafelton/git-flair/flair/bin:/private/groups/brookslab/cafelton/git-flair/flair/src/flair:$PATH"

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--genomealignedbam',
                        help='Sorted and indexed bam file aligned to the genome')
    parser.add_argument('-g', '--genome', type=str,
                        help='FastA of reference genome, can be minimap2 indexed')
    parser.add_argument('-o', '--output', default='flair',
                        help='output file name base for FLAIR isoforms (default: flair)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (4)')
    parser.add_argument('-f', '--gtf', default='',
                        help='GTF annotation file, used for renaming FLAIR isoforms '
                             'to annotated isoforms and adjusting TSS/TESs')

    mutexc = parser.add_mutually_exclusive_group(required=False)
    mutexc.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
                                               'Use this option if you aligned your short-reads with STAR, '
                                               'STAR will automatically output this file')
    mutexc.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from short-read alignment with junctions_from_sam)')
    parser.add_argument('--junction_support', type=int, default=1,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support.')
    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (15)')

    parser.add_argument('-s', '--support', type=float, default=3.0,
                        help='''minimum number of supporting reads for an isoform (3)''')
    parser.add_argument('--stringent', default=False, action='store_true',
                        help='''specify if all supporting reads need to be full-length
                (spanning 25 bp of the first and last exons)''')
    parser.add_argument('--check_splice', default=False, action='store_true',
                        help='''enforce coverage of 4 out of 6 bp around each splice site and no
                insertions greater than 3 bp at the splice site''')
    # parser.add_argument('--quality', type=int, default=0,
    #                     help='minimum MAPQ of read assignment to an isoform (1)')
    parser.add_argument('-w', '--end_window', type=int, default=100,
                        help='window size for comparing TSS/TES (100)')

    parser.add_argument('--noaligntoannot', default=False, action='store_true',
                        help='''related to old annotation_reliant, now specify if you don't want
                        an initial alignment to the annotated sequences and only want transcript
                        detection from the genomic alignment.
                         Will be slightly faster but less accurate if the annotation is good''')
    parser.add_argument('-n', '--no_redundant', default='none',
                        help='''For each unique splice junction chain, report options include:
            none--best TSSs/TESs chosen for each unique set of splice junctions;
            longest--single TSS/TES chosen to maximize length;
            best_only--single most supported TSS/TES used in conjunction chosen (none)''')
    parser.add_argument('--max_ends', type=int, default=2,
                        help='maximum number of TSS/TES picked per isoform (2)')
    parser.add_argument('--filter', default='default',
                        help='''Report options include:
            nosubset--any isoforms that are a proper set of another isoform are removed;
            default--subset isoforms are removed based on support;
            comprehensive--default set + all subset isoforms;
            ginormous--comprehensive set + single exon subset isoforms''')

    parser.add_argument('--parallelmode', default='auto:1GB',
                        help='''parallelization mode. Default: "auto:1GB" This indicates an automatic threshold where
                            if the file is less than 1GB, parallelization is done by chromosome, but if it's larger,
                            parallelization is done by region of non-overlapping reads. Other modes: bychrom, byregion,
                            auto:xGB - for setting the auto threshold, it must be in units of GB.''')
    parser.add_argument('--predictCDS', default=False, action='store_true',
                        help='specify if you want to predict the CDS of the final isoforms. '
                             'Will be output in the final bed file but not the gtf file. '
                             'Productivity annotation is also added in the name field, '
                             'which is detailed further in the predictProductivity documentation')
    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''specify if intermediate and temporary files are to be kept for debugging.
                Intermediate files include: promoter-supported reads file,
                read assignments to firstpass isoforms''')
    parser.add_argument('--keep_sup', default=False, action='store_true',
                        help='''specify if you want to keep supplementary alignments to define isoforms''')

    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        parser.error('No arguments passed. Please provide bam file, genome, and gtf')

    args, unknown = parser.parse_known_args()
    if unknown:
        logging.info(f'unrecognized arguments: {" ".join(unknown)}')

    check_file_paths(args)
    args = add_preset_args(args)
    # args.quality = '0' if args.trust_ends else args.quality
    # if args.mm2_args:
    #     args.mm2_args = args.mm2_args.split(',')
    return args


def check_file_paths(args):
    if not args.genomealignedbam:
        raise FlairInputDataError(f'Please include the --genomealignedbam option')
    if not args.genome:
        raise FlairInputDataError(f'Please include the --genome option\n')
    if not os.path.exists(args.genomealignedbam):
        raise FlairInputDataError(f'Aligned reads file path does not exist: {args.genomealignedbam}')
    if not os.path.exists(args.genome):
        raise FlairInputDataError(f'Genome file path does not exist: {args.genome}')
    if not (args.parallelmode in {'bychrom', 'byregion'}
            or (args.parallelmode[:5] == 'auto:'
                and ((args.parallelmode[-2:] == 'GB' and args.parallelmode[5:-2].replace(".", "").isnumeric())
                     or args.parallelmode[5:].replace(".", "").isnumeric()))):
        raise FlairInputDataError(
            f'parallelmode {args.parallelmode} is not in an allowed format. See docs for allowed formats')


def add_preset_args(args):
    args.mm2_args = []
    args.quality = 0
    args.trust_ends = False
    args.remove_internal_priming = False
    args.isoformtss = True
    return args


def makecorrecttempdir():
    tempDirName = str(uuid.uuid4())
    try:
        current_directory = os.getcwd()
        tempDir = os.path.join(current_directory, tempDirName)
        os.mkdir(tempDir)
    except OSError:
        raise OSError(f"Creation of the directory {tempDirName} failed")
    return tempDir + '/'


def bin_search(query, data):
    """ Query is a coordinate interval. Binary search for the query in sorted data,
        which is a list of coordinates. Finishes when an overlapping value of query and
        data exists and returns the index in data. """
    i = int(round(len(data) / 2))  # binary search prep
    lower, upper = 0, len(data)
    while True:
        if upper - lower < 2:  # stop condition but not necessarily found
            break
        if data[i][1] < query[0]:
            lower = i
            i = int(round((i + upper) / 2))
        elif data[i][0] > query[1]:
            upper = i
            i = int(round((lower + i) / 2))
        else:  # found
            break
    return i


def generateKnownSSDatabase(args, tempDir):
    # Convert gtf to bed and split by chromosome.
    juncs, chromosomes, knownSS = dict(), set(), dict()  # initialize juncs for adding to db

    if args.gtf:
        juncs, chromosomes, knownSS = gtfToSSBed(args.gtf, knownSS, False, False, False)

    # Do the same for the other juncs file.
    if args.junction_tab or args.junction_bed:
        if args.junction_tab:
            shortread, type = args.junction_tab, 'tab'
        else:
            shortread, type = args.junction_bed, 'bed'
        juncs, chromosomes, addFlag = addOtherJuncs(juncs, type, shortread, args.junction_support, chromosomes,
                                                    False, knownSS, False, False)
        if addFlag == False:
            logging.info(f'ERROR Added no extra junctions from {shortread}')

    # added to allow annotations not to be used.
    if len(list(juncs.keys())) < 1:
        raise FlairInputDataError("No junctions from GTF or junctionsBed to correct with")

    annotationFiles = dict()
    for chrom in chromosomes:
        annotationFiles[chrom] = os.path.join(tempDir, "%s_known_juncs.bed" % chrom)
        with open(os.path.join(tempDir, "%s_known_juncs.bed" % chrom), "w") as out:
            if chrom in juncs:
                data = juncs[chrom]
                sortedData = sorted(list(data.keys()), key=lambda item: item[0])
                for k in sortedData:
                    annotation = data[k]
                    c1, c2, strand = k
                    print(chrom, c1, c2, annotation, ".", strand, sep="\t", file=out)
    return chromosomes, annotationFiles


def correctsingleread(bedread, intervalTree, junctionBoundaryDict):
    # FIXME: this was copied from ssPrep.py::correctReads() and modified
    juncs = bedread.juncs
    strand = bedread.strand
    c1Type, c2Type = ("donor", "acceptor") if strand == "+" else ("acceptor", "donor")

    newJuncs = list()
    ssStrands = set()

    for x in juncs:
        c1, c2 = x[0], x[1]
        if c1 not in junctionBoundaryDict:
            junctionBoundaryDict = ssCorrect(c1, strand, c1Type, intervalTree, junctionBoundaryDict, False)
        if c2 not in junctionBoundaryDict:
            junctionBoundaryDict = ssCorrect(c2, strand, c2Type, intervalTree, junctionBoundaryDict, False)

        c1Corr = junctionBoundaryDict[c1].ssCorr.coord
        c2Corr = junctionBoundaryDict[c2].ssCorr.coord
        # don't allow junctions outside or near the ends of the reads
        ends_slop = 8
        if not ((bedread.start + ends_slop) <= c1Corr < (bedread.end - ends_slop)):
            return None
        if not ((bedread.start + ends_slop) <= c2Corr < (bedread.end - ends_slop)):
            return None

        ssTypes = [junctionBoundaryDict[c1].ssCorr.ssType, junctionBoundaryDict[c2].ssCorr.ssType]

        ssStrands.add(junctionBoundaryDict[c1].ssCorr.strand)
        ssStrands.add(junctionBoundaryDict[c2].ssCorr.strand)

        if None in ssTypes:  # or ssTypes[0] == ssTypes[1]: # Either two donors or two acceptors or both none.
            return None
        newJuncs.append((c1Corr, c2Corr))

    blocks, sizes, starts = juncsToBed12(bedread.name, bedread.refchrom, bedread.start, bedread.end, newJuncs)
    if blocks is None:
        return None  # tmp until BED construction bugs are fixed

    # 0 length exons, remove them.
    if min(sizes) == 0:
        return None

    bedread.juncs = newJuncs
    bedread.esizes = sizes
    bedread.estarts = starts
    bedread.setexons()
    return bedread


def getrgb(name, strand, junclen):
    if name[:4] == 'ENST':
        return '3,28,252'
    elif junclen == 0:
        return "99,99,99"
    elif strand == '+':
        return "27,158,119"
    else:
        return "217,95,2"


def getexonsfromjuncs(juncs, start, end):
    if len(juncs) == 0:
        estarts = [0]
        esizes = [end - start]
    else:
        estarts = [0] + [x[1] - start for x in juncs]
        esizes = [juncs[0][0] - start] + [juncs[i + 1][0] - juncs[i][1] for i in range(len(juncs) - 1)] + [
            end - juncs[-1][1]]
    return estarts, esizes


class BedRead(object):
    def __init__(self):
        self.name = None

    def generate_from_cigar(self, alignstart, is_reverse, cigartuples, readname, referencename, qualscore,
                            juncDirection):
        # alignstart, is_reverse, cigartuples, readname, referencename, qualscore, juncDirection
        positiveTxn = "27,158,119"  # green
        negativeTxn = "217,95,2"  # orange
        unknownTxn = "99,99,99"
        refpos = alignstart
        intronblocks = []
        hasmatch = False
        for block in cigartuples:
            if block[0] == 3 and hasmatch:  # intron, pay attention
                intronblocks.append([refpos, refpos + block[1]])
                refpos += block[1]
            elif block[0] in {0, 7, 8, 2}:  # consumes reference
                refpos += block[1]
                if block[0] in {0, 7, 8}:
                    hasmatch = True  # match
        # dirtowrite = '-' if is_reverse else '+'
        # chr1  476363  497259  ENST00000455464.7_ENSG00000237094.12    1000    -       476363  497259  0       3       582,169,151,    0,8676,20745,
        esizes, estarts = intronChainToestarts(intronblocks, alignstart, refpos)
        if juncDirection not in {'+', '-'}:
            juncDirection = "-" if is_reverse else "+"

        junctions = []
        for i in range(len(estarts) - 1):
            junctions.append((alignstart + estarts[i] + esizes[i], alignstart + estarts[i + 1]))

        self.refchrom = referencename
        self.start = alignstart
        self.end = refpos
        self.name = readname
        self.score = min(qualscore, 1000)
        self.strand = juncDirection
        self.blockcount = len(intronblocks) + 1
        self.esizes = esizes
        self.estarts = estarts
        self.juncs = tuple(junctions)
        self.setexons()

    def setexons(self):
        self.exons = [(self.start + self.estarts[i], self.start + self.estarts[i] + self.esizes[i]) for i in
                      range(len(self.estarts))]

    def getsequence(self, genome):
        exons = self.exons
        if self.strand == '-':
            exons = exons[::-1]
        exonseq = []
        for i in range(len(exons)):
            thisexonseq = genome.fetch(self.refchrom, exons[i][0], exons[i][1])
            if self.strand == '-':
                thisexonseq = revcomp(thisexonseq)
            exonseq.append(thisexonseq)
        return ''.join(exonseq)

    def generatefromvals(self, chrom, start, end, name, score, strand, juncs):
        self.refchrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.juncs = juncs
        self.estarts, self.esizes = getexonsfromjuncs(juncs, start, end)
        self.setexons()

    def getbedline(self):
        rgbcolor = getrgb(self.name, self.strand, len(self.juncs))
        bedline = [self.refchrom, self.start, self.end, self.name, self.score, self.strand,
                   self.start, self.end, rgbcolor, len(self.estarts), ','.join([str(x) for x in self.esizes]),
                   ','.join([str(x) for x in self.estarts])]
        bedline = [str(x) for x in bedline]
        return bedline


class AnnotData(object):
    def __init__(self):
        self.transcripttoexons = {}
        self.transcripttoinfo = {}
        self.alltranscripts = []
        self.juncstotranscript = {}
        self.junctogene = {}
        self.allannotse = []
        self.genetoannotjuncs = {}
        self.genetostrand = {}

    def returndata(self):
        return self.juncstotranscript, self.junctogene, self.allannotse, self.genetoannotjuncs, self.genetostrand, self.transcripttoexons, self.alltranscripts

def generate_region_dict(allregions):
    chromtoregions, regionstoannotdata = {}, {}
    for chrom, rstart, rend in allregions:
        if chrom not in chromtoregions:
            chromtoregions[chrom] = []
        chromtoregions[chrom].append((rstart, rend))
        regionstoannotdata[(chrom, rstart, rend)] = AnnotData()
    return chromtoregions, regionstoannotdata

def get_tname_to_exons(gtf):
    allchromtotranscripttoexons = {}
    for line in open(gtf):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]
            if ty == 'transcript' or ty == 'exon':
                if chrom not in allchromtotranscripttoexons:
                    allchromtotranscripttoexons[chrom] = {}
                this_transcript = line[8][line[8].find('transcript_id') + 15:]
                this_transcript = this_transcript[:this_transcript.find('"')]
                this_gene = line[8].split('gene_id "')[1].split('"')[0].replace('_', '-')
                if (this_transcript, this_gene) not in allchromtotranscripttoexons[chrom]:
                    allchromtotranscripttoexons[chrom][(this_transcript, this_gene)] = [(None, None, strand), []]

                if ty == 'transcript':
                    allchromtotranscripttoexons[chrom][(this_transcript, this_gene)][0] = (start, end, strand)
                elif ty == 'exon':
                    allchromtotranscripttoexons[chrom][(this_transcript, this_gene)][1].append((start, end))
    return allchromtotranscripttoexons

def get_annot_tstart_tend(tinfo):
    tstart, tend, strand = tinfo[0]
    if tstart == None:
        tstart = min([x[0] for x in tinfo[1]])
        tend = max([x[1] for x in tinfo[1]])
    return tstart, tend, strand

def save_transcript_annot_to_region(transcript, gene, thisregion, regionstoannotdata, tstart, tend, strand, texons):
    sortedexons = sorted(texons)
    regionstoannotdata[thisregion].transcripttoexons[(transcript, gene)] = tuple(sortedexons)
    juncs = []
    for i in range(len(sortedexons) - 1):
        juncs.append((sortedexons[i][1], sortedexons[i + 1][0]))
    regionstoannotdata[thisregion].alltranscripts.append((transcript, gene, strand))
    if gene not in regionstoannotdata[thisregion].genetostrand:
        regionstoannotdata[thisregion].genetostrand[gene] = strand
    if len(juncs) == 0:
        regionstoannotdata[thisregion].allannotse.append((tstart, tend, strand, gene))
    else:
        regionstoannotdata[thisregion].juncstotranscript[tuple(juncs)] = (transcript, gene)
        if gene not in regionstoannotdata[thisregion].genetoannotjuncs:
            regionstoannotdata[thisregion].genetoannotjuncs[gene] = set()
        for j in juncs:
            if j not in regionstoannotdata[thisregion].junctogene:
                regionstoannotdata[thisregion].junctogene[j] = set()
            regionstoannotdata[thisregion].junctogene[j].add((transcript, gene))
            regionstoannotdata[thisregion].genetoannotjuncs[gene].add(j)
    regionstoannotdata[thisregion].allannotse = sorted(regionstoannotdata[thisregion].allannotse)
    return regionstoannotdata

def get_annot_for_chrom(chromregions, rchrom, regionstoannotdata, chrom_transcripttoexons):
    for transcript, gene in chrom_transcripttoexons:
        tinfo = chrom_transcripttoexons[(transcript, gene)]
        tstart, tend, strand = get_annot_tstart_tend(tinfo)
        for rstart, rend in chromregions:
            if rstart < tstart < rend or rstart < tend < rend:
                thisregion = (rchrom, rstart, rend)
                regionstoannotdata = save_transcript_annot_to_region(transcript, gene, thisregion, regionstoannotdata, tstart, tend, strand, tinfo[1])


    return regionstoannotdata


def getannotinfo(gtf, allregions):
    chromtoregions, regionstoannotdata = generate_region_dict(allregions)
    allchromtotranscripttoexons = get_tname_to_exons(gtf)

    for rchrom in allchromtotranscripttoexons:
        if rchrom in chromtoregions: # only get annot for regions that exist in reads
            regionstoannotdata = get_annot_for_chrom(chromtoregions[rchrom], rchrom, regionstoannotdata, allchromtotranscripttoexons[rchrom])
    return regionstoannotdata


def getcountsamcommand(args, refbed, outputname, mapfile, isannot):
    # count sam transcripts ; the dash at the end means STDIN
    count_cmd = ['filter_transcriptome_align.py', '--sam', '-',
                 '-o', outputname, '-t', 1, ###feeding 1 thread in because this is already multithreaded here
                 '--quality', str(args.quality), '-w', str(args.end_window)]
    if mapfile:
        count_cmd += ['--generate_map', mapfile]
    if args.stringent or isannot:
        count_cmd += ['--stringent']
    if args.check_splice:
        count_cmd += ['--check_splice']
    if args.check_splice or args.stringent or isannot:
        count_cmd += ['-i', refbed]  # annotated isoform bed file
    if args.trust_ends:
        count_cmd += ['--trust_ends']
    if args.remove_internal_priming:
        count_cmd += ['--remove_internal_priming', '--intprimingthreshold', str(args.intprimingthreshold),
                      '--intprimingfracAs', str(args.intprimingfracAs), '--transcriptomefasta', args.transcriptfasta]
    if args.remove_internal_priming and isannot:
        count_cmd += ['--permissive_last_exons']
    return tuple(count_cmd)


def transcriptomealignandcount(args, inputreads, alignfasta, refbed, outputname, mapfile, isannot):
    # minimap (results are piped into count_sam_transcripts.py)
    # '--split-prefix', 'minimap2transcriptomeindex', doesn't work with MD tag
    if type(inputreads) == str:
        inputreads = [inputreads]
    mm2_cmd = tuple(
        ['minimap2', '-a', '-t', str(args.threads), '-N', '4', '--MD'] + args.mm2_args + [alignfasta] + inputreads)
    # FIXME add in step to filter out chimeric reads here
    # FIXME really need to go in and check on how count_sam_transcripts is working
    count_cmd = getcountsamcommand(args, refbed, outputname, mapfile, isannot)
    pipettor.run([mm2_cmd, count_cmd])

# START METHODS TO EDIT - HARRISON

def getbestends(currgroup, end_window):
    bestends = []
    if len(currgroup) > int(end_window):
        allstarts = Counter([x[0] for x in currgroup])
        allends = Counter([x[1] for x in currgroup])
        for start1, end1, strand1, name1 in currgroup:
            weightedscore = allstarts[start1] + allends[end1]
            bestends.append((weightedscore, start1, end1, strand1, name1))
    else:
        for start1, end1, strand1, name1 in currgroup:
            score, weightedscore = 0, 0
            for start2, end2, strand2, name2 in currgroup:
                if abs(start1 - start2) <= end_window and abs(end1 - end2) <= end_window:
                    score += 2
                    weightedscore += ((end_window - abs(start1 - start2)) / end_window) + \
                                     ((end_window - abs(end1 - end2)) / end_window)
            bestends.append((weightedscore, start1, end1, strand1, name1))
    bestends.sort(reverse=True)
    # DO I WANT TO ADD CORRECTION TO NEARBY ANNOTATED TSS/TTS????
    return bestends[0]


def combinefinalends(currgroup):
    if len(currgroup) == 1:
        return currgroup[0]
    else:
        currgroup.sort(key=lambda x: x[2])  # sort by marker
        allreads = [y for x in currgroup for y in x[-1]]
        if currgroup[0] != 'a':  # if no annotated iso, sort further
            currgroup.sort(key=lambda x: len(x[-1]), reverse=True)
        bestiso = currgroup[0]
        bestiso[-1] = allreads
        return bestiso


def groupreadsbyends(readinfos, sortindex, end_window):
    sortedends = sorted(readinfos, key=lambda x:x[sortindex])
    newgroups, group = [], []
    lastedge = 0
    for isoinfo in sortedends:
        edge = isoinfo[sortindex]
        if edge-lastedge <= end_window:
            group.append(isoinfo)
        else:
            if len(group) > 0: newgroups.append(group)
            group = [isoinfo]
        lastedge = edge
    if len(group) > 0:
        newgroups.append(group)
    return newgroups

# MAIN METHOD - CALLS OTHERS IN GROUP
# readends is a list containing elements with: (read.start, read.end, read.strand, read.name)
# If the reads are spliced, the group will contain only the info for reads with a shared splice junction
# if the reads are unspliced, the group will contain info for all unspliced reads in a given chromosome/region, depending on how you're parallelizing
# The output is a list containing elements with: [weighted.end.score, group.start, group.end, group.strand, representative.read.name, [list of all read names in group]]
# weighted.end.score (represents how many reads have ends similar to this exact position)
# You can rewrite this to redo the grouping method, but please maintain the inputs and outputs.
def collapseendgroups(end_window, readends, dogetbestends=True):
    startgroups = groupreadsbyends(readends, 0, end_window)
    allendgroups, isoendgroups = [], []
    for startgroup in startgroups:
        allendgroups.extend(groupreadsbyends(startgroup, 1, end_window))
    for endgroup in allendgroups:
        if dogetbestends:
            isoendgroups.append(list(getbestends(endgroup, end_window)) + [[x[3] for x in endgroup]])
        else:
            isoendgroups.append(combinefinalends(endgroup))
    return isoendgroups

# END METHODS TO EDIT





def identify_good_match_to_annot(args, tempprefix, thischrom, annottranscripttoexons, alltranscripts, genome):
    goodaligntoannot, firstpasssingleexons, supannottranscripttojuncs = [], set(), {}
    if not args.noaligntoannot and len(alltranscripts) > 0:
        transcripttostrand = {}
        with open(tempprefix + '.annotated_transcripts.bed', 'w') as annotbed, open(
                tempprefix + '.annotated_transcripts.fa', 'w') as annotfa:
            for transcript, gene, strand in alltranscripts:
                transcripttostrand[(transcript, gene)] = strand
                exons = annottranscripttoexons[(transcript, gene)]
                start, end = exons[0][0], exons[-1][1]
                estarts = [x[0] - start for x in exons]
                esizes = [x[1] - x[0] for x in exons]
                bedline = [thischrom, start, end, transcript + '_' + gene, '.', strand, start, end, '0', len(exons),
                           ','.join([str(x) for x in esizes]), ','.join([str(x) for x in estarts])]
                if strand == '-':
                    exons = exons[::-1]
                exonseq = []
                for i in range(len(exons)):
                    thisexonseq = genome.fetch(thischrom, exons[i][0], exons[i][1])
                    if strand == '-':
                        thisexonseq = revcomp(thisexonseq)
                    exonseq.append(thisexonseq)
                annotbed.write('\t'.join([str(x) for x in bedline]) + '\n')
                annotfa.write('>' + transcript + '_' + gene + '\n')
                annotfa.write(''.join(exonseq) + '\n')
        transcriptomealignandcount(args, tempprefix + '.reads.fasta',
                                   tempprefix + '.annotated_transcripts.fa',
                                   tempprefix + '.annotated_transcripts.bed',
                                   tempprefix + '.matchannot.counts.tsv',
                                   tempprefix + '.matchannot.read.map.txt', True)
        with open(tempprefix + '.matchannot.bed', 'w') as annotbed:
            for line in open(tempprefix + '.matchannot.read.map.txt'):
                striso, reads = line.rstrip().split('\t', 1)
                reads = reads.split(',')
                goodaligntoannot.extend(reads)
                if len(reads) >= args.support:
                    transcript = '_'.join(striso.split('_')[:-1])
                    gene = striso.split('_')[-1]
                    if (transcript, gene) in annottranscripttoexons:
                        exons = annottranscripttoexons[(transcript, gene)]
                        start, end = exons[0][0], exons[-1][1]
                        estarts = [x[0] - start for x in exons]
                        esizes = [x[1] - x[0] for x in exons]
                        strand = transcripttostrand[(transcript, gene)]
                        bedline = [thischrom, start, end, transcript + '_' + gene, len(reads), strand, start, end, '0',
                                   len(exons), ','.join([str(x) for x in esizes]), ','.join([str(x) for x in estarts])]
                        annotbed.write('\t'.join([str(x) for x in bedline]) + '\n')
                        firstpasssingleexons.update(set(exons))
                        annotjuncs = tuple([(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)])
                        supannottranscripttojuncs[(transcript, gene)] = (len(reads), annotjuncs)
    else:
        with open(tempprefix + '.matchannot.counts.tsv', 'w') as countsout, open(
                tempprefix + '.matchannot.read.map.txt', 'w') as mapout, open(tempprefix + '.matchannot.bed',
                                                                              'w') as annotbed:
            pass
    goodaligntoannot = set(goodaligntoannot)
    return goodaligntoannot, firstpasssingleexons, supannottranscripttojuncs


def filtercorrectgroupreads(args, tempprefix, rchrom, rstart, rend, samfile, goodaligntoannot, intervalTree,
                            junctionBoundaryDict):
    sjtoends = {}
    shortchromfasta = open(tempprefix + 'reads.notannotmatch.fasta', 'w')
    c = 0
    for read in samfile.fetch(rchrom, int(rstart), int(rend)):
        if not read.is_secondary and (not read.is_supplementary or args.keep_sup):
            if read.query_name not in goodaligntoannot:
                c += 1
                shortchromfasta.write('>' + read.query_name + '\n')
                shortchromfasta.write(read.get_forward_sequence() + '\n')
                if read.mapping_quality >= args.quality:
                    junc_strand = inferMM2JuncStrand(read)
                    if junc_strand == 'ambig':
                       junc_strand = "-" if read.is_reverse else "+"
                    bedread = BedRead()
                    bedread.generate_from_cigar(read.reference_start, read.is_reverse, read.cigartuples,
                                                read.query_name,
                                                read.reference_name, read.mapping_quality, junc_strand)
                    correctedread = correctsingleread(bedread, intervalTree, junctionBoundaryDict)
                    if correctedread:
                        junckey = tuple(sorted(correctedread.juncs))
                        if junckey not in sjtoends:
                            sjtoends[junckey] = []
                        sjtoends[junckey].append((correctedread.start, correctedread.end, correctedread.strand, correctedread.name))
    shortchromfasta.close()
    return sjtoends


def groupfirstpasssingleexon(goodendswithsupreads):
    lastend, furthergroups, thisgroup = 0, [], []
    goodendswithsupreads.sort(key=lambda x: [x[1], x[2]])
    for weightedscore, start, end, strand, name, endsupport in goodendswithsupreads:
        if (start >= lastend or (lastend - start) / (end - start) > 0.5) and len(thisgroup) > 0:
            furthergroups.append(thisgroup)
            thisgroup = []
        thisgroup.append([weightedscore, start, end, strand, name, endsupport])
        if end > lastend:
            lastend = end
    if len(thisgroup) > 0:
        furthergroups.append(thisgroup)
    logging.info(f'single exon groups: {len(furthergroups)}')
    return furthergroups


def filterendsbyredundantandsupport(args, goodendswithsupreads):
    bestends = []
    for i in range(len(goodendswithsupreads)):
        goodendswithsupreads[i][-1] = len(goodendswithsupreads[i][-1])  # last val is read support
    goodendswithsupreads.sort(key=lambda x: [x[0], x[2] - x[1]],
                              reverse=True)  # first by weighted score, then by length
    juncsupport = sum([x[-1] for x in goodendswithsupreads])
    if juncsupport >= args.support:
        if args.no_redundant == 'none':
            if goodendswithsupreads[0][-1] < args.support:
                goodendswithsupreads = [goodendswithsupreads[0]]
                goodendswithsupreads[0][-1] = juncsupport
            else:
                goodendswithsupreads = [x for x in goodendswithsupreads if x[-1] >= args.support]
                goodendswithsupreads = goodendswithsupreads[:args.max_ends]  # select only top most supported ends
            for theseends in goodendswithsupreads:
                bestends.append(theseends)
        else:
            # best_only uses the default sorting
            if args.no_redundant == 'longest':
                goodendswithsupreads.sort(reverse=True, key=lambda x: x[2] - x[1])
            thisbest = goodendswithsupreads[0]
            thisbest[-1] = juncsupport  # all reads for junc are counted as support
            bestends.append(thisbest)
    return bestends


def processjuncstofirstpassisos(args, tempprefix, thischrom, sjtoends, firstpasssingleexons):
    firstpassunfiltered, firstpassjunctoname = {}, {}
    with open(tempprefix + '.firstpass.unfiltered.bed', 'w') as isoout, open(tempprefix + '.firstpass.reallyunfiltered.bed', 'w') as isoout2:
        for juncs in sjtoends:
            goodendswithsupreads = collapseendgroups(args.end_window, sjtoends[juncs])
            for bestscore, beststart, bestend, beststrand, bestname, thesereads in goodendswithsupreads:
                thisscore = len(thesereads)
                thisiso = BedRead()
                thisiso.generatefromvals(thischrom, beststart, bestend, bestname, thisscore, beststrand, juncs)
                isoout2.write('\t'.join(thisiso.getbedline()) + '\n')
            if juncs == ():
                bestends = [x[:-1] + [len(x[-1])] for x in goodendswithsupreads if len(x[-1]) >= args.support]
            else: bestends = filterendsbyredundantandsupport(args, goodendswithsupreads)
            for bestscore, beststart, bestend, beststrand, bestname, thisscore in bestends:
                thisiso = BedRead()
                thisiso.generatefromvals(thischrom, beststart, bestend, bestname, thisscore, beststrand, juncs)
                firstpassunfiltered[bestname] = thisiso
                isoout.write('\t'.join(thisiso.getbedline()) + '\n')
                if juncs == ():
                    firstpasssingleexons.add((thisiso.exons[0][0], thisiso.exons[0][1], thisiso.name))
                else:
                    for j in juncs:
                        if j not in firstpassjunctoname:
                            firstpassjunctoname[j] = set()
                        firstpassjunctoname[j].add(bestname)
                    firstpasssingleexons.update(thisiso.exons)

    firstpasssingleexons = sorted(list(firstpasssingleexons))
    return firstpassunfiltered, firstpassjunctoname, firstpasssingleexons


def filtersplicediso(args, thisiso, firstpassjunctoname, firstpassunfiltered, junctogene, supannottranscripttojuncs,
                     annottranscripttoexons):
    isoswithsimilarjuncs = set()
    for j in thisiso.juncs:
        isoswithsimilarjuncs.update(firstpassjunctoname[j])
        if j in junctogene:
            isoswithsimilarjuncs.update(junctogene[j])  # annot isos
    issubset = [0, 0]  # first exon is a subset, last exon is a subset
    firstexon, lastexon = thisiso.exons[0], thisiso.exons[-1]
    superset_support = []
    for otherisoname in isoswithsimilarjuncs:
        if otherisoname != thisiso.name:
            if type(otherisoname) == tuple:
                if otherisoname in supannottranscripttojuncs:
                    otherisoscore, otherisojuncs = supannottranscripttojuncs[otherisoname]
                    otherisoexons = annottranscripttoexons[otherisoname]
                else:
                    continue
            else:
                otheriso = firstpassunfiltered[otherisoname]
                otherisoscore, otherisojuncs, otherisoexons = otheriso.score, otheriso.juncs, otheriso.exons
            if len(thisiso.juncs) < len(otherisojuncs):
                thisisostrjuncs, otherisostrjuncs = str(thisiso.juncs)[1:-1].rstrip(','), str(otherisojuncs)[1:-1]
                if thisisostrjuncs in otherisostrjuncs:  # if junctions are subset
                    # check whether first + last exon overlap
                    for i in range(len(otherisoexons)):

                        otherexon = otherisoexons[i]
                        if i == 0 or i == len(otherisoexons) - 1:  # is first or last exon of other transcript
                            if firstexon[1] == otherexon[1] or lastexon[0] == otherexon[0]:
                                if firstexon[1] == otherexon[1]:
                                    issubset[0] = 1
                                elif lastexon[0] == otherexon[0]:
                                    issubset[1] = 1
                                superset_support.append(otherisoscore)
                        else:
                            if firstexon[1] == otherexon[1] and firstexon[0] >= otherexon[0] - 10:
                                issubset[0] = 1
                                superset_support.append(otherisoscore)
                            if lastexon[0] == otherexon[0] and lastexon[1] <= otherexon[1] + 10:
                                issubset[1] = 1
                                superset_support.append(otherisoscore)
    if sum(issubset) < 2:  # both first and last exon have to overlap
        return True
    elif args.filter != 'nosubset':
        if thisiso.score >= args.support and thisiso.score > max(superset_support) * 1.2:
            return True
    return False


def filtersingleexoniso(args, groupediso, currgroup, firstpassunfiltered):
    thisiso = firstpassunfiltered[groupediso[2]]
    exprcomp = []
    iscontained = False
    for compiso in currgroup:
        if compiso != groupediso:
            if compiso[0] - 10 <= groupediso[0] and groupediso[1] <= compiso[1] + 10:
                if len(compiso) == 2 or args.filter == 'nosubset':  # is exon from spliced transcript
                    iscontained = True
                    break  # filter out
                else:  # is other single exon - check relative expression
                    otherscore = firstpassunfiltered[compiso[2]].score
                    thisscore = thisiso.score
                    if thisscore >= args.support and otherscore * 1.2 < thisscore:
                        exprcomp.append(True)
                    else:
                        exprcomp.append(False)
    if not iscontained and all(exprcomp):
        return True
    else:
        return False


def filtersingleexongroup(args, currgroup, firstpassunfiltered, firstpass):
    for groupediso in currgroup:
        if len(groupediso) == 3:  # is single exon with name
            if filtersingleexoniso(args, groupediso, currgroup, firstpassunfiltered):
                firstpass[groupediso[2]] = firstpassunfiltered[groupediso[2]]
    return firstpass


def filterallsingleexon(args, firstpasssingleexons, firstpassunfiltered, firstpass):
    lastend = 0
    currgroup = []
    for isoinfo in firstpasssingleexons:
        start, end = isoinfo[0], isoinfo[1]
        if start < lastend:
            currgroup.append(isoinfo)
        else:
            if len(currgroup) > 0:
                firstpass = filtersingleexongroup(args, currgroup, firstpassunfiltered, firstpass)
            currgroup = [isoinfo]
        if end > lastend:
            lastend = end
    if len(currgroup) > 0:
        firstpass = filtersingleexongroup(args, currgroup, firstpassunfiltered, firstpass)

    return firstpass


def filterfirstpassisos(args, firstpassunfiltered, firstpassjunctoname, firstpasssingleexons, junctogene,
                        supannottranscripttojuncs, annottranscripttoexons):
    if args.filter == 'ginormous':
        firstpass = firstpassunfiltered
    else:
        firstpass = {}
        for isoname in firstpassunfiltered:
            thisiso = firstpassunfiltered[isoname]
            if thisiso.juncs != ():
                if args.filter == 'comprehensive':
                    firstpass[isoname] = thisiso
                else:
                    passesfiltering = filtersplicediso(args, thisiso, firstpassjunctoname, firstpassunfiltered,
                                                       junctogene, supannottranscripttojuncs, annottranscripttoexons)
                    if passesfiltering:
                        firstpass[isoname] = thisiso
        # HANDLE SINGLE EXONS SEPARATELY - group first - one traversal of list
        firstpass = filterallsingleexon(args, firstpasssingleexons, firstpassunfiltered, firstpass)

    return firstpass


def combinetempfilesbysuffix(args, tempprefixes, suffixes):
    for filesuffix in suffixes:
        with open(args.output + filesuffix, 'wb') as allannotcounts:
            for tempprefix in tempprefixes:
                with open(tempprefix + filesuffix, 'rb') as fd:
                    shutil.copyfileobj(fd, allannotcounts, 1024 * 1024 * 10)


def getsplicedisogenehits(thisiso, junctogene, genetoannotjuncs):
    gene_hits = {}
    if thisiso.juncs != ():
        for j in thisiso.juncs:
            if j in junctogene:
                for transcript, gene in junctogene[j]:
                    if gene not in gene_hits:
                        gene_hits[gene] = [0, -1 * len(genetoannotjuncs[gene])]
                    gene_hits[gene][0] += 1
    return gene_hits


def getsingleexongenehits(thisiso, allannotse):
    gene_hits = {}
    thisexon = thisiso.exons[0]
    index = bin_search(thisexon, allannotse)
    for annotexoninfo in allannotse[index - 2:index + 2]:  # start, end, strand, gene
        overlap = min(thisexon[1], annotexoninfo[1]) - max(thisexon[0], annotexoninfo[0])
        if overlap > 0:
            # base coverage of long-read isoform by the annotated isoform
            thisisofrac = float(overlap) / (thisexon[1] - thisexon[0])
            # base coverage of the annotated isoform by the long-read isoform
            annotfrac = float(overlap) / (annotexoninfo[1] - annotexoninfo[0])
            if thisisofrac > 0.5 and annotfrac > 0.8:
                if annotexoninfo[3] not in gene_hits or thisisofrac > gene_hits[annotexoninfo[3]][0]:
                    gene_hits[annotexoninfo[3]] = [thisisofrac, annotfrac]
    return gene_hits


def getgenenamesandwritefirstpass(tempprefix, thischrom, firstpass, juncstotranscript, junctogene, allannotse,
                                  genetoannotjuncs, genetostrand, genome):
    with open(tempprefix + '.firstpass.bed', 'w') as isoout, open(tempprefix + '.firstpass.fa', 'w') as seqout:
        # THIS IS WHERE WE CAN GET GENES AND ADJUST NAMES
        annotnametousedcounts = {}
        for isoname in firstpass:
            thisiso = firstpass[isoname]
            # Adjust name based on annotation
            thistranscript, thisgene = thisiso.name, None#thischrom.replace('_', '-') + ':' + str(round(thisiso.start, -3))
            if thisiso.juncs != () and thisiso.juncs in juncstotranscript:
                thistranscript, thisgene = juncstotranscript[thisiso.juncs]
                if thistranscript in annotnametousedcounts:
                    annotnametousedcounts[thistranscript] += 1
                    thistranscript = thistranscript + '-endvar' + str(annotnametousedcounts[thistranscript])
                else:
                    annotnametousedcounts[thistranscript] = 1
            else:
                if thisiso.juncs != ():
                    gene_hits = getsplicedisogenehits(thisiso, junctogene, genetoannotjuncs)
                else:
                    gene_hits = getsingleexongenehits(thisiso, allannotse)  # single exon

                if gene_hits:
                    sortedgenes = sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)
                    thisgene = sortedgenes[0][0]
            if thisgene:
                thisiso.strand = genetostrand[thisgene]
            else:
                thisgene = thischrom.replace('_', '-') + ':' + str(round(thisiso.start, -3))
            thisiso.name = thistranscript + '_' + thisgene
            isoout.write('\t'.join(thisiso.getbedline()) + '\n')
            seqout.write('>' + thisiso.name + '\n')
            seqout.write(thisiso.getsequence(genome) + '\n')


def getisogenefromname(isogene):
    iso = '_'.join(isogene.split('_')[:-1])
    gene = isogene.split('_')[-1]
    return iso, gene


def processdetectedisos(args, mapfile, bedfile, marker, genetojuncstoends):
    ogisotoreads = {}
    for line in open(mapfile):
        isogene, reads = line.rstrip().split('\t', 1)
        reads = reads.split(',')
        support = len(reads)
        if support >= args.support:
            iso, gene = getisogenefromname(isogene)
            isoid = (marker, iso)
            ogisotoreads[isoid] = reads

    for line in open(bedfile):
        line = line.rstrip().split('\t')
        chrom, start, end, name, score, strand = line[:6]
        iso, gene = getisogenefromname(name)
        isoid = (marker, iso)
        if isoid in ogisotoreads:
            start, end = int(start), int(end)
            estarts = [int(x) for x in line[11].rstrip(',').split(',')]
            esizes = [int(x) for x in line[10].rstrip(',').split(',')]
            juncs = []
            for i in range(len(estarts) - 1):
                juncs.append((start + estarts[i] + esizes[i], start + estarts[i + 1]))
            junckey = (chrom, strand, tuple(juncs))
            if gene not in genetojuncstoends:
                genetojuncstoends[gene] = {}
            if junckey not in genetojuncstoends[gene]:
                genetojuncstoends[gene][junckey] = []
            genetojuncstoends[gene][junckey].append([start, end, isoid, ogisotoreads[isoid]])

    return genetojuncstoends


def revcomp(seq):
    compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 
                'R':'Y', 'Y':'R','K':'M','M':'K','S':'S','W':'W', 'B':'V','V':'B','D':'H','H':'D'}
    seq = seq.upper()
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq[::-1])


def getbedgtfoutfrominfo(endinfo, chrom, strand, juncs, gene, genome):
    start, end, isoid, readnames = endinfo
    score = min(len(readnames), 1000)
    marker, iso = isoid
    estarts, esizes = getexonsfromjuncs(juncs, start, end)
    bedline = [chrom, start, end, iso + '_' + gene, score, strand, start, end,
               getrgb(iso, strand, juncs), len(estarts), ','.join([str(x) for x in esizes]),
               ','.join([str(x) for x in estarts])]
    gtflines = []
    gtflines.append([chrom, 'FLAIR', 'transcript', start + 1, end, score, strand, '.',
                     'gene_id "' + gene + '"; transcript_id "' + iso + '";'])
    exons = [(start + estarts[i], start + estarts[i] + esizes[i]) for i in range(len(estarts))]
    if strand == '-':
        exons = exons[::-1]
    exonseq = []
    for i in range(len(exons)):
        gtflines.append([chrom, 'FLAIR', 'exon', exons[i][0] + 1, exons[i][1], score, strand, '.',
                         'gene_id "' + gene + '"; transcript_id "' + iso + '"; exon_number ' + str(i + 1)])
        thisexonseq = genome.fetch(chrom, exons[i][0], exons[i][1])
        if strand == '-':
            thisexonseq = revcomp(thisexonseq)
        exonseq.append(thisexonseq)
    return '\t'.join([str(x) for x in bedline]) + '\n', gtflines, ''.join(exonseq)


def combineannotnovelwriteout(args, genetojuncstoends, genome):
    with open(args.output + '.isoforms.bed', 'w') as isoout, open(args.output + '.isoform.read.map.txt', 'w') as mapout, open(
            args.output + '.isoforms.gtf', 'w') as gtfout, open(args.output + '.isoforms.fa', 'w') as seqout, open(args.output + '.isoform.counts.txt', 'w') as countsout:
        for gene in genetojuncstoends:
            gtflines, tstarts, tends = [], [], []
            for chrom, strand, juncs in genetojuncstoends[gene]:
                endslist = genetojuncstoends[gene][(chrom, strand, juncs)]
                endslist = collapseendgroups(args.end_window, endslist, False)
                # FIXME could try accounting for all reads assigned to isoforms - assign them to closest ends
                # not sure how much of an issue this is
                if juncs != ():
                    if args.no_redundant == 'best_only':
                        endslist.sort(key=lambda x: [len(x[-1]), x[1] - x[0]], reverse=True)
                        endslist = [endslist[0]]
                    elif args.no_redundant == 'longest':
                        endslist.sort(key=lambda x: [x[1] - x[0]], reverse=True)
                        endslist = [endslist[0]]
                    else:
                        endslist.sort(key=lambda x: [len(x[-1]), x[1] - x[0]], reverse=True)
                        endslist = endslist[:args.max_ends]
                nametousedcounts = {}
                for isoinfo in endslist:
                    marker, iso = isoinfo[2]
                    iso = iso.split('-endvar')[0]
                    if iso in nametousedcounts:
                        nametousedcounts[iso] += 1
                        iso = iso + '-endvar' + str(nametousedcounts[iso])
                    else:
                        nametousedcounts[iso] = 1
                    isoinfo[2] = (marker, iso)
                    bedline, gtffortranscript, tseq = getbedgtfoutfrominfo(isoinfo, chrom, strand, juncs, gene, genome)
                    isoout.write(bedline)
                    tstarts.append(isoinfo[0])
                    tends.append(isoinfo[1])
                    gtflines.extend(gtffortranscript)
                    mapout.write(iso + '_' + gene + '\t' + ','.join(isoinfo[3]) + '\n')
                    countsout.write(iso + '_' + gene + '\t' + str(len(isoinfo[3])) + '\n')
                    seqout.write('>' + iso + '_' + gene + '\n')
                    seqout.write(tseq + '\n')
            gtflines.insert(0, [chrom, 'FLAIR', 'gene', min(tstarts) + 1, max(tends), '.', gtflines[0][6], '.',
                                'gene_id "' + gene + '";'])
            for g in gtflines:
                gtfout.write('\t'.join([str(x) for x in g]) + '\n')

def decide_parallel_mode(parallel_option, genomealignedbam):
    ###parallel_option format already validated in option validation method
    if parallel_option in {'bychrom', 'byregion'}: return parallel_option
    else:
        if parallel_option[-2:] == 'GB': threshold = float(parallel_option[5:-2])
        else: threshold = float(parallel_option[5:])
        filesizeGB = os.path.getsize(genomealignedbam) / 1e+9
        if filesizeGB > threshold: return 'byregion'
        else: return 'bychrom'

def runcollapsebychrom(listofargs):
    args, tempprefix, splicesiteannot_chrom, juncstotranscript, junctogene, allannotse, genetoannotjuncs, genetostrand, annottranscripttoexons, allannottranscripts = listofargs
    # first extract reads for chrom as fasta
    tempsplit = tempprefix.split('/')[-1].split('-')
    rchrom, rstart, rend = '-'.join(tempsplit[:-2]), tempsplit[-2], tempsplit[-1]
    pipettor.run([('samtools', 'view', '-h', args.genomealignedbam, rchrom + ':' + rstart + '-' + rend),
                  ('samtools', 'fasta', '-')],
                 stdout=open(tempprefix + '.reads.fasta', 'w'))
    # then align reads to transcriptome and run count_sam_transcripts
    genome = pysam.FastaFile(args.genome)
    goodaligntoannot, firstpasssingleexons, supannottranscripttojuncs = identify_good_match_to_annot(args, tempprefix,
                                                                                                 rchrom,
                                                                                                 annottranscripttoexons,
                                                                                                 allannottranscripts,
                                                                                                 genome)
    # load splice junctions for chrom
    intervalTree, junctionBoundaryDict = buildIntervalTree(splicesiteannot_chrom, args.ss_window, rchrom, False)

    samfile = pysam.AlignmentFile(args.genomealignedbam, 'rb')
    sjtoends = filtercorrectgroupreads(args, tempprefix, rchrom, rstart, rend, samfile, goodaligntoannot, intervalTree,
                                       junctionBoundaryDict)
    samfile.close()

    firstpassunfiltered, firstpassjunctoname, firstpasssingleexons = processjuncstofirstpassisos(args, tempprefix,
                                                                                                 rchrom, sjtoends,
                                                                                                 firstpasssingleexons)

    firstpass = filterfirstpassisos(args, firstpassunfiltered, firstpassjunctoname, firstpasssingleexons,
                                    junctogene, supannottranscripttojuncs, annottranscripttoexons)
    temptoremove = [tempprefix + '.reads.fasta', tempprefix + 'reads.notannotmatch.fasta']
    if not args.noaligntoannot and len(allannottranscripts) > 0:
        temptoremove.extend([tempprefix + '.annotated_transcripts.bed', tempprefix + '.annotated_transcripts.fa'])
    if len(firstpass.keys()) > 0:
        getgenenamesandwritefirstpass(tempprefix, rchrom, firstpass, juncstotranscript, junctogene,
                                      allannotse, genetoannotjuncs, genetostrand, genome)
        transcriptomealignandcount(args, tempprefix + 'reads.notannotmatch.fasta',
                                   tempprefix + '.firstpass.fa',
                                   tempprefix + '.firstpass.bed',
                                   tempprefix + '.novelisos.counts.tsv',
                                   tempprefix + '.novelisos.read.map.txt', False)
        temptoremove.extend([tempprefix + '.firstpass.fa'])
    else:
        with open(tempprefix + '.firstpass.fa', 'w') as faout, open(tempprefix + '.firstpass.bed', 'w') as bedout, \
                open(tempprefix + '.novelisos.counts.tsv', 'w') as countsout, open(
            tempprefix + '.novelisos.read.map.txt', 'w') as mapout:
            pass
    genome.close()
    if not args.keep_intermediate:
        for f in temptoremove:
            os.remove(f)


def collapsefrombam():
    args = get_args()
    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)
    allchrom = genome.references
    logging.info('making temp dir')
    tempDir = makecorrecttempdir()
    logging.info('Getting regions')
    allregions = []
    if decide_parallel_mode(args.parallelmode, args.genomealignedbam) == 'bychrom':
        for chrom in allchrom:
            chromsize = genome.get_reference_length(chrom)
            allregions.append((chrom, 0, chromsize))
    else:
        pipettor.run(['flair_partition', '--min_partition_items', '1000', '--threads', str(args.threads), '--bam=' + args.genomealignedbam,
                      tempDir + 'regions.bed'])
        for line in open(tempDir + 'regions.bed'):
            line = line.rstrip().split('\t')
            chrom, start, end = line[0], int(line[1]), int(line[2])
            allregions.append((chrom, start, end))
    logging.info(f'Number of regions {len(allregions)}')
    logging.info('Generating splice site database')
    knownchromosomes, annotationFiles = generateKnownSSDatabase(args, tempDir)

    regionstoannotdata = {}
    if args.gtf:
        logging.info('Extracting annotation from GTF')
        regionstoannotdata = getannotinfo(args.gtf, allregions)
    logging.info('splitting by chunk')
    chunkcmds = []
    tempprefixes = []
    for rchrom, rstart, rend in allregions:
        if rchrom in knownchromosomes:
            juncstotranscript, junctogene, allannotse, genetoannotjuncs, genetostrand, annottranscripttoexons, allannottranscripts = {}, {}, [], {}, {}, {}, []
            if args.gtf:
                juncstotranscript, junctogene, allannotse, genetoannotjuncs, genetostrand, annottranscripttoexons, allannottranscripts = \
                regionstoannotdata[(rchrom, rstart, rend)].returndata()

            splicesiteannot_chrom = annotationFiles[rchrom]
            tempprefix = tempDir + '-'.join([rchrom, str(rstart), str(rend)])
            chunkcmds.append([args, tempprefix, splicesiteannot_chrom, juncstotranscript,
                              junctogene, allannotse, genetoannotjuncs, genetostrand, annottranscripttoexons,
                              allannottranscripts])
            tempprefixes.append(tempprefix)
    mp.set_start_method('fork')
    p = mp.Pool(args.threads)
    childErrs = set()
    c = 1
    for i in p.imap(runcollapsebychrom, chunkcmds):
        logging.info(f'\rdone running chunk {c} of {len(chunkcmds)}')
        childErrs.add(i)
        c += 1
    p.close()
    p.join()
    if len(childErrs) > 1:
        raise ValueError(childErrs)

    if not args.noaligntoannot:
        combinetempfilesbysuffix(args, tempprefixes,
                                 ['.matchannot.counts.tsv', '.matchannot.read.map.txt', '.matchannot.bed'])
    combinetempfilesbysuffix(args, tempprefixes,
                             ['.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed', '.firstpass.bed', '.novelisos.counts.tsv',
                              '.novelisos.read.map.txt'])

    if not args.keep_intermediate:
        shutil.rmtree(tempDir)

    if not args.noaligntoannot:
        genetojuncstoends = processdetectedisos(args, args.output + '.matchannot.read.map.txt',
                                            args.output + '.matchannot.bed', 'a', {})
    else: genetojuncstoends = {}
    genetojuncstoends = processdetectedisos(args, args.output + '.novelisos.read.map.txt',
                                            args.output + '.firstpass.bed', 'n',
                                            genetojuncstoends)

    combineannotnovelwriteout(args, genetojuncstoends, genome)
    if not args.keep_intermediate:
        files_to_remove = ['.firstpass.reallyunfiltered.bed',
                           '.firstpass.unfiltered.bed',
                           '.firstpass.bed',
                           '.novelisos.counts.tsv',
                           '.novelisos.read.map.txt']
        if not args.noaligntoannot:
            files_to_remove += ['.matchannot.bed',
                           '.matchannot.counts.tsv',
                           '.matchannot.read.map.txt']
        for f in files_to_remove:
            os.remove(args.output + f)

    if args.predictCDS:
        prodcmd = ('predictProductivity',
                   '-i', args.output+'.isoforms.bed',
                   '-o', args.output + '.isoforms.CDS',
                   '--gtf', args.gtf,
                   '--genome_fasta', args.genome,
                   '--longestORF')
        pipettor.run([prodcmd])
        # os.rename(args.output + '.isoforms.CDS.bed', args.output + '.isoforms.bed')
        os.remove(args.output + '.isoforms.CDS.info.tsv')
    genome.close()


if __name__ == "__main__":
    collapsefrombam()
