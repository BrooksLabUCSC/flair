#! /usr/bin/env python3

import sys, argparse, os, pipettor, glob, uuid, shutil
import pysam
from flair_align import inferMM2JuncStrand, intronChainToestarts
from ssUtils import addOtherJuncs, gtfToSSBed
from ssPrep import buildIntervalTree, ssCorrect, juncsToBed12
from multiprocessing import Pool
import time


# export PATH="/private/groups/brookslab/cafelton/git-flair/flair/bin:/private/groups/brookslab/cafelton/git-flair/flair/src/flair:$PATH"

def getargs():
    parser = argparse.ArgumentParser(description='flair-collapse parse options',
                                     usage='''python3 flair_straightfrombam.py -g genome.fa -b reads.genomealigned.bam [options]''')
    parser.add_argument('-b', '--genomealignedbam',
                        help='Optional: sorted and indexed bam file (or files) aligned to the genome. Only use this if you have already aligned reads to the genome for some other purpose')
    parser.add_argument('-g', '--genome', type=str,
                        help='FastA of reference genome, can be minimap2 indexed')
    parser.add_argument('-o', '--output', default='flair.collapse',
                        help='output file name base for FLAIR isoforms (default: flair.collapse)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (4)')
    parser.add_argument('-f', '--gtf', default='',
                        help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')

    parser.add_argument('-j', '--shortread', type=str, default='',
                        help='bed format splice junctions from short-read sequencing')
    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (15)')

    parser.add_argument('-s', '--support', type=float, default=3.0,
                        help='''minimum number of supporting reads for an isoform;
                if s < 1, it will be treated as a percentage of expression of the gene (3)''')
    parser.add_argument('--stringent', default=False, action='store_true',
                        help='''specify if all supporting reads need to be full-length
                (80%% coverage and spanning 25 bp of the first and last exons)''')
    parser.add_argument('--check_splice', default=False, action='store_true',
                        help='''enforce coverage of 4 out of 6 bp around each splice site and no
                insertions greater than 3 bp at the splice site''')
    parser.add_argument('--quality', type=int, default=0,
                        help='minimum MAPQ of read assignment to an isoform (1)')

    parser.add_argument('--noaligntoannot', default=False, action='store_true',
                        help='''related to old annotation_reliant, now specify if you don't want an initial alignment to the annotated sequences and only want transcript detection from the genomic alignment. Will be slightly faster but less accurate if the annotation is good''')
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

    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        sys.exit(1)

    args, unknown = parser.parse_known_args()
    if unknown:
        sys.stderr.write('unrecognized arguments: {}\n'.format(' '.join(unknown)))

    args = checkfilepaths(args)
    # args.quality = '0' if args.trust_ends else args.quality
    # if args.mm2_args:
    #     args.mm2_args = args.mm2_args.split(',')
    return args


def checkfilepaths(args):
    if not args.genomealignedbam:
        sys.stderr.write(f'Please include the --genomealignedbam option\n')
        sys.exit(1)
    if not args.genome:
        sys.stderr.write(f'Please include the --genome option\n')
        sys.exit(1)
    if not os.path.exists(args.genomealignedbam):
        sys.stderr.write(f'Aligned reads file path does not exist: {args.genomealignedbam}\n')
        sys.exit(1)
    if not os.path.exists(args.genome):
        sys.stderr.write('Genome file path does not exist: {}\n'.format(args.genome))
        sys.exit(1)
    return args


def makecorrecttempdir():
    tempDirName = str(uuid.uuid4())
    try:
        current_directory = os.getcwd()
        tempDir = os.path.join(current_directory, tempDirName)
        os.mkdir(tempDir)
    except OSError:
        print("Creation of the directory %s failed" % tempDirName, file=sys.stderr)
        sys.exit(1)
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

    sys.stderr = open(os.path.join(tempDir, 'splicelog.txt'), 'w')
    if args.gtf:
        juncs, chromosomes, knownSS = gtfToSSBed(args.gtf, knownSS, False, False, False)

    # Do the same for the other juncs file.
    if args.shortread:
        juncs, chromosomes, addFlag = addOtherJuncs(juncs, args.shortread, chromosomes, args.genome,
                                                    False, knownSS, False, False)
        if addFlag == False:
            sys.stderr.write('\nERROR Added no extra junctions from {}\n\n'.format(args.shortread))
            sys.exit(1)
    knownSS = dict()

    sys.stderr = sys.__stderr__

    # added to allow annotations not to be used.
    if len(list(juncs.keys())) < 1:
        print("No junctions from GTF or junctionsBed to correct with. Exiting...", file=sys.stderr)
        sys.exit(1)

    annotationFiles = dict()
    chrom2 = set()
    for chrom, data in juncs.items():
        if len(data) > 0:
            chrom2.add(chrom)
            annotationFiles[chrom] = os.path.join(tempDir, "%s_known_juncs.bed" % chrom)
            with open(os.path.join(tempDir, "%s_known_juncs.bed" % chrom), "w") as out:
                sortedData = sorted(list(data.keys()), key=lambda item: item[0])
                for k in sortedData:
                    annotation = data[k]
                    c1, c2, strand = k
                    print(chrom, c1, c2, annotation, ".", strand, sep="\t", file=out)
    return chrom2, annotationFiles


def correctsingleread(bedread, intervalTree, junctionBoundaryDict):
    juncs = bedread.juncs
    strand = bedread.strand
    c1Type, c2Type = ("donor", "acceptor") if strand == "+" else ("acceptor", "donor")

    newJuncs = list()
    ssStrands = set()
    novelSS = False

    for x in juncs:
        c1, c2 = x[0], x[1]
        if c1 not in junctionBoundaryDict:
            junctionBoundaryDict = ssCorrect(c1, strand, c1Type, intervalTree, junctionBoundaryDict, False)
        if c2 not in junctionBoundaryDict:
            junctionBoundaryDict = ssCorrect(c2, strand, c2Type, intervalTree, junctionBoundaryDict, False)

        # c1Obj, c2Obj = junctionBoundaryDict[c1], junctionBoundaryDict[c2] unused

        c1Corr = junctionBoundaryDict[c1].ssCorr.coord
        c2Corr = junctionBoundaryDict[c2].ssCorr.coord

        ssTypes = [junctionBoundaryDict[c1].ssCorr.ssType, junctionBoundaryDict[c2].ssCorr.ssType]

        ssStrands.add(junctionBoundaryDict[c1].ssCorr.strand)
        ssStrands.add(junctionBoundaryDict[c2].ssCorr.strand)

        if None in ssTypes:  # or ssTypes[0] == ssTypes[1]: # Either two donors or two acceptors or both none.
            novelSS = True
        newJuncs.append((c1Corr, c2Corr))

    blocks, sizes, starts = juncsToBed12(bedread.start, bedread.end, newJuncs)

    # 0 length exons, remove them.
    if min(sizes) == 0: novelSS = True

    if novelSS:
        return None
    else:
        bedread.juncs = newJuncs
        bedread.esizes = sizes
        bedread.estarts = starts
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
                if block[0] in {0, 7, 8}: hasmatch = True  # match
        # dirtowrite = '-' if is_reverse else '+'
        # chr1	476363	497259	ENST00000455464.7_ENSG00000237094.12	1000	-	476363	497259	0	3	582,169,151,	0,8676,20745,
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
        self.score = qualscore
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
        if self.strand == '-': exons = exons[::-1]
        exonseq = []
        for i in range(len(exons)):
            thisexonseq = genome.fetch(self.refchrom, exons[i][0], exons[i][1])
            if self.strand == '-': thisexonseq = revcomp(thisexonseq)
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

    def returndata(self):
        return self.juncstotranscript, self.junctogene, self.allannotse, self.genetoannotjuncs, self.transcripttoexons, self.alltranscripts


def getannotinfo(gtf, allregions):
    chromtoregions, regionstoannotdata = {}, {}
    for chrom, rstart, rend in allregions:
        if chrom not in chromtoregions: chromtoregions[chrom] = []
        chromtoregions[chrom].append((rstart, rend))
        regionstoannotdata[(chrom, rstart, rend)] = AnnotData()
    allchromtotranscripttoexons = {}
    for line in open(gtf):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]
            if ty == 'transcript' or ty == 'exon':
                if chrom not in allchromtotranscripttoexons: allchromtotranscripttoexons[chrom] = {}
                this_transcript = line[8][line[8].find('transcript_id') + 15:]
                this_transcript = this_transcript[:this_transcript.find('"')]
                this_gene = line[8].split('gene_id "')[1].split('"')[0].replace('_', '-')
                if (this_transcript, this_gene) not in allchromtotranscripttoexons[chrom]:
                    allchromtotranscripttoexons[chrom][(this_transcript, this_gene)] = [None, []]

                if ty == 'transcript':
                    allchromtotranscripttoexons[chrom][(this_transcript, this_gene)][0] = (start, end, strand)
                elif ty == 'exon':
                    allchromtotranscripttoexons[chrom][(this_transcript, this_gene)][1].append((start, end))
    for rchrom in allchromtotranscripttoexons:
        for transcript, gene in allchromtotranscripttoexons[rchrom]:
            tstart, tend, strand = allchromtotranscripttoexons[rchrom][(transcript, gene)][0]
            for rstart, rend in chromtoregions[rchrom]:
                if rstart < tstart < rend or rstart < tend < rend:
                    thisregion = (rchrom, rstart, rend)
                    sortedexons = sorted(allchromtotranscripttoexons[rchrom][(transcript, gene)][1])
                    regionstoannotdata[thisregion].transcripttoexons[(transcript, gene)] = tuple(sortedexons)
                    juncs = []
                    for i in range(len(sortedexons) - 1):
                        juncs.append((sortedexons[i][1], sortedexons[i + 1][0]))
                    regionstoannotdata[thisregion].alltranscripts.append((transcript, gene, strand))
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
    return regionstoannotdata  # juncstotranscript, junctogene, allannotse, genetoannotjuncs, transcripttoexons, chromtotranscript


def getcountsamcommand(args, refbed, outputname, mapfile, isannot):
    # count sam transcripts ; the dash at the end means STDIN
    count_cmd = ['count_sam_transcripts.py', '--sam', '-',
                 '-o', outputname, '-t', str(args.threads),
                 '--quality', str(args.quality), '-w', str(args.end_window)]
    if mapfile:
        count_cmd += ['--generate_map', mapfile]
    if args.stringent or isannot:
        count_cmd += ['--stringent']
    if args.check_splice:  # or isannot:
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
    ##'--split-prefix', 'minimap2transcriptomeindex', doesn't work with MD tag
    if type(inputreads) == str: inputreads = [inputreads]
    mm2_cmd = tuple(
        ['minimap2', '-a', '-t', str(args.threads), '-N', '4', '--MD'] + args.mm2_args + [alignfasta] + inputreads)
    ###FIXME add in step to filter out chimeric reads here
    ###FIXME really need to go in and check on how count_sam_transcripts is working
    count_cmd = getcountsamcommand(args, refbed, outputname, mapfile, isannot)
    # print(' '.join(mm2_cmd))
    # print(' '.join(count_cmd))
    # sys.stderr.write('Aligning and counting supporting reads for transcripts\n')
    pipettor.run([mm2_cmd, count_cmd])


def getbestends(currgroup):
    bestends = []
    for start1, end1, strand1, name1 in currgroup:
        score, weightedscore = 0, 0
        for start2, end2, strand2, name2 in currgroup:
            if abs(start1 - start2) <= args.end_window and abs(end1 - end2) <= args.end_window:
                score += 2
                weightedscore += ((args.end_window - abs(start1 - start2)) / args.end_window) + \
                                 ((args.end_window - abs(end1 - end2)) / args.end_window)
        bestends.append((weightedscore, start1, end1, strand1, name1))
    bestends.sort(reverse=True)
    ###DO I WANT TO ADD CORRECTION TO NEARBY ANNOTATED TSS/TTS????
    return bestends[0]


def combinefinalends(currgroup):
    if len(currgroup) == 1:
        return currgroup[0]
    else:
        currgroup.sort(key=lambda x: x[2])  ##sort by marker
        allreads = [y for x in currgroup for y in x[-1]]
        if currgroup[0] != 'a':  ##if no annotated iso, sort further
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
    if len(group) > 0: newgroups.append(group)
    return newgroups

def collapseendgroups(readends, dogetbestends=True):
    # print([x[:2] for x in readends])
    startgroups = groupreadsbyends(readends, 0, args.end_window)
    allendgroups, isoendgroups = [], []
    for startgroup in startgroups:
        allendgroups.extend(groupreadsbyends(startgroup, 1, args.end_window))
    for endgroup in allendgroups:
        if dogetbestends:
            isoendgroups.append(list(getbestends(endgroup)) + [[x[3] for x in endgroup]])
        else:
            isoendgroups.append(combinefinalends(endgroup))
    # print([x[:2] for x in isoendgroups])
    return isoendgroups




def addpresetargs(args):
    args.mm2_args = []
    args.end_window = 100
    args.trust_ends = False
    args.remove_internal_priming = False
    args.isoformtss = True
    return args


def identifygoodmatchtoannot(args, tempprefix, thischrom, annottranscripttoexons, alltranscripts, genome):
    goodaligntoannot, firstpasssingleexons, supannottranscripttojuncs = [], set(), {}
    # if args.transcriptfasta:
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
                if strand == '-': exons = exons[::-1]
                exonseq = []
                for i in range(len(exons)):
                    thisexonseq = genome.fetch(thischrom, exons[i][0], exons[i][1])
                    if strand == '-': thisexonseq = revcomp(thisexonseq)
                    exonseq.append(thisexonseq)
                annotbed.write('\t'.join([str(x) for x in bedline]) + '\n')
                annotfa.write('>' + transcript + '_' + gene + '\n')
                annotfa.write(''.join(exonseq) + '\n')
        transcriptomealignandcount(args, tempprefix + '.reads.fasta',
                                   tempprefix + '.annotated_transcripts.fa',  # args.transcriptfasta,
                                   tempprefix + '.annotated_transcripts.bed',  # args.annotated_bed,
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
    for read in samfile.fetch(rchrom, int(rstart), int(rend)):
        if not read.is_secondary and not read.is_supplementary:
            if read.query_name not in goodaligntoannot:
                shortchromfasta.write('>' + read.query_name + '\n')
                shortchromfasta.write(read.get_forward_sequence() + '\n')
                if read.mapping_quality >= args.quality:
                    juncstrand = inferMM2JuncStrand(read)
                    bedread = BedRead()
                    bedread.generate_from_cigar(read.reference_start, read.is_reverse, read.cigartuples,
                                                read.query_name,
                                                read.reference_name, read.mapping_quality, juncstrand)
                    correctedread = correctsingleread(bedread, intervalTree, junctionBoundaryDict)
                    if correctedread:
                        junckey = tuple(sorted(correctedread.juncs))
                        if junckey not in sjtoends: sjtoends[junckey] = []
                        sjtoends[junckey].append(
                            (correctedread.start, correctedread.end, correctedread.strand, correctedread.name))
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
        if end > lastend: lastend = end
    if len(thisgroup) > 0: furthergroups.append(thisgroup)
    return furthergroups


def filterendsbyredundantandsupport(args, goodendswithsupreads):
    bestends, nopass = [], []
    for i in range(len(goodendswithsupreads)):
        goodendswithsupreads[i][-1] = len(goodendswithsupreads[i][-1])  ##last val is read support
    goodendswithsupreads.sort(key=lambda x: [x[0], x[2] - x[1]],
                              reverse=True)  ##first by weighted score, then by length
    juncsupport = sum([x[-1] for x in goodendswithsupreads])
    if juncsupport >= args.support:
        if args.no_redundant == 'none':
            if goodendswithsupreads[0][-1] < args.support:
                goodendswithsupreads = [goodendswithsupreads[0]]
                goodendswithsupreads[0][-1] = juncsupport
            else:
                goodendswithsupreads = [x for x in goodendswithsupreads if x[-1] >= args.support]
                goodendswithsupreads = goodendswithsupreads[:args.max_ends]  ###select only top most supported ends
            for theseends in goodendswithsupreads:
                bestends.append(theseends)
        else:
            # if args.no_redundant == 'best_only': ###this uses the default sorting
            if args.no_redundant == 'longest':
                goodendswithsupreads.sort(reverse=True, key=lambda x: x[2] - x[1])
            # bestscore, beststart, bestend, beststrand, bestname, endsupport = goodendstosupreads[0]
            thisbest = goodendswithsupreads[0]
            thisbest[-1] = juncsupport  ####all reads for junc are counted as support
            bestends.append(thisbest)
    else:
        thisbest = goodendswithsupreads[0]
        thisbest[-1] = juncsupport
        nopass.append(thisbest)
    return bestends, nopass


def processjuncstofirstpassisos(args, tempprefix, thischrom, sjtoends, firstpasssingleexons):
    firstpassunfiltered, firstpassjunctoname = {}, {}
    with open(tempprefix + '.firstpass.unfiltered.bed', 'w') as isoout:
        for juncs in sjtoends:
            goodendswithsupreads = collapseendgroups(sjtoends[juncs])
            ###single exon - collapse overlapping intervals if need to pick best for each region
            if juncs == ():
                groupedgoodendswithsupreads = groupfirstpasssingleexon(goodendswithsupreads)
            else:
                groupedgoodendswithsupreads = [goodendswithsupreads]

            for goodendswithsupreads in groupedgoodendswithsupreads:
                bestends, nopass = filterendsbyredundantandsupport(args, goodendswithsupreads)
                for bestscore, beststart, bestend, beststrand, bestname, thisscore in bestends:
                    thisiso = BedRead()
                    thisiso.generatefromvals(thischrom, beststart, bestend, bestname, thisscore, beststrand, juncs)
                    firstpassunfiltered[bestname] = thisiso
                    isoout.write('\t'.join(thisiso.getbedline()) + '\n')
                    if juncs == ():
                        firstpasssingleexons.add((thisiso.exons[0][0], thisiso.exons[0][1], thisiso.name))
                    else:
                        for j in juncs:
                            if j not in firstpassjunctoname: firstpassjunctoname[j] = set()
                            firstpassjunctoname[j].add(bestname)
                        firstpasssingleexons.update(thisiso.exons)
                for bestscore, beststart, bestend, beststrand, bestname, thisscore in nopass:
                    thisiso = BedRead()
                    thisiso.generatefromvals(thischrom, beststart, bestend, bestname, thisscore, beststrand, juncs)
                    isoout.write('\t'.join(thisiso.getbedline()) + '\n')
    firstpasssingleexons = sorted(list(firstpasssingleexons))
    return firstpassunfiltered, firstpassjunctoname, firstpasssingleexons


def filtersplicediso(args, thisiso, firstpassjunctoname, firstpassunfiltered, junctogene, supannottranscripttojuncs,
                     annottranscripttoexons):
    isoswithsimilarjuncs = set()
    for j in thisiso.juncs:
        isoswithsimilarjuncs.update(firstpassjunctoname[j])
        if j in junctogene:
            isoswithsimilarjuncs.update(junctogene[j])  ##annot isos
    issubset = [0, 0]  # first exon is a subset, last exon is a subset
    firstexon, lastexon = thisiso.exons[0], thisiso.exons[-1]
    # firstexonlen, lastexonlen = firstexon[1]-firstexon[0], lastexon[1]-lastexon[0]
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
                if thisisostrjuncs in otherisostrjuncs:  ###if junctions are subset
                    ###check whether first + last exon overlap
                    for i in range(len(otherisoexons)):

                        otherexon = otherisoexons[i]
                        if i == 0 or i == len(otherisoexons) - 1:  ##is first or last exon of other transcript
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
    # if thisiso.name == 'm54284U_201026_025447/141558985/ccs': print(thisiso.name, issubset, thisiso.score, superset_support)
    if sum(issubset) < 2:  ###both first and last exon have to overlap
        return True  # firstpass[isoname] = thisiso
    elif args.filter != 'nosubset':
        if thisiso.score >= args.support and thisiso.score > max(superset_support) * 1.2:
            return True  # firstpass[isoname] = thisiso
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
                    break  ##filter out
                else:  ### is other single exon - check relative expression
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
    # firstpass[thisiso.name] = thisiso


def filtersingleexongroup(args, currgroup, firstpassunfiltered, firstpass):
    for groupediso in currgroup:
        if len(groupediso) == 3:  ##is single exon with name
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
        if end > lastend: lastend = end
    if len(currgroup) > 0: firstpass = filtersingleexongroup(args, currgroup, firstpassunfiltered, firstpass)

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
                    if passesfiltering: firstpass[isoname] = thisiso
        ###HANDLE SINGLE EXONS SEPARATELY - group first - one traversal of list
        firstpass = filterallsingleexon(args, firstpasssingleexons, firstpassunfiltered, firstpass)

    return firstpass


def combinetempfilesbysuffix(args, tempprefixes, suffixes):
    for filesuffix in suffixes:
        with open(args.output + filesuffix, 'wb') as allannotcounts:
            for tempprefix in tempprefixes:  # + ['notwellaligned']:
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
    for annotexoninfo in allannotse[index - 2:index + 2]:  ##start, end, strand, gene
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
                                  genetoannotjuncs, genome):
    with open(tempprefix + '.firstpass.bed', 'w') as isoout, open(tempprefix + '.firstpass.fa', 'w') as seqout:
        ####THIS IS WHERE WE CAN GET GENES AND ADJUST NAMES
        annotnametousedcounts = {}
        for isoname in firstpass:
            thisiso = firstpass[isoname]
            ####Adjust name based on annotation
            thistranscript, thisgene = thisiso.name, thischrom.replace('_', '-') + ':' + str(round(thisiso.start, -3))
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
                    gene_hits = getsingleexongenehits(thisiso, allannotse)  ##single exon

                if gene_hits:
                    sortedgenes = sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)
                    thisgene = sortedgenes[0][0]
            thisiso.name = thistranscript + '_' + thisgene
            isoout.write('\t'.join(thisiso.getbedline()) + '\n')
            seqout.write('>' + thisiso.name + '\n')
            seqout.write(thisiso.getsequence(genome) + '\n')


def getisogenefromname(isogene):
    iso = '_'.join(isogene.split('_')[:-1])
    gene = isogene.split('_')[-1]
    return iso, gene


def processdetectedisos(args, mapfile, bedfile, marker, genetojuncstoends):  # ogisotoreads, , genetoisoorigins):
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
                # genetoisoorigins[gene] = set()
            # genetoisoorigins[gene].add(marker)
            if junckey not in genetojuncstoends[gene]: genetojuncstoends[gene][junckey] = []
            genetojuncstoends[gene][junckey].append([start, end, isoid, ogisotoreads[isoid]])

    return genetojuncstoends  # , genetoisoorigins, ogisotoreads,


def revcomp(seq):
    compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    seq = seq.upper()
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq[::-1])


def getbedgtfoutfrominfo(endinfo, chrom, strand, juncs, gene, genome):
    start, end, isoid, readnames = endinfo
    score = len(readnames)
    marker, iso = isoid
    estarts, esizes = getexonsfromjuncs(juncs, start, end)
    bedline = [chrom, start, end, iso + '_' + gene, score, strand, start, end,
               getrgb(iso, strand, juncs), len(estarts), ','.join([str(x) for x in esizes]),
               ','.join([str(x) for x in estarts])]
    gtflines = []
    gtflines.append([chrom, 'FLAIR', 'transcript', start + 1, end, score, strand, '.',
                     'gene_id "' + gene + '"; transcript_id "' + iso + '";'])
    exons = [(start + estarts[i], start + estarts[i] + esizes[i]) for i in range(len(estarts))]
    if strand == '-': exons = exons[::-1]
    exonseq = []
    for i in range(len(exons)):
        gtflines.append([chrom, 'FLAIR', 'exon', exons[i][0] + 1, exons[i][1], score, strand, '.',
                         'gene_id "' + gene + '"; transcript_id "' + iso + '"; exon_number ' + str(i + 1)])
        thisexonseq = genome.fetch(chrom, exons[i][0], exons[i][1])
        if strand == '-': thisexonseq = revcomp(thisexonseq)
        exonseq.append(thisexonseq)
    return '\t'.join([str(x) for x in bedline]) + '\n', gtflines, ''.join(exonseq)


def combineannotnovelwriteout(args, genetojuncstoends, genome):
    with open(args.output + '.isoforms.bed', 'w') as isoout, open(args.output + '.read.map.txt', 'w') as mapout, open(
            args.output + '.isoforms.gtf', 'w') as gtfout, open(args.output + '.isoforms.fa', 'w') as seqout:
        for gene in genetojuncstoends:
            gtflines, tstarts, tends = [], [], []
            for chrom, strand, juncs in genetojuncstoends[gene]:
                endslist = genetojuncstoends[gene][(chrom, strand, juncs)]
                endslist = collapseendgroups(endslist, False)

                # start, end, isoid, ogisotoreads[isoid]
                ###TODO could try accounting for all reads assigned to isoforms - assign them to closest ends
                ###not sure how much of an issue this is
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
                    seqout.write('>' + iso + '_' + gene + '\n')
                    seqout.write(tseq + '\n')
            gtflines.insert(0, [chrom, 'FLAIR', 'gene', min(tstarts) + 1, max(tends), '.', gtflines[0][6], '.',
                                'gene_id "' + gene + '";'])
            for g in gtflines:
                gtfout.write('\t'.join([str(x) for x in g]) + '\n')


def runcollapsebychrom(listofargs):
    args, tempprefix, splicesiteannot_chrom, juncstotranscript, junctogene, allannotse, genetoannotjuncs, annottranscripttoexons, allannottranscripts = listofargs
    ###first extract reads for chrom as fasta
    tempsplit = tempprefix.split('/')[-1].split('-')
    rchrom, rstart, rend = '-'.join(tempsplit[:-2]), tempsplit[-2], tempsplit[-1]
    pipettor.run([('samtools', 'view', '-h', args.genomealignedbam, rchrom + ':' + rstart + '-' + rend),
                  ('samtools', 'fasta', '-')],
                 stdout=open(tempprefix + '.reads.fasta', 'w'))
    ###then align reads to transcriptome and run count_sam_transcripts
    genome = pysam.FastaFile(args.genome)
    goodaligntoannot, firstpasssingleexons, supannottranscripttojuncs = identifygoodmatchtoannot(args, tempprefix,
                                                                                                 rchrom,
                                                                                                 annottranscripttoexons,
                                                                                                 allannottranscripts,
                                                                                                 genome)
    ###load splice junctions for chrom
    intervalTree, junctionBoundaryDict = buildIntervalTree(splicesiteannot_chrom, args.ss_window, rchrom, False)

    print('processing reads into firstpass transcripts for', rchrom, rstart, rend, file=sys.stderr)
    samfile = pysam.AlignmentFile(args.genomealignedbam, 'rb')
    sjtoends = filtercorrectgroupreads(args, tempprefix, rchrom, rstart, rend, samfile, goodaligntoannot, intervalTree,
                                       junctionBoundaryDict)
    samfile.close()
    print('intron chains:', len(sjtoends.keys()), file=sys.stderr)

    firstpassunfiltered, firstpassjunctoname, firstpasssingleexons = processjuncstofirstpassisos(args, tempprefix,
                                                                                                 rchrom, sjtoends,
                                                                                                 firstpasssingleexons)

    firstpass = filterfirstpassisos(args, firstpassunfiltered, firstpassjunctoname, firstpasssingleexons,
                                    junctogene, supannottranscripttojuncs, annottranscripttoexons)
    print('firstpass isos: ', len(firstpass.keys()), file=sys.stderr)
    temptoremove = ['rm', tempprefix + '.reads.fasta', tempprefix + 'reads.notannotmatch.fasta']
    if not args.noaligntoannot and len(allannottranscripts) > 0:
        temptoremove.extend([tempprefix + '.annotated_transcripts.bed', tempprefix + '.annotated_transcripts.fa'])
    if len(firstpass.keys()) > 0:
        getgenenamesandwritefirstpass(tempprefix, rchrom, firstpass, juncstotranscript, junctogene,
                                      allannotse, genetoannotjuncs, genome)
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
    pipettor.run([tuple(temptoremove)])


def collapse(args):
    args = addpresetargs(args)
    genome = pysam.FastaFile(args.genome)
    allchrom = genome.references
    tempDir = makecorrecttempdir()
    sys.stderr.write('Getting regions\n')
    t1 = time.time()
    allregions = []
    if os.path.getsize(args.genomealignedbam) < 1e+9:  # less than 1G
        for chrom in allchrom:
            chromsize = genome.get_reference_length(chrom)
            allregions.append((chrom, 0, chromsize))
    else:
        pipettor.run([('bedtools', 'bamtobed', '-i', args.genomealignedbam),
                      ('bedPartition', '-minPartitionItems=1000', '-parallel=' + str(args.threads), '/dev/stdin',
                       tempDir + 'regions.bed')])
        for line in open(tempDir + 'regions.bed'):
            line = line.rstrip().split('\t')
            chrom, start, end = line[0], int(line[1]), int(line[2])
            allregions.append((chrom, start, end))
    sys.stderr.write('Number of regions ' + str(len(allregions)) + '\n')
    sys.stderr.write('Generating splice site database\n')
    knownchromosomes, annotationFiles = generateKnownSSDatabase(args, tempDir)

    regionstoannotdata = {}
    if args.gtf:
        sys.stderr.write('Extracting annotation from GTF\n')
        regionstoannotdata = getannotinfo(args.gtf, allregions)
    sys.stderr.write('splitting by chunk\n')
    chunkcmds = []
    tempprefixes = []
    for rchrom, rstart, rend in allregions:
        if rchrom in knownchromosomes:
            juncstotranscript, junctogene, allannotse, genetoannotjuncs, annottranscripttoexons, allannottranscripts = {}, {}, [], {}, {}, []
            if args.gtf:
                juncstotranscript, junctogene, allannotse, genetoannotjuncs, annottranscripttoexons, allannottranscripts = \
                regionstoannotdata[(rchrom, rstart, rend)].returndata()

            splicesiteannot_chrom = annotationFiles[rchrom]
            tempprefix = tempDir + '-'.join([rchrom, str(rstart), str(rend)])
            # runcollapsebychrom([args, tempprefix, splicesiteannot_chrom, juncstotranscript,
            #                    junctogene, allannotse, genetoannotjuncs, annottranscripttoexons,
            #                    allannottranscripts])
            chunkcmds.append([args, tempprefix, splicesiteannot_chrom, juncstotranscript,
                              junctogene, allannotse, genetoannotjuncs, annottranscripttoexons,
                              allannottranscripts])
            tempprefixes.append(tempprefix)
    t2 = time.time()
    print('region overhead', t2 - t1)
    sys.stderr.write('running by chunk\n')
    p = Pool(args.threads)
    childErrs = set()
    for i in p.imap(runcollapsebychrom, chunkcmds):
        childErrs.add(i)
    p.close()
    p.join()
    if len(childErrs) > 1:
        print(childErrs, file=sys.stderr)
        sys.exit(1)

    if not args.noaligntoannot:
        combinetempfilesbysuffix(args, tempprefixes,
                                 ['.matchannot.counts.tsv', '.matchannot.read.map.txt', '.matchannot.bed'])
    combinetempfilesbysuffix(args, tempprefixes,
                             ['.firstpass.unfiltered.bed', '.firstpass.bed', '.novelisos.counts.tsv',
                              '.novelisos.read.map.txt'])

    shutil.rmtree(tempDir)

    genetojuncstoends = processdetectedisos(args, args.output + '.matchannot.read.map.txt',
                                            args.output + '.matchannot.bed', 'a', {})  # , {})
    genetojuncstoends = processdetectedisos(args, args.output + '.novelisos.read.map.txt',
                                            args.output + '.firstpass.bed', 'n',
                                            genetojuncstoends)  # , genetoisoorigins)

    combineannotnovelwriteout(args, genetojuncstoends, genome)
    genome.close()


if __name__ == "__main__":
    args = getargs()
    collapse(args)

# export PATH="/private/groups/brookslab/cafelton/git-flair/flair/bin:/private/groups/brookslab/cafelton/git-flair/flair/src/flair:$PATH"
##cd /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test
##minimap2 -ax splice -t 12 GRCh38.chr22.genome.fa testflairstraightfrombam/sample1.badread5X.fastq | samtools view -hb - | samtools sort - > testflairstraightfrombam/sample1.badread5X.genomealigned.bam; samtools index testflairstraightfrombam/sample1.badread5X.genomealigned.bam
###python3 /private/groups/brookslab/cafelton/git-flair/flair/src/flair/flair_straightfrombam.py -g GRCh38.chr22.genome.fa -f gencode.v38.annotation.chr22.gtf -b testflairstraightfrombam/sample1.badread5X.genomealigned.bam -o testflairstraightfrombam/031225 --no_redundant longest

# time python3 /private/groups/brookslab/cafelton/git-flair/flair/src/flair/flair_straightfrombam.py -g /private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa -f /private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf -b WTC11.ENCFF370NFS.aligned.bam --transcriptfasta generate -o 031225_NFS_straightfrombam_longest --no_redundant longest; time python3 /private/groups/brookslab/cafelton/git-flair/flair/src/flair/flair_straightfrombam.py -g /private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa -f /private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf -b WTC11.ENCFF370NFS.aligned.bam --transcriptfasta generate -o 031225_NFS_straightfrombam_bestonly --no_redundant best_only; time python3 /private/groups/brookslab/cafelton/git-flair/flair/src/flair/flair_straightfrombam.py -g /private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa -f /private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf -b WTC11.ENCFF370NFS.aligned.bam --transcriptfasta generate -o 031225_NFS_straightfrombam_allends
###add short-read juncs + rerun
##1608504 1607899 487995 0


##flair correct -g GRCh38.chr22.genome.fa -f gencode.v38.annotation.chr22.gtf --shortread /private/groups/brookslab/cafelton/lrgasp-wtc11/WTC11_all.SJ.out.tab -q WTC11.ENCFF370NFS.chr22.genomealigned.bed -o wtc11-chr22-ogflair/031325
###flair collapse -g GRCh38.chr22.genome.fa -f gencode.v38.annotation.chr22.gtf --quality 0 --isoformtss --stringent --check_splice --annotation_reliant generate --no_gtf_end_adjustment --keep_intermediate -q wtc11-chr22-ogflair/031325_all_corrected.bed -r WTC11.ENCFF370NFS.chr22.genomealigned.fasta -o wtc11-chr22-ogflair/031325
##time python3 /private/groups/brookslab/cafelton/git-flair/flair/src/flair/flair_straightfrombam.py -g GRCh38.chr22.genome.fa -f gencode.v38.annotation.chr22.gtf --shortread /private/groups/brookslab/cafelton/lrgasp-wtc11/WTC11_all.SJ.out.tab --stringent --check_splice -b WTC11.ENCFF370NFS.chr22.genomealigned.bam -o wtc11-chr22-straightfrombam/040125

"""

python /private/groups/brookslab/cafelton/bin/lrgasp-challenge-1-evaluation-main/sqanti3_lrgasp.challenge1.py /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/wtc11-chr22-ogflair/031325.isoforms.gtf /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/gencode.v38.annotation.chr22.gtf /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/GRCh38.chr22.genome.fa --gtf --cage_peak /private/groups/brookslab/cafelton/lrgasp-wtc11/refTSS.human.bed --polyA_motif_list /private/groups/brookslab/cafelton/lrgasp-wtc11/polyA_list.txt --polyA_peak /private/groups/brookslab/cafelton/lrgasp-wtc11/WTC11_all_polyApeaks.bed -c /private/groups/brookslab/cafelton/lrgasp-wtc11/WTC11_all.SJ.out.tab -o 031725
cd ../../wtc11-chr22-straightfrombam/squanti
python /private/groups/brookslab/cafelton/bin/lrgasp-challenge-1-evaluation-main/sqanti3_lrgasp.challenge1.py /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/wtc11-chr22-straightfrombam/031325.isoforms.gtf /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/gencode.v38.annotation.chr22.gtf /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/GRCh38.chr22.genome.fa --gtf --cage_peak /private/groups/brookslab/cafelton/lrgasp-wtc11/refTSS.human.bed --polyA_motif_list /private/groups/brookslab/cafelton/lrgasp-wtc11/polyA_list.txt --polyA_peak /private/groups/brookslab/cafelton/lrgasp-wtc11/WTC11_all_polyApeaks.bed -c /private/groups/brookslab/cafelton/lrgasp-wtc11/WTC11_all.SJ.out.tab -o 031725

python /private/groups/brookslab/cafelton/bin/lrgasp-challenge-1-evaluation-main/sqanti3_lrgasp.challenge1.py /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/wtc11-chr22-straightfrombam/040125.isoforms.gtf /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/gencode.v38.annotation.chr22.gtf /private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/GRCh38.chr22.genome.fa --gtf --cage_peak /private/groups/brookslab/cafelton/lrgasp-wtc11/refTSS.human.bed --polyA_motif_list /private/groups/brookslab/cafelton/lrgasp-wtc11/polyA_list.txt --polyA_peak /private/groups/brookslab/cafelton/lrgasp-wtc11/WTC11_all_polyApeaks.bed -c /private/groups/brookslab/cafelton/lrgasp-wtc11/WTC11_all.SJ.out.tab -o 040125

"""
