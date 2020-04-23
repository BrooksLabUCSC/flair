#!/usr/bin/env python3

import sys, os, re, argparse, textwrap, subprocess, shutil, random
sys.path.append('/pod/home/jeltje/lib')
from Bio import SeqIO
from collections import defaultdict
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

This program corrects splices in (nanopore read) alignments (psl format), using a whole genome annotation and a junctions file (genepred format).
The junctions file should be derived from a short read alignment from the same sample, and can be created using
junctionsFromSam.py, which is based on 
https://raw.githubusercontent.com/anbrooks/juncBASE/master/preProcess_getASEventReadCounts.py

(NOTE TO SELF: junctions are reported but not used for any decisions yet. This could be optional.)
NOTE TO USER: If you don't have a junctions file, just use the genome annotation file as input for both the -j and -a parameters.

Before looking for splices, small gaps in the alignment are merged. The gap size can be changed; default is 30.
Splices are corrected if they are within a 'wiggle' distance of the intron start or intron end in the genome annotation.

For any novel splices, consensus sequences (GT/AG etc) are checked using the genome sequence. If no consensus is found, output junctions have '.' in the strand field.

The mRnaToGenes program must be in $PATH

OUTPUTS:
    novel.txt        contains a list of novel junctions with enough supporting reads
    novelsplices.bed contains the novel junctions in bed format
    junctions.bed    contains all junctions found in the query file. The score field contains the number of alignments with this junction.
    notfound.txt     is a list of query donor and acceptor sites that could not be found within the allowed wiggle distance
    multihit.txt     is a list of query donor and acceptor sites that had multiple hits within the allowed wiggle distance
    corrected.gp     contains all query annotations, splice corrected where possible

Notes: Novel junctions are only reported if they are identical in at least 3 annotations (by default). This means that it is possible to miss junctions for which alignments are close but not identical. To see all novel junctions, set --novelthreshold to 1


        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('-a','--annotations', type=str, required=True, help="genepred format genome annotation")
group.add_argument('-j','--junctions', type=str, required=True, help="genepred format junctions from SAM file")
group.add_argument('-g','--genofasta', type=str, required=True, help="genome fasta file")
group.add_argument('-q','--query', type=str, required=True, help="psl format alignment to be corrected")
group.add_argument('-w', '--wiggle', type=int, required = True, help="wiggle room for annotated splice junction")
group.add_argument('-o', '--outdir', type=str, default = '.', help="output directory")
parser.add_argument('-m', '--mergesize', type=int, default=30,  help="merge genome alignment gaps of this size (30)")
parser.add_argument('-n', '--novelthreshold', type=int, default=3,  help="report any novel junctions that are confirmed by at least this many reads")
parser.add_argument('-r', '--reportnovel', dest='reportnovel', action='store_true', help='report novel junctions with n supporting reads even if the novel junction is within wiggle room of an annotated junction')
parser.add_argument('-mode', '--ambiguousmode', dest='expand', type=str, default='leave',  \
    help="[leave/correct/expand/closest] leave: leave splice site uncorrected if ambiguous, correct: correct splice site to ambiguous site with probability proportional to the frequency of that site's usage, expand: report all possible combinations of corrections for ambiguous reads, closest: correct ambiguous splice site to the closer annotated site")
parser.set_defaults(expand=False)

class GpHit(object):
    """
    holds one genepred format alignment; extracts intron locations
    """
    def __init__(self, inline):
        fields = inline.split('\t')
        # self.printBase = ('\t').join(fields[:8])
        self.fields = fields
        self.qName, self.chrom, self.strand = map(str, fields[:3])
        self.genoStart, self.genoEnd, self.cdsStart, self.cdsEnd, self.blockCount = map(int, fields[3:8])
        self.chromStarts = fields[8].strip(',').split(',')
        self.chromEnds = fields[9].strip(',').split(',')
        self.lefts = map(int, self.chromEnds[:-1])
        self.rights = map(int, self.chromStarts[1:])
        self.makeIntrons() 

    def makeIntrons(self):
        self.introns = []
        for i in xrange(len(self.lefts)):
            intron = ('{}-{}').format(self.lefts[i], self.rights[i])
            self.introns.append(intron)

    def doPrint(self, ohandle):
        """
        Print current alignment in genePred format
        """
        self.printBase = '\t'.join([self.qName] + self.fields[1:8])
        if self.blockCount == 1:
            ohandle.write('{}\t{},\t{},\n'.format(self.printBase, self.chromStarts[0], self.chromEnds[0]))
            return
        chromStarts = str(self.genoStart)
        for i in self.rights:
            chromStarts += ',' + str(i)
        chromEnds = '' 
        for i in self.lefts:
            chromEnds += str(i) + ','
        chromEnds += str(self.genoEnd)
        ohandle.write('{}\t{},\t{},\n'.format(self.printBase, chromStarts, chromEnds))
        
def die_roll(olist):
    """ 
    olist is a list of tuples, each tuple containing the outcome name and the frequency of that outcome
    """
    total = sum([tup[1] for tup in olist])
    r = random.randint(0, total)
    count = 0
    outcome = olist[0][0]
    for t in olist:
        count += t[1]
        outcome = t[0]
        if count > r:
            return outcome
    return outcome


def correctCoord(chrom, coord, jtree, wiggle, outnf, outmh, side='left', anntree='', chr_seq='', strand='.'):
    """
    Compare coord (integer) to intervals in jtree. If exact integer is found, return. If integer falls within range, 
    return middle of range. If integer falls within multiple ranges, return original integer and print warning.
    if integer is not found, return original and print warning

    If there are multiple hits coordinates within the wiggle room, len(found) > 1. Depending on args.expand, this
    function will return either all the found coordinates or return a single coordinate, picked by a weighted die if
    the coordinates are from short read junction data based on the coordinate's frequency in the short read data.

    If there are no coordinates supported by short read junctions, then correctcoord will search within annotated
    junctions if provided. 
    """
    global correct
    global incorrect
    global unable
    global novels
    global uncorrected
    perfect, found = overlaps(jtree, coord, wiggle)
    if perfect:
        correct += 1
        if args.expand == 'expand':
            return [coord]
        return coord
    if len(found) >1:
        outmh.write("{}:{}\n".format(chrom, coord))
        olist = []  # list of possible coordinate outcomes
        founddist = []  # list of distances from coord to annotated coords
        if args.expand == 'expand':
            return found
        for f in found:
            founddist += [(f, abs(f - coord))]
            if anntree != '':  # search within short read supported junctions
                if side == 'left' and f in jfreq_lefts:  # short read support
                    olist += [(f, jfreq_lefts[f])]
                elif side == 'right' and f in jfreq_rights:
                    olist += [(f, jfreq_rights[f])]
            else:
                if side == 'left' and f in annot_lefts:
                    olist += [(f, annot_lefts[f])]
                elif side == 'right' and f in annot_rights:
                    olist += [(f, annot_rights[f])]
        incorrect += 1
        if args.expand == 'expand':
            return [olist]
        elif args.expand == 'leave':
            return coord
        elif args.expand == 'closest':
            founddist.sort(key = lambda x: x[1])
            if founddist[0][1] != founddist[1][1]:  # if the coord is closer to one of the found sites
                return founddist[0][0]
            return die_roll(founddist[:2])
        return die_roll(olist)  # 'correct'
    if len(found) == 0:
        if anntree != '':
            return correctCoord(chrom, coord, anntree, wiggle, outnf, outmh, side, '', chr_seq)  # search annotated junctions (anntree as jtree)
        if args.expand == 'closest':
            if novelValidity(chr_seq, coord, side,strand):
                novels += 1
                return coord
            unable += 1
            return False
        outnf.write("{}:{}\n".format(chrom, coord))
        uncorrected += 1
        unable += 1
        if args.expand == 'expand':
            return [coord]
        return coord
    if args.expand =='expand':
        return [found[0]]
    return found[0]



class Junctions():
    """
        Parse genepred format junctions 
    """
    def __init__(self):
        self.lefts = IntervalTree()
        self.rights = IntervalTree()
        self.introns = set()
        self.istrands = dict()
    def add(self, line, wiggle, jfreq_lefts, jfreq_rights):
        """
           Add starts and ends as ranges in an intervaltree (faster for querying)
        """
        fields = line.strip().split('\t');

        chromStarts = fields[8].strip(',').split(',')
        rightset = map(int, chromStarts[1:])
        chromEnds = fields[9].strip(',').split(',')
        leftset = map(int, chromEnds[:-1])
        strand = fields[2]
        for i in xrange(len(leftset)):
            intron = ('{}-{}').format(leftset[i], rightset[i])
            self.iadd(strand, intron)
            self.lefts[leftset[i] - wiggle: leftset[i] + wiggle +1] = intron
            self.rights[rightset[i] - wiggle: rightset[i] + wiggle +1] = intron
            self.introns.add(intron)
            if len(fields) < 11:
                annot_lefts[leftset[i]] = 1
                annot_rights[rightset[i]] = 1
                continue
            jfreq_lefts[leftset[i]] = int(fields[10])
            jfreq_rights[rightset[i]] = int(fields[10])

    def iadd(self, strand, intron):
        """
        Keep track of intron associated strand
        """
        if intron in self.istrands:
            if self.istrands[intron] != strand:
                # print >>sys.stderr, "WARNING, strand ambiguity for annotated junction", intron,  self.istrands[intron], strand
                self.istrands[intron] = '.'
        else:
            self.istrands[intron] = strand


def splitByChr(genepred, outdir):
    """
    Split input genepred file into one file per chromosome, named chrN.gp
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    chrList = []
    curChr = None
    foundchrs = set() 
    
    with open(genepred, 'r') as gp:
        for line in gp:
            chrom =line.split('\t')[1]
            if not curChr == chrom:
                if curChr:
                    with open(os.path.join(outdir, curChr + '.gp'), 'a') as o:
                        o.write(''.join(chrList))
                chrList = []
                curChr = chrom
                foundchrs.add(chrom)
            chrList.append(line)  
        with open(os.path.join(outdir, curChr + '.gp'), 'a') as o:
            o.write(''.join(chrList))
    return foundchrs

def overlaps(tree, coord, wiggle):
   """
    Returns the junction(s) that overlap
   """
   found = set()
   jcts = tree.search(coord)
   for i in jcts:
        if i.begin + wiggle == coord:
            return True, found
        found.add(i.begin + wiggle)  # why is this here..... ?
   return False, list(found)

def which(prog):
    """See if a program exists on the system and return the path"""
    cmd = ["which",prog]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0: 
        print >>sys.stderr, 'ERROR, cannot find program {}, please install and/or add to $PATH'.format(prog)
        sys.exit(1)
    return res

def listCount(mylist, wiggle, cutoff):
    """
    Count occurrences of coordinate in list within wiggle distance, print if more than cutoff occurrences
    """
    mylist = sorted(mylist)
    cur = 0
    hits = []
    for i in mylist:
        if cur + wiggle >= i:
            hits.append(i)
        else:
            if len(hits) >= cutoff:
                print >>sys.stderr, hits
            hits = []
        cur = i

def makeBed(chr, chr_seq, intron, score, outfile, strandlist):
    """
    Create bed format intron output with a blocksize of 20
    If the intron is not found in existing annotation, or if strand info conflicts, 
    infers strand from intron start and end sequence 
    If the strand cannot be determined, the strand field is output as a period
    """
    start, end = map(int, intron.split('-'))
    strand = '.'
    if intron in strandlist:
        strand = strandlist[intron]
    if strand == '.':
        strand = disambiguateJcnStr(chr_seq, start, end)

    cstart = start -20
    cend = end +20
    bedline = '{chr}\t{cstart}\t{cend}\t{chr}:{intron}\t{score}\t{strand}\t{cstart}\t{cend}\t0,0,255\t2\t20,20\t0,{bstart},\n'.format(chr=chr, 
         strand=strand, cstart=cstart, cend=cend, intron=intron, score=score, bstart=end-cstart)
    outfile.write(bedline)
    return bedline


def disambiguateJcnStr(chr_seq, start, end):
    """
    Will use splice site sequence to infer strand
    If no strand can be determined, returns '.'
    This function is from  disambiguate_junctions.py by Angela Brooks
    """

    strand = '.'

    # extract the sequence from the chromosome
    intron_seq = chr_seq[start-1:end].upper()

    if intron_seq.startswith("GT") and intron_seq.endswith("AG"):
        strand = "+"
    elif intron_seq.startswith("CT") and intron_seq.endswith("AC"):
        strand = "-"
    # Other common splice site sequence
    elif intron_seq.startswith("GC") and intron_seq.endswith("AG"):
        strand = "+"
    elif intron_seq.startswith("CT") and intron_seq.endswith("GC"):
        strand = "-"
    # minor spliceosome
    elif intron_seq.startswith("AT") and intron_seq.endswith("AC"):
        strand = "+"
    elif intron_seq.startswith("GT") and intron_seq.endswith("AT"):
        strand = "-"
    # Priority to 5' splice site since there is more information
    elif intron_seq.startswith("GT"):
        strand = "+"
    elif intron_seq.endswith("AC"):
        strand = "-"
    elif intron_seq.endswith("AG"):
        strand = "+"
    elif intron_seq.startswith("CT"):
        strand = "-"
    # else:
        # print >>sys.stderr, "Cannot find strand for {}-{}".format(start, end)
    return strand

def novelValidity(chr_seq, start,side, strand,toprint='', chrom=''):
    """
    Will use splice site sequence to infer validity of splice motif
    """
    # extract the sequence from the chromosome
    intron_seq = chr_seq[start-2:start+2].upper()
    if side == 'left':
        if strand == '+' and intron_seq.endswith("GT"):
            return True
        elif strand == '-' and intron_seq.endswith('CT'):
            return True
        elif strand == '.' and (intron_seq.endswith("GT") or intron_seq.endswith("CT")):
            return True
        return False
    if strand == '+' and intron_seq.startswith("AG"):
        return True
    elif strand == '-' and intron_seq.startswith('AC'):
        return True
    elif strand == '.' and (intron_seq.endswith("AG") or intron_seq.startswith("AC")):
        return True
    return False

def findChr(records, this_chr):
    """
    Will identify chr sequence, then call disambiguateJcnStr
    This function is from  disambiguate_junctions.py by Angela Brooks
    """
    chr_seq = None
    if this_chr in records:
        chr_seq = records[this_chr].seq
    elif this_chr.lstrip("chr") in records:
        chr_seq = records[unFormatChr(this_chr)].seq
    else:
        print "Cannot find %s in genome sequence." % this_chr
        return 'BAD'
    return chr_seq

# Main
args = parser.parse_args()

# Sanity check
if args.mergesize <= args.wiggle:
    print >>sys.stderr, "ERROR: mergesize (currently {}) must be larger than wiggle (currently: {}), please correct and rerun".format(args.mergesize, args.wiggle)
    sys.exit(1)

# Setup temporary output location
tmpdir = '/tmp/spliceCorr_' + str(os.getpid())
if os.path.exists(tmpdir):
    shutil.rmtree(tmpdir)
os.makedirs(tmpdir)

# Convert the genome alignment psl to genepred - doing this first because it is most likely to fail
mrnaToGene = which('mrnaToGene')
msize = '-insertMergeSize={}'.format(args.mergesize)
queryname = args.query[args.query.rfind('/')+1:]
gpgeno = '{}/{}.gp'.format(tmpdir, queryname)
cmd = [mrnaToGene, '-noCds', msize, args.query, gpgeno]
try:
    subprocess.check_call(cmd)
except subprocess.CalledProcessError:
    print >>sys.stderr, "Please have mrnaToGene in your path, exiting."
    sys.exit(1)

# split all input files into chromosomes
workdir = args.outdir
if not os.path.exists(workdir):
    os.makedirs(workdir)
jdir = os.path.join(tmpdir, 'junctions')
adir = os.path.join(tmpdir, 'annots')
qdir = os.path.join(tmpdir, 'aligns')
splitByChr(args.junctions, jdir)
splitByChr(args.annotations, adir)
chroms = splitByChr(gpgeno, qdir)

# read genome fasta
try:
    records = SeqIO.index(args.genofasta, "fasta")
except:
    print >>sys.stderr, "ERROR: Could not open genome sequence."
    sys.exit(1)

# Per-chromosome analysis
expanded = {}
uncorrected = 0
correct, novels = 0, 0
incorrect, unable = 0, 0
with open(os.path.join(workdir, 'corrected.gp'), 'w') as outgp, \
  open(os.path.join(workdir, 'novelsplices.bed'), 'w') as splicebed, \
  open(os.path.join(workdir, 'junctions.bed'), 'w') as outbed, \
  open(os.path.join(workdir, 'novel.txt'), 'w') as outnv, \
  open(os.path.join(workdir, 'novel_uncorrected.txt'),'w') as outnvu, \
  open(os.path.join(workdir, 'notfound.txt'), 'w') as outnf, \
  open(os.path.join(workdir, 'multihit.txt'), 'w') as outmh:
    outnv.write("junction\talignment\tRNASeq\tannotation\n")
    for c in chroms:
        print >>sys.stderr, c
        cors_ann = Junctions()
        cors_short = Junctions()
        junctionIntrons = []
        annotIntrons = []
        jfreq_lefts = {}  # junction frequencies for smaller coordinate of a junction from short read
        jfreq_rights = {}
        annot_lefts = {}
        annot_rights = {}
        # read in junctions file
        try:
            temp = open(os.path.join(jdir, c + '.gp'), 'r')
            temp = open(os.path.join(adir, c + '.gp'), 'r')
        except:
            # print >>sys.stderr, "Chrom not found in junction file."
            continue
        with open(os.path.join(jdir, c + '.gp'), 'r') as f:
            for line in f:
                cors_short.add(line, args.wiggle, jfreq_lefts, jfreq_rights)
                hit = GpHit(line.strip())
                junctionIntrons.extend(hit.introns)
        # add the annotation file
        with open(os.path.join(adir, c + '.gp'), 'r') as f:
            for line in f:
                cors_ann.add(line, args.wiggle, annot_lefts, annot_rights)
                hit = GpHit(line.strip())
                annotIntrons.extend(hit.introns)
        # read alignment file and correct intronstarts and intron ends
        if args.annotations == args.junctions:
            jfreq_lefts = annot_lefts
            jfreq_rights = annot_rights
        expansions_lefts = []
        expansions_rights = []
        alignIntrons = []
        allIntrons = []
        chr_seq = findChr(records, c) 
        if chr_seq == 'BAD':  # i also added this - alison
            continue  

        with open(os.path.join(qdir, c + '.gp'), 'r') as f:
            for line in f:
                hit = GpHit(line.strip())
                newlefts = []
                expansions_lefts = []  # a list of lists of newlefts
                inferredstrand = hit.strand
                if not hit.lefts:
                    hit.doPrint(outgp)
                    continue
                for i in hit.lefts:
                    cc = correctCoord(c, i, cors_short.lefts, args.wiggle, outnf, outmh, 'left', cors_ann.lefts, chr_seq, inferredstrand)
                    if args.expand == 'expand':
                        if not expansions_lefts:
                            expansions_lefts = [[e] for e in cc]
                            continue
                        temp = []
                        for left in cc:
                            temp += [e + [left] for e in expansions_lefts]
                        expansions_lefts = temp
                    else:
                        if not cc:
                            break
                        newlefts.append(cc)
                if not cc:
                    continue
                newrights = []
                expansions_rights = []
                for i in hit.rights:
                    cc = correctCoord(c, i, cors_short.rights, args.wiggle, outnf, outmh, 'right', cors_ann.rights, chr_seq, inferredstrand)
                    if args.expand == 'expand':
                        if not expansions_rights:
                            expansions_rights = [[e] for e in cc]
                            continue
                        temp = []
                        for right in cc:
                            temp += [e + [right] for e in expansions_rights]
                        expansions_rights = temp
                    else:
                        if not cc:
                            break
                        newrights.append(cc)
                if not cc:
                    continue

                if args.reportnovel:
                    hit.makeIntrons()  # all uncorrected introns
                    allIntrons.extend(hit.introns)

                hit.lefts = newlefts
                hit.rights = newrights

                qindex = 0
                if args.expand == 'expand':
                    hit.qName += '_' + str(qindex)
                    for leftset in expansions_lefts:
                        for rightset in expansions_rights:
                            hit.lefts = leftset
                            hit.rights = rightset
                            hit.makeIntrons()
                            alignIntrons.extend(hit.introns)
                            hit.doPrint(outgp)
                            qindex += 1
                            hit.qName = hit.qName[:hit.qName.rfind('_')] + str(qindex)
                    if qindex not in expanded:
                        expanded[qindex] = 0
                    expanded[qindex] += 1
                else:
                    hit.makeIntrons()
                    alignIntrons.extend(hit.introns)
                    hit.doPrint(outgp)
        

        # now for every unique junction in the alignment, count the number of supporting reads and check if the junction is present in the annotation and in the junctions file
        # Also print out every junction in bed format and check for consensus splice sites.
        for i in (set(alignIntrons)):
            alcount = alignIntrons.count(i) 
            jucount = junctionIntrons.count(i)
            ancount = annotIntrons.count(i)
            bedline = makeBed(c, chr_seq, i, alcount, outbed, cors_short.istrands)
            if not ancount and alcount >= args.novelthreshold:
               outnv.write("{}:{}\t{}\t{}\t{}\n".format(c, i, alignIntrons.count(i), junctionIntrons.count(i), annotIntrons.count(i)))
               splicebed.write(bedline)

        if args.reportnovel:
            for i in set(allIntrons):
                alcount = allIntrons.count(i) 
                jucount = junctionIntrons.count(i)
                ancount = annotIntrons.count(i)
                bedline = makeBed(c, chr_seq, i, alcount, outbed, cors_short.istrands)
                if not ancount and alcount >= args.novelthreshold:
                   outnvu.write("{}:{}\t{}\t{}\t{}\n".format(c, i, allIntrons.count(i), junctionIntrons.count(i), annotIntrons.count(i)))
                   splicebed.write(bedline)

if args.expand == 'expand':
    print(expanded)
# sys.stderr.write('# splice sites with no short read/gtf annotation:{}\n'.format(uncorrected))
shutil.rmtree(tmpdir)

# print('perfect splice sites: {}, corrected splice sites: {}, novel splice support: {}, unable to be corrected: {}'.format(correct, incorrect, novels,unable))

# print textwrap.dedent('''\


# SUCCESS

# These output files were created:
#     novel.txt        contains a list of novel junctions with enough supporting reads
#     novelsplices.bed contains novel junctions
#     junctions.bed    contains all junctions found in the query file. The score field contains the number of alignments with this junction.
#     notfound.txt     is a list of query donor and acceptor sites that could not be found within the allowed wiggle distance
#     multihit.txt     is a list of query donor and acceptor sites that had multiple hits within the allowed wiggle distance
#     corrected.gp     contains all query annotations, splice corrected where possible

# Notes: 
#    The junctions in junctions.bed are NOT NECESSARILY correct splices. They just represent all remaining alignment gaps after splice correction.
#    Novel junctions are only reported if they are identical in at least 3 annotations. This means that it is possible to miss junctions for which alignments are close but not identical. To see all novel junctions, set --novelthreshold to 1
#    If the bed file contains a period in the strand field, the junction is likely an artifact (for example a misalignment of exon ends or an unfilled gap in an exon)

# PROGRAM FINISHED

#         ''')

