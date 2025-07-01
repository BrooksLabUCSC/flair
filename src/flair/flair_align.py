#! /usr/bin/env python3

import sys
import argparse
import os
import pipettor
import pysam
from flair import remove_internal_priming

def intronChainToestarts(ichain, start, end):
    esizes, estarts = [], [0,]
    for i in ichain:
        esizes.append(i[0] - (start + estarts[-1]))
        estarts.append(i[1] - start)
    esizes.append(end - (start + estarts[-1]))
    return esizes, estarts


def inferMM2JuncStrand(read):
    # minimap gives junction strand denoted as 'ts'
    # the sign corresponds to the alignment orientation, where + agrees and - disagrees
    orientation = read.flag
    try:
        juncDir = read.get_tag('ts')
    except:
        juncDir = None

    # Try to resolve strand by looking for polyA
    if not juncDir:
        left, right = read.cigar[0], read.cigar[-1]
        s1, s2 = read.seq[:50], read.seq[-50:]
        # pa = str()
        if ("T" * 10 in s1 and left[0] == 4 and left[1] >= 10) and (
                        "A" * 10 in s2 and right[0] == 4 and right[1] >= 10):
            # probably internal priming
            juncDir = "ambig"

        elif ("T" * 10 in s1 and left[0] == 4 and left[1] >= 10):
            # maps to positive strand but has a rev comp polyA
            juncDir = "-" if orientation == 16 else "+"
        # print("anti")
        # pa = "ppa"
        elif ("A" * 10 in s2 and right[0] == 4 and right[1] >= 10):
            # maps to positive strand but has a sense polyA
            juncDir = "+" if orientation == 16 else "-"
        # print("sense")
        # pa = "ppa"
        else:
            # no polyA or polyT. Fragment?
            juncDir = "ambig"
        # pa = "nan"

    else:
        if orientation == 0 and juncDir == "+":
            juncDir = "+"
        elif orientation == 0 and juncDir == "-":
            juncDir = "-"
        elif orientation == 16 and juncDir == "+":
            juncDir = "-"
        elif orientation == 16 and juncDir == "-":
            juncDir = "+"
    return juncDir

def bed_from_cigar(alignstart, is_reverse, cigartuples, readname, referencename, qualscore, juncDirection):
    positiveTxn = "27,158,119"  # green
    negativeTxn = "217,95,2"  # orange
    unknownTxn = "99,99,99"
    refpos = alignstart
    intronblocks = []
    hasmatch = False
    for block in cigartuples:
        if block[0] == 3 and hasmatch: #intron, pay attention
            intronblocks.append([refpos, refpos + block[1]])
            refpos += block[1]
        elif block[0] in {0, 7, 8, 2}:  # consumes reference
            refpos += block[1]
            if block[0] in {0, 7, 8}: hasmatch = True#match
    # dirtowrite = '-' if is_reverse else '+'
    #chr1   476363  497259  ENST00000455464.7_ENSG00000237094.12    1000    -       476363  497259  0       3       582,169,151,    0,8676,20745,
    esizes, estarts = intronChainToestarts(intronblocks,alignstart, refpos)
    rgbcolor = unknownTxn
    if juncDirection == "+": rgbcolor = positiveTxn
    elif juncDirection == "-": rgbcolor = negativeTxn
    else:
        juncDirection = "-" if is_reverse else "+"
    outline = [referencename, str(alignstart), str(refpos), readname, str(qualscore), juncDirection, str(alignstart), str(refpos), rgbcolor, str(len(intronblocks) + 1), ','.join([str(x) for x in esizes]) + ',', ','.join([str(x) for x in estarts]) + ',']
    return outline

def doalignment(args):
    # minimap
    mm2_cmd = ['minimap2', '-ax', 'splice', '-s', str(args.minfragmentsize),
               '-G', args.maxintronlen, '-t', str(args.threads)]
    if args.nvrna:
        mm2_cmd += ['-uf', '-k14']
    if args.junction_bed:
        mm2_cmd += ['--junc-bed', args.junction_bed]
    mm2_cmd += ['-secondary=no']
    if args.mm_index:
        mm2_cmd += [args.mm_index]
    else:
        mm2_cmd += [args.genome]
    mm2_cmd += args.reads

    samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', args.output + '_unfiltered.bam', '-')
    samtools_index_cmd = ('samtools', 'index', args.output + '_unfiltered.bam')
    pipettor.run([mm2_cmd, samtools_sort_cmd])
    pipettor.run([samtools_index_cmd])

def dofiltering(args, inbam, filterreadmap=None):
    samfile = pysam.AlignmentFile(inbam, 'rb')#args.output + '_unfiltered.bam', "rb")
    withsup = None
    if args.filtertype == 'separate':
        withsup = pysam.AlignmentFile(args.output + '_chimeric.bam', "wb", template=samfile)
    trash = pysam.AlignmentFile(args.output + '_trash.bam', "wb", template=samfile)
    fullyaligned = pysam.AlignmentFile(args.output + '.bam', "wb", template=samfile)
    outbed = open(args.output + '.bed', 'w')
    readstoremove = set()
    annottranscriptends, genome = None, None
    if args.remove_internal_priming:
        annottranscriptends = remove_internal_priming.getannotends(args.gtf)
        genome = pysam.FastaFile(args.genome)
    if filterreadmap:
        filteredout = pysam.AlignmentFile(args.output + '_filteredout.bam', "wb", template=samfile)
        for line in open(filterreadmap):
            reads = line.rstrip().split('\t', 1)[1].split(',')
            readstoremove.update(set(reads))
    for read in samfile.fetch():
        if read.is_mapped and not read.is_secondary:
            mapq = read.mapping_quality
            if mapq < args.quality:
                trash.write(read)
                continue
            if filterreadmap:
                if read.query_name in readstoremove:
                    filteredout.write(read)
                    continue
            if args.remove_singleexon:
                if read.get_cigar_stats()[0][3] == 0:
                    trash.write(read)
                    continue
            if args.remove_internal_priming:
                # print(read.query_name)
                notinternalpriming = remove_internal_priming.removeinternalpriming(read.reference_name, read.reference_start, read.reference_end, read.is_reverse,
                                                          genome, annottranscriptends, None, args.intprimingthreshold, args.intprimingfracAs)
                # print(read.query_name, notinternalpriming)
                if not notinternalpriming:
                    trash.write(read)
                    continue
            if read.has_tag('SA'):
                if args.filtertype == 'separate':
                    withsup.write(read)
                elif args.filtertype == 'removesup':
                    trash.write(read)
                elif read.is_supplementary:
                    trash.write(read)
                else:
                    fullyaligned.write(read)
                    juncstrand = inferMM2JuncStrand(read)
                    bedline = bed_from_cigar(read.reference_start, read.is_reverse, read.cigartuples,
                                                                     read.query_name, read.reference_name, mapq, juncstrand)
                    outbed.write('\t'.join(bedline) + '\n')
            else:
                fullyaligned.write(read)
                juncstrand = inferMM2JuncStrand(read)
                bedline = bed_from_cigar(read.reference_start, read.is_reverse, read.cigartuples, read.query_name,
                                                                 read.reference_name, mapq, juncstrand)
                outbed.write('\t'.join(bedline) + '\n')
        else:
            trash.write(read)
    samfile.close()
    if filterreadmap:
        filteredout.close()
        pysam.index(args.output + '_filteredout.bam')
    if args.filtertype == 'separate':
        withsup.close()
        pysam.index(args.output + '_chimeric.bam')
    trash.close()
    fullyaligned.close()
    outbed.close()
    pysam.index(args.output + '_trash.bam')
    pysam.index(args.output + '.bam')


# note to self: select interpreter for conda
def align():
    parser = argparse.ArgumentParser(description='flair-align parse options',
            usage='flair_align -g genome.fa -r <reads.fq>|<reads.fa> [options]')
    required = parser.add_argument_group('required named arguments')
    atleastone = parser.add_argument_group('Either one of the following arguments is required')
    required.add_argument('-r', '--reads', nargs='+', type=str, required=True,
            help='FastA/FastQ files of raw reads')
    atleastone.add_argument('-g', '--genome', type=str,
            help='FastA of reference genome, can be minimap2 indexed')
    atleastone.add_argument('--mm_index', type=str, default='',
            help='minimap2 index .mmi file')
    parser.add_argument('-o', '--output', default='flair.aligned',
            help='output file name base (default: flair.aligned)')
    parser.add_argument('-t', '--threads', type=int, default=4,
            help='minimap2 number of threads (4)')
    parser.add_argument('--junction_bed', default='',
            help='annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment')
    parser.add_argument('--nvrna', action='store_true', default=False,
            help='specify this flag to use native-RNA specific alignment parameters for minimap2')
    parser.add_argument('--quality', type=int, default=0,
            help='minimum MAPQ of read alignment to the genome (0)')
    parser.add_argument('--minfragmentsize', type=int, default=80,
                                            help='minimum size of alignment kept, used in minimap -s. More important when doing downstream fusion detection')
    parser.add_argument('--maxintronlen', default='200k',
                                            help='maximum intron length in genomic alignment. Longer can help recover more novel isoforms with long introns')
    parser.add_argument('--filtertype', type=str, default='removesup',
    help='method of filtering chimeric alignments (potential fusion reads). Options: removesup (default), separate (required for downstream work with fusions), keepsup (keeps supplementary alignments for isoform detection, does not allow gene fusion detection)')
    parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
            help='''Suppress minimap progress statements from being printed''')
    parser.add_argument('--remove_internal_priming', default=False, action='store_true',
                                       help='specify if want to remove reads with internal priming')
    parser.add_argument('-f', '--gtf', type=str,
            help='reference annotation, only used if --remove_internal_priming is specified, recommended if so')
    parser.add_argument('--intprimingthreshold', type=int, default=12,
                                            help='number of bases that are at leas 75%% As required to call read as internal priming')
    parser.add_argument('--intprimingfracAs', type=float, default=0.75,
                                            help='number of bases that are at least 75%% As required to call read as internal priming')
    parser.add_argument('--remove_singleexon', default=False, action='store_true',
                                            help='specify if want to remove unspliced reads')

    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        sys.exit(1)
    if 'align' in sys.argv:
        sys.argv.remove('align')
    args = parser.parse_args()

    # do we have multiple inputs?
    if ',' in args.reads[0]:
        args.reads = args.reads[0].split(',')
    for rfile in args.reads:
        if not os.path.exists(rfile):
            raise ValueError(f'Error: read file does not exist: {rfile}')

    if args.filtertype not in {'keepsup', 'removesup', 'separate'}:
        raise ValueError(f'filtertype {args.filtertype} is not valid. Use keepsup, removesup, or separate')
    doalignment(args)

    ##run filtering
    dofiltering(args, args.output + '_unfiltered.bam')

    pipettor.run([('rm', args.output + '_unfiltered.bam', args.output + '_unfiltered.bam.bai')])
    bedout = args.output + '.bed'
    # bam2Bed12(bamout, bedout, args.keep_supplementary)
    return bedout

if __name__ == "__main__":
    align()
