#! /usr/bin/env python3

import argparse
import os
import pipettor
import pysam
import logging
from flair.pycbio.sys import cli
from flair import FlairInputDataError

FILTER_KEEPSUP = 'keepsup'
FILTER_REMOVESUP = 'removesup'
FILTER_SEPARATE = 'separate'
FILTERS = (FILTER_KEEPSUP, FILTER_REMOVESUP, FILTER_SEPARATE)

def parse_args():
    desc = "FLAIR align outputs an unfiltered bam file and a filtered bam file for use in the downstream pipeline"
    parser = argparse.ArgumentParser(description=desc)

    reads = parser.add_argument_group('required named arguments')
    reads.add_argument('-r', '--reads', nargs='+', type=str, required=True,
                       help='FASTA/FASTQ file(s) of raw reads, either space or comma separated')
    genome = parser.add_argument_group('Either one of the following arguments is required')
    genome.add_argument('-g', '--genome', type=str,
                        help='FASTA of reference genome, can be minimap2 indexed')
    genome.add_argument('--mm_index', type=str, default='',
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
    parser.add_argument('--filtertype', type=str, choices=FILTERS, default=FILTER_REMOVESUP,
                        help='method of filtering chimeric alignments (potential fusion reads). Options: removesup (default), separate (required for downstream work with fusions), keepsup (keeps supplementary alignments for isoform detection, does not allow gene fusion detection)')  # noqa: E501
    parser.add_argument('--minfragmentsize', type=int, default=80,
                        help='minimum size of alignment kept, used in minimap -s. More important when doing downstream fusion detection')
    parser.add_argument('--maxintronlen', default='200k',
                        help='maximum intron length in genomic alignment. Longer can help recover more novel isoforms with long introns')
    parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
                        help='''Suppress minimap progress statements from being printed''')
    args = cli.parseArgsWithLogging(parser)

    reads = []
    for rfiles in args.reads:
        for rfile in rfiles.split(','):
            if not os.path.exists(rfile):
                raise FlairInputDataError(f'Error: read file does not exist: {rfile}')
            reads.append(rfile)
    args.reads = reads
    return args

def doalignment(args):
    # minimap
    mm2_cmd = ['minimap2', '-ax', 'splice', '-s', str(args.minfragmentsize),
               '-G', args.maxintronlen, '--MD', '-t', str(args.threads)]
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

    samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', args.output + '.bam', '-')
    samtools_index_cmd = ('samtools', 'index', args.output + '.bam')
    pipettor.run([mm2_cmd, samtools_sort_cmd])
    pipettor.run([samtools_index_cmd])

def dofiltering(args, inbam):
    samfile = pysam.AlignmentFile(inbam, 'rb')
    outbam = pysam.AlignmentFile(args.output + '.filtered.bam', "wb", template=samfile)
    withsup = None
    if args.filtertype == FILTER_SEPARATE:
        withsup = pysam.AlignmentFile(args.output + '_chimeric.bam', "wb", template=samfile)
    totalalignments, mappednotsec, supplementary, primary = 0, 0, 0, 0
    dropped_unmapped_secondary, dropped_quality, dropped_supplementary = 0, 0, 0
    for read in samfile.fetch():
        totalalignments += 1
        if read.is_mapped and not read.is_secondary:
            mappednotsec += 1
            if read.mapping_quality < args.quality:
                dropped_quality += 1
                logging.debug(f"read dropped: low quality ({read.mapping_quality} < {args.quality}): {read.query_name}")
            elif read.is_supplementary:
                supplementary += 1
                if args.filtertype == FILTER_SEPARATE:
                    withsup.write(read)
                elif args.filtertype == FILTER_KEEPSUP:
                    outbam.write(read)
                else:
                    dropped_supplementary += 1
                    logging.debug(f"read dropped: supplementary removed: {read.query_name}")
            else:
                primary += 1
                if read.has_tag('SA') and args.filtertype == FILTER_SEPARATE:
                    withsup.write(read)
                else:
                    outbam.write(read)
        else:
            dropped_unmapped_secondary += 1
            logging.debug(f"read dropped: unmapped or secondary: {read.query_name}")
    logging.info(f'total alignments in bam file (includes unaligned reads): {totalalignments}')
    logging.info(f'total non-secondary alignments: {mappednotsec}')
    logging.info(f'total primary alignments with quality >= {args.quality}: {primary}')
    logging.info(f'total supplementary alignments with quality >= {args.quality}: {supplementary}')
    logging.info(f'reads dropped: unmapped or secondary: {dropped_unmapped_secondary}')
    logging.info(f'reads dropped: quality < {args.quality}: {dropped_quality}')
    logging.info(f'reads dropped: supplementary removed: {dropped_supplementary}')
    samfile.close()
    outbam.close()
    pysam.index(args.output + '.filtered.bam')
    if withsup is not None:
        withsup.close()
        pysam.index(args.output + '_chimeric.bam')


def align():
    args = parse_args()
    doalignment(args)
    dofiltering(args, args.output + '.bam')

def main():
    align()


if __name__ == "__main__":
    main()
