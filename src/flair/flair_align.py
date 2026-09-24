#! /usr/bin/env python3

import os
import pipettor
import pysam
import logging
from flair import FlairInputDataError

FILTER_KEEPSUP = 'keepsup'
FILTER_REMOVESUP = 'removesup'
FILTER_SEPARATE = 'separate'
FILTERS = (FILTER_KEEPSUP, FILTER_REMOVESUP, FILTER_SEPARATE)

###
# command line
###
def add_subparser(subparsers):
    desc = "FLAIR align outputs an unfiltered bam file and a filtered bam file for use in the downstream pipeline"
    parser = subparsers.add_parser('align', help="Align reads to the genome",
                                   description=desc)
    reads = parser.add_argument_group('required named arguments')
    reads.add_argument('-r', '--reads', nargs='+', type=str, required=True,
                       help='FASTA/FASTQ file(s) of raw reads, either space or comma separated')
    genome = parser.add_mutually_exclusive_group(required=True)
    genome.add_argument('-g', '--genome', type=str,
                        help='FASTA of reference genome, can be minimap2 indexed')
    genome.add_argument('--mm_index', type=str, default='',
                        help='minimap2 index .mmi file')
    parser.add_argument('-o', '--output', default='flair.aligned',
                        help='output file name base (default: %(default)s)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (default: %(default)s)')
    parser.add_argument('--junction_bed', default='',
                        help='annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment')
    parser.add_argument('--native_rna', action='store_true',
                        help='use native-RNA specific alignment parameters for minimap2')
    parser.add_argument('--quality', type=int, default=0,
                        help='minimum MAPQ of read alignment to the genome (default: %(default)s)')
    parser.add_argument('--filter_type', type=str, choices=FILTERS, default=FILTER_REMOVESUP,
                        help='method of filtering chimeric alignments (potential fusion reads); '
                        'separate is required for downstream work with fusions, keepsup keeps '
                        'supplementary alignments for isoform detection and does not allow gene '
                        'fusion detection (default: %(default)s)')
    parser.add_argument('--min_fragment_size', type=int, default=80,
                        help='minimum size of alignment kept, used in minimap -s. More important '
                        'when doing downstream fusion detection (default: %(default)s)')
    parser.add_argument('--max_intron_len', default='200k',
                        help='maximum intron length in genomic alignment. Longer can help recover '
                        'more novel isoforms with long introns (default: %(default)s)')
    parser.set_defaults(entry=align_cmd)

def align_cmd(args):
    align(reads=split_reads_files(args.reads), genome=args.genome, mm_index=args.mm_index,
          output=args.output, threads=args.threads, junction_bed=args.junction_bed,
          native_rna=args.native_rna, quality=args.quality, filter_type=args.filter_type,
          min_fragment_size=args.min_fragment_size, max_intron_len=args.max_intron_len)

def split_reads_files(reads_arg):
    "reads files may be comma-separated as well as repeated; check them up front"
    reads = []
    for rfiles in reads_arg:
        for rfile in rfiles.split(','):
            if not os.path.exists(rfile):
                raise FlairInputDataError(f'read file does not exist: {rfile}')
            reads.append(rfile)
    return reads

###
# alignment
###
def doalignment(*, reads, genome, mm_index, output, threads, junction_bed, native_rna,
                min_fragment_size, max_intron_len):
    # minimap
    mm2_cmd = ['minimap2', '-ax', 'splice', '-s', str(min_fragment_size),
               '-G', max_intron_len, '--MD', '-t', str(threads)]
    if native_rna:
        mm2_cmd += ['-uf', '-k14']
    if junction_bed:
        mm2_cmd += ['--junc-bed', junction_bed]
    mm2_cmd += ['--secondary=no']
    if mm_index:
        mm2_cmd += [mm_index]
    else:
        mm2_cmd += [genome]
    mm2_cmd += reads

    samtools_sort_cmd = ('samtools', 'sort', '-@', str(threads), '-o', output + '.bam', '-')
    samtools_index_cmd = ('samtools', 'index', output + '.bam')
    pipettor.run([mm2_cmd, samtools_sort_cmd])
    pipettor.run([samtools_index_cmd])

def dofiltering(inbam, *, output, quality, filter_type):
    samfile = pysam.AlignmentFile(inbam, 'rb')
    outbam = pysam.AlignmentFile(output + '.filtered.bam', "wb", template=samfile)
    withsup = None
    if filter_type == FILTER_SEPARATE:
        withsup = pysam.AlignmentFile(output + '_chimeric.bam', "wb", template=samfile)
    totalalignments, mappednotsec, supplementary, primary = 0, 0, 0, 0
    dropped_unmapped_secondary, dropped_quality, dropped_supplementary = 0, 0, 0
    for read in samfile.fetch():
        totalalignments += 1
        if read.is_mapped and not read.is_secondary:
            mappednotsec += 1
            if read.mapping_quality < quality:
                dropped_quality += 1
                logging.debug(f"read dropped: low quality ({read.mapping_quality} < {quality}): {read.query_name}")
            elif read.is_supplementary:
                supplementary += 1
                if filter_type == FILTER_SEPARATE:
                    withsup.write(read)
                elif filter_type == FILTER_KEEPSUP:
                    outbam.write(read)
                else:
                    dropped_supplementary += 1
                    logging.debug(f"read dropped: supplementary removed: {read.query_name}")
            else:
                primary += 1
                if read.has_tag('SA') and filter_type == FILTER_SEPARATE:
                    withsup.write(read)
                else:
                    outbam.write(read)
        else:
            dropped_unmapped_secondary += 1
            logging.debug(f"read dropped: unmapped or secondary: {read.query_name}")
    # fetch() with no region walks the index, which never yields an unplaced unmapped
    # read, so this count cannot include them
    logging.info(f'total alignments in bam file: {totalalignments}')
    logging.info(f'total non-secondary alignments: {mappednotsec}')
    logging.info(f'total primary alignments with quality >= {quality}: {primary}')
    logging.info(f'total supplementary alignments with quality >= {quality}: {supplementary}')
    logging.info(f'reads dropped: unmapped or secondary: {dropped_unmapped_secondary}')
    logging.info(f'reads dropped: quality < {quality}: {dropped_quality}')
    logging.info(f'reads dropped: supplementary removed: {dropped_supplementary}')
    samfile.close()
    outbam.close()
    pysam.index(output + '.filtered.bam')
    if withsup is not None:
        withsup.close()
        pysam.index(output + '_chimeric.bam')


def align(*, reads, genome, mm_index, output, threads, junction_bed, native_rna,
          quality, filter_type, min_fragment_size, max_intron_len):
    doalignment(reads=reads, genome=genome, mm_index=mm_index, output=output,
                threads=threads, junction_bed=junction_bed, native_rna=native_rna,
                min_fragment_size=min_fragment_size, max_intron_len=max_intron_len)
    dofiltering(output + '.bam', output=output, quality=quality, filter_type=filter_type)
