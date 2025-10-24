#! /usr/bin/env python3

import os
import sys
import re
import argparse
import pipettor
import os
import numpy as np
import codecs
import tempfile
import time
import logging
import pipettor, pysam

os.environ['OPENBLAS_NUM_THREADS'] = '1'

def quantify(isoform_sequences=''):
    parser = argparse.ArgumentParser(description='takes in many long-read RNA-seq reads files and quantifies them against a single transcriptome. A stringent, full-read-match-based approach')
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-r', '--reads_manifest', action='store', dest='r', type=str,
                    required=True, help='Tab delimited file containing sample id, condition, batch, reads.fq')
    if not isoform_sequences:
        required.add_argument('-i', '--isoforms', action='store', dest='i',
                type=str, required=True, help='FastA of FLAIR collapsed isoforms')
    parser.add_argument('-o', '--output', type=str, action='store', dest='o', default='flair.quantify',
            help='''output file name base for FLAIR quantify (default: flair.quantify)''')
    parser.add_argument('-t', '--threads', type=int,
            action='store', dest='t', default=4, help='minimap2 number of threads (4)')
    parser.add_argument('--temp_dir', default='', action='store', dest='temp_dir',
            help='''directory to put temporary files. use './" to indicate current directory
            (default: python tempfile directory)''')
    parser.add_argument('--sample_id_only', default=False, action='store_true', dest='sample_id_only',
            help='''only use sample id in output header''')
    parser.add_argument('--tpm', action='store_true', dest='tpm', default=False,
            help='Convert counts matrix to transcripts per million and output as a separate file named <output>.tpm.tsv')
    parser.add_argument('--quality', type=int, action='store', dest='quality', default=0,
            help='''minimum MAPQ of read assignment to an isoform (0)''')
    parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends',
            help='specify if reads are generated from a long read method with minimal fragmentation')
    parser.add_argument('--generate_map', default=False, action='store_true', dest='generate_map',
            help='''create read-to-isoform assignment files for each sample (default: not specified)''')
    parser.add_argument('--isoform_bed', '--isoformbed', default='', type=str, action='store', dest='isoforms',
            help='''isoform .bed file, must be specified if --stringent or check_splice is specified''')
    parser.add_argument('--stringent', default=False, action='store_true', dest='stringent',
            help='''Supporting reads must cover 80 percent of their isoform and extend at least 25 nt into the
            first and last exons. If those exons are themselves shorter than 25 nt, the requirement becomes
            'must start within 4 nt from the start" or "must end within 4 nt from the end" ''')
    parser.add_argument('--check_splice', default=False, action='store_true', dest='check_splice',
            help='''enforce coverage of 4 out of 6 bp around each splice site and no
            insertions greater than 3 bp at the splice site''')
    parser.add_argument('--output_bam', default=False, action='store_true', dest='output_bam',
                                            help='whether to output bam file of reads aligned to correct isoforms')
    args = parser.parse_args()
    if isoform_sequences:
        args.i = isoform_sequences
        args.o += '.counts_matrix.tsv'
    if (args.stringent or args.check_splice):
        if not args.isoforms:
            raise Exception('Please specify isoform models as .bed file using --isoform_bed')
        elif not os.path.exists(args.isoforms):
            raise Exception('Isoform models bed file path does not exist: ' + args.isoforms)
        elif args.isoforms.endswith('.psl'):
            raise Exception('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.')
    if not os.path.exists(args.i):
        raise Exception('Isoform sequences fasta file path does not exist: ' + args.i)

    samData = list()
    with codecs.open(args.r, 'r', encoding='utf-8', errors='ignore') as lines:
        for line in lines:
            cols = line.rstrip().split('\t')
            if len(cols) < 4:
                raise Exception(f'Expected 4 columns in tab-delimited manifest.tsv, got {len(cols)}. Exiting.')

            sample, group, batch, readFile = cols
            if args.sample_id_only is False:
                if '_' in sample or '_' in group or '_' in batch:
                    raise Exception(f'Please do not use underscores in the id, condition, or batch fields of {args.r}.')

            readFileRoot = tempfile.NamedTemporaryFile().name
            if args.temp_dir != '':
                if not os.path.isdir(args.temp_dir):
                    pipettor.run([('mkdir', '-p', args.temp_dir)])
                readFileRoot = args.temp_dir + '/' + readFileRoot[readFileRoot.rfind('/')+1:]

            if not os.path.exists(readFile):
                raise Exception('Query file path does not exist: {}'.format(readFile))

            samData.append(cols + [readFileRoot + '.sam'])
    logging.info(f'Writing temporary files with prefixes similar to {readFileRoot}')

    for num, data in enumerate(samData, 0):
        sample, group, batch, readFile, samOut = data
        logging.info(f'Aligning and quantifying isoforms for sample {sample}_{batch}, {num+1}/{len(samData)}')
        mm2_command = ('minimap2', '--MD', '-a', '-N', '4', '-t', str(args.t), args.i, readFile)

        count_cmd = ['filter_transcriptome_align.py', '-s', '-',
                     '-o', samOut+'.counts.txt', '-t', str(args.t), '--quality', str(args.quality)]
        if args.trust_ends:
            count_cmd += ['--trust_ends']
        if args.stringent:
            count_cmd += ['--stringent']
        if args.check_splice:
            count_cmd += ['--check_splice']
        if args.check_splice or args.stringent:
            count_cmd += ['-i', args.isoforms]
        if args.generate_map:
            count_cmd += ['--generate_map', args.o+'.'+sample+'.'+group+'.isoform.read.map.txt']
        if args.output_bam:
            count_cmd += ['--output_bam', args.o+'.'+sample+'.'+group+'.flair.aligned.bam']


        pipettor.run([mm2_command, tuple(count_cmd)])

    logging.info(f'Writing counts to {args.o}.counts.tsv')
    countData = dict()
    for num, data in enumerate(samData):
        sample, group, batch, readFile, samOut = data
        for line in open(samOut+'.counts.txt'):
            line = line.rstrip().split('\t')
            iso, numreads = line[0], line[1]
            if iso not in countData:
                countData[iso] = np.zeros(len(samData))
            countData[iso][num] = numreads



    countMatrix = open(args.o+'.counts.tsv', 'w')

    if args.sample_id_only:
        countMatrix.write('\t'.join(['ID']+[x[0] for x in samData])+'\n')
    else:
        countMatrix.write('ids\t%s\n' % '\t'.join(['_'.join(x[:3]) for x in samData]))

    features = sorted(list(countData.keys()))
    for f in features:
        countMatrix.write('%s\t%s\n' % (f, '\t'.join(str(int(x)) for x in countData[f])))

    countMatrix.close()

    if args.tpm:
        pipettor.run([('counts_to_tpm.py', args.o+'.counts.tsv', args.o+'.tpm.tsv')])
    return args.o+'.counts.tsv'

if __name__ == '__main__':
    # FIXME: need proper error handling
    sys.exit(quantify())
