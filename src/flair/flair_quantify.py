#! /usr/bin/env python3

import os
import sys
import re
import argparse
import subprocess
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import codecs
import tempfile
import time
import pipettor, pysam


def quantify(isoform_sequences=''):
    parser = argparse.ArgumentParser(description='flair-quantify parse options',
            usage='flair quantify -r reads_manifest.tsv -i isoforms.fa [options]')
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
                    subprocess.check_call(['mkdir', '-p', args.temp_dir])
                readFileRoot = args.temp_dir + '/' + readFileRoot[readFileRoot.rfind('/')+1:]

            if not os.path.exists(readFile):
                raise Exception('Query file path does not exist: {}'.format(readFile))

            samData.append(cols + [readFileRoot + '.sam'])
    sys.stderr.write('Writing temporary files with prefixes similar to {}\n'.format(readFileRoot))

    for num, data in enumerate(samData, 0):
        sample, group, batch, readFile, samOut = data
        sys.stderr.write('Aligning sample %s_%s, %s/%s \n' % (sample, batch, num+1, len(samData)))
        mm2_command = ['minimap2', '--MD', '-a', '-N', '4', '-t', str(args.t), args.i, readFile]

        # TODO: Replace this with proper try/except Exception as ex
        try:
            if subprocess.call(mm2_command, stdout=open(samOut, 'w'),
                               stderr=open(samOut+'.mm2_stderr.txt', 'w')):
                raise Exception('Check {} file'.format(samOut+'.mm2_stderr.txt'))
        except:
            raise Exception('''Possible minimap2 error, please check that all file, directory and executable paths exist''')
        subprocess.check_call(['rm', samOut+'.mm2_stderr.txt'])
        sys.stderr.flush()



        sys.stderr.write('Quantifying isoforms for sample %s_%s: %s/%s \n' % (sample, batch, num+1, len(samData)))

        count_cmd = ['filter_transcriptome_align.py', '-s', samOut,
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


        subprocess.check_call(count_cmd)
        sys.stderr.flush()

        # if args.output_bam:
        #     sys.stderr.write('Filtering bam reads for sample %s_%s: %s/%s \r' % (sample, batch, num+1, len(samData)))
        #     readToIso = {}
        #     for line in open(args.o+'.'+sample+'.'+group+'.isoform.read.map.txt'):
        #         line = line.rstrip().split('\t', 1)
        #         for r in line[1].split(','):
        #             readToIso[r] = line[0]
        #
        #     newsam = open(samOut.split('.sam')[0] + '-filtered.sam', 'w')
        #     nametoseq = {}
        #     impsecondary = []
        #     for line in open(samOut):
        #         if line[0] == '@': newsam.write(line)
        #         else:
        #             line = line.split('\t')
        #             read, iso, flag = line[0], line[2], line[1]
        #             if flag == '0' or flag == '16':
        #                 nametoseq[read] = [line[9], line[10]]
        #                 if read in readToIso and readToIso[read] == iso:
        #                     # if line[4] == '0': line[4] = '60'
        #                     newsam.write('\t'.join(line))
        #             elif read in readToIso and readToIso[read] == iso:
        #                 impsecondary.append(line)
        #     for line in impsecondary:
        #         line[9] = nametoseq[line[0]][0]
        #         line[10] = nametoseq[line[0]][1]
        #         line[5] = line[5].replace('H', 'S')
        #         # if line[4] == '0': line[4] = '60'
        #         newsam.write('\t'.join(line))
        #     newsam.close()
        #     subprocess.check_call(['samtools', 'sort', '-@', str(args.t), samOut.split('.sam')[0] + '-filtered.sam', '-o', args.o+'.'+sample+'.'+group+'.flair.aligned.bam'])
        #     subprocess.check_call(['samtools', 'index', args.o+'.'+sample+'.'+group+'.flair.aligned.bam'])
        #
        # subprocess.check_call(['rm', samOut])

    sys.stderr.write('Writing counts to {} \n'.format(args.o + '.counts.tsv'))
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
        countMatrix.write('%s\t%s\n' % (f, '\t'.join(str(x) for x in countData[f])))

    countMatrix.close()
    sys.stderr.flush()
    sys.stderr.write('\n')

    if args.tpm:
        subprocess.check_call(['counts_to_tpm.py', args.o+'.counts.tsv', args.o+'.tpm.tsv'])
    return args.o+'.counts.tsv'

if __name__ == '__main__':
    # FIXME: need proper error handling
    sys.exit(quantify())
