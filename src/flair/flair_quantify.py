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
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1,
		help='''minimum MAPQ of read assignment to an isoform (1)''')
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
	args, unknown = parser.parse_known_args()
	if unknown:
		sys.stderr.write('Quantify unrecognized arguments: {}\n'.format(' '.join(unknown)))

	if isoform_sequences:
		args.i = isoform_sequences
		args.o += '.counts_matrix.tsv'
	if (args.stringent or args.check_splice):
		if not args.isoforms:
			sys.stderr.write('Please specify isoform models as .bed file using --isoform_bed\n')
			return 1
		elif not os.path.exists(args.isoforms):
			sys.stderr.write('Isoform models bed file path does not exist\n')
			return 1
		elif args.isoforms.endswith('.psl'):
			sys.stderr.write('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.\n')
			return 1
	if not os.path.exists(args.i):
		sys.stderr.write('Isoform sequences fasta file path does not exist\n')
		return 1

	samData = list()
	with codecs.open(args.r, 'r', encoding='utf-8', errors='ignore') as lines:
		for line in lines:
			cols = line.rstrip().split('\t')
			if len(cols) < 4:
				sys.stderr.write(f'Expected 4 columns in tab-delimited manifest.tsv, got {len(cols)}. Exiting.\n')
				return 1

			sample, group, batch, readFile = cols
			if args.sample_id_only is False:
				if '_' in sample or '_' in group or '_' in batch:
					sys.stderr.write(f'Please do not use underscores in the id, condition, or batch fields of {args.r}. Exiting. \n')
					return 1

			readFileRoot = tempfile.NamedTemporaryFile().name
			if args.temp_dir != '':
				if not os.path.isdir(args.temp_dir):
					if subprocess.call(['mkdir', '-p', args.temp_dir]):
						sys.stderr.write('Could not make temporary directory {}\n'.format(args.temp_dir))
						return 1
				readFileRoot = args.temp_dir + '/' + readFileRoot[readFileRoot.rfind('/')+1:]

			if not os.path.exists(cols[3]):
				sys.stderr.write('Query file path does not exist: {}\n'.format(cols[3]))
				return 1
			samData.append(cols + [readFileRoot + '.sam'])
	sys.stderr.write('Writing temporary files with prefixes similar to {}\n'.format(readFileRoot))

	for num, sample in enumerate(samData, 0):
		sys.stderr.write('Step 1/3. Aligning sample %s_%s, %s/%s \n' % (sample[0], sample[2], num+1, len(samData)))
		mm2_command = ['minimap2', '--MD', '-a', '-N', '4', '-t', str(args.t), args.i, sample[-2]]

		# TODO: Replace this with proper try/except Exception as ex
		try:
			# samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', sample[-1], '-')
			# samtools_index_cmd = ('samtools', 'index', sample[-1])
			# pipettor.run([mm2_cmd, samtools_sort_cmd])
			# pipettor.run([samtools_index_cmd])
			if subprocess.call(mm2_command, stdout=open(sample[-1], 'w'),
				stderr=open(sample[-1]+'.mm2_stderr.txt', 'w')):
				sys.stderr.write('Check {} file\n'.format(sample[-1]+'.mm2_stderr.txt'))
				return 1
		except:
			sys.stderr.write('''Possible minimap2 error, please check that all file, directory,
				and executable paths exist\n''')
			return 1
		subprocess.check_call(['rm', sample[-1]+'.mm2_stderr.txt'])
		sys.stderr.flush()

	countData = dict()
	for num, data in enumerate(samData):
		sample, group, batch, readFile, samOut = data
		sys.stderr.write('Step 2/3. Quantifying isoforms for sample %s_%s: %s/%s \n' % (sample, batch, num+1, len(samData)))

		count_cmd = ['count_sam_transcripts.py', '-s', samOut,
			'-o', samOut+'.counts.txt', '-t', str(args.t), '--quality', str(args.quality)]
		if args.trust_ends:
			count_cmd += ['--trust_ends']
		if args.stringent:
			count_cmd += ['--stringent']
		if args.check_splice:
			count_cmd += ['--check_splice']
		if args.check_splice or args.stringent:
			count_cmd += ['-i', args.isoforms]
		if args.generate_map or args.output_bam:
			count_cmd += ['--generate_map', args.o+'.'+sample+'.'+group+'.isoform.read.map.txt']

		subprocess.check_call(count_cmd)
		for line in open(samOut+'.counts.txt'):
			line = line.rstrip().split('\t')
			iso, numreads = line[0], line[1]
			if iso not in countData:
				countData[iso] = np.zeros(len(samData))
			countData[iso][num] = numreads
		sys.stderr.flush()

		if args.output_bam:
			sys.stderr.write('Step 2.5/3. Filtering bam reads for sample %s_%s: %s/%s \r' % (sample, batch, num+1, len(samData)))
			readToIso = {}
			for line in open(args.o+'.'+sample+'.'+group+'.isoform.read.map.txt'):
				line = line.rstrip().split('\t', 1)
				for r in line[1].split(','):
					readToIso[r] = line[0]

			newsam = open(samOut.split('.sam')[0] + '-filtered.sam', 'w')
			nametoseq = {}
			impsecondary = []
			for line in open(samOut):
				if line[0] == '@': newsam.write(line)
				else:
					line = line.split('\t')
					read, iso, flag = line[0], line[2], line[1]
					if flag == '0' or flag == '16':
						nametoseq[read] = line[9]
						if read in readToIso and readToIso[read] == iso: newsam.write('\t'.join(line))
					elif read in readToIso and readToIso[read] == iso:
						impsecondary.append(line)
			for line in impsecondary:
				line[9] = nametoseq[line[0]]
				line[5] = line[5].replace('H', 'S')
				newsam.write('\t'.join(line))
			newsam.close()
			if subprocess.call(['samtools', 'sort', '-@', str(args.t), samOut.split('.sam')[0] + '-filtered.sam', '-o', args.o+'.'+sample+'.'+group+'.flair.aligned.bam']):
				sys.stderr.write(f'Samtools issue with sorting minimap2 sam\n')
				return 1
			# subprocess.check_call(['rm', samOut.split('.sam')[0] + '-filtered.sam'])

			subprocess.check_call(['samtools', 'index', args.o+'.'+sample+'.'+group+'.flair.aligned.bam'])

		subprocess.check_call(['rm', samOut])

	sys.stderr.write('Step 3/3. Writing counts to {} \n'.format(args.o+'.counts.tsv'))
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
	status = quantify()
	os.exit(status)
