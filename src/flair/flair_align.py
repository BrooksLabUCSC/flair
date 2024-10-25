#! /usr/bin/env python3

import sys
import argparse
import os
import pipettor
import pysam
os.environ['OPENBLAS_NUM_THREADS'] = '1'
# from bam2Bed12 import bam2Bed12

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
	#chr1	476363	497259	ENST00000455464.7_ENSG00000237094.12	1000	-	476363	497259	0	3	582,169,151,	0,8676,20745,
	esizes, estarts = intronChainToestarts(intronblocks,alignstart, refpos)
	rgbcolor = unknownTxn
	if juncDirection == "+": rgbcolor = positiveTxn
	elif juncDirection == "-": rgbcolor = negativeTxn
	else:
		juncDirection = "-" if is_reverse else "+"
	outline = [referencename, str(alignstart), str(refpos), readname, str(qualscore), juncDirection, str(alignstart), str(refpos), rgbcolor, str(len(intronblocks) + 1), ','.join([str(x) for x in esizes]) + ',', ','.join([str(x) for x in estarts]) + ',']
	return outline

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
	parser.add_argument('--quality', type=int, default=1,
		help='minimum MAPQ of read alignment to the genome (1)')
	parser.add_argument('--minfragmentsize', type=int, default=80,
						help='minimum size of alignment kept, used in minimap -s. More important when doing downstream fusion detection')
	parser.add_argument('-f', '--filtertype', type=str, default='keepsup',
        help='method of filtering chimeric alignments (potential fusion reads). Options: removesup (default), separate (required for downstream work with fusions), keepsup (keeps supplementary alignments for isoform detection, does not allow gene fusion detection)')
	parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
		help='''Suppress minimap progress statements from being printed''')
	
	no_arguments_passed = len(sys.argv) == 1
	if no_arguments_passed:
		parser.print_help()
		sys.exit(1)
	if 'align' in sys.argv:
		sys.argv.remove('align')
	args, unknown = parser.parse_known_args()
	if unknown and not args.quiet:
		sys.stderr.write('Align unrecognized arguments: {}\n'.format(' '.join(unknown)))

	# do we have multiple inputs?
	if ',' in args.reads[0]:
		args.reads = args.reads[0].split(',')
	for rfile in args.reads:
		if not os.path.exists(rfile):
			sys.stderr.write(f'Check that read file {rfile} exists\n')
			sys.exit(1)

	# define outputs
	bamout = args.output+'.bam'
	bedout = args.output+'.bed'

	# minimap
	mm2_cmd = ['minimap2', '-ax', 'splice', '-s', str(args.minfragmentsize), '-t', str(args.threads), args.genome]+args.reads
	if args.mm_index:
		mm2_cmd[5] = args.mm_index
	if args.nvrna:
		mm2_cmd[3:3] = ['-uf', '-k14']
	if args.junction_bed:
		mm2_cmd[3:3] = ['--junc-bed', args.junction_bed]
	# if str(args.keep_supplementary) != '0':
	# 	mm2_cmd[3:3] = ['-N', str(args.keep_supplementary)]
	# else:
	mm2_cmd[3:3] = ['--secondary=no']
	mm2_cmd = tuple(mm2_cmd)
	
	# samtools; the dash at the end means STDIN
	# samtools_filter_cmd = ('samtools', 'view', '-q', str(args.quality), '-h', '-')
	samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', args.output+'_unfiltered.bam', '-')
	samtools_index_cmd = ('samtools', 'index', args.output+'_unfiltered.bam')
	if not args.quiet:
		print('flair align minimap cmd:', mm2_cmd, file=sys.stderr)
		# print('flair align samtools filter cmd:', samtools_filter_cmd, file=sys.stderr)
		print('flair align samtools sort cmd:', samtools_sort_cmd, file=sys.stderr)
		print('flair align samtools index cmd:', samtools_index_cmd, file=sys.stderr)
	# pipettor.run([mm2_cmd, samtools_filter_cmd, samtools_sort_cmd])
	pipettor.run([mm2_cmd, samtools_sort_cmd])
	pipettor.run([samtools_index_cmd])

	##run filtering
	samfile = pysam.AlignmentFile(args.output + '_unfiltered.bam', "rb")
	if args.filtertype == 'separate':
		withsup = pysam.AlignmentFile(args.output + '_chimeric.bam', "wb", template=samfile)
	trash = pysam.AlignmentFile(args.output + '_trash.bam', "wb", template=samfile)
	fullyaligned = pysam.AlignmentFile(args.output + '.bam', "wb", template=samfile)
	outbed = open(bedout, 'w')
	for read in samfile.fetch():
		if read.is_mapped and not read.is_secondary:
			if read.has_tag('SA'):
				if args.filtertype == 'separate':
					withsup.write(read)
				elif args.filtertype == 'removesup':
					trash.write(read)
				elif read.is_supplementary:
					trash.write(read)
				else:
					mapq = read.mapping_quality
					if mapq > args.quality:
						fullyaligned.write(read)
						juncstrand = inferMM2JuncStrand(read)
						bedline = bed_from_cigar(read.reference_start, read.is_reverse, read.cigartuples,
												 read.query_name, read.reference_name, mapq, juncstrand)
						outbed.write('\t'.join(bedline) + '\n')
					else:
						trash.write(read)
			else:
				mapq = read.mapping_quality
				if mapq > args.quality:
					fullyaligned.write(read)
					juncstrand = inferMM2JuncStrand(read)
					bedline = bed_from_cigar(read.reference_start, read.is_reverse, read.cigartuples, read.query_name,
											 read.reference_name, mapq, juncstrand)
					outbed.write('\t'.join(bedline) + '\n')
				else:
					trash.write(read)
		else:
			trash.write(read)
	samfile.close()
	if args.filtertype == 'separate':
		withsup.close()
		pysam.index(args.output + '_chimeric.bam')
	trash.close()
	fullyaligned.close()
	outbed.close()
	pysam.index(args.output + '_trash.bam')
	pysam.index(args.output + '.bam')

	pipettor.run([('rm', args.output + '_unfiltered.bam', args.output + '_unfiltered.bam.bai')])

	# bam2Bed12(bamout, bedout, args.keep_supplementary)
	return bedout

if __name__ == "__main__":
	align()

