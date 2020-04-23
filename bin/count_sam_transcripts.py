#!/usr/bin/env python3

import sys, argparse, re, csv, math, os
from multiprocessing import Pool
from collections import Counter
from collections import namedtuple

parser = argparse.ArgumentParser(description='for counting transcript abundances after aligning reads' +\
			' to transcripts; for multiple mappers, only the best alignment for each read is used', \
			usage='python -s samfile -o outputfile')
required = parser.add_argument_group('required named arguments')
required.add_argument('-s', '--sam', type=str, default='', required=True, action='store', \
	dest='sam', help='sam file')
required.add_argument('-o', '--output', default=0.25, type=str, required=True, action='store', \
	dest='output', help='output file name')
parser.add_argument('-i', '--isoforms', type=str, action='store', dest='isoforms', default='', \
	help='specify isoforms.bed or .psl file if --stringent is specified')
parser.add_argument('-w', '--window', type=int, action='store', dest='w', default=10, \
	help='number of bases for determining which end is best match (10)')
parser.add_argument('--stringent', default='', action='store_true', dest='stringent', \
	help='only count if read alignment passes stringent criteria')
parser.add_argument('-t', '--threads', default=4, type=int, action='store', dest='t', \
	help='number of threads to use')
parser.add_argument('--quality', default=1, type=int, action='store', dest='quality', \
	help='minimum quality threshold to consider if ends are to be trusted (1)')
parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends', \
	help='specify if reads are generated from a long read method with minimal fragmentation')
parser.add_argument('--generate_map', default='', action='store', type=str, dest='generate_map', \
	help='''specify an output path for a txt file of which isoform each read is assigned to''')

args = parser.parse_args()
try:
	sam = open(args.sam)
	outfilename = args.output
	if args.stringent:
		isbed = args.isoforms[-3:].lower() != 'psl'
		isoform_file = open(args.isoforms)
except:
	sys.stderr.write('Sam file issue and/or if --stringent was specified, so must a valid isoforms file\n')
	sys.exit(1)

def is_stringent(tname, blocksizes, blockstarts):
	''' The additional stringent criteria applied if --stringent is specified. Read must cover 80% of its
	corresponding transcript and the read must extend at least 25 nt into the first and last exons of the tn.'''
	read_start, read_end = blockstarts[0], blockstarts[-1]+blocksizes[-1]
	isoform_length = transcript_lengths[tname]
	first_blocksize, last_blocksize = terminal_exons[tname]['left_size'], terminal_exons[tname]['right_size']

	right_coverage = left_coverage = False
	if len(blocksizes) == 1:  # single exon transcript
		if read_start < 25 and read_end > isoform_length - 25:
			right_coverage = left_coverage = True
	else:
		if first_blocksize < 25:  # if the terminal exon is sufficiently short, relax by 2 bp
			if read_start < 2:
				left_coverage = True
		elif read_start <= (first_blocksize - 25):
			left_coverage = True

		if last_blocksize < 25:
			if (isoform_length - read_end) < 2:
				right_coverage = True
		if (isoform_length-last_blocksize + 25) <= read_end:
			right_coverage = True

	# coverage = proportion of bases of the isoform that the read covers
	coverage = sum(blocksizes) / isoform_length

	return right_coverage and left_coverage and coverage > 0.8

def count_transcripts_for_reads(read_names):
	counts = {}  # isoform to counts
	isoform_read = {}  # isoform and their supporting read ids 
	for r in read_names:
		transcripts = reads[r]  # all potential transcripts this read is assigned to

		if not args.stringent and len(set(transcripts)) == 1:  # only one possible transcript
			for t in transcripts:
				quality = reads[r][t].mapq
				if quality >= args.quality:  # if transcript passes mapq threshold
					if t not in counts:
						counts[t] = 0
					counts[t] += 1
					if args.generate_map:
						if t not in isoform_read:
							isoform_read[t] = []
						isoform_read[t] += [r]
					break
			continue

		# based on sorted mapq, if the read aligns much better to one transcript than the other candidates
		ordered_transcripts = sorted(transcripts.items(), key=lambda x: transcripts[x[0]].mapq)
		if (not args.stringent and ordered_transcripts[-1][1].mapq > ordered_transcripts[-2][1].mapq \
			and ordered_transcripts[-1][1].mapq > 7) or \
			(not args.stringent and not args.trust_ends and ordered_transcripts[-1][1].mapq >= args.quality):
			t = ordered_transcripts[-1][0]
			if t not in counts:
				counts[t] = 0
			counts[t] += 1
			if args.generate_map:
				if t not in isoform_read:
					isoform_read[t] = []
				isoform_read[t] += [r]
			continue

		if not args.trust_ends and not args.stringent:  # multiple mapper, can't assign
			continue

		# read in cigar info for stringent/trust_ends modes into transcript_coverage dict
		# {transcript: (sum_matches, unmapped_left, unmapped_right, softclip_left, softclip_right)}
		transcript_coverage = {}
		for t, t_info in ordered_transcripts:

			cigar, pos = t_info.cigar, t_info.startpos

			relstart = 0
			blocksizes, relblockstarts = [], []

			matches = re.findall('([0-9]+)([A-Z])', cigar)
			num, op = int(matches[0][0]), matches[0][1]
			if op == 'H':  # check for H and then S at beginning of cigar 
				matches = matches[1:]
				relstart += num

			num, op = int(matches[0][0]), matches[0][1]
			if op == 'S':
				softclip_left = num
				matches = matches[1:]
			else:
				softclip_left = 0

			for m in matches:  # does not check for H and S 
				num, op = int(m[0]), m[1]
				if op == 'M':  # consumes reference
					blocksizes += [num]
					relblockstarts += [relstart]
					relstart += num
				elif op == 'D':  # consumes reference
					relstart += num
				elif op == 'N':  # consumes reference
					relstart += num

			num, op = int(matches[-1][0]), matches[-1][1]
			if op == 'H':  # check for H and S at the end of cigar 
				matches = matches[:-1]
			num, op = int(matches[-1][0]), matches[-1][1]
			if op == 'S':
				softclip_right = num
			else:
				softclip_right = 0

			blockstarts = [pos + s for s in relblockstarts]
			read_left = blockstarts[0]
			read_right = blockstarts[-1] + blocksizes[-1]
			if args.stringent and is_stringent(t, blocksizes, blockstarts) \
				or args.trust_ends:
				transcript_coverage[t] = (sum(blocksizes), read_left, \
					transcript_lengths[t] - read_right, softclip_left, softclip_right)

		if args.stringent and not transcript_coverage:  # no transcripts passed stringent criteria
			continue

		# sort by highest number of base matches
		ranked_transcripts = sorted(transcript_coverage.items(), key=lambda x: x[1][0])
		best_t, best_t_info = ranked_transcripts[0]
		for t, t_info in ranked_transcripts[1:]:
			if t_info[0] + args.w < best_t_info[0]:
				continue
			# does this transcript have less softclipping than the best transcript and
			# does this transcript reach the ends better than the best transcript
			if t_info[3]+t_info[4] <= best_t_info[3]+best_t_info[4] and \
			t_info[1]+t_info[2] <= best_t_info[1]+best_t_info[2]:
				best_t, best_t_info = t, t_info
		if best_t not in counts:
			counts[best_t] = 0
		counts[best_t] += 1

		if args.generate_map:
			if best_t not in isoform_read:
				isoform_read[best_t] = []
			isoform_read[best_t] += [r]

	if args.generate_map:
		return counts, isoform_read
	else:
		return counts, None

terminal_exons = {}  # first last exons of the firstpass flair reference transcriptome
if args.stringent:
	for line in isoform_file:
		line = line.rstrip().split('\t')
		if isbed:
			name, left, right = line[3], int(line[1]), int(line[2])
			blocksizes = [int(n) for n in line[10].split(',')[:-1]]
		else:
			name, left, right = line[9], int(line[11]), int(line[12])
			blocksizes = [int(n) for n in line[18].split(',')[:-1]]
		terminal_exons[name] = {}
		terminal_exons[name]['left_size'] = blocksizes[0]
		terminal_exons[name]['right_size'] = blocksizes[-1]
		terminal_exons[name]['left_pos'] = left
		terminal_exons[name]['right_pos'] = right

aln = namedtuple('aln', 'cigar mapq startpos')
transcript_lengths = {}
reads = {}
for line in sam:
	if line.startswith('@'):
		if line.startswith('@SQ'):
			name = line[line.find('SN:') + 3:].split()[0]
			length = float(line[line.find('LN:') + 3:].split()[0])
			transcript_lengths[name] = length
		continue
	line = line.rstrip().split('\t')
	read, transcript, cigar, quality, pos = line[0], line[2], line[5], int(line[4]), int(line[3]) - 1
	if transcript == '*':
		continue
	if read not in reads:
		reads[read] = {}
	reads[read][transcript] = {}
	reads[read][transcript] = aln(cigar=cigar, mapq=quality, startpos=pos)

if __name__ == '__main__':

	grouped_reads = []
	allread_names = list(reads.keys())
	groupsize = int(math.ceil(len(allread_names)/args.t))  # cast to int bc python 2.7 needs it
	if groupsize == 0:
		grouped_reads = [allread_names]
	else:
		i = 0
		for group in range(args.t):
			new_i = i + groupsize
			grouped_reads += [allread_names[i:new_i]]
			i = new_i

	p = Pool(args.t)
	counts = p.map(count_transcripts_for_reads, grouped_reads)
	p.terminate()

	merged_counts = Counter(counts[0][0])
	for res in counts[1:]:
		merged_counts += Counter(res[0])


	with open(outfilename, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
		for t in merged_counts:
			writer.writerow([t, merged_counts[t]])

	if args.generate_map:
		merged_map = {}
		for res in counts:
			for i in res[1]:
				if i not in merged_map:
					merged_map[i] = []
				merged_map[i] += res[1][i]

		with open(args.generate_map, 'wt') as outfile:
			writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
			for i in merged_map:
				writer.writerow([i, ','.join(merged_map[i])])



