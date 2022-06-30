#!/usr/bin/env python3

import sys, argparse, re, csv, math, os
from multiprocessing import Pool
from collections import Counter
from collections import namedtuple

parser = argparse.ArgumentParser(description='''for counting transcript abundances after 
	aligning reads to transcripts; for multiple mappers, only the best alignment 
	for each read is used, usage=python -s samfile -o outputfile''')
required = parser.add_argument_group('required named arguments')
required.add_argument('-s', '--sam', type=argparse.FileType('r'), required=True,
	dest='sam', help='sam file or - for STDIN')
required.add_argument('-o', '--output', default='counts.txt', type=str, required=True, action='store',
	dest='output', help='output file name')
parser.add_argument('-i', '--isoforms', type=str, action='store', dest='isoforms', default='',
	help='specify isoforms.bed or .psl file if --stringent and/or --check_splice is specified')
parser.add_argument('--stringent', default=False, action='store_true', dest='stringent',
	help='only count if read alignment passes stringent criteria')
parser.add_argument('--check_splice', default=False, action='store_true', dest='check_splice',
	help='''enforce coverage of 4 out of 6 bp around each splice site and no 
	insertions greater than 3 bp at the splice site''')
parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends',
	help='specify if reads are generated from a long read method with minimal fragmentation')
parser.add_argument('-t', '--threads', default=4, type=int, action='store', dest='t',
	help='number of threads to use')
parser.add_argument('-w', '--window', type=int, action='store', dest='w', default=10,
	help='number of bases for determining which end is best match (10)')
parser.add_argument('--quality', default=1, type=int, action='store', dest='quality',
	help='minimum quality threshold to consider if ends are to be trusted (1)')
parser.add_argument('--generate_map', default='', action='store', type=str, dest='generate_map',
	help='''specify an output path for a txt file of which isoform each read is assigned to''')
parser.add_argument('--fusion_dist', default='', action='store', dest='fusion_dist',
	help='''minimium distance between separate read alignments on the same chromosome to be
	considered a fusion, otherwise no reads will be assumed to be fusions''')
parser.add_argument('--minimal_input', default='', action='store_true', dest='minimal_input',
	help='''input file is not actually a sam file, but only the info necessary''')
args = parser.parse_args()

sam = args.sam

if args.stringent or args.fusion_dist or args.check_splice:
	if not os.path.exists(args.isoforms):
		sys.stderr.write('A valid isoforms bed file needs to be specified: {}\n'.format(args.isoforms))
		sys.exit(1)
	isbed = args.isoforms[-3:].lower() != 'psl'
	isoform_file = open(args.isoforms)

if args.fusion_dist:
	args.trust_ends = True

outfilename = args.output

min_insertion_len = 3
ss_window = 3
num_match_in_ss_window = 4
trust_ends_window = 50
large_indel_tolerance = 25

def is_stringent(tname, blocksizes, blockstarts):
	''' The additional criteria applied if --stringent is specified: Read must cover 80% of its
	corresponding transcript and the read must extend at least 25 nt into the first and last exons of the tn.'''
	isoform_length = transcript_lengths[tname]
	if sum(blocksizes) / isoform_length < 0.8:
		return False

	read_start, read_end = blockstarts[0], blockstarts[-1]+blocksizes[-1]
	first_blocksize, last_blocksize = terminal_exons[tname]['left_size'], terminal_exons[tname]['right_size']

	# covers enough bases into the first and last exons
	right_coverage = left_coverage = False
	if len(blocksizes) == 1:  # single exon transcript
		if read_start < 25 and read_end > isoform_length - 25:
			right_coverage = left_coverage = True
	else:
		if first_blocksize < 25:
		# if the terminal exon is sufficiently short, relax by 5 bp bc drna misses bases on 5' end
			if read_start < 5:
				left_coverage = True
		elif args.trust_ends:
			left_coverage = read_start <= trust_ends_window
		elif read_start <= (first_blocksize - 25):
			left_coverage = True
		else:
			return False

		if last_blocksize < 25:
			if (isoform_length - read_end) < 5:
				right_coverage = True
		elif args.trust_ends:
			right_coverage = (isoform_length - read_end) < trust_ends_window
		elif (isoform_length-last_blocksize + 25) <= read_end:
			right_coverage = True
		else:
			return False

	return right_coverage and left_coverage

def check_splice(tname, pos, covered_pos, insertion_pos):
	''' Returns True if all splice sites are covered. Splice site is considered covered if
	the number of matches according to CIGAR string  num_match_in_ss_window '''
	i = 0  # exon index
	for ss in splice_sites[tname]:
		if sum(covered_pos[ss-pos-3:ss-pos+3]) <= num_match_in_ss_window:
			return False
	return True

def are_far(transcript_1, transcript_2):
	t1_left, t1_right = terminal_exons[transcript_1]['left_pos'], terminal_exons[transcript_1]['right_pos']
	t2_left, t2_right = terminal_exons[transcript_2]['left_pos'], terminal_exons[transcript_2]['right_pos']
	return t1_left > t2_right+int(args.fusion_dist) and t1_right > t2_right or \
	t1_right < t2_left - int(args.fusion_dist) and t1_left < t2_left


def count_transcripts_for_reads(read_names):
	counts = {}  # isoform to counts
	isoform_read = {}  # isoform and their supporting read ids
	for r in read_names:
		transcripts = reads[r]  # all potential transcripts this read is assigned to

		if not args.stringent and not args.check_splice:
			if len(set(transcripts)) == 1:  # only one possible transcript
				for t in transcripts:
					quality = reads[r][t].mapq
					if t not in counts:
						counts[t] = 0
					counts[t] += 1
					if args.generate_map:
						if t not in isoform_read:
							isoform_read[t] = []
						isoform_read[t] += [r]

				continue
			# based on sorted mapq, if the read aligns much better to one transcript than the other candidates
			ordered_transcripts = sorted(transcripts.items(), key=lambda x: transcripts[x[0]].mapq)
			if (ordered_transcripts[-1][1].mapq > ordered_transcripts[-2][1].mapq \
				and ordered_transcripts[-1][1].mapq > 7) or \
				(not args.trust_ends and ordered_transcripts[-1][1].mapq >= args.quality):
				t = ordered_transcripts[-1][0]
				if t not in counts:
					counts[t] = 0
				counts[t] += 1
				if args.generate_map:
					if t not in isoform_read:
						isoform_read[t] = []
					isoform_read[t] += [r]
				continue

		if not args.trust_ends and not args.stringent and not args.check_splice:
		# multiple mapper, can't assign
			continue

		# read in cigar info for stringent/trust_ends modes into transcript_coverage dict
		# {transcript: (sum_matches, unmapped_left, unmapped_right, softclip_left, softclip_right)}
		transcript_coverage = {}
		for t, t_info in transcripts.items():
			covered_pos = []  # 0 if transcript pos isn't covered with M, 1 if pos is covered
			insertion_pos = set()
			cigar, pos = t_info.cigar, t_info.startpos

			relstart = 0
			blocksizes, relblockstarts = [], []

			matches = re.findall('([0-9]+)([A-Z])', cigar)
			num, op = int(matches[0][0]), matches[0][1]
			if op == 'H':  # check for H and then S at beginning of cigar 
				matches = matches[1:]
				# relstart += num

			num, op = int(matches[0][0]), matches[0][1]
			if op == 'S':
				softclip_left = num
				matches = matches[1:]
			else:
				softclip_left = 0

			indel_detected = False
			for m in matches:  # does not check for H and S 
				num, op = int(m[0]), m[1]
				if op == 'M':  # consumes reference
					blocksizes += [num]
					relblockstarts += [relstart]
					relstart += num
					covered_pos += [1] * num
				elif op == 'D':  # consumes reference
					relstart += num
					covered_pos += [0] * num
					if num > large_indel_tolerance:
						large_indel = True
						break
				elif op == 'N':  # consumes reference
					relstart += num
					covered_pos += [0] * num
				elif op == 'I' and num >= min_insertion_len:  # note insertions for check_splice
					if args.check_splice:
						if relstart+pos in splice_sites[t] or relstart+pos+1 in splice_sites[t] \
						or relstart+pos-1 in splice_sites[t]:
							indel_detected = True
							break
						insertion_pos.add(relstart)
					if num > large_indel_tolerance:
						indel_detected = True
						break

			if indel_detected:
				continue

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
			# a read can pass stringent assignment to a transcript but may not pass splice checking,
			# and vice-versa. if both stringent and check_splice are not specified, check each
			# individually. if only trust_ends is supplied, the read-transcript assignment does not
			# need to pass but the assignment with the highest coverage will be prioritized
			if (args.check_splice and args.stringent and check_splice(t, pos, covered_pos, insertion_pos) and \
				is_stringent(t, blocksizes, blockstarts)) \
			or (args.stringent and not args.check_splice and is_stringent(t, blocksizes, blockstarts) or \
			args.check_splice and not args.stringent and check_splice(t, pos, covered_pos, insertion_pos)) or \
			(not args.stringent and not args.check_splice and args.trust_ends):
				transcript_coverage[t] = (sum(blocksizes), read_left, \
					transcript_lengths[t] - read_right, softclip_left, softclip_right, transcripts[t].mapq)

		if not transcript_coverage:  # no transcripts passed stringent and/or check_splice criteria
			continue
		# sort by highest number of base matches
		ranked_transcripts = sorted(transcript_coverage.items(), key=lambda x: x[1][0], reverse=True)
		# then sort by highest mapq
		ranked_transcripts = sorted(ranked_transcripts, key=lambda x: x[1][5], reverse=True)

		best_t, best_t_info = ranked_transcripts[0]

		# reminder that best_t_info is sum(blocksizes), left pos, right pos, softclip_left, softclip_right)
		for t, t_info in ranked_transcripts[1:]:
			if not args.check_splice and t_info[0] + args.w < best_t_info[0]:
				continue
			# does this transcript have less softclipping than the best transcript and
			# does this transcript reach the ends better than the best transcript
			if t_info[3]+t_info[4] <= best_t_info[3]+best_t_info[4] and \
			t_info[1]+t_info[2] <= best_t_info[1]+best_t_info[2]:
				best_t, best_t_info = t, t_info
		if best_t not in counts:
			counts[best_t] = 0
		counts[best_t] += 1

		if args.fusion_dist:  # looking for the second best isoform match
			second_t, second_t_info = best_t, (0, sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize)
			for t, t_info in ranked_transcripts[1:]:
				if t_info[0] + args.w < best_t_info[0]:
					continue
				# make sure that this transcript isn't already the best transcript
				if t == best_t:
					continue

				if t_info[3]+t_info[4] <= second_t_info[3]+second_t_info[4] and \
				t_info[1]+t_info[2] <= second_t_info[1]+second_t_info[2] and \
				(terminal_exons[t]['c'] != terminal_exons[second_t]['c'] or \
					are_far(t, second_t)):
					second_t, second_t_info = t, t_info
			if second_t != best_t:  # if a second best match was found
				if second_t not in counts:
					counts[second_t] = 0
				counts[second_t] += 1


		if args.generate_map:
			if best_t not in isoform_read:
				isoform_read[best_t] = []
			isoform_read[best_t] += [r]

			if args.fusion_dist and second_t != best_t:
				if second_t not in isoform_read:
					isoform_read[second_t] = []
				isoform_read[second_t] += [r]

	if args.generate_map:
		return counts, isoform_read
	else:
		return counts, None

terminal_exons = {}  # first last exons of the firstpass flair reference transcriptome
splice_sites = {}  # splice site positions related to transcript start site
if args.stringent or args.check_splice or args.fusion_dist:
	for line in isoform_file:
		line = line.rstrip().split('\t')
		if isbed:
			name, left, right, chrom = line[3], int(line[1]), int(line[2]), line[0]
			blocksizes = [int(n) for n in line[10].split(',')[:-1]]
		else:
			name, left, right, chrom = line[9], int(line[11]), int(line[12]), line[13]
			blocksizes = [int(n) for n in line[18].split(',')[:-1]]
		terminal_exons[name] = {}
		terminal_exons[name]['left_size'] = blocksizes[0]
		terminal_exons[name]['right_size'] = blocksizes[-1]
		terminal_exons[name]['left_pos'] = left
		terminal_exons[name]['right_pos'] = right
		terminal_exons[name]['c'] = chrom
		if args.check_splice:
			splice_sites[name] = set()
			last_pos = 0
			for b in blocksizes[:-1]:
				last_pos += b
				splice_sites[name].add(last_pos)

aln = namedtuple('aln', 'cigar mapq startpos')

transcript_lengths = {}
reads = {}
for line in sam:
	if line.startswith('@'):
		if line.startswith('@SQ'):
			name = line[line.find('SN:') + 3:].split()[0]
			if name.count('|') > 2:
				name = name[:name.find('|')]
			length = float(line[line.find('LN:') + 3:].split()[0])
			transcript_lengths[name] = length
		continue
	line = line.rstrip().split('\t')
	if args.minimal_input:
		read, transcript, cigar, quality, pos = line[0], line[1], line[2], int(line[3]), int(line[4]) - 1
	else:
		read, transcript, cigar, quality, pos = line[0], line[2], line[5], int(line[4]), int(line[3]) - 1
	if quality < args.quality:
		continue
	if transcript == '*':
		continue
	if read not in reads:
		reads[read] = {}
	if transcript.count('|')>2:
		transcript = transcript[:transcript.find('|')]
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
			for t in res[1]:
				if t not in merged_map:
					merged_map[t] = []
				merged_map[t] += res[1][t]

		with open(args.generate_map, 'wt') as outfile:
			writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
			for t in merged_map:
				if t in merged_counts:
					writer.writerow([t, ','.join(merged_map[t])])



