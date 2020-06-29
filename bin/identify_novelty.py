import sys, csv, os

try:
	psl = open(sys.argv[1])
	isbed = sys.argv[1][-3:].lower() != 'psl'
	gtf = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py isoforms.psl|bed .gtf outfile_with_novel_categories.psl|bed \n')
	sys.exit(1)

def get_junctions(line):
	junctions = set()
	starts = [int(n) + 1 for n in line[20].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[18].split(',')[:-1]]  # for indexing pupropses
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_junctions_bed12(line):
	junctions = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart + 1 for n in line[11].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_exons_psl(line):
	exons = set()
	starts = [int(n) + 1 for n in line[20].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[18].split(',')[:-1]]
	for e in range(len(starts)):
		exons.add((starts[e], starts[e]+sizes[e]))
	return exons

def get_exons_bed12(line):
	exons = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	for e in range(len(starts)):
		exons.add((starts[e], starts[e]+sizes[e]))
	return exons

def overlap(coord0, coord1, tol=1):
	coord0, coord1 = sorted([coord0, coord1], key = lambda x: x[0])
	return (coord0[0] <= coord1[0] and coord1[0] <= coord0[1] - tol)

def bin_search(query, data):
	""" Query is a coordinate interval. Binary search for the query in sorted data, 
	which is a list of coordinates. Finishes when an overlapping value of query and 
	data exists and returns the index in data. """
	i = int(round(len(data)/2))  # binary search prep
	lower, upper = 0, len(data)
	while True:
		if upper - lower < 2:  # stop condition but not necessarily found
			break
		if data[i][1] < query[0]:
			lower = i
			i = int(round((i + upper)/2))
		elif data[i][0] > query[1]:
			upper = i
			i = int(round((lower + i)/2))
		else:  # found
			break
	return i

prev_transcript, prev_chrom = '', ''
all_juncs = {}
all_exons = {}
junc_to_trans = {}
trans_to_juncs = {}

for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'exon':
		continue
	if chrom not in trans_to_juncs:
		trans_to_juncs[chrom] = {}
		junc_to_trans[chrom] = {}
		all_juncs[chrom] = set()
		all_exons[chrom] = set()
	this_transcript = line[8][line[8].find('transcript_id')+15:]
	this_transcript = this_transcript[:this_transcript.find('"')]


	if this_transcript != prev_transcript:
		if prev_transcript:
			trans_to_juncs[prev_chrom][prev_transcript] = junctions
			all_juncs[chrom].update(junctions)
		junctions = set()
		prev_transcript = this_transcript
		prev_chrom = chrom
	elif strand == '-' and end < prev_start:
		junctions.add((end, prev_start))
		if (end, prev_start) not in junc_to_trans[chrom]:
			junc_to_trans[chrom][(end, prev_start)] = set()
		junc_to_trans[chrom][(end, prev_start)].add(this_transcript)
	else:
		junctions.add((prev_end, start))
		if (prev_end, start) not in junc_to_trans[chrom]:
			junc_to_trans[chrom][(prev_end, start)] = set()
		junc_to_trans[chrom][(prev_end, start)].add(this_transcript)
	all_exons[chrom].add((start, end))
	prev_start = start
	prev_end = end
trans_to_juncs[chrom][this_transcript] = junctions
all_juncs[chrom].update(junctions)

for chrom in all_exons:
	 all_exons[chrom] = sorted(list(all_exons[chrom]), key=lambda x: x[0])

novel, total = 0, 0
seenjunctions = {}
transcript_counts = {}
novelexpr = 0
annotexpr = 0
nic = 0
net = 0  # novel exon total
nloci = 0 # novel loci
linenum = 0
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for line in psl:
		# if linenum % 10000 == 0:
		# 	print(linenum, nic, net, nloci, novel)
		linenum += 1

		line = line.rstrip().split('\t')
		if isbed:
			chrom, name, start, end, num_exons = line[0], line[3], int(line[1]), int(line[2]), int(line[9])
		else:
			chrom, name, start, end, num_exons = line[13], line[9], int(line[15]), int(line[16]), int(line[17])
		if chrom not in trans_to_juncs or num_exons == 1:
			continue

		if isbed:
			junctions = get_junctions_bed12(line)
			exons = get_exons_bed12(line)
		else:
			junctions = get_junctions(line)
			exons = get_exons_psl(line)

		# if chrom not in seenjunctions:
		# 	seenjunctions[chrom] = []
		# elif junctions in seenjunctions[chrom]:
		# 	continue
		# seenjunctions[chrom] += [junctions]

		# total += 1
		# transcript = ''
		# name = line[9]
		ic = True  # in category. meaning all junctions are in gencode
		possible_transcript_matches = set()

		for j in junctions:
			if j not in all_juncs[chrom]:
				ic = False
			if j in junc_to_trans[chrom]:
				possible_transcript_matches.update(junc_to_trans[chrom][j])

		subset_iso = False  # truncated version of an annotated isoform
		transcript_match = ''
		# print(chrom, possible_transcript_matches)
		for t in possible_transcript_matches:
			annot_juncs = trans_to_juncs[chrom][t]
			if junctions == annot_juncs:
				transcript_match = t
				break
			elif str(sorted(list(junctions)))[1:-1] in str(sorted(list(annot_juncs)))[1:-1]:
				subset_iso = True

		novelexon = 0
		if not subset_iso and not ic:
			for e in exons:  # this isoform
				if e not in all_exons[chrom]:
					overlaps = False
					# for e_ in all_exons[chrom]:
					# 	if overlap(e, e_,0):
					# 		overlaps = True
					# 		break
					i = bin_search(e, all_exons[chrom])
					for e_ in all_exons[chrom][i-2:i+2]:
						if overlap(e, e_, 0):
							overlaps = True
							break
					if not overlaps:
						novelexon += 1

		# if transcript not in transcript_counts:
		# 	transcript_counts[transcript] = 0
		# else:
		# 	transcript_counts[transcript] += 1

		if not transcript_match:
			# novel += 1
			# try:
			# 	novelexpr += int(line[-1])
			# except:
			# 	novelexpr += 1
			if subset_iso:
				writer.writerow(line + ['novel_subset'])
			elif ic:
				writer.writerow(line + ['novel_in_category'])
				# nic += 1
			elif novelexon:
				if novelexon < num_exons:
					writer.writerow(line + ['novel_exon'])
					# net += 1
				else:
					# nloci += 1
					writer.writerow(line + ['novel_locus'])
			else:
				writer.writerow(line + ['novel_junction'])

		else:  # annotated transcript identified
			# if '-' in name[-3:]:
			# 	name = name[:name.rfind('-')]
			# if transcript_counts[transcript] == 0:
			# 	if '_' in name:
			# 		line[9] = transcript + '_' + name[name.rfind('_')+1:]
			# 	else:
			# 		line[9] = transcript + '_' + name
			# else:
			# 	if '_' in name:
			# 		line[9] = transcript + '_' + name[name.rfind('_')+1:] \
			# 		 + '-' + str(transcript_counts[transcript])			
			# 	else:
			# 		line[9] = transcript + '_' + name + '-' + str(transcript_counts[transcript])
			writer.writerow(line + ['annotated'])
			# try:
			# 	annotexpr += int(line[-1])
			# except:
			# 	annotexpr += 1
# sys.stderr.write('{} out of {} isoforms are novel\n'.format(novel, total))
# sys.stderr.write('{}/{} reads are assigned to novel ones\n'.format(novelexpr, novelexpr+annotexpr))
# sys.stderr.write('{} novel in category\n'.format(nic))
# sys.stderr.write('{} novel loci\n'.format(nloci))
# sys.stderr.write('{} novel with at least one completely novel exon\n'.format(net))




