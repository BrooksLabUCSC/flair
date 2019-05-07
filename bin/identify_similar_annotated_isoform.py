import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py psl gtf isos_matched.psl \n')
	sys.stderr.write('renames if there is an exact splice junction chain match, \
		disregards differences in TSS/TES\n')
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
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

prev_transcript = ''
# annotated_juncs = {}  # chrom: [(junctions, transcript_name), ... ]
junc_to_tn = {}  # chrom: {intron: [transcript_names], ... }
tn_to_juncs = {}  # chrom: {transcript_name: (junction1, junction2), ... }


for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'exon':
		continue
	if chrom not in junc_to_tn:
		# annotated_juncs[chrom] = []
		junc_to_tn[chrom] = {}
		tn_to_juncs[chrom] = {}
	this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('gene_type')-3]  # p specific to gencode v24
	this_transcript = line[8][line[8].find('transcript_id')+15:]
	this_transcript = this_transcript[:this_transcript.find('"')]

	if this_transcript != prev_transcript:
		if prev_transcript:
			# annotated_juncs[chrom] += [(junctions, prev_transcript)]
			tn_to_juncs[chrom][prev_transcript] = junctions
			for j in junctions:
				if j not in junc_to_tn[chrom]:
					junc_to_tn[chrom][j] = set()
				junc_to_tn[chrom][j].add(prev_transcript)


		junctions = set()
		prev_transcript = this_transcript
	elif strand == '-':
		junctions.add((end, prev_start))
	else:
		junctions.add((prev_end, start))
	prev_start = start
	prev_end = end
# annotated_juncs[chrom] += [(junctions, this_transcript)]
tn_to_juncs[chrom][this_transcript] = junctions
for j in junctions:
	junc_to_tn[chrom][j] = this_transcript

novel, total = 0, 0
seenjunctions = {}
transcript_counts = {}
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])
		if chrom not in junc_to_tn:
			continue
		junctions = get_junctions(line)

		total += 1
		subset = False
		transcript = ''

		if junctions:
			matches = set()
			for j in junctions:
				if j in junc_to_tn[chrom]:
					matches.update(junc_to_tn[chrom][j])
			for t in matches:
				if tn_to_juncs[chrom][t] == junctions:
					transcript = t  # annotated transcript identified
					break


		if not transcript:
			novel += 1
			writer.writerow(line)
		else:
			if transcript not in transcript_counts:
				transcript_counts[transcript] = 0
			else:
				transcript_counts[transcript] += 1
			if '-' in name[-3:]:
				name = name[:name.rfind('-')]
			if transcript_counts[transcript] == 0:
				if '_' in name:
					line[9] = transcript + '_' + name[name.rfind('_')+1:]
				else:
					line[9] = transcript + '_' + name
			else:
				if '_' in name:
					line[9] = transcript + '_' + name[name.rfind('_')+1:] \
					 + '-' + str(transcript_counts[transcript])			
				else:
					line[9] = transcript + '_' + name + '-' + str(transcript_counts[transcript])
			writer.writerow(line)
			line[9] = name
# sys.stderr.write('{} out of {} isoforms have novel splice junction chains\n'.format(novel, total))

