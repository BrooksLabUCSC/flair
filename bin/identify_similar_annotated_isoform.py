import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py psl gtf isos_matched.psl \n')
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
annotated_juncs = {}

for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'exon':
		continue
	if chrom not in annotated_juncs:
		annotated_juncs[chrom] = []
	this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('gene_type')-3]  # p specific to gencode v24
	this_transcript = line[8][line[8].find('transcript_id')+15:]
	this_transcript = this_transcript[:this_transcript.find('"')]

	if this_transcript != prev_transcript:
		if prev_transcript:
			annotated_juncs[chrom] += [(junctions, prev_transcript)]
		junctions = set()
		prev_transcript = this_transcript
	elif strand == '-':
		junctions.add((end, prev_start))
	else:
		junctions.add((prev_end, start))
	prev_start = start
	prev_end = end

annotated_juncs[chrom] += [junctions]
novel, total = 0, 0
seenjunctions = {}
transcript_counts = {}
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])
		if chrom not in annotated_juncs:
			continue
		junctions = get_junctions(line)
		# if chrom not in seenjunctions:
		# 	seenjunctions[chrom] = []
		# elif junctions in seenjunctions[chrom]:
		# 	continue
		# seenjunctions[chrom] += [junctions]
		total += 1
		subset = False
		transcript = ''
		name = line[9]
		for j,t in annotated_juncs[chrom]:
			if junctions == j:
				transcript = t
				break
		if transcript not in transcript_counts:
			transcript_counts[transcript] = 0
		else:
			transcript_counts[transcript] += 1

		if not transcript:
			novel += 1
			writer.writerow(line)
		else:  # annotated transcript identified
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
sys.stderr.write('{} out of {} isoforms are novel\n'.format(novel, total))

