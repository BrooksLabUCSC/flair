import sys, csv
import scipy.stats as sps

try:
	psl = open(sys.argv[1])
	isbed = sys.argv[1][-3:] == 'bed'
	outfilename3 = sys.argv[2]
	outfilename5 = sys.argv[3]
	if len(sys.argv) > 4:
		counts_tsv = open(sys.argv[4])
	else:
		counts_tsv = ''
	simple = len(sys.argv) > 5
	wiggle = 10  # whatever wiggle used in sscorrect...
except:
	sys.stderr.write('usage: script.py isoforms.psl|.bed alt_acceptor.drim alt_donor.drim [counts_matrix] [simple]\n')
	sys.stderr.write('currently just outputs drimseq input\n')
	sys.exit(1)

def get_junctions_psl(starts, sizes):
	junctions = [] #set()
	for b in range(len(starts)-1):
		# junctions.add((starts[b]+sizes[b], starts[b+1]))
		junctions += [(starts[b]+sizes[b], starts[b+1], starts[b], starts[b+1]+sizes[b+1])]
	return junctions

def pslreader(psl, sample_names, iso_counts, fiveprimeon=False):
	junctiondict = {}
	for line in psl:
		line = line.rstrip().split('\t')

		if isbed:
			chrom, name, strand, start, end = line[0], line[3], line[5], int(line[1]), int(line[2])
		else:
			chrom, name, strand, start, end = line[13], line[9], line[8], int(line[15]), int(line[16])
		if '_' in name:
			name = name[:name.rfind('_')]

		if iso_counts and name not in iso_counts:
			continue

		chrom = strand + chrom  # stranded comparisons

		if isbed:
			starts = [int(n) + start for n in line[11].split(',')[:-1]]
			sizes = [int(n) for n in line[10].split(',')[:-1]]
		else:
			sizes = [int(n) for n in line[18].split(',')[:-1]]
			starts = [int(n) for n in line[20].split(',')[:-1]]
		junctions = get_junctions_psl(starts, sizes)
		if chrom not in junctiondict:
			junctiondict[chrom] = {}
		for j_index in range(len(junctions)):
			j = junctions[j_index]
			if fiveprimeon:
				fiveprime, threeprime = j[1], j[0]
				exon_end, exon_start = j[2], j[3]
			else:
				fiveprime, threeprime = j[0], j[1]
				exon_end, exon_start = j[3], j[2]
			if strand == '-':
				fiveprime, threeprime = threeprime, fiveprime
				exon_end = exon_start
			if fiveprime not in junctiondict[chrom]:
				junctiondict[chrom][fiveprime] = {}  # 5' end anchor
			if threeprime not in junctiondict[chrom][fiveprime]:
				junctiondict[chrom][fiveprime][threeprime] = {}
				junctiondict[chrom][fiveprime][threeprime]['counts'] = [0]*len(sample_names)
				junctiondict[chrom][fiveprime][threeprime]['isos'] = []# isoform list for this junction
				junctiondict[chrom][fiveprime][threeprime]['exon_end'] = exon_end  # for detecting exon skipping
			elif (not fiveprimeon and strand == '+') or (fiveprimeon and strand == '-'):
				if exon_end < junctiondict[chrom][fiveprime][threeprime]['exon_end']: # pick shorter exon end
					junctiondict[chrom][fiveprime][threeprime]['exon_end'] = exon_end
			else:
				if exon_end > junctiondict[chrom][fiveprime][threeprime]['exon_end']:
					junctiondict[chrom][fiveprime][threeprime]['exon_end'] = exon_end	
			junctiondict[chrom][fiveprime][threeprime]['isos'] += [name]
			for c in range(len(sample_names)): 
				junctiondict[chrom][fiveprime][threeprime]['counts'][c] += iso_counts[name][c]
	return junctiondict

def find_altss(alljuncs, writer, fiveprimeon=False):
	""" If fiveprimeon is True, then alternative 5' SS will be reported instead. """
	n = 0  # linenum
	for chrom in alljuncs:
		for fiveprime in alljuncs[chrom]:
			if len(alljuncs[chrom][fiveprime]) == 1:  # if there is only one 3' end, not an alt 3' junction
				continue

			if simple and len(alljuncs[chrom][fiveprime]) > 2:
				continue
			all_tp = list(alljuncs[chrom][fiveprime].keys())

			for tp1 in all_tp:  # tp1 = 3' SS
				exon_end = alljuncs[chrom][fiveprime][tp1]['exon_end']
				for tp2 in all_tp: # tp2 is also a 3' SS for the same 5' anchor as tp1
					if tp1 == tp2:
						continue
					if abs(tp2-tp1) < wiggle:
						continue
					if abs(fiveprime - tp1) > abs(fiveprime - tp2):  # this will get tested when tp2 is tp1
						continue
					inclusion = tp1
					exclusion = tp2
					strand = chrom[0]

					if (not fiveprimeon and strand == '+') or (fiveprimeon and strand == '-'):
						if tp2 > exon_end:  # exon skipping. tp2 does not overlap tp1's exon
							continue
					else:
						if tp2 < exon_end:  # exon skipping for alt SS upstream of anchor 
							continue

					chrom_clean = chrom[1:]
					event = chrom_clean+':'+str(fiveprime)+'-'+str(inclusion)+'_'+chrom_clean+':'+str(fiveprime)+'-'+str(exclusion)
					writer.writerow(['exclusion_'+str(fiveprime)+'_'+str(n), event] + alljuncs[chrom][fiveprime][exclusion]['counts'] +\
					[','.join(alljuncs[chrom][fiveprime][exclusion]['isos'])] )
					writer.writerow(['inclusion_'+str(fiveprime)+'_'+str(n), event] + alljuncs[chrom][fiveprime][inclusion]['counts'] +\
					[','.join(alljuncs[chrom][fiveprime][inclusion]['isos'])] )
					n += 1

iso_counts = {}
sample_names = []
if counts_tsv:
	sample_names = counts_tsv.readline().rstrip().split('\t')[1:]
	for line in counts_tsv:
		line = line.rstrip().split('\t')
		iso = line[0][:line[0].rfind('_')]
		iso_counts[iso] = [float(x) for x in line[1:]]

alljuncs = pslreader(psl, sample_names, iso_counts)

with open(outfilename3, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	# writer.writerow(['seqname','gene_id','event_id','alternative_transcripts','total_transcripts']) # ioe
	writer.writerow(['feature_id', 'gene_id']+sample_names+['isoform_ids'])
	find_altss(alljuncs, writer)


alljuncs = pslreader(open(sys.argv[1]), sample_names, iso_counts, fiveprimeon=True)  # reopening
with open(outfilename5, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	# writer.writerow(['seqname','gene_id','event_id','alternative_transcripts','total_transcripts'])  # ioe
	writer.writerow(['feature_id', 'gene_id']+sample_names+['isoform_ids'])

	find_altss(alljuncs, writer, fiveprimeon=True)
