#!/usr/bin/env python3
import sys, csv, os

try:
	psl = open(sys.argv[1])
	isgtf = sys.argv[1][-3:] == 'gtf'
	isbed = sys.argv[1][-3:] == 'bed'
	outfilenamebase = sys.argv[2]
	if len(sys.argv) > 3:
		counts_tsv = open(sys.argv[3])
	else:
		counts_tsv = ''
	wiggle = 10  # minimum distance apart for alt SS to be tested
except:
	sys.stderr.write('usage: call_diffsplice_events.py .psl|.bed out.tsv [counts_tsv]\n')
	sys.exit(1)

def overlap(coords0, coords1):
	return coords1[0] >= coords0[0] and coords1[0] <= coords0[1] or \
		coords1[1] >= coords0[0] and coords1[1] <= coords0[1]

def get_junctions_psl(starts, sizes):
	junctions = []
	for b in range(len(starts)-1):
		junctions += [(starts[b]+sizes[b], starts[b+1], starts[b], starts[b+1]+sizes[b+1])]
	return junctions

def parse_iso_id(iso_gene):
	if '_' not in iso_gene:
		return iso_gene
	if '_chr' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_chr')]
	elif '_XM' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_XM')]
	elif '_XR' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_XR')]
	elif '_NM' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_NM')]
	elif '_NR' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_NR')]
	elif '_R2_' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_R2_')]
	else:
		iso = iso_gene[:iso_gene.rfind('_')]
	return iso

def update_altsplice_dict(jdict, fiveprime, threeprime, exon_start, exon_end, sample_names,\
	iso_counts, search_threeprime=True):
	if fiveprime not in jdict[chrom]:
		jdict[chrom][fiveprime] = {}  # 5' end anchor if search_threeprime
	if threeprime not in jdict[chrom][fiveprime]:
		jdict[chrom][fiveprime][threeprime] = {}
		jdict[chrom][fiveprime][threeprime]['counts'] = [0]*len(sample_names)
		jdict[chrom][fiveprime][threeprime]['isos'] = []# isoform list for this junction
		jdict[chrom][fiveprime][threeprime]['exon_end'] = exon_end  # for detecting exon skipping
	elif (search_threeprime and strand == '+') or (not search_threeprime and strand == '-'):
		if exon_end < jdict[chrom][fiveprime][threeprime]['exon_end']: # pick shorter exon end
			jdict[chrom][fiveprime][threeprime]['exon_end'] = exon_end
	else:
		if exon_end > jdict[chrom][fiveprime][threeprime]['exon_end']:
			jdict[chrom][fiveprime][threeprime]['exon_end'] = exon_end
	jdict[chrom][fiveprime][threeprime]['isos'] += [name]
	for c in range(len(sample_names)): 
		jdict[chrom][fiveprime][threeprime]['counts'][c] += iso_counts[name][c]
	return jdict

def find_altss(alljuncs, writer, search_threeprime=True):
	""" If fiveprimeon is True, then alternative 5' SS will be reported instead. """
	for chrom in alljuncs:
		for fiveprime in alljuncs[chrom]:
			if len(alljuncs[chrom][fiveprime]) == 1:  # if there is only one 3' end, not an alt 3' junction
				continue

			all_tp = list(alljuncs[chrom][fiveprime].keys())

			n = 0
			for tp1 in all_tp:  # tp1 = three prime SS number 1
				exon_end = alljuncs[chrom][fiveprime][tp1]['exon_end']
				for tp2 in all_tp: # tp2 is also a 3' SS with the same 5' anchor as tp1
					if tp1 == tp2 or abs(tp2-tp1) < wiggle or abs(fiveprime - tp1) > abs(fiveprime - tp2):
						# two sites are the same, too close together, or have already been tested in another order
						continue
					inclusion = tp1
					exclusion = tp2
					strand, chrom_clean = chrom[0], chrom[1:]

					if (search_threeprime and strand == '+') or (not search_threeprime and strand == '-'):
						if tp2 > exon_end:  # exon skipping. tp2 does not overlap tp1's exon
							continue
					elif tp2 < exon_end:  # exon skipping for alt SS upstream of anchor 
							continue

					feature_suffix = chrom_clean+':'+str(fiveprime) if n == 0 else  chrom_clean+':'+str(fiveprime)+'-'+str(n)
					event = chrom_clean+':'+str(fiveprime)+'-'+str(inclusion)+'_'+chrom_clean+':'+str(fiveprime)+'-'+str(exclusion)

					writer.writerow(['inclusion_'+feature_suffix, event] + alljuncs[chrom][fiveprime][inclusion]['counts'] +\
					[','.join(alljuncs[chrom][fiveprime][inclusion]['isos'])] )
					writer.writerow(['exclusion_'+feature_suffix, event] + alljuncs[chrom][fiveprime][exclusion]['counts'] +\
					[','.join(alljuncs[chrom][fiveprime][exclusion]['isos'])] )
					n += 1

iso_counts = {}
sample_names = []
if counts_tsv:
	sample_names = counts_tsv.readline().rstrip().split('\t')[1:]
	for line in counts_tsv:
		line = line.rstrip().split('\t')
		iso = parse_iso_id(line[0])
		iso_counts[iso] = [float(x) for x in line[1:]]

isoforms = {}  # ir detection
ir_junctions = {}  # ir detection
a3_junctions = {}  # alt 3' ss detection
a5_junctions = {}  # alt 5' ss detection
for line in psl:
	line = line.rstrip().split('\t')

	if isbed:
		chrom, name, start, end, strand = line[0], line[3], int(line[1]), int(line[2]), line[5]
	else:
		chrom, name, start, end, strand = line[13], line[9], int(line[15]), int(line[16]), line[8]

	name = parse_iso_id(name)

	if iso_counts and name not in iso_counts:
		continue

	if isbed:
		blockstarts = [int(n) + start for n in line[11].split(',')[:-1]]
		blocksizes = [int(n) for n in line[10].split(',')[:-1]]
	else:
		blocksizes = [int(x) for x in line[18].split(',')[:-1]]
		blockstarts = [int(x) for x in line[20].split(',')[:-1]]

	chrom = strand + chrom  # stranded comparisons
	if chrom not in isoforms:
		isoforms[chrom] = {}
		ir_junctions[chrom] = {}
		a3_junctions[chrom] = {}
		a5_junctions[chrom] = {}

	isoforms[chrom][name] = {}
	isoforms[chrom][name]['sizes']	= blocksizes
	isoforms[chrom][name]['starts']	= blockstarts
	isoforms[chrom][name]['range']	= start, end

	these_jcns = get_junctions_psl(blockstarts, blocksizes)
	for j_index in range(len(these_jcns)):
		j =  these_jcns[j_index]
		fiveprime, threeprime = j[0], j[1]
		exon_end, exon_start = j[3], j[2]

		if strand == '-':
			fiveprime, threeprime = threeprime, fiveprime
			exon_end, exon_start = exon_start, exon_end

		a3_junctions = update_altsplice_dict(a3_junctions, fiveprime, threeprime, \
			exon_start, exon_end, sample_names, iso_counts)
		a5_junctions = update_altsplice_dict(a5_junctions, threeprime, fiveprime, \
			exon_end, exon_start, sample_names, iso_counts, search_threeprime=False)

		j = (j[0], j[1])  # IR junctions do not need the flanking exon info from get_junctions_psl
		if j not in ir_junctions[chrom]:  # ir detection
			ir_junctions[chrom][j] = {}
			ir_junctions[chrom][j]['exclusion'] = {}
			ir_junctions[chrom][j]['inclusion'] = {}
			ir_junctions[chrom][j]['exclusion']['counts'] = [0]*len(sample_names)
			ir_junctions[chrom][j]['inclusion']['counts'] = [0]*len(sample_names)
			ir_junctions[chrom][j]['exclusion']['isos'] = []
			ir_junctions[chrom][j]['inclusion']['isos'] = []
		ir_junctions[chrom][j]['exclusion']['isos'] += [name]
		for c in range(len(sample_names)):
			ir_junctions[chrom][j]['exclusion']['counts'][c] += iso_counts[name][c]

with open(outfilenamebase+'.alt3.events.quant.tsv', 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	writer.writerow(['feature_id', 'coordinate']+sample_names+['isoform_ids'])
	find_altss(a3_junctions, writer)

with open(outfilenamebase+'.alt5.events.quant.tsv', 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	writer.writerow(['feature_id', 'coordinate']+sample_names+['isoform_ids'])
	find_altss(a5_junctions, writer, search_threeprime=False)

seen_j = set()
with open(outfilenamebase + '.ir.events.quant.tsv', 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	writer.writerow(['feature_id', 'coordinate']+sample_names+['isoform_ids'])
	for chrom in ir_junctions:
		for j in ir_junctions[chrom]:
			for iname in isoforms[chrom]:  # compare with all other isoforms to find IR
				if iname in ir_junctions[chrom][j]['exclusion']['isos']:  # is an exclusion isoform
					continue
				start, end = isoforms[chrom][iname]['range']
				if start > j[1] or end < j[0]:  # isoform boundaries do not overlap junction
					continue
				starts, sizes = isoforms[chrom][iname]['starts'], isoforms[chrom][iname]['sizes']
				for start, size in zip(starts[1:], sizes[1:]):
					estart, eend = start, start+size  # exon start, exon end
					if estart < j[0] and eend > j[1]:  # retention
						ir_junctions[chrom][j]['inclusion']['isos'] += [iname]
						for c in range(len(sample_names)):
							ir_junctions[chrom][j]['inclusion']['counts'][c] += iso_counts[iname][c]

		for j in ir_junctions[chrom]:
			incounts = ir_junctions[chrom][j]['inclusion']['counts']
			if sum(incounts) == 0:
				continue
			if not sample_names:
				ir_junctions[chrom][j]['exclusion']['counts'] = ir_junctions[chrom][j]['inclusion']['counts'] = []

			chrom_clean = chrom[1:]
			event = chrom_clean+':'+str(j[0])+'-'+str(j[1])
			writer.writerow(['inclusion_'+event+chrom[0], event] + ir_junctions[chrom][j]['inclusion']['counts'] +\
			[','.join(ir_junctions[chrom][j]['inclusion']['isos'])] )
			writer.writerow(['exclusion_'+event+chrom[0], event] + ir_junctions[chrom][j]['exclusion']['counts'] +\
			[','.join(ir_junctions[chrom][j]['exclusion']['isos'])] )
		ir_junctions[chrom] = None



