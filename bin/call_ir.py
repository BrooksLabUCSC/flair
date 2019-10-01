import sys, csv

try:
	psl = open(sys.argv[1])
	isgtf = sys.argv[1][-3:] == 'gtf'
	isbed = sys.argv[1][-3:] == 'bed'
	outfilename = sys.argv[2]
	if len(sys.argv) > 3:
		counts_tsv = open(sys.argv[3])
	else:
		counts_tsv = ''
except:
	sys.stderr.write('usage: script.py .psl|.bed out.tsv [counts_tsv]\n')
	sys.exit(1)

def overlap(coords0, coords1):
	return coords1[0] >= coords0[0] and coords1[0] <= coords0[1] or \
		coords1[1] >= coords0[0] and coords1[1] <= coords0[1]

def get_junctions_psl(starts, sizes):
	junctions = set()
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

iso_counts = {}
sample_names = []
if counts_tsv:
	sample_names = counts_tsv.readline().rstrip().split('\t')[1:]
	for line in counts_tsv:
		line = line.rstrip().split('\t')
		iso = line[0][:line[0].rfind('_')]
		iso_counts[iso] = [float(x) for x in line[1:]]

isoforms = {}
junctions = {}

for line in psl:
	line = line.rstrip().split('\t')
	if isbed:
		chrom, name, start, end, strand = line[0], line[3], int(line[1]), int(line[2]), line[5]
	else:
		chrom, name, start, end, strand = line[13], line[9], int(line[15]), int(line[16]), line[8]
	if '_' in name:
		name = name[:name.rfind('_')]
	if iso_counts and name not in iso_counts:
		continue
	if isbed:
		blockstarts = [int(n) + start for n in line[11].split(',')[:-1]]
		blocksizes = [int(n) for n in line[10].split(',')[:-1]]
	else:
		blocksizes = [int(x) for x in line[18].split(',')[:-1]]
		blockstarts = [int(x) for x in line[20].split(',')[:-1]]
	if chrom not in isoforms:
		isoforms[chrom] = {}
		junctions[chrom] = {}
	if strand not in isoforms[chrom]:
		isoforms[chrom][strand] = {}
		junctions[chrom][strand] = {}

	isoforms[chrom][strand][name] = {}
	isoforms[chrom][strand][name]['entry'] = line
	isoforms[chrom][strand][name]['sizes'] = blocksizes
	isoforms[chrom][strand][name]['starts'] = blockstarts
	isoforms[chrom][strand][name]['range'] = start, end

	these_jcns = get_junctions_psl(blockstarts, blocksizes)
	for j in these_jcns:
		if j not in junctions[chrom][strand]:
			junctions[chrom][strand][j] = {}
			junctions[chrom][strand][j]['exclusion'] = {}
			junctions[chrom][strand][j]['inclusion'] = {}
			junctions[chrom][strand][j]['exclusion']['counts'] = [0]*len(sample_names)
			junctions[chrom][strand][j]['inclusion']['counts'] = [0]*len(sample_names)
			junctions[chrom][strand][j]['exclusion']['isos'] = []
			junctions[chrom][strand][j]['inclusion']['isos'] = []

		junctions[chrom][strand][j]['exclusion']['isos'] += [name]
		for c in range(len(sample_names)):
			junctions[chrom][strand][j]['exclusion']['counts'][c] += iso_counts[name][c]

introncoords = set()
allcoords = set()
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	writer.writerow(['feature_id', 'gene_id']+sample_names+['isoform_ids'])
	n = 0
	for chrom in junctions:
		for strand in junctions[chrom]:
			for j in junctions[chrom][strand]:
				for iname in isoforms[chrom][strand]:  # compare with all other isoforms to find IR
					if iname in junctions[chrom][strand][j]['exclusion']['isos']:  # is an exclusion isoform
						continue
					start, end = isoforms[chrom][strand][iname]['range']
					if start > j[1] or end < j[0]:  # isoform boundaries do not overlap junction
						continue
					starts, sizes = isoforms[chrom][strand][iname]['starts'], isoforms[chrom][strand][iname]['sizes']
					for start, size in zip(starts[1:], sizes[1:]):
						estart, eend = start, start+size  # exon start, exon end
						if estart < j[0] and eend > j[1]:  # retention
							junctions[chrom][strand][j]['inclusion']['isos'] += [iname]
							for c in range(len(sample_names)):
								junctions[chrom][strand][j]['inclusion']['counts'][c] += iso_counts[iname][c]

			for j in junctions[chrom][strand]:
				incounts = junctions[chrom][strand][j]['inclusion']['counts']
				if sum(incounts) == 0:
					continue
				if not sample_names:
					junctions[chrom][j][strand]['exclusion']['counts'] = junctions[chrom][strand][j]['inclusion']['counts'] = []

				event = chrom+':'+str(j[0])+'-'+str(j[1])
				writer.writerow(['exclusion_'+str(n), event] + junctions[chrom][strand][j]['exclusion']['counts'] +\
				[','.join(junctions[chrom][strand][j]['exclusion']['isos'])] )
				writer.writerow(['inclusion_'+str(n), event] + junctions[chrom][strand][j]['inclusion']['counts'] +\
				[','.join(junctions[chrom][strand][j]['inclusion']['isos'])] )
				n += 1	
			junctions[chrom][strand] = 0
