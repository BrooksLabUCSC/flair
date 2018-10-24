import sys, csv

try:
	psl = open(sys.argv[1])
	mode = sys.argv[2]
	if sys.argv[2] != 'ginormous' and sys.argv[2] != 'comprehensive':
		mode = 'default'  # all typos will report default
	pslout = sys.argv[3]
except:
	print('usage: script.py collapsed.psl [default/comprehensive/ginormous] filtered.psl')
	sys.exit()

def get_junctions_psl(starts, sizes):
	junctions = set()
	if len(starts) != 1:
		for b in range(len(starts)-1):
			junctions.add((starts[b]+sizes[b], starts[b+1]))
		return junctions

def get_junctions_bed12(line):
	junctions = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_exons_psl(starts, sizes):
	exons = []
	for e in range(len(starts)):
		exons += [(starts[e], starts[e]+sizes[e])]
	return exons

def overlap(coords0, coords1, tol=1):
	coords0, coords1 = sorted([coords0, coords1], key = lambda x: x[0])
	return (coords0[0] < coords1[0] and coords1[0] < coords0[1] - tol)

def exon_overlap(coords0, coords1, left=True, tol=1):
	len0, len1 = coords0[1] - coords0[0], coords1[1] - coords1[0]
	if left:
		return coords0[1] == coords1[1] and (len0 - tol) < len1
	return coords0[0] == coords1[0] and (len0 - tol) < len1

def contained(coords0, coords1):
	return coords1[0] <= coords0[0] and coords1[1] >= coords0[1]  # complete coverage of coords0 by coords1

def bin_search(query, data):
	""" Query is coordinates. Binary search for the query in data, which is a list of coordinates.
	Finishes when an overlapping value of query and data exists and returns the index in data. """
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
	
isoforms = {}
allevents = {}
jcn_to_name = {}
if sys.argv[1][-3:] == 'bed':
	for line in psl:
		line = line.rstrip().split()

else:
	for line in psl:
		line = line.rstrip().split()
		chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])
		sizes = [int(n) for n in line[18].split(',')[:-1]]
		starts = [int(n) for n in line[20].split(',')[:-1]]
		junctions = get_junctions_psl(starts, sizes)
		exons = get_exons_psl(starts, sizes)
		jname = str(sorted(list(junctions)))
		if chrom not in isoforms:
			isoforms[chrom] = {}
			allevents[chrom] = {}
			allevents[chrom]['alljunctions'] = set()
			allevents[chrom]['allexons'] = set()
			jcn_to_name[chrom] = {}
		# if '-' in name[-4:]:
		# 	name = name[:name.rfind('-')]
		if name not in isoforms[chrom]:
			isoforms[chrom][name] = {}
			isoforms[chrom][name]['line'] = [line]
			isoforms[chrom][name]['jset'] = junctions
			isoforms[chrom][name]['jname'] = jname
			isoforms[chrom][name]['exons'] = exons
			allevents[chrom]['alljunctions'].update(junctions)
			allevents[chrom]['allexons'].update(exons)
			for j in junctions:
				if j not in jcn_to_name[chrom]:
					jcn_to_name[chrom][j] = set()
				jcn_to_name[chrom][j].add(name)

keepisoforms = []
for chrom in isoforms:
	sys.stderr.write(chrom+'\n')
	alljunctions = sorted(list(allevents[chrom]['alljunctions']), key=lambda x: x[0])
	allexons = sorted(list(allevents[chrom]['allexons']), key=lambda x: x[0])
	for n in isoforms[chrom]:
		if mode == 'ginormous':  # currently includes all
			keepisoforms += isoforms[chrom][n]['line']
			continue
		if len(isoforms[chrom][n]['jset']) >= 1:  # multi-exon isoforms
			if mode == 'comprehensive':  # all spliced subsets allowed
				keepisoforms += isoforms[chrom][n]['line']
				continue
			similar_isos = set()
			for j in isoforms[chrom][n]['jset']:
				similar_isos.update(jcn_to_name[chrom][j])

			issubset = [0, 0]  # first exon is a subset, last exon is a subset
			subset_juncs = set()  # if true, will be populated with introns of other isoforms
			exons = sorted(list(isoforms[chrom][n]['exons']),key=lambda x: x[0])
			first_exon, last_exon = exons[0], exons[1]
			for n_ in similar_isos:  # compare with all similar
				if n == n_:
					continue
				elif isoforms[chrom][n]['jname'][1:-1] in isoforms[chrom][n_]['jname'] and \
				len(isoforms[chrom][n]['jname']) < len(isoforms[chrom][n_]['jname']):
					for exon in list(isoforms[chrom][n_]['exons']):
						if exon_overlap(first_exon, exon, tol=10):
							issubset[0] = 1
						elif exon_overlap(last_exon, exon, left=False, tol=10):
							issubset[1] = 1
					subset_juncs.update(isoforms[chrom][n_]['jset'])
					
			if sum(issubset) != 2:  # see if first/last exon overlaps a junction 
				ir = False
				for exon in sorted(list(isoforms[chrom][n]['exons'])):
					for junction in sorted(list(subset_juncs)):
						if overlap(exon, junction, 10):
							ir = True
							break
				if ir:  # if there is intron retention, it is not considered a subset isoform
					keepisoforms += isoforms[chrom][n]['line']
		else:  # single exon isoforms
			exon = list(isoforms[chrom][n]['exons'])[0]
			i = bin_search(exon, alljunctions)
			iscontained = False 
			for j in alljunctions[i-2:i+2]:
				if contained(exon, j):
					iscontained = True
					break
			if iscontained:  # boolean: exon is completely contained within a splice junction
				keepisoforms += isoforms[chrom][n]['line']
				continue

			i = bin_search(exon, allexons)
			iscontained = False  # boolean: exon completely contained in another exon
			for e in allexons[i-2:i+2]:
				if contained(exon, e):
					iscontained = True
					break
			if not iscontained:
				keepisoforms += isoforms[chrom][n]['line']


with open(pslout, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for iso in keepisoforms:
		writer.writerow(iso)
