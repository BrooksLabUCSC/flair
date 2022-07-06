#!/usr/bin/env python3
import sys, csv, math, os

try:
	psl = open(sys.argv[1])
	isbed = sys.argv[1][-3:].lower() != 'psl' 
	mode = sys.argv[2] if sys.argv[2] in ['ginormous', 'comprehensive', 'nosubset'] else 'default' # all typos will report default
	pslout = sys.argv[3]
	tol = 100 if len(sys.argv) <= 4 else int(sys.argv[4])
	keep_extra_column = len(sys.argv) > 5
except:
	sys.stderr.write('usage: filter_collapsed_isoforms.py collapsed.psl (nosubset/default/comprehensive/ginormous)' + \
		' filtered.psl [tolerance] [keep_extra_column]\n')
	sys.exit(1)

def get_junctions(starts, sizes):
	junctions = set()
	if len(starts) != 1:
		for b in range(len(starts)-1):
			junctions.add((starts[b]+sizes[b], starts[b+1]))
		return junctions

def get_exons(starts, sizes):
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

def contained(coords0, coords1):  # complete coverage of coords0 by coords1
	return coords1[0] <= coords0[0] and coords1[1] >= coords0[1]

def bin_search_left(query, data):
	""" Query is a coordinate interval. Approximate binary search for the left coordinate of 
	the query in data sorted by the right coordinate. Finishes when the first interval in data with
	a right coordinate that is greater than the query's left coordinate is found. """
	lower, upper, i = 0, len(data), int(math.floor(len(data)/2))  # binary search prep
	if upper == 0:
		return set()
	while True:
		if lower == i or upper == i:
			if i - 1 > 0 and data[i-1][1] >= query[0]:
				lower = i = i - 1
			else:
				if data[i][1] < query[0]:
					i = i + 1
				break
		elif data[i][1] < query[0]:  # i is too low of an index
			lower = i
			i = int(math.floor((lower + upper)/2.))  # floor because of python 2 and 3 differences
		else: # i is too high of an index
			upper = i
			i = int(math.floor((lower + upper)/2.))
	return data[i:min(i+40, len(data))]

def bin_search_right(query, data):
	""" Query is a coordinate interval. Approximate binary search for the left coordinate of 
	the query in data sorted by the right coordinate. Finishes when the first interval in data with
	a left coordinate that is greater than the query's right coordinate is found. Kept separate for
	my own readability. """
	lower, upper, i = 0, len(data), int(math.floor(len(data)/2))
	if upper == 0:
		return set()
	while True:
		if lower == i or upper == i:
			if i + 1 < len(data) and data[i+1][0] <= query[1]:
				upper = i = i + 1
			else:
				break
		elif data[i][0] < query[1]:
			lower = i
			i = int(math.floor((lower + upper)/2.))
		else:
			upper = i
			i = int(math.floor((lower + upper)/2.))
	return data[max(0, i-40):i+1]

isoforms, allevents, jcn_to_name = {}, {}, {}
for line in psl:
	line = line.rstrip().split()
	try:
		line[-1] = float(line[-1])
		support = True
	except:
		sys.stderr.write('Expecting an extra column of integers\n')
		sys.exit(1)
	if isbed:
		chrom, name, chrstart = line[0], line[3], int(line[1])
		sizes = [int(n) for n in line[10].split(',')[:-1]]
		starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	else:
		chrom, name = line[13], line[9]
		sizes = [int(n) for n in line[18].split(',')[:-1]]
		starts = [int(n) for n in line[20].split(',')[:-1]]
	junctions = get_junctions(starts, sizes)
	exons = get_exons(starts, sizes)
	if chrom not in isoforms:
		isoforms[chrom] = {}
		allevents[chrom] = {}
		allevents[chrom]['allexons'] = set()
		jcn_to_name[chrom] = {}
		allevents[chrom]['all_se_exons'] = set()
	isoforms[chrom][name] = {}
	isoforms[chrom][name]['line'] = line
	isoforms[chrom][name]['exons'] = exons
	if junctions:  # multi-exon isoform
		isoforms[chrom][name]['jset'] = junctions
		allevents[chrom]['allexons'].update(exons)
		jname = str(sorted(list(junctions)))
		isoforms[chrom][name]['jname'] = jname
		for j in junctions:
			if j not in jcn_to_name[chrom]:
				jcn_to_name[chrom][j] = set()
			jcn_to_name[chrom][j].add(name)
	else:
		exon = exons[0]
		allevents[chrom]['all_se_exons'].add((exon[0], exon[1], line[-1]))

keepisoforms = []
for chrom in isoforms:
	allexons_left = sorted(list(allevents[chrom]['allexons']), key=lambda x: x[0])
	allexons_right = sorted(allexons_left, key=lambda x: x[1])
	all_se_exons_left = sorted(list(allevents[chrom]['all_se_exons']), key=lambda x: x[0])
	all_se_exons_right = sorted(all_se_exons_left, key=lambda x: x[1])
	for n in isoforms[chrom]:
		if mode == 'ginormous':  # currently includes all
			keepisoforms += [isoforms[chrom][n]['line']]
			continue
		if 'jname' in isoforms[chrom][n]:  # multi-exon isoforms
			if mode == 'comprehensive':  # all spliced subsets allowed
				keepisoforms += [isoforms[chrom][n]['line']]
				continue
			similar_isos = set()
			for j in isoforms[chrom][n]['jset']:
				similar_isos.update(jcn_to_name[chrom][j])

			issubset = [0, 0]  # first exon is a subset, last exon is a subset
			exons = isoforms[chrom][n]['exons']
			first_exon, last_exon = exons[0], exons[-1]
			superset_support = []
			for n_ in similar_isos:  # compare with all similar
				if n == n_:
					continue
				elif isoforms[chrom][n]['jname'][1:-1] in isoforms[chrom][n_]['jname'] and \
				len(isoforms[chrom][n]['jname']) < len(isoforms[chrom][n_]['jname']):  # is subset
					for exon in isoforms[chrom][n_]['exons']:
						if exon_overlap(first_exon, exon, tol=tol):
							superset_support += [isoforms[chrom][n_]['line'][-1]]
							issubset[0] = 1
						if exon_overlap(last_exon, exon, left=False, tol=tol):
							superset_support += [isoforms[chrom][n_]['line'][-1]]
							issubset[1] = 1
			if sum(issubset) == 0:  # not a subset
				keepisoforms += [isoforms[chrom][n]['line']]
			elif mode == 'nosubset' or len(exons) < 4:  # is a subset and will be removed
				continue
			elif isoforms[chrom][n]['line'][-1] > 3 and \
				isoforms[chrom][n]['line'][-1] > (sum(superset_support)/len(superset_support))*1.2:
					keepisoforms += [isoforms[chrom][n]['line']]
		else:  # single exon isoforms
			exon = isoforms[chrom][n]['exons'][0]
			iscontained = False  # is exon contained in another exon (e) of a multi-exon transcript?
			b = bin_search_left(exon, allexons_right)
			c = bin_search_right(exon, allexons_left)
			for e in set(b).intersection(set(c)):
				if contained(exon, (e[0]-tol, e[1]+tol)):
					iscontained = True
					break
			if iscontained:  # if so, then filter out
				continue

			contained_expr = False  # exon is contained in another exon but is more highly expressed
			b = bin_search_left(exon, all_se_exons_right)
			c = bin_search_right(exon, all_se_exons_left)
			for e in set(b).intersection(set(c)):  # is exon contained in another single exon transcript e?
				if (e[1] - e[0]) <= (exon[1] - exon[0]):  # e must be larger than exon to remove exon
					continue
				if contained(exon, (e[0]-tol, e[1]+tol)):
					iscontained = True
					if mode != 'nosubset' and isoforms[chrom][n]['line'][-1] > 3 and \
						e[2] * 1.2 < isoforms[chrom][n]['line'][-1]:
						contained_expr = True  # keep exon
			if not iscontained or contained_expr:
				keepisoforms += [isoforms[chrom][n]['line']]

with open(pslout, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for iso in keepisoforms:
		if keep_extra_column:
			writer.writerow(iso)
		else:
			writer.writerow(iso[:-1])
