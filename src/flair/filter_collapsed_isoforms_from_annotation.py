#!/usr/bin/env python3
import sys, csv, math, os, argparse

parser = argparse.ArgumentParser(description='options',
	usage='python filter_collapsed_isoforms_from_annotation.py -i in.bed|psl -a annotated.bed -o out [options]')
required = parser.add_argument_group('required named arguments')
required.add_argument('-i', action='store', dest='i', type=str, required=True,
	help='input bed/psl of collapsed isoforms to be filtered')
required.add_argument('-a', action='store', dest='a', type=str, required=True,
	help='annotated .bed of isoforms with read support')
required.add_argument('-o', action='store', dest='o', type=str, required=True,
	help='output file, .bed or .psl matching the input file')
parser.add_argument('--new_map', action='store', dest='new_map', type=str, default='',
	help='output annotated map file for isos there were merged')
parser.add_argument('--map_i', action='store', dest='map_i', type=str, required=False, default='', 
	help='''isoform-read map file for the annotated isoforms to retain subset isoforms with 
	sufficient read support. must also specify --map_a.''')
parser.add_argument('--map_a', action='store', dest='map_a', type=str, required=False, default='', 
	help='''isoform-read map file for the flair-collapse isoforms)''')
parser.add_argument('-w', action='store', dest='w', type=int, required=False, default=100, 
	help='''number of extra basepairs on a terminal exon for a subset isoform to be kept (default=100)''')
parser.add_argument('-s', '--support', default=3, action='store', dest='s', type=float,
	help='minimum number of supporting reads for an isoform (3)')

args = parser.parse_args()

try:
	psl = open(args.i)
	annotated = open(args.a)
except:
	sys.stderr.write('Issue with opening files (filter_collapsed_isoforms_from_annotation.py)\n')
	sys.exit(1)

isbed = args.i[-3:].lower() != 'psl'
tol = args.w

def split_iso_gene(iso_gene):
    if '_' not in iso_gene:
    	return iso_gene, 'NA'
    elif '_chr' in iso_gene:
        splitchar = '_chr'
    elif '_XM' in iso_gene:
        splitchar = '_XM'
    elif '_XR' in iso_gene:
        splitchar = '_XR'
    elif '_NM' in iso_gene:
        splitchar = '_NM'
    elif '_NR' in iso_gene:
        splitchar = '_NR'
    elif '_R2_' in iso_gene:
        splitchar = '_R2_'
    elif '_NC_' in iso_gene:
        splitchar = '_NC_'
    else:
        splitchar = '_'
    iso = iso_gene[:iso_gene.rfind(splitchar)]
    gene = iso_gene[iso_gene.rfind(splitchar)+1:]
    return iso, gene

def get_info(line, isbed):
	if isbed:
		chrom, name, chrstart = line[0], line[3], int(line[1])
		sizes = [int(n) for n in line[10].split(',')[:-1]]
		starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	else:
		chrom, name = line[13], line[9]
		sizes = [int(n) for n in line[18].split(',')[:-1]]
		starts = [int(n) for n in line[20].split(',')[:-1]]
	return chrom, name, sizes, starts

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

iso_support = {}
annotated_iso_read_map = {}
flair_iso_read_map = {}
if args.map_a:
	for line in open(args.map_a):
		iso, reads = line.rstrip().split('\t')
		num_supporting = len(reads.split(','))
		if num_supporting < int(args.s):
			continue
		iso_support[iso] = num_supporting
		annotated_iso_read_map[iso] = reads
	for line in open(args.map_i):
		iso, reads = line.rstrip().split('\t')
		# if '_' in iso:
		# 	iso, gene = split_iso_gene(iso)
		flair_iso_read_map[iso] = reads
		if iso in iso_support:
			iso_support[iso] += len(reads.split(','))
		else:
			iso_support[iso] = len(reads.split(','))

keep_isoforms = []
isoforms, allevents, jcn_to_name, all_iso_info = {}, {}, {}, {}
for line in annotated:
	line = line.rstrip().split()
	chrom, name, sizes, starts = get_info(line, isbed)
	keep_isoforms += [name]
	junctions = get_junctions(starts, sizes)
	exons = get_exons(starts, sizes)
	if chrom not in isoforms:
		isoforms[chrom] = {}
		allevents[chrom] = {}
		allevents[chrom]['allexons'] = set()
		jcn_to_name[chrom] = {}
		allevents[chrom]['all_se_exons'] = set()
	isoforms[chrom][name] = {}
	isoforms[chrom][name]['exons'] = exons
	all_iso_info[name] = line
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
		allevents[chrom]['all_se_exons'].add((exon[0], exon[1], iso_support[name]))

in_annotation = 0
for line in psl:
	line = line.rstrip().split()
	chrom, name, sizes, starts = get_info(line, isbed)
	# if split_iso_gene(name)[0] in keep_isoforms:
	# 	in_annotation += 1
	# 	continue
	junctions = get_junctions(starts, sizes)
	exons = get_exons(starts, sizes)
	if chrom not in isoforms:
		isoforms[chrom] = {}
		allevents[chrom] = {}
		allevents[chrom]['allexons'] = set()
		jcn_to_name[chrom] = {}
		allevents[chrom]['all_se_exons'] = set()
	isoforms[chrom][name] = {}
	isoforms[chrom][name]['exons'] = exons
	all_iso_info[name] = line
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
		allevents[chrom]['all_se_exons'].add((exon[0], exon[1], iso_support[name]))

for chrom in isoforms:
	allexons_left = sorted(list(allevents[chrom]['allexons']), key=lambda x: x[0])
	allexons_right = sorted(allexons_left, key=lambda x: x[1])
	all_se_exons_left = sorted(list(allevents[chrom]['all_se_exons']), key=lambda x: x[0])
	all_se_exons_right = sorted(all_se_exons_left, key=lambda x: x[1])
	for n in isoforms[chrom]:  # n is isoform name
		if n in keep_isoforms:  # is an annotated isoform and will automatically be kept
			continue

		if 'jname' in isoforms[chrom][n]:  # multi-exon isoforms

			similar_isos = set()  # similar isos share any splice jcns
			for j in isoforms[chrom][n]['jset']:
				similar_isos.update(jcn_to_name[chrom][j])

			issubset = [False, False]  # first exon is a subset, last exon is a subset
			exons = isoforms[chrom][n]['exons']
			first_exon, last_exon = exons[0], exons[-1]
			superset_iso = set()
			# if '373715' in n:
				# print(n)
			too_similar_to_annotated = False
			for n_ in similar_isos:
				if n == n_:
					continue
				# if isoform n and isoform n_ have the exact same splice junction chain and
				# the ends of n are w/in 50 bp of the ends of n_, they are too close, do not keep
				# if n_ is already in keep_isos
				if isoforms[chrom][n]['jname'][1:-1] in isoforms[chrom][n_]['jname'] and \
				len(isoforms[chrom][n]['jname']) == len(isoforms[chrom][n_]['jname']) and \
				abs(isoforms[chrom][n_]['exons'][0][0]-first_exon[0]) < args.w and \
				abs(isoforms[chrom][n_]['exons'][-1][1]-last_exon[1]) < args.w and \
				n_ in keep_isoforms:
					too_similar_to_annotated = True
					if n_ not in annotated_iso_read_map:
						print(n, 'too similar to', n_)
						print(isoforms[chrom][n_]['exons'][0])
					else:
						annotated_iso_read_map[n_] += ','+flair_iso_read_map[n]
					break
				# are all the junctions in this isoform n in isoform n_ and
				# are there the same number or more splice junctions in n_
				if isoforms[chrom][n]['jname'][1:-1] in isoforms[chrom][n_]['jname'] and \
				len(isoforms[chrom][n]['jname']) <= len(isoforms[chrom][n_]['jname']):
					for exon in isoforms[chrom][n_]['exons']:
						if exon_overlap(first_exon, exon, tol=args.w):
							issubset[0] = True
							superset_iso.add(n_)
						elif exon_overlap(last_exon, exon, left=False, tol=args.w):
							issubset[1] = True
							superset_iso.add(n_)

			if not too_similar_to_annotated:
				if not issubset[0] and not issubset[1]:  # not a subset
					keep_isoforms += [n]
				elif iso_support[n] > 3:  
					superset_support = []
					for n_ in superset_iso:
						superset_support += [iso_support[n_]]
					# subset isoform will be kept if it has higher expression of the average of its superset isoforms
					if iso_support[n] > (sum(superset_support)/len(superset_support))*0.5:
						keep_isoforms += [n]

		else:  # single exon isoforms

			exon = isoforms[chrom][n]['exons'][0]
			iscontained = False  # is exon contained in another exon (e) of a multi-exon transcript?
			b = bin_search_left(exon, allexons_right)
			c = bin_search_right(exon, allexons_left)
			for e in set(b).intersection(set(c)):
				if contained(exon, (e[0]-args.w, e[1]+args.w)):
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
				if contained(exon, (e[0]-args.w, e[1]+args.w)):
					iscontained = True
					# if iso_support[n] > 3 and \
					# 	e[2] * 1.2 < iso_support[n]:
					# 	contained_expr = True  # keep exon
			if not iscontained or contained_expr:
				keep_isoforms += [n]


# used_names = set()
# name_remapping = {}
# with open(args.o, 'wt') as outfile:  # renamed isoform bed
# 	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)

# 	for iso in keep_isoforms:  # isoform bed
# 		new_name, gene = split_iso_gene(iso)
# 		if gene != 'NA':
# 			if new_name in used_names:
# 				if '-' in new_name[-3:]:
# 					sys.stderr.write('Issue, exiting {}\n'.format(iso))
# 					sys.exit(1)
# 				new_name += '-0'
# 		used_names.add(new_name)

# 		line = all_iso_info[iso]
# 		if isbed:
# 			line[3] = new_name + '_' + gene
# 		else:
# 			line[9] = new_name + '_' + gene
# 		writer.writerow(line)
# 		name_remapping[iso] = new_name

# if args.new_map:
# 	with open(args.new_map, 'wt') as outfile:
# 		writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
# 		# writer.writerow(['transcript_id', 'read_ids'])
# 		for iso in annotated_iso_read_map:
# 			writer.writerow([name_remapping[iso], annotated_iso_read_map[iso]])

# 		for iso in keep_isoforms:
# 			if iso not in annotated_iso_read_map:  # flair collapse iso that is in keep_isoforms
# 				writer.writerow([name_remapping[iso], flair_iso_read_map[iso]])


with open(args.o, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for iso in keep_isoforms:
		writer.writerow(all_iso_info[iso])

if args.new_map:
	with open(args.new_map, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
		for iso in annotated_iso_read_map:
			writer.writerow([iso, annotated_iso_read_map[iso]])
		for iso in keep_isoforms:
			if iso not in annotated_iso_read_map:  # flair collapse iso that is in keep_isoforms
				writer.writerow([iso, flair_iso_read_map[iso]])





