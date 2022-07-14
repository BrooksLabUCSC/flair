#!/usr/bin/env python3

import sys, csv, argparse, math, os
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='collapse parse options', \
			usage='python collapse_isoforms_precise.py -q <query.psl>/<query.bed> [options]')
required = parser.add_argument_group('required named arguments')
required.add_argument('-q', '--query', type=str, default='', required=True, action='store', \
	dest='q', help='BED12 or PSL file of aligned/corrected reads. PSL should end in .psl')
parser.add_argument('-o', '--output', type=str, action='store', \
	dest='o', default='', help='specify output file, should agree with query file type')
parser.add_argument('-w', '--window', default=100, type=int, \
	action='store', dest='w', help='window size for grouping TSS/TES (100)')
parser.add_argument('-s', '--support', default=0.25, type=float, action='store', \
	dest='s', help='minimum proportion(s<1) or number of supporting reads(s>=1) for an isoform (0.1)')
parser.add_argument('-f', '--gtf', default='', type=str, \
	action='store', dest='f', help='GTF annotation file for selecting annotated TSS/TES')
parser.add_argument('-m', '--max_results', default=2, type=int, action='store', dest='m', \
	help='maximum number of novel TSS or TES picked per isoform unless --no_redundant is specified (2)')
parser.add_argument('-t', '--threads', default=2, type=int, \
	action='store', dest='t', help='number of threads to use (2)')
parser.add_argument('-n', '--no_redundant', default='none', action='store', dest='n', \
	help='For each unique splice junction chain, report options include: \
	none: multiple best TSSs/TESs chosen for each unique set of splice junctions, see M; \
	longest: TSS/TES chosen to maximize length; \
	best_only: single best TSS/TES used in conjunction chosen; \
	longest/best_only override max_results argument immediately before output \
	resulting in one isoform per unique set of splice junctions (default: none)')
parser.add_argument('-c', '--clean', default=False, action='store_true', dest='c', \
	help='Specify this to not append read support to the end of each entry (default: not specified)')
parser.add_argument('-i', '--isoformtss', default=False, action='store_true', dest='i', \
	help='when specified, TSS/TES for each isoform will be determined from supporting reads \
	for individual isoforms (default: not specified, determined at the gene level)')
parser.add_argument('-nosplice', default='chrM', action='store', dest='nosplice', \
	help='Comma separated list of chromosomes that should not have spliced isoforms (default: chrM)')
parser.add_argument('--quiet', default=False, \
	action='store_true', dest='quiet', help='suppress output to stderr')
args = parser.parse_args()

try:
	max_results, window, minsupport, psl = args.m, args.w, args.s, open(args.q)
	if args.f:
		gtf = open(args.f)
except:
	sys.stderr.write('Make sure all files (query, GTF) have valid paths and can be opened\n')
	sys.exit(1)

bed = args.q[-3:].lower() != 'psl'
if args.o:
	if args.o[-3:].lower() != args.q[-3:].lower():
		sys.stderr.write('Make sure input and output file extensions agree\n')
		sys.exit(1)
else:  # default output name
	args.o = args.q[:-3]+'collapsed.bed' if bed else args.q[:-3]+'collapsed.psl'

def get_junctions(line):
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	if len(starts) == 1:
		return 0, 0
	junctions = set()
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions, (starts[0]+sizes[0], starts[-1])

def get_junctions_bed12(line):
	chrstart = int(line[1])
	starts = [int(n) + chrstart + 1 for n in line[11].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return 0, 0
	junctions = set()
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions, (starts[0]+sizes[0], starts[-1])

def overlap(coord0, coord1, tol=0):
	if coord1[0] < coord0[0]:
		coord0, coord1 = coord1, coord0
	isoverlap = coord0[0] <= coord1[0] and coord1[0] <= coord0[1] - tol
	if isoverlap:
		return isoverlap, (min(coord0[1], coord1[1]) - coord1[0])/(coord0[1] - coord0[0])# + \
			#(coord0[1] - coord1[0])/(coord1[1] - coord1[0]) / 2
	return isoverlap, 0

def bin_search(query, data):
	""" Query is a coordinate interval. Approximate binary search for the query in sorted data, 
	which is a list of coordinates. Finishes when the closest overlapping value of query and 
	data is found and returns the index in data. """
	i = int(math.floor(len(data)/2))  # binary search prep
	lower, upper = 0, len(data)
	if not upper:
		return -1
	tried = set()
	rightfound = ''  # null value in place of 0, which is a valid value for rightfound
	while not (data[i][0] <= query[0] and data[i][1] >= query[0]):  # query left coordinate not found in data yet
		if data[i][0] <= query[1] and data[i][1] >= query[1]:  # query right found, will keep looking for left 
			rightfound = i
		if data[i][1] < query[0]:  # i is too low of an index
			lower = i
			i = int(math.floor((lower + upper)/2.))
		else:  # i is too high of an index
			upper = i
			i = int(math.floor((lower + upper)/2.))

		if i in tried or i == upper:
			if data[i][0] >= query[0] and data[i][1] <= query[1]: # data interval sandwiched inside query
				break
			elif i + 1 < len(data) and data[i+1][0] > query[0] and data[i+1][1] < query[1]:  # data can be incremented
				i = i + 1
			else:
				i = rightfound if rightfound != '' else -1
			break
		tried.add(i)
	return i

def interval_insert(query, data):
	i = int(round(len(data)/2))
	lower, upper = 0, len(data)
	if not upper:
		return [query]
	while True:
		if data[i][0] <= query[0]:  # i is too low of an index
			lower = i
			i = int(math.floor((lower + upper)/2.))
			if i == lower:
				break
		elif data[i][0] >= query[0]:  # i is too high of an index
			upper = i
			i = int(math.floor((lower + upper)/2.))
			if i == upper:
				i = -1
				break
	data = data[:i+1] + [query] + data[i+1:]
	return data

def add_se(sedict, tss, tes, line, sedictkeys, support=1):
	loci = {}  # altered loci
	c = bin_search((tss, tes), sedictkeys)
	if c != -1:  # found
		for coord in sedictkeys[c:]:  # coordinates of a locus of single-exon reads
			if overlap((tss, tes), coord):
				if tss not in sedict[coord]['tss']:
					sedict[coord]['tss'][tss] = 0
				if tss not in sedict[coord]['tss_tes']:
					sedict[coord]['tss_tes'][tss] = {}
				if tes not in sedict[coord]['tss_tes'][tss]:
					sedict[coord]['tss_tes'][tss][tes] = 0
				if not loci:
					sedict[coord]['tss'][tss] += support  # add tss/tes data to existing locus
					sedict[coord]['tss_tes'][tss][tes] += support
					newcoord = tss if tss < coord[0] else coord[0], tes if tes > coord[1] else coord[1]
					loci[newcoord] = sedict.pop(coord)
				else:
					newlocus = sedict.pop(coord)
					oldlocus = list(loci.keys())[0]  # previous locus the tss/tes overlapped with
					if tss in newlocus['tss']:
						newlocus['tss'][tss] += support
						if tes in newlocus['tss_tes'][tss]:
							newlocus['tss_tes'][tss][tes] += support
						else:
							newlocus['tss_tes'][tss][tes] = loci[oldlocus]['tss_tes'][tss][tes] + support
					loci[oldlocus]['tss_tes'].update(newlocus['tss_tes'])  # combine previous and current loci
					loci[oldlocus]['tss'].update(newlocus['tss'])
					newcoord = oldlocus[0] if oldlocus[0] < coord[0] else coord[0],\
						oldlocus[1] if oldlocus[1] > coord[1] else coord[1]  # combined locus name
					loci[newcoord] = loci.pop(oldlocus)  # add all info under combined locus name
			else:
				break
		sedict.update(loci)
		sedictkeys = sorted(sedict.keys(), key=lambda x: x[0])
	else:  # define new locus
		locus = (tss,tes)
		sedict[locus] = {}
		sedict[locus]['tss'] = {}
		sedict[locus]['tss'][tss] = support
		sedict[locus]['tss_tes'] = {}
		sedict[locus]['tss_tes'][tss] = {}
		sedict[locus]['tss_tes'][tss][tes] = support
		sedict[locus]['line'] = line
		sedict[locus]['bounds'] = (tss, tes)
		sedictkeys = interval_insert((tss,tes),sedictkeys)
	return sedict, sedictkeys

def run_add_se(chrom):
	sedict, sekeys = {}, {}
	sedict[chrom] = {}
	sekeys[chrom] = []  # sorted keys for the singleexon dictionary
	for se in all_se_by_chrom[chrom]:
		sedict[chrom], sekeys[chrom] = add_se(sedict[chrom], se[0], se[1], \
			all_se_by_chrom[chrom][se]['line'], sekeys[chrom], all_se_by_chrom[chrom][se]['count'])
	return sedict

def iterative_add_se(sedict, chrom, group, se):
	tss, tes = se[0], se[1]
	support = all_se_by_chrom[chrom][se]['count']
	if group not in sedict[chrom]:
		sedict[chrom][group] = {}
		sedict[chrom][group]['tss'] = {}
		sedict[chrom][group]['tss_tes'] = {}
		sedict[chrom][group]['bounds'] = [float(tss), float(tes), support]
		sedict[chrom][group]['line'] = all_se_by_chrom[chrom][se]['line']
		sedict[chrom][group]['strand'] = {}
	else:  # calculate mean tss and tes of this group
		if tss > sedict[chrom][group]['bounds'][1] or tes > sedict[chrom][group]['bounds'][1] + window:
			return sedict, False
		sedict[chrom][group]['bounds'][2] += support
		sedict[chrom][group]['bounds'][0] = sedict[chrom][group]['bounds'][0] + \
			support*(tss - sedict[chrom][group]['bounds'][0])/sedict[chrom][group]['bounds'][2]
		sedict[chrom][group]['bounds'][1] = sedict[chrom][group]['bounds'][1] + \
			support*(tes - sedict[chrom][group]['bounds'][1])/sedict[chrom][group]['bounds'][2]

	strand = all_se_by_chrom[chrom][se]['line'][5] if bed else all_se_by_chrom[chrom][se]['line'][8]
	if strand not in sedict[chrom][group]['strand']:
		sedict[chrom][group]['strand'][strand] = 0
	sedict[chrom][group]['strand'][strand] += 1 

	if tss not in sedict[chrom][group]['tss']:
		sedict[chrom][group]['tss'][tss] = 0
		sedict[chrom][group]['tss_tes'][tss] = {}
	sedict[chrom][group]['tss'][tss] += support

	if tes not in sedict[chrom][group]['tss_tes'][tss]:
		sedict[chrom][group]['tss_tes'][tss][tes] = 0
	sedict[chrom][group]['tss_tes'][tss][tes] += support
	return sedict, True

def run_iterative_add_se(chrom):  # add single exon genes iteratively, assumes that reads are sorted
	group = 0  # group number
	sedict = {}
	sedict[chrom] = {}
	se_ordered = all_se_by_chrom[chrom].keys()
	se_ordered = sorted(se_ordered, key = lambda x: x[1])
	se_ordered = sorted(se_ordered, key = lambda x: x[0])
	sedict, added = iterative_add_se(sedict, chrom, group, se_ordered[0])
	for se in se_ordered[1:]:
		overlapped_loci = []
		overlapped_intervals = []
		things = []
		for g in reversed(range(max(0, group-6), group+1)):
			isoverlap, coverage = overlap(se, sedict[chrom][g]['bounds'])
			if isoverlap:
				overlapped_loci += [(g, coverage)]
			else:
				break

		overlapped_loci = sorted(overlapped_loci, key = lambda x: x[1], reverse=True)
		added = False
		for loci in overlapped_loci:
			g = loci[0]
			sedict, added = iterative_add_se(sedict, chrom, g, se)
			if added:
				break
		if not added:
			group += 1
			sedict, added = iterative_add_se(sedict, chrom, group, se)
	return sedict

def find_best_tss(sites, finding_tss, remove_used):
	nearby = dict.fromkeys(sites, 0)  # key: site, value: number of supporting reads in window
	wnearby = dict.fromkeys(sites, 0)  # key: site, value: weighted number of supporting reads in window
	bestsite = [0, 0, 0, 0]  # TSS position, wnearby value, nearby value, number of supporting reads
	if finding_tss:
		orderedsites = sorted(list(sites.keys()))  # this prioritizes the longer in case of ties
	else:
		orderedsites = sorted(list(sites.keys()), reverse=True)
	for s in orderedsites:  # calculate the auxiliary support this site within window
		for s_ in sites:
			if s == s_:
				nearby[s] += sites[s]
				wnearby[s] += sites[s]
			elif abs(s - s_) < window:
				wnearby[s] += (window - abs(s - s_))/float(window * sites[s_])  # downweighted by distance from this site
				nearby[s] += sites[s_]
		if wnearby[s] > bestsite[1]:  # update bestsite if this site has more supporting reads
			bestsite = [s, wnearby[s], nearby[s], sites[s]]

	if remove_used:
		for s in list(sites.keys()):  # remove reads supporting the bestsite to find alternative TSSs
			if abs(s - bestsite[0]) <= window:
				sites.pop(s)
	return sites, bestsite

def find_tsss(sites, total, finding_tss=True, max_results=2, chrom='', junccoord='', remove_used=True):
	""" Finds the best TSSs within Sites. If find_tss is False, some 
	assumptions are changed to search specifically for TESs. I also assume that the correct 
	splice site will be the more represented, so measures to filter out degraded reads are 
	recommended. """
	remaining = float(sum(list(sites.values())))  # number isoforms with these junctions
	found_tss = []  # TSSs found
	used_annotated = set()
	avg = remaining/ (len(sites))
	while ((minsupport < 1 and remaining/total > minsupport) or remaining >= minsupport) and \
			len(found_tss) < max_results:  
		sites, bestsite = find_best_tss(sites, finding_tss, remove_used)
		newremaining = sum(list(sites.values()))
		used = remaining - newremaining
		remaining = newremaining
		if len(found_tss) >= 1 and (args.n == 'best_only' or \
			(minsupport < 1 and (used/total) < minsupport or bestsite[3] < 4 or used < avg)):
		# second+ end site called stringently
			break
		closest_annotated = 1e15  # just a large number
		if annot_tss and chrom in annot_tss:  # args.f supplied
			for t in range(bestsite[0]-window, bestsite[0]+window):
				if (finding_tss and t in annot_tss[chrom] or not finding_tss and t in annot_tes[chrom]) \
					and abs(t - bestsite[0]) < abs(closest_annotated - bestsite[0]) and \
					(junccoord and (finding_tss and t < junccoord[0] or not finding_tss and t > junccoord[1])):
					closest_annotated = t
		if closest_annotated < 1e15:
			if len(found_tss) >= max_results:
				break
			found_tss += [[closest_annotated] + bestsite[1:4]]
		else:
			found_tss += [bestsite]
	return found_tss

def find_best_sites(sites_tss_all, sites_tes_all, junccoord, chrom='', max_results=max_results):
	""" sites_tss_all = {tss: count}
	sites_tes_all = {tss: {tes: count}}
	specific_tes = {tes: count} for a specific set of tss within given window
	junccoord is coordinate of the isoform's first splice site and last splice site"""
	total = float(sum(list(sites_tss_all.values())))  # total number of reads for this splice junction chain
	# if junccoord in [(36178823, 36180248), (36179025, 36180248)]:
	# 	print('--', junccoord)
	# 	print(sites_tss_all)
	found_tss = find_tsss(sites_tss_all, total, finding_tss=True, max_results=max_results, \
		chrom=chrom, junccoord=junccoord)
	# if junccoord in [(36178823, 36180248), (36179025, 36180248)]:
	# 	print(found_tss)
	if not found_tss:
		return ''

	ends = []
	for tss in found_tss:
		specific_tes = {}  # the specific end sites associated with this tss
		for tss_ in sites_tes_all:
			if abs(tss_ - tss[0]) <= window:
				for tes in sites_tes_all[tss_]:
					if tes not in specific_tes:
						specific_tes[tes] = 0
					specific_tes[tes] += sites_tes_all[tss_][tes]

		found_tes = find_tsss(specific_tes, total, finding_tss=False, max_results=max_results, chrom=chrom, \
			junccoord=junccoord)
		for tes in found_tes:
			ends += [[tss[0], tes[0], max(tss[3], tes[3]), tss[3], tes[3]]]
	return ends

def run_se_collapse(chrom):
	senames = {}
	towrite = []
	used_ends = {}
	all_se_starts = {}
	all_se_ends = {}
	for locus in singleexon[chrom]:
		line = list(singleexon[chrom][locus]['line'])
		locus_info = singleexon[chrom][locus]
		ends = find_best_sites(locus_info['tss'], locus_info['tss_tes'], \
			locus_info['bounds'], chrom, max_results=1)
		name = line[3] if bed else line[9]
		if name not in senames:
			senames[name] = 0
		strand = sorted(singleexon[chrom][locus]['strand'].items(),key=lambda x: x[1])[-1][1]
		for tss, tes, support, tsscount, tescount in ends:
			if (tss, tes) not in used_ends:
				used_ends[(tss, tes)] = len(towrite)
			else:
				if not args.c:
					towrite[used_ends[(tss, tes)]][-1] += support
				continue
			if tss not in all_se_starts:
				all_se_starts[tss] = 0
			if tes not in all_se_ends:
				all_se_ends[tes] = 0
			all_se_starts[tss] += support
			all_se_ends[tes] += support

			i = senames[name]
			if not bed:
				edited_line = edit_line(line, tss, tes, tes-tss)
				edited_line[5] = strand
			else:
				edited_line = edit_line_bed12(line, tss, tes, tes-tss)
				edited_line[8] = strand
			if not args.c:
				edited_line += [support]
			if i >= 1:  # to avoid redundant names for isoforms from the same general locus
				if '_' in name:
					newname = name[:name.rfind('_')]+'-'+str(i)+name[name.rfind('_'):]
				else:
					newname = name+'-'+str(i)
				if not bed:
					edited_line[9] = newname
				else:
					edited_line[3] = newname
			senames[name] += 1
			towrite += [edited_line]

	if args.i:
		return towrite
	new_towrite = []
	used_ends = {}
	for line in towrite:
		line = list(line)
		if bed:
			tss, tes = int(line[1]), int(line[2])
		else:
			tss, tes = int(line[15]), int(line[16])
		tss_support = all_se_starts[tss]  # current tss
		for t in reversed(range(tss-window, tss+window)):
			if t in all_se_starts:
				t_support = all_se_starts[t]  # comparison tss
				if t_support > tss_support:
					tss, tss_support = t, t_support
		tes_support = all_se_ends[tes]  # current tes
		for t in range(tes-window, tes+window):
			if t in all_se_ends:
				t_support = all_se_ends[t]  # comparison tes
				if t_support > tes_support:
					tes, tes_support = t, t_support
		if (tss, tes) not in used_ends:
			used_ends[(tss, tes)] = len(new_towrite)
		else:
			if not args.c:
				new_towrite[used_ends[(tss, tes)]][-1] += line[-1]  # support
			continue
		if bed:
			newline = edit_line_bed12(line, tss, tes)
		else:
			newline = edit_line(line, tss, tes)
		new_towrite += [newline]
	return new_towrite

	return towrite

def run_find_best_sites(chrom):
	allends = {}  # counts of all TSSs/TESs by chromosome
	allends[chrom] = {}
	allends[chrom]['tss'] = {}
	allends[chrom]['tes'] = {}
	towrite = {}  # isoforms to be written
	towrite[chrom] = {}
	for jset in isoforms[chrom]:  # unique splice junction chain
		towrite[chrom][jset] = []
		line = list(isoforms[chrom][jset]['line'])
		junccoord = isoforms[chrom][jset]['junccoord']
		jset_ends = find_best_sites(isoforms[chrom][jset]['tss'], \
			isoforms[chrom][jset]['tss_tes'], junccoord, chrom)

		if args.i and args.n == 'longest':
			tss_sorted = sorted(jset_ends, key=lambda x: x[0])[0]  # smallest left coord
			tes_sorted = sorted(jset_ends, key=lambda x: x[1])[-1]  # largest right coord 
			if not args.c:
				line += [(tss_sorted[3]+tes_sorted[3])/2.]
			tss, tes = tss_sorted[0], tes_sorted[1]
			if not bed:
				towrite[chrom][jset] += [edit_line(line, tss, tes)]
			else:
				towrite[chrom][jset] += [edit_line_bed12(line, tss, tes)]
			continue
		if args.i and args.n == 'best_only':
			tss, tes = jset_ends[0][:2]
			if not args.c:
				line += [jset_ends[0][3]]
			if not bed:
				towrite[chrom][jset] += [edit_line(line, tss, tes)]
			else:
				towrite[chrom][jset] += [edit_line_bed12(line, tss, tes)]
			continue  # 1 isoform per unique set of junctions

		i = 0
		name = line[9] if not bed else line[3]
		for tss, tes, support, tsscount, tescount in jset_ends:
			if not bed:
				edited_line = edit_line(line, tss, tes)
			else:
				edited_line = edit_line_bed12(line, tss, tes)
			if not args.c:
				edited_line += [support]
			if i >= 1:  # to avoid redundant names for isoforms with the same junctions
				if '_' in name:
					newname = name[:name.rfind('_')]+'-'+str(i)+name[name.rfind('_'):]
				else:
					newname = name+'-'+str(i)
				if not bed:
					edited_line[9] = newname
				else:
					edited_line[3] = newname
			if args.i:
				towrite[chrom][jset] += [edited_line]
			else:  # all isoforms will go through another pass to homogenize ends
				towrite[chrom][jset] += [edited_line + [junccoord]]
				if tss not in allends[chrom]['tss']:
					allends[chrom]['tss'][tss] = 0
				allends[chrom]['tss'][tss] += tsscount
				if tes not in allends[chrom]['tes']:
					allends[chrom]['tes'][tes] = 0
				allends[chrom]['tes'][tes] = tescount
			if args.n != 'longest':
				i += 1
	if args.i:
		return towrite

	new_towrite = {}  # another pass through all isoforms, moving tss/tes to be more uniform within a gene
	new_towrite[chrom] = {}
	for jset in towrite[chrom]:
		new_towrite[chrom][jset] = []
		jsetends = set()
		endpair = [1e15, 0, 0]  # only applies if no_redundant was specified
		for line in towrite[chrom][jset]:  # adjust isoform TSS/TES using allends dict
			line = list(line)
			junccoord = line[-1]
			if bed:
				tss, tes = int(line[1]), int(line[2])
			else:
				tss, tes = int(line[15]), int(line[16])
			tss_support = allends[chrom]['tss'][tss]  # current tss
			for t in reversed(range(tss-window, tss+window)):
				if t in allends[chrom]['tss']:
					# if 'e68436c4' in line[3]:
					# 	print(t, allends[chrom]['tss'][t], tss_support, tss)
					t_support = allends[chrom]['tss'][t]  # comparison tss
					if (t_support > tss_support and t < junccoord[0]) or \
					(t_support == tss_support and t < tss):
						tss, support = t, t_support

			tes_support = allends[chrom]['tes'][tes]  # current tes
			for t in range(tes-window, tes+window):
				if t in allends[chrom]['tes']:
					t_support = allends[chrom]['tes'][t]  # comparison tes
					if (t_support > tes_support and t > junccoord[1]) or \
					(t_support == tes_support and t > tes):
						tes, tes_support = t, t_support

			if (tss,tes) not in jsetends:
				jsetends.add((tss,tes))
				if args.n == 'best_only':
					endpair = [tss, tes, support] if endpair[2] < tss_support else endpair
				elif args.n == 'longest':
					endpair[0] = tss if tss < endpair[0] else endpair[0]
					endpair[1] = tes if tes > endpair[1] else endpair[1]
				else:
					if bed:
						newline = edit_line_bed12(line, tss, tes)
					else:
						newline = edit_line(line, tss, tes)
					new_towrite[chrom][jset] += [newline[:-1]]
		if endpair[0] != 1e15:  # write best_only or longest option isoforms
			if bed:
				newline = edit_line_bed12(line, endpair[0], endpair[1])
			else:
				newline = edit_line(line, endpair[0], endpair[1])
			new_towrite[chrom][jset] += [newline[:-1]]
	return new_towrite

def edit_line(line, tss, tes, blocksize=''):
	line = list(line)
	if blocksize:  # single exon transcript
		line[18] = str(blocksize) + ','
		line[20] = str(tss) + ','
		line[15] = tss
		line[16] = tes
		return line
	bsizes = [int(x) for x in line[18].split(',')[:-1]]
	bstarts = [int(x) for x in line[20].split(',')[:-1]]
	tstart = int(line[15])  # current chrom start
	tend = int(line[16])
	bsizes[0] += tstart - tss
	bsizes[-1] += tes - tend
	bstarts[0] = tss
	line[15] = tss
	line[16] = tes
	line[18] = ','.join([str(x) for x in bsizes])+','
	line[20] = ','.join([str(x) for x in bstarts])+','
	return line

def edit_line_bed12(line, tss, tes, blocksize=''):
	line = list(line)
	if blocksize:
		line[10] = str(blocksize) + ','
		line[11] = '0,'
		line[1] = line[6] = tss
		line[2] = line[7] = tes
		return line
	bsizes = [int(x) for x in line[10].split(',')[:-1]]
	bstarts = [int(x) for x in line[11].split(',')[:-1]]
	tstart = int(line[1])  # current chrom start
	tend = int(line[2])
	bsizes[0] += tstart - tss
	bsizes[-1] += tes - tend
	bstarts = [x + int(line[1]) for x in bstarts]
	bstarts[0] = tss
	line[1] = line[6] = tss
	line[2] = line[7] = tes
	line[10] = ','.join([str(x) for x in bsizes])+','
	line[11] = ','.join([str(x - int(line[1])) for x in bstarts])+','
	return line


annot_tss = {}  # annotated left terminal sites per chromosome from GTF
annot_tes = {}  # annotated right terminal sites
if args.f:
	for line in gtf:
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end = line[0], line[2], int(line[3]) - 1, int(line[4])
		if ty == 'transcript':
			if chrom not in annot_tss:
				annot_tss[chrom] = set()
				annot_tes[chrom] = set()
			annot_tss[chrom].add(start)
			annot_tes[chrom].add(end)
	if not args.quiet:
		sys.stderr.write('Annotated ends extracted from GTF\n')

isoforms = {}  # spliced isoforms
all_se_by_chrom = {}  # single-exon reads grouped by exact start and end site
nosplice_chroms = set(args.nosplice.split(','))
for line in psl:
	line = tuple(line.rstrip().split('\t'))
	if bed:
		chrom, tss, tes = line[0], int(line[1]), int(line[2])
		junctions, junccoord = get_junctions_bed12(line)
	else:
		chrom, tss, tes = line[13], int(line[15]), int(line[16])
		junctions, junccoord =  get_junctions(line)

	if not junctions:  # single-exon isoforms
		if chrom not in all_se_by_chrom:
			all_se_by_chrom[chrom] = {}
		if (tss, tes) not in all_se_by_chrom[chrom]:
			all_se_by_chrom[chrom][(tss, tes)] = {}
			all_se_by_chrom[chrom][(tss, tes)]['count'] = 0
			all_se_by_chrom[chrom][(tss, tes)]['line'] = line
		all_se_by_chrom[chrom][(tss, tes)]['count'] += 1
		continue
	elif chrom in nosplice_chroms:
		continue

	junctions = str(sorted(list(junctions)))  # splice junction chain
	# if '(36179025, 36179262)' in junctions:
	# 	print(tss, junccoord, line[3])

	if chrom not in isoforms:
		isoforms[chrom] = {}
	if junctions not in isoforms[chrom]:
		isoforms[chrom][junctions] = {}
		isoforms[chrom][junctions]['tss'] = {}  # ends stored for each unique sjc
		isoforms[chrom][junctions]['tss_tes'] = {}
		isoforms[chrom][junctions]['line'] = line
		isoforms[chrom][junctions]['junccoord'] = junccoord

	if tss not in isoforms[chrom][junctions]['tss']:
		isoforms[chrom][junctions]['tss'][tss] = 0  # tss usage count
		isoforms[chrom][junctions]['tss_tes'][tss] = {}  # TESs for this TSS
	isoforms[chrom][junctions]['tss'][tss] += 1

	if tes not in isoforms[chrom][junctions]['tss_tes'][tss]:
		isoforms[chrom][junctions]['tss_tes'][tss][tes] = 0
	isoforms[chrom][junctions]['tss_tes'][tss][tes] += 1

if not args.quiet:
	sys.stderr.write('Read data extracted\n')

chrom_names = []  # sorted by descending total number of unique single-exon isoforms
for chrom in all_se_by_chrom:
	chrom_names += [(chrom, len(all_se_by_chrom[chrom]))]
chrom_names = [chrom for chrom,num in sorted(chrom_names, key=lambda x: x[1], reverse=True)]
if __name__ == '__main__':
	p = Pool(args.t)
	res = p.map(run_iterative_add_se, chrom_names)
	p.terminate()
all_se_by_chrom = None

singleexon = {}  # single-exon isoforms
for r in res:  # combine results
	singleexon.update(r)

if not args.quiet:
	sys.stderr.write('Single-exon genes grouped, collapsing\n')

with open(args.o, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)

	if __name__ == '__main__':
		p = Pool(args.t)
		res = p.map(run_se_collapse, chrom_names)
		p.terminate()
	singleexon = None

	for r in res:
		for edited_line in r:
			writer.writerow(edited_line)

	if __name__ == '__main__':
		p = Pool(args.t)
		res = p.map(run_find_best_sites, list(isoforms.keys()))
		p.terminate()
	isoforms = None

	for towrite in res:
		for chrom in towrite:
			for jset in towrite[chrom]:
				for iso in towrite[chrom][jset]:
					writer.writerow(iso)
