import sys, csv, argparse

parser = argparse.ArgumentParser(description='collapse parse options', \
			usage='python collapse_isoforms_precise.py -q <query.psl>/<query.bed> [options]')
required = parser.add_argument_group('required named arguments')
required.add_argument('-q', '--query', type=str, default='', required=True, \
	action='store', dest='q', help='BED12 or PSL file of aligned/corrected reads. PSL should end in .psl')
parser.add_argument('-o', '--output', type=str, \
	action='store', dest='o', default='', help='specify output file, psl or bed12')
parser.add_argument('-w', '--window', default=20, type=int, \
	action='store', dest='w', help='window size for grouping TSS/TES (20)')
parser.add_argument('-s', '--support', default=1, type=int, \
	action='store', dest='s', help='minimum number of supporting reads for an isoform (1)')
parser.add_argument('-f', '--gtf', default='', type=str, \
	action='store', dest='f', help='GTF annotation file for selecting annotated TSS/TES')
parser.add_argument('-m', '--max_results', default=2, type=int, \
	action='store', dest='m', help='maximum number of novel TSS or TES picked per isoform (2)')
parser.add_argument('-n', '--no_redundant', default='none', \
	action='store', dest='n', help='Report options include: \
	none: best TSSs/TESs chosen for each unique set of splice junctions; \
	longest: TSS/TES chosen to maximize length; \
	best_only: best TSS/TES used in conjunction chosen; \
	longest/best_only override max_results argument immediately before output \
	resulting in one isoform per unique set of splice junctions (default: not specified)')
parser.add_argument('-i', '--isoformtss', default=False, \
	action='store_true', dest='i', help='when specified, TSS/TES for each isoform will be determined \
	from supporting reads for individual isoforms, not from the gene level (default: not specified)')
parser.add_argument('--quiet', default=False, \
	action='store_true', dest='quiet', help='suppress output to stderr')
args = parser.parse_args()

try:
	max_results, window, minsupport, psl = args.m, args.w, args.s, open(args.q)
	if args.f:
		gtf = open(args.f)
except:
	sys.stderr.write('Make sure all files (query, GTF) have valid paths and can be opened\n')
	sys.exit()
bed = args.q[-3:] != 'psl'
pslout = True
if args.o:
	if args.o[-3:] != 'psl':
		pslout = False
else:  # default output name
	args.o = args.q[:-3]+'collapsed.bed12' if bed else args.q[:-3]+'collapsed.psl'

def get_junctions(line):
	junctions = set()
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions, starts[0], starts[-1]+sizes[-1], (starts[0]+sizes[0], starts[-1])

def get_junctions_bed12(line):
	junctions = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions, (starts[0]+sizes[0], starts[-1])

def get_start_end(line):
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	return starts[0], starts[-1]+sizes[-1]

def find_best_tss(sites, total, finding_tss):
	nearby = dict.fromkeys(sites, 0)  # key site, value number of supporting reads window
	wnearby = dict.fromkeys(sites, 0)  # key site, value weighted number of supporting reads window
	bestsite = (0, 0, 0)  # TSS position, number of supporting reads, weighted supporting reads
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
			bestsite = (s, wnearby[s], nearby[s], sites[s])

	for s in list(sites.keys()):  # remove reads supporting the bestsite to find alternative TSSs
		if abs(s - bestsite[0]) <= window:
			sites.pop(s)
	return sites, bestsite

def find_tsss(sites, finding_tss=True, max_results=2, chrom='', junccoord=''):
	""" Finds the best TSSs within Sites. If find_tss is False, some 
	assumptions are changed to search specifically for TESs. I also assume that the correct 
	splice site will be the more represented, so measures to filter out degraded reads are 
	recommended. """
	remaining = total = float(sum(list(sites.values())))
	found_tss = []  # TSSs found will be ordered by descending importance
	novel_tss = 0  # number of novel TSS found
	while ((minsupport < 1 and minsupport > 0.5 and remaining/total > minsupport) or \
	remaining >= minsupport):  # stop cases: bestsite encompasses >50% of all sites OR 
	# no other site would pass the supporting read threshold AND no more results are to be reported
		sites, bestsite = find_best_tss(sites, total, finding_tss)
		remaining = sum(list(sites.values()))
		if bestsite == (0, 0, 0, 0):
			break
		closest_annotated = 1e15  # just a large number
		if annotends:  # args.f supplied
			for t in range(bestsite[0]-window, bestsite[0]+window):
				if t in annotends[chrom] and t - bestsite[0] < closest_annotated - bestsite[0]:
					closest_annotated = t
			if finding_tss and closest_annotated >= junccoord[0] or \
			not finding_tss and closest_annotated <= junccoord[1]:
				closest_annotated = 1e15  # annotated site was invalid
		if closest_annotated < 1e15 and closest_annotated not in found_tss:
			found_tss += [(closest_annotated, bestsite[1], bestsite[2], bestsite[3])]
		else:
			if novel_tss >= max_results:  # limit to max_results number of novel end sites
				continue
			if len(found_tss) > max_results:  # if annotated end site(s) have been identified, stop searching
				break
			found_tss += [bestsite]
			novel_tss += 1

	return found_tss  # port over the code that prioritizes annotated start sites

def find_best_sites(sites_tss_all, sites_tes_all, junccoord, chrom=''):
	""" sites_tss_all = {tss: count}
	sites_tes_all = {tss: {tes: count}}
	specific_tes = {tes: count} for a specific set of tss within given window
	junccoord is coordinate of the isoform's first splice site and last splice site"""
	total = float(sum(list(sites_tss_all.values())))  # number isoforms with these junctions
	found_tss = find_tsss(sites_tss_all, finding_tss=True, max_results=max_results, \
		chrom=chrom, junccoord=junccoord)
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

		found_tes = find_tsss(specific_tes, finding_tss=False, max_results=max_results, chrom=chrom, \
			junccoord=junccoord)
		for tes in found_tes:
			ends += [(tss[0], tes[0], tes[2], tss[3], tes[3])]
	return ends

def single_exon_pairs(sedict):
	""" """
	found_tss = find_tsss(sedict, finding_tss=True, max_results=1)
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
		found_tes = find_tsss(specific_tes, finding_tss=False)
		for tes in found_tes:
			ends += [(tss[0], tes[0], tes[2])]

	tss_nearby = dict.fromkeys(sedict['tss'], 0)
	tes_nearby = dict.fromkeys(sedict['tes'], 0)
	for tss in sedict['tss']:
		for tss_ in sedict['tss']:
			tss_nearby[tss] += sedict['tss'][tss_]
	for tes in sedict['tes']:
		for tes_ in sedict['tes']:
			tes_nearby[tes] += sedict['tes'][tes_]
	alltss = sorted(tss_nearby.keys())
	alltes = sorted(tes_nearby.keys())
	tssi, tesi = 0, 0
	tss_stack = []
	pairs = []
	while tssi < len(alltss) and tesi < len(alltes):
		if tes_nearby[alltes[tesi]] < 3:
			tssi += 1
			continue
		if alltes[tesi] < alltss[tssi]:
			pairs += []
			tss_stack = [alltss[tssi]]
		else:
			tss_stack += [alltss[tssi]]
	return pairs

def edit_line(line, tss, tes, blocksize=''):
	if blocksize:
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
	if blocksize:
		line[10] = str(blocksize) + ','
		line[11] = '0,'
		line[1] = line[6] = tss
		line[2] = line[7] = tes
		return line

	bsizes = [int(x) for x in line[10].split(',')[:-1]]

	oldbsizes = bsizes

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

annotends = {}  # all annotated TSS/TES sites per chromosome
if args.f:
	for line in gtf:
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end = line[0], line[2], int(line[3]), int(line[4])
		if ty == 'transcript':
			if chrom not in annotends:
				annotends[chrom] = set()
			annotends[chrom].add(start)
			annotends[chrom].add(end)

isoforms, singleexon = {}, {}  # spliced isoforms and single exon isoforms
for line in psl:
	line = tuple(line.rstrip().split('\t'))
	if bed:
		chrom = line[0]
		tss, tes = int(line[1]), int(line[2])
		junctions, junccoord = get_junctions_bed12(line)
	else:
		chrom = line[13]
		junctions, tss, tes, junccoord = get_junctions(line)

	if not junctions:  # single-exon isoforms 
		if chrom not in singleexon:
			singleexon[chrom] = {}
			singleexon[chrom]['tss'] = {}
			singleexon[chrom]['tss_tes'] = {}
			singleexon[chrom]['line'] = {}
		if tss not in singleexon[chrom]['tss']:
			singleexon[chrom]['tss'][tss] = 0
			singleexon[chrom]['tss_tes'][tss] = {}
		singleexon[chrom]['tss'][tss] += 1
		if tes not in singleexon[chrom]['tss_tes'][tss]:
			singleexon[chrom]['tss_tes'][tss][tes] = 0
		singleexon[chrom]['tss_tes'][tss][tes] += 1
		continue

	junctions = str(sorted(list(junctions)))  # hashable but still unique

	if chrom not in isoforms:
		isoforms[chrom] = {}
	if junctions not in isoforms[chrom]:
		isoforms[chrom][junctions] = {}
		isoforms[chrom][junctions]['tss'] = {}
		isoforms[chrom][junctions]['line'] = line
		isoforms[chrom][junctions]['junccoord'] = junccoord
		isoforms[chrom][junctions]['tss_tes'] = {}
	if tss not in isoforms[chrom][junctions]['tss']:
		isoforms[chrom][junctions]['tss'][tss] = 0  #[0, []] tss usage count, list of ends
		isoforms[chrom][junctions]['tss_tes'][tss] = {}
	isoforms[chrom][junctions]['tss'][tss] += 1
	if tes not in isoforms[chrom][junctions]['tss_tes'][tss]:
		isoforms[chrom][junctions]['tss_tes'][tss][tes] = 0
	isoforms[chrom][junctions]['tss_tes'][tss][tes] += 1

if not args.quiet:
	sys.stderr.write('Read data extracted, collapsing\n')

with open(args.o, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	allends = {}  # counts of TSSs/TESs by chromosome
	writtenisos = {}
	for chrom in isoforms:
		if chrom not in allends:
			allends[chrom] = {}
			allends[chrom]['tss'] = {}
			allends[chrom]['tes'] = {}
			writtenisos[chrom] = {}
		for jset in isoforms[chrom]:
			if jset not in writtenisos[chrom]:
				writtenisos[chrom][jset] = []
			line = isoforms[chrom][jset]['line']
			junccoord = isoforms[chrom][jset]['junccoord']
			ends = find_best_sites(isoforms[chrom][jset]['tss'], isoforms[chrom][jset]['tss_tes'], \
				junccoord, chrom)
			name = line[9] if pslout else line[3]

			if args.i and args.n == 'longest':
				tss = sorted(ends, key=lambda x: x[0])[0][0]
				tes = sorted(ends, key=lambda x: x[1])[-1][1]
				if pslout:
					edited_line = edit_line(list(line), tss, tes)
				else:
					edited_line = edit_line_bed12(list(line), tss, tes)
				writer.writerow(edited_line)
				continue  # 1 isoform per unique set of junctions
			if args.i and args.n == 'best_only':
				tss, tes = sorted(ends, key=lambda x: x[2])[-1][:2]
				if pslout:
					edited_line = edit_line(list(line), tss, tes)
				else:
					edited_line = edit_line_bed12(list(line), tss, tes)
				writer.writerow(edited_line)
				continue

			i = 0
			for tss, tes, support, tsscount, tescount in ends:
				if pslout:
					edited_line = edit_line(list(line), tss, tes)
				else:
					edited_line = edit_line_bed12(list(line), tss, tes)
				if i >= 1:
					if pslout:  # to avoid redundant names for isoforms with the same junctions
						edited_line[9] = name+'-'+str(i)	
					else:
						edited_line[3] = name+'-'+str(i)
				if args.i:
					writer.writerow(edited_line)
				else:  # all isoforms will go through a second pass to homogenize novel ends
					if args.n == 'best_only':
						writtenisos[chrom][jset] += [edited_line + [support, junccoord]]
					else:
						writtenisos[chrom][jset] += [edited_line + [junccoord]]
					allends[chrom]['tss'][tss] = tsscount
					allends[chrom]['tes'][tes] = tescount
				i += 1
 
	# senames = {}
	# for chrom in singleexon:
	# 	pairs = single_exon_pairs(singleexon[chrom])
	# 	i = 0
	# 	for p in pairs:
	# 		tss, tes, support = p
	# 		if pslout:
	# 			line = edit_line(list(line), tss, tes)
	# 		else:
	# 			line = edit_line_bed12(list(line), tss, tes)
	# 		writer.writerow(line + [support])

	if args.i:
		sys.exit()

	isoforms = 0
	for chrom in writtenisos:
		for jset in writtenisos[chrom]:
			jsetends = set()
			endpair = [1e15, 0, 0]  # only applies if no_redundant was specified
			for line in writtenisos[chrom][jset]:  # adjust isoform TSS/TES using allends dict
				if bed:
					tss, tes = int(line[1]), int(line[2])
				else:
					tss, tes = get_start_end(line)
				junccoord = line[-1]
				for t in range(tss-window, tss+window):
					if t in allends[chrom]['tss']:
						t_support = allends[chrom]['tss'][t]  # comparison tss
						tss_support = allends[chrom]['tss'][tss]  # current best tss
						if (t_support > tss_support and t < junccoord[0]) or \
						(t_support == tss_support and t < tss):
							tss = t
							support = line[-2]

				for t in range(tes-window, tes+window):
					if t in allends[chrom]['tes']:
						t_support = allends[chrom]['tes'][t]  # comparison tes
						tes_support = allends[chrom]['tes'][tes]  # current best tes
						if (t_support > tes_support and t > junccoord[1]) or \
						(t_support == tes_support and t > tes):
							tes = t
							support = line[-2]  # only used for best_only option

				if (tss,tes) not in jsetends:
					jsetends.add((tss,tes))
					if args.n == 'best_only':
						endpair = [tss, tes, support] if endpair[2] < support else endpair
					elif args.n == 'longest':
						endpair[0] = tss if tss < endpair[0] else endpair[0]
						endpair[1] = tes if tes > endpair[1] else endpair[1]
					else:
						if bed:
							newline = edit_line_bed12(list(line), tss, tes)
						else:
							newline = edit_line(list(line), tss, tes)
						writer.writerow(newline[:-1])

			if endpair[0] != 1e15:
				if bed:
					newline = edit_line_bed12(list(line), endpair[0], endpair[1])
				else:
					newline = edit_line(list(line), endpair[0], endpair[1])
				if args.n == 'best_only':
					writer.writerow(newline[:-2])
				else:
					writer.writerow(newline[:-1])		





