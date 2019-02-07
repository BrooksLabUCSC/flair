import sys, csv, argparse, math
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='collapse parse options', \
			usage='python collapse_isoforms_precise.py -q <query.psl>/<query.bed> [options]')
required = parser.add_argument_group('required named arguments')
required.add_argument('-q', '--query', type=str, default='', required=True, action='store', \
	dest='q', help='BED12 or PSL file of aligned/corrected reads. PSL should end in .psl')
parser.add_argument('-o', '--output', type=str, action='store', \
	dest='o', default='', help='specify output file, should agree with query file type')
parser.add_argument('-w', '--window', default=100, type=int, \
	action='store', dest='w', help='window size for grouping TSS/TES (20)')
parser.add_argument('-s', '--support', default=0.25, type=float, action='store', \
	dest='s', help='minimum proportion(s<1) or number of supporting reads(s>=1) for an isoform (0.25)')
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
	sys.exit()
bed = args.q[-3:].lower() != 'psl'
pslout = True
if args.o:
	if args.o[-3:].lower() != 'psl':
		pslout = False
else:  # default output name
	args.o = args.q[:-3]+'collapsed.bed12' if bed else args.q[:-3]+'collapsed.psl'

def get_junctions(line):
	junctions = set()
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	if len(starts) == 1:
		return 0, starts[0], starts[-1]+sizes[-1], 0
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions, starts[0], starts[-1]+sizes[-1], (starts[0]+sizes[0], starts[-1])

def get_junctions_bed12(line):
	junctions = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return 0, 0
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions, (starts[0]+sizes[0], starts[-1])

def get_start_end(line):
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	return starts[0], starts[-1]+sizes[-1]

def overlap(coord0, coord1, tol=0):
	coord0, coord1 = sorted([coord0, coord1], key = lambda x: x[0])
	return (coord0[0] <= coord1[0] and coord1[0] <= coord0[1] - tol)

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

def find_best_tss(sites, finding_tss):
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

def find_tsss(sites, total, finding_tss=True, max_results=2, chrom='', junccoord=''):
	""" Finds the best TSSs within Sites. If find_tss is False, some 
	assumptions are changed to search specifically for TESs. I also assume that the correct 
	splice site will be the more represented, so measures to filter out degraded reads are 
	recommended. """
	remaining = float(sum(list(sites.values())))  # number isoforms with these junctions
	found_tss = []  # TSSs found will be ordered by descending importance
	novel_tss = 0  # number of novel TSS found
	while ((minsupport < 1 and remaining/total > minsupport) or remaining >= minsupport):  
		sites, bestsite = find_best_tss(sites, finding_tss)
		used = remaining - sum(list(sites.values()))
		if minsupport < 1 and (used/(total) < minsupport or used < 3) and len(found_tss) >= 1:
		# found at minimum the best site, and this site did not surpass minsupport percentage of reads
			break
		remaining = sum(list(sites.values()))
		closest_annotated = 1e15  # just a large number
		if annotends:  # args.f supplied
			for t in range(bestsite[0]-window, bestsite[0]+window):
				if t in annotends[chrom] and (t - bestsite[0]) < (closest_annotated - bestsite[0]):
					closest_annotated = t
			if junccoord and (finding_tss and closest_annotated >= junccoord[0] or \
			not finding_tss and closest_annotated <= junccoord[1]):
				closest_annotated = 1e15  # annotated site was invalid
		if closest_annotated < 1e15 and closest_annotated not in found_tss:
			found_tss += [(closest_annotated, bestsite[1], bestsite[2], bestsite[3], bestsite[0])]
		else:
			if novel_tss >= max_results:  # limit to max_results number of novel end sites
				continue
			if len(found_tss) > max_results:  # if annotated end site(s) have been identified, stop searching
				break
			found_tss += [bestsite]
			novel_tss += 1
	return found_tss

def find_best_sites(sites_tss_all, sites_tes_all, junccoord, chrom='', max_results=max_results):
	""" sites_tss_all = {tss: count}
	sites_tes_all = {tss: {tes: count}}
	specific_tes = {tes: count} for a specific set of tss within given window
	junccoord is coordinate of the isoform's first splice site and last splice site"""
	total = float(sum(list(sites_tss_all.values())))  # total number of reads for this splice junction chain
	found_tss = find_tsss(sites_tss_all, total, finding_tss=True, max_results=max_results, \
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

		found_tes = find_tsss(specific_tes, total, finding_tss=False, max_results=max_results, chrom=chrom, \
			junccoord=junccoord)
		for tes in found_tes:
			ends += [(tss[0], tes[0], tes[2], tss[3], tes[3])]
	return ends

def run_se_collapse(chrom):
	senames = {}
	towrite = []
	for locus in singleexon[chrom]:
		line = singleexon[chrom][locus]['line']
		locus_info = singleexon[chrom][locus]
		ends = find_best_sites(locus_info['tss'], locus_info['tss_tes'], \
			locus_info['bounds'], chrom, max_results=max_results)
		name = line[9] if pslout else line[3]
		if name not in senames:
			senames[name] = 0
		for tss, tes, support, tsscount, tescount in ends:
			i = senames[name]
			senames[name] += 1
			if pslout:
				edited_line = edit_line(list(line), tss, tes, tes-tss)
			else:
				edited_line = edit_line_bed12(list(line), tss, tes, tes-tss)
			if i >= 1:  # to avoid redundant names for isoforms from the same general locus
				if pslout:
					edited_line[9] = name+'-'+str(i)	
				else:
					edited_line[3] = name+'-'+str(i)
			towrite += [edited_line]
	return towrite

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

annotends = {}  # all annotated TSS/TES sites per chromosome from GTF
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
		chrom = line[13]
		junctions, tss, tes, junccoord =  get_junctions(line)

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
	res = p.map(run_add_se, chrom_names)
	p.close()
all_se_by_chrom = None
singleexon = {}  # single-exon isoforms
for r in res:
	singleexon.update(r)

if not args.quiet:
	sys.stderr.write('Single-exon genes grouped, collapsing\n')

with open(args.o, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')

	if __name__ == '__main__':
		p = Pool(args.t)
		res = p.map(run_se_collapse, chrom_names)
		p.terminate()

	for r in res:
		for edited_line in r:
			writer.writerow(edited_line)
	singleexon = None

	allends = {}  # counts of all TSSs/TESs by chromosome
	towrite = {}  # isoforms to be written
	for chrom in isoforms:
		if chrom not in allends:
			allends[chrom] = {}
			allends[chrom]['tss'] = {}
			allends[chrom]['tes'] = {}
			towrite[chrom] = {}
		for jset in isoforms[chrom]:
			if jset not in towrite[chrom]:
				towrite[chrom][jset] = []
			line = isoforms[chrom][jset]['line']
			junccoord = isoforms[chrom][jset]['junccoord']
			ends = find_best_sites(isoforms[chrom][jset]['tss'], \
				isoforms[chrom][jset]['tss_tes'], junccoord, chrom)

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
			name = line[9] if pslout else line[3]
			for tss, tes, support, tsscount, tescount in ends:
				if pslout:
					edited_line = edit_line(list(line), tss, tes)
				else:
					edited_line = edit_line_bed12(list(line), tss, tes)
				if i >= 1:  # to avoid redundant names for isoforms with the same junctions
					if pslout:
						edited_line[9] = name+'-'+str(i)	
					else:
						edited_line[3] = name+'-'+str(i)
				if args.i:
					writer.writerow(edited_line)
				else:  # all isoforms will go through a second pass to homogenize ends
					if args.n == 'best_only':
						towrite[chrom][jset] += [edited_line + [support, junccoord]]
					else:  # longest
						towrite[chrom][jset] += [edited_line + [junccoord]]
					if tss not in allends[chrom]['tss']:
						allends[chrom]['tss'][tss] = 0
					allends[chrom]['tss'][tss] += tsscount
					if tes not in allends[chrom]['tes']:
						allends[chrom]['tes'][tes] = 0
					allends[chrom]['tes'][tes] = tescount
				i += 1

	if args.i:
		sys.exit()
	isoforms = None

	for tes in allends['chr6']['tes']:
		if tes < (36569936 + window) and tes > (36569936 - window):
			print(tes, allends['chr6']['tes'][tes])

	for chrom in towrite:
		for jset in towrite[chrom]:
			jsetends = set()
			endpair = [1e15, 0, 0]  # only applies if no_redundant was specified
			for line in towrite[chrom][jset]:  # adjust isoform TSS/TES using allends dict
				if bed:
					tss, tes = int(line[1]), int(line[2])
				else:
					tss, tes = get_start_end(line)
				junccoord = line[-1]
				tss_support = allends[chrom]['tss'][tss]  # current tss
				for t in range(tss+window, tss-window, -1):
					if t in allends[chrom]['tss']:
						t_support = allends[chrom]['tss'][t]  # comparison tss
						if (t_support > tss_support and t < junccoord[0]) or \
						(t_support == tss_support and t < tss):
							tss = t
							support = t_support

				tes_support = allends[chrom]['tes'][tes]  # current tes
				oldtes = tes
				if tes == 36570029:
					print(tes_support)
				for t in range(tes-window, tes+window):
					if t in allends[chrom]['tes']:
						t_support = allends[chrom]['tes'][t]  # comparison tes
						if oldtes == 36570029:
							print(t_support, t)
						if (t_support > tes_support and t > junccoord[1]) or \
						(t_support == tes_support and t > tes):
							tes = t
							tes_support = t_support  # only used for best_only option
				if oldtes == 36570029:
					print(tes, tes_support)
				if tes == 36570029:
					print('this weirdo tes got corrected to 36570029', oldtes, line)


				if (tss,tes) not in jsetends:
					jsetends.add((tss,tes))
					if args.n == 'best_only':
						endpair = [tss, tes, support] if endpair[2] < tss_support else endpair
					elif args.n == 'longest':
						endpair[0] = tss if tss < endpair[0] else endpair[0]
						endpair[1] = tes if tes > endpair[1] else endpair[1]
					else:
						if bed:
							newline = edit_line_bed12(list(line), tss, tes)
						else:
							newline = edit_line(list(line), tss, tes)
						writer.writerow(newline[:-1])
			if endpair[0] != 1e15:  # write best_only or longest option isoforms
				if bed:
					newline = edit_line_bed12(list(line), endpair[0], endpair[1])
				else:
					newline = edit_line(list(line), endpair[0], endpair[1])
				if args.n == 'best_only':
					writer.writerow(newline[:-2])
				else:
					writer.writerow(newline[:-1])		

