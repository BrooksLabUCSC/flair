#!/usr/bin/env python3

import sys
import csv
import argparse
import math
import os
from multiprocessing import Pool
from bed import Bed, BedBlock

parser = argparse.ArgumentParser(description='collapse parse options',
			usage='python collapse_isoforms_precise.py -q <query.psl>/<query.bed> [options]')
required = parser.add_argument_group('required named arguments')
required.add_argument('-q', '--query', type=str, required=True, 
	help='BED12 or PSL file of aligned/corrected reads. PSL should end in .psl')
parser.add_argument('-o', '--output', type=str, 
	help='specify output file, should agree with query file type')
parser.add_argument('-w', '--window', default=100, type=int,
	 help='window size for grouping TSS/TES (100)')
parser.add_argument('-s', '--support', default=0.25, type=float, 
	help='minimum proportion(s<1) or number of supporting reads(s>=1) for an isoform (0.1)')
parser.add_argument('-f', '--gtf', type=str,
	 help='GTF annotation file for selecting annotated TSS/TES')
parser.add_argument('-m', '--max_results', default=2, type=int, 
	help='maximum number of novel TSS or TES picked per isoform unless --no_redundant is specified (2)')
parser.add_argument('-t', '--threads', default=2, type=int,
	 help='number of threads to use (2)')
parser.add_argument('-n', '--no_redundant', default='none',
	help='For each unique splice junction chain, report options include: \
	none: multiple best TSSs/TESs chosen for each unique set of splice junctions, see M; \
	longest: TSS/TES chosen to maximize length; \
	best_only: single best TSS/TES used in conjunction chosen; \
	longest/best_only override max_results argument immediately before output \
	resulting in one isoform per unique set of splice junctions (default: none)')
parser.add_argument('-c', '--clean', default=False, action='store_true',
	help='Specify this to not append read support to the end of each entry (default: not specified)')
parser.add_argument('-i', '--isoformtss', default=False, action='store_true',
	help='when specified, TSS/TES for each isoform will be determined from supporting reads \
	for individual isoforms (default: not specified, determined at the gene level)')
parser.add_argument('--nosplice', default='chrM',
	help='Comma separated list of chromosomes that should not have spliced isoforms (default: chrM)')
parser.add_argument('--quiet', default=False, action='store_true',
	help='suppress output to stderr')

args = parser.parse_args()

if args.output:
	if args.output[-3:].lower() != args.query[-3:].lower():
		sys.stderr.write('Make sure input and output file extensions agree\n')
		sys.exit(1)
else:  # default output name
	args.output = args.query[:-3]+'collapsed.bed'

# This renaming of arguments is in preparation of turning this program into a function
queryfile=args.query
max_results=args.max_results
threads=args.threads
window=args.window
isoformtss=args.isoformtss
minsupport=args.support
gtfname=args.gtf
clean=args.clean
nosplice=args.nosplice
no_redundant=args.no_redundant
outputfname=args.output
quiet=args.quiet

#def collapse_isoforms_precise(queryfile, max_results=2, window=100, threads=2, clean=False,
#	minsupport=0.25, gtfname=None, no_redundant='none', nosplice='chrM', isoformtss=False, 
#	outputfname=None, quiet=False):

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

	strand = all_se_by_chrom[chrom][se]['line'][5]
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
	se_ordered = sorted(se_ordered, key=lambda x: x[1])
	se_ordered = sorted(se_ordered, key=lambda x: x[0])
	sedict, added = iterative_add_se(sedict, chrom, group, se_ordered[0])
	for se in se_ordered[1:]:
		overlapped_loci = []
		#overlapped_intervals = [] # unused
		for g in reversed(range(max(0, group-6), group+1)):
			isoverlap, coverage = overlap(se, sedict[chrom][g]['bounds'])
			if isoverlap:
				overlapped_loci += [(g, coverage)]
			else:
				break

		overlapped_loci = sorted(overlapped_loci, key=lambda x: x[1], reverse=True)
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
	#used_annotated = set() # unused
	avg = remaining / (len(sites))
	while ((minsupport < 1 and remaining/total > minsupport) or remaining >= minsupport) and \
			len(found_tss) < max_results:
		sites, bestsite = find_best_tss(sites, finding_tss, remove_used)
		newremaining = sum(list(sites.values()))
		used = remaining - newremaining
		remaining = newremaining
		if len(found_tss) >= 1 and (no_redundant == 'best_only' or
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
	found_tss = find_tsss(sites_tss_all, total, finding_tss=True, max_results=max_results,
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

		found_tes = find_tsss(specific_tes, total, finding_tss=False, max_results=max_results, chrom=chrom,
			junccoord=junccoord)
		for tes in found_tes:
			ends += [[tss[0], tes[0], max(tss[3], tes[3]), tss[3], tes[3]]]
	return ends


def run_se_collapse(chrom):
	'''Collapse single exon isoforms'''
	senames = {}
	towrite = []
	used_ends = {}
	all_se_starts = {}
	all_se_ends = {}
	bedObjs = []
	for locus in singleexon[chrom]:
		line = list(singleexon[chrom][locus]['line'])
		locus_info = singleexon[chrom][locus]
		ends = find_best_sites(locus_info['tss'], locus_info['tss_tes'],
			locus_info['bounds'], chrom, max_results=1)
		name = line[3]
		if name not in senames:
			senames[name] = 0
		strand = sorted(singleexon[chrom][locus]['strand'].items(),key=lambda x: x[1])[-1][1]
		for tss, tes, support, tsscount, tescount in ends:
			if (tss, tes) not in used_ends:
				used_ends[(tss, tes)] = len(towrite)
			else:
				if not clean:
					towrite[used_ends[(tss, tes)]][-1] += support
				continue
			if tss not in all_se_starts:
				all_se_starts[tss] = 0
			if tes not in all_se_ends:
				all_se_ends[tes] = 0
			all_se_starts[tss] += support
			all_se_ends[tes] += support

			i = senames[name]
		#	print('from single exon, sending', tes, tss, tes-tss, file=sys.stderr)
			edited_line = edit_line_bed12(line, tss, tes, tes-tss)
			edited_line[8] = strand
			if not clean:
				edited_line += [support]
			if i >= 1:  # to avoid redundant names for isoforms from the same general locus
				if '_' in name:
					newname = name[:name.rfind('_')]+'-'+str(i)+name[name.rfind('_'):]
				else:
					newname = name+'-'+str(i)
				edited_line[3] = newname
			senames[name] += 1
			towrite += [edited_line]

	if isoformtss:
		return towrite
	new_towrite = []
	used_ends = {}
	for line in towrite:
		line = list(line)
		tss, tes = int(line[1]), int(line[2])
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
			if not clean:
				new_towrite[used_ends[(tss, tes)]][-1] += line[-1]  # support
			continue
	#	print('from towrite, sending', tss, tes, file=sys.stderr)
		newline = edit_line_bed12(line, tss, tes)
		new_towrite += [newline]
	return new_towrite


def run_find_best_sites(chrom):

	allends = {}  # counts of all TSSs/TESs by chromosome
	allends[chrom] = {}
	allends[chrom]['tss'] = {}
	allends[chrom]['tes'] = {}
	towrite = {}  # isoforms to be written
	towrite[chrom] = {}
	bedObjs = []
	for jset in isoforms[chrom]:  # unique splice junction chain
		towrite[chrom][jset] = []
		line = list(isoforms[chrom][jset]['line'])
		junccoord = isoforms[chrom][jset]['junccoord']
		jset_ends = find_best_sites(isoforms[chrom][jset]['tss'],
			isoforms[chrom][jset]['tss_tes'], junccoord, chrom)

		if isoformtss and no_redundant == 'longest':
			tss_sorted = sorted(jset_ends, key=lambda x: x[0])[0]  # smallest left coord
			tes_sorted = sorted(jset_ends, key=lambda x: x[1])[-1]  # largest right coord
			if not clean:
				line += [(tss_sorted[3]+tes_sorted[3])/2.]
			tss, tes = tss_sorted[0], tes_sorted[1]
			print('from jset longest', tss, tes, file=sys.stderr)
			towrite[chrom][jset] += [edit_line_bed12(line, tss, tes)]
			# TODO: remove this
			ding = edit_line_bed12(line, tss, tes)
			if len(ding) > 12:
				bedObjs.append(Bed.parse(ding[:13], numStdCols=12))
			else:
				bedObjs.append(Bed.parse(ding))
			continue
		if isoformtss and no_redundant == 'best_only':
			tss, tes = jset_ends[0][:2]
			if not clean:
				line += [jset_ends[0][3]]
			print('from jset best_only', tss, tes, file=sys.stderr)
			towrite[chrom][jset] += [edit_line_bed12(line, tss, tes)]
			# TODO: remove this
			ding = edit_line_bed12(line, tss, tes)
			if len(ding) > 12:
				bedObjs.append(Bed.parse(ding[:13], numStdCols=12))
			else:
				bedObjs.append(Bed.parse(ding))
			continue  # 1 isoform per unique set of junctions

		i = 0
		name = line[3]
		for tss, tes, support, tsscount, tescount in jset_ends:
		#	print('from jset overall, sending', tss, tes, file=sys.stderr)
			edited_line = edit_line_bed12(line, tss, tes)
			if not clean:
				edited_line += [support]
			if i >= 1:  # to avoid redundant names for isoforms with the same junctions
				if '_' in name:
					newname = name[:name.rfind('_')]+'-'+str(i)+name[name.rfind('_'):]
				else:
					newname = name+'-'+str(i)
				edited_line[3] = newname
			if isoformtss:
				towrite[chrom][jset] += [edited_line]
				# TODO: remove this
				if clean:
					bedObjs.append(Bed.parse(edited_line[:12]))
				else:
					bedObjs.append(Bed.parse(edited_line[:13], numStdCols=12))

			else:  # all isoforms will go through another pass to homogenize ends
				towrite[chrom][jset] += [edited_line + [junccoord]]
				# TODO: remove this
				if clean:
					bedObjs.append(Bed.parse(edited_line[:12]))
				else:
					bedObjs.append(Bed.parse(edited_line[:13], numStdCols=12))
				if tss not in allends[chrom]['tss']:
					allends[chrom]['tss'][tss] = 0
				allends[chrom]['tss'][tss] += tsscount
				if tes not in allends[chrom]['tes']:
					allends[chrom]['tes'][tes] = 0
				allends[chrom]['tes'][tes] = tescount
			if no_redundant != 'longest':
				i += 1
	if isoformtss:
		return bedObjs
		#return towrite

	bedObjs = []
	new_towrite = {}  # another pass through all isoforms, moving tss/tes to be more uniform within a gene
	new_towrite[chrom] = {}
	for jset in towrite[chrom]:
		new_towrite[chrom][jset] = []
		jsetends = set()
		endpair = [1e15, 0, 0]  # only applies if no_redundant was specified
		for line in towrite[chrom][jset]:  # adjust isoform TSS/TES using allends dict
			line = list(line)
			junccoord = line[-1]
			tss, tes = int(line[1]), int(line[2])
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
				if no_redundant == 'best_only':
					endpair = [tss, tes, support] if endpair[2] < tss_support else endpair
				elif no_redundant == 'longest':
					endpair[0] = tss if tss < endpair[0] else endpair[0]
					endpair[1] = tes if tes > endpair[1] else endpair[1]
				else:
		#			print('from unfound, sending', tss, tes, file=sys.stderr)
					newline = edit_line_bed12(line, tss, tes)
					if clean:
						bedObjs.append(Bed.parse(newline[:12]))
					else:
						bedObjs.append(Bed.parse(newline[:13], numStdCols=12))
					new_towrite[chrom][jset] += [newline[:-1]]
		if endpair[0] != 1e15:  # write best_only or longest option isoforms
			print('endpair', tss, tes, file=sys.stderr)
			newline = edit_line_bed12(line, endpair[0], endpair[1])
			new_towrite[chrom][jset] += [newline[:-1]]
			if clean:
				bedObjs.append(Bed.parse(newline[:12]))
			else:
				bedObjs.append(Bed.parse(newline[:13], numStdCols=12))
	return bedObjs

def bed_format(bedfields, start, end, singleExon=False):
	if singleExon is True:
		bedObj = Bed( [BedBlock(start, end)])


def edit_line_bed12(line, tss, tes, blocksize=''):
	# for single exons this is mostly a a bed line, but multiexon has a bunch of extra fields.
	# the code as is pretends that everything is ORF, which isn't right.
	line = list(line)
	if blocksize:
		bedObj = Bed(line[0], int(line[1]), int(line[2]), name=line[3], strand=line[5], 
			blocks=[BedBlock(tss, tes)], thickStart=int(tss), thickEnd=int(tes), numStdCols=12) 
		
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
if gtfname:
	gtf = open(gtfname, 'r')
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
	if not quiet:
		sys.stderr.write('Annotated ends extracted from GTF\n')

isoforms = {}  # spliced isoforms
all_se_by_chrom = {}  # single-exon reads grouped by exact start and end site
nosplice_chroms = set(nosplice.split(','))
query = open(queryfile, 'r')
for line in query:
	line = tuple(line.rstrip().split('\t'))
	chrom, tss, tes = line[0], int(line[1]), int(line[2])
	junctions, junccoord = get_junctions_bed12(line)

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

if not quiet:
	sys.stderr.write('Read data extracted\n')

chrom_names = []  # sorted by descending total number of unique single-exon isoforms
for chrom in all_se_by_chrom:
	chrom_names += [(chrom, len(all_se_by_chrom[chrom]))]
chrom_names = [chrom for chrom,num in sorted(chrom_names, key=lambda x: x[1], reverse=True)]
res = []
#if __name__ == '__main__':
#	p = Pool(threads)
#	res = p.map(run_iterative_add_se, chrom_names)
#	p.close()
#	p.join()
for chrom in chrom_names:
	res.append(run_iterative_add_se(chrom))

all_se_by_chrom = None

singleexon = {}  # single-exon isoforms
for r in res:  # combine results
	singleexon.update(r)

if not quiet:
	sys.stderr.write('Single-exon genes grouped, collapsing\n')

# print all the results to a bed file. 
# nb: this does not output a sorted file
with open(outputfname, 'wt') as outfile:

	# single exon genes 
	res = []
	for chrom in chrom_names:
		outlines = run_se_collapse(chrom)
		for edited_line in outlines:
			if clean:
				bedObject = Bed.parse(edited_line[:12])
			else:
				bedObject = Bed.parse(edited_line[:13], numStdCols=12)
			bedObject.write(outfile)

	# multiexon genes
	res = []
	for chrom in isoforms.keys():
		bedObjs = run_find_best_sites(chrom)
		[b.write(outfile) for b in bedObjs]
	
