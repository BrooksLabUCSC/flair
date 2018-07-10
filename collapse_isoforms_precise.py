import sys, csv

try:
	psl = open(sys.argv[1])
	if sys.argv[1][-3:] == 'bed':
		bed = True
	else:
		bed = False
	maxnum = int(sys.argv[2])
	minsupport = float(sys.argv[3])
	outfilename = sys.argv[4]
except:
	sys.stderr.write('usage: script.py psl max_tes/tss_threshold min_support outfilename \n')
	sys.stderr.write('usage: script.py bed max_tes/tss_threshold min_support outfilename \n')
	sys.stderr.write('if min_support < 1 it will be used as a percentage \n')
	sys.exit(1)

def get_junctions(line):
	junctions = set()
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		# junctions += [(starts[b]+sizes[b], starts[b+1])]
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_junctions_bed12(line):
	junctions = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_start_end(line):
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	return starts[0], starts[-1]+sizes[-1]

def find_best_site(sites, find_tss=True):  # sites is a dictionary of key site and value frequency
	total = float(sum(list(sites.values())))
	if total < minsupport:
		return ''
	nearby = dict.fromkeys(sites, 0)  # site, distance away
	bestsite = (0, 0)
	for s in sites:  # calculate number of reads supporting this site within maxnum
		for s_ in sites:
			if abs(s - s_) <= maxnum:
				nearby[s] += sites[s_]
		if nearby[s] > bestsite[1]:
			bestsite = (s, nearby[s])
	# return [bestsite]  # no alternative sites allowed. sorry :c 
	if minsupport < 1 and minsupport > 0.5 and bestsite[1]/total >= minsupport:
		return [bestsite]
	elif total - bestsite[1] < minsupport:  # best
		return [bestsite]
	for s in sites:  # remove bestsite reads
		if abs(s - bestsite[0]) <= maxnum:
			nearby[s] = 0
	other_sites = [bestsite]
	for s in sites:  # find other sites with sufficient coverage
		close = False
		if nearby[s] < minsupport:
			continue
		worseothersite = ''
		for os in other_sites:
			if abs(os[0] - s) <= maxnum * 1.5:
				if nearby[s] > os[1]:
					worseothersite = os
				close = True
				break
		if worseothersite:
			other_sites.remove(worseothersite)
			other_sites += [(s, nearby[s])]  # replacing 
		if not close:
			other_sites += [(s, nearby[s])]  # addingdd AQother_sites = sorted(other_sites, key=lambda x: x[1], reverse=True)
	other_sites = sorted(other_sites, key=lambda x: x[1], reverse=True)	
	if len(other_sites) > 1:
		# print('--')
		# print(other_sites)
		if minsupport < 1:
			other_sites2 = [o for o in other_sites if o[1]/total >= minsupport]
			if not other_sites2:
				if find_tss:
					other_sites2 = sorted(other_sites, key=lambda x: x[0])
				else:
					other_sites2 = sorted(other_sites, key=lambda x: x[0], reverse=True)
				other_sites2 = [other_sites2[0]]
		else:
			other_sites2 = [o for o in other_sites if o[1] >= minsupport]
			bestsite = ''
			# print(other_sites2)
			if not other_sites2:
				if find_tss:
					other_sites2 = sorted(other_sites, key=lambda x: x[0])
				else:
					other_sites2 = sorted(other_sites, key=lambda x: x[0], reverse=True)
				other_sites2 = [other_sites2[0]]
			elif len(other_sites2) > 2:
				if other_sites2[0][1] > other_sites2[1][1]:
					bestsite = other_sites2[0]
					other_sites2 = other_sites2[1:]
			# print(other_sites2, bestsite)
			if len(other_sites2) >= 2:
				if other_sites2[0][1] > other_sites2[1][1]:  # second best
					if bestsite:
						other_sites2 = [bestsite, other_sites2[0]]
					else:
						other_sites2 = [other_sites2[0]]
				elif other_sites2[0][1] == other_sites2[1][1]:  # no clear best, take the furthest
					if find_tss:
						other_sites2 = sorted(other_sites2, key=lambda x: x[0])
					else:
						other_sites2 = sorted(other_sites2, key=lambda x: x[0], reverse=True)
					if bestsite and bestsite != other_sites2[0]:
						other_sites2 = [bestsite, other_sites2[0]]
					else:
						other_sites2 = [other_sites2[0]]
				else:
					other_sites2 = other_sites2[:2]
		# print(other_sites2)
		other_sites = other_sites2
	return other_sites

def find_best_site2(sites_tss_all, sites_tes_all):
	""" sites is a dict of key start site and value frequency, sites_tes are for
	associated ends and their frequencies
	sites_tes_all = {tss: {tes: freq}}"""
	sites_tss = find_best_site(sites_tss_all, find_tss = True)
	if not sites_tss:
		return ''
	sepairs = []
	for tss in sites_tss:
		specific_tes = {}  # the specific end sites associated with this tss start
		for tss_ in sites_tes_all:
			if abs(tss_ - tss[0]) <= maxnum:
				for tes in sites_tes_all[tss_]:
					if tes not in specific_tes:
						specific_tes[tes] = 0
					specific_tes[tes] += sites_tes_all[tss_][tes]
		sites_tes = find_best_site(specific_tes, find_tss = False)
		# if len(sites_tss) > 1:
		# 	print(tss, specific_tes, tes)
		for tes in sites_tes:
			sepairs += [(tss[0], tes[0], tes[1])]
	return sepairs

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
	if bsizes[0] < 0:
		print(tss, tstart, tend, line)
		return ''
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
		line[1] = tss
		line[2] = tes
		return line

	bsizes = [int(x) for x in line[10].split(',')[:-1]]
	bstarts = [int(x) for x in line[11].split(',')[:-1]]
	tstart = int(line[1])  # current chrom start
	tend = int(line[2])
	bsizes[0] += tstart - tss
	bsizes[-1] += tes - tend
	if bsizes[0] < 0:
		print(tss, tstart, tend, line)
		return ''
	bstarts[0] = tss
	line[1] = tss
	line[2] = tes
	line[10] = ','.join([str(x) for x in bsizes])+','
	line[11] = ','.join([str(x) - int(line[1]) for x in bstarts])+','
	return line

def single_exon_pairs(sedict):  # incomplete
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
	tssi = 0
	tesi = 0
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

isoforms = {}
n = 0
singleexon = {}
for line in psl:
	line = tuple(line.rstrip().split('\t'))
	if bed:
		chrom = line[0]
		tss, tes = int(line[1]), int(line[2])
		junctions = get_junctions_bed12(line)
	else:
		chrom = line[13]
		tss, tes = get_start_end(line)
		junctions = get_junctions(line)

	if not junctions:
		junctions = 'se'  # all single exons will be considered together
		if chrom not in singleexon:
			singleexon[chrom] = {}
			singleexon[chrom]['tss'] = {}
			singleexon[chrom]['tes'] = {}
		if tss not in singleexon[chrom]['tss']:
			singleexon[chrom]['tss'][tss] = 0
		singleexon[chrom]['tss'][tss] += 1
		if tes not in singleexon[chrom]['tes']:
			singleexon[chrom]['tes'][tes] = 0
		singleexon[chrom]['tes'][tes] += 1
		continue

	junctions = str(sorted(list(junctions)))  # hashable but still unique
	if chrom not in isoforms:
		isoforms[chrom] = {}
	if junctions not in isoforms[chrom]:
		isoforms[chrom][junctions] = {}
		isoforms[chrom][junctions]['tss'] = {}
		isoforms[chrom][junctions]['tes'] = {}
		isoforms[chrom][junctions]['line'] = line
		isoforms[chrom][junctions]['tss_tes'] = {}
	if tss not in isoforms[chrom][junctions]['tss']:
		isoforms[chrom][junctions]['tss'][tss] = 0#[0, []]  # tss usage count, list of ends
		isoforms[chrom][junctions]['tss_tes'][tss] = {}
	isoforms[chrom][junctions]['tss'][tss] += 1
	if tes not in isoforms[chrom][junctions]['tss_tes'][tss]:
		isoforms[chrom][junctions]['tss_tes'][tss][tes] = 0
	isoforms[chrom][junctions]['tss_tes'][tss][tes] += 1
	if tes not in isoforms[chrom][junctions]['tes']:
		isoforms[chrom][junctions]['tes'][tes] = 0
	isoforms[chrom][junctions]['tes'][tes] += 1

dist = {}
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in isoforms:
		for jset in isoforms[chrom]:
			sepairs = find_best_site2(isoforms[chrom][jset]['tss'], isoforms[chrom][jset]['tss_tes'])
			# if tss:
			# 	print(isoforms[chrom][jset], tss)
			# 	sys.exit()
			# tes = find_best_site2(isoforms[chrom][jset]['tes'])
			if len(sepairs) not in dist:
				dist[len(sepairs)] = 0
			dist[len(sepairs)] += 1
			if not tes:
				continue
			line = isoforms[chrom][jset]['line']
			# if len(tss[0]) and len(tes[0]) == 1:
			# 	line = edit_line(line, tss[0], tes[0])
			# 	writer.writerow(line + [max(tss[1], tes[1])])
			# 	if tss[1] != tes[1]:
					# print(tss[1], tes[1],isoforms[chrom][jset]['tss'], isoforms[chrom][jset]['tes'] )
			# else:
			# if len(tes) > 1 and len(tss) > 1:
			# 	print(line, tes, tss)
			i = 0

			for tss_,tes_,support in sepairs:
				# if tes_[0] < 0 or tss_[0] < 0:
				# 	print(tes_, tss_)
				# print(tss_, tes_)
				if bed:
					templine = edit_line_bed12(list(line), tss_, tes_)
				else:
					templine = edit_line(list(line), tss_, tes_)
				if not templine:
					continue
				if i >= 1:
					if bed:
						templine[3] = templine[3]+'-'+str(i)
					else:
						templine[9] = templine[9]+'-'+str(i)
					writer.writerow(templine + [support])
				else:
					writer.writerow(templine + [support])
				i += 1
	# for chrom in singleexon:
	# 	pairs = single_exon_pairs(singleexon[chrom])
	# 	for p in pairs:
	# 		tss, tes, freq = p
	# 		templine = edit_line(list(line), tss, tes, tes-tss)
	# 		if not templine:
	# 			continue
	# 		writer.writerow(templine + [freq])
print(dist)






