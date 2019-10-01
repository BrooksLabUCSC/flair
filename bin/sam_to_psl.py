import sys, csv, re

try:
	sam = open(sys.argv[1])
	outfilename = sys.argv[2]
	if len(sys.argv) > 3:
		if sys.argv[3] == 'quick':
			quick = True
			chromsizefile = ''
		else:
			quick = False
			chromsizefile = sys.argv[3]
	else:
		chromsizefile = ''

except:
	sys.stderr.write('usage: script.py samfile outpsl [chromsizefile]\n')
	sys.stderr.write('written for minimap sams\n')
	sys.exit(1)

if chromsizefile:
	chromsizes = {}
	for line in open(chromsizefile):
		line = line.rstrip().split('\t')
		chromsizes[line[0]] = line[1]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in sam:
		if line.startswith('@'):
			continue
		line = line.split('\t')
		tname = line[2]
		if tname == '*':
			continue
		qname, flag, pos, cigar, seq, qual = line[0], int(line[1]), int(line[3]) - 1, line[5], line[9], line[10]
		matches = re.findall('([0-9]+)([A-Z])', cigar)
		matchlen = mismatches = relstart = qstart = qconsumed = 0
		blocksizes, relblockstarts, qstarts = [], [], []
		tend = pos
		qnuminsert = 0
		qbaseinsert = 0
		tnuminsert = 0
		tbaseinsert = 0  # deletion
		qsize_backup = 0
		for m in matches:
			num, op = int(m[0]), m[1]
			if op == 'M':  # consumes reference
				blocksizes += [num]
				relblockstarts += [relstart]
				qstarts += [qconsumed]
				relstart += num
				matchlen += num
				tend += num
				qconsumed += num
				qsize_backup += num
			elif op == 'D':  # consumes reference
				relstart += num
				mismatches += num
				tend += num
				qnuminsert += num
			elif op == 'I':
				qconsumed += num
				tbaseinsert += num
				tnuminsert += 1
				qsize_backup += num
			elif op == 'N':  # consumes reference
				tend += num
				relstart += num
			elif op == 'S':
				if not qstart and not matchlen:
					qstart = num
				qsize_backup += num
			elif op == 'H':
				if not qstart and not matchlen:
					qstart = num
					relstart += num
				qsize_backup += num  # technically does not consume q but useful when comparing a read's secondary alignments
			else:
				sys.stderr.write(op + ' unrecognized\n')
				sys.exit(1)

		blockstarts = ','.join([str(pos + s) for s in relblockstarts]) + ','
		blocksizes = ','.join([str(s) for s in blocksizes]) + ','
		if quick:
			writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, 0, qname, 0, 0, 0, \
				tname, 0, 0, 0, 0, blocksizes, 0, blockstarts])
			continue

		blockcount = len(relblockstarts)
		relblockstarts = ','.join([str(s) for s in relblockstarts]) + ','
		qstarts = ','.join([str(qstart + s) for s in qstarts]) + ','

		qend = qconsumed + qstart
		ncount = seq.count('N')
		qsize = len(seq)
		qsize = qsize_backup
		if chromsizefile:
			tsize = chromsizes[tname]  # chromosome length
		else:
			tsize = 0
		tstart = pos
		strand = '-' if flag & 0x10 else '+'  # flag&0x10 is 1 when the strand is -
		mismatches = qbaseinsert = qnuminsert = tnuminsert = tbaseinsert = 0
		writer.writerow([matchlen, mismatches, 0, ncount, qnuminsert, qbaseinsert, \
			tnuminsert, tbaseinsert, strand, qname, qsize, qstart, qend, \
			tname, tsize, tstart, tend, blockcount, blocksizes, qstarts, blockstarts])

