import sys, csv, re

try:
	sam = open(sys.argv[1])
	chromsizefile = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py samfile chromsizes outpsl\n')
	sys.stderr.write('written for minimap sams\n')
	sys.exit(1)

def cigar_to_blocks(matches):  # parses cigar string matches, returns columns 19 and 20 of psl
	return 

chromsizes = {}
for line in chromsizefile:
	line = line.rstrip().split('\t')
	chromsizes[line[0]] = line[1]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in sam:
		if line.startswith('@'):
			continue
		line = line.rstrip().split('\t')
		qname, flag, tname, pos, cigar, seq, qual = line[0], int(line[1]), line[2], int(line[3]), line[5], line[9], line[10]
		if tname == '*':
			continue
		pos = pos - 1
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
					sys.stderr.write(op + '\n')
		qend = qconsumed + qstart
		ncount = seq.count('N')
		qsize = len(seq)
		# if qsize == 1:
		qsize = qsize_backup
		tsize = chromsizes[tname]  # chromosome length
		tstart = pos
		strand = '-' if flag & 0x10 else '+'  # flag&0x10 is 1 when the strand is -
		blockstarts = [str(pos + s) for s in relblockstarts]
		blockcount = len(blockstarts)
		qstarts = ','.join([str(qstart + s) for s in qstarts]) + ','
		blocksizes = ','.join([str(s) for s in blocksizes]) + ','
		relblockstarts = ','.join([str(s) for s in relblockstarts]) + ','
		blockstarts = ','.join(blockstarts) + ','
		mismatches = qbaseinsert = qnuminsert = tnuminsert = tbaseinsert = 0
		writer.writerow([matchlen, mismatches, 0, ncount, qnuminsert, qbaseinsert, \
			tnuminsert, tbaseinsert, strand, qname, qsize, qstart, qend, \
			tname, tsize, tstart, tend, blockcount, blocksizes, qstarts, blockstarts])