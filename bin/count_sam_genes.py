import sys, csv

try:
	sam = open(sys.argv[1])
	outfilename = sys.argv[2]
	if len(sys.argv) > 3:
		tlengths = open(sys.argv[3])
	else:
		tlengths = ''
except:
	sys.stderr.write('usage: script.py samfile outfilename [transcript_lengths]\n')
	sys.stderr.write('generates a file to see if there is a length bias\n')
	sys.stderr.write('does not count multiple mappers, consider filtering for primary alignments first\n')
	sys.exit(1)

reads = {}
for line in sam:
	if line.startswith('@'):
		continue
	line = line.rstrip().split('\t')
	read, transcript, cigar = line[0], line[2], line[5]
	if read not in reads:
		reads[read] = {}
		reads[read]['cigar'] = cigar
		reads[read]['transcript'] = [transcript]
	else:
		reads[read]['transcript'] += [transcript]

transcripts = {}  # counts
for r in reads:
	t = reads[r]['transcript']
	if len(t) == 1:
		if t[0] not in transcripts:
			transcripts[t[0]] = {}
			transcripts[t[0]]['counts'] = 0
			transcripts[t[0]]['length'] = 0
		transcripts[t[0]]['counts'] += 1
	# else:  # how to deal with multiple mappers? ie more than one transcript for a read
	# 	for i in range(len(t)):
	# 		if t[i] not in transcripts:
	# 			transcripts[t[i]] = {}
	# 			transcripts[t[i]]['counts'] = 0
	# 			transcripts[t[i]]['length'] = 0
	# 		transcripts[t[i]]['counts'] += 1./len(t)

if tlengths:
	for line in tlengths:
		line = line.rstrip().split('\t')
		if line[0] not in transcripts:
			continue
		transcripts[line[0]]['length'] = line[1]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for t in transcripts:
		if transcripts[t]['length'] == 0:
			writer.writerow([t, transcripts[t]['counts']])
		else:
			writer.writerow([t, transcripts[t]['counts'], transcripts[t]['length']])
