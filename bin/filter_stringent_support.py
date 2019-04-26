import sys, csv

try:
	isoforms = open(sys.argv[1])
	alignment = open(sys.argv[2])
	minsupport = int(sys.argv[3])
	outfilename = sys.argv[4]
	if len(sys.argv) > 5:
		outfilename2 = sys.argv[5]
	else:
		outfilename2 = ''
except:
	sys.stderr.write('usage: script.py alignment.sam.psl isoforms.psl minsupport out_isoforms.psl [out_assignments.txt]\n')
	sys.exit(1)

isoform_info = {}
for line in isoforms:
	line = line.rstrip().split('\t')
	blocksizes = [float(n) for n in line[18].split(',')[:-1]]
	isoform_info[line[9]] = [sum(blocksizes), blocksizes[0], blocksizes[-1], line]

iso_read = {}  # isoform-read assignments for reads that span 25bp of the first and last exon
for line in alignment:  # reads aligned to the isoforms sam-turned-psl
	line = line.rstrip().split('\t')
	read, isoform = line[9], line[13]  # names

	blocksizes = [int(n) for n in line[18].split(',')[:-1]]
	blockstarts = [int(n) for n in line[20].split(',')[:-1]]
	read_start, read_end = blockstarts[0], blockstarts[-1]+blocksizes[-1]

	info = isoform_info[isoform]
	isoform_length, first_blocksize, last_blocksize = info[0:3]

	right_coverage = left_coverage = False
	if first_blocksize < 25:
		if read_start < 2:
			left_coverage = True
	elif read_start <= (first_blocksize - 25):
		left_coverage = True

	if last_blocksize < 25:
		if (isoform_length - read_end) < 2:
			right_coverage = True
	if (isoform_length-last_blocksize + 25) <= read_end:
		right_coverage = True

	coverage = sum(blocksizes) / isoform_length
	# coverage = proportion of bases of the isoform that the read covers

	if isoform not in iso_read:
		iso_read[isoform] = []
	if right_coverage and left_coverage and coverage > 0.8:  
		iso_read[isoform] += [[read, isoform, coverage]]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for iso in iso_read:
		supporting = iso_read[iso]  # supporting reads
		if len(supporting) >= minsupport:
			writer.writerow(isoform_info[iso][3])

if outfilename2:  # map file
	with open(outfilename2, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t')
		for iso in iso_read:
			supporting = iso_read[iso]
			if len(supporting) >= minsupport:
				for s in supporting:
					writer.writerow(s)
