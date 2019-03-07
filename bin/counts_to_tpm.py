import sys, csv

try:
	sizefile = open(sys.argv[1])
	counts_matrix = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py fasta.sizes counts_matrix.tsv count_matrix.tpm.tsv\n')
	sys.stderr.write('convenience script for obtaining a file of isoform sizes: bin/fasta_seq_lengths.py\n')
	sys.exit()

sizes = {}
for line in sizefile:
	line = line.rstrip().split('\t')
	sizes[line[0]] = float(line[1])

header = counts_matrix.readline().split('\t')
num_samples = len(header[1:])
matrix_data = [header]
all_rpk = [0] * num_samples
for line in counts_matrix:
	line = line.rstrip().split('\t')
	isoform_id = line[0]
	rpk = [float(count)/sizes[isoform_id] for count in line[1:]]
	for n in range(num_samples):
		all_rpk[n] += rpk[n]
	matrix_data += [[isoform_id] + rpk]

all_rpk = [rpk/1e6 for rpk in all_rpk]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	writer.writerow(matrix_data[0])
	for line in matrix_data[1:]:
		for n in range(num_samples):
			line[n+1] = line[n+1]/all_rpk[n]
		writer.writerow(line)