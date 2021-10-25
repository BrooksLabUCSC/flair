#!/usr/bin/env python3
import sys, csv, os, argparse, pysam, subprocess

parser = argparse.ArgumentParser(description='options',
	usage='python script.py psl|bed genome.fa outfilename [options]')
parser.add_argument('psl', type=str,
	action='store', help='isoforms in psl or bed format')
parser.add_argument('genome', type=str,
	action='store', help='genomic sequence')
parser.add_argument('outfilename', type=str,
	action='store', help='Name of output file')
parser.add_argument('--isoform_haplotypes', action='store', dest='isoform_haplotypes',
	type=str, help='isoform haplotype assignments')
parser.add_argument('-v', '--vcf', action='store', dest='vcf',
	type=str, help='vcf file')
parser.add_argument('--vcf_out', action='store', dest='vcf_out', default='',
	type=str, help='vcf output file name')
# parser.add_argument('--models_out', action='store', dest='iso_out', default='', type=str,
# 	help='isoform psl/bed out, will contain additional isoforms created from unphased variants')
args = parser.parse_args()

fastq = args.outfilename[-2:].lower() in ['fq', 'fastq']
isbed = args.psl[-3:].lower() != 'psl'

if args.vcf and not (args.vcf and args.isoform_haplotypes):
	sys.stderr.write('Must provide both vcf and haplotype information if vcf is provided\n')
	sys.exit(1)

psldata = {}
for line in open(args.psl):  # or bed
	line = line.rstrip().split('\t')
	chrom = line[0] if isbed else line[13]
	if chrom not in psldata:
		psldata[chrom] = []
	psldata[chrom] += [line]

haplotype = {}  # isoform to haplotype
if args.isoform_haplotypes:
	for line in open(args.isoform_haplotypes):
		line = line.rstrip().split('\t')
		if line[1] != 'NA':
			haplotype[line[0]] = [int(hp) for hp in line[1].split(',')]

if args.vcf:
	vcf = pysam.VariantFile(args.vcf, 'r')
	try:
		vcf.fetch(chrom)
	except ValueError:
		if args.vcf[-3:] != '.gz':
			subprocess.call(['bgzip', '-c', args.vcf], stdout=open(args.vcf+'.gz', 'w'))
			args.vcf = args.vcf+'.gz'
		subprocess.call(['tabix', '-fp', 'vcf', args.vcf])
		vcf = pysam.VariantFile(args.vcf, 'r')
	used_variants = {}
	variant_string_to_record = {}

def get_sequence(entry, seq):
	if isbed:
		start = int(entry[1])
		blockstarts = [int(n) + start for n in entry[11].split(',')[:-1]]
		blocksizes = [int(n) for n in entry[10].split(',')[:-1]]
		strand = entry[5]
	else:
		blocksizes = [int(n) for n in entry[18].split(',')[:-1]]
		blockstarts = [int(n) for n in entry[20].split(',')[:-1]]
		strand = entry[8]
	pulledseq = ''
	for block in range(len(blockstarts)):
		pulledseq += seq[blockstarts[block]:blockstarts[block]+blocksizes[block]]
	if strand == '-':
		pulledseq = revcomp(pulledseq)
	return pulledseq

def get_sequence_with_variants(entry, seq, name):
	if isbed:
		start = int(entry[1])
		blockstarts = [int(n) + start for n in entry[11].split(',')[:-1]]
		blocksizes = [int(n) for n in entry[10].split(',')[:-1]]
		strand = entry[5]
		name = entry[3]
	else:
		blocksizes = [int(n) for n in entry[18].split(',')[:-1]]
		blockstarts = [int(n) for n in entry[20].split(',')[:-1]]
		strand = entry[8]
		name = entry[9]

	if chrom not in vcf.header.contigs:
		variants = []
	else:
		variants = vcf.fetch(chrom, blockstarts[0], blockstarts[-1]+blocksizes[-1],reopen=True)

	# get variants for this haplotype
	v_to_add = []
	for v in variants:
		sample_name = list(v.samples)[0]
		variant_ps = v.samples[sample_name]['PS']
		variant_gt = v.samples[sample_name]['GT']
		variant_ac = v.info['AC']
		if variant_gt == (1,1):
			v_to_add += [v]

		elif name not in haplotype or variant_ps not in haplotype[name]:

			continue
		else:
			v_to_add += [v]
	v_to_add.reverse()  # add variants starting at the downstream coordinate

	pulledseq = ''
	for block in range(len(blockstarts)):
		exon_seq = seq[blockstarts[block]:blockstarts[block]+blocksizes[block]]

		for v in v_to_add:
			if v.pos > blockstarts[block] and v.pos < blockstarts[block]+blocksizes[block]:
				if v.ref != exon_seq[v.pos-blockstarts[block]-1]:
					print('VCF ref {} does not match genome ref base {}'.format(v.ref, 
						exon_seq[v.pos-blockstarts[block] - 2:v.pos-blockstarts[block] + 2]))
				exon_seq = exon_seq[:v.pos-blockstarts[block]-1] + v.alts[0] + exon_seq[v.pos-blockstarts[block]:]
				
				vstring = str(v)
				if vstring not in variant_string_to_record:
					variant_string_to_record[vstring] = v


					used_variants[vstring] = set()

				used_variants[vstring].add(name)
		pulledseq += exon_seq
	if strand == '-':
		pulledseq = revcomp(pulledseq)
	return pulledseq

revcomp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y',
'Y':'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}

def revcomp(seq):
	rev_seq = ''
	for i in reversed(range(len(seq))):
		rev_seq += revcomp_dict[seq[i]]
	return rev_seq

with open(args.outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	seq, chrom = '', ''
	ignore = False
	for line in open(args.genome):
		line = line.rstrip()
		if line.startswith('>'):
			if not chrom:
				chrom = line.split()[0][1:]
				continue
			if chrom in psldata:  # or bed
				for entry in psldata[chrom]:
					name = entry[3] if isbed else entry[9]
					if fastq:
						writer.writerow(['@' + name])
					else:
						writer.writerow(['>' + name])
					if args.vcf:
						pulledseq = get_sequence_with_variants(entry, seq, name)
					else:
						pulledseq = get_sequence(entry, seq)
					writer.writerow([pulledseq])
					if fastq:
						writer.writerow(['+'])
						writer.writerow(['@'*len(pulledseq)])
			chrom = line.split()[0][1:]
			ignore = chrom not in psldata
			seq = ''
		elif not ignore:
			seq += line.upper()

	if chrom in psldata:  # last chromosome
		for entry in psldata[chrom]:
			name = entry[3] if isbed else entry[9]
			if fastq:
				writer.writerow(['@'+name])
			else:
				writer.writerow(['>'+name])
			pulledseq = get_sequence(entry,seq)
			writer.writerow([pulledseq])
			if fastq:
				writer.writerow(['+'])
				writer.writerow(['@'*len(pulledseq)])
if args.vcf:
	header = vcf.header
	header.add_meta('FORMAT', items=[('ID',"ISO"), ('Number',1), ('Type','String'),
		('Description','Isoforms')])
	if not args.vcf_out:
		args.vcf_out = args.vcf[:-3]+'used_variants.vcf'
	vcf_outfile = pysam.VariantFile(args.vcf_out, 'w', header=vcf.header)
	for v in used_variants:
		vline = variant_string_to_record[v]
		vline.samples[list(vcf.header.samples)[0]]['ISO'] = ','.join(used_variants[v])

		vcf_outfile.write(vline)

