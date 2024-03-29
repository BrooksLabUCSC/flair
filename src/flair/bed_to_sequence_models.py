#!/usr/bin/env python3
import sys, csv, os, argparse, pysam, subprocess
from collections import namedtuple

# PLEASE NOTE: the models section does not work, keeping this file in case we want to fix this later
# models are not called in the flair pipeline

def main():
	parser = argparse.ArgumentParser(description='options',
		usage='python script.py psl|bed genome.fa outfilename [options]')
	parser.add_argument('psl', type=str,
		help='isoforms in psl or bed format')
	parser.add_argument('genome', type=str,
		help='genomic sequence')
	parser.add_argument('outfilename', type=str,
		help='Name of output file')
	parser.add_argument('-v', '--vcf', type=str, 
		help='vcf file with flair phased transcripts')
	# longshot phased arguments
	parser.add_argument('--isoform_haplotypes', type=str, 
		help='isoform haplotype assignments')
	parser.add_argument('--vcf_out',  default='', type=str, 
		help='''vcf output file name. Vcf output can only be created if the input is vcf and 
		isoform_haplotypes have been provided''')
	parser.add_argument('--models_out', default='', type=str,
		help='isoform psl/bed output file name, will contain additional isoforms created from unphased variants')
	args = parser.parse_args()

	fastq, isbed = False, False
	if args.outfilename[-2:].lower() in ['fq', 'fastq']:
		fastq = True
	if args.psl[-3:].lower() != 'psl':
		isbed = True
	bed_to_sequence(args.outfilename, args.genome, args.psl, args.models_out, args.vcf, 
	args.isoform_haplotypes, isbed=isbed, vcf_out=args.vcf_out, fastq=fastq)


def revcomp(seq):
	revcomp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y',
		'Y':'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
	rev_seq = ''
	for i in reversed(range(len(seq))):
		rev_seq += revcomp_dict[seq[i]]
	return rev_seq

def split_iso_gene(iso_gene):
    if '_' not in iso_gene:
    	return iso_gene, 'NA'
    elif '_chr' in iso_gene:
        splitchar = '_chr'
    elif '_XM' in iso_gene:
        splitchar = '_XM'
    elif '_XR' in iso_gene:
        splitchar = '_XR'
    elif '_NM' in iso_gene:
        splitchar = '_NM'
    elif '_NR' in iso_gene:
        splitchar = '_NR'
    elif '_R2_' in iso_gene:
        splitchar = '_R2_'
    elif '_NC_' in iso_gene:
        splitchar = '_NC_'
    else:
        splitchar = '_'
    iso = iso_gene[:iso_gene.rfind(splitchar)]
    gene = iso_gene[iso_gene.rfind(splitchar)+1:]
    return iso, gene


def get_sequence(entry, seq, isbed=False):
	if isbed:
		start = int(entry[1])
		blockstarts = [int(n) + start for n in entry[11].split(',')[:-1]]
		blocksizes = [int(n) for n in entry[10].split(',')[:-1]]
		strand = entry[5]
	else:
		blocksizes = [int(n) for n in entry[18].split(',')[:-1]]
		blockstarts = [int(n) for n in entry[20].split(',')[:-1]]
		strand = entry[8]
	pulled_seq = ''
	for block in range(len(blockstarts)):
		pulled_seq += seq[blockstarts[block]:blockstarts[block]+blocksizes[block]]
	if strand == '-':
		pulled_seq = revcomp(pulled_seq)
	return pulled_seq

def add_variants_to_seq(variant_list, no_variant_sequence, starts, sizes, strand = '+', chrom='chr1', 
	iso_name=False):
	pulled_seq = ''

	for block in range(len(starts)):
		exon_seq = no_variant_sequence[starts[block]:starts[block]+sizes[block]]
		for v in variant_list:
			if v.pos > starts[block] and v.pos < starts[block]+sizes[block]:
				if v.ref != exon_seq[v.pos-starts[block]-1]:
					print('VCF ref {} does not match genome ref base {}, at {}:{}'.format(v.ref, 
						exon_seq[v.pos-starts[block] - 2:v.pos-starts[block] + 2], v.chrom, v.pos))
				exon_seq = exon_seq[:v.pos-starts[block]-1] + v.alts[0] + exon_seq[v.pos-starts[block]:]

				if iso_name:
					vstring = str(v)
					if vstring not in variant_string_to_record:
						variant_string_to_record[vstring] = v


						used_variants[vstring] = set()

					used_variants[vstring].add(iso_name)

		pulled_seq += exon_seq

	return pulled_seq


def get_sequence_with_variants(entry, seq, pvcf, isbed=False, models_out=False, isoform_haplotypes=False):
	''' Entry is the isoform model line, seq is the genomic sequence for this chromosome'''
	
	unphased_variant_support = 3
	psldata = dict()

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

	# get variants for this haplotype
	if not isoform_haplotypes:
		if name not in iso_to_variants: # TODO: iso_to_variants never gets filled
			pulled_seq = get_sequence(entry, seq, isbed)
			return pulled_seq
		v_to_add = iso_to_variants[name]
		v_to_add.reverse()
	else:
		if chrom not in pvcf.header.contigs:
			variants = []
		else:
			variants = pvcf.fetch(chrom, blockstarts[0], blockstarts[-1]+blocksizes[-1],reopen=True)
		v_to_add = []
		v_to_add_alt = []
		for v in variants:
			sample_name = list(v.samples)[0]
			variant_ps = v.samples[sample_name]['PS']
			variant_gt = v.samples[sample_name]['GT']
			variant_ac = v.info['AC']
			if variant_gt == (1,1):
				v_to_add += [v]
				v_to_add_alt += [v]

			elif models_out and variant_gt == (0,1) and not variant_ps and \
			variant_ac[0] > unphased_variant_support and variant_ac[1] > unphased_variant_support:
				print(entry[0], v.pos, variant_ps, variant_ac)
				v_to_add_alt += [v]

			elif name not in haplotype or variant_ps not in haplotype[name]:

				continue
			else:
				v_to_add += [v]
				v_to_add_alt += [v]
		v_to_add.reverse()  # add variants starting from the end in case of indels 
		v_to_add_alt.reverse()

		if len(v_to_add_alt) != len(v_to_add):
			iso, gene = split_iso_gene(name)
			name_ref, name_alt = '>' + iso+':0_'+gene, '>' + iso+':1_'+gene
		else:
			name_ref = name
	iso_name = False
	if isoform_haplotypes:
		iso_name = name_ref
	pulled_seq = add_variants_to_seq(v_to_add, seq, blockstarts, blocksizes, strand, entry[0], iso_name=iso_name)

	if strand == '-':
		pulled_seq = revcomp(pulled_seq)
	if isoform_haplotypes:
		if not models_out or len(v_to_add_alt) == len(v_to_add):
			return pulled_seq

		alt_seq = add_variants_to_seq(v_to_add_alt, seq, blockstarts, blocksizes, iso_name=name_alt)
		if strand == '-':
			alt_seq = revcomp(alt_seq)
		return name_ref, pulled_seq, name_alt, alt_seq # TODO: This CANNOT be right
	else:
		return pulled_seq


def write_sequences(writer, chrom, pvcf, isbed, psldata, seq, models_out=False, isoform_haplotypes=False, fastq=False):

	models = []
	for entry in psldata[chrom]:

		if isoform_haplotypes:
			pulled_seqs = get_sequence_with_variants(entry, seq, pvcf, isbed, models_out, isoform_haplotypes)
			writer.writerow(['>' + name])
			writer.writerow([pulled_seq])   # TODO: This CANNOT be right

		else:
			name = entry[3] if isbed else entry[9]
			if fastq:
				writer.writerow(['@' + name])
			else:
				writer.writerow(['>' + name])
			pulled_seq = get_sequence(entry, seq, isbed)

			writer.writerow([pulled_seq])
			if fastq:
				writer.writerow(['+'])
				writer.writerow(['@'*len(pulled_seq)])
	return models


def bed_to_sequence(outfilename, genome, pslfile, models_out=False, inputvcf=False, 
	isoform_haplotypes=False, isbed=False, vcf_out=False, fastq=False):

	if inputvcf and not (inputvcf and isoform_haplotypes):
		sys.stderr.write('Must provide both vcf and haplotype information if vcf is provided\n')
		sys.exit(1)
	if vcf_out and not inputvcf:
		sys.stderr.write('Not going to write isoform models without input vcf\n')
		vcf_out = False
	print(inputvcf, models_out, vcf_out, isbed)
	if not inputvcf and models_out or not isbed and models_out:
		sys.stderr.write('Not going to write isoform models without vcf or in psl format\n')
		models_out = False
	psldata = {}
	for line in open(pslfile):  # or bed
		line = line.rstrip().split('\t')
		chrom = line[0] if isbed else line[13]
		if chrom not in psldata:
			psldata[chrom] = []
		psldata[chrom] += [line]

	haplotype = {}  # isoform to haplotype
	if isoform_haplotypes:
		for line in open(isoform_haplotypes):
			line = line.rstrip().split('\t')
			if line[1] != 'NA':
				haplotype[line[0]] = [int(hp) for hp in line[1].split(',')]

	pvcf = False
	if inputvcf:
		pvcf = pysam.VariantFile(inputvcf, 'r')
		try:
			pvcf.fetch(chrom)
		except ValueError:
			if inputvcf[-3:] != '.gz':
				subprocess.check_call(['bgzip', '-c', inputvcf], stdout=open(inputvcf+'.gz', 'w'))
				inputvcf = inputvcf+'.gz'
			subprocess.check_call(['tabix', '-fp', 'vcf', inputvcf])
			pvcf = pysam.VariantFile(inputvcf, 'r')
		used_variants = {}
		variant_string_to_record = {}

	with open(outfilename, 'wt') as outfile1:
		writer1 = csv.writer(outfile1, delimiter='\t', lineterminator=os.linesep)
		seq, chrom = '', ''
		ignore = False
		models_to_write = []

		for line in open(genome):
			line = line.rstrip()
			if line.startswith('>'):
				if not chrom:
					chrom = line.split()[0][1:]
					continue
				if chrom in psldata:  # or bed				
					write_sequences(writer1, chrom, pvcf, isbed, psldata, seq, models_out, isoform_haplotypes, fastq)

				chrom = line.split()[0][1:]
				ignore = chrom not in psldata
				seq = ''
			elif not ignore:
				seq += line.upper()

		if chrom in psldata:  # last chromosome
			write_sequences(writer1, chrom, pvcf, isbed, psldata, seq, isoform_haplotypes, fastq)

	if models_out:
		with open(models_out, 'wt') as outfile2:
			writer2 = csv.writer(outfile2, delimiter='\t', lineterminator=os.linesep)
			for entry in models_to_write:
				writer2.writerow(entry)

	if pvcf and isoform_haplotypes:
		header = pvcf.header
		header.add_meta('FORMAT', items=[('ID',"ISO"), ('Number',1), ('Type','String'),
			('Description','Isoforms')])
		if not vcf_out:
			vcf_out = inputvcf[:-3]+'used_variants.vcf'
		vcf_outfile = pysam.VariantFile(vcf_out, 'w', header=pvcf.header)
		for v in used_variants:
			vline = variant_string_to_record[v]
			vline.samples[list(pvcf.header.samples)[0]]['ISO'] = ','.join(used_variants[v])

			vcf_outfile.write(vline)

if __name__ == "__main__":
    main()
