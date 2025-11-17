#!/usr/bin/env python3
import sys, csv, os, argparse, pysam, subprocess
import pysam
import logging
from flair.flair_transcriptome import revcomp
from flair import FlairInputDataError

def main():
    parser = argparse.ArgumentParser(description='options')
    parser.add_argument('bed', type=str,
            action='store', help='isoforms in bed format')
    parser.add_argument('genome', type=str,
            action='store', help='genomic sequence')
    parser.add_argument('outfilename', type=str,
            action='store', help='Name of output file')
    parser.add_argument('-v', '--vcf', action='store', dest='vcf',
            type=str, help='vcf file with flair phased transcripts')
    # longshot phased arguments
    parser.add_argument('--isoform_haplotypes', action='store', dest='isoform_haplotypes',
            type=str, help='isoform haplotype assignments')
    parser.add_argument('--vcf_out', action='store', dest='vcf_out', default='',
            type=str, help='vcf output file name')

    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        parser.error("No arguments passed, please provide genome and bed isoform files")
    args = parser.parse_args()

    if args.vcf and not (args.vcf and args.isoform_haplotypes):
        raise FlairInputDataError('Must provide both vcf and haplotype information if vcf is provided')


    bed_to_sequence(query=args.bed, genome=args.genome, outfilename=args.outfilename,
             isoform_haplotypes=args.isoform_haplotypes, vcfinput=args.vcf,
             vcf_out=args.vcf_out)


# NOTE: using functions inside bed_to_sequence because all of them rely on 'global' variables
# This really should be rewritten.
def bed_to_sequence(query, genome, outfilename, isoform_haplotypes=False, vcfinput=False,
                    vcf_out=False):
    used_variants = dict()
    variant_string_to_record = dict()
    fastq = outfilename[-2:].lower() in ['fq', 'fastq']
    unphased_variant_support = 3

    beddata = {}
    for line in open(query):  # or bed
        line = line.rstrip().split('\t')
        chrom = line[0]
        if chrom not in beddata:
            beddata[chrom] = []
        beddata[chrom] += [line]

    haplotype = {}  # isoform to haplotype
    if isoform_haplotypes:
        for line in open(isoform_haplotypes):
            line = line.rstrip().split('\t')
            if line[1] != 'NA':
                haplotype[line[0]] = [int(hp) for hp in line[1].split(',')]

    vcf = False
    if vcfinput:
        vcf = pysam.VariantFile(vcfinput, 'r')
        try:
            vcf.fetch(chrom)
        # TODO: check for gz beforehand
        except ValueError:
            if vcfinput[-3:] != '.gz':
                subprocess.check_call(['bgzip', '-c', vcfinput], stdout=open(vcfinput+'.gz', 'w'))
                vcfinput = vcfinput+'.gz'
            subprocess.check_call(['tabix', '-fp', 'vcf', vcfinput])
            vcf = pysam.VariantFile(vcfinput, 'r')

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


    def get_sequence(entry, seq):
        start = int(entry[1])
        blockstarts = [int(n) + start for n in entry[11].rstrip(',').split(',')]
        blocksizes = [int(n) for n in entry[10].rstrip(',').split(',')]
        strand = entry[5]
        pulled_seq = ''
        for block in range(len(blockstarts)):
            pulled_seq += seq[blockstarts[block]:blockstarts[block]+blocksizes[block]]
        if strand == '-':
            pulled_seq = revcomp(pulled_seq)
        return pulled_seq

    def add_variants_to_seq(variant_list, no_variant_sequence, starts, sizes, strand = '+', chrom='chr1', iso_name=''):
        pulled_seq = ''

        for block in range(len(starts)):
            exon_seq = no_variant_sequence[starts[block]:starts[block]+sizes[block]]
            for v in variant_list:
                if v.pos > starts[block] and v.pos < starts[block]+sizes[block]:
                    if v.ref != exon_seq[v.pos-starts[block]-1]:
                        print('VCF ref {} does not match genome ref base {}, at {}:{}'.format(v.ref,
                                exon_seq[v.pos-starts[block] - 2:v.pos-starts[block] + 2], v.chrom, v.pos))
                    exon_seq = exon_seq[:v.pos-starts[block]-1] + v.alts[0] + exon_seq[v.pos-starts[block]:]

                    if isoform_haplotypes:
                        vstring = str(v)
                        if vstring not in variant_string_to_record:
                            variant_string_to_record[vstring] = v


                            used_variants[vstring] = set()

                        used_variants[vstring].add(iso_name)

            pulled_seq += exon_seq

        return pulled_seq


    def get_sequence_with_variants(entry, seq):
        ''' Entry is the isoform model line, seq is the genomic sequence for this chromosome'''
        start = int(entry[1])
        blockstarts = [int(n) + start for n in entry[11].split(',')[:-1]]
        blocksizes = [int(n) for n in entry[10].split(',')[:-1]]
        strand = entry[5]
        name = entry[3]


        # get variants for this haplotype
        if isoform_haplotypes and name in haplotype:
            if chrom not in vcf.header.contigs:
                variants = []
            else:
                variants = vcf.fetch(chrom, blockstarts[0], blockstarts[-1]+blocksizes[-1],reopen=True)
            v_to_add = []
            v_to_add_alt = []
            for v in variants:
                sample_name = list(v.samples)[0]
                variant_ps = v.samples[sample_name]['PS']
                variant_gt = v.samples[sample_name]['GT']
                if variant_gt == (1,1):
                    v_to_add.append(v)
                    v_to_add_alt.append(v)

                elif variant_ps in haplotype[name]:
                    if variant_gt == (0,1):
                        v_to_add.append(v)
                    elif variant_gt == (1,0):
                        v_to_add_alt.append(v)
            v_to_add.reverse()  # add variants starting from the end in case of indels
            v_to_add_alt.reverse()

            iso, gene = split_iso_gene(name)
            iso = ':'.join(iso.split(':')[:-1])
            name_ref, name_alt = '>' + iso+':0|1_'+gene, '>' + iso+':1|0_'+gene
            names = [name_ref, name_alt]
            pulled_seq = [add_variants_to_seq(v_to_add, seq, blockstarts, blocksizes, strand, entry[0], iso_name=name_ref),
                          add_variants_to_seq(v_to_add_alt, seq, blockstarts, blocksizes, strand, entry[0],
                                              iso_name=name_alt)]
        else:
            pulled_seq = [seq]
            names = [name]

        if strand == '-':
            pulled_seq = [revcomp(x) for x in pulled_seq]

        return names, pulled_seq


    def write_sequences(beddata_chrom, seq):

        seenisos = set()
        models = []
        for entry in beddata_chrom:

            name = entry[3]
            if vcfinput:
                iso, gene = split_iso_gene(name)
                iso = ':'.join(iso.split(':')[:-1])
                if iso not in seenisos:
                    names, pulled_seq = get_sequence_with_variants(entry, seq)
                    for i in range(len(pulled_seq)):
                        n,p = names[i], pulled_seq[i]
                        writer.writerow(['>' + n])
                        writer.writerow([p])
                    seenisos.add(iso)

            else:
                if fastq:
                    writer.writerow(['@' + name])
                else:
                    writer.writerow(['>' + name])
                pulled_seq = get_sequence(entry, seq)

                writer.writerow([pulled_seq])
                if fastq:
                    writer.writerow(['+'])
                    writer.writerow(['@'*len(pulled_seq)])
        return models


    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        seq, chrom = '', ''



        genome = pysam.FastaFile(genome)
        for chrom in beddata:
            write_sequences(beddata[chrom], genome.fetch(chrom).upper())



    if vcf and isoform_haplotypes:
        header = vcf.header
        header.add_meta('FORMAT', items=[('ID',"ISO"), ('Number',1), ('Type','String'),
                ('Description','Isoforms')])
        if not vcf_out:
            vcf_out = vcf[:-3]+'used_variants.vcf'
        vcf_outfile = pysam.VariantFile(vcf_out, 'w', header=vcf.header)
        for v in used_variants:
            vline = variant_string_to_record[v]
            vline.samples[list(vcf.header.samples)[0]]['ISO'] = ','.join(used_variants[v])

            vcf_outfile.write(vline)

if __name__ == "__main__":
    main()
