#! /usr/bin/env python3
import sys, csv,  os, argparse, pysam, gzip
from multiprocessing import Pool
import ast

parser = argparse.ArgumentParser(description='options',
    usage='''python assign_variants_to_transcripts.py -i reads.bam -v variants.vcf
    --map_out out.map --bed_out out.bed [options] > vcf.out ''')
required = parser.add_argument_group('required named arguments')
required.add_argument('--bam', action='store', dest='bam', type=str, required=True,
    help='input bam of aligned reads to be assigned')
required.add_argument('-i', action='store', dest='i', type=str, required=True,
    help='input bed of isoform models')
required.add_argument('-v', '--vcf', action='store', dest='vcf', type=str, required = True,
    help='called variants vcf file')
required.add_argument('--map', action='store', dest='map', type=str, required=True,
    help='''isoform-read map file''')
required.add_argument('--bed_out', action='store', dest='bed_out', type=str, required=True,
  help='''updated isoform file for the flair-collapse isoforms''')
required.add_argument('--map_out', action='store', dest='map_out', type=str, required=True,
  help='''updated isoform-read map file for the flair-collapse isoforms''')
parser.add_argument('-s', '--support', default=3, action='store', dest='support', type=float,
  help='minimum absolute number of supporting reads supporting a haplotype-specific isoform (3)')
parser.add_argument('--support_fraction', default=0.05, action='store', dest='support_fraction', type=float,
  help='''minimum proportion of supporting reads supporting a haplotype-specific isoform in
  a gene (0.05)''')
parser.add_argument('-t', '--threads', default=4, action='store', dest='t', type=int,
  help='number of threads (4)')
args = parser.parse_args()


isbed = args.i[-3:].lower() != 'psl'
isoform_to_chrom = {}

bed_lines = {}
for line in open(args.i):
    line = line.rstrip().split('\t')
    if isbed:
        chrom, name = line[0], line[3]
    else:
        chrom, name = line[13], line[9]
    isoform_to_chrom[name] = chrom

    bed_lines[name] = line

## building an isoform to read map dictionary separated by chromosome for parallelizability
isoform_to_reads = {}
read_to_isoform = {}
chrom_names = set()
for line in open(args.map):
    name, reads = line.rstrip().split('\t')

    ## iso was filtered out for low coverage but is still present in map file
    if name not in isoform_to_chrom:
        continue
    chrom = isoform_to_chrom[name]
    chrom_names.add(chrom)

    if chrom not in isoform_to_reads:
        isoform_to_reads[chrom] = {}
        read_to_isoform[chrom] = {}
    if name not in isoform_to_reads[chrom]:
        isoform_to_reads[chrom][name] = set()
    for r in reads.split(','):
        isoform_to_reads[chrom][name].add(r)
        read_to_isoform[chrom][r] = name


if args.vcf[-3:] == '.gz':
    vcf = gzip.open(args.vcf, 'rt')
else:
    vcf = open(args.vcf)

all_variants = {}
for line in vcf:
    if line.startswith('#'):
        continue
    line = line.rstrip().split('\t')
    chrom, pos = line[0], int(line[1])
    ref, alt = line[3], line[4]
    this_pos = chrom + ':' + str(pos)

    if len(ref) > 1 or len(alt) > 1:
        continue
    if chrom not in all_variants:
        all_variants[chrom] = {}
    all_variants[chrom][pos] = (ref, alt)


def determine_phase_sets(chrom):
    read_variants = {}
    new_isos = {}
    new_map = {}
    variant_to_isos = {}
    ## cataloging amount of support for each variant in the vcf on this chromosome
    if chrom in all_variants:
        samfile = pysam.AlignmentFile(args.bam, "rb" )
        for pileupcolumn in samfile.pileup(chrom, min_base_quality=3):
            this_chrom = pileupcolumn.reference_name

            pos = pileupcolumn.pos + 1
            chrom_pos = chrom+':'+str(pos)
            if pos not in all_variants[chrom]:
                continue

            for pileupread in pileupcolumn.pileups:
                read_name = pileupread.alignment.query_name

                if not pileupread.is_del and not pileupread.is_refskip:
                    base =  pileupread.alignment.query_sequence[pileupread.query_position]

                    ## at this position, does this read contain the non-reference allele?
                    if base == all_variants[chrom][pileupcolumn.pos+1][1]:
                        if read_name not in read_variants:
                            read_variants[read_name] = set()
                        read_variants[read_name].add( (pos, 'A'))

    for isoform in isoform_to_reads[chrom]:
        phase_sets = {}
        for read in isoform_to_reads[chrom][isoform]:

            ## check if each supporting read for this isoform contains variants
            if read in read_variants:
                this_set = str(sorted(read_variants[read]))
                if this_set not in phase_sets:
                    phase_sets[this_set] = []
                phase_sets[this_set] += [read]

        phase_set_index = 0
        phased_read_names = set()
        ## check if enough reads contain the same set of variants
        for phase_set, reads in sorted(phase_sets.items(), key=lambda x: len(x[1]), reverse=True):
            if len(reads) < args.support or \
            float(len(reads)) / len(isoform_to_reads[chrom][isoform]) < args.support_fraction:
                break

            phase_set_index += 1
            new_iso_name = isoform + '_PS'+str(phase_set_index)
            new_isos[new_iso_name] = (isoform, phase_set)  # old name, contained variants
            new_map[new_iso_name] = reads
            phased_read_names.update(reads)

            for variant in ast.literal_eval(phase_set):
                chrom_pos, status = variant
                if status == 'A':
                    if chrom_pos not in variant_to_isos:
                        variant_to_isos[chrom_pos] = set()
                    variant_to_isos[chrom_pos].add(new_iso_name)

        ## check if there are enough reads to support the reference bases-only version of isoform
        if len(isoform_to_reads[chrom][isoform]) - len(phased_read_names) >= args.support:
            new_isos[isoform] = 'NA'
            new_map[isoform] = isoform_to_reads[chrom][isoform] - phased_read_names

    res_dict = {}
    res_dict[(chrom, 'map')] = new_map
    res_dict[(chrom, 'vcf')] = variant_to_isos
    res_dict[(chrom, 'bed')] = new_isos
    return res_dict


if __name__ == '__main__':
    p = Pool(args.t)
    res = p.map(determine_phase_sets, chrom_names)
    p.terminate()


with open(args.map_out, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    for chrom_res in res:
        for chrom, dict_name in chrom_res:
            if dict_name == 'map':
                for iso in chrom_res[(chrom, dict_name)]:
                    writer.writerow([iso, ','.join(chrom_res[(chrom, dict_name)][iso])])

with open(args.bed_out, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    for chrom_res in res:
        for chrom, dict_name in chrom_res:
            if dict_name == 'bed':
                for iso in chrom_res[(chrom, dict_name)]:
                    if chrom_res[(chrom, dict_name)][iso] == 'NA':
                        writer.writerow(bed_lines[iso])
                    else:
                        no_var_iso_name = chrom_res[(chrom, dict_name)][iso][0] 
                        line = bed_lines[no_var_iso_name]
                        line[3] = iso
                        writer.writerow(line)

if args.vcf[-3:] == '.gz':
    vcf = gzip.open(args.vcf, 'rt')
else:
    vcf = open(args.vcf)

for line in vcf:
    line = line.rstrip().split('\t')
    if line[0].startswith('#'):
        print('\t'.join(line))
    else:
        chrom = line[0]
        pos = int(line[1])

        for chrom_res in res:
            if (chrom, 'vcf') in chrom_res:
                if pos not in chrom_res[(chrom, 'vcf')]:
                    continue
                line[-1] += ':'+','.join(chrom_res[(chrom, 'vcf')][pos])
                print('\t'.join(line))
