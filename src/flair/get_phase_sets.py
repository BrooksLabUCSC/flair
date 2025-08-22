#!/usr/bin/env python3
import sys
import csv
import os
import argparse
import pysam
from flair import FlairInputDataError

def main():
    parser = argparse.ArgumentParser(description='options')
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-i', '--isoforms',
            type=str, required=True, help='isoforms in bed format')
    required.add_argument('-m', '--isoform_reads_map',
            type=str, required=True, help='file mapping of supporting reads to isoforms')
    required.add_argument('-b', '--bam',
            type=str, required=True, help='genomic alignment of reads with additional haplotype tags')
    required.add_argument('-o', '--output',
            type=str, required=True, help='output mapping of isoforms to phase sets')
    required.add_argument('--out_iso', action='store', dest='outiso',
            type=str, required=True, help='output new set of isoforms with variants')
    #parser.add_argument('-t', '--threads',
    #       type=int, required=False, default=4, help='number of threads (default=4)')
    parser.add_argument('--comprehensive', action='store_true', dest='comprehensive',
            required=False, default=False, help='always create an isoform model without variants')
    # parser.add_argument('-of', action='store_true', dest='of', \
    #       required=False, help='''Specify this argument to force overwriting of files in
    #       an existing output directory''')
    args = parser.parse_args()
    get_phase_sets(isoforms=args.isoforms, isoform_reads_map=args.isoform_reads_map,
            bam=args.bam, output=args.output, outiso=args.outiso, comprehensive=args.comprehensive)


def split_iso_gene(iso_gene):
    if '_chr' in iso_gene:
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
    elif '_' not in iso_gene:
        return iso_gene, 'NA'
    else:
        splitchar = '_'
    iso = iso_gene[:iso_gene.rfind(splitchar)]
    gene = iso_gene[iso_gene.rfind(splitchar)+1:]
    return iso, gene

def get_phase_sets(isoforms, isoform_reads_map, bam, output, outiso, outmap=None, comprehensive=False):
    isoform_model = {}
    for line in open(isoforms):
        line = line.rstrip().split('\t')
        isoform_model[line[3]] = line

    read_isoform = {}
    phase_sets = {}
    for line in open(isoform_reads_map):
        line = line.rstrip().split('\t')
        reads, isoform = line[1].split(','), line[0]
        phase_sets[isoform] = {}
        for r in reads:
            read_isoform[r] = isoform

    bam = pysam.AlignmentFile(bam, 'rb')
    num_unassigned, num_reads_in_bam = 0, 0.
    for read in bam.fetch():
        if not read.cigarstring:
            continue
        num_reads_in_bam += 1
        if read.query_name not in read_isoform:
            num_unassigned += 1
            continue

        try:
            ps_tag = str(read.get_tag('PS'))
            hp_tag = str(read.get_tag('HP'))
        except KeyError:
            ps_tag = ''
            hp_tag = ''
        isoform = read_isoform[read.query_name]
        if (ps_tag, hp_tag) not in phase_sets[isoform]:
            phase_sets[isoform][(ps_tag, hp_tag)] = set()
        phase_sets[isoform][(ps_tag,hp_tag)].add(read.query_name)

    if num_reads_in_bam == 0:
        raise FlairInputDataError(f'No reads in bam {bam} (get_phase_sets.py)')
    elif (num_reads_in_bam - num_unassigned)/len(read_isoform) < 0.9:
        raise FlairInputDataError(f'{num_unassigned} out of {num_reads_in_bam} reads not assigned to an isoform, '
                         f'{len(read_isoform)} reads in map')

    bam.close()

    #reads_isoforms = None  unused
    if outmap: outmap = open(outmap, 'wt')
    with open(output, 'wt') as outfile, open(outiso, 'wt') as outfile_iso:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        writer_iso = csv.writer(outfile_iso, delimiter='\t', lineterminator=os.linesep)
        for i in phase_sets:
            total_supporting_reads = float(sum([len(ps[1]) for ps in phase_sets[i].items()]))
            if ('', '') in phase_sets[i]:
                phase_sets[i].pop(('', ''))

            # sort all phase sets for this iso by number of supporting reads
            pss = sorted(phase_sets[i].items(), key=lambda ps: ps[1], reverse=True)
            iso, gene = split_iso_gene(i)

            tot_hap_reads = sum([len(x[1]) for x in pss])

            if pss and tot_hap_reads/total_supporting_reads >= 0.1:
                for j in range(len(pss)):
                    ps_tag, hp_tag = pss[j][0]
                    hap_name = iso+'-PS:'+ps_tag+':'+hp_tag+'_'+gene
                    num_reads = len(pss[j][1])
                    writer.writerow([hap_name, ps_tag, hp_tag, num_reads, total_supporting_reads])
                    writer_iso.writerow(isoform_model[i][:3] + [hap_name] + isoform_model[i][4:])
                    if outmap: outmap.write(hap_name + '\t' + ','.join(pss[j][1]) + '\n')
            if not pss or tot_hap_reads/total_supporting_reads < 0.1 or comprehensive: ##not high enough variant frequency or all isoforms requested
                hap_name = iso+'-PS:NA'+'_'+gene
                writer.writerow([hap_name, 'NA'])
                writer_iso.writerow(isoform_model[i][:3] + [hap_name] + isoform_model[i][4:])


if __name__ == "__main__":
    main()
