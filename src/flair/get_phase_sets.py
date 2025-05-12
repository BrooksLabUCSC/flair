#!/usr/bin/env python3
import sys
import csv
import os
import argparse
import pysam

def main():
    parser = argparse.ArgumentParser(description='options',
            usage='python get_phase_sets.py -i isoforms.bed -m isoform_read_map.txt -b reads.bam -o out [options]')
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

def get_phase_sets(isoforms, isoform_reads_map, bam, output, outiso, comprehensive=False):
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
            phase_sets[isoform][(ps_tag, hp_tag)] = 0
        phase_sets[isoform][(ps_tag,hp_tag)] += 1

    if num_reads_in_bam == 0:
        sys.stderr.write('No reads in bam {} (get_phase_sets.py)\n'.format(bam))
        sys.exit(1)
    elif (num_reads_in_bam - num_unassigned)/len(read_isoform) < 0.9:
        sys.stderr.write('{} out of {} reads not assigned to an isoform, {} reads in map\n'.format(num_unassigned, num_reads_in_bam, len(read_isoform)))
        sys.exit(1)

    bam.close()

    #reads_isoforms = None  unused

    with open(output, 'wt') as outfile, open(outiso, 'wt') as outfile_iso:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        writer_iso = csv.writer(outfile_iso, delimiter='\t', lineterminator=os.linesep)
        for i in phase_sets:
            total_supporting_reads = float(sum([ps[1] for ps in phase_sets[i].items()]))

            if ('', '') in phase_sets[i]:
                phase_sets[i].pop(('', ''))

            # sort all phase sets for this iso by number of supporting reads
            pss = sorted(phase_sets[i].items(), key=lambda ps: ps[1], reverse=True)

            if not pss:  # no variants
                writer.writerow([i, 'NA'])
                writer_iso.writerow(isoform_model[i])
                continue

            ps_tag, hp_tag = pss[0][0]
            num_reads = pss[0][1]

            if num_reads / total_supporting_reads > 0.1:
                iso, gene = split_iso_gene(i)
                var1 = iso+'-PS:'+ps_tag+':'+hp_tag+'_'+gene

                if num_reads / total_supporting_reads > 0.8 or comprehensive:  # create a new isoform name
                    var2 = iso+'-PS:NA'+'_'+gene
                    writer.writerow([var1, ps_tag, hp_tag, num_reads, total_supporting_reads])
                    writer.writerow([var2, 'NA'])
                    writer_iso.writerow(isoform_model[i][:3] + [var1] + isoform_model[i][4:])
                    writer_iso.writerow(isoform_model[i][:3] + [var2] + isoform_model[i][4:])

                else:  # only create an isoform with variant
                    writer.writerow([i, ps_tag, hp_tag, num_reads, total_supporting_reads])
                    writer_iso.writerow(isoform_model[i][:3] + [var1] + isoform_model[i][4:])
            elif len(pss) > 1 and abs(num_reads / total_supporting_reads - 0.5) < 0.1 \
                    and abs(pss[1][1] / total_supporting_reads - 0.5) < 0.1:  # two haplotypes
                iso, gene = split_iso_gene(i)
                var1 = iso+'-PS:'+pss[0][0]+'_'+gene
                var2 = iso+'-PS:'+pss[1][0]+'_'+gene
                writer.writerow([var1, pss[0][0]+','+pss[1][0]])
                writer.writerow([var2, pss[0][0]+','+pss[1][0]])
                writer_iso.writerow(isoform_model[i][:3] + [var1] + isoform_model[i][4:])
                writer_iso.writerow(isoform_model[i][:3] + [var2] + isoform_model[i][4:])
            else:  # variants not high enough in frequency
                # print(isoform_model[i][:4], pss, '.', total_supporting_reads)
                writer.writerow([i, 'NA'])
                writer_iso.writerow(isoform_model[i])

if __name__ == "__main__":
    main()
