#!/usr/bin/env python3
import os
import argparse
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import matplotlib
matplotlib.use("Agg")
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.colors as mcolors
from flair import FlairInputDataError

def parse_args():
    desc = '''The script will produce two images, one of the isoform models and another of the usage proportions.
    The most highly expressed isoforms across all the samples will be plotted.
    The minor isoforms are aggregated into a gray bar. You can toggle min_reads or
    color_palette to plot more isoforms.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('isoforms', type=str, action='store',
                        help='isoforms in bed format')
    parser.add_argument('counts_matrix', type=str, action='store',
                        help='isoform counts')
    parser.add_argument('gene_name', type=str, action='store',
                        help='''Name of gene, must correspond with the gene names in the
    isoform and counts matrix files''')
    parser.add_argument('-o', action='store', dest='o', type=str,
                        help='prefix used for output files (default=gene_name)')
    parser.add_argument('--min_reads', action='store', dest='min_reads', type=int, default=6,
                        help='minimum number of total supporting reads for an isoform to be visualized (default=6)')
    parser.add_argument('-v', '--vcf', action='store', dest='vcf',
                        help='''VCF containing the isoform names that include each variant
                        in the last sample column''')
    parser.add_argument('--palette', action='store', dest='palette',
                        help='''provide a palette file if you would like to visualize more than 7
                        isoforms at once or change the palette used. each line contains a hex color for each isoform''')
    return parser.parse_args()

hex_colors = ['#ba748a', '#3498db', "#34495e"]
name_colors = ['xkcd:windows blue', 'xkcd:faded green', 'xkcd:dusty purple', 'xkcd:amber']
gray = "xkcd:greyish"
reverse_complement = {'C': 'G', 'G': 'C', 'A': 'T', 'T': 'A'}

def parse_bed(bedfh, names=False, plotany=False, keepiso=set()):
    info = []
    usednames = []
    lowbound, upbound = 1e9, 0
    for line in bedfh:
        line = line.rstrip().split('\t')
        name, start, end = line[3], int(line[1]), int(line[2])
        blocksizes = [int(n) for n in line[10].split(',')[:-1]]
        blockstarts = [int(n) + start for n in line[11].split(',')[:-1]]

        if '_' in name[-4:]:  # for isoforms with productivity appended to the name
            flag = name[name.rfind('_') + 1:]
            if flag not in ['PRO', 'PTC']:
                flag = 'Z'+flag  # ngo and nstop are alphabetically after pro and ptc
            name = name[:name.rfind('_')]
        elif '_PTC' in name:
            flag = 'PTC'
            name = name.replace('_PTC', '')
        elif '_PRO' in name:
            flag = 'PRO'
            name = name.replace('_PRO', '')
        elif '_NST' in name:
            flag = 'Z'
            name = name.replace('_NST', '')
        elif '_NGO' in name:
            flag = 'Z'
            name = name.replace('_NGO', '')

        else:
            flag = 'PRO'

        if name not in keepiso and name[:name.rfind('_')] not in keepiso:
            continue

        strand = line[5]
        lowbound = min(lowbound, start)
        upbound = max(upbound, end)
        info += [[blocksizes, blockstarts, keepiso[name], flag, name]]  # exon sizes, exon starts, color, prod, iso
        usednames += [name]
    upbound += 100  # padding up and downstream of isoforms
    lowbound -= 100
    for i in range(len(info)):
        info[i][1] = [n - lowbound for n in info[i][1]]

    try:
        if names:
            return info, usednames
        else:
            return info, lowbound, upbound, strand, usednames
    except NameError:
        raise NameError('Check that the gene name is present in the bed file')
    except Exception as e:
        raise e


def pack(data, rev=True, color=False, tosort=True):
    starts = [max(d[1]) for d in data] if rev else [min(d[1]) for d in data] # sort by right or left end
    if tosort:
        data = [d for (s,d) in sorted(zip(starts, data))]
    else:
        data = [d for s, d in zip(starts, data)]
    packed = [[data[0]]]
    ends = [max(data[0][1]) + data[0][0][data[0][1].index(max(data[0][1]))]]
    for i in range(1, len(data)):
        min_start = min(data[i][1])
        end = max(data[i][1]) + data[i][0][data[i][1].index(max(data[i][1]))]
        pos = -1
        for j in range(len(packed)):
            if ends[j] + 5 < min_start:  # added the plus 5 for spacing between reads
                pos = j  # pack this read with the read at position j
                break
        if pos >= 0:
            packed[pos] += [data[i]]
            ends[pos] = end
        else:
            packed += [[data[i]]]
            ends += [end]
    return packed


def plot_blocks(data, panel, names, iso_to_variant, upper, lower, strand, base_colors, height=.5, l=0.8):
    panel.set_xlim(1, upper - lower + 2)
    if strand == '-':  # flip axes so that the isoforms are plotted 5' -> 3'
        panel.set_xlim(upper - lower + 2, 1)
    panel.set_ylim(-.6, len(data) * 2-.4)

    panel.tick_params(axis='both', which='both',
                       bottom=False, labelbottom=False,
                       left=False, labelleft=False,
                       right=False, labelright=False,
                       top=False, labeltop=False)

    di = 0  # data index
    #ni = 0  # name index unused
    ew = 1.5  # edgewidth on all blocks
    for i in range(0, len(data)*2, 2):  # each line
        read = data[di]

        for j in range(len(read)):  # plot isoform block
            line = read[j]
            sizes, starts, color, flag, iso_name = line[0], line[1], line[2], line[3], line[4][:line[4].rfind('_')]

            iso_start = data[di][0+j][1][0]
            iso_end = data[di][0+j][1][-1]+data[di][0+j][0][-1]
            if strand == '+':
                reverse_text_alignment = iso_start/(upper - lower) > 0.81
                x = iso_start if not reverse_text_alignment else iso_end
            else:
                reverse_text_alignment = iso_end/(upper - lower) < 0.19
                x = iso_end if not reverse_text_alignment else iso_start

            if reverse_text_alignment:  # plot iso name
                panel.text(x, i-height/2 + 1, iso_name, fontsize=8, ha='right', va='center')
            else:
                panel.text(x, i-height/2 + 1, iso_name, fontsize=8, ha='left', va='center')

            for k in range(len(sizes)):  # each block of each read
                if flag == 'PRO':
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2],
                        sizes[k], height, facecolor=color, linewidth=ew, edgecolor=color, zorder=10)
                    panel.add_patch(rectangle)
                elif flag == 'PTC':
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2],
                        sizes[k], height, facecolor='none', edgecolor=color,
                        linewidth=ew, zorder=10,hatch='////')
                    panel.add_patch(rectangle)
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2],
                        sizes[k], height, facecolor=color, linewidth=0, zorder=10,alpha=0.2)
                    panel.add_patch(rectangle)
                else:
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2],
                        sizes[k], height, facecolor='none', edgecolor=color, linewidth=ew, zorder=10,alpha=1)
                    panel.add_patch(rectangle)
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2],
                        sizes[k], height, facecolor=color, linewidth=0, zorder=10,alpha=0.2)
                    panel.add_patch(rectangle)
                if k > 0:
                    panel.plot([starts[k-1]+sizes[k-1], starts[k]], [i]*2, 'k-', lw=l)

            if iso_name in iso_to_variant:
                for var in iso_to_variant[iso_name]:
                    if strand != '+':
                        basecol = base_colors[reverse_complement[var[3]]]
                    else:
                        basecol = base_colors[var[3]]
                    rectangle = mplpatches.Rectangle([var[1], i - height/2-.2],
                        (upper-lower)/500, height+.4, facecolor=basecol, linewidth=0, zorder=10,alpha=.8)
                    panel.add_patch(rectangle)

        di += 1

def plot_isoform_usage(args):
    args = parse_args()
    bedfh = open(args.isoforms)
    counts_matrix = open(args.counts_matrix)
    if not args.o:
        args.o = args.gene_name

    color_palette = name_colors + hex_colors
    base_colors = {'C': color_palette[0], 'A': color_palette[1], 'G': color_palette[4], 'T': color_palette[5]}
    if args.palette:
        color_palette = []
        for line in open(args.palette):
            color_palette += line.rstrip()

    keepiso = {}  # isoforms that they have a sufficient proportion of reads mapping to them
    sample_ids = counts_matrix.readline().rstrip().split('\t')[1:]
    proportions = []
    totals = [0]*len(sample_ids)

    figwidth = 1 + len(sample_ids)*7/10
    plt.figure(figsize=(figwidth,6))  # proportion usage figure
    figstart = 0.7/figwidth

    panel = plt.axes([figstart, 0.11, 1-figstart-0.02, 0.88], frameon=False)  # plotting the proportion of expression
    panel.tick_params(axis='both',which='both',
                       bottom=True, labelbottom=True,
                       left=True, labelleft=True,
                       right=False, labelright=False,
                       top=False, labeltop=False, labelsize=8)

    gray_bar = [['lowexpr']+[0]*len(sample_ids)+[gray]]  # the minor isoform bar is gray

    for line in counts_matrix:
        line = line.rstrip().split('\t')
        if args.gene_name not in line[0]:
            continue
        counts = [float(x) for x in line[1:]]
        if all(x < args.min_reads for x in counts):
            for i in range(len(sample_ids)):
                gray_bar[0][i+1] += counts[i]   # add to gray bar bc expression is too low
            for i in range(len(sample_ids)):
                totals[i] += counts[i]
            continue
        proportions += [[line[0]]+counts+[sum(counts)]]

    colori = 0
    proportions = sorted(proportions, key=lambda x: x[-1], reverse=True)  # sort by expression
    proportions_color = []
    for i in range(len(proportions)):
        p = proportions[i]
        counts = p[1:-1]
        for i in range(len(sample_ids)):
            totals[i] += counts[i]
        if colori == len(color_palette):  # add to gray bar bc colors ran out
            for i in range(len(sample_ids)):
                gray_bar[0][i+1] += counts[i]
            continue
        proportions_color += [p+[color_palette[colori]]]
        keepiso[p[0]] = color_palette[colori]
        colori += 1

    proportions = [gray_bar[0]] + proportions_color

    if len(proportions) == 1:
        raise FlairInputDataError('''Needs more than 1 isoform with sufficient representation, check gene_name in
            your counts file, then try toggling min_reads''')

    proportions = sorted(proportions, key=lambda x:x[1])[::-1]
    heights = [0]*len(sample_ids)
    for iso in proportions:
        for i in range(len(sample_ids)):
            if totals[i] == 0:
                continue
            percentage = iso[1+i]/totals[i]*100
            rectangle = mplpatches.Rectangle([i+1-0.4, heights[i]],
                        width=0.9, height=percentage, facecolor=iso[-1], linewidth=0)
            panel.add_patch(rectangle)
            if percentage >= 7.5:  # add usage percentage
                panel.text(i+1+.07, heights[i] + percentage/2, str(round(percentage,1))+'%',
                 fontsize=10, ha='center',va='center',color='white')
            if percentage >= 12.25:  # add read num
                if iso[1+i] == 1:
                    panel.text(i+1.45, heights[i]+1.75, str(int(iso[1+i]))+' read',
                     fontsize=6, ha='right',va='center',color='white')
                else:
                    panel.text(i+1.44, heights[i]+1.75, str(int(iso[1+i]))+' reads',
                     fontsize=6, ha='right',va='center',color='white')

            heights[i] += percentage

    xlim = len(sample_ids) + 0.5  # guide lines every 20%
    panel.plot([-1, xlim], [0,     0], 'r--', lw=.75, alpha=.3, zorder=0)
    panel.plot([-1, xlim], [20,   20], 'r--', lw=.75, alpha=.3, zorder=0)
    panel.plot([-1, xlim], [40,   40], 'r--', lw=.75, alpha=.3, zorder=0)
    panel.plot([-1, xlim], [60,   60], 'r--', lw=.75, alpha=.3, zorder=0)
    panel.plot([-1, xlim], [80,   80], 'r--', lw=.75, alpha=.3, zorder=0)
    panel.plot([-1, xlim], [100, 100], 'r--', lw=.75, alpha=.3, zorder=0)

    panel.set_xticks(np.arange(1,xlim))
    panel.set_xticklabels(sample_ids, rotation=20, ha='right')
    panel.set_xlim(0.5, xlim)
    panel.set_ylim(0, 100)
    panel.set_ylabel('Percent Usage', fontsize=12)
    # plt.savefig(args.o+'_proportion.pdf', transparent=True, dpi=600)  # uncomment to output as pdf
    plt.savefig(args.o+'_usage.png', dpi=600)

    # plotting isoform structures
    fig_0 = plt.figure(figsize=(9, 3))

    panel = plt.axes([0.005, 0.015, .99, 0.97], frameon=True)  # annotation


    isoforms, lower, upper, strand, names = parse_bed(bedfh, keepiso=keepiso)
    isoforms = sorted(isoforms,key=lambda x: x[3], reverse=True)  # sort by productivity
    packed = pack(isoforms, rev=False, tosort=False)

    iso_to_variant = {}
    if args.vcf:
        for line in open(args.vcf):
            if line.startswith('#'):
                continue
            line = line.rstrip().split('\t')
            if args.gene_name not in line[-1]:
                continue
            variant_isos = ''.join(line[-1].split(':')[5:]).split(',')
            for k in keepiso:
                if ''.join(k.split(':')) in variant_isos:
                    iso_name = k[:k.rfind('_')]
                    if iso_name not in iso_to_variant:
                        iso_to_variant[iso_name] = []
                    iso_to_variant[iso_name] += [(line[0], int(line[1])-lower, line[3], line[4])]


    plot_blocks(packed, panel, names, iso_to_variant, upper, lower, strand, base_colors, l=1)

    # plt.savefig(args.o+'_bars.pdf', transparent=True, dpi=600)  # uncomment to output as pdf
    plt.savefig(args.o+'_isoforms.png', dpi=600)

def main():
    plot_isoform_usage(parse_args())

if __name__ == "__main__":
    main()
