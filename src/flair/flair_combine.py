#! /usr/bin/env python3

import argparse
import logging
from flair.bed_to_gtf import bed_to_gtf
from flair.pycbio.hgdata.bed import BedReader
from statistics import mode

def bedReadToIntronChain(bed):
    introns = []
    for i in range(len(bed.blocks) - 1):
        introns.append((bed.blocks[i].end, bed.blocks[i + 1].start))
    # strand is not accounted for here, all intron chains will be left to right
    return tuple(introns)

def intronChainToestarts(ichain, start, end):
    esizes, estarts = [], [0,]
    for i in ichain:
        esizes.append(i[0] - (start + estarts[-1]))
        estarts.append(i[1] - start)
    esizes.append(end - (start + estarts[-1]))
    return esizes, estarts

def getbestends(isodata):
    bestiso = (None, None, None, None, 0)
    for info in isodata:
        if info[4] > bestiso[4]:
            bestiso = info
        elif info[4] == bestiso[4]:
            if info[1] - info[0] > bestiso[1] - bestiso[0]:
                bestiso = info
    return bestiso

def combineIsos(isolist, endwindow):
    isolist.sort()
    isoendgroups = {}
    laststart, lastend = 0, 0
    currgroup = []
    for isoinfo in isolist:
        start, end = isoinfo[0], isoinfo[1]
        if start - laststart <= endwindow and end - lastend <= endwindow:
            currgroup.append(isoinfo)
        else:
            if len(currgroup) > 0:
                isoendgroups[getbestends(currgroup)] = currgroup
            currgroup = [isoinfo]
        laststart, lastend = start, end
    if len(currgroup) > 0:
        isoendgroups[getbestends(currgroup)] = currgroup
    return isoendgroups


def cleanisoname(isoname):
    # removes PAR_Y from end of isoform IDs
    # this is deprecated in new gencode annot, but required for backwards compatibility
    # the _ are disruptive downstream
    return ''.join(isoname.split('_PAR_Y'))


def combine():  # noqa: C901 - FIXME: reduce complexity
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--manifest', required=True, type=str,
                        help="path to manifest files that points to transcriptomes to combine. Each line of file should be tab separated with sample name, sample type (isoform or fusionisoform), path/to/isoforms.bed, path/to/isoforms.fa, path/to/combined.isoform.read.map.txt."  # noqa: E501
                             " fa and read.map.txt files are not required, although if .fa files are not provided for each sample a .fa output will not be generated")
    parser.add_argument('-o', '--output_prefix', default='flair.combined.isoforms',
                        help="path to collapsed_output.bed file. default: 'collapsed_flairomes'")
    parser.add_argument('-w', '--endwindow', type=int, default=200,
                        help="window for comparing ends of isoforms with the same intron chain. Default:200bp")
    parser.add_argument('-p', '--minpercentusage', type=int, default=10,
                        help="minimum percent usage required in one sample to keep isoform in combined transcriptome. Default:10")
    parser.add_argument('-c', '--convert_gtf', action='store_true',
                        help="[optional] whether to convert the combined transcriptome bed file to gtf")
    parser.add_argument('-s', '--include_se', action='store_true',
                        help='whether to include single exon isoforms. Default: dont include')
    parser.add_argument('--end_filter', default='longest',
                        help='type of filtering transcript ends. Options: longest(default), usage, or none, or a number for the maximum amount of ends allowed for a single splice junction chain')
    parser.add_argument('--min_reads', type=int, default=3,
                        help='min reads from all samples to call isoform')

    args = parser.parse_args()
    manifest = args.manifest
    outprefix = args.output_prefix
    endwindow = int(args.endwindow)
    minpercentusage = int(args.minpercentusage) / 100.

    logging.info('parsing manifest')
    # FIXME: remove fasta file input, add genome input, generate reference by getting sequence from genome
    bedfiles, mapfiles, samples, fafiles = [], [], [], []
    for line in open(manifest):
        line = line.rstrip().split('\t')
        # print(line)
        if not (3 <= len(line) <= 5):
            raise Exception(f'Expected between 3 to 5 columns in manifest, got {len(line)} in {manifest}')
        samples.append(line[0] + '__' + line[1])
        bedfiles.append(line[2])
        if len(line) > 3:
            fafiles.append(line[3])
        else:
            fafiles.append('')           # FIXME: switch to None
        if len(line) > 4:
            mapfiles.append(line[4])
        else:
            mapfiles.append('')

    # all samples have fasta file, so fasta
    generatefa = all([len(x) > 0 for x in fafiles])

    logging.info('loading isoforms from individual samples')
    intronchaintoisos = {}
    sampletoseq = {}
    for i in range(len(samples)):
        sample = samples[i]
        isfusion = sample.split('__')[-1] == 'fusionisoform'
        genetoreads, isotoreads = {}, {}
        if mapfiles[i] != '':
            for line in open(mapfiles[i]):
                iso, reads = line.rstrip().split('\t', 1)
                gene = iso.split('_')[-1]
                numreads = len(reads.split(','))
                if gene not in genetoreads:
                    genetoreads[gene] = numreads
                else:
                    genetoreads[gene] += numreads
                isotoreads[iso] = numreads
        if isfusion:
            fnametoinfo = {}
            for bed in BedReader(bedfiles[i], fixScores=True):
                chr, start, end, strand, isoname = bed.chrom, bed.chromStart, bed.chromEnd, bed.strand, bed.name
                isoname = '_'.join(isoname.split('_')[1:])
                if bed.blockCount > 1:
                    ichain = bedReadToIntronChain(bed)
                else:
                    ichain = chr + '-' + str(int(round(start, -4)))
                if isoname not in fnametoinfo:
                    fnametoinfo[isoname] = []
                fnametoinfo[isoname].append([chr, strand, start, end, ichain])
            for isoname in fnametoinfo:
                gene = isoname.split('_')[-1]
                loci = fnametoinfo[isoname]
                if loci[0][1] == '+':
                    start = loci[0][2]
                    loci[0][2] = None
                else:
                    start = loci[0][3]
                    loci[0][3] = None

                if loci[-1][1] == '+':
                    end = loci[-1][3]
                    loci[-1][3] = None
                else:
                    end = loci[-1][2]
                    loci[-1][2] = None
                ichainid = tuple([tuple(x) for x in loci])
                if mapfiles[i] != '':
                    isousage = isotoreads[isoname] / genetoreads[gene]
                    isocounts = isotoreads[isoname]
                else:
                    isousage, isocounts = 1, 0
                if ichainid not in intronchaintoisos:
                    intronchaintoisos[ichainid] = []
                isoname = cleanisoname(isoname)
                intronchaintoisos[ichainid].append((start, end, sample, isoname, isousage, isocounts))
        else:  # not loading fusion reads
            for bed in BedReader(bedfiles[i], fixScores=True):
                chr, start, end, strand, isoname = bed.chrom, bed.chromStart, bed.chromEnd, bed.strand, bed.name
                gene = isoname.split('_')[-1]
                ichain = None
                if bed.blockCount > 1:  # removing single exon isoforms, may want to add this as a user input option later - although how am I handling single exon isoforms? Are they all getting stored as the same empty intron chain? that seems bad
                    ichain = bedReadToIntronChain(bed)
                elif args.include_se:
                    ichain = chr + '-' + str(int(round(start, -4)))
                if ichain:
                    ichainid = (chr, strand, ichain)
                    if mapfiles[i] != '':
                        isousage = isotoreads[isoname] / genetoreads[gene]
                        isocounts = isotoreads[isoname]
                    else:
                        isousage, isocounts = 1, 0
                    if ichainid not in intronchaintoisos:
                        intronchaintoisos[ichainid] = []
                    isoname = cleanisoname(isoname)
                    intronchaintoisos[ichainid].append((start, end, sample, isoname, isousage, isocounts))

        if generatefa:
            last = None
            sampletoseq[sample] = {}
            for line in open(fafiles[i]):
                if line[0] == '>':
                    last = line[1:].rstrip()
                else:
                    last = cleanisoname(last)
                    sampletoseq[sample][last] = line.rstrip()

    logging.info('combining isoforms')
    finalisostosupport = {}
    isocount = 1
    outbed, outcounts = open(outprefix + '.bed', 'w'), open(outprefix + '.counts.tsv', 'w')
    if generatefa:
        outfa = open(outprefix + '.fa', 'w')
    outmap = open(outprefix + '.isoform.map.txt', 'w')
    # FIXME Need to remove gene from ichain!!
    isomap = {}
    for ichainid in intronchaintoisos:
        # chr, strand, gene, ichain = ichainid
        collapsedIsos = combineIsos(intronchaintoisos[ichainid], endwindow)
        isse = isinstance(ichainid[-1], str)
        isfusion = isinstance(ichainid[0], tuple)
        maxintronchainusage = 0
        totintronchaincounts = 0
        ichainendscount = 1
        ends_for_sorting = []
        for start, end, sample, isoname, isousage, isocounts in collapsedIsos:
            maxisousage = max([x[4] for x in collapsedIsos[(start, end, sample, isoname, isousage, isocounts)]])
            totintronchaincounts += sum([x[5] for x in collapsedIsos[(start, end, sample, isoname, isousage, isocounts)]])
            if maxisousage > maxintronchainusage:
                maxintronchainusage = maxisousage
            if args.end_filter == 'longest':
                ends_for_sorting.append((abs(end - start), start, end, sample, isoname, isousage, isocounts))
            else:
                ends_for_sorting.append((maxisousage, start, end, sample, isoname, isousage, isocounts))
        ends_for_sorting.sort(reverse=True)
        if args.end_filter == 'none' or (maxintronchainusage > minpercentusage and (totintronchaincounts > args.min_reads or totintronchaincounts == 0)):  # the only way for the tot counts to be 0 is if there's no map files provided, allow that
            if args.end_filter in {'longest', 'usage'}:
                ends_for_sorting = [ends_for_sorting[0]]
            elif args.end_filter.is_numeric():
                ends_for_sorting = ends_for_sorting[:int(args.end_filter)]
            for _, start, end, sample, isoname, isousage, isocounts in ends_for_sorting:
                theseisos = collapsedIsos[(start, end, sample, isoname, isousage, isocounts)]
                theseisos.sort(key=lambda x: x[1] - x[0], reverse=True)  # longest first
                maxisousage = max([x[4] for x in theseisos])
                totisocounts = sum([x[5] for x in theseisos])
                if ichainendscount == 1 or (args.end_filter.isnumeric() and (totisocounts > int(args.min_reads) or totisocounts == 0)):
                    if isfusion:
                        outgene = mode([x[3].split('_')[-1] for x in theseisos])
                        outname = 'flairiso' + str(isocount) + '-' + str(ichainendscount) + '_' + outgene
                    else:
                        outname = None
                        # this is for prioritizing annotated transcript names above unannotated transcript names
                        # FIXME breaks if annotation is not gencode/ensembl
                        for i in theseisos:
                            if i[3][:4] == 'ENST' and len(i[3].split('ENSG')[0]) < 25 and len(i[3].split('ENSG')) == 2:
                                outname = str(isocount) + '-' + str(ichainendscount) + '_' + i[3]
                                break
                        if not outname:
                            outgene = None
                            for i in theseisos:
                                if len(i[3].split('ENSG')) > 1:
                                    outgene = 'ENSG' + i[3].split('ENSG')[-1]
                                if not outgene or outgene[:4] != 'ENSG':
                                    if len(i[3].split('chr')) > 1:
                                        outgene = 'chr' + i[3].split('chr')[-1]
                            if not outgene:
                                outgene = mode([x[3].split('_')[-1] for x in theseisos])
                            outname = 'flairiso' + str(isocount) + '-' + str(ichainendscount) + '_' + outgene

                    # output bed line
                    if isfusion:
                        ichainid = [list(x) for x in ichainid]
                        # ichain id is: [chr, strand, start, end, ichain]
                        if ichainid[0][1] == '+':
                            ichainid[0][2] = start
                        else:
                            ichainid[0][3] = start
                        if ichainid[-1][1] == '+':
                            ichainid[-1][3] = end
                        else:
                            ichainid[-1][2] = end
                        for gindex in range(len(ichainid)):
                            chr, strand, fstart, fend, ichain = ichainid[gindex]
                            if isinstance(ichain, str):
                                esizes, estarts = [fend - fstart], [0]
                            else:
                                esizes, estarts = intronChainToestarts(ichain, fstart, fend)
                            outbed.write('\t'.join([chr, str(fstart), str(fend),
                                                    'fusiongene' + str(gindex + 1) + '_' + outname, '1000', strand,
                                                    str(fstart), str(fend), '0', str(len(esizes)),
                                                    ','.join([str(x) for x in esizes]) + ',',
                                                    ','.join([str(x) for x in estarts]) + ',']) + '\n')
                    else:
                        chr, strand, ichain = ichainid
                        if isse:
                            esizes, estarts = [end - start], [0]
                        else:
                            esizes, estarts = intronChainToestarts(ichain, start, end)
                        outbed.write(
                            '\t'.join([chr, str(start), str(end), outname, '1000', strand, str(start), str(end), '0',
                                       str(len(esizes)), ','.join([str(x) for x in esizes]) + ',',
                                       ','.join([str(x) for x in estarts]) + ',']) + '\n')

                    # output sequence
                    if generatefa:
                        isoseq = sampletoseq[sample][isoname]
                        outfa.write('>' + outname + '\n' + isoseq + '\n')
                else:
                    outgene = None
                    for i in theseisos:
                        if i[3].split('_')[-1][:4] == 'ENSG':
                            outgene = i[3].split('_')[-1]
                    if not outgene:
                        outgene = mode([x[3].split('_')[-1] for x in theseisos])
                    outname = 'lowexpiso_' + outgene

                if outname not in isomap:
                    isomap[outname] = []
                isomap[outname].extend([x[2] + '..' + x[3] for x in theseisos])
                # get counts
                if outname not in finalisostosupport:
                    finalisostosupport[outname] = {x: 0 for x in samples}
                for isoinfo in theseisos:
                    finalisostosupport[outname][isoinfo[2]] += isoinfo[5]
                ichainendscount += 1
            isocount += 1
        else:
            for start, end, sample, isoname, isousage, isocounts in collapsedIsos:
                theseisos = collapsedIsos[(start, end, sample, isoname, isousage, isocounts)]
                outgene = None
                for i in theseisos:
                    if i[3].split('_')[-1][:4] == 'ENSG':
                        outgene = i[3].split('_')[-1]
                if not outgene:
                    outgene = mode([x[3].split('_')[-1] for x in theseisos])
                outname = 'lowexpiso_' + outgene
                if outname not in isomap:
                    isomap[outname] = []
                for info in collapsedIsos:
                    theseisos = collapsedIsos[info]
                    isomap[outname].extend([x[2] + '..' + x[3] for x in theseisos])
    outbed.close()
    if generatefa:
        outfa.close()

    for newiso in isomap:
        outmap.write(newiso + '\t' + '\t'.join(isomap[newiso]) + '\n')

    outcounts.write('\t'.join(['ids'] + samples) + '\n')
    for name in finalisostosupport:
        outline = [name]
        for s in samples:
            outline.append(str(finalisostosupport[name][s]))
        outcounts.write('\t'.join(outline) + '\n')
    outcounts.close()
    outmap.close()

    if args.convert_gtf:
        bed_to_gtf(query=outprefix + '.bed', outputfile=outprefix + '.gtf')


def main():
    combine()


if __name__ == "__main__":
    main()
