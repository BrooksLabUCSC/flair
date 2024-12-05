#! /usr/bin/env python3

import sys
import argparse
import os, glob
import pipettor
import pysam

os.environ['OPENBLAS_NUM_THREADS'] = '1'



##NEED: reads manifest, transcriptome fasta, transcriptome bed??

def getvariants():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--manifest', required=True, type=str,
                        help="path to manifest files that points to sample names + bam files aligned to transcriptome. Each line of file should be tab separated.")
    parser.add_argument('-o', '--output_prefix', default='flair',
                        help="path to collapsed_output.bed file. default: 'flair'")
    parser.add_argument('-i', '--isoforms',
                        help="path to transcriptome fasta file")
    parser.add_argument('-b', '--bedisoforms',
                        help="path to transcriptome bed file")
    parser.add_argument('-g', '--genome',
                          type=str, required=True, help='FastA of reference genome')
    parser.add_argument('-f', '--gtf',
                        type=str, required=True,
                        help='GTF annotation file')
    # parser.add_argument('-w', '--endwindow', type=int, default=200,
    #                     help="window for comparing ends of isoforms with the same intron chain. Default:200bp")
    # parser.add_argument('-p', '--minpercentusage', type=int, default=10,
    #                     help="minimum percent usage required in one sample to keep isoform in combined transcriptome. Default:10")
    # parser.add_argument('-c', '--convert_gtf', action='store_true',
    #                     help="[optional] whether to convert the combined transcriptome bed file to gtf")
    # parser.add_argument('-s', '--include_se', action='store_true',
    #                     help='whether to include single exon isoforms. Default: dont include')
    # parser.add_argument('-f', '--filter', default='usageandlongest',
    #                     help='type of filtering. Options: usageandlongest(default), usageonly, none')


    args = parser.parse_args()


    # pipettor.run([('samtools', 'faidx', args.isoforms)])

    samples, bamfiles = [], []
    for line in open(args.manifest):
        sample, bamfile = line.rstrip().split('\t')
        samples.append(sample)
        bamfiles.append(bamfile)

        # ##run longshot on each aligned bam file
        # longshotcmd = ('longshot', '--bam', bamfile, '--ref', args.isoforms, '--out', sample + '.flairaligned.vcf', '-F', '-P', '0.00001', '-q', '10', '--output-ref'), '-d' '/private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/longshotdebug')
        # pipettor.run([longshotcmd])
        # # break

    # ##combine vcf files for all samples
    # out = open(args.output_prefix + '.flairaligned.longshot.vcf', 'w')
    # vartoalt = {}
    # hasheader = False
    # for sample in samples:
    #     for line in open(sample + '.flairaligned.vcf'):
    #         if line[0] != '#':
    #             line = line.rstrip().split('\t')
    #             refinfo, alt = tuple(line[:4]), line[4]
    #             if refinfo not in vartoalt: vartoalt[refinfo] = set()
    #             vartoalt[refinfo].add(alt)
    #         elif not hasheader:
    #             out.write(line)
    #     hasheader = True
    # for r in vartoalt:
    #     alts = ','.join(vartoalt[r])
    #     l = [r[0], int(r[1]), r[2], r[3], alts]
    #     out.write('\t'.join([str(x) for x in l]) + '\n')
    # out.close()

    # ###load vcf variants - could skip this and instead load variants dict directly when combining vcfs above
    # ###for longshot variant calling, can force it to output ‘genotype’ snvs, but it refuses to output alt allele
	# ###current workaround is to make the longshot vcf report positions only, but that doesn’t seem good long-term
    # vcfvars = {}
    # for line in open(args.output_prefix + '.flairaligned.longshot.vcf'):
    #     if line[0] != '#':
    #         line = line.split('\t', 6)
    #         chr, pos, id, ref, alt = line[:5]
    #         alt = alt.split(',')
    #         pos = int(pos)
    #         if chr not in vcfvars: vcfvars[chr] = set()
    #         for a in alt:
    #             if len(a) == len(ref):  ##substitution
    #                 # vcfvars[chr].add((pos, 'S', a))
    #                 vcfvars[chr].add((pos, 'S'))
    #             elif len(a) < len(ref):  # deletion
    #                 # vcfvars[chr].add((pos, 'D', len(ref) - len(a)))
    #                 vcfvars[chr].add((pos, 'D'))
    #             else:
    #                 # vcfvars[chr].add((pos, 'I', len(a) - len(ref)))
    #                 vcfvars[chr].add((pos, 'I'))
    #
    #
    # ####this could become really unwieldy with many or large files, but for now I'm trying to not output a text file with this data, since that's also unwieldy
    # isotomutstringtoreads = {}
    # for sindex in range(len(samples)):
    #     sample, bamfile = samples[sindex], bamfiles[sindex]
    #     samfile = pysam.AlignmentFile(bamfile, 'rb')
    #     c = 0
    #     for s in samfile:
    #         if not s.is_secondary and not s.is_supplementary and s.is_mapped:
    #             c += 1
    #             if c % 100000 == 0: print(c, 'reads checked')
    #             thesemutations = []
    #             cigar = s.cigartuples
    #             alignstart = s.reference_start
    #             ref, quer = 0, 0
    #             readseq = s.query_sequence
    #             for block in cigar:
    #                 if block[0] in {0, 7, 8}:  # match, consumes both
    #                     ref += block[1]
    #                     quer += block[1]
    #                 elif block[0] in {1, 4}:  # consumes query ###1 is insertion
    #                     if block[0] == 1: thesemutations.append((ref + alignstart + 1, 'I', readseq[quer:quer+block[1]])) ##insertion has the inserted base sequence
    #                     quer += block[1]
    #                 elif block[0] in {2, 3}:  # consumes reference ##2 is deletion
    #                     if block[0] == 2: thesemutations.append((ref + alignstart + 1, 'D', block[1])) ###deletion has a number
    #                     ref += block[1]
    #             alignedbases = s.get_aligned_pairs(with_seq=True, matches_only=True)
    #             for i in alignedbases:
    #                 if i[2] and i[2].islower():
    #                     thesemutations.append((i[1] + 1, 'S', readseq[i[0]]))
    #             filteredmuts = []
    #             for m in thesemutations:
    #                 if (s.reference_name in vcfvars and m[:2] in vcfvars[s.reference_name]) or (
    #                         m[1] == 'D' and m[2] >= 10) or (m[1] == 'I' and len(m[2]) >= 10): filteredmuts.append(m)
    #             isoname = s.reference_name
    #             if isoname not in isotomutstringtoreads: isotomutstringtoreads[isoname] = {}
    #             thesemutations = frozenset(filteredmuts)
    #             if thesemutations not in isotomutstringtoreads[isoname]: isotomutstringtoreads[isoname][thesemutations] = set()
    #             isotomutstringtoreads[isoname][thesemutations].add(sample + '__' + s.query_name)
    #     samfile.close()
    # print('done parsing reads')
    #
    # isotoseq = {}
    # isoname = None
    # for line in open(args.isoforms):
    #     if line[0] == '>': isoname = line[1:].rstrip()
    #     else: isotoseq[isoname] = line.rstrip()
    #
    # isomutsout = open(args.output_prefix + '.isoswithvars.counts.tsv', 'w')
    # isomutsfa = open(args.output_prefix + '.isoswithvars.fa', 'w')
    # isomutsout.write('\t'.join(['isoname__muts'] + samples) + '\n')
    # filtervalue = 0.2
    # unassigned = frozenset()
    # c = 0
    # for iso in isotomutstringtoreads:
    #     ##for each isoform, get the count of total reads supporting each individual variant
    #     ##also get total reads for isoform
    #     ##for each individual mutation, figure out which ones pass threshold of support
    #     ##go back through, make filtered set of mutstrings that only keeps mutstrings that pass support
    #     c += 1
    #     totreads = 0
    #     mutcounts = {}
    #     mutsetstoreads = isotomutstringtoreads[iso]
    #
    #     for mutset in mutsetstoreads:
    #         reads = isotomutstringtoreads[iso][mutset]
    #         rcount = len(reads)
    #         totreads += rcount
    #         for m in mutset:
    #             if m not in mutcounts: mutcounts[m] = 0
    #             mutcounts[m] += rcount
    #
    #     goodsupport = {}
    #     updatedfiltervalue = filtervalue
    #     while updatedfiltervalue < 1 and len(mutsetstoreads.keys()) > 0:
    #         goodmuts = set()
    #         for i in mutcounts:
    #             if mutcounts[i] > 2 and mutcounts[i] / totreads > updatedfiltervalue: goodmuts.add(i)
    #         filtmutsets = {}
    #         filtsettoog = {}
    #         for mutset in mutsetstoreads:
    #             filtset = mutset & goodmuts
    #             if filtset not in filtmutsets:
    #                 filtmutsets[filtset] = set()
    #                 filtsettoog[filtset] = set()
    #             filtmutsets[filtset].update(mutsetstoreads[mutset])
    #             filtsettoog[filtset].add(mutset)
    #         for mutset in filtmutsets:
    #             mutsupport = len(filtmutsets[mutset])
    #             if mutset in goodsupport:
    #                 goodsupport[mutset].update(filtmutsets[mutset])
    #                 for ogmutset in filtsettoog[mutset]:
    #                     mutsetstoreads.pop(ogmutset)
    #             elif mutsupport > 1 and mutsupport / totreads > filtervalue:
    #                 goodsupport[mutset] = filtmutsets[mutset]  # mutsupport
    #                 for ogmutset in filtsettoog[mutset]:
    #                     mutsetstoreads.pop(ogmutset)
    #         updatedfiltervalue += 0.1
    #     if len(mutsetstoreads.keys()) > 0:
    #         goodsupport[unassigned] = set() ##unassigned is set to no muts currently
    #         for mutset in mutsetstoreads:
    #             goodsupport[unassigned].update(mutsetstoreads[mutset])
    #     varperisocount = 1
    #     for mutpos in goodsupport:
    #         mutinfo = [list(x) for x in sorted(list(mutpos))]
    #         ##getting read sequence
    #         for x in range(len(mutinfo)):
    #             mutinfo[x][0] = mutinfo[x][0] - 1
    #         mutinfo.sort(reverse=True)
    #         thisseq = list(isotoseq[iso])
    #         realmutpos = []
    #         for pos, muttype, alt in mutinfo:
    #             if muttype == 'S':
    #                 thisseq[pos] = alt
    #             elif muttype == 'D':
    #                 count = alt
    #                 thisseq = thisseq[:pos] + thisseq[pos + count:]
    #                 realmutpos = [(x[0] - count, x[1], x[2]) for x in realmutpos]
    #             elif muttype == 'I':
    #                 count = len(alt)
    #                 thisseq = thisseq[:pos] + list(alt) + thisseq[pos:]
    #                 realmutpos = [(x[0] + count, x[1], x[2]) for x in realmutpos]
    #             realmutpos.append((pos, muttype, alt))
    #         realmutpos.sort()
    #
    #         if len(realmutpos) == 0: muttotext = 'nomuts'
    #         # elif mutpos == unassigned: muttotext = 'unassigned'
    #         else:
    #             muttotext = ','.join(['.'.join([str(y) for y in x]) for x in realmutpos])
    #         outname = str(varperisocount) + '-' + iso + '__' + muttotext
    #
    #         ###goal: output counts table and fa file
    #         stor = {s:0 for s in samples}
    #         for r in goodsupport[mutpos]:
    #             s = r.split('__')[0]
    #             stor[s] += 1
    #         outcounts = [outname] + [str(stor[s]) for s in samples]
    #         isomutsout.write('\t'.join(outcounts) + '\n')
    #
    #         ##writing out fasta line
    #         # isomutsfa.write('>' + outname + '\n')
    #         isomutsfa.write('>' + str(varperisocount) + '-' + iso + ' ' + muttotext + '\n')
    #         isomutsfa.write(''.join(thisseq) + '\n')
    #         varperisocount += 1
    # isomutsout.close()
    # isomutsfa.close()
    # print('done generating variant-aware transcript models, aligning to genome for visualization')
    #
    # ##align isoform models to the genome
    # ###align to transcriptome with --secondary=no
    # mm2_cmd = ('minimap2', '-ax', 'splice', '-s', '40', '-t', '4', '--secondary=no', args.genome, args.output_prefix + '.isoswithvars.fa')
    # samtools_sort_cmd = ('samtools', 'sort', '-o', args.output_prefix + '.isoswithvars.bam', '-')
    # samtools_index_cmd = ('samtools', 'index', args.output_prefix + '.isoswithvars.bam')
    # pipettor.run([mm2_cmd, samtools_sort_cmd])
    # pipettor.run([samtools_index_cmd])

    ###get productivity for transcripts without variants
    # prodcmd = ('predictProductivity', '-i', args.bedisoforms, '-o', args.output_prefix + '.isoforms.productivity', '--gtf', args.gtf, '--genome_fasta', args.genome, '--longestORF')
    # pipettor.run([prodcmd])

    # path = os.path.dirname(os.path.realpath(__file__)) + '/'
    # ##adjust productivity prediction to account for variants
    # prodcmd = ('python3', path + 'predict_aaseq_withvar.py', args.output_prefix + '.isoforms.productivity.info.tsv',
    #             args.output_prefix + '.isoswithvars.fa', args.output_prefix + '.isoswithvars.productivity.info.tsv')
    # pipettor.run([prodcmd])

    ##future desire - get proper bed file with orfs from predicting productivity from variant-aware transcripts

    ##need to get counts table for productivity and aaseq
    isotoprodaaseq, aaseqtoinfo, genetoaaseqtocounts, genetoprodtocounts = {}, {}, {}, {}
    aaseqcount = 1
    for line in open(args.output_prefix + '.isoswithvars.productivity.info.tsv'):
        line = line.rstrip().split('\t')
        isoname, prod = line[0], line[1]
        if len(line) == 2: aaseq = ''
        else: aaseq = line[-1]
        gname = isoname.split('_')[-1]
        if prod != 'PRO': aaseq = 'NOTPRO'
        else:
            if (gname, aaseq) not in aaseqtoinfo:
                aaseqtoinfo[(gname, aaseq)] = [aaseqcount, []]
                aaseqcount += 1
            aaseqtoinfo[(gname, aaseq)][1].append(isoname)
            aaseq = 'aaseq' + str(aaseqtoinfo[(gname, aaseq)][0])
        isotoprodaaseq[isoname] = (gname, prod, aaseq)

    for line in open(args.output_prefix + '.isoswithvars.counts.tsv'):
        line = line.rstrip().split('\t')
        if line[0] != 'isoname__muts':
            isoname = '__'.join(line[0].split('__')[:-1])
            counts = [int(x) for x in line[1:]]
            if isoname in isotoprodaaseq:
                gname, prod, aaseq = isotoprodaaseq[isoname]
                if gname not in genetoprodtocounts:
                    genetoprodtocounts[gname] = {}
                    genetoaaseqtocounts[gname] = {}
                if prod not in genetoprodtocounts[gname]: genetoprodtocounts[gname][prod] = counts
                else: genetoprodtocounts[gname][prod] = [genetoprodtocounts[gname][prod][x] + counts[x] for x in range(len(counts))]
                if aaseq not in genetoaaseqtocounts[gname]: genetoaaseqtocounts[gname][aaseq] = counts
                else: genetoaaseqtocounts[gname][aaseq] = [genetoaaseqtocounts[gname][aaseq][x] + counts[x] for x in range(len(counts))]
        else:
            countssamples = line[1:]

    out = open(args.output_prefix + '.aaseq.key.tsv', 'w')
    for gname, aaseq in aaseqtoinfo:
        aaseqid, tlist = aaseqtoinfo[(gname, aaseq)]
        out.write('\t'.join([gname, 'aaseq' + str(aaseqid), ','.join(tlist), aaseq]) + '\n')
    out.close()

    out = open(args.output_prefix + '.aaseq.counts.tsv', 'w')
    out.write('\t'.join(['aaseqid_gene'] + countssamples) + '\n')
    for gname in genetoaaseqtocounts:
        for id in genetoaaseqtocounts[gname]:
            counts = genetoaaseqtocounts[gname][id]
            out.write('\t'.join([id + '_' + gname] + [str(x) for x in counts]) + '\n')
    out.close()

    prodtypes = ['PRO', 'PTC', 'NST', 'NGO']
    out = open(args.output_prefix + '.productivity.counts.tsv', 'w')
    out.write('\t'.join(['productivity_gene'] + countssamples) + '\n')
    for gname in genetoprodtocounts:
        for p in prodtypes:
            if p in genetoprodtocounts[gname]:
                counts = genetoprodtocounts[gname][p]
                out.write('\t'.join([p + '_' + gname] + [str(x) for x in counts]) + '\n')
    out.close()



if __name__ == "__main__":
    getvariants()



