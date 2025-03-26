#! /usr/bin/env python3

import sys
import argparse
import os, glob, math
import pipettor
import pysam
import pybedtools
import numpy as np

os.environ['OPENBLAS_NUM_THREADS'] = '1'
compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

def getStarts(gtf):
    starts = list()
    with open(gtf) as lines:
        for l in lines:
            if l[0] == "#": continue
            cols = l.rstrip().split("\t")
            chrom, c1, c2, strand = cols[0], int(cols[3])-1, int(cols[4]), cols[6]
            if cols[2] == "start_codon":
                gene = cols[8][cols[8].find('gene_id')+len('gene_id')+2:]
                gene = gene[:gene.find('"')]
                # gene = re.search("(ENSG[^\.]+)", cols[-1]).group(1)

                starts.append((chrom,c1,c2,gene,".",strand))
    if (len(starts)) == 0:
        sys.stderr.write('ERROR, no start codons were found in', gtf)
        sys.exit(1)
    return starts

class Isoform(object):
    '''
    Object to handle isoform related data.

    attributes:

    methods:

    '''

    def __init__(self, name=None):
        self.name = name
        # self.pro = "UNK"
        self.chrom = ""

        self.seqvariants = {}
        self.exons     = set()
        self.starts    = set()
        # self.orfs      = list()
        self.exonSizes = list()
        self.ptcpoint = 0
        self.startcodons = []

class SeqVar(object):
    def __init__(self, name=None, seq=None):
        self.name = name
        self.pro = "UNK"
        self.sequence = seq
        self.orfs = list()

def getStartRelPos(genomicStartPos,exon, exons, isoObj):
    '''
    is handed a genomic position, the exon it occurs in, all exons,
    and returns the position relative to all exons
    '''
    exonNum = exons.index(exon)
    isoObj.exonSizes = [x[1]-x[0] for x in exons]
    relativeStart = None
    # First get start position relative to transcript sequence.
    if isoObj.strand == "+":
        relativeStart = genomicStartPos - exons[exonNum][0] + sum([x for x in isoObj.exonSizes[:exonNum]])
    elif isoObj.strand == "-":
        relativeStart = len(isoObj.sequence) - (genomicStartPos - exons[exonNum][0] + sum([x for x in isoObj.exonSizes[:exonNum]])) - 3

    return relativeStart

def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

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

    outprefix = '/'.join(args.output_prefix.split('/')[:-1])
    if len(outprefix) > 0: outprefix = outprefix + '/'

    ###REMOVE, TESTING
    args.output_prefix = args.output_prefix.split('/')[-1]
    pipettor.run([('samtools', 'faidx', args.isoforms)])

    samples, bamfiles = [], []
    for line in open(args.manifest):
        sample, bamfile = line.rstrip().split('\t')
        samples.append(sample)
        bamfiles.append(bamfile)

        ##run longshot on each aligned bam file
        print('running longshot on', sample)
        # longshotcmd = ('longshot', '--bam', bamfile, '--ref', args.isoforms, '--out', outprefix +  sample + '.flairaligned.vcf', '-F', '-P', '0.00001', '-q', '10', '--output-ref')#, '-d' '/private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/longshotdebug')
        longshotcmd = ('longshot', '--bam', bamfile, '--ref', args.isoforms, '--out', args.output_prefix +  sample + '.flairaligned.vcf', '-F', '-P', '0.00001', '-q', '10', '--output-ref')#, '-d' '/private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/longshotdebug')
        pipettor.run([longshotcmd])
        # ## break

    print(args.output_prefix, outprefix)
    ##combine vcf files for all samples
    out = open(args.output_prefix + '.flairaligned.longshot.vcf', 'w')
    vartoalt = {}
    hasheader = False
    for sample in samples:
        # for line in open(outprefix +  sample + '.flairaligned.vcf'):
        for line in open(args.output_prefix + sample + '.flairaligned.vcf'):
            if line[0] != '#':
                line = line.rstrip().split('\t')
                refinfo, alt = tuple(line[:4]), line[4]
                if refinfo not in vartoalt: vartoalt[refinfo] = set()
                vartoalt[refinfo].add(alt)
            elif not hasheader:
                out.write(line)
        hasheader = True
    for r in vartoalt:
        alts = ','.join(vartoalt[r])
        l = [r[0], int(r[1]), r[2], r[3], alts]
        out.write('\t'.join([str(x) for x in l]) + '\n')
    out.close()
    print('done combining vcf variants')


    isotoblocks = {}
    genetoiso = {}
    for line in open(args.bedisoforms):
        line = line.rstrip().split('\t')
        thischr, iso, dir, start, esizes, estarts = line[0], line[3], line[5], int(line[1]), [int(x) for x in line[10].split(',')[:-1]], [int(x) for x in line[11].split(',')[:-1]]
        exonblocks = []  ##block is gstart, tstart, len
        if dir == '-':
            esizes = esizes[::-1]
            estarts = estarts[::-1]
        currtstart = 0
        for i in range(len(esizes)):
            exonblocks.append((currtstart, estarts[i] + start, esizes[i], thischr, dir))
            currtstart += esizes[i]
        exonblocks.sort()

        if iso[:10] == 'fusiongene': iso = '_'.join(iso.split('_')[1:])
        if iso not in isotoblocks: isotoblocks[iso] = exonblocks
        else:
            tstartadj = isotoblocks[iso][-1][0] + isotoblocks[iso][-1][2]
            exonblocks = [(x[0] + tstartadj,) + x[1:] for x in exonblocks]
            isotoblocks[iso].extend(exonblocks)

        gene = iso.split('_')[-1]
        if gene not in genetoiso: genetoiso[gene] = set()
        genetoiso[gene].add(iso)
    print('done parsing iso blocks')

    ###load vcf variants - could skip this and instead load variants dict directly when combining vcfs above
    ###for longshot variant calling, can force it to output ‘genotype’ snvs, but it refuses to output alt allele
	###current workaround is to make the longshot vcf report positions only, but that doesn’t seem good long-term
    vcfvars = {}
    for line in open(args.output_prefix + '.flairaligned.longshot.vcf'):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            if len(line) == 4: line.append('')
            iso, tpos, id, ref, alt = line[:5]
            alt = alt.split(',')
            tpos = int(tpos)

            for a in alt:
                if len(a) == len(ref):  thistype='S'
                elif len(a) < len(ref): thistype='D'
                else: thistype='i'

                blocks = isotoblocks[iso]
                genomePos = None
                for tstart, gstart, bsize, thischr, dir in blocks:
                    if (dir == '+' and tstart <= tpos < tstart + bsize) or (
                            dir == '-' and tstart <= tpos - 1 < tstart + bsize):
                        if dir == '+':
                            genomePos = gstart + (tpos - tstart)
                        else:
                            genomePos = gstart + ((tstart + bsize) - (tpos - 1))
                        break
                if genomePos:
                    thisgene = iso.split('_')[-1]
                    for iso2 in genetoiso[thisgene]:
                        blocks = isotoblocks[iso2]
                        tpos2 = None
                        for tstart, gstart, bsize, thischr, dir in blocks:
                            if gstart <= genomePos < gstart + bsize:
                                if dir == '+':
                                    tpos2 = tstart + (genomePos-gstart)
                                else:
                                    tpos2 = (tstart+bsize+1) - (genomePos-gstart)
                                break
                        if tpos2:
                            if iso2 not in vcfvars: vcfvars[iso2] = set()
                            vcfvars[iso2].add((tpos2, thistype))
                            # print(iso2, tpos2)





    isotoseq = {}
    isoname = None
    for line in open(args.isoforms):
        if line[0] == '>': isoname = line[1:].rstrip()
        else: isotoseq[isoname] = line.rstrip()

    # isotomutstringtoreads = {}
    if not os.path.exists(args.output_prefix + '_tempvarfiles'): os.makedirs(args.output_prefix + '_tempvarfiles')
    genenames = set()
    ####this could become really unwieldy with many or large files, but for now I'm trying to not output a text file with this data, since that's also unwieldy
    for sindex in range(len(samples)):
        sample, bamfile = samples[sindex], bamfiles[sindex]
        samfile = pysam.AlignmentFile(bamfile, 'rb')
        c = 0
        for s in samfile:
            if s.is_mapped: #and not s.is_supplementary: ##not s.is_secondary and ###trusting FLAIR quantify to only output one alignment per read
                c += 1
                if c % 100000 == 0: print(c, 'reads checked')
                thesemutations = []
                cigar = s.cigartuples
                alignstart = s.reference_start
                ref, quer = 0, 0
                readseq = s.query_sequence
                for block in cigar:
                    if block[0] in {0, 7, 8}:  # match, consumes both
                        ref += block[1]
                        quer += block[1]
                    elif block[0] in {1, 4}:  # consumes query ###1 is insertion
                        if block[0] == 1: thesemutations.append((ref + alignstart + 1, 'I', readseq[quer:quer+block[1]])) ##insertion has the inserted base sequence
                        quer += block[1]
                    elif block[0] in {2, 3}:  # consumes reference ##2 is deletion
                        if block[0] == 2: thesemutations.append((ref + alignstart + 1, 'D', block[1])) ###deletion has a number
                        ref += block[1]
                alignedbases = s.get_aligned_pairs(with_seq=True, matches_only=True)
                for i in alignedbases:
                    if i[2] and i[2].islower():
                        thesemutations.append((i[1] + 1, 'S', readseq[i[0]]))
                filteredmuts = []
                for m in thesemutations:
                    if (s.reference_name in vcfvars and m[:2] in vcfvars[s.reference_name]) or (
                            m[1] == 'D' and m[2] >= 10) or (m[1] == 'I' and len(m[2]) >= 10): filteredmuts.append(m)
                isoname = s.reference_name
                genename = isoname.split('_')[-1].split('.')[0][:-2]
                genenames.add(genename)
                filteredmuts = ['.'.join([str(z) for z in x]) for x in filteredmuts]
                with open(args.output_prefix + '_tempvarfiles/' + genename + '.txt', 'a') as tempvarout:
                    tempvarout.write('\t'.join([isoname, str(sindex) + '__' + s.query_name, ','.join(filteredmuts)]) + '\n')

                # if isoname not in isotomutstringtoreads: isotomutstringtoreads[isoname] = {}
                # thesemutations = frozenset(filteredmuts)
                # if thesemutations not in isotomutstringtoreads[isoname]: isotomutstringtoreads[isoname][thesemutations] = set()
                # isotomutstringtoreads[isoname][thesemutations].add(str(sindex) + '__' + s.query_name)
        samfile.close()
        print('done parsing reads for', sample)

    genenames = set()
    for f in os.listdir(args.output_prefix + '_tempvarfiles/'):
        if f[0] != '.' and 'processed' not in f:
            genenames.add(f.split('.txt')[0])
    # genenames = {'ENSG000001307_copy',}



    genomepostosupport = {}
    for gene in genenames:
        alreadyprocessed = False
        isotomutstringtoreads = {}
        # print(gene)
        with open(args.output_prefix + '_tempvarfiles/' + gene + '.txt', 'r') as mutstringfile:
            for line in mutstringfile:
                isoname, readname, mutstring = line.rstrip('\n').split('\t')
                if '__' in mutstring:
                    alreadyprocessed = True
                    break
                if mutstring == '': thesemutations = tuple()
                else:
                    thesemuts = [x.split('.') for x in mutstring.split(',')]
                    thesemutations = []
                    for m in thesemuts:
                        m[0] = int(m[0])
                        thesemutations.append(tuple(m))
                    thesemutations = tuple(thesemutations)
                if isoname not in isotomutstringtoreads: isotomutstringtoreads[isoname] = {}
                if thesemutations not in isotomutstringtoreads[isoname]: isotomutstringtoreads[isoname][thesemutations] = set()
                isotomutstringtoreads[isoname][thesemutations].add(readname)
        if alreadyprocessed: continue
        ###go through mutstrings for each transcript
        ###convert all vars to genomic positions
        ###get support for each var at the genomic level
        with open(args.output_prefix + '_tempvarfiles/' + gene + '_processed.txt', 'w') as mutstringfile2:
            for iso in isotomutstringtoreads:
                # mutsettogenomepos = {}
                mutsetstoreads = isotomutstringtoreads[iso]
                blocks = isotoblocks[iso]
                for mutset in mutsetstoreads:
                    rcount = len(isotomutstringtoreads[iso][mutset])
                    newmutpos = []
                    for m in mutset: #pos, type, str
                        tpos = m[0]
                        genomeset = None
                        for tstart, gstart, bsize, thischr, dir in blocks:
                            if (dir == '+' and tstart <= tpos < tstart + bsize) or (dir == '-' and tstart <= tpos-1 < tstart + bsize):
                                if dir == '+': genomePos = gstart + (tpos - tstart)
                                else: genomePos = gstart + ((tstart + bsize) - (tpos-1))
                                genomeset = (thischr, dir, str(genomePos), str(tpos), m[1], m[2])
                                break
                        # print(genomeset, m, blocks,)
                        if genomeset:
                            if genomeset not in genomepostosupport: genomepostosupport[genomeset] = rcount
                            else: genomepostosupport[genomeset] += rcount
                            newmutpos.append('..'.join(genomeset))
                    mutstringfile2.write('\t'.join([iso, ','.join(newmutpos), ','.join(isotomutstringtoreads[iso][mutset])]) + '\n')
    print('done preprocessing iso muts')

    isovarcounts, isovarfa = {}, {}
    muttexttogenomeinfo = {}
    for gene in genenames:
        isotomutstringtoreads = {}
        # print(gene)
        with open(args.output_prefix + '_tempvarfiles/' + gene + '_processed.txt', 'r') as mutstringfile3:
            for line in mutstringfile3:
                iso, genomemutpos, readnames = line.rstrip('\n').split('\t')
                if iso not in isotomutstringtoreads: isotomutstringtoreads[iso] = {}
                if genomemutpos == '': genomemutpos = frozenset()
                else: genomemutpos = frozenset((tuple(x.split('..')) for x in genomemutpos.split(',')))
                isotomutstringtoreads[iso][genomemutpos] = set(readnames.split(','))
        for iso in isotomutstringtoreads:
            ##for each isoform, get the count of total reads supporting each individual variant
            ##also get total reads for isoform
            ##for each individual mutation, figure out which ones pass threshold of support
            mutsetstoreads = isotomutstringtoreads[iso]

            allmuts = set()
            for mutset in mutsetstoreads:
                allmuts.update(mutset)
            goodmuts = set()
            for m in allmuts:
                # genomeset = mutsettogenomepos[(iso, m)]
                if genomepostosupport[m] >= 3: goodmuts.add(m)

            goodsupport = {}
            for mutset in mutsetstoreads:
                filtset = mutset & goodmuts
                if filtset not in goodsupport: goodsupport[filtset] = isotomutstringtoreads[iso][mutset]
                else: goodsupport[filtset].update(isotomutstringtoreads[iso][mutset])

            for mutpos in goodsupport:
                # print(mutpos, len(goodsupport[mutpos]))


                # genometext = []
                # for m in mutpos:
                #     # genomeset = mutsettogenomepos[(iso, m)]
                #     genometext.append('.'.join(m))

                mutinfo = sorted([x[:3] + (int(x[3])-1, ) + x[4:] for x in mutpos], reverse=True, key=lambda x:x[3])
                ##getting read sequence
                # for x in range(len(mutinfo)):
                #     mutinfo[x] = list(mutinfo[x])
                #     mutinfo[x][0] = mutinfo[x][0] - 1

                # mutinfo.sort(reverse=True)
                thisseq = list(isotoseq[iso])
                realmutpos = []
                for thischr, strand, genomepos, pos, muttype, alt in mutinfo:
                    # print(len(thisseq), thischr, strand, pos, muttype, alt)
                    if muttype == 'S':
                        thisseq[pos] = alt
                    elif muttype == 'D':
                        count = int(alt)
                        thisseq = thisseq[:pos] + thisseq[pos + count:]
                        realmutpos = [x[:3] + (x[3]-count, ) + x[4:] for x in realmutpos]
                    elif muttype == 'I':
                        count = len(alt)
                        thisseq = thisseq[:pos] + list(alt) + thisseq[pos:]
                        realmutpos = [x[:3] + (x[3]+count, ) + x[4:] for x in realmutpos]
                    realmutpos.append((thischr, strand, genomepos, pos, muttype, alt))
                realmutpos.sort()

                if len(realmutpos) == 0: muttotext = 'nomuts'
                else:
                    muttotext = ','.join(['..'.join([str(y) for y in x]) for x in realmutpos])
                # genometext = ','.join(genometext)

                outname = iso + '__' + muttotext
                # muttexttogenomeinfo[outname] = genometext

                if iso not in isovarcounts: isovarcounts[iso] = {}
                if outname not in isovarcounts[iso]:
                    isovarcounts[iso][outname] = [0 for s in samples]
                    isovarfa[outname] = ''.join(thisseq)

                # isovarcounts[iso][outname][sindex] = len(goodsupport[mutpos])
                for r in goodsupport[mutpos]:
                    s = int(r.split('__')[0])
                    isovarcounts[iso][outname][s] += 1


    isomutsout = open(args.output_prefix + '.isoswithvars.counts.tsv', 'w')
    isomutsfa = open(args.output_prefix + '.isoswithvars.fa', 'w')
    isomutsout.write('\t'.join(['isoname', 'varsontranscript', 'varsongenome'] + samples) + '\n')
    genomevarsout = open(args.output_prefix + '.isovars.genomicpos.bed', 'w')
    for gset in genomepostosupport:
        if genomepostosupport[gset] >= 3:
            #genomeset = (thischr, dir, str(genomePos), str(tpos), m[1], m[2])
            # genomevarsout.write('\t'.join([gset[0], str(gset[2]-1), str(gset[2]), '.'.join([str(x) for x in gset[3:]]), '0', gset[1]]) + '\n')
            genomevarsout.write('\t'.join([gset[0], str(int(gset[2])-1), gset[2], str(genomepostosupport[gset]) + '-' + '..'.join(gset[3:]), '0', gset[1]]) + '\n')

    for iso in isovarcounts:
        varperisocount = 1
        for isovar in isovarcounts[iso]:
            # genometext = muttexttogenomeinfo[isovar]
            i, m = isovar.split('__')
            counts = isovarcounts[iso][isovar]
            outseq = isovarfa[isovar]
            outname = str(varperisocount) + '-' + iso
            # isomutsout.write('\t'.join([outname, m, genometext] + [str(x) for x in counts]) + '\n')
            isomutsout.write('\t'.join([outname, m] + [str(x) for x in counts]) + '\n')
            isomutsfa.write('>' + outname + ' ' + m + '\n')
            isomutsfa.write(outseq + '\n')
            varperisocount += 1

    isomutsout.close()
    isomutsfa.close()
    genomevarsout.close()



    # print('done generating variant-aware transcript models, aligning to genome for visualization')
    # #
    # ##align isoform models to the genome
    # ###align to transcriptome with --secondary=no
    # mm2_cmd = ('minimap2', '-ax', 'splice', '-s', '40', '-t', '12', '--secondary=no', args.genome, args.output_prefix + '.isoswithvars.fa')
    # samtools_sort_cmd = ('samtools', 'sort', '-o', args.output_prefix + '.isoswithvars.bam', '-')
    # samtools_index_cmd = ('samtools', 'index', args.output_prefix + '.isoswithvars.bam')
    # pipettor.run([mm2_cmd, samtools_sort_cmd])
    # pipettor.run([samtools_index_cmd])




    # ##get productivity for transcripts without variants
    # print('getting prod')
    path = os.path.dirname(os.path.realpath(__file__)) + '/'
    prodcmd = (path + '../../bin/predictProductivity-Feb25-2', '-i', args.bedisoforms, '-o', args.output_prefix + '.isoforms.productivity', '--gtf', args.gtf, '--genome_fasta', args.genome, '--longestORF')
    pipettor.run([prodcmd])
    # print(' '.join(prodcmd))

    ##adjust productivity prediction to account for variants
    prodcmd = ('python3', path + 'predict_aaseq_withvar.py', args.output_prefix + '.isoforms.productivity.info.tsv',
                args.output_prefix + '.isoswithvars.fa', args.output_prefix + '.isoswithvars.productivity.info.tsv')
    pipettor.run([prodcmd])


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
        # if line[0] != 'isoname__muts':
        if line[0] != 'isoname':
            isoname = line[0]#'__'.join(line[0].split('__')[:-1])
            counts = [int(x) for x in line[2:]]#[2:]]
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















    ##future desire - get proper bed file with orfs from predicting productivity from variant-aware transcripts

    ###load reference gtf, get transcript pos of annotated start codons
    ##load transcript bed, get transcript pos of ptc point
    ##for productivity, now 3 categories - PRO, PST, NGO (no good orf)
    ##load cpat transcripts with no orf, output as NGO
    ###load cpat reults get transcripts to possible ORFs, filter out those <0.35

    # #cpat -g isofusionvars082824.isoswithvars.fa -o isofusionvars082824.isoswithvars.cpat.tsv -d Human_logitModel.RData -x Human_Hexamer.tsv
    # logitmodelpath = 'Human_logitModel.RData'
    # hexamerpath = 'Human_Hexamer.tsv'
    #
    # # pipettor.run([('cpat', '-g', args.output_prefix + '.isoswithvars.fa', '-o', args.output_prefix + '.isoswithvars.cpat', '-d', logitmodelpath, '-x', hexamerpath)])
    # pipettor.run([('cpat', '-g', args.output_prefix + '.flair.combined.isoforms.fa', '-o', args.output_prefix + '.isoforms.cpat', '-d', logitmodelpath, '-x', hexamerpath)])
    # print('done with cpat')
    #
    # ###get positions of ATGs on transcript
    # isoDict = {}
    # isoname = None
    # # for line in open(args.output_prefix + '.isoswithvars.fa'):
    # for line in open(args.output_prefix + '.flair.combined.isoforms.fa'):
    #     if line[0] == '>':
    #         isoname = line[1:].rstrip().split(' ')[0].upper()
    #     else:
    #         seq = line.rstrip()
    #         bedisoname = isoname#'-'.join(isoname.split('-')[1:])
    #         if bedisoname not in isoDict:
    #             isoDict[bedisoname] = Isoform(bedisoname)
    #         isoDict[bedisoname].seqvariants[isoname] = SeqVar(isoname, seq)
    #         for i in range(len(seq) - 3):
    #             codon = seq[i:i + 3]
    #             if codon == 'ATG': isoDict[bedisoname].startcodons.append(i)
    #
    # ###this whole thing is error prone if there are any insertions or deletions. Need a file that can tell me how indels affect exons
    # ##maybe from fasta file, get insertion pos and len, this is part of Isoform
    # ##figure out which exon indel is in, add/subtract from this exon len based on indel
    # for line in open(args.bedisoforms):
    #     line = line.rstrip().split('\t')
    #     isoname = line[3].upper()
    #     if isoname in isoDict:
    #         esizes = [int(x) for x in line[-2].rstrip(',').split(',')]
    #         isostart = int(line[1])
    #         estarts = [int(x) for x in line[-1].rstrip(',').split(',')]
    #         if line[-3] == '1': ptcpointont = 0
    #         else:
    #             if esizes[-2] < 55: ptcpointont = sum(esizes[:-2])
    #             else: ptcpointont = sum(esizes[:-1]) - 55
    #         exons = [(isostart + estarts[x], isostart + estarts[x] + esizes[x]) for x in range(len(estarts))]
    #         isoDict[isoname].ptcpoint = ptcpointont
    #         isoDict[isoname].strand = line[5]
    #         isoDict[isoname].chrom = line[0]
    #         isoDict[isoname].exons = exons
    #
    # annotstarts = getStarts(args.gtf)
    # bt = pybedtools.BedTool(args.bedisoforms)
    # b6 = bt.bed6()
    # st = pybedtools.BedTool(annotstarts)
    # bt_st = b6.intersect(st, s=True, split=True, wao=True)
    # for intersection in bt_st:
    #     read = intersection[3].upper()
    #     if read in isoDict:
    #         # iso,gene = read.split("_")
    #         overlap = intersection[-1]
    #         goStart = int(intersection[-6])
    #         exonCoord = (int(intersection[1]), int(intersection[2]))
    #         if overlap != "3": #check for full overlap of start codon with transcript
    #             continue
    #         exons = isoDict[read].exons
    #         exonNum = exons.index(exonCoord)
    #         exonSizes = [x[1] - x[0] for x in exons]
    #         # First get start position relative to transcript sequence. ##NOTE This might be incorrect for transcripts that have insertions or deletions, which could impact sequence length beyond what is represented in the bed file
    #         relativeStart = None
    #         if isoDict[read].strand == "+":
    #             relativeStart = (goStart - exons[exonNum][0]) + sum([x for x in exonSizes[:exonNum]])
    #         elif isoDict[read].strand == "-":
    #             relativeStart = (sum(exonSizes) - (goStart - exons[exonNum][0]) + sum([x for x in exonSizes[:exonNum]])) - 3
    #         isoDict[read].starts.add(relativeStart)
    #
    #
    # ###currently have ptc point and annotated start pos
    # ##for each transcript:
    # ##for each orf:
    # ##get whether has proper stop codon ##not sure this is necessary, I think CPAT checks for this
    # ##get whether stops before PTC point
    # ##get whether start codon matches annotation
    # ##adjust ctat score based on number of upstream ATGs
    # ##if no orfs have good score and proper start: NGO
    # ##if transcripts are all PTC:
    # ##  rank by isannot then updated ctat_score, output best as PTC
    # ##if have non-PTC transcripts:
    # #   rank by isannot then updated ctat_score, output best as PRO
    # atg_growth = 0.5
    # atg_shift = 10
    #
    # # for line in open(args.output_prefix + '.isoswithvars.cpat.no_ORF.txt'):
    # for line in open(args.output_prefix + '.isoforms.cpat.no_ORF.txt'):
    #     isoname = line.rstrip().upper()
    #     bedisoname = isoname#'-'.join(isoname.split('-')[1:])
    #     isoDict[bedisoname].seqvariants[isoname].pro = 'NGO'
    # # for line in open(args.output_prefix + '.isoswithvars.cpat.ORF_prob.tsv'):
    # for line in open(args.output_prefix + '.isoforms.cpat.ORF_prob.tsv'):
    #     line = line.rstrip().split('\t')
    #     if line[0] != 'ID':
    #         isoid, isolen, orfstart, orfend, orflen, codingprob = line[0], int(line[1]), int(line[4]), int(line[5]), int(line[6]), float(line[9])
    #         isoid = '_'.join(isoid.split('_')[:-2])
    #         bedisoname = isoid#'-'.join(isoid.split('-')[1:])
    #         orfstart -= 1
    #         if codingprob > 0.35:
    #             hasstopcodon = orfend < len(isoDict[bedisoname].seqvariants[isoid].sequence)-1
    #             isPTC = orfend < isoDict[bedisoname].ptcpoint
    #             hasannotstart = orfstart in isoDict[bedisoname].starts
    #             numpriorstarts = sum([1 if x < orfstart else 0 for x in isoDict[bedisoname].startcodons])
    #             atgscore = 1 - 1/( 1+ np.exp(-atg_growth*(numpriorstarts - atg_shift)))
    #             orfscore = codingprob * atgscore
    #             iscoding = not isPTC
    #             isoDict[bedisoname].seqvariants[isoid].orfs.append((hasstopcodon, iscoding, hasannotstart, orfscore, orfstart, orfend, orflen))
    #
    # # out = open(args.output_prefix + '.isoswithvars.orfs.tsv', 'w')
    # out = open(args.output_prefix + '.isoforms.orfs.tsv', 'w')
    # for bedisoname in isoDict:
    #     for isoid in isoDict[bedisoname].seqvariants:
    #         orfs = sorted(isoDict[bedisoname].seqvariants[isoid].orfs, reverse=True)
    #         if len(orfs) == 0: isoDict[bedisoname].seqvariants[isoid].pro = 'NGO'
    #         else:
    #             bestorf = orfs[0]
    #             if not bestorf[0]: isoDict[bedisoname].seqvariants[isoid].pro = 'NST'
    #             elif not bestorf[1]: ###if the bestORF is PTC
    #                 isoDict[bedisoname].seqvariants[isoid].pro = 'PTC'
    #             else:
    #                 isoDict[bedisoname].seqvariants[isoid].pro = 'PRO'
    #                 isoDict[bedisoname].seqvariants[isoid].orfs = orfs
    #         outline = [isoid, isoDict[bedisoname].seqvariants[isoid].pro]
    #         if isoDict[bedisoname].seqvariants[isoid].pro == 'PRO':
    #             bestorf = isoDict[bedisoname].seqvariants[isoid].orfs[0]
    #             orfstart, orfend = bestorf[4], bestorf[5]
    #             orfseq = isoDict[bedisoname].seqvariants[isoid].sequence[orfstart:orfend]
    #             aaseq = translate(orfseq)
    #             # print(isoid, aaseq[-3:], isoDict[bedisoname].seqvariants[isoid].sequence[orfend-3:orfend+6], orfend, len(isoDict[bedisoname].seqvariants[isoid].sequence))
    #             outline.extend([str(orfstart), str(orfend), aaseq])
    #         out.write('\t'.join(outline) + '\n')
    # out.close()











if __name__ == "__main__":
    getvariants()



