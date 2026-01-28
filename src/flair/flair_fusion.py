#! /usr/bin/env python3

import sys
import argparse
import os, glob
import pipettor
import pysam
import logging
import gzip
from flair.gtf_to_bed import gtf_to_bed
from flair.bed_to_sequence import bed_to_sequence
from flair.convert_synthetic_to_genome_bed import convert_synthetic_isos
from flair import transcriptomic_chimeras
from flair import genomic_chimeras
from collections import defaultdict
from flair import FlairInputDataError


def def_value():
    return set()

def report_nofusions(outputprefix):
    logging.info('no fusions detected. Exiting')
    for file in [outputprefix + '.fusions.isoforms.bed', outputprefix + '.fusions.isoforms.fa']:
        f = open(file, 'w')
        f.close()

def detectfusions():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-g', '--genome',
                          type=str, required=True, help='FastA of reference genome')
    parser.add_argument('-f', '--gtf',
                        type=str, required=True, help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
    required.add_argument('-r', '--reads', nargs='+',
                          type=str, required=True, help='FastA/FastQ files of raw reads, can specify multiple files')
    required.add_argument('-b', '--genomechimbam',
                          type=str, required=True, help='bam file of chimeric reads from genomic alignment from flair align')
    parser.add_argument('--transcriptchimbam',
                        help='Optional: bam file of chimeric reads from transcriptomic alignment. If not provided, this will be made for you')
    parser.add_argument('-o', '--output', default='flair.fusion',
                        help='output file name base for FLAIR isoforms (default: flair.collapse)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (4)')
    parser.add_argument('--minfragmentsize', type=int, default=40,
                        help='minimum size of alignment kept, used in minimap -s. More important when doing downstream fusion detection')
    parser.add_argument('-s', '--support', type=float, default=3.0,
                        help='''minimum number of supporting reads for a fusion (3)''')
    parser.add_argument('--maxloci', type=int, default=2,
                        help='''max loci detected in fusion. Set higher for detection of 3-gene+ fusions''')
    path = os.path.dirname(os.path.realpath(__file__)) + '/'

    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        parser.error("No arguments passed. Please provide a bam file, reads, genome, and annotation file")

    args = parser.parse_args()

    if ',' in args.reads[0]:
        args.reads = args.reads[0].split(',')
    for rfile in args.reads:
        if not os.path.exists(rfile):
            raise FlairInputDataError(f'Read file path does not exist: {rfile}')

    if not os.path.exists(args.genome):
        raise FlairInputDataError(f'Genome file path does not exist: {args.genome}')
    if not os.path.exists(args.gtf):
        if not args.gtf:
            raise FlairInputDataError('Please specify annotated gtf with -f ')
        else:
            raise FlairInputDataError('GTF file path does not exist')


    # if args.annotated_fa == 'generate':
    # get transcript sequences
    args.annotated_bed = args.output + '.annotated_transcripts.bed'
    gtf_to_bed(args.annotated_bed, args.gtf, include_gene=True)
    args.annotated_fa = args.output + '.annotated_transcripts.fa'
    bed_to_sequence(query=args.output + '.annotated_transcripts.bed', genome=args.genome,
                    outfilename=args.annotated_fa)

    ####NEED TO REMEMBER THAT FUSION DETECTION RELIES ON HAVING PROPERLY STRANDED READS - need to add stranding step and/or better documentation on this

    ###Processing the gtf file so many times is really inefficient, how can we resolve this??

    ###align to transcriptome with --secondary=no
    if not args.transcriptchimbam:
        mm2_cmd = ['minimap2', '-a', '-s', str(args.minfragmentsize), '-t', str(args.threads), '--secondary=no',
                   args.annotated_fa] + args.reads
        mm2_cmd = tuple(mm2_cmd)

        # samtools; the dash at the end means STDIN
        samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', args.output + '_unfilteredtranscriptome.bam', '-')
        samtools_index_cmd = ('samtools', 'index', args.output + '_unfilteredtranscriptome.bam')
        pipettor.run([mm2_cmd, samtools_sort_cmd])
        pipettor.run([samtools_index_cmd])

        ##filter transcriptome alignment to chimeric only and remove the rest

        ##run filtering
        samfile = pysam.AlignmentFile(args.output + '_unfilteredtranscriptome.bam', "rb")
        withsup = pysam.AlignmentFile(args.output + '_transcriptomechimeric.bam', "wb", template=samfile)
        for read in samfile.fetch():
            if read.is_mapped and not read.is_secondary:
                if read.has_tag('SA'):
                    withsup.write(read)
        samfile.close()
        withsup.close()
        pysam.index(args.output + '_transcriptomechimeric.bam')

        pipettor.run([('rm', args.output + '_unfilteredtranscriptome.bam', args.output + '_unfilteredtranscriptome.bam.bai')])
        args.transcriptchimbam = args.output + '_transcriptomechimeric.bam'
        print('aligned to transcriptome')


    geneannot, genetoinfo, annot, genetoexons = {}, {}, {}, {}
    genetoinfo = {}
    for line in open(args.gtf):
        if line.startswith('#'):
            continue
        line = line.rstrip().split('\t', 8)
        chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]
        if ty in {'gene', 'exon', 'transcript'}:
            gene_id = line[8].split('gene_id "')[1].split('"')[0]
            gene_id = gene_id.replace('_', '-')
            if ty == 'gene':
                genetoinfo[gene_id] = [chrom, start, end, strand, []]

            elif ty == 'exon':
                transcript_id = line[8].split('transcript_id "')[1].split('"')[0]
                if gene_id not in genetoexons: genetoexons[gene_id] = {}
                if transcript_id not in genetoexons[gene_id]: genetoexons[gene_id][transcript_id] = []
                genetoexons[gene_id][transcript_id].append((start, end))
                if chrom not in annot: annot[chrom] = defaultdict(def_value)
                # for i in range(round(start, -1), round(end, -1), 10):
                for i in range(round(start, -2), round(end, -2) + 1, 100):
                    annot[chrom][i].add((gene_id, strand))
            elif ty == 'transcript':
                end5 = start if strand == '+' else end
                genetoinfo[gene_id][-1].append(end5)
    print('read gtf')
    intronLocs, intronToGenome = {}, {}
    for g in genetoexons:
        chrom, start, end, strand, end5s = genetoinfo[g]
        for t in genetoexons[g]:
            myexons = sorted(genetoexons[g][t])
            first, last = myexons[0], myexons[-1]
            mylocs = [[0, myexons[0][0] - 500, myexons[0][0]]]  ##add start of transcript
            runningtot = 0
            for i in range(len(myexons) - 1):
                runningtot += myexons[i][1] - myexons[i][0]  ##add size of last exon
                mylocs.append([runningtot, myexons[i][1], myexons[i + 1][0]])  ##add intron
            runningtot += myexons[-1][1] - myexons[-1][0]
            mylocs.append([runningtot, myexons[-1][1], myexons[-1][1] + 500])
            if strand == '-':
                mylocs = [[runningtot - mylocs[x][0], mylocs[x][1], mylocs[x][2]] for x in range(len(mylocs))]
            intronLocs[t] = sorted([x[0] for x in mylocs])
            intronToGenome[t] = {x[0]: (x[1], x[2]) for x in mylocs}
    print('loaded annot')

    tchim = transcriptomic_chimeras.idTranscriptomicChimeras(args.transcriptchimbam, genetoinfo, intronLocs, intronToGenome, args.support, maxloci=args.maxloci, reqdisttostart=15000)
    print('done getting fusions from transcriptome')
    # print(tchim)
    combchim = genomic_chimeras.idGenomicChimeras(args.genomechimbam, annot, geneannot, genetoinfo, args.support, maxloci=args.maxloci, reqdisttostart=15000)
    # print(combchim)
    print('done gettting fusions from genome')
    #
    for f in tchim:
        if f not in combchim:
            combchim[f] = tchim[f]
        else:
            combchim[f]['reads'] = combchim[f]['reads'] | tchim[f]['reads']
            combchim[f]['disttostart'].extend(tchim[f]['disttostart'])
            combchim[f]['qdist'].extend(tchim[f]['qdist'])
            genes = f.split('__')
            for g in genes:
                gends, tends = combchim[f][g][1:], tchim[f][g][1:]
                if gends[0] < gends[1]:
                    combchim[f][g][1] = min(gends[0], tends[0])
                    combchim[f][g][2] = max(gends[1], tends[1])
                else:
                    combchim[f][g][1] = max(gends[0], tends[0])
                    combchim[f][g][2] = min(gends[1], tends[1])

    bedout = open(args.output + '.prelimfusions.bed', 'w')
    fusionreads = set()
    for f in combchim:
        fusionreads.update(combchim[f]['reads'])
        genes = f.split('__')
        for g in genes:
            startdiststr = ','.join([str(x) for x in combchim[f]['disttostart']])
            qdiststr = ','.join([str(x) for x in combchim[f]['qdist']])
            bedname = g + '__' + '--'.join(genes) #+ '__' + startdiststr + '__' + qdiststr
            strand = '+' if combchim[f][g][1] < combchim[f][g][2] else '-'
            bedpos = [str(x) for x in sorted(combchim[f][g][1:])]

            bedline = [combchim[f][g][0]] + bedpos + [bedname, str(len(combchim[f]['reads'])), strand, startdiststr]#, qdiststr] #+ [str(x) for x in combchim[f][g][1:]] + ['255,0,0']
            bedout.write('\t'.join(bedline) + '\n')

    bedout.close()


    temp = args.reads[0].split('.')
    if temp[-1] == 'gz': temp = temp[:-1]
    freadsname = args.output + '.chimreads.' + temp[-1]


    freads = open(freadsname, 'w')
    for file in args.reads:
        last = False
        if file.split('.')[-1] == 'gz':
            readsfile = gzip.open(file, 'rt')
            file = file[:-3]
        else: readsfile = open(file)

        if file.split('.')[-1] == 'fasta' or file.split('.')[-1] == 'fa':
            for line in readsfile:
                if line[0] == '>':
                    readname = line[1:].rstrip().split()[0]
                    if readname in fusionreads: last = True
                    else: last = False
                if last: freads.write(line)
        else:
            linecount = 0
            for line in readsfile:
                if linecount % 4 == 0:
                    readname = line[1:].rstrip().split()[0]
                    if readname in fusionreads: last = True
                    else: last = False
                if last: freads.write(line)
                linecount += 1
    freads.close()
    print('done processing fusion reads')

    if os.path.getsize(args.output + '.prelimfusions.bed') == 0:
        report_nofusions(args.output)
        return


    makesynthcommand = ['python3', path + 'make_synthetic_fusion_reference.py', '-a', args.gtf, '-g', args.genome,
                        '-o', args.output, '-c', args.output + '.prelimfusions.bed']
    pipettor.run([makesynthcommand])
    if os.path.getsize(args.output + '-syntheticFusionGenome.fa') == 0:
        report_nofusions(args.output)
        return

    faidxcommand = ['samtools', 'faidx', args.output + '-syntheticFusionGenome.fa']
    mm2_cmd = ['minimap2', '-ax', 'splice', '-s', str(args.minfragmentsize), '-t', str(args.threads), '-un',
               '--secondary=no', '-G', '1000k', args.output + '-syntheticFusionGenome.fa', freadsname]
    samtools_filter_cmd = ('samtools', 'view', '-F', '2048', '-hb', '-')
    samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', args.output + '.syntheticAligned.bam', '-')
    samtools_index_cmd = ('samtools', 'index', args.output + '.syntheticAligned.bam')
    # bamtobedcmd = ('bamToBed', '-bed12', '-i', args.output + '.syntheticAligned.bam')
    bamtobedcmd = ('bedtools', 'bamtobed', '-bed12', '-i', args.output + '.syntheticAligned.bam')
    getsscommand = ['python3', path + 'synthetic_splice_sites.py', args.output + '.syntheticAligned.bed',
                        args.output + '-syntheticReferenceAnno.gtf', args.output + '.syntheticAligned.SJ.bed', args.output + '-syntheticBreakpointLoc.bed', '8', '2', args.output + '-syntheticFusionGenome.fa']#'15', '2']
    ##NOT ADDING GTF ANNOT TO correct or collapse - I think this will save time down the line
    correctcommand = ['python3', path + 'flair_correct.py', '-t', args.threads, '-q', args.output + '.syntheticAligned.bed',
                      '-f', args.output + '-syntheticReferenceAnno.gtf',
                      '--output', args.output + '.syntheticAligned.flair', '--junction_bed', args.output + '.syntheticAligned.SJ.bed', '--ss_window', '8']
    collapsecommand = ['python3', path + 'flair_collapse.py', '-t', args.threads, '-q', args.output + '.syntheticAligned.flair_all_corrected.bed',
                      '-g', args.output + '-syntheticFusionGenome.fa', #'-f', args.output + '-syntheticReferenceAnno.gtf',
                      '--output', args.output + '.syntheticAligned.flair', '-r', freadsname, '--end_window', '300', #'--stringent', '--check_splice', #'--annotation_reliant', 'generate',
                      '--generate_map', '--quality', '0', '--support', '2', '--fusion_breakpoints', args.output + '-syntheticBreakpointLoc.bed', '--allow_paralogs']


    ##currently need to run correct and collapse as subprocess because they expect specific args, need to fix this at some point I think
    # print(makesynthcommand)

    pipettor.run([faidxcommand])
    print('synth genome made')
    pipettor.run([mm2_cmd, samtools_filter_cmd, samtools_sort_cmd])
    print('done aligning to synthetic fusion genome')
    pipettor.run([samtools_index_cmd])
    pipettor.run([bamtobedcmd], stdout=args.output + '.syntheticAligned.bed')
    pipettor.run([getsscommand])
    print('done getting ss')
    pipettor.run([correctcommand])
    print('done with correct')
    pipettor.run([collapsecommand])

    ##clean up isoform/gene names for args.output + '.combined.isoform.read.map.txt', '.isoforms.bed', '.isoforms.fa'
    oldnametonewname = {}
    out = open(args.output + '.syntheticAligned.isoforms.bed', 'w')
    c = 0
    for line in open(args.output + '.syntheticAligned.flair.isoforms.bed'):
        c +=1
        line = line.rstrip().split('\t')
        newname = 'fusioniso' + str(c) + '_' + line[0]
        oldnametonewname[line[3]] = newname
        line[3] = newname
        out.write('\t'.join(line) + '\n')
    out.close()
    out = open(args.output + '.syntheticAligned.isoform.read.map.txt', 'w')
    for line in open(args.output + '.syntheticAligned.flair.isoform.read.map.txt'):
        line = line.split('\t', 1)
        line[0] = oldnametonewname[line[0]]
        out.write('\t'.join(line))
    out.close()

    maxpromiscuity = 4
    convert_synthetic_isos(args.gtf, args.output + '.syntheticAligned.isoforms.bed',
                      args.output + '.syntheticAligned.isoform.read.map.txt', freadsname,
                      args.output + '-syntheticBreakpointLoc.bed', args.output + '.fusions.isoforms.bed', os.path.realpath(__file__).split('flair_fusion')[0] + 'dgd_Hsa_all_v71.tsv', maxpromiscuity)
    goodisos = set()
    for line in open(args.output + '.fusions.isoforms.bed'):
        line = line.rstrip().split('\t')
        goodisos.add('_'.join(line[3].split('_')[1:]))

    out = open(args.output + '.fusions.isoforms.fa', 'w')
    good = False
    for line in open(args.output + '.syntheticAligned.flair.isoforms.fa'):
        if line[0] == '>':
            oldname = line[1:].rstrip()
            if oldnametonewname[oldname] in goodisos: good = True
            else: good = False
            if good: out.write('>' + oldnametonewname[oldname] + '\n')
        elif good: out.write(line)
    out.close()

    os.rename(args.output + '.syntheticAligned.isoform.read.map.txt', args.output + '.fusion.isoform.read.map.txt')

    #removing extra FLAIR files
    for filename in glob.glob(args.output + '.syntheticAligned.flair*'):
        os.remove(filename)






if __name__ == '__main__':
    detectfusions()
