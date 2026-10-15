import argparse
import pysam
from flair import FlairInputDataError
from flair.isoform_data import get_reverse_complement
from flair.gtf_io import gtf_record_parser, gtf_write_row, GtfAttrsSet
from flair.pycbio.hgdata.bed import BedReader

# Makes a synthetic genome of fusion contigs, each the loci of a preliminary fusion
# joined in fusion order, with an annotation of the fusion genes' transcripts on them
# and the breakpoint positions between the loci.


def read_breakpoints(chimbp):
    """From the preliminary fusions bed: fusion -> its loci (gene, chrom, 5' bound, 3' bound)
    in order, gene -> [chrom, start, end, strand] over all its loci, and gene -> an empty
    transcript dict"""
    allBP, fgenes, transcripts = {}, {}, {}
    # ['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'breakpointCoord', 'outerEdgeCoord', 'readSupport']
    for bed in BedReader(chimbp, numStdCols=6, fixScores=True):
        chrom, name, strand = bed.chrom, bed.name, bed.strand
        start, end = max(0, bed.chromStart), max(0, bed.chromEnd)
        gene, fusion = name.split('__')
        fusion = tuple(fusion.split('--'))
        if strand == '+':
            allBP.setdefault(fusion, []).append((gene, chrom, start, end))
        else:
            allBP.setdefault(fusion, []).append((gene, chrom, end, start))
        transcripts.setdefault(gene, {})
        if gene not in fgenes:
            fgenes[gene] = [chrom, start, end, strand]
        else:
            fgenes[gene][1] = min(start, fgenes[gene][1])
            fgenes[gene][2] = max(end, fgenes[gene][2])
    return allBP, fgenes, transcripts


class ChromSizes:
    "chrom sizes of the genome, an error for a sequence the genome doesn't have"
    def __init__(self, genome, genome_path, anno_path):
        self.sizes = dict(zip(genome.references, genome.lengths))
        self.genome_path, self.anno_path = genome_path, anno_path

    def get(self, chrom):
        size = self.sizes.get(chrom)
        if size is None:
            raise FlairInputDataError(f"sequence '{chrom}' in annotation {self.anno_path} is not in genome {self.genome_path}; "
                                      "supply an annotation and genome from the same assembly, using the same sequence names")
        return size


def read_transcripts(anno, fgenes, transcripts, chrom_sizes):
    """Fill in the fusion genes' extents from the annotation, padded by 500 bases, and their
    transcripts' exons and start codons, in genomic order"""
    for rec in gtf_record_parser(anno, include_features={'gene', 'exon', 'start_codon'}, attrs=GtfAttrsSet.ALL):
        genename = rec.gene_id
        if genename not in fgenes:
            continue
        if rec.feature == 'gene':
            # learned that can't assume that transcript appears in anno only once - two diff ENSG can have same hugo name
            # the 500 base padding is held within the sequence; a gene near either end has less room than that
            fgenes[genename] = (rec.chrom,
                                min([max(0, rec.start - 500), fgenes[genename][1]]),
                                max([min(chrom_sizes.get(rec.chrom), rec.end + 500), fgenes[genename][2]]),
                                rec.strand)
        elif rec.feature in ('exon', 'start_codon'):
            transcripts[genename].setdefault(rec.transcript_id, []).append((rec.start, rec.end, rec.feature))
    # in genomic order whatever the order of the GTF: minus strand exons were put in
    # genomic order by reversing their file order, which only works for an annotation
    # listing them in transcript order, as GENCODE does.  A position-sorted GTF gave
    # transcripts ending before they started, and flair transcriptome failed on them
    for genename in transcripts:
        for tname in transcripts[genename]:
            transcripts[genename][tname].sort()


def locus_bounds(order, numloci, gene, leftbound, rightbound, fgenes, chrom_size):
    """A locus's 5' and 3' bounds on its chromosome, the 5' gene's taken to the gene's 5'
    end and the 3' gene's to its 3' end"""
    # the breakpoint BED is padded by 1000 and clamped only at the low end, so a
    # bound can run past the chromosome.  genome.fetch then returns a short
    # sequence while startLoc advances by the nominal distance, which shifts every
    # later locus and the breakpoint offsets against the FASTA actually written
    leftbound = max(0, min(leftbound, chrom_size))
    rightbound = max(0, min(rightbound, chrom_size))
    if order == 0:  # 5' gene, correct end to 5' end of gene
        leftbound = fgenes[gene][1] if leftbound < rightbound else fgenes[gene][2]
    elif order == numloci - 1:
        rightbound = fgenes[gene][2] if leftbound < rightbound else fgenes[gene][1]
    return leftbound, rightbound


def synthetic_exons(exons, order, leftbound, rightbound, startLoc):
    """A transcript's exons in a locus, in fusion contig coordinates; start codons only for
    the 5' gene"""
    chunk = []
    if leftbound < rightbound:  # positive strand transcript
        for start, end, feature in exons:
            if start >= leftbound and end <= rightbound and (order == 0 or feature == 'exon'):
                chunk.append(((start - leftbound) + startLoc, (end - leftbound) + startLoc, feature))
    else:  # negative strand transcript
        for start, end, feature in reversed(exons):
            if start >= rightbound and end <= leftbound and (order == 0 or feature == 'exon'):
                chunk.append(((leftbound - end) + startLoc, (leftbound - start) + startLoc, feature))
    return chunk


def build_fusion(loci, fgenes, transcripts, genome, chrom_sizes):
    """A fusion contig's locus labels, sequences, start positions and, per locus, transcript
    -> its exons there"""
    labels, sequence, allstartloc, allisochunks = [], [], [], []
    startLoc = 0
    for order, (gene, thisChr, leftbound, rightbound) in enumerate(loci):
        leftbound, rightbound = locus_bounds(order, len(loci), gene, leftbound, rightbound, fgenes,
                                             chrom_sizes.get(thisChr))
        if leftbound < rightbound:
            sequence.append(genome.fetch(thisChr, leftbound, rightbound))
        else:
            sequence.append(get_reverse_complement(genome.fetch(thisChr, rightbound, leftbound)))
        labels.append('..'.join([str(x) for x in [gene, thisChr, leftbound, rightbound]]))
        allisochunks.append({tname: synthetic_exons(exons, order, leftbound, rightbound, startLoc)
                             for tname, exons in transcripts[gene].items()})
        allstartloc.append(startLoc)
        startLoc = startLoc + abs(rightbound - leftbound)
    return labels, sequence, allstartloc, allisochunks


def multi_exon_isoforms(allisochunks):
    "[transcript, exons] of each locus's distinct transcripts with more than one exon there"
    finalisochunks = []
    for isochunks in allisochunks:
        seen = []
        for iso, chunk in isochunks.items():
            if chunk != [] and len([x for x in chunk if x[2] == 'exon']) > 1 and chunk not in seen:
                finalisochunks.append([iso, chunk])
                seen.append(chunk)
    return finalisochunks


def write_fusion(outs, fusion, labels, sequence, allstartloc, finalisochunks, isocount):
    "write a fusion contig, its breakpoints and annotation; returns the next transcript number"
    out, annoOut, bpOut = outs
    fusionname = '--'.join(labels)
    geneid = '--'.join(fusion)
    fusionchrname = '--'.join(['..'.join(x.split('*')) for x in fusion])
    out.write('>' + fusionchrname + '\n')
    out.write(''.join(sequence) + '\n')
    for s in range(1, len(fusion)):
        bpOut.write('\t'.join([fusionchrname, str(allstartloc[s]), str(allstartloc[s]), 'breakpoint-' + str(s) + '--' + fusionname]) + '\n')
    gtf_write_row(annoOut, fusionchrname, 'SYNTHFUSION', 'gene', 0, len(''.join(sequence)), None, '+', None,
                  gene_id=geneid)
    for iso, isoexons in finalisochunks:
        transcriptid = str(isocount) + '-' + iso
        isocount += 1
        gtf_write_row(annoOut, fusionchrname, 'SYNTHFUSION', 'transcript', isoexons[0][0], isoexons[-1][1], None, '+', None,
                      gene_id=geneid, transcript_id=transcriptid)
        for exon in isoexons:
            gtf_write_row(annoOut, fusionchrname, 'SYNTHFUSION', exon[2], exon[0], exon[1], None, '+', None,
                          gene_id=geneid, transcript_id=transcriptid)
    return isocount


def main():
    parser = argparse.ArgumentParser(description='make synthetic fusion reference')
    parser.add_argument('-c', '--chimbp', action='store', help='bed file of fusion breakpoints')
    parser.add_argument('-g', '--genome', action='store', dest='g', help='path to genome')
    parser.add_argument('-a', '--anno', action='store', dest='a', help='path to anno.gtf')
    parser.add_argument('-o', '--output', action='store', help='output file prefix')
    args = parser.parse_args()
    prefix = args.output

    allBP, fgenes, transcripts = read_breakpoints(args.chimbp)
    genome = pysam.FastaFile(args.g)
    chrom_sizes = ChromSizes(genome, args.g, args.a)
    read_transcripts(args.a, fgenes, transcripts, chrom_sizes)

    with open(prefix + '-syntheticFusionGenome.fa', 'w') as out, \
            open(prefix + '-syntheticReferenceAnno.gtf', 'w') as annoOut, \
            open(prefix + '-syntheticBreakpointLoc.bed', 'w') as bpOut:
        isocount = 1
        for fusion in allBP:
            labels, sequence, allstartloc, allisochunks = build_fusion(allBP[fusion], fgenes, transcripts, genome, chrom_sizes)
            isocount = write_fusion((out, annoOut, bpOut), fusion, labels, sequence, allstartloc,
                                    multi_exon_isoforms(allisochunks), isocount)


if __name__ == '__main__':
    main()
