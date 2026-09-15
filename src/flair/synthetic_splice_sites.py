"""Call splice junctions from reads aligned to a synthetic fusion reference.

Splice sites that are not close to an annotated junction are clustered within
sjwiggle and snapped to the best-supported position in their cluster.
"""
import argparse
from collections import Counter
import pysam
from flair.gtf_io import gtf_record_parser, GtfAttrsSet
from flair.pycbio.hgdata.bed import BedReader

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignedbed', help="BED12 of reads aligned to the synthetic reference")
    parser.add_argument('referencegtf', help="annotation on the synthetic reference")
    parser.add_argument('outfile', help="BED6 of called splice junctions")
    parser.add_argument('refbp', help="BED3 of the fusion breakpoint on each synthetic chromosome")
    parser.add_argument('sjwiggle', type=int, help="splice sites this far apart are one site")
    parser.add_argument('readcov', type=int, help="reads needed to support a splice site")
    parser.add_argument('refgenome', help="FASTA of the synthetic reference")
    return parser.parse_args()


def read_annot_splice_juncs(referencegtffile):
    "annotated exon bounds per chromosome, 1-based to match the read introns"
    juncs = {}
    for rec in gtf_record_parser(referencegtffile, include_features={'exon'}, attrs=GtfAttrsSet.FLAIR):
        juncs.setdefault(rec.chrom, set()).add((rec.start + 1, rec.end))
    return juncs


def read_breakpoints(refbpfile):
    "fusion breakpoint position of each synthetic chromosome"
    return {bed.chrom: bed.chromStart for bed in BedReader(refbpfile, numStdCols=3)}


def _bed_introns(bed):
    "introns of a read as (chrom, start, end), with end 1-based"
    return [(bed.chrom, bed.blocks[i].end, bed.blocks[i + 1].start + 1)
            for i in range(len(bed.blocks) - 1)]


def count_read_introns(alignedbedfile):
    "read count and the strands those reads were on, for each intron"
    introns = {}
    for bed in BedReader(alignedbedfile, fixScores=True):
        for intron in _bed_introns(bed):
            count, strands = introns.get(intron, (0, Counter()))
            strands[bed.strand] += 1
            introns[intron] = (count + 1, strands)
    return introns


def find_motif_strand(genome, chrom, start, end):
    "strand implied by the GT/AG or CT/AC motifs, None when neither matches"
    donor = genome.fetch(chrom, start, start + 2)
    acceptor = genome.fetch(chrom, end - 3, end - 1)
    if (donor == 'GT') and (acceptor == 'AG'):
        return '+'
    elif (donor == 'CT') and (acceptor == 'AC'):
        return '-'
    else:
        return None


def _intron_strand(motif_strand, read_strands):
    "the motif strand when there is one, otherwise the commonest read strand"
    if motif_strand is not None:
        return motif_strand
    return read_strands.most_common(1)[0][0]


def near_annot_junc(annot_juncs, chrom, start, end, sjwiggle):
    "is the intron within sjwiggle of an annotated junction at both ends"
    for sj in annot_juncs.get(chrom, ()):
        if (abs(start - sj[0]) <= sjwiggle) and (abs(end - sj[1]) <= sjwiggle):
            return True
    return False


def _novel_intron_strand(genome, breakpoints, annot_juncs, intron, strands, sjwiggle):
    """Strand of an intron worth calling: one with a splice motif or crossing the
    breakpoint, and not already annotated.  None when the intron is not worth calling."""
    chrom, start, end = intron
    spans_breakpoint = start <= breakpoints[chrom] <= end
    motif_strand = find_motif_strand(genome, chrom, start, end)
    if not (spans_breakpoint or (motif_strand is not None)):
        return None
    elif near_annot_junc(annot_juncs, chrom, start, end, sjwiggle):
        return None
    else:
        return _intron_strand(motif_strand, strands)


def collect_novel_splice_sites(read_introns, annot_juncs, breakpoints, genome, sjwiggle):
    """Unannotated splice-site positions per chromosome, each repeated once per
    supporting read, along with the introns they came from."""
    chr_to_novel_ss = {}
    novel_introns = {}
    for intron, (readcount, strands) in read_introns.items():
        strand = _novel_intron_strand(genome, breakpoints, annot_juncs, intron, strands, sjwiggle)
        if strand is not None:
            positions = chr_to_novel_ss.setdefault(intron[0], [])
            positions.extend([intron[1]] * readcount)
            positions.extend([intron[2]] * readcount)
            novel_introns[intron] = (readcount, strand)
    return chr_to_novel_ss, novel_introns


def _cluster(positions, sjwiggle):
    "group sorted positions that are within sjwiggle of their predecessor"
    group = []
    prev = None
    for pos in positions:
        if (prev is not None) and (pos - prev > sjwiggle):
            yield group
            group = []
        group.append(pos)
        prev = pos
    if len(group) > 0:
        yield group


def _consensus_positions(group, sjwiggle, readcov):
    "the best-supported positions of a cluster, kept more than sjwiggle apart"
    chosen = []
    for pos, count in Counter(group).most_common():
        if (count >= readcov) and (count >= len(group) / 10):
            if not any(fp - sjwiggle <= pos <= fp + sjwiggle for fp in chosen):
                chosen.append(pos)
    return chosen


def consensus_splice_sites(chr_to_novel_ss, sjwiggle, readcov):
    "novel splice-site positions collapsed to the consensus positions of their cluster"
    chr_to_good_ss = {}
    for chrom, positions in chr_to_novel_ss.items():
        chr_to_good_ss[chrom] = set()
        for group in _cluster(sorted(positions), sjwiggle):
            if len(group) >= readcov:
                chr_to_good_ss[chrom].update(_consensus_positions(group, sjwiggle, readcov))
    return chr_to_good_ss


def _snap(pos, good_ss, sjwiggle):
    "nearest consensus splice site within sjwiggle, or pos when there is none"
    near = [sj for sj in good_ss if abs(pos - sj) <= sjwiggle]
    return min(near, key=lambda sj: abs(pos - sj)) if len(near) > 0 else pos


def snap_introns(novel_introns, chr_to_good_ss, sjwiggle):
    "read support per junction after both ends are snapped to the consensus sites"
    support = {}
    for (chrom, start, end), (readcount, strand) in novel_introns.items():
        good_ss = chr_to_good_ss.get(chrom, ())
        junc = (chrom, _snap(start, good_ss, sjwiggle), _snap(end, good_ss, sjwiggle) - 1, strand)
        support[junc] = support.get(junc, 0) + readcount
    return support


def write_splice_juncs(outfilename, splice_junc_support, min_support=2):
    "junctions with at least min_support reads, as BED6"
    with open(outfilename, 'w') as out:
        for (chrom, start, end, strand), support in splice_junc_support.items():
            if support >= min_support:
                out.write('\t'.join([chrom, str(start), str(end), '.', str(support), strand]) + '\n')


def synthetic_splice_sites(alignedbed, referencegtf, outfile, refbp, sjwiggle, readcov, refgenome):
    genome = pysam.FastaFile(refgenome)
    annot_juncs = read_annot_splice_juncs(referencegtf)
    breakpoints = read_breakpoints(refbp)
    read_introns = count_read_introns(alignedbed)
    chr_to_novel_ss, novel_introns = collect_novel_splice_sites(read_introns, annot_juncs, breakpoints,
                                                                genome, sjwiggle)
    chr_to_good_ss = consensus_splice_sites(chr_to_novel_ss, sjwiggle, readcov)
    write_splice_juncs(outfile, snap_introns(novel_introns, chr_to_good_ss, sjwiggle))


def main():
    args = parse_args()
    synthetic_splice_sites(args.alignedbed, args.referencegtf, args.outfile, args.refbp,
                           args.sjwiggle, args.readcov, args.refgenome)


if __name__ == '__main__':
    main()
