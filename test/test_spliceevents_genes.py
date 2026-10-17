"""
flair spliceevents: the gene each read's junction chain is counted for, and the
regions reads are corrected in
"""
import pysam
from flair import SeqRange
from flair.isoform_data import ReadRec
from flair.read_processing import add_corrected_read_to_groups
from flair.flair_spliceevents import group_juncs_by_annot_gene, combine_regions, correction_regions

# on chr1 +: gene A, one transcript of exons 100-200, 300-400 and 500-600, and gene B,
# exons 1000-1100, 1200-1300 and 1400-1500 in one transcript, and 1000-1100 and
# 1400-1500 in another
A_JUNCS = ((200, 300), (400, 500))
B_JUNCS = ((1100, 1200), (1300, 1400))
B_SKIP = ((1100, 1400),)
SJC_TO_GENE = {A_JUNCS: 'A', B_JUNCS: 'B', B_SKIP: 'B'}
JUNC_TO_GENE = {j: gene for juncs, gene in SJC_TO_GENE.items() for j in juncs}
GENE_TO_JUNCS = {'A': set(A_JUNCS), 'B': set(B_JUNCS) | set(B_SKIP)}
GENE_TO_EXONS = {'A': ((100, 200), (300, 400), (500, 600)), 'B': ((1000, 1100), (1200, 1300), (1400, 1500))}


def _genes(*reads):
    """group_juncs_by_annot_gene of reads given as (junctions, start, end): the genes'
    chains, the chains of no gene, and the number of single-exon reads"""
    sj_to_ends = {}
    for n, (juncs, start, end) in enumerate(reads):
        add_corrected_read_to_groups(ReadRec('chr1', '+', juncs, start, end, f'r{n}'), sj_to_ends)
    genes, nogene, se_reads = group_juncs_by_annot_gene(sj_to_ends, SJC_TO_GENE, JUNC_TO_GENE, GENE_TO_EXONS, GENE_TO_JUNCS)
    return {gene: sorted(chains) for gene, chains in genes.items()}, sorted(nogene), len(se_reads)


def test_an_annotated_chain_is_its_transcripts_gene():
    assert _genes((A_JUNCS, 100, 600), (B_SKIP, 1000, 1500)) == ({'A': [A_JUNCS], 'B': [B_SKIP]}, [], 0)


def test_a_new_chain_is_the_gene_of_most_of_its_junctions():
    # one of A's junctions and a new one
    chain = ((200, 300), (400, 550))
    assert _genes((chain, 100, 700))[0] == {'A': [chain]}


def test_a_chain_of_no_annotated_junction_is_the_gene_whose_exons_it_overlaps_most():
    # a new intron inside A's exons
    chain = ((330, 370),)
    assert _genes((chain, 300, 400))[0] == {'A': [chain]}


def test_a_chain_outside_all_genes_is_of_no_gene():
    chain = ((2200, 2300),)
    assert _genes((chain, 2100, 2400)) == ({}, [chain], 0)


def test_single_exon_reads_are_counted_apart():
    assert _genes(((), 100, 200), (A_JUNCS, 100, 600)) == ({'A': [A_JUNCS]}, [], 1)


def test_overlapping_and_touching_regions_are_merged():
    regions = [SeqRange('chr1', 500, 900), SeqRange('chr1', 100, 300), SeqRange('chr1', 300, 400),
               SeqRange('chr1', 200, 250), SeqRange('chr2', 100, 200)]
    assert combine_regions(regions) == [SeqRange('chr1', 100, 400), SeqRange('chr1', 500, 900), SeqRange('chr2', 100, 200)]


def _bam(path, reads):
    "a sorted, indexed bam of reads given as (start, cigar), on chr1"
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 100000}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for n, (start, cigar) in enumerate(sorted(reads)):
            read = pysam.AlignedSegment(bam.header)
            read.query_name, read.reference_id, read.reference_start = f'r{n}', 0, start
            read.cigarstring, read.mapping_quality = cigar, 60
            read.query_sequence = 'A' * read.query_length
            bam.write(read)
    pysam.index(str(path))
    return str(path)


def test_correction_regions_reach_over_the_reads_overlapping_a_region(tmp_path):
    # a read from 900 spliced to 5100-5200 crosses the region 1000-5000's end, in one sample
    bam1 = _bam(tmp_path / 's1.bam', [(900, '100M4100N100M'), (2000, '100M')])
    bam2 = _bam(tmp_path / 's2.bam', [(1500, '100M')])
    report = [SeqRange('chr1', 1000, 5000)]
    regions, reports = correction_regions(report, [bam1, bam2])
    assert regions == [SeqRange('chr1', 900, 5200)]
    assert reports == {('chr1', 900, 5200): report}


def test_regions_reaching_over_each_other_are_corrected_together(tmp_path):
    # one read crosses from the first region into the second
    bam = _bam(tmp_path / 's1.bam', [(1900, '200M800N200M')])
    report = [SeqRange('chr1', 1000, 2000), SeqRange('chr1', 2500, 3500)]
    regions, reports = correction_regions(report, [bam])
    assert regions == [SeqRange('chr1', 1000, 3500)]
    assert reports == {('chr1', 1000, 3500): report}
