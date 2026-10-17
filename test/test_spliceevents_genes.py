"""
flair spliceevents: the gene each read's junction chain is counted for, and the
regions reads are corrected in
"""
import pysam
from flair import SeqRange
from flair.annotation_data import annot_data_from_gtf
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.isoform_data import Isoform, Junc
from flair.flair_spliceevents import (annotated_genes, novel_genes_of_chains, gene_junctions, group_chains, combine_regions,
                                      correction_regions, ReadEnds, ReadChain, read_read_chains, write_read_chains)

# on chr1, 0-based starts: genes A and B on +, B with a transcript skipping its middle
# exon, and D and E on -
GTF = [('A', 'TA', '+', [(100, 200), (300, 400), (500, 600)]),
       ('B', 'TB', '+', [(1000, 1100), (1200, 1300), (1400, 1500)]),
       ('B', 'TB2', '+', [(1000, 1100), (1400, 1500)]),
       ('D', 'TD', '-', [(5000, 5100), (5200, 5300)]),
       ('E', 'TE', '-', [(6000, 6100), (6200, 6300)])]
A_JUNCS = ((200, 300), (400, 500))
B_SKIP = ((1100, 1400),)


def _annots(tmp_path):
    gtf = tmp_path / 'annotation.gtf'
    with open(gtf, 'w') as fh:
        for gene, transcript, strand, exons in GTF:
            attrs = f'gene_id "{gene}"; transcript_id "{transcript}"; gene_name "{gene}";'
            fh.write(f'chr1\ttest\ttranscript\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t{strand}\t.\t{attrs}\n')
            for start, end in exons:
                fh.write(f'chr1\ttest\texon\t{start + 1}\t{end}\t.\t{strand}\t.\t{attrs}\n')
    gtf_data = gtf_data_parser(str(gtf), attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)
    return annot_data_from_gtf(gtf_data, SeqRange('chr1', 0, 10 ** 6))


def _genes(tmp_path, juncs, strand='+'):
    return annotated_genes(Isoform('chr1', strand, tuple(Junc(*j) for j in juncs)), _annots(tmp_path))


def test_an_annotated_chain_is_its_transcripts_gene(tmp_path):
    assert _genes(tmp_path, A_JUNCS) == (('A',), '+')
    assert _genes(tmp_path, B_SKIP) == (('B',), '+')


def test_a_new_chain_sharing_junctions_is_their_gene(tmp_path):
    # one of A's junctions and a new one
    assert _genes(tmp_path, ((200, 300), (400, 550))) == (('A',), '+')


def test_a_readthrough_chain_is_its_genes_in_the_order_it_reads_them(tmp_path):
    assert _genes(tmp_path, ((200, 300), (400, 500), (600, 1000), (1100, 1200))) == (('A', 'B'), '+')
    # on -, E, the higher, is read first
    assert _genes(tmp_path, ((5100, 5200), (5300, 6000), (6100, 6200)), '-') == (('E', 'D'), '-')


def test_a_chain_of_unknown_strand_gets_its_genes_strand(tmp_path):
    assert _genes(tmp_path, A_JUNCS, '.') == (('A',), '+')
    # and one of no gene has none, to be dropped
    assert _genes(tmp_path, ((2200, 2300),), '.') == ((), None)


def test_a_chain_outside_all_genes_is_of_no_annotated_gene(tmp_path):
    assert _genes(tmp_path, ((2200, 2300),)) == ((), '+')


def test_chains_of_no_annotated_gene_get_the_same_novel_gene_in_every_sample(tmp_path):
    # chains sharing the splice site at 2300, one in each sample
    s1 = [ReadChain((), ((2200, 2300),), '+', [ReadEnds(2100, 2400)] * 3), ReadChain(('A',), A_JUNCS, '+', [ReadEnds(100, 600)])]
    s2 = [ReadChain((), ((2300, 2400),), '+', [ReadEnds(2250, 2500)] * 2)]
    novel_gene = novel_genes_of_chains('chr1', [s1, s2])
    novel = 'chr1:2100-2500:+'
    assert novel_gene == {('+', ((2200, 2300),)): novel, ('+', ((2300, 2400),)): novel}
    # a novel gene's junctions are its reads', an annotated one's the annotation's
    assert gene_junctions([s1, s2], novel_gene, _annots(tmp_path)) == {novel: {(2200, 2300), (2300, 2400)}, 'A': set(A_JUNCS)}


def test_a_readthroughs_reads_count_with_its_5_gene_for_starts_and_its_3_gene_for_ends():
    readthrough = ReadChain(('A', 'B'), ((200, 300), (600, 1000)), '+', [ReadEnds(150, 1300)])
    a = ReadChain(('A',), A_JUNCS, '+', [ReadEnds(100, 600)])
    views, strands, cluster_genes, _ = group_chains('chr1', ['s1'], [[readthrough, a]], {})
    assert views['five'] == [{'A': {readthrough.juncs: readthrough.ends, A_JUNCS: a.ends}}]
    assert views['three'] == [{'B': {readthrough.juncs: readthrough.ends}, 'A': {A_JUNCS: a.ends}}]
    # the two chains share A's donor at 200: one cluster, of both genes
    assert list(cluster_genes.values()) == ['A,B']


def test_a_gene_nested_in_anothers_intron_is_a_cluster_of_its_own():
    outer = ReadChain(('O',), ((1000, 5000),), '+', [ReadEnds(900, 5100)])
    nested = ReadChain(('N',), ((2000, 2100), (2200, 2300)), '+', [ReadEnds(1900, 2400)])
    _, _, cluster_genes, _ = group_chains('chr1', ['s1'], [[outer, nested]], {})
    assert sorted(cluster_genes.values()) == ['N', 'O']


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


def test_read_chains_written_and_read_back_are_the_same_reads(tmp_path):
    chains = [ReadChain(('A',), A_JUNCS, '+', [ReadEnds(100, 600), ReadEnds(120, 580)]),
              ReadChain(('A', 'B'), ((200, 300), (600, 1000)), '+', [ReadEnds(150, 1300)]),
              ReadChain((), ((2200, 2300),), '-', [ReadEnds(2100, 2400)])]
    path = str(tmp_path / 's1.read_chains.tsv')
    write_read_chains(path, 'chr1', chains)
    assert open(path).read().splitlines()[1] == 'A,B\tchr1\t200-300,600-1000\t+\t150-1300'
    assert read_read_chains(path) == chains
