"""
flair spliceevents: the gene each read's junction chain is counted for, and the
regions reads are corrected in
"""
import pysam
from flair import SeqRange
from flair.annotation_data import annot_data_from_gtf
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.isoform_data import Isoform, Junc
from flair.spliceevents_events import SplicingEvent, extract_terminal_info
from flair.spliceevents_genes import novel_genes_of_chains, gene_junctions, group_chains
from flair.spliceevents_reads import (annotated_genes, ReadEnds, ReadChain, read_read_chains, write_read_chains,
                                      annotated_terminal_exons, annotated_terminal_splice_sites, AnnotatedFragmentExons,
                                      AnnotatedTerminalExon, fragment_chains)
from flair.spliceevents_regions import (combine_regions, correction_regions, partition_regions, gene_clusters,
                                        intergenic_cuts, event_in_regions)

# on chr1, 0-based starts: genes A and B on +, B with a transcript skipping its middle
# exon, and D and E on -
GTF = [('A', 'TA', '+', [(100, 200), (300, 400), (500, 600)]),
       ('B', 'TB', '+', [(1000, 1100), (1200, 1300), (1400, 1500)]),
       ('B', 'TB2', '+', [(1000, 1100), (1400, 1500)]),
       ('D', 'TD', '-', [(5000, 5100), (5200, 5300)]),
       ('E', 'TE', '-', [(6000, 6100), (6200, 6300)])]
A_JUNCS = ((200, 300), (400, 500))
B_SKIP = ((1100, 1400),)


def _gtf_data(tmp_path, gtf_transcripts=GTF, chrom='chr1'):
    gtf = tmp_path / 'annotation.gtf'
    with open(gtf, 'w') as fh:
        for gene, transcript, strand, exons in gtf_transcripts:
            attrs = f'gene_id "{gene}"; transcript_id "{transcript}"; gene_name "{gene}";'
            fh.write(f'{chrom}\ttest\ttranscript\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t{strand}\t.\t{attrs}\n')
            for start, end in exons:
                fh.write(f'{chrom}\ttest\texon\t{start + 1}\t{end}\t.\t{strand}\t.\t{attrs}\n')
    return gtf_data_parser(str(gtf), attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)


def _annots(tmp_path, gtf_transcripts=GTF):
    return annot_data_from_gtf(_gtf_data(tmp_path, gtf_transcripts), SeqRange('chr1', 0, 10 ** 6))


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


def test_an_annotated_gene_of_one_exon_has_its_reads_junctions(tmp_path):
    # as HEXIM1, whose reads splice inside its one exon: without junctions it had
    # no span, and its starts and ends were reported only in a region from 0
    annots = _annots(tmp_path, GTF + [('S', 'TS', '+', [(3000, 4000)])])
    chains = [ReadChain(('S',), ((3200, 3500),), '+', [ReadEnds(3000, 4000)] * 3)]
    assert gene_junctions([chains], {}, annots) == {'S': {(3200, 3500)}}


def test_a_readthroughs_reads_count_with_its_5_gene_for_starts_and_its_3_gene_for_ends():
    readthrough = ReadChain(('A', 'B'), ((200, 300), (600, 1000)), '+', [ReadEnds(150, 1300)])
    a = ReadChain(('A',), A_JUNCS, '+', [ReadEnds(100, 600)])
    views, strands, chain_genes, _ = group_chains('chr1', ['s1'], [[readthrough, a]], {})
    assert views['five'] == [{'A': {readthrough.juncs: readthrough.ends, A_JUNCS: a.ends}}]
    assert views['three'] == [{'B': {readthrough.juncs: readthrough.ends}, 'A': {A_JUNCS: a.ends}}]
    # the two chains share A's donor at 200: one cluster
    assert [sorted(chains) for chains in views['cluster'][0].values()] == [sorted([A_JUNCS, readthrough.juncs])]
    assert chain_genes[('+', readthrough.juncs)] == ('A', 'B')


def test_a_gene_nested_in_anothers_intron_is_a_cluster_of_its_own():
    outer = ReadChain(('O',), ((1000, 5000),), '+', [ReadEnds(900, 5100)])
    nested = ReadChain(('N',), ((2000, 2100), (2200, 2300)), '+', [ReadEnds(1900, 2400)])
    views, _, _, _ = group_chains('chr1', ['s1'], [[outer, nested]], {})
    assert sorted(list(chains) for chains in views['cluster'][0].values()) == [[outer.juncs], [nested.juncs]]


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


def test_regions_are_cut_in_the_gaps_between_genes(tmp_path):
    transcripts = _gtf_data(tmp_path).transcripts
    region = SeqRange('chr1', 0, 10000)
    # A 100-600 (3 exons), B 1000-1500 (5, of two transcripts), D 5000-5300 (2), E 6000-6300 (2)
    assert gene_clusters(region, transcripts) == [[100, 600, 3], [1000, 1500, 5], [5000, 5300, 2], [6000, 6300, 2]]
    assert intergenic_cuts(region, transcripts) == [800, 3250, 5650]
    # pieces of at least 5 exons, but the last
    assert intergenic_cuts(region, transcripts, min_exons=5) == [3250]
    # genes cut by the region are clipped to it
    assert gene_clusters(SeqRange('chr1', 1200, 5100), transcripts) == [[1200, 1500, 5], [5000, 5100, 2]]


def test_pieces_own_the_reads_starting_in_them(tmp_path):
    # a read from 1400 to 3600 crosses the cut at 3250, and is only the second piece's
    bam = _bam(tmp_path / 's1.bam', [(200, '300M'), (1100, '100M'), (1400, '100M2000N100M'), (6000, '100M')])
    region = SeqRange('chr1', 0, 10000)
    pieces = partition_regions([region], [bam], {('chr1', 0, 10000): [800, 3250, 5650]})
    # each corrected over its piece and its reads; the third, 3250-5650, has none
    assert pieces == [(SeqRange('chr1', 0, 800), SeqRange('chr1', 0, 800), region, 1),
                      (SeqRange('chr1', 800, 3600), SeqRange('chr1', 800, 3250), region, 2),
                      (SeqRange('chr1', 5650, 10000), SeqRange('chr1', 5650, 10000), region, 1)]


def test_the_first_piece_owns_the_reads_starting_before_the_region(tmp_path):
    bam = _bam(tmp_path / 's1.bam', [(900, '200M')])
    region = SeqRange('chr1', 1000, 5000)
    pieces = partition_regions([region], [bam], {('chr1', 1000, 5000): [3000]})
    assert pieces == [(SeqRange('chr1', 900, 3000), SeqRange('chr1', 0, 3000), region, 1)]


def test_a_cut_regions_events_are_reported_by_the_piece_they_start_in():
    event = SplicingEvent('ir|chr1:3001-3500(+)', 'ir', 'chr1', '+', {}, 'A')
    event.add('spliced', {}, [(3000, 3500)], ())
    region = [SeqRange('chr1', 0, 10000)]
    assert event_in_regions(event, region, None)
    assert event_in_regions(event, region, None, SeqRange('chr1', 800, 3250))
    assert not event_in_regions(event, region, None, SeqRange('chr1', 3250, 5650))


def test_read_chains_written_and_read_back_are_the_same_reads(tmp_path):
    chains = [ReadChain(('A',), A_JUNCS, '+', [ReadEnds(100, 600), ReadEnds(120, 580)]),
              ReadChain(('A', 'B'), ((200, 300), (600, 1000)), '+', [ReadEnds(150, 1300)]),
              ReadChain((), ((2200, 2300),), '-', [ReadEnds(2100, 2400)])]
    path = str(tmp_path / 's1.read_chains.tsv')
    write_read_chains(path, 'chr1', chains)
    assert open(path).read().splitlines()[1] == 'A,B\tchr1\t200-300,600-1000\t+\t150-1300\t'
    assert read_read_chains(path) == chains


def test_first_and_last_exons_and_their_splice_sites_are_by_chromosome(tmp_path):
    terminal = annotated_terminal_exons(_gtf_data(tmp_path).transcripts)
    assert set(terminal) == {'chr1'}
    assert AnnotatedTerminalExon('first', '+', 100, 200, 'A') in terminal['chr1']
    assert AnnotatedTerminalExon('last', '+', 500, 600, 'A') in terminal['chr1']
    # on -, the first exon is the higher
    assert AnnotatedTerminalExon('first', '-', 5200, 5300, 'D') in terminal['chr1']
    assert AnnotatedTerminalExon('last', '-', 5000, 5100, 'D') in terminal['chr1']
    afe_ss, ale_ss = annotated_terminal_splice_sites(terminal)
    assert afe_ss['chr1']['+'] == [200, 1100] and afe_ss['chr1']['-'] == [5200, 6200]
    assert ale_ss['chr1']['+'] == [500, 1400] and ale_ss['chr1']['-'] == [5100, 6100]
    # another chromosome's sites are its own
    other = annotated_terminal_exons(_gtf_data(tmp_path, chrom='chr2').transcripts)
    assert annotated_terminal_splice_sites(other)[0].keys() == {'chr2'}


# A's last exon is 500-600 and first 100-200 on +; D's last 5000-5100 and first
# 5200-5300 on -; and an antisense gene, AS, overlapping B's last exon, 1400-1500
FRAGMENT_GTF = GTF + [('AS', 'TAS', '-', [(1450, 1550), (1700, 1800)]), ('S', 'TS', '+', [(8000, 9000)])]


def _fragment_exons(tmp_path):
    gtf_data = _gtf_data(tmp_path, FRAGMENT_GTF)
    annots = annot_data_from_gtf(gtf_data, SeqRange('chr1', 0, 10 ** 6))
    return AnnotatedFragmentExons(annotated_terminal_exons(gtf_data.transcripts)['chr1'], annots)


def test_an_unspliced_read_in_a_last_exon_is_a_3_fragment(tmp_path):
    exons = _fragment_exons(tmp_path)
    assert exons.fragment(520, 650) == ('last', 'A', '+')
    # its start up to 100 bases into the intron before, its end past the transcript's
    assert exons.fragment(420, 580) == ('last', 'A', '+')
    assert exons.fragment(380, 580) is None
    # on -, by its end
    assert exons.fragment(4900, 5050) == ('last', 'D', '-')
    assert exons.fragment(4900, 5250) is None


def test_an_unspliced_read_in_a_first_exon_is_a_5_fragment(tmp_path):
    exons = _fragment_exons(tmp_path)
    assert exons.fragment(50, 180) == ('first', 'A', '+')
    assert exons.fragment(5220, 5400) == ('first', 'D', '-')


def test_an_unspliced_read_elsewhere_is_not_a_fragment(tmp_path):
    exons = _fragment_exons(tmp_path)
    # in an internal exon
    assert exons.fragment(320, 380) is None
    # intronic
    assert exons.fragment(420, 480) is None
    # a single-exon gene's
    assert exons.fragment(8100, 8200) is None


def test_a_fragment_is_of_the_gene_sharing_the_most_bases_with_it(tmp_path):
    exons = _fragment_exons(tmp_path)
    # in B's last exon, 1400-1500, and the antisense gene's, 1450-1550: 60 bases of
    # B's and 30 of AS's
    assert exons.fragment(1420, 1480) == ('last', 'B', '+')
    # 50 of each: no one gene
    assert exons.fragment(1450, 1500) is None
    # 90 of AS's last exon, which on - ends at its splice site, 1550
    assert exons.fragment(1460, 1560) == ('last', 'AS', '-')


def test_a_fragments_trusted_strand_must_be_its_genes(tmp_path):
    exons = _fragment_exons(tmp_path)
    assert exons.fragment(520, 650, '+') == ('last', 'A', '+')
    assert exons.fragment(520, 650, '-') is None
    # only the genes of the trusted strand count, so a tie between strands isn't one
    assert exons.fragment(1450, 1500, '+') == ('last', 'B', '+')
    assert exons.fragment(1450, 1500, '-') == ('last', 'AS', '-')


class _Read:
    def __init__(self, start, end, strand='+'):
        self.start, self.end, self.strand = start, end, strand


def test_fragments_are_a_chain_per_gene_and_kind(tmp_path):
    reads = [_Read(520, 650), _Read(530, 640, '-'), _Read(50, 180), _Read(320, 380)]
    chains = fragment_chains(reads, _fragment_exons(tmp_path), trust_strand=False)
    assert chains == [ReadChain(('A',), (), '+', [ReadEnds(50, 180)], 'first'),
                      ReadChain(('A',), (), '+', [ReadEnds(520, 650), ReadEnds(530, 640)], 'last')]
    # with a trusted strand, the read on - isn't A's
    chains = fragment_chains(reads, _fragment_exons(tmp_path), trust_strand=True)
    assert [len(chain.ends) for chain in chains] == [1, 1]


def test_fragments_count_only_for_their_genes_starts_or_ends():
    a = ReadChain(('A',), A_JUNCS, '+', [ReadEnds(100, 600)])
    last = ReadChain(('A',), (), '+', [ReadEnds(520, 650)] * 2, 'last')
    first = ReadChain(('A',), (), '+', [ReadEnds(50, 180)], 'first')
    views, _, chain_genes, coverage = group_chains('chr1', ['s1'], [[a, last, first]], {})
    assert [list(chains) for chains in views['cluster'][0].values()] == [[A_JUNCS]]
    assert views['three'][0]['A'] == {A_JUNCS: a.ends, (): last.ends}
    assert views['five'][0]['A'] == {A_JUNCS: a.ends, (): first.ends}
    assert coverage['+'].covering((300, 400)) == {'s1': 1}
    alljuncs, ale, ends = extract_terminal_info(['s1'], views['three'], 'A', '+', 'ale')
    assert sorted(ends) == [(600, 's1'), (650, 's1'), (650, 's1')]
    assert list(ale) == [500] and alljuncs.keys() == set(A_JUNCS)
    _, _, starts = extract_terminal_info(['s1'], views['five'], 'A', '+', 'afe')
    assert sorted(starts) == [(50, 's1'), (100, 's1')]


def test_fragments_are_written_and_read_back_and_old_files_still_read(tmp_path):
    chains = [ReadChain(('A',), A_JUNCS, '+', [ReadEnds(100, 600)]),
              ReadChain(('A',), (), '+', [ReadEnds(520, 650)], 'last')]
    path = str(tmp_path / 's1.read_chains.tsv')
    write_read_chains(path, 'chr1', chains)
    assert read_read_chains(path) == chains
    with open(path, 'w') as fh:
        fh.write('A\tchr1\t200-300,400-500\t+\t100-600\n')
    assert read_read_chains(path) == chains[:1]
