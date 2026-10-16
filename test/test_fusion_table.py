"""
OUTPUT.fusions.tsv: per final fusion, its preliminary and final reads, isoforms, and best
supported breakpoint and isoform, most final reads first.
"""
from flair.fusion_table import breakpoint, isoform_breakpoints, fusion_rows, write_table, FUSION_TABLE_COLUMNS

# pieces are (position in fusion, chrom, start, end, strand) in 0-based half-open coordinates
ERGIC3 = (0, 'chr20', 35542078, 35548865, '+')
KRAS = (1, 'chr12', 25205245, 25227412, '-')


def test_breakpoint_is_the_last_and_first_bases_1_based():
    "ERGIC3 (+) ends at its piece's last base, KRAS (-) starts at its piece's highest base"
    assert breakpoint(ERGIC3, KRAS) == 'chr20:35548865(+)--chr12:25227412(-)'
    minus5 = (0, 'chr17', 64502925, 64503571, '-')
    plus3 = (1, 'chr1', 1000, 2000, '+')
    assert breakpoint(minus5, plus3) == 'chr17:64502926(-)--chr1:1001(+)'


def test_each_junction_in_fusion_order():
    ddx5 = (1, 'chr17', 64502925, 64503571, '-')
    kras = (2, 'chr12', 25205245, 25245395, '-')
    first = (0, 'chr20', 35542020, 35542941, '+')
    assert isoform_breakpoints([ERGIC3, KRAS]) == 'chr20:35548865(+)--chr12:25227412(-)'
    assert isoform_breakpoints([kras, first, ddx5]) == \
        'chr20:35542941(+)--chr17:64503571(-),chr17:64502926(-)--chr12:25245395(-)'


def test_rows_top_isoform_and_breakpoint_and_order():
    "two isoforms of one fusion share a breakpoint, and a smaller fusion sorts after it"
    other5, other3 = (0, 'chr1', 100, 500, '+'), (1, 'chr2', 800, 900, '+')
    other3_alt = (1, 'chr2', 850, 900, '+')
    iso_pieces = {'iso1': [ERGIC3, KRAS], 'iso2': [ERGIC3, KRAS],
                  'iso3': [other5, other3], 'iso4': [other5, other3_alt], 'iso5': [other5, other3_alt]}
    iso_contig = {'iso1': 'G1--G2', 'iso2': 'G1--G2', 'iso3': 'G3--chr2:800-900',
                  'iso4': 'G3--chr2:800-900', 'iso5': 'G3--chr2:800-900'}
    iso_reads = {'iso1': 17, 'iso2': 11, 'iso3': 5, 'iso4': 3, 'iso5': 3}
    prelim = {'G1--G2': 32, 'G3--chr2:800-900': 12}
    names = {'G1': 'ERGIC3', 'G2': 'KRAS', 'G3': 'GENE3'}
    iso_genes = {iso: tuple(contig.split('--')) for iso, contig in iso_contig.items()}
    rows = fusion_rows(prelim, iso_contig, iso_genes, iso_reads, iso_pieces, names)
    assert rows == [
        ['ERGIC3--KRAS', 'G1--G2', 32, 28, 2, 'chr20:35548865(+)--chr12:25227412(-)', 28, 'iso1', 17],
        # iso4 and iso5 share a breakpoint with 6 reads, more than iso3's 5
        ['GENE3--chr2:800-900', 'G3--chr2:800-900', 12, 11, 3, 'chr1:500(+)--chr2:851(+)', 6, 'iso3', 5]]


def test_empty_table_has_the_header(tmp_path):
    tsv = tmp_path / 'f.fusions.tsv'
    write_table([], str(tsv))
    assert tsv.read_text() == '\t'.join(FUSION_TABLE_COLUMNS) + '\n'


def test_rows_are_by_final_genes_with_their_contigs_prelim_reads():
    "isoforms of two contigs given the same genes are one fusion, with both contigs' reads"
    iso_pieces = {'iso1': [ERGIC3, KRAS], 'iso2': [ERGIC3, KRAS]}
    iso_contig = {'iso1': 'G1--G2', 'iso2': 'G1b--G2'}
    iso_genes = {'iso1': ('G1', 'G2'), 'iso2': ('G1', 'G2')}
    rows = fusion_rows({'G1--G2': 20, 'G1b--G2': 10}, iso_contig, iso_genes, {'iso1': 17, 'iso2': 11}, iso_pieces,
                       {'G1': 'ERGIC3', 'G2': 'KRAS'})
    assert [row[:5] for row in rows] == [['ERGIC3--KRAS', 'G1--G2', 30, 28, 2]]
