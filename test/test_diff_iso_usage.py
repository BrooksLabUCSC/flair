"""
diff_iso_usage of the counts of flair spliceevents and variantquant, besides a
quantify counts matrix (test-diff-iso-usage)
"""
from flair import splicing_events_tsv, var_counts_tsv
from flair.counts_matrix_tsv import CountsRow, write_counts_matrix
from flair.diff_iso_usage import add_adjusted_pvalues, diff_iso_usage
from flair.iso_usage_tsv import IsoUsageReader

EVENT = 'es|chr1:101-200(+)|G1'


def _usage(path):
    with IsoUsageReader(path) as reader:
        return [(r.geneID, r.isoID, r['this_iso_s1_count'], r['other_isos_s1_count'], r['s1_PSI'])
                for r in reader]


def test_an_events_other_row_counts_toward_its_alternatives_but_is_not_tested(tmp_path):
    counts = str(tmp_path / 'x.spliceevents.counts.tsv')
    with open(counts, 'w') as fh:
        fh.write('\t'.join(list(splicing_events_tsv.EVENT_COLUMNS) + ['s1', 's2']) + '\n')
        for alternative, s1, s2 in (('inc', 6, 1), ('exc', 2, 7), ('other', 2, 0)):
            fh.write('\t'.join([EVENT, alternative, 'es', 'G1'] + [''] * 8 + [str(s1), str(s2)]) + '\n')
    out = str(tmp_path / 'usage.tsv')
    diff_iso_usage(counts, 's1', 's2', out)
    assert _usage(out) == [(EVENT, EVENT + '|inc', '6', '4', '0.6'),
                           (EVENT, EVENT + '|exc', '2', '8', '0.2')]


def test_an_alleles_share_of_its_sites_reads(tmp_path):
    counts = str(tmp_path / 'x.var.counts.tsv')
    with var_counts_tsv.var_counts_writer(counts, ['s1', 's2']) as writer:
        writer.writeAllele('chr1:5', 'A', 'G', ('.', '.', '.'), 'G1', [(9, 1), (2, 8)])
    out = str(tmp_path / 'usage.tsv')
    diff_iso_usage(counts, 's1', 's2', out)
    assert _usage(out) == [('chr1:5:A>G', 'chr1:5:A>G|yes', '1', '9', '0.1'),
                           ('chr1:5:A>G', 'chr1:5:A>G|no', '9', '1', '0.9')]


def _gene_matrix(tmp_path):
    "one gene of samples s1, s2: i1 and i2 switch, i3 of 1 read, i4 steady"
    counts = str(tmp_path / 'x.counts.tsv')
    write_counts_matrix(counts, ['s1', 's2'], [CountsRow('G1', 'i1', [20, 2]), CountsRow('G1', 'i2', [2, 20]),
                                               CountsRow('G1', 'i3', [1, 0]), CountsRow('G1', 'i4', [5, 5])])
    return counts


def _rows(path):
    with IsoUsageReader(path) as reader:
        return {r.isoID: r for r in reader}


def test_an_isoform_of_fewer_reads_than_min_expression_is_not_tested(tmp_path):
    out = str(tmp_path / 'usage.tsv')
    diff_iso_usage(_gene_matrix(tmp_path), 's1', 's2', out, min_expression=5)
    rows = _rows(out)
    # i3 is written, with its PSIs, its p-values NA; its read still counts for the others
    assert (rows['i3'].fisher_pval, rows['i3'].fisher_adj_pval, rows['i3'].s1_PSI) == ('NA', 'NA', '0.036')
    assert rows['i1']['other_isos_s1_count'] == '8'
    assert all(rows[iso].fisher_adj_pval != 'NA' for iso in ('i1', 'i2', 'i4'))


def test_each_isoforms_other_counts_are_its_genes_total_less_its_own(tmp_path):
    out = str(tmp_path / 'usage.tsv')
    diff_iso_usage(_gene_matrix(tmp_path), 's1', 's2', out)
    rows = _rows(out)
    # the gene's 28 and 27 reads in s1 and s2
    assert {iso: (r['other_isos_s1_count'], r['other_isos_s2_count']) for iso, r in rows.items()} == {
        'i1': ('8', '25'), 'i2': ('26', '7'), 'i3': ('27', '27'), 'i4': ('23', '22')}


def test_pvalues_are_adjusted_over_the_rows_tested():
    # Benjamini-Hochberg of 0.01, 0.04, 0.03 over three tests: p * 3 / rank, then
    # the minimum of those at or above it; the untested row stays NA
    rows = add_adjusted_pvalues([['G1', 'i1', 0.01], ['G1', 'i2', 'NA'], ['G2', 'i3', 0.04], ['G2', 'i4', 0.03]])
    assert [row[3] if row[3] == 'NA' else round(row[3], 6) for row in rows] == [0.03, 'NA', 0.04, 0.04]
    assert add_adjusted_pvalues([['G1', 'i1', 'NA']]) == [['G1', 'i1', 'NA', 'NA']]
