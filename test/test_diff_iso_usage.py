"""
diff_iso_usage of the counts of flair spliceevents and variantquant, besides a
quantify counts matrix (test-diff-iso-usage)
"""
from flair import splicing_events_tsv, var_counts_tsv
from flair.diff_iso_usage import diff_iso_usage
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
    assert _usage(out) == [(EVENT, EVENT + '|inc', '6.0', '4.0', '0.6'),
                           (EVENT, EVENT + '|exc', '2.0', '8.0', '0.2')]


def test_an_alleles_share_of_its_sites_reads(tmp_path):
    counts = str(tmp_path / 'x.var.counts.tsv')
    with var_counts_tsv.var_counts_writer(counts, ['s1', 's2']) as writer:
        writer.writeAllele('chr1:5', 'A', 'G', ('.', '.', '.'), 'G1', [(9, 1), (2, 8)])
    out = str(tmp_path / 'usage.tsv')
    diff_iso_usage(counts, 's1', 's2', out)
    assert _usage(out) == [('chr1:5:A>G', 'chr1:5:A>G|yes', '1.0', '9.0', '0.1'),
                           ('chr1:5:A>G', 'chr1:5:A>G|no', '9.0', '1.0', '0.9')]
