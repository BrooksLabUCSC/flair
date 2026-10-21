"""
Round trip tests for the TSV format modules: what a writer writes, its reader
reads back as the same values.  Most of these readers have no caller yet, so
this is the only thing that runs them.
"""
import pytest
from flair.aaseq_pred_tsv import AaSeqPredReader, AaSeqPredWriter
from flair.aaseq_tsv import AaSeqReader, write_aaseqs
from flair.allele_group_counts_tsv import AlleleGroupCountsReader, allele_group_counts_writer
from flair.combined_map_tsv import CombinedMapReader, CombinedMapWriter
from flair.deseq2_counts_tsv import Deseq2CountsReader, write_deseq2_counts
from flair.diffsplice_fishers_tsv import DiffspliceFishersReader, diffsplice_fishers_writer
from flair.drimseq_counts_tsv import DrimSeqCountsReader, drimseq_counts_writer
from flair.es_events_tsv import EsEventsReader, EsEventsWriter
from flair.event_quant_tsv import EventQuantReader, columns as event_quant_columns, event_quant_writer
from flair.formula_matrix_tsv import FormulaMatrixReader, FormulaRow, write_formula_matrix
from flair.indel_vars_tsv import IndelVarsReader, IndelVarsWriter
from flair.iso_allele_counts_tsv import IsoAlleleCountsReader, iso_allele_counts_writer
from flair.iso_usage_change_tsv import IsoUsageChangeReader, iso_usage_change_writer
from flair.iso_usage_tsv import IsoUsageReader, iso_usage_writer
from flair.marked_isoforms_bed import MarkedIsoformsReader, MarkedIsoformsWriter
from flair.read_ends_tsv import ReadEndsReader, ReadEndsWriter
from flair.read_map_tsv import ReadMapReader, ReadMapWriter
from flair.retained_introns_tsv import RetainedIntronsReader, RetainedIntronsWriter
from flair.seq_lengths_tsv import (SeqLengthHistogramReader, SeqLengthHistogramWriter,
                                   SeqLengthsReader, SeqLengthsWriter)
from flair.splicing_events_tsv import SplicingEventsReader, SplicingEventsWriter
from flair.splicing_totals_tsv import SplicingTotalsReader, SplicingTotalsWriter
from flair.read_chains_tsv import ReadChainsReader, read_chains_writer
from flair.splicing_outliers_tsv import SplicingOutliersReader, SplicingOutliersWriter
from flair.tpm_tsv import CountsRow, TpmReader, tpm_writer
from flair.transcript_counts_tsv import TranscriptCountsReader, TranscriptCountsWriter
from flair.var_counts_tsv import VarCountsReader, var_counts_writer

def chained_messages(ex):
    "every message on an exception chain, for checking why a row was refused"
    messages = []
    while ex is not None:
        messages.append(str(ex))
        ex = ex.__cause__
    return '\n'.join(messages)

def tsv_path(tmp_path, name):
    return str(tmp_path / name)

def written(writer_class, path, rows, **kwargs):
    with writer_class(path, **kwargs) as writer:
        for row in rows:
            writer.writeRow(row)
    return path

def test_read_map_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'reads.tsv')
    with ReadMapWriter(path) as writer:
        writer.writeReads('iso1', ['readA', 'readB'])
        writer.writeReads('iso2', [])
    rows = list(ReadMapReader(path))
    assert [(r.name, r.reads) for r in rows] == [('iso1', ['readA', 'readB']),
                                                 ('iso2', [])]

def test_aaseq_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'x.aaseq.tsv')
    write_aaseqs(path, {'MKV': 'FLP00000001', 'MKW': 'FLP00000002'})
    rows = list(AaSeqReader(path))
    assert [(r.aaseq_id, r.aaseq) for r in rows] == [('FLP00000001', 'MKV'),
                                                     ('FLP00000002', 'MKW')]

def test_aaseq_pred_round_trip(tmp_path):
    path = written(AaSeqPredWriter, tsv_path(tmp_path, 'pred.tsv'),
                   [('iso1', 'PRO', ['5utr', '3utr'], 'MKV'), ('iso2', 'PTC', [], 'M')])
    rows = list(AaSeqPredReader(path))
    assert [r.utr_vars for r in rows] == [['5utr', '3utr'], []]

def test_allele_group_counts_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'ag.counts.tsv')
    with allele_group_counts_writer(path, with_normal=True) as writer:
        writer.writeRow(('ps1', 'ag1', '7', '2', 'yes'))
    with AlleleGroupCountsReader(path) as reader:
        assert reader.with_normal
        assert [r.somatic for r in reader] == ['yes']

def test_allele_group_counts_without_normal(tmp_path):
    path = tsv_path(tmp_path, 'ag.counts.tsv')
    with allele_group_counts_writer(path, with_normal=False) as writer:
        writer.writeRow(('ps1', 'ag1', '7'))
    with AlleleGroupCountsReader(path) as reader:
        assert not reader.with_normal

def test_iso_allele_counts_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'ia.counts.tsv')
    with iso_allele_counts_writer(path, with_normal=False) as writer:
        writer.writeRow(('gene1', 'iso1', 'ps1', 'ag1', 'iso1|ag1', 'FLP00000001', '5'))
    with IsoAlleleCountsReader(path) as reader:
        assert not reader.with_normal
        assert [r.aaseq_id for r in reader] == ['FLP00000001']

def test_combined_map_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'map.tsv')
    with CombinedMapWriter(path) as writer:
        writer.writeSources('FLT00000001', ['s1,s2:FLT00000003', 's3:FLT00000009'])
    rows = list(CombinedMapReader(path))
    assert rows[0].source_ids == ['s1,s2:FLT00000003', 's3:FLT00000009']

def test_deseq2_counts_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'ds2.tsv')
    write_deseq2_counts(path, ['s1', 's2'], ['gene1', 'gene2'], [[1, 2], [3, 4]])
    with Deseq2CountsReader(path) as reader:
        assert reader.sample_columns == ['s1', 's2']
        assert [(r.feature_id, r.s1, r.s2) for r in reader] == [('gene1', 1, 2),
                                                                ('gene2', 3, 4)]

def test_drimseq_counts_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'drim.tsv')
    with drimseq_counts_writer(path, ['s1']) as writer:
        writer.writeIsoform('gene1', 'iso1', [5])
        writer.writeIsoform('gene1', 'iso2', [6])
    with DrimSeqCountsReader(path) as reader:
        assert reader.sample_columns == ['s1']
        assert [(r.irow, r.feature_id, r.s1) for r in reader] == [(0, 'iso1', 5),
                                                                  (1, 'iso2', 6)]

def test_formula_matrix_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'formula.tsv')
    rows = [FormulaRow('s1', 'ctl', 'b1'), FormulaRow('s2', 'test', 'b2')]
    write_formula_matrix(path, rows, with_batch=True)
    with FormulaMatrixReader(path) as reader:
        assert reader.with_batch
        assert list(reader) == rows

def test_formula_matrix_without_batch(tmp_path):
    path = tsv_path(tmp_path, 'formula.tsv')
    write_formula_matrix(path, [FormulaRow('s1', 'ctl', 'b1')], with_batch=False)
    with FormulaMatrixReader(path) as reader:
        assert not reader.with_batch
        assert list(reader) == [FormulaRow('s1', 'ctl', None)]

def test_es_events_round_trip(tmp_path):
    path = written(EsEventsWriter, tsv_path(tmp_path, 'es.tsv'),
                   [('chr1:100-200', '+', 2, 1, ['isoA', 'isoB'], ['isoC'])])
    rows = list(EsEventsReader(path))
    assert (rows[0].inclusion_isos, rows[0].exclusion_isos) == (['isoA', 'isoB'], ['isoC'])

def test_event_quant_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'events.quant.tsv')
    with event_quant_writer(path, ['s1', 's2']) as writer:
        writer.writeSide('inclusion', 'chr1:100', 'chr1:100-200', [3, 4], ['isoA'])
    with EventQuantReader(path) as reader:
        assert reader.sample_columns == ['s1', 's2']
        rows = list(reader)
    assert rows[0].feature_id == 'inclusion_chr1:100'
    assert rows[0].isoform_ids == ['isoA']

def test_event_quant_rejects_other_columns(tmp_path):
    path = tsv_path(tmp_path, 'bad.tsv')
    with open(path, 'w') as fh:
        fh.write('a\tb\tc\n')
    with pytest.raises(Exception, match="event quant TSV must have the columns"):
        EventQuantReader(path)

def test_diffsplice_fishers_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'fishers.tsv')
    columns = event_quant_columns(['s1', 's2'])
    with diffsplice_fishers_writer(path, columns, 's1', 's2') as writer:
        writer.writeRow(['inclusion_e', 'e', '1', '2', ['isoA', 'isoB'], 0.5])
    with DiffspliceFishersReader(path) as reader:
        assert reader.pval_column == 's1-s2_pval'
        assert [r.isoform_ids for r in reader] == [['isoA', 'isoB']]

def test_indel_vars_round_trip(tmp_path):
    path = written(IndelVarsWriter, tsv_path(tmp_path, 'indels.tsv'),
                   [('chr1', 100, 'insertion', 'A', 'AT', 3, 10, ['dn'])])
    assert [r.filters for r in IndelVarsReader(path)] == [['dn']]

def test_iso_usage_change_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'usage.tsv')
    with iso_usage_change_writer(path) as writer:
        writer.writeColumns(gene_id='gene1', isoform_id='iso1',
                            delta_usage=0.25, adj_pval=0.01)
    rows = list(IsoUsageChangeReader(path))
    assert (rows[0].delta_usage, rows[0].adj_pval) == (0.25, 0.01)

def test_iso_usage_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'iso_usage.tsv')
    with iso_usage_writer(path, 'sA', 'sB') as writer:
        writer.writeRow(['gene1', 'iso1', 'NA', 1, 2, 3, 4, 'NA', 'NA', 'NA'])
    rows = list(IsoUsageReader(path))
    assert rows[0].geneID == 'gene1'
    assert rows[0].sA_PSI == 'NA'

def test_marked_isoforms_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'marked.bed')
    bed_row = ['chr1', '100', '200', 'iso1', '0', '+']
    with MarkedIsoformsWriter(path, len(bed_row)) as writer:
        writer.writeIsoform(bed_row, True)
        writer.writeIsoform(bed_row, False)
    rows = list(MarkedIsoformsReader(path))
    assert [r.retains_intron for r in rows] == [1, 0]

def test_read_ends_round_trip(tmp_path):
    path = written(ReadEndsWriter, tsv_path(tmp_path, 'ends.tsv'),
                   [('readA', 'iso1', 0, 5, 10, 2, 7, 12),
                    ('readB', 'iso2', None, None, 3, None, None, 4)])
    rows = list(ReadEndsReader(path))
    assert rows[0].start_sj_index == 0
    assert rows[1].start_sj_index is None
    assert rows[1].end_tend_dist == 4

def test_retained_introns_round_trip(tmp_path):
    path = written(RetainedIntronsWriter, tsv_path(tmp_path, 'introns.tsv'),
                   [('chr1', 100, 200, '+')])
    assert [(r.start, r.end) for r in RetainedIntronsReader(path)] == [(100, 200)]

def test_seq_lengths_round_trip(tmp_path):
    path = written(SeqLengthsWriter, tsv_path(tmp_path, 'len.tsv'),
                   [('seqA', 100), ('seqB', 200)])
    assert [(r.name, r.length) for r in SeqLengthsReader(path)] == [('seqA', 100),
                                                                    ('seqB', 200)]

def test_seq_length_histogram_round_trip(tmp_path):
    path = written(SeqLengthHistogramWriter, tsv_path(tmp_path, 'hist.tsv'),
                   [(100, 2), (200, 1)])
    assert [(r.length, r.count) for r in SeqLengthHistogramReader(path)] == [(100, 2),
                                                                             (200, 1)]

def test_splicing_events_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'events.tsv')
    with SplicingEventsWriter(path, ['s1', 's2']) as writer:
        writer.writeRow(['es|chr1:201-300(+)|G1', 'inclusion', 'es', ['G1'], ['NAME1'], ['G2', 'G3'], [],
                         ['chr1:101-200', 'chr1:301-400'], [], ['chr1:101-400'], ['chr1:201-300'], [],
                         '3', ''])
    with SplicingEventsReader(path) as reader:
        assert reader.sample_columns == ['s1', 's2']
        rows = list(reader)
    assert rows[0].junctions_included == ['chr1:101-200', 'chr1:301-400']
    assert rows[0].junctions_excluded == []
    assert (rows[0].additional_genes, rows[0].s1, rows[0].s2) == (['G2', 'G3'], '3', '')

def test_splicing_totals_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'totals.tsv')
    with SplicingTotalsWriter(path, ['s1', 's2']) as writer:
        writer.writeRow(['es|chr1:201-300(+)|G1', 'es', ['G1'], ['NAME1'], [], [], 'event_reads', 40, 0])
    with SplicingTotalsReader(path) as reader:
        assert reader.sample_columns == ['s1', 's2']
        rows = list(reader)
    assert (rows[0].total, rows[0].gene, rows[0].s1, rows[0].s2) == ('event_reads', ['G1'], 40, 0)

def test_splicing_outliers_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'outliers.tsv')
    with SplicingOutliersWriter(path) as writer:
        writer.writeRow(['es|chr1:201-300(+)|G1', 'inclusion', 'es', ['G1'], ['NAME1'], [], [], 's1',
                         0.5, 0.1, 0.9, 3, [7, 20], 0.4, 4.0])
    rows = list(SplicingOutliersReader(path))
    assert rows[0]['event_reads;total_locus_reads'] == [7, 20]
    assert (rows[0].gene, rows[0].sample) == (['G1'], 's1')

def test_read_chains_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'read_chains.tsv')
    with read_chains_writer(path) as writer:
        writer.writeRow([['A', 'B'], 'chr1', ((200, 300), (600, 1000)), '+', ((150, 1300),), ''])
        writer.writeRow([[], 'chr1', (), '-', ((10, 50),), 'last'])
    rows = list(ReadChainsReader(path))
    assert (rows[0].genes, rows[0].juncs, rows[0].ends) == (['A', 'B'], ((200, 300), (600, 1000)), ((150, 1300),))
    assert (rows[1].genes, rows[1].juncs, rows[1].fragment) == ([], (), 'last')

def test_transcript_counts_round_trip(tmp_path):
    path = written(TranscriptCountsWriter, tsv_path(tmp_path, 'counts.tsv'),
                   [('iso1', 5), ('iso2', 0)])
    assert [(r.transcript, r.count) for r in TranscriptCountsReader(path)] == [('iso1', 5),
                                                                               ('iso2', 0)]

def test_tpm_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'x.tpm.tsv')
    with tpm_writer(path, ['s1', 's2']) as writer:
        writer.writeCountsRow(CountsRow('gene1', 'iso1', [1.5, 2.5]))
    with TpmReader(path) as reader:
        assert [r.counts for r in reader] == [[1.5, 2.5]]

def test_var_counts_round_trip(tmp_path):
    path = tsv_path(tmp_path, 'var.counts.tsv')
    with var_counts_writer(path, ['s1', 's2']) as writer:
        writer.writeRow(['chr1:100', 'A', 'G', '+', '.', '.', 'gene1', (7, 3), (1, 0)])
    with VarCountsReader(path) as reader:
        assert reader.sample_columns == ['s1', 's2']
        rows = list(reader)
    assert (rows[0].alt, rows[0].strand, rows[0].s1, rows[0].s2) == ('G', '+', [7, 3], [1, 0])
