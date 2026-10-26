"""Isoforms whose usage within their gene changed, from the diffexp t-test on
counts normalized by gene.

Columns, in order:

  gene_id        the gene
  isoform_id     the isoform
  delta_usage    change in the isoform's share of its gene, condition_b minus
                 condition_a, between -1 and 1
  adj_pval       Holm-Sidak adjusted p-value; only rows below 0.05 are written

Written by diffexp as isoforms_sig_exp_change_norm_by_gene.tsv.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('gene_id', 'isoform_id', 'delta_usage', 'adj_pval')
TYPE_MAP = {'gene_id': str, 'isoform_id': str,
            'delta_usage': float, 'adj_pval': float}

class IsoUsageChangeReader(TsvReader):
    def __init__(self, iso_usage_change_tsv):
        super().__init__(iso_usage_change_tsv, typeMap=TYPE_MAP)

class IsoUsageChangeWriter(TsvWriter):
    def __init__(self, iso_usage_change_tsv):
        super().__init__(iso_usage_change_tsv, columns=COLUMNS, typeMap=TYPE_MAP)

@contextmanager
def iso_usage_change_writer(iso_usage_change_tsv):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(iso_usage_change_tsv)
    with fileOps.AtomicFileCreate(iso_usage_change_tsv) as tmp_tsv:
        with IsoUsageChangeWriter(tmp_tsv) as writer:
            yield writer
