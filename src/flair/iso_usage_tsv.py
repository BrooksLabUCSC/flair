"""Per-isoform usage compared between two samples, from diff_iso_usage.

Columns, in order:

  geneID                       the gene
  isoID                        the isoform
  fisher_pval                  Fisher's exact p-value for this isoform against
                               the rest of its gene, NA when either sample has
                               no reads to test
  this_iso_<sample1>_count     reads on this isoform in each sample
  this_iso_<sample2>_count
  other_isos_<sample1>_count   reads on the gene's other isoforms
  other_isos_<sample2>_count
  <sample1>_PSI                this isoform's share of its gene, each sample
  <sample2>_PSI
  delta_PSI                    sample2 PSI minus sample1 PSI

The sample names are in the column names, so the sign of delta_PSI can be read
from the file without knowing the order the samples were given in.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter

def columns(colname1, colname2):
    return ('geneID', 'isoID', 'fisher_pval',
            f'this_iso_{colname1}_count', f'this_iso_{colname2}_count',
            f'other_isos_{colname1}_count', f'other_isos_{colname2}_count',
            f'{colname1}_PSI', f'{colname2}_PSI', 'delta_PSI')

class IsoUsageReader(TsvReader):
    "values stay strings: a p-value or a PSI can be the string NA"
    def __init__(self, iso_usage_tsv):
        super().__init__(iso_usage_tsv, defaultColType=str)

class IsoUsageWriter(TsvWriter):
    def __init__(self, iso_usage_tsv, colname1, colname2):
        # every value is written as given: counts, p-values and the NA that marks
        # an isoform with nothing to test
        super().__init__(iso_usage_tsv, columns=columns(colname1, colname2),
                         defaultColType=str)

@contextmanager
def iso_usage_writer(iso_usage_tsv, colname1, colname2):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(iso_usage_tsv)
    with fileOps.AtomicFileCreate(iso_usage_tsv) as tmp_tsv:
        with IsoUsageWriter(tmp_tsv, colname1, colname2) as writer:
            yield writer
