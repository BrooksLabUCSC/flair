"""The counts matrix handed to DRIMSeq, one row per isoform.

Columns, in order:

  irow         a serial number over the rows of this file.  DRIMSeq indexes its
               data frame by the first column and refuses duplicate values, which
               the isoform ids can carry when one isoform appears twice
  gene_id      the gene whose isoforms DRIMSeq compares against each other
  feature_id   the isoform
  <sample>     one integer column per sample, named as in the formula matrix

Written by diffexp into its working directory; not an output anyone reads
directly.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter

ID_COLUMNS = ('irow', 'gene_id', 'feature_id')

class DrimSeqCountsReader(TsvReader):
    def __init__(self, drimseq_counts_tsv):
        super().__init__(drimseq_counts_tsv, defaultColType=int,
                         typeMap={'gene_id': str, 'feature_id': str})

    @property
    def sample_columns(self):
        return self.columns[len(ID_COLUMNS):]

class DrimSeqCountsWriter(TsvWriter):
    def __init__(self, drimseq_counts_tsv, sample_columns):
        super().__init__(drimseq_counts_tsv,
                         columns=list(ID_COLUMNS) + list(sample_columns),
                         defaultColType=str)
        self.irow = 0

    def writeIsoform(self, gene_id, feature_id, counts):
        "rows are numbered in the order they are written"
        self.writeRow([self.irow, gene_id, feature_id] + list(counts))
        self.irow += 1

@contextmanager
def drimseq_counts_writer(drimseq_counts_tsv, sample_columns):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(drimseq_counts_tsv)
    with fileOps.AtomicFileCreate(drimseq_counts_tsv) as tmp_tsv:
        with DrimSeqCountsWriter(tmp_tsv, sample_columns) as writer:
            yield writer
