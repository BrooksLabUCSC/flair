"""The formula matrix: the design handed to the R differential tests.

Columns, in order:

  sample_id    names one counts column of the matrix being tested
  condition    the condition that sample belongs to
  batch        the batch it belongs to, written only when batch is modelled

R reads this with row.names = 1, so sample_id is the row name and must be unique.
DESeq2 chooses its design by whether a batch column is present, condition + batch
when it is and condition alone when it is not, so the column is written only when
the batch is to be modelled.  DRIMSeq is told by an option instead and ignores the
difference.

This is the same shape as a sample info TSV, but it describes one comparison
rather than the samples of a counts matrix: the sample ids are the ones the
comparison uses, and batch may be absent.
"""
from contextlib import contextmanager
from collections import namedtuple
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('sample_id', 'condition')
COLUMNS_WITH_BATCH = ('sample_id', 'condition', 'batch')

class FormulaRow(namedtuple('FormulaRow', COLUMNS_WITH_BATCH)):
    "one sample of the comparison; batch is None when it is not modelled"
    __slots__ = ()

class FormulaMatrixReader(TsvReader):
    "batch is there only when it was modelled; the row then carries None for it"
    def __init__(self, formula_matrix_tsv):
        super().__init__(formula_matrix_tsv, defaultColType=str)

    @property
    def with_batch(self):
        return 'batch' in self.columns

    def __iter__(self):
        for row in super().__iter__():
            yield FormulaRow(row.sample_id, row.condition,
                             row.batch if self.with_batch else None)

class FormulaMatrixWriter(TsvWriter):
    def __init__(self, formula_matrix_tsv, *, with_batch):
        self.with_batch = with_batch
        super().__init__(formula_matrix_tsv,
                         columns=COLUMNS_WITH_BATCH if with_batch else COLUMNS,
                         defaultColType=str)

    def writeFormulaRow(self, row):
        if self.with_batch:
            self.writeRow((row.sample_id, row.condition, row.batch))
        else:
            self.writeRow((row.sample_id, row.condition))

@contextmanager
def formula_matrix_writer(formula_matrix_tsv, *, with_batch):
    "write a formula matrix atomically; the path appears only once complete"
    fileOps.ensureFileDir(formula_matrix_tsv)
    with fileOps.AtomicFileCreate(formula_matrix_tsv) as tmp_tsv:
        with FormulaMatrixWriter(tmp_tsv, with_batch=with_batch) as writer:
            yield writer

def write_formula_matrix(formula_matrix_tsv, rows, *, with_batch):
    with formula_matrix_writer(formula_matrix_tsv, with_batch=with_batch) as writer:
        for row in rows:
            writer.writeFormulaRow(row)
