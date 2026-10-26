"""Transcripts per million, one row per isoform.

Columns, in order:

  gene_id      the gene the isoform belongs to
  isoform_id   the isoform
  <sample>     one column per sample, named by its sample id, holding TPM

The same shape as a counts matrix, with the counts replaced by TPM, so it is read
and written through counts_matrix_tsv with float values.  flair quantify --tpm
writes it beside the counts matrix.
"""
from flair.counts_matrix_tsv import (CountsMatrixReader, CountsRow,  # noqa: F401
                                     counts_matrix_writer, read_counts_rows,
                                     read_sample_columns)

class TpmReader(CountsMatrixReader):
    "the counts matrix reader with float values"
    def __init__(self, tpm_tsv):
        super().__init__(tpm_tsv, countType=float)

def tpm_writer(tpm_tsv, sample_columns):
    "write a TPM matrix atomically; values are floats rather than counts"
    return counts_matrix_writer(tpm_tsv, sample_columns, countType=float)

def read_tpm_rows(tpm_tsv):
    "every row of the matrix, in file order"
    with TpmReader(tpm_tsv) as reader:
        return list(reader)
