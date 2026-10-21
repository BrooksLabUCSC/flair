"""The length of every sequence in a FASTA, and the histogram of those lengths.

Sequence lengths:

  name         the FASTA record name, the whole line after '>'
  length       the number of residues

Histogram:

  length       a length that at least one sequence has
  count        how many sequences have it, in increasing length order

Both written by fasta_seq_lengths.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

LENGTH_COLUMNS = ('name', 'length')
HISTOGRAM_COLUMNS = ('length', 'count')

class SeqLengthsReader(TsvReader):
    def __init__(self, seq_lengths_tsv):
        super().__init__(seq_lengths_tsv,
                         typeMap={'name': str, 'length': int})

class SeqLengthHistogramReader(TsvReader):
    def __init__(self, histogram_tsv):
        super().__init__(histogram_tsv, defaultColType=int)

class SeqLengthsWriter(TsvWriter):
    def __init__(self, seq_lengths_tsv):
        super().__init__(seq_lengths_tsv, columns=LENGTH_COLUMNS,
                         typeMap={'name': str, 'length': int})

class SeqLengthHistogramWriter(TsvWriter):
    def __init__(self, histogram_tsv):
        super().__init__(histogram_tsv, columns=HISTOGRAM_COLUMNS,
                         defaultColType=int)
