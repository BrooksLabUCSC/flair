"""Introns that an isoform retains, one row per distinct intron.

Columns, in order:

  chrom        the chromosome
  start        the intron start, zero based
  end          the intron end
  strand       the strand of the isoform retaining it

Written by mark_intron_retention beside its marked BED.  An intron appears once
however many isoforms retain it.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('chrom', 'start', 'end', 'strand')
TYPE_MAP = {'chrom': str, 'start': int, 'end': int, 'strand': str}

class RetainedIntronsReader(TsvReader):
    def __init__(self, retained_introns_tsv):
        super().__init__(retained_introns_tsv, typeMap=TYPE_MAP)

class RetainedIntronsWriter(TsvWriter):
    def __init__(self, retained_introns_tsv):
        super().__init__(retained_introns_tsv, columns=COLUMNS, typeMap=TYPE_MAP)
