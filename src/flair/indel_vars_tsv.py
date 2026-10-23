"""Insertions and deletions called by identify_vars, one row per variant.

Columns, in order:

  chrom        the chromosome
  pos          the position, one based as a VCF writes it
  indel_type   insertion or deletion
  ref_seq      the reference sequence
  var_seq      the variant sequence
  indel_reads  reads carrying the variant
  tot_cov      reads covering the position
  filters      filter tags, dn marking a dense cluster of variants around this
               one; comma separated in the file and a list on the row
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType

COLUMNS = ('chrom', 'pos', 'indel_type', 'ref_seq', 'var_seq', 'indel_reads',
           'tot_cov', 'filters')
TYPE_MAP = {'chrom': str, 'pos': int, 'indel_type': str, 'ref_seq': str,
            'var_seq': str, 'indel_reads': int, 'tot_cov': int,
            'filters': commaListType}

class IndelVarsReader(TsvReader):
    def __init__(self, indel_vars_tsv):
        super().__init__(indel_vars_tsv, typeMap=TYPE_MAP)

class IndelVarsWriter(TsvWriter):
    def __init__(self, indel_vars_tsv, *, outFh=None):
        super().__init__(indel_vars_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)
