"""VCF variants grouped by the region they fall in, written by
annotate_group_vcf_vars.

Columns, in order:

  chrom        the chromosome
  region       the region the variant falls in
  pos          the position, one based as the VCF writes it
  ref          the reference allele
  alts         the alternate alleles, comma separated in the file and a list
               on the row
  name         the variant name from the VCF
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType

COLUMNS = ('chrom', 'region', 'pos', 'ref', 'alts', 'name')

class RegionVarsReader(TsvReader):
    def __init__(self, region_vars_tsv):
        super().__init__(region_vars_tsv, defaultColType=str,
                         typeMap={'alts': commaListType})

class RegionVarsWriter(TsvWriter):
    def __init__(self, region_vars_tsv, *, outFh=None):
        super().__init__(region_vars_tsv, columns=COLUMNS, defaultColType=str,
                         typeMap={'alts': commaListType},
                         outFh=outFh)
