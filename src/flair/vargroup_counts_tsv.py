"""Reads grouped by the set of variants they carry, from flair variantquant.

Columns, in order:

  chrom            the reference the reads aligned to
  gene             the gene the group was assigned to, the commonest among the
                   reads' variants
  num_reads        how many reads carry this set of positions
  num_positions    how many variant positions the set has
  mods_per_read    how many of the positions each read carries, one per read;
                   comma separated in the file and a list on the row
  varpos           the variant positions, likewise

Positions are one based, as the VCF they came from writes them.


"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaIntListType, commaListType

COLUMNS = ('chrom', 'gene', 'num_reads', 'num_positions', 'mods_per_read', 'varpos')
TYPE_MAP = {'chrom': str, 'gene': str, 'num_reads': int, 'num_positions': int,
            'mods_per_read': commaIntListType, 'varpos': commaListType}

class VarGroupCountsReader(TsvReader):
    def __init__(self, vargroup_counts_tsv):
        super().__init__(vargroup_counts_tsv, typeMap=TYPE_MAP)

class VarGroupCountsWriter(TsvWriter):
    def __init__(self, vargroup_counts_tsv, *, outFh=None):
        super().__init__(vargroup_counts_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)

@contextmanager
def vargroup_counts_writer(vargroup_counts_tsv):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(vargroup_counts_tsv)
    with fileOps.AtomicFileCreate(vargroup_counts_tsv) as tmp_tsv:
        with VarGroupCountsWriter(tmp_tsv) as writer:
            yield writer
