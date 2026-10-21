"""Read counts of each allele group, from flair alleles.

Columns, in order:

  phase_set        the phase set the group belongs to
  allele_group     the group, named by the variants it carries
  tumor_counts     reads supporting it in the tumor BAM
  normal_counts    reads supporting it in the normal BAM, with --norm_bam
  somatic          yes when the normal support is below --read_support, so the
                   group looks somatic; written only with --norm_bam

The last two columns exist only when a normal BAM was given, since without one
there is nothing to call somatic against.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('phase_set', 'allele_group', 'tumor_counts')
NORMAL_COLUMNS = ('normal_counts', 'somatic')

def columns(with_normal):
    return list(COLUMNS) + (list(NORMAL_COLUMNS) if with_normal else [])

class AlleleGroupCountsReader(TsvReader):
    "the normal columns are there only when the run had a normal BAM"
    def __init__(self, allele_group_counts_tsv):
        super().__init__(allele_group_counts_tsv, defaultColType=str)

    @property
    def with_normal(self):
        return NORMAL_COLUMNS[0] in self.columns

class AlleleGroupCountsWriter(TsvWriter):
    def __init__(self, allele_group_counts_tsv, *, with_normal):
        super().__init__(allele_group_counts_tsv, columns=columns(with_normal),
                         defaultColType=str)

@contextmanager
def allele_group_counts_writer(allele_group_counts_tsv, *, with_normal):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(allele_group_counts_tsv)
    with fileOps.AtomicFileCreate(allele_group_counts_tsv) as tmp_tsv:
        with AlleleGroupCountsWriter(tmp_tsv, with_normal=with_normal) as writer:
            yield writer
