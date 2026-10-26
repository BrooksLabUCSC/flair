"""Read counts of each isoform and allele group combination, from flair
isoalleles.

Columns, in order:

  gene                      the gene, empty on the rows that count an allele
                            group across isoforms
  source_isoform            the isoform before allele labelling
  phase_set                 the phase set
  allele_group              the allele group
  allele_labeled_isoform    the isoform as labelled with its allele
  aaseq_id                  the protein it translates to, empty when it has none
  tumor_counts              reads supporting it in the tumor BAM
  normal_counts             reads in the normal BAM, written only when one was
                            given

See aaseq_tsv for the sequence behind aaseq_id.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('gene', 'source_isoform', 'phase_set', 'allele_group',
           'allele_labeled_isoform', 'aaseq_id', 'tumor_counts')
NORMAL_COLUMN = 'normal_counts'

def columns(with_normal):
    return list(COLUMNS) + ([NORMAL_COLUMN] if with_normal else [])

class IsoAlleleCountsReader(TsvReader):
    "the normal counts column is there only when the run had a normal BAM"
    def __init__(self, iso_allele_counts_tsv):
        super().__init__(iso_allele_counts_tsv, defaultColType=str)

    @property
    def with_normal(self):
        return NORMAL_COLUMN in self.columns

class IsoAlleleCountsWriter(TsvWriter):
    def __init__(self, iso_allele_counts_tsv, *, with_normal, outFh=None):
        super().__init__(iso_allele_counts_tsv, columns=columns(with_normal),
                         defaultColType=str, outFh=outFh)

@contextmanager
def iso_allele_counts_writer(iso_allele_counts_tsv, *, with_normal):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(iso_allele_counts_tsv)
    with fileOps.AtomicFileCreate(iso_allele_counts_tsv) as tmp_tsv:
        with IsoAlleleCountsWriter(tmp_tsv, with_normal=with_normal) as writer:
            yield writer
