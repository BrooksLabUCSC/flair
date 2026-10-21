"""The amino acid sequence behind each FLP id.

Columns, in order:

  aaseq_id     the FLP id used in the aaseq_id column of an isoforms BED
  aaseq        the protein sequence itself

One row per distinct sequence; isoforms that translate alike share an id, so the
sequence is written once however many isoforms carry it.  flair transcriptome
writes <output>.aaseq.tsv and flair alleles writes
<output>.aaseq.allelegroups.tsv.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('aaseq_id', 'aaseq')

class AaSeqReader(TsvReader):
    def __init__(self, aaseq_tsv):
        super().__init__(aaseq_tsv, defaultColType=str)

class AaSeqWriter(TsvWriter):
    def __init__(self, aaseq_tsv):
        super().__init__(aaseq_tsv, columns=COLUMNS, defaultColType=str)

@contextmanager
def aaseq_writer(aaseq_tsv):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(aaseq_tsv)
    with fileOps.AtomicFileCreate(aaseq_tsv) as tmp_tsv:
        with AaSeqWriter(tmp_tsv) as writer:
            yield writer

def write_aaseqs(aaseq_tsv, aaseq_to_id):
    "aaseq_to_id maps each sequence to its id, as the callers build it"
    with aaseq_writer(aaseq_tsv) as writer:
        for aaseq, aaseq_id in aaseq_to_id.items():
            writer.writeRow((aaseq_id, aaseq))
