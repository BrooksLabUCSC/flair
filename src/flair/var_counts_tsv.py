"""Reads showing each allele of each site, from flair variantquant, a row per allele of
a site and group of its reads.

Columns, in order:

  varpos       the site, chrom:position
  ref          the reference base, the genome's
  alt          the allele counted, one of the site's alts; a site of several alts
               has a row each
  strand       the strand of the reads, with --trust_strand, otherwise '.'
  phase_set    the phase set of the reads, with --haplotypes, otherwise '.'
  haplotype    their haplotype in it, with --haplotypes, otherwise '.'; the reads
               of no haplotype have a row of '.' too
  gene         the gene of the exons over the site, the first of the sorted list
               when there are more than one, empty when there are none
  <sample>     one column per sample, holding the reads covering the site not
               showing the alt and those showing it, ';' separated in the file and a
               pair on the row, so the two add up to the site's coverage whatever its
               other alts

Positions are one based, as the VCF they came from writes them.  Every run has
these columns, whether or not it splits the reads by strand or haplotype.

An allele's rows are written when a read of some sample shows the alt, or always
with --output_all, each when some sample has at least --min_coverage reads of its
group covering the site.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import read_header_columns, semicolonIntListType

ID_COLUMNS = ('varpos', 'ref', 'alt', 'strand', 'phase_set', 'haplotype', 'gene')

def columns(sample_columns):
    return list(ID_COLUMNS) + list(sample_columns)

class VarCountsReader(TsvReader):
    "every sample column holds the reads not showing and showing the alt, as a pair"
    def __init__(self, var_counts_tsv):
        sample_columns = read_header_columns(var_counts_tsv)[len(ID_COLUMNS):]
        super().__init__(var_counts_tsv, defaultColType=str,
                         typeMap={col: semicolonIntListType for col in sample_columns})

    @property
    def sample_columns(self):
        return self.columns[len(ID_COLUMNS):]

class VarCountsWriter(TsvWriter):
    def __init__(self, var_counts_tsv, sample_columns, *, outFh=None):
        super().__init__(var_counts_tsv, columns=columns(sample_columns),
                         defaultColType=str, outFh=outFh,
                         typeMap={col: semicolonIntListType for col in sample_columns})

@contextmanager
def var_counts_writer(var_counts_tsv, sample_columns):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(var_counts_tsv)
    with fileOps.AtomicFileCreate(var_counts_tsv) as tmp_tsv:
        with VarCountsWriter(tmp_tsv, sample_columns) as writer:
            yield writer
