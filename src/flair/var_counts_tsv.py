"""Reads showing each allele of each site, from flair variantquant, two rows per allele
of a site and group of its reads: the reads showing the alt, and those not.

Columns, in order:

  var_name     the allele and group of reads the two rows count, as
               chrom:position:ref>alt, then (strand) with --trust_strand, then
               |phase_set|haplotype for the reads of a phase set or haplotype
  is_var       'yes' for the row of the reads showing the alt, 'no' for the row of
               those covering the site not showing it, so the two add up to the
               site's coverage whatever its other alts
  varpos       the site, chrom:position
  ref          the reference base, the genome's
  alt          the allele counted, one of the site's alts; a site of several alts
               has rows each
  strand       the strand of the reads, with --trust_strand, otherwise '.'
  phase_set    the phase set of the reads, with --haplotypes, otherwise '.'
  haplotype    their haplotype in it, with --haplotypes, otherwise '.'; the reads
               of no haplotype have rows of '.' too
  gene         the gene of the exons over the site, the first of the sorted list
               when there are more than one, empty when there are none
  <sample>     one integer column per sample, the reads of the row

var_name and is_var are a group and a feature, as the gene_id and isoform_id of a
counts matrix, so flair diffexp and diff_iso_usage test a change in an allele's
share of its site's reads between conditions (see counts_matrix_tsv).

Positions are one based, as the VCF they came from writes them.  Every run has
these columns, whether or not it splits the reads by strand or haplotype.

An allele's rows are written when a read of some sample shows the alt, or always
with --output_all, each pair when some sample has at least --min_coverage reads of
its group covering the site.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import read_header_columns

ID_COLUMNS = ('var_name', 'is_var', 'varpos', 'ref', 'alt', 'strand', 'phase_set', 'haplotype', 'gene')
# is_var of the rows of the reads showing the alt and of those not showing it
IS_VAR, NOT_VAR = 'yes', 'no'
# a strand, phase set or haplotype the reads were not split by
NO_GROUP = '.'

def var_name(varpos, ref, alt, strand, phase_set, haplotype):
    "the name of an allele of a site and a group of its reads, as var_name"
    name = f'{varpos}:{ref}>{alt}'
    if strand != NO_GROUP:
        name += f'({strand})'
    if (phase_set, haplotype) != (NO_GROUP, NO_GROUP):
        name += f'|{phase_set}|{haplotype}'
    return name

def columns(sample_columns):
    return list(ID_COLUMNS) + list(sample_columns)

class VarCountsReader(TsvReader):
    def __init__(self, var_counts_tsv):
        sample_columns = read_header_columns(var_counts_tsv)[len(ID_COLUMNS):]
        super().__init__(var_counts_tsv, defaultColType=str,
                         typeMap={col: int for col in sample_columns})

    @property
    def sample_columns(self):
        return self.columns[len(ID_COLUMNS):]

class VarCountsWriter(TsvWriter):
    def __init__(self, var_counts_tsv, sample_columns, *, outFh=None):
        super().__init__(var_counts_tsv, columns=columns(sample_columns),
                         defaultColType=str, outFh=outFh,
                         typeMap={col: int for col in sample_columns})

    def writeAllele(self, varpos, ref, alt, group, gene, sample_counts):
        """the two rows of an allele and group of reads, (strand, phase_set,
        haplotype): sample_counts are per sample the reads covering the site not
        showing the alt and those showing it"""
        name = var_name(varpos, ref, alt, *group)
        self.writeRow([name, IS_VAR, varpos, ref, alt, *group, gene] + [alt_reads for _, alt_reads in sample_counts])
        self.writeRow([name, NOT_VAR, varpos, ref, alt, *group, gene] + [not_alt for not_alt, _ in sample_counts])

@contextmanager
def var_counts_writer(var_counts_tsv, sample_columns):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(var_counts_tsv)
    with fileOps.AtomicFileCreate(var_counts_tsv) as tmp_tsv:
        with VarCountsWriter(tmp_tsv, sample_columns) as writer:
            yield writer
