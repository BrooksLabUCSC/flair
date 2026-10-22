"""The tables of flair isoalleles, OUTPUT.isoalleles.<table>.tsv, one per entry of
TABLES.  Its counts are a counts matrix (counts_matrix_tsv), its proteins an aaseq
table (aaseq_tsv), and its sequences and isoforms a FASTA and a flair BED.

A group is an isoform's reads of one haplotype (H1, H2, untagged) and, in a pair,
one somatic set (H2-S1), or of an unknown set (H2-S?); see isoallele_genes.  A
variant is written chrom:pos:ref>alt, pos 1-based as in a VCF.

sequences: one row per sequence, an isoform's group with the variants applied
  gene_id                       the gene
  isoform_id                    the isoform
  allele_group                  the group, as H2-S1; the sequence is isoform|group
  variants                      the variants applied, comma separated
  terminal_exon_somatic_added   the terminal-exon somatic variants also applied
  productivity                  PRO, PTC, NST or NGO, as flair transcriptome calls it
  productivity_note             why there is no ORF or no stop, if so
  aaseq_id                      the protein's id: the isoform's own, or FLPV and a hash
  protein_change                the changes from the isoform's own protein
  somatic_protein_change        the changes from the same haplotype's protein
                                without the somatic variants
  cds                           the CDS, chrom:start-end, 0-based half open
  phase_sets                    the phase sets of the isoform's phased hets

variants: one row per variant used
  gene_id, isoform_ids          the gene, and the isoforms it is in, comma separated
  variant                       the variant
  kind                          germline homozygous, germline phased het, germline
                                unphased het (placed by reads, or not placed), somatic
  haplotype_or_sets             the haplotype of a het, or the somatic sets carrying it
  detail                        a het's phase set, a somatic variant's normal reads

somatic_sets: one row per somatic set of a gene, of a pair
  gene_id, somatic_set          the gene, and the set, S1 and so on, by reads, most first
  variants                      the somatic variants the set's reads carry
  reads                         the reads, in every sample, carrying exactly them

excluded: the tumor's calls left out of a pair as somatic, too few normal reads
  gene_id, variant, normal_reads, isoform_ids, reason

warnings: one row per warning
  gene_id, isoform_id, warning, detail

terminal_sites: per isoform, its terminal-exon somatic variants, per group
  gene_id, isoform_id, variant, allele_group
  reads_reaching, reads_carrying the group's reads that inform the site, and carry it
  decision                      added, dropped, or that the site defines sets

genes: per gene and sample, how many of its spliced reads the read map assigns
  gene_id, sample
  spliced_reads                 primary reads of the gene's strand with an intron,
                                over one of its exons
  assigned_spliced_reads        those assigned to one of its isoforms
  assigned_share                their share, empty for no spliced reads
  assigned_reads                every read assigned to one of its isoforms
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter

TABLES = {
    'sequences': ('gene_id', 'isoform_id', 'allele_group', 'variants', 'terminal_exon_somatic_added',
                  'productivity', 'productivity_note', 'aaseq_id', 'protein_change', 'somatic_protein_change',
                  'cds', 'phase_sets'),
    'variants': ('gene_id', 'isoform_ids', 'variant', 'kind', 'haplotype_or_sets', 'detail'),
    'somatic_sets': ('gene_id', 'somatic_set', 'variants', 'reads'),
    'excluded': ('gene_id', 'variant', 'normal_reads', 'isoform_ids', 'reason'),
    'warnings': ('gene_id', 'isoform_id', 'warning', 'detail'),
    'terminal_sites': ('gene_id', 'isoform_id', 'variant', 'allele_group', 'reads_reaching', 'reads_carrying',
                       'decision'),
    'genes': ('gene_id', 'sample', 'spliced_reads', 'assigned_spliced_reads', 'assigned_share', 'assigned_reads'),
}


def table_path(output, table):
    return f'{output}.isoalleles.{table}.tsv'


class IsoallelesTableReader(TsvReader):
    "values stay strings: counts and shares can be empty"
    def __init__(self, path, table):
        super().__init__(path, defaultColType=str)
        if tuple(self.columns) != TABLES[table]:
            raise ValueError(f'{path}: expected the columns of the isoalleles {table} table')


class IsoallelesTableWriter(TsvWriter):
    def __init__(self, path, table, *, outFh=None):
        super().__init__(path, columns=TABLES[table], defaultColType=str, outFh=outFh)


@contextmanager
def isoalleles_table_writer(path, table):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(path)
    with fileOps.AtomicFileCreate(path) as tmp_tsv:
        with IsoallelesTableWriter(tmp_tsv, table) as writer:
            yield writer
