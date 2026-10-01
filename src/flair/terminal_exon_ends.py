"""Terminal exon end normalization, shared by flair transcriptome and flair quantify.

The first and last exons of each spliced isoform in a gene are extended to the
furthest end seen among the gene's isoforms whose terminal exon has a splice site
within NORM_END_SS_WINDOW of theirs, plus NORM_END_EXTRA_LEN of padding so reads
that run a little past the end still align.  Single-exon isoforms are not
normalized.  The padded ends are clamped to the chromosome.  Everything here works
on plain coordinates, so it does not depend on how a caller represents exons.
"""
from flair import FlairInputDataError

# padding added beyond the furthest terminal exon end
NORM_END_EXTRA_LEN = 100
# terminal exons with splice sites this close are grouped together
NORM_END_SS_WINDOW = 50


class GeneTerminalExonEnds:
    """Furthest terminal exon ends of one gene, grouped by the splice site at the
    inner edge of the terminal exon."""
    def __init__(self, gene_id):
        self.gene_id = gene_id
        self.first_exon_starts = {}  # first exon splice site (its end) -> min start
        self.last_exon_ends = {}  # last exon splice site (its start) -> max end

    def add(self, first_start, first_ss, last_ss, last_end):
        if first_ss not in self.first_exon_starts:
            self.first_exon_starts[first_ss] = first_start
        else:
            self.first_exon_starts[first_ss] = min(self.first_exon_starts[first_ss], first_start)
        if last_ss not in self.last_exon_ends:
            self.last_exon_ends[last_ss] = last_end
        else:
            self.last_exon_ends[last_ss] = max(self.last_exon_ends[last_ss], last_end)

    def first_exon_start(self, first_ss):
        "min start of first exons with a splice site within the window of first_ss"
        return min(self.first_exon_starts[i]
                   for i in range(first_ss - NORM_END_SS_WINDOW, first_ss + NORM_END_SS_WINDOW)
                   if i in self.first_exon_starts)

    def last_exon_end(self, last_ss):
        "max end of last exons with a splice site within the window of last_ss"
        return max(self.last_exon_ends[i]
                   for i in range(last_ss - NORM_END_SS_WINDOW, last_ss + NORM_END_SS_WINDOW)
                   if i in self.last_exon_ends)

    def furthest_ends(self, first_ss, last_ss):
        """return the unpadded (start, end) for a spliced isoform whose first exon ends
        at first_ss and whose last exon starts at last_ss"""
        return self.first_exon_start(first_ss), self.last_exon_end(last_ss)

    def normalized_ends(self, first_ss, last_ss, chrom_len):
        """return the furthest_ends padded by NORM_END_EXTRA_LEN and clamped to
        [0, chrom_len]"""
        start, end = self.furthest_ends(first_ss, last_ss)
        return max(start - NORM_END_EXTRA_LEN, 0), min(end + NORM_END_EXTRA_LEN, chrom_len)


class TerminalExonEnds:
    """Furthest terminal exon ends by gene."""
    def __init__(self):
        # FIXME: this is temporary.  The code groups by (gene_id, strand)
        # for reasons that are suspected to be bugs in stranding.
        self._by_gene = {}  # (gene_id, strand) -> GeneTerminalExonEnds
        self._gene_id_to_strand = {}

    def _obtain(self, gene_id, strand):
        "get current entry or create a new one"
        gene_key = (gene_id, strand)
        gene_entry = self._by_gene.get(gene_key)
        if gene_entry is None:
            gene_entry = GeneTerminalExonEnds(gene_id)
            self._by_gene[gene_key] = gene_entry

        existing_strand = self._gene_id_to_strand.get(gene_id)
        if existing_strand is None:
            self._gene_id_to_strand[gene_id] = strand
        elif strand != existing_strand:
            raise FlairInputDataError(f"gene id '{gene_id}' has transcripts on both strands, "
                                      f"'{existing_strand}' and '{strand}'; give each strand its own "
                                      "gene id")
        return gene_entry

    def add_transcript(self, gene_id, strand, exon_bounds):
        """add a transcript given its exons as (start, end) pairs in genomic order;
        single-exon transcripts only record the gene"""
        gene_entry = self._obtain(gene_id, strand)
        if len(exon_bounds) > 1:
            (first_start, first_ss), (last_ss, last_end) = exon_bounds[0], exon_bounds[-1]
            gene_entry.add(first_start, first_ss, last_ss, last_end)

    def fetch(self, gene_id, strand) -> GeneTerminalExonEnds:
        """return entry or error"""
        return self._by_gene[(gene_id, strand)]

    def furthest_ends(self, gene_id, strand, exon_bounds):
        """return the unpadded (start, end) for a spliced transcript in the gene, given
        its exons as (start, end) pairs in genomic order"""
        return self.fetch(gene_id, strand).furthest_ends(exon_bounds[0][1], exon_bounds[-1][0])

    def normalized_ends(self, gene_id, strand, exon_bounds, chrom_len):
        """return the padded, clamped (start, end) for a spliced transcript in the gene,
        given its exons as (start, end) pairs in genomic order"""
        return self.fetch(gene_id, strand).normalized_ends(exon_bounds[0][1], exon_bounds[-1][0], chrom_len)
