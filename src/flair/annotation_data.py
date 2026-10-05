"""External gene annotation data loaded from GTF.

Provides indexed lookups for junction-to-gene mapping, gene strand,
transcript exon structure, and single-exon gene tracking.  Used by
flair_transcriptome and flair_spliceevents for read correction, gene
assignment, and isoform filtering.
"""

from bisect import bisect_left
from flair.isoform_data import Exon, exons_to_juncs


class AnnotData(object):
    def __init__(self):
        # map of (transcript_id, gene_id) -> tuple of Exon
        self.transcript_to_exons = {}

        # list of (transcript_id, gene_id, strand)
        self.transcripts = []

        # map of junction chain tuple -> (transcript_id, gene_id)
        self.juncchain_to_transcript = {}

        # map of Junc -> set of (transcript_id, gene_id)
        self.junc_to_gene = {}

        # map of splice site position (a junction start or end) -> set of gene_id
        self.splice_site_to_genes = {}

        # single-exon annotations by strand: {'+': [], '-': []}
        # each entry is Exon(start, end, gene_id), sorted for binary search
        # FIXME: rename once it is figured out how this works in get_single_exon_gene_overlaps
        self.all_annot_SE = {'+': [], '-': []}

        # map of strand -> gene_id -> set of Exon
        # FIXME: why is strand needed here
        self.spliced_exons = {'+': {}, '-': {}}

        # map of gene_id -> set of Junc
        self.gene_to_annot_juncs = {}

        # map of gene_id -> strand
        self.gene_to_strand = {}

        # map of gene_id -> tuple of sorted exon coordinate tuples
        # union of all exons across all transcripts in the gene
        self.gene_to_exons = {}

        # map of junction chain tuple -> gene_id
        self.sjc_to_gene = {}

        # map of "transcript_id_gene_id" -> junction chain tuple
        # key format matches BED name field from gtf_to_bed --include_gene
        self.transcript_to_sjc = {}

        # map of Junc -> gene_id (single value, last gene seen wins)
        # used by spliceevents for simple junction-to-gene lookup
        self.junc_to_gene_id = {}

        self.gene_to_cds_starts = {}

        self.transcript_to_nmd_except = {}

        self.start_codon_count = 0

        # annotated transcripts' genomic start and end, for keeping subset isoforms
        # that match them: strand -> sorted list of (start, end), only for
        # transcripts whose tags don't say an end wasn't found
        self.transcript_ends = {'+': [], '-': []}
        # each side's confirmed ends on its own, for a subset isoform truncated on
        # one side: strand -> sorted list of genomic starts, and of ends
        self.confirmed_starts = {'+': [], '-': []}
        self.confirmed_ends = {'+': [], '-': []}

    def has_transcript_ends(self, strand, start, end, window):
        """does an annotated transcript on strand, with both ends confirmed, start
        and end within window of start and end"""
        ends = self.transcript_ends[strand]
        i = bisect_left(ends, (start - window,))
        while i < len(ends) and ends[i][0] <= start + window:
            if abs(ends[i][1] - end) <= window:
                return True
            i += 1
        return False

    def has_transcript_end(self, strand, side, pos, window):
        """does an annotated transcript on strand have a confirmed genomic start
        (side 0) or end (side 1) within window of pos"""
        positions = (self.confirmed_starts if side == 0 else self.confirmed_ends)[strand]
        i = bisect_left(positions, pos - window)
        return i < len(positions) and positions[i] <= pos + window


def annot_data_from_gtf(gtf_data, region):
    """Build AnnotData for a region from a pre-partitioned GtfData object."""
    annots = AnnotData()
    if gtf_data is None:
        return annots
    region_map = {region: annots}
    for trans in gtf_data.transcripts:
        if len(trans.exons) > 0:
            _process_transcript(annots, region, region_map, trans)
    # finalize gene_to_exons as sorted tuples
    for gene_id in annots.gene_to_exons:
        annots.gene_to_exons[gene_id] = tuple(sorted(annots.gene_to_exons[gene_id]))
    # once, not once per transcript, which made annotation loading quadratic.  The
    # binary search over these needs them sorted
    for se_strand in ('+', '-'):
        annots.all_annot_SE[se_strand] = sorted(annots.all_annot_SE[se_strand])
        annots.transcript_ends[se_strand].sort()
        annots.confirmed_starts[se_strand].sort()
        annots.confirmed_ends[se_strand].sort()
    return annots

def _process_transcript(annots, region, region_map, trans):
    exons = [Exon(exon.start, exon.end) for exon in trans.exons]
    sorted_exons = sorted(exons)
    t_start = sorted_exons[0].start
    t_end = sorted_exons[-1].end
    _save_transcript_annot(trans.transcript_id, trans.gene_id, region,
                           region_map, t_start, t_end, trans.strand, sorted_exons,
                           trans.attrs['tag'], trans.start_codon)

def _save_cds_starts(gene_id, start_codon, strand, annots):
    if gene_id not in annots.gene_to_cds_starts:
        annots.gene_to_cds_starts[gene_id] = set()
    if start_codon is not None:
        annots.start_codon_count += 1
        if strand == '+':
            annots.gene_to_cds_starts[gene_id].add(start_codon.start)
        else:
            annots.gene_to_cds_starts[gene_id].add(start_codon.end)

def _save_spliced_transcript_info(gene_id, t_exons, juncs, transcript_id, strand, annots):
    if gene_id not in annots.spliced_exons[strand]:
        annots.spliced_exons[strand][gene_id] = set()
    annots.spliced_exons[strand][gene_id].update(set(t_exons))
    annots.juncchain_to_transcript[tuple(juncs)] = (transcript_id, gene_id)
    annots.sjc_to_gene[tuple(juncs)] = gene_id
    annots.transcript_to_sjc[f"{transcript_id}_{gene_id}"] = tuple(juncs)
    if gene_id not in annots.gene_to_annot_juncs:
        annots.gene_to_annot_juncs[gene_id] = set()
    for j in juncs:
        if j not in annots.junc_to_gene:
            annots.junc_to_gene[j] = set()
        annots.junc_to_gene[j].add((transcript_id, gene_id))
        annots.junc_to_gene_id[j] = gene_id
        annots.gene_to_annot_juncs[gene_id].add(j)
        for site in (j.start, j.end):
            if site not in annots.splice_site_to_genes:
                annots.splice_site_to_genes[site] = set()
            annots.splice_site_to_genes[site].add(gene_id)


# GENCODE tags for a transcript end that couldn't be confirmed, the 5' (start)
# and 3' (end) of the mRNA
_MRNA_START_NOT_FOUND_TAG = 'mRNA_start_NF'
_MRNA_END_NOT_FOUND_TAG = 'mRNA_end_NF'


def _save_transcript_ends(annots, strand, t_start, t_end, transcript_tags):
    """record the transcript's ends, unless its tags say either wasn't found;
    a subset isoform is only kept for matching both ends of one transcript"""
    five_prime_found = _MRNA_START_NOT_FOUND_TAG not in transcript_tags
    three_prime_found = _MRNA_END_NOT_FOUND_TAG not in transcript_tags
    if five_prime_found and three_prime_found:
        annots.transcript_ends[strand].append((t_start, t_end))
    # the genomic start is the 5' end on +, the 3' end on -
    if (five_prime_found if strand == '+' else three_prime_found):
        annots.confirmed_starts[strand].append(t_start)
    if (three_prime_found if strand == '+' else five_prime_found):
        annots.confirmed_ends[strand].append(t_end)

def _save_transcript_annot(transcript_id, gene_id, region, region_map, t_start, t_end,
                           strand, t_exons, transcript_tags, start_codon):
    # region is a SeqRegion object
    annots = region_map[region]

    _save_cds_starts(gene_id, start_codon, strand, annots)
    annots.transcript_to_nmd_except[transcript_id] = False
    if 'NMD_exception' in transcript_tags:
        annots.transcript_to_nmd_except[transcript_id] = True

    _save_transcript_ends(annots, strand, t_start, t_end, transcript_tags)
    annots.transcript_to_exons[(transcript_id, gene_id)] = tuple(t_exons)
    juncs = exons_to_juncs(t_exons)
    annots.transcripts.append((transcript_id, gene_id, strand))
    if gene_id not in annots.gene_to_strand:
        annots.gene_to_strand[gene_id] = strand
    # accumulate exons per gene (as coordinate tuples for spliceevents compatibility)
    exon_coords = set((e.start, e.end) for e in t_exons)
    if gene_id not in annots.gene_to_exons:
        annots.gene_to_exons[gene_id] = exon_coords
    else:
        annots.gene_to_exons[gene_id].update(exon_coords)
    if len(juncs) == 0:
        annots.all_annot_SE[strand].append(Exon(t_start, t_end, gene_id))
    else:
        _save_spliced_transcript_info(gene_id, t_exons, juncs, transcript_id, strand, annots)
