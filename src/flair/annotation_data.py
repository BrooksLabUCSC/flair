"""External gene annotation data loaded from GTF.

Provides indexed lookups for junction-to-gene mapping, gene strand,
transcript exon structure, and single-exon gene tracking.  Used by
flair_transcriptome and flair_spliceevents for read correction, gene
assignment, and isoform filtering.
"""

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
    return annots

def _process_transcript(annots, region, region_map, trans):
    trans.gene_id = trans.gene_id.replace('_', '-')
    exons = [Exon(exon.start, exon.end) for exon in trans.exons]
    sorted_exons = sorted(exons)
    t_start = sorted_exons[0].start
    t_end = sorted_exons[-1].end
    _save_transcript_annot(trans.transcript_id, trans.gene_id, region,
                           region_map, t_start, t_end, trans.strand, sorted_exons)


def _save_transcript_annot(transcript_id, gene_id, region, region_map, t_start, t_end, strand, t_exons):
    # region is a SeqRegion object
    annots = region_map[region]
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
    for strand in ['+', '-']:
        annots.all_annot_SE[strand] = sorted(annots.all_annot_SE[strand])  # FIXME: make set? Colette note: needs to be sorted for binary search later
