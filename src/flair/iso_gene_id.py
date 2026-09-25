"""The composite isoform_gene id that flair writes as a BED name and as a counts
matrix row id.

flair_transcriptome builds it as transcript_id + '_' + gene_id, and '_' occurs in
gene ids of its own accord, so the id cannot be split on the separator alone.  The
split below looks for the accession prefixes flair has met, longest match last,
and falls back to the final '_'.  It is a heuristic and it is wrong for a gene id
that neither starts with one of those prefixes nor is free of underscores, such as
a GENCODE PAR id or a gene name.

There were three copies of this before, in bed_to_gtf, diff_iso_usage and es_as,
and they had already diverged: two had an '_R2_' branch and one did not.
"""

# accession prefixes a gene id may start with, each preceded by '_' in the composite id
_GENE_ID_PREFIXES = ('_chr', '_XM', '_XR', '_NM', '_NR', '_R2_')


def split_iso_gene(iso_gene):
    "the transcript and gene halves of a composite id"
    for prefix in _GENE_ID_PREFIXES:
        if prefix in iso_gene:
            at = iso_gene.rfind(prefix)
            return iso_gene[:at], iso_gene[at + 1:]
    at = iso_gene.rfind('_')
    return iso_gene[:at], iso_gene[at + 1:]


def parse_gene_id(iso_gene):
    "the gene half of a composite id"
    return split_iso_gene(iso_gene)[1]
