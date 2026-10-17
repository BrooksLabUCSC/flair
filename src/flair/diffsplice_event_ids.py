"""The ids naming an alternative splicing event in the diffsplice
events.quant.tsv files, used in both the feature_id and the coordinate column.

The strand is part of every id.  DRIMSeq requires feature_id to be unique and
takes coordinate as the gene_id it groups by, and an event at the same
coordinates on the two strands is two events, not one.  The form,
chrom:start-end(strand), is the one flair_spliceevents and intron_support use.
"""

def site_id(chrom, strand, site):
    "a splice site as chrom:pos(strand)"
    return f'{chrom}:{site}({strand})'

def junction_id(chrom, strand, start, end):
    "a junction, or any other interval, as chrom:start-end(strand)"
    return f'{chrom}:{start}-{end}({strand})'
