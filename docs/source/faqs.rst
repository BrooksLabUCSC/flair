FAQ
===

1. Where did flair align, flair correct and flair collapse go?
--------------------------------------------------------------
``flair transcriptome`` does what correct and collapse did, starting from a sorted,
indexed BAM of the reads aligned to the genome, which you make with minimap2 yourself
in place of ``flair align``.  :ref:`migration-3.0-3.1-label` lists every option that
was renamed, split or removed, and what replaces it.

2. flair transcriptome uses too much memory, what can I do?
-----------------------------------------------------------
A run holds the parsed annotation once, and each worker holds the reads of the cluster
of overlapping reads it is building transcripts for.

- lower ``-t``, ``--threads``: the run forks that many workers, each working on a
  cluster of its own
- give ``--region`` or ``--region_bed`` to build transcripts in the regions you care
  about: only their clusters are run, and only their chromosomes' annotation is parsed

Splitting the input and running the pieces separately is no longer anything to do: a
run already works one cluster at a time, reading that cluster's reads from the BAM.

3. flair transcriptome reports single exon isoforms inside introns of genes.
-----------------------------------------------------------------------------
Most of those reads are oligo-dT primed on a run of A in the genome rather than on a
transcript's poly(A) tail.  Such reads are removed by default, and
``--keep_internally_primed`` keeps them, for libraries whose reads were not primed on
poly(A).

What survives still has to clear two bars:

- ``--single_exon_support`` reads at the same clustered ends, 3 by default
- a poly(A) consensus across the group giving it a strand, a group that gets none
  being dropped; a single exon read's alignment carries no strand, unless
  ``--trust_strand`` or ``--directRNA`` says the library is directional

``--annotation_reliant`` is gone.  Reads are aligned to the annotation's transcripts
by default when a GTF is given with ``-f``, and ``--no_align_to_annot`` turns that
off; the annotation is what names a gene, not what restricts the output to one.
