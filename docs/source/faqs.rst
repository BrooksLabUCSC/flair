FAQ
===

1. Where did flair align, flair correct and flair collapse go?
--------------------------------------------------------------
``flair transcriptome`` does what correct and collapse did, starting from a sorted,
indexed BAM of the reads aligned to the genome, made with minimap2 in place of
``flair align``.  :ref:`migration-3.0-3.1-label` lists every option that
was renamed, split or removed, and what replaces it.

2. flair transcriptome uses too much memory, what can I do?
-----------------------------------------------------------
A run keeps the parsed annotation once, and each worker keeps the reads of the cluster
of overlapping reads it is building transcripts for.

- lower ``-t``, ``--threads``: the run forks that many workers, each working on a
  cluster of its own
- give ``--region`` or ``--region_bed`` to build transcripts in the regions of
  interest: only their clusters are run, and only their chromosomes' annotation is
  parsed

Splitting the input and running the pieces separately is no longer needed: a run
already works one cluster at a time, reading that cluster's reads from the BAM.

3. flair transcriptome reports single exon isoforms inside introns of genes.
-----------------------------------------------------------------------------
Most of those reads are oligo-dT primed on the genome's own A rather than on a
transcript's poly(A) tail: 60% A in the 20 genome bases after the read's 3' end with
no poly(A) tail clipped there, or 80% A in its last 20 aligned bases.  Such reads are
removed by default, and ``--keep_internally_primed`` keeps them, for libraries whose
reads were not primed on poly(A).

A read that survives must also pass two tests:

- ``--single_exon_support`` reads at the same clustered ends, 3 by default
- a poly(A) consensus across the group giving it a strand, a group that gets none
  being dropped; a single exon read's alignment has no strand, unless
  ``--trust_strand`` or ``--directRNA`` states that the library is directional

``--annotation_reliant`` is gone.  Reads are aligned to the annotation's transcripts
by default when a GTF is given with ``-f``, and ``--no_align_to_annot`` turns that
off; the annotation is what names a gene, not what restricts the output to one.

4. Do my reads need unique ids?
--------------------------------
Within each BAM, yes: an id must name one read, and FLAIR stops with an error naming a
repeated one.  Between samples ids may repeat, since each sample's BAM is counted on its
own.  Running ``flair transcriptome`` on each replicate and combining the transcriptomes
requires nothing of the ids; merging replicates into one BAM first means labelling
each run's reads.  :ref:`alignment-label` has the rule and a command that labels them.

5. Which file from pbmm2 or dorado do I give FLAIR?
----------------------------------------------------
The reads aligned to the genome, coordinate sorted, with its index beside it, which is
``pbmm2 align --sort`` output or what ``dorado aligner`` writes to ``--output-dir``.  The
basecaller's reads file contains no alignments, and a transcriptome aligned BAM has a
transcript's coordinates rather than the genome's.  :ref:`alignment-label` names the
files and the options.
