.. _ids-label:

FLAIR gene and isoform identifiers
==================================

FLAIR assigns its own identifiers to the genes, isoforms and protein sequences it
defines.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Identifier
     - Example
     - Description
   * - gene id
     - ``FLG00000001``
     - Identifies a gene, annotated or novel.
   * - isoform id
     - ``FLT00000001``
     - Identifies one isoform of one gene.
   * - amino acid sequence id
     - ``FLP00000001``
     - Identifies one protein sequence. Isoforms that translate to the same
       sequence share an id; an isoform with no coding prediction has none.

Each is a three letter prefix and eight digits.

The numbers are assigned in file order and are unique within one FLAIR run, not
across runs. Two runs over the same reads can give the same gene or isoform
different numbers, so an id is only meaningful alongside the output it was written
with. ``flair combine`` numbers afresh, so an id in a combined transcriptome is
unrelated to the id the same isoform had in the per-sample transcriptome it came
from.

A fusion isoform has an isoform id and a gene id of its own. The genes it joins keep
their own gene ids, listed in transcript order.

One id per column
-----------------

No FLAIR file puts an isoform and its gene in one field. The gene is a column of its
own: ``gene_id`` in an isoforms BED, ``gene_id`` beside ``isoform_id`` in a counts
matrix, and the column ``gtf_to_bed --include_gene`` adds. An id is therefore read,
never parsed.

The annotated gene or transcript an isoform was matched to, when it was matched to
one, is likewise its own column rather than part of the id.

Releases before FLAIR 3 joined the two with an underscore,
``ENST00000225792.10_ENSG00000108654.15``. Splitting that needs a guess, since gene
ids contain underscores of their own, so FLAIR no longer writes or reads it. The one
exception is ``bed_to_gtf`` given a BED that FLAIR did not write, which has no column
to read and so must guess.
