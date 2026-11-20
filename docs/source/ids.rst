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

Derived identifiers
-------------------

Three things FLAIR names are named for what they are made of, so their ids do have
parts, with a separator that the id's own component ids cannot contain:

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - What
     - Form
     - Where
   * - an allele of an isoform
     - ``FLT00000001|H1``
     - the counts matrix, BED, FASTA and read map of ``flair isoalleles``
   * - a splicing event
     - ``es|chr1:100-200(+)|ENSG00000188243.13``
     - the ``event`` column and the BED name of ``flair spliceevents``, whose BED
       name adds ``|alternative``
   * - a fusion
     - ``ERGIC3--KRAS``
     - the ``fusion`` and ``gene_ids`` columns of ``flair fusion``'s table

A row of a ``flair spliceevents`` counts table is named ``event|alternative`` where
a counts matrix names an isoform, since that is what the row counts.

Releases before FLAIR 3 joined the two with an underscore,
``ENST00000225792.10_ENSG00000108654.15``. Splitting that needs a guess, since gene
ids contain underscores of their own, so FLAIR no longer writes or reads it.
``bed_to_gtf`` rejects a BED that FLAIR did not write, which has no gene column,
unless given ``--force``, which uses each isoform name as its gene id too.
