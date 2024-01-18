FAQ
===

1. Flair collapse uses too much memory, what can I do?
------------------------------------------------------
Flair's memory requirements increase with larger input files.
If your bed file is over 1 Gigabyte, consider splitting it by chromosome
and then running separately on each file.

2. Flair collapse reports a lot of single exon isoforms inside introns of genes.
--------------------------------------------------------------------------------
Flair tries to find different isoforms of a gene, which means that it tends to have trouble with single exon transcripts that do not overlap exons of known genes. Since it has no way of knowing that these are not separate genes (especially because single exon reads often have no clear orientation) it reports them all so you don't lose information.
The best way to deal with this is to use --annotation_reliant. This restricts Flair to only genes present in the input gtf.

 
