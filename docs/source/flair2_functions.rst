FLAIR2 capabilities
==================

FLAIR2 has an updated isoform detection algorithm and an added feature of variant-aware isoform detection.


Performance increases
-------------------------------

FLAIR2 run with the ``--annotation_reliant`` argument invokes an alignment of the reads to an annotated transcriptome first, followed by novel isoform detection. This can be run with or without ``--check_splice``, which enforces higher quality matching specifically around each splice site for read-to-isoform assignment steps.

.. code:: sh


   flair collapse --check_splice --annotation_reliant generate -f annotation.gtf -g genome.fa -r reads.fa -q corrected.bed [options]

If you are running collapse with the same transcript reference multiple times, you can specify the previously generated transcript sequence file to the ``--annotation_reliant`` argument instead.

Variant integration
-------------------------------

FLAIR has two modalities for phasing variants to discover variant-aware transcript models. The first uses phasing information from longshot, which is comprised of a phase set determined for each read and a set of variants corresponding to each phase set. For the second modality, FLAIR can approach phasing variants that is agnostic to ploidy, which may be worthy of exploration if working with RNA edits and potential cancer-related aneuploidies: 1) given variant calls, FLAIR tabulates the most frequent combinations of variants present in each isoform from the supporting read sequences; 2) from the isoform-defining collapse step, FLAIR generates a set of reads assigned to each isoform; so 3) isoforms that have sufficient read support for a collection of mismatches are determined. This latter method of phasing focuses solely on frequency of groups of mismatches that co-occur within reads and does not use ploidy information to refine haplotypes, allowing for the generation of multiple haplotypes within a gene and transcript model. 


Longshot
~~~~~~~~~~~~

Longshot provides phased read outputs, which can be supplied to flair-collapse via ``--longshot_vcf`` and ``--longshot_bam``. The outputs of collapse are the following: 1) isoform models as a bed file, 2) the subset of variants from the longshot vcf that were used, and 3) isoform sequences with variants as a fasta file. The isoform models and variants can be viewed by aligning the isoform sequences and using IGV or other visualization tools.
.. code:: sh


   longshot --bam flair.aligned.bam --ref genome.fa --out longshot.vcf --out_bam longshot.bam --min_allele_qual 3 -F
   samtools index longshot.bam

   flair collapse -r reads.fa -q corrected.bed -g genome.fa --longshot_vcf longshot.vcf --longshot_bam flair.longshot.bam [options]

   minimap2 -ax splice --secondary=no genome.fa flair.collapse.isoforms.fa > flair.collapse.isoforms.fa.sam
   samtools sort flair.collapse.isoforms.fa.sam -o flair.collapse.isoforms.fa.bam
   samtools index flair.collapse.isoforms.fa.bam


Any vcf
~~~~~~~~~~~~

FLAIR2 can also take a vcf agnostic to the variant caller and spike variants in given any isoform model file. If enough supporting reads for an individual isoform contain the same pattern of variants, then FLAIR will create an additional isoform with _PSX appended to the name. Flair-collapse needs to be run with ``--generate_map``.
.. code:: sh


   flair collapse -r reads.fa -q corrected.bed -g genome.fa --generate_map [options]

   assign_variants_to_transcripts --bam flair.aligned.bam -i flair.collapse.isoforms.bed -v variants.vcf --map flair.collapse.isoform.read.map.txt --bed_out out.bed --map_out out.map > out.vcf 

   psl_to_sequence out.bed ~/bl/Brooks_LRGASP/reference/lrgasp_grcm39_sirvs.fasta out.fa -v out.vcf

   minimap2 -ax splice --secondary=no genome.fa out.fa > out.fa.sam
   samtools sort out.fa.sam -o out.fa.bam
   samtools index out.fa.bam


