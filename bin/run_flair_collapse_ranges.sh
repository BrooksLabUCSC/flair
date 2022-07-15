GENOME_FA=genome.fa
CORRECTED_BED=all_corrected.bed
ANNOTATION_GTF=annotation.gtf
READS_BAM=aligned.bam
TEMP_DIR=flair_temp
OUTFILE_NAME=flair.collapse

# partition the bed into independent regions
bedPartition $CORRECTED_BED $CORRECTED_BED.ranges.bed

# tabix index the corrected reads bed file
sort -k1,1 -k2,2n $CORRECTED_BED > $CORRECTED_BED.sorted.bed
bgzip $CORRECTED_BED.sorted.bed
tabix -f --preset bed --zero-based $CORRECTED_BED.sorted.bed.gz

# run flair collapse for each independent region
parallel -a $CORRECTED_BED.ranges.bed python ../flair.py collapse -q $CORRECTED_BED.sorted.bed.gz -g $GENOME_FA -r $READS_BAM -f $ANNOTATION_GTF --temp_dir $TEMP_DIR -o $TEMP_DIR/temp --quiet --range

# concatenate results from all regions
cat $TEMP_DIR/temp*isoforms.bed > $OUTFILE_NAME.isoforms.bed
cat $TEMP_DIR/temp*isoforms.fa > $OUTFILE_NAME.isoforms.fa
cat $TEMP_DIR/temp*isoforms.gtf > $OUTFILE_NAME.isoforms.btf