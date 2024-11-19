import sys, subprocess, os





###DO THIS IN TERMINAL FIRST
#export PATH="<PATHTOFLAIR>/bin:<PATHTOFLAIR>/src/flair:$PATH"
flair_path = '/private/groups/brookslab/cafelton/git-flair/flair/'
output_dir = '/private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/'
reference_genome = '/private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa'
reference_annot = '/private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf'
outprefix = 'isofusionvars082824'
sampleprefix = '082824-'

for i in range(1,4):
    outprefix = output_dir + outprefix + '-' + str(i) + '_flair/'
    infa = output_dir + outprefix + '-' + str(i) + '-10x.fa'
    if not os.path.exists(outprefix):
        os.mkdir(outprefix)
    aligncommand = 'python3 ' + flair_path + 'flair.py align -r ' + infa + \
                   ' -t 12 -g ' + reference_genome + ' -f separate --minfragmentsize 40 -o ' + \
                    outprefix + 'flair.align'
    # print(aligncommand)
    subprocess.call(aligncommand, shell=True)
    correctcommand = 'python3 ' + flair_path + 'flair.py correct -q ' + outprefix + 'flair.align.bed' + \
                     ' -t 12 -g ' + reference_genome + ' ' + \
                     '-f ' + reference_annot + ' -o ' + \
                      outprefix + 'flair'
    # print(correctcommand)
    subprocess.call(correctcommand, shell=True)
    collapsecommand = 'python3 ' + flair_path + 'flair.py collapse -q ' + outprefix + 'flair_all_corrected.bed' + \
                     ' -t 12 -g ' + reference_genome + ' ' + \
                     '-f ' + reference_annot + ' -o ' + \
                    outprefix + 'flair.collapse -r ' + infa + ' --generate_map --annotation_reliant generate ' + \
                    '--stringent --check_splice --quality 0 --isoformtss -n longest'
    # print(collapsecommand)
    subprocess.call(collapsecommand, shell=True)

    fusioncommand = 'python3 ' + flair_path + 'src/flair/flair_detectfusions.py' + \
                         ' -t 12 -g ' + reference_genome + ' ' + \
                     '-f ' + reference_annot + ' -o ' + \
                    outprefix + 'flair.fusion -r ' + infa + ' -b ' + outprefix + 'flair.align_chimeric.bam' + \
                    ' --annotated_fa ' + outprefix + 'flair.collapse.annotated_transcripts.fa'
    # print(fusioncommand)
    subprocess.call(fusioncommand, shell=True)

out = open(output_dir + outprefix + '_combinemanifest.tsv', 'w')
for i in range(1,4):
    outprefix = output_dir + outprefix + '-' + str(i) + '_flair/'
    #WTC11.ENCFF245IPA	/private/groups/brookslab/tswon/flair_transcriptome_optimization/samples/WTC11.ENCFF245IPA-withshortread.isoforms.bed	/private/groups/brookslab/tswon/flair_transcriptome_optimization/samples/WTC11.ENCFF245IPA-withshortread.combined.isoform.read.map.txt
    outline = [sampleprefix + str(i), 'isoform', outprefix + 'flair.collapse.isoforms.bed', outprefix + 'flair.collapse.isoforms.fa', outprefix + 'flair.collapse.combined.isoform.read.map.txt']
    out.write('\t'.join(outline) + '\n')
    outline = [sampleprefix + str(i), 'fusionisoform', outprefix + 'flair.fusion.fusions.isoforms.bed', outprefix + 'flair.fusion.syntheticAligned.isoforms.fa', outprefix + 'flair.fusion.syntheticAligned.isoform.read.map.txt']
    out.write('\t'.join(outline) + '\n')
out.close()

combinecommand = 'python3 ' + flair_path + 'flair.py combine -m ' + output_dir + outprefix + '_combinemanifest.tsv -o ' + output_dir + outprefix + '.flair.combined.isoforms'
subprocess.call(combinecommand, shell=True)

out = open(output_dir + outprefix + '_quantifymanifest.tsv', 'w')
for i in range(1,4):
    outprefix = output_dir + outprefix + '-' + str(i) + '_flair/'
    infa = output_dir + outprefix + '-' + str(i) + '-10x.fa'
    outline = [sampleprefix + str(i), 'c1', 'b1', infa]
    out.write('\t'.join(outline) + '\n')
out.close()

quantifycommand = 'python3 ' + flair_path + 'flair.py quantify ' + \
                  '--reads_manifest ' + output_dir + outprefix + '_quantifymanifest.tsv ' + \
                  '-i ' + output_dir + outprefix + '.flair.combined.isoforms.fa ' + \
                  '-o ' + output_dir + outprefix + '.flair.quantify ' + \
                  '--output_bam --quality 0 --threads 12'
##Don't use this version
# quantifycommand = 'python3 ' + flair_path + 'flair.py quantify ' + \
#                   '--reads_manifest /private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/isofusionvars082824_quantifymanifest.tsv ' + \
#                   '-i /private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/isofusionvars082824.flair.combined.isoforms.fa ' + \
#                   '--isoform_bed /private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/isofusionvars082824.flair.combined.isoforms.bed ' + \
#                   '-o /private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/isofusionvars082824.flair.quantify ' + \
#                   '--generate_map --output_bam --quality 0 --stringent --check_splice --threads 12'
subprocess.call(quantifycommand, shell=True)


out = open(output_dir + outprefix + '_variantsmanifest.tsv', 'w')
for i in range(1,4):
    outprefix = output_dir + outprefix + '-' + str(i) + '_flair/'
    infa = output_dir + outprefix + '.flair.quantify.' + sampleprefix + str(i) + '.c1.flair.aligned.bam'
    outline = [sampleprefix + str(i), infa]
    out.write('\t'.join(outline) + '\n')
out.close()

variantscommand = 'python3 ' + flair_path + 'src/flair/flair_variants.py -m ' + output_dir + outprefix + '_variantsmanifest.tsv ' + \
                  '-i ' + output_dir + outprefix + '.flair.combined.isoforms.fa ' + \
                  '-b ' + output_dir + outprefix + '.flair.combined.isoforms.bed ' + \
                  '-o ' + output_dir + outprefix + ' ' + \
                  '-g ' + reference_genome + ' ' + \
                  '-f  ' + reference_annot + ' '
subprocess.call(variantscommand, shell=True)
