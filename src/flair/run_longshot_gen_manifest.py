import sys
import os
import pipettor


manifest = sys.argv[1]
refisoforms = sys.argv[2]
outprefix = sys.argv[3]
outmanifestfile = sys.argv[4]

outprefix = outprefix.split('/')[-1]

if not os.path.exists(refisoforms + '.fai'):
    pipettor.run([('samtools', 'faidx', refisoforms)])

outmanifest = open(outmanifestfile, 'w')
for line in open(manifest):
    sample, bamfile = line.rstrip().split('\t')

    ##run longshot on each aligned bam file
    print('running longshot on', sample)
    longshotcmd = ('longshot', '--bam', bamfile, '--ref', refisoforms, '--out', outprefix + '.' +  sample + '.flairaligned.vcf', '-F', '-P', '0.00001', '-q', '10', '--output-ref'  )  # , '-d' '/private/groups/brookslab/cafelton/testflairanyvcf/simisofusionvars/longshotdebug')
    pipettor.run([longshotcmd])
    outmanifest.write('\t'.join([sample, bamfile, outprefix + '.' +  sample + '.flairaligned.vcf']) + '\n')
outmanifest.close()
