import pysam
import sys
from flair.isoform_data import get_reverse_complement
from flair.pycbio.hgdata.bed import BedReader


def get_iso_info_from_bed(file):
    isotoinfo = {}
    for bed in BedReader(file, fixScores=True):
        chrom, start, end, iso, strand = bed.chrom, bed.chromStart, bed.chromEnd, bed.name, bed.strand
        esizes = [len(blk) for blk in bed.blocks]
        estarts = [blk.start - start for blk in bed.blocks]
        isizes = [estarts[i + 1] - (esizes[i] + estarts[i]) for i in range(len(esizes) - 1)]
        isotoinfo[iso] = {'chrom': chrom, 'strand': strand, 'start': start, 'end': end, 'esizes': esizes, 'introns': isizes}
    return isotoinfo

def convert_start_pos(isoinfo, isostart):
    esizes = isoinfo['esizes']
    introns = isoinfo['introns']
    if isoinfo['strand'] == '-':
        esizes, introns = esizes[::-1], introns[::-1]
    for i in range(len(esizes)):
        if isostart < sum(esizes[:i + 1]):
            if isoinfo['strand'] == '+':
                return isoinfo['start'] + isostart + sum(introns[:i])
            else:
                return isoinfo['end'] - (isostart + sum(introns[:i]))  # this is now the end position of the transcript on the genome


def add_introns_to_block(exonbounds, introns, tpos, blocktype, blocklen, newcigar):
    for e in range(len(exonbounds)):
        if tpos <= exonbounds[e] < tpos + blocklen:
            if exonbounds[e] - tpos > 0:
                newcigar.append((blocktype, exonbounds[e] - tpos))
            newcigar.append((3, introns[e]))  # insert intron
            blocklen -= exonbounds[e] - tpos
            tpos += exonbounds[e] - tpos
    tpos += blocklen
    return exonbounds, tpos, blocktype, blocklen, newcigar

def convert_cigar(isoinfo, cigar, startpos):
    tpos = startpos
    newcigar = []
    esizes = isoinfo['esizes']
    introns = isoinfo['introns']
    if isoinfo['strand'] == '-':
        esizes, introns = esizes[::-1], introns[::-1]
    exonbounds = [sum(esizes[:i + 1]) for i in range(len(esizes) - 1)]  # excludes last exon, that boundary is not intronic
    for blocktype, blocklen in cigar:
        if blocktype in {0, 2, 3, 7, 8}:  # consumes reference
            # check if this crosses an exon boundary
            exonbounds, tpos, blocktype, blocklen, newcigar = add_introns_to_block(exonbounds, introns, tpos, blocktype, blocklen, newcigar)

        if blocklen > 0:  # insertions always written, but don't alter tpos. Also write any remaining part of block
            newcigar.append((blocktype, blocklen))
    if isoinfo['strand'] == '-':
        newcigar = newcigar[::-1]  # on the negative strand, must reverse cigar string
    return newcigar

def generate_new_segment(s=None, newname=None, newseq=None, newquals=None, newflag=None, newtags=None, newrefid=None, newrefstart=None, newcigart=None, newmapq=None):
    a = pysam.AlignedSegment()
    a.query_name = newname if newname is not None else s.query_name
    a.query_sequence = newseq if newseq is not None else s.query_sequence
    a.flag = newflag if newflag is not None else s.flag
    a.mapping_quality = newmapq if newmapq is not None else s.mapping_quality
    a.query_qualities = newquals if newquals is not None else s.query_qualities
    a.tags = newtags if newtags is not None else s.tags
    a.reference_id = newrefid if newrefid is not None else s.reference_id
    a.reference_start = newrefstart if newrefstart is not None else s.reference_start
    a.cigartuples = newcigart if newcigart is not None else s.cigartuples
    return a


isoformbed = sys.argv[1]
transcriptomebam = sys.argv[2]
genomesamplebam = sys.argv[3]
outprefix = sys.argv[4]


isotoinfo = get_iso_info_from_bed(isoformbed)


samfile = pysam.AlignmentFile(transcriptomebam, 'rb')
headerfile = pysam.AlignmentFile(genomesamplebam, 'rb')
outfile = pysam.AlignmentFile(f'{outprefix}.genomelift.bam', 'wb', template=headerfile)

outchroms = outfile.references

# CURRENTLY NONE OF THIS WORKS FOR NEGATIVE STRAND READS
for s in samfile:
    if s.is_mapped:
        thisisoinfo = isotoinfo[s.reference_name]
        newrefid = outchroms.index(thisisoinfo['chrom'])
        newrefstart = convert_start_pos(thisisoinfo, s.reference_start)
        newcigart = convert_cigar(thisisoinfo, s.cigartuples, s.reference_start)
        if thisisoinfo['strand'] == '+':
            a = generate_new_segment(s, newrefid=newrefid, newrefstart=newrefstart, newcigart=newcigart)
        else:
            newcigart = convert_cigar(thisisoinfo, s.cigartuples, s.reference_start)
            totalaligneddist = sum([x[1] for x in newcigart if x[0] not in {1, 4, 5}])
            newrefstart = newrefstart - totalaligneddist  # was read end, correct to read start
            newflag = s.flag - 16 if s.is_reverse else s.flag + 16
            newseq = get_reverse_complement(s.query_sequence) if s.query_sequence else None
            newquals = s.query_qualities[::-1] if s.query_qualities else None
            a = generate_new_segment(s, newrefid=newrefid, newrefstart=newrefstart, newcigart=newcigart,
                                     newflag=newflag, newquals=newquals, newseq=newseq)

        outfile.write(a)


headerfile.close()
samfile.close()
outfile.close()

# pysam.sort('-o', f'{sample}.flairomealigned.mygenomelift.sorted.bam', f'{sample}.flairomealigned.mygenomelift.bam')
# pysam.index(f'{sample}.flairomealigned.mygenomelift.sorted.bam')
pysam.sort('-o', f'{outprefix}.genomelift.sorted.bam', f'{outprefix}.genomelift.bam')
pysam.index(f'{outprefix}.genomelift.sorted.bam')
