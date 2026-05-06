#!/usr/bin/env python3

# FIXME: this is convoluted, this reuses count_sam_transcripts command line an functions
# these should be in a libary

import pysam
import logging
import shutil
from flair.count_sam_transcripts import parse_args, check_args, read_isoforms_bed, IsoAln, get_best_transcript, write_output
from flair.remove_internal_priming import removeinternalpriming
import multiprocessing as mp
from flair import FlairInputDataError
from flair.io_utils import make_temp_dir


def generate_alignment_obj_for_read(args, genome, transcript_to_exons, transcriptaligns, header):
    filteredtranscriptaligns = {}
    for alignment in transcriptaligns:
        transcript = alignment.reference_name
        if args.remove_internal_priming:
            intprim_annot = transcript_to_exons if args.permissive_last_exons else None
            not_internal_priming = removeinternalpriming(alignment.reference_name,
                                                         alignment.reference_start,
                                                         alignment.reference_end, False,
                                                         genome, None, intprim_annot,
                                                         args.intprimingthreshold,
                                                         args.intprimingfracAs)
        else:
            not_internal_priming = True
        if not not_internal_priming:
            logging.debug(f"read dropped: internal priming on {transcript}: {alignment.query_name}")
        else:
            pos = alignment.reference_start
            try:
                alignscore = alignment.get_tag('AS')
                mdtag = alignment.get_tag('MD')
            except KeyError as ex:
                raise FlairInputDataError(f"Missing AS or MD tag in alignment of '{alignment.query_name}' in '{args.sam}'") from ex
            cigar = alignment.cigartuples
            tlen = header.get_reference_length(transcript)
            filteredtranscriptaligns[transcript] = IsoAln(transcript, pos, cigar, tlen, alignscore, mdtag)
    return filteredtranscriptaligns

def process_read_chunk(chunkinfo):  # noqa: C901 - FIXME: reduce complexity
    chunkindex, readstoaligns, temp_dir, info, args, headeroutfilename, clippingdata = chunkinfo
    genome = None
    if args.remove_internal_priming:
        genome = pysam.FastaFile(args.transcriptomefasta)

    headerfile = pysam.AlignmentFile(headeroutfilename, 'rb')
    if args.output_bam:
        temp_out_file = pysam.AlignmentFile(temp_dir + 'readChunk' + str(chunkindex) + '.bam', 'wb', template=headerfile)
    results = []

    for readname in readstoaligns:
        transcriptaligns = [pysam.AlignedSegment.fromstring(x, headerfile.header) for x in readstoaligns[readname]]
        filteredtranscriptaligns = generate_alignment_obj_for_read(args, genome, info.transcript_to_exons, transcriptaligns, headerfile.header)
        finaltnames = []
        this_clipping = clippingdata[readname] if readname in clippingdata else None
        if len(filteredtranscriptaligns) > 0:
            assignedts = get_best_transcript(filteredtranscriptaligns, info, this_clipping,
                                             stringent=args.stringent, check_splice=args.check_splice,
                                             fusion_breakpoints=args.fusion_breakpoints,
                                             allow_UTR_indels=args.allow_UTR_indels,
                                             trimmedreads=args.trimmedreads,
                                             soft_clipping_buffer=args.soft_clipping_buffer,
                                             output_endpos=args.output_endpos,
                                             trust_ends=args.trust_ends, end_norm_dist=args.end_norm_dist)
            if not assignedts:
                logging.debug(f"read dropped: no passing transcript assignment: {readname}")
            else:
                for assignedt, gtstart, gtend in assignedts:
                    finaltnames.append(assignedt)
                    results.append((readname, assignedt, gtstart, gtend))
        else:
            logging.debug(f"read dropped: all transcript alignments filtered (internal priming or other): {readname}")
        if args.output_bam and len(finaltnames) > 0:
            readseq, readerr = None, None
            for alignment in transcriptaligns:
                if not alignment.is_secondary:  # supplementary already filtered earlier
                    readseq = alignment.query_sequence
                    readerr = alignment.query_qualities

            for alignment in transcriptaligns:
                if alignment.reference_name in finaltnames:
                    alignment.mapping_quality = 60
                    if alignment.is_secondary:
                        alignment.query_sequence = readseq
                        alignment.query_qualities = readerr
                        if alignment.is_reverse:
                            alignment.flag = 16
                        else:
                            alignment.flag = 0
                        alignment.cigarstring = alignment.cigarstring.replace('H', 'S')
                    temp_out_file.write(alignment)
    if genome:
        genome.close()
    if args.output_bam:
        temp_out_file.close()
    return results


def bam_to_read_aligns(samfile, chunksize, temp_dir, info,  # noqa: C901 - FIXME: reduce complexity
                       args, headeroutfilename, readstoclipping):
    lastname = None
    lastaligns = []
    readchunk = {}
    chunkindex = 1
    clippingdata = {}
    for read in samfile:
        readname = read.query_name
        if readname != lastname:
            if len(readchunk) == chunksize:
                logging.info(f'\rstarting chunk {chunkindex}')
                yield (chunkindex, readchunk, temp_dir, info, args, headeroutfilename, clippingdata)
                readchunk = {}
                clippingdata = {}
                chunkindex += 1
            if len(lastaligns) > 0:
                readchunk[lastname] = lastaligns
                if lastname in readstoclipping:
                    clippingdata[lastname] = readstoclipping[lastname]
                lastaligns = []
            lastname = readname
        # removing supplementary alignments has the biggest effect on the output. I think we should do it though
        if read.is_mapped and not read.is_supplementary and read.mapping_quality >= args.quality:
            lastaligns.append(read.to_string())
        else:
            if not read.is_mapped:
                logging.debug(f"read dropped: unmapped: {read.query_name}")
            elif read.is_supplementary:
                logging.debug(f"read dropped: supplementary alignment: {read.query_name}")
            else:
                logging.debug(f"read dropped: low quality ({read.mapping_quality} < {args.quality}): {read.query_name}")
    if len(lastaligns) > 0:
        readchunk[lastname] = lastaligns
        if lastname in readstoclipping:
            clippingdata[lastname] = readstoclipping[lastname]
    if len(readchunk) > 0:
        logging.info(f'\rstarting chunk {chunkindex}')
        yield (chunkindex, readchunk, temp_dir, info, args, headeroutfilename, clippingdata)

def process_alignments(args, info):  # noqa: C901 - FIXME: reduce complexity
    logging.info('processing alignments')
    samfile = pysam.AlignmentFile(args.sam, 'r')
    prefix = args.output.split('.txt')[0]
    temp_dir = make_temp_dir(prefix)
    headeroutfilename = temp_dir + 'headerfile.bam'
    hfile = pysam.AlignmentFile(headeroutfilename, 'wb', template=samfile)
    hfile.close()
    pysam.index(temp_dir + 'headerfile.bam')

    readstoclipping = {}
    if args.trimmedreads:
        for line in open(args.trimmedreads):
            rname, clipping = line.rstrip().split('\t')
            readstoclipping[rname] = int(clipping)

    chunksize = 1000

    chunkresults = []
    mp.set_start_method('fork', force=True)

    args.sam = ''   # required to pass args to multiprocessing

    # write method to yield chunks
    # for chunk in chunkyielder

    with mp.Pool(args.threads) as p:
        for r in p.imap_unordered(process_read_chunk, bam_to_read_aligns(samfile, chunksize, temp_dir, info,
                                                                         args, headeroutfilename, readstoclipping)):
            chunkresults.append(r)
    logging.info('starting to combine temp files')

    if args.output_endpos:
        endout = open(args.output_endpos, 'w')

    transcript_to_reads = {}
    for i in range(len(chunkresults)):
        for read, transcript, gtstart, gtend in chunkresults[i]:
            if transcript not in transcript_to_reads:
                transcript_to_reads[transcript] = []
            transcript_to_reads[transcript].append((read, gtstart, gtend))
            if args.output_endpos:
                endout.write('\t'.join([str(x) for x in [read, transcript, gtstart[0], gtstart[1], gtend[0], gtend[1]]]) + '\n')
    if args.output_endpos:
        endout.close()
    write_output(args, transcript_to_reads)

    if args.output_bam:
        outfile = pysam.AlignmentFile(temp_dir + 'combined_unsorted.bam', 'wb', template=samfile)
        for i in range(len(chunkresults)):
            tempfile = pysam.AlignmentFile(temp_dir + 'readChunk' + str(i + 1) + '.bam', 'rb')
            for read in tempfile:
                outfile.write(read)
            tempfile.close()
        outfile.close()
        pysam.sort('-o', args.output_bam, temp_dir + 'combined_unsorted.bam')
        pysam.index(args.output_bam)
    shutil.rmtree(temp_dir)


if __name__ == '__main__':
    logging.info('processing annotation')

    args = parse_args()
    args = check_args(args)
    info = read_isoforms_bed(
        isoforms=args.isoforms, stringent=args.stringent, check_splice=args.check_splice,
        fusion_dist=args.fusion_dist, fusion_breakpoints=args.fusion_breakpoints,
        output_endpos=args.output_endpos, unique_bound=args.unique_bound)
    process_alignments(args, info)
