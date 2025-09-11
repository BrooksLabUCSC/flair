#!/usr/bin/env python3

import os
import sys
import argparse
import pysam
import logging
import shutil
from flair.count_sam_transcripts import *
from flair.flair_transcriptome import makecorrecttempdir
import multiprocessing as mp
from time import sleep
from flair import FlairInputDataError


def generate_alignment_obj_for_read(args, genome, transcript_to_exons, transcriptaligns, header):
    filteredtranscriptaligns = {}
    for alignment in transcriptaligns:
        transcript = alignment.reference_name
        if args.remove_internal_priming:
            intprimannot = transcript_to_exons if args.permissive_last_exons else None
            notinternalpriming = remove_internal_priming.removeinternalpriming(alignment.reference_name,
                                                                               alignment.reference_start,
                                                                               alignment.reference_end, False,
                                                                               genome, None, intprimannot,
                                                                               args.intprimingthreshold,
                                                                               args.intprimingfracAs)
        else:
            notinternalpriming = True
        if transcript == 'HISEQ:1287:HKCG7BCX3:1:1102:8019:82885': print(notinternalpriming)
        if notinternalpriming:
            pos = alignment.reference_start
            try:
                alignscore = alignment.get_tag('AS')
                mdtag = alignment.get_tag('MD')
            except KeyError as ex:
                raise Exception(f"Missing AS or MD tag in alignment of '{alignment.query_name}' in '{args.sam}'") from ex
            cigar = alignment.cigartuples
            tlen = header.get_reference_length(transcript)
            filteredtranscriptaligns[transcript] = IsoAln(transcript, pos, cigar, tlen, alignscore, mdtag)
    return filteredtranscriptaligns



def process_read_chunk(chunkinfo):
    chunkindex, readstoaligns, tempDir, transcript_to_exons, transcript_to_bp_ss_index, args, headeroutfilename, clippingdata, transcript_to_genomic_ends = chunkinfo
    genome = None
    if args.remove_internal_priming:
        genome = pysam.FastaFile(args.transcriptomefasta)

    headerfile = pysam.AlignmentFile(headeroutfilename, 'rb')
    if args.output_bam:
        tempOutFile = pysam.AlignmentFile(tempDir + 'readChunk' + str(chunkindex) + '.bam', 'wb', template=headerfile)
    results = []

    for readname in readstoaligns:
        transcriptaligns = [pysam.AlignedSegment.fromstring(x, headerfile.header) for x in readstoaligns[readname]]
        filteredtranscriptaligns = generate_alignment_obj_for_read(args, genome, transcript_to_exons, transcriptaligns, headerfile.header)
        finaltnames = []
        thisclipping = clippingdata[readname] if readname in clippingdata else None
        if len(filteredtranscriptaligns) > 0:
            # print(readname)
            assignedts = get_best_transcript(filteredtranscriptaligns, args, transcript_to_exons, transcript_to_bp_ss_index, thisclipping, transcript_to_genomic_ends)
            if assignedts:
                for assignedt, gtstart, gtend in assignedts:
                    finaltnames.append(assignedt)
                    results.append((readname, assignedt, gtstart, gtend))
        if args.output_bam and len(finaltnames) > 0:
            readseq, readerr = None, None
            for alignment in transcriptaligns:
                if not alignment.is_secondary: #supplementary already filtered earlier
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
                    tempOutFile.write(alignment)
    if genome:
        genome.close()
    if args.output_bam:
        tempOutFile.close()
    return results


def report_thread_error(error):
    raise ValueError(error)


def bam_to_read_aligns(samfile, chunksize, tempDir, transcript_to_exons, transcript_to_bp_ss_index,
                                            args, headeroutfilename, readstoclipping, transcript_to_genomic_ends):
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
                yield (chunkindex, readchunk, tempDir, transcript_to_exons, transcript_to_bp_ss_index, args, headeroutfilename, clippingdata, transcript_to_genomic_ends)
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
    if len(lastaligns) > 0:
        readchunk[lastname] = lastaligns
        if lastname in readstoclipping:
            clippingdata[lastname] = readstoclipping[lastname]
    if len(readchunk) > 0:
        logging.info(f'\rstarting chunk {chunkindex}')
        yield (chunkindex, readchunk, tempDir, transcript_to_exons, transcript_to_bp_ss_index, args, headeroutfilename, clippingdata, transcript_to_genomic_ends)

def process_alignments(args, transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends):
    logging.info('processing alignments')
    samfile = pysam.AlignmentFile(args.sam, 'r')
    tempDir = makecorrecttempdir()
    headeroutfilename = tempDir + 'headerfile.bam'
    hfile = pysam.AlignmentFile(headeroutfilename, 'wb', template=samfile)
    hfile.close()
    pysam.index(tempDir + 'headerfile.bam')

    readstoclipping = {}
    if args.trimmedreads:
        for line in open(args.trimmedreads):
            rname, clipping = line.rstrip().split('\t')
            readstoclipping[rname] = int(clipping)

    chunksize = 1000

    chunkresults = []
    mp.set_start_method('fork')
    p = mp.Pool(args.threads)

    args.sam = '' # required to pass args to multiprocessing

    # write method to yield chunks
    # for chunk in chunkyielder

    for r in p.imap_unordered(process_read_chunk, bam_to_read_aligns(samfile, chunksize, tempDir, transcript_to_exons,
                                                                transcript_to_bp_ss_index, args, headeroutfilename, readstoclipping, transcript_to_genomic_ends)):
        chunkresults.append(r)

    p.close()
    p.join()
    logging.info('starting to combine temp files')

    if args.output_endpos:
        endout = open(args.output_endpos, 'w')

    transcripttoreads = {}
    for i in range(len(chunkresults)):
        for read, transcript, gtstart, gtend in chunkresults[i]:
            if transcript not in transcripttoreads: transcripttoreads[transcript] = []
            transcripttoreads[transcript].append(read)
            if args.output_endpos:
                endout.write('\t'.join([str(x) for x in [read, transcript, gtstart, gtend]]) + '\n')
    if args.output_endpos:
        endout.close()
    write_output(args, transcripttoreads)

    if args.output_bam:
        outfile = pysam.AlignmentFile(tempDir + 'combined_unsorted.bam', 'wb', template=samfile)
        for i in range(len(chunkresults)):
            tempfile = pysam.AlignmentFile(tempDir + 'readChunk' + str(i+1) + '.bam', 'rb')
            for read in tempfile:
                outfile.write(read)
            tempfile.close()
        outfile.close()
        pysam.sort('-o', args.output_bam, tempDir + 'combined_unsorted.bam')
        pysam.index(args.output_bam)
    shutil.rmtree(tempDir)

if __name__ == '__main__':
    logging.info('processing annotation')
    args = parse_args()
    args = check_args(args)
    transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends = get_annot_info(args)
    process_alignments(args, transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends)


