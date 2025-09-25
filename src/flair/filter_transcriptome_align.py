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


def generate_alignment_obj_for_read(args, genome, transcripttoexons, transcriptaligns, header):
    filteredtranscriptaligns = {}
    for alignment in transcriptaligns:
        transcript = alignment.reference_name
        if args.remove_internal_priming:
            intprimannot = transcripttoexons if args.permissive_last_exons else None
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
    chunkindex, readstoaligns, tempDir, transcripttoexons, transcripttobpssindex, args, headeroutfilename = chunkinfo
    genome = None
    if args.remove_internal_priming:
        genome = pysam.FastaFile(args.transcriptomefasta)

    headerfile = pysam.AlignmentFile(headeroutfilename, 'rb')
    if args.output_bam:
        tempOutFile = pysam.AlignmentFile(tempDir + 'readChunk' + str(chunkindex) + '.bam', 'wb', template=headerfile)

    results = []

    for readname in readstoaligns:
        transcriptaligns = [pysam.AlignedSegment.fromstring(x, headerfile.header) for x in readstoaligns[readname]]
        filteredtranscriptaligns = generate_alignment_obj_for_read(args, genome, transcripttoexons, transcriptaligns, headerfile.header)
        finaltnames = []
        if len(filteredtranscriptaligns) > 0:
            assignedts = getbesttranscript(filteredtranscriptaligns, args, transcripttoexons, transcripttobpssindex)
            if assignedts:
                for assignedt in assignedts:
                    finaltnames.append(assignedt)
                    results.append((readname, assignedt))
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


def bam_to_read_aligns(samfile, chunksize, tempDir, transcripttoexons, transcripttobpssindex,
                                            args, headeroutfilename):
    lastname = None
    lastaligns = []
    readchunk = {}
    chunkindex = 1
    for read in samfile:
        readname = read.query_name
        if readname != lastname:
            if len(readchunk) == chunksize:
                logging.info(f'\rstarting chunk {chunkindex}')
                yield (chunkindex, readchunk, tempDir, transcripttoexons, transcripttobpssindex, args, headeroutfilename)
                readchunk = {}
                chunkindex += 1
            if len(lastaligns) > 0:
                readchunk[lastname] = lastaligns
                lastaligns = []
            lastname = readname
        # removing supplementary alignments has the biggest effect on the output. I think we should do it though
        if read.is_mapped and not read.is_supplementary and read.mapping_quality >= args.quality:
            lastaligns.append(read.to_string())
    if len(lastaligns) > 0:
        readchunk[lastname] = lastaligns
    if len(readchunk) > 0:
        logging.info(f'\rstarting chunk {chunkindex}')
        yield (chunkindex, readchunk, tempDir, transcripttoexons, transcripttobpssindex, args, headeroutfilename)

def process_alignments(args, transcripttoexons, transcripttobpssindex):
    logging.info('processing alignments')
    samfile = pysam.AlignmentFile(args.sam, 'r')
    tempDir = makecorrecttempdir()
    headeroutfilename = tempDir + 'headerfile.bam'
    hfile = pysam.AlignmentFile(headeroutfilename, 'wb', template=samfile)
    hfile.close()
    pysam.index(tempDir + 'headerfile.bam')

    chunksize = 1000

    chunkresults = []
    mp.set_start_method('fork')
    p = mp.Pool(args.threads)

    args.sam = '' # required to pass args to multiprocessing

    # write method to yield chunks
    # for chunk in chunkyielder

    for r in p.imap_unordered(process_read_chunk, bam_to_read_aligns(samfile, chunksize, tempDir, transcripttoexons,
                                                                transcripttobpssindex, args, headeroutfilename)):
        chunkresults.append(r)

    p.close()
    p.join()
    logging.info('starting to combine temp files')

    transcripttoreads = {}
    for i in range(len(chunkresults)):
        for read, transcript in chunkresults[i]:
            if transcript not in transcripttoreads:
                transcripttoreads[transcript] = []
            transcripttoreads[transcript].append(read)
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
    args = parseargs()
    args = checkargs(args)
    transcripttoexons, transcripttobpssindex = getannotinfo(args)
    process_alignments(args, transcripttoexons, transcripttobpssindex)
