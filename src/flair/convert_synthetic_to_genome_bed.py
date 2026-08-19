import logging
from flair.pycbio.hgdata.bed import BedBlock, BedReader
from flair.flair_bed import FlairBed

def get_iso_to_reads(readmapfile):
    isoreadsup = {}
    for line in open(readmapfile):
        line = line.rstrip().split('\t')
        reads = line[1].split(',')
        isoreadsup[line[0]] = set(reads)
    return isoreadsup

def get_synth_info(breakpointfile):
    synthchrtoinfo = {}
    for bed in BedReader(breakpointfile, numStdCols=4):
        synthchrtoinfo[bed.chrom] = '--'.join(bed.name.split('--')[1:])
    return synthchrtoinfo

def get_paralog_ref(paralogfile):
    grouptogenes = {}
    genetoparalogs = {}
    # 5  244     2       69321074        69338940        1       ENSG00000205572 SERF1B  small EDRK-rich factor 1B (centromeric) [Source:HGNC Symbol;Acc:10756]
    for line in open(paralogfile):
        line = line.split('\t')
        group = line[1]
        gname = line[6]
        if group not in grouptogenes:
            grouptogenes[group] = set()
        grouptogenes[group].add(gname)

    for group in grouptogenes:
        gset = grouptogenes[group]
        for g in gset:
            genetoparalogs[g] = group
    return genetoparalogs

def get_gene_name_conv(annotgtf):
    from flair.gtf_io import gtf_record_parser, GtfAttrsSet
    genetoname = {}
    for rec in gtf_record_parser(annotgtf, include_features={'gene'}, attrs=GtfAttrsSet.ALL):
        geneid = rec.gene_id
        genename = rec.gene_name if rec.gene_name else geneid
        genetoname[geneid.split('.')[0]] = genename
    return genetoname

def identify_promiscuous_genes(isoformsbed, genetoparalogs):
    from flair.pycbio.hgdata.bed import BedReader
    locustopartners = {}
    for bed in BedReader(isoformsbed, fixScores=True):
        fnames = set(bed.chrom.split('--'))
        fnames = {x.split('.')[0] if x[:3] != 'chr' else x.split('-')[0] + '-' + str(round(int(x.split('-')[1]), -6))
                  for x in fnames}
        for i in fnames:
            other = fnames - {i, }
            newother = frozenset([genetoparalogs[g] if g in genetoparalogs else g for g in other])
            if i not in locustopartners:
                locustopartners[i] = set()
            locustopartners[i].add(newother)
    return locustopartners


def identify_fusion_problems(fgenes, locustopartners, maxpromiscuity, genetoname, genetoparalogs, genomic_chroms, isosup):

    ispromiscuous, areparalogs, areig, allnames, gcount = False, [], [], set(), []
    for g in fgenes:
        gcount.append(list(fgenes).count(g))
        if len(locustopartners[g]) > maxpromiscuity:
            ispromiscuous = True
        if g in genetoname:
            allnames.add(genetoname[g][:4])
            if genetoname[g][:3] in {'IGK', 'IGH', 'IGL', 'IGF'}:
                areig.append(True)
            else:
                areig.append(False)
        else:
            allnames.add(g)
            # if g[:3] == 'chr':
            #     areig.append(True)
            # else:
            areig.append(False)
        if g in genetoparalogs:
            other = set(fgenes) - {g, }
            other = {genetoparalogs[g2] if g2 in genetoparalogs else g2 for g2 in other}
            if genetoparalogs[g] in other or len(other) == 0:
                areparalogs.append(True)
            else:
                areparalogs.append(False)
        else:
            areparalogs.append(False)
    overall = (all([x == 1 for x in gcount]) and isosup >= 1 and not ispromiscuous
               and any(x is False for x in areparalogs) and all(x is False for x in areig)
               and len(allnames) > 1 and 'chrM' not in genomic_chroms
               and all([x[:3] == 'chr' for x in genomic_chroms]))
    # if overall:
    #     print(fgenes, 'prom', ispromiscuous, 'para', areparalogs, 'ig', areig, 'uniquenames', allnames, 'goodchroms', 'chrM' not in genomic_chroms and all([x[:3] == 'chr' for x in genomic_chroms]), 'gcount', gcount)
    return overall

def get_locus_bounds(synthinfo):
    locuslen = [abs(x[3] - x[2]) for x in synthinfo]
    locusbounds = []
    laststart = 0
    for i in range(len(locuslen)):
        locusbounds.append((laststart, locuslen[i] + laststart))
        laststart += locuslen[i]
    return locusbounds

def separate_exons_by_locus(esizes, estarts, numloci, locusbounds, start, thickStart, thickEnd):
    starts = [None for x in range(numloci)]
    exonindexes = [[] for x in range(numloci)]
    thickEdges = [[None, None] for x in range(numloci)]
    for order in range(numloci):
        if locusbounds[order][0] < thickStart < locusbounds[order][1]:
            thickEdges[order][0] = thickStart
        if locusbounds[order][0] < thickEnd < locusbounds[order][1]:
            thickEdges[order][1] = thickEnd
    for i in range(len(esizes)):
        thisstart, thisend = start + estarts[i], start + estarts[i] + esizes[i]
        for order in range(numloci):
            if locusbounds[order][0] <= thisstart and thisend <= locusbounds[order][1]:
                if starts[order] is None:
                    starts[order] = estarts[i]  # thisstart #- locusbounds[order][0]
                exonindexes[order].append(i)
    return starts, exonindexes, thickEdges

def convert_to_genomic_coords(numloci, synthinfo, exonindexes, starts, locusbounds, esizes, estarts, start, iso, hasThickEdges, thickEdges):
    genomicbounds = []
    bedinfo = []
    for order in range(numloci):
        genename, genomicchr, leftbound, rightbound = synthinfo[order]
        locusesizes = [esizes[i] for i in exonindexes[order]]
        locusestarts = [estarts[i] for i in exonindexes[order]]
        locusestarts = [x - starts[order] for x in locusestarts]
        totlen = locusestarts[-1] + locusesizes[-1]
        locusdir = '+' if leftbound < rightbound else '-'
        # name = 'gene' + str(order + 1) + '_' + iso
        thickStart, thickEnd = None, None
        if leftbound > rightbound:  # reverse direction
            temp = []
            for i in range(len(locusestarts) - 1, -1, -1):
                temp.append(totlen - (locusestarts[i] + locusesizes[i]))
            locusestarts = temp
            locusesizes = locusesizes[::-1]
            chromStart = leftbound - (((start + starts[order]) - locusbounds[order][0]) + totlen)
            if thickEdges[order][0] is not None:
                thickEnd = leftbound - (thickEdges[order][0] - locusbounds[order][0])
            if thickEdges[order][1] is not None:
                thickStart = leftbound - (thickEdges[order][1] - locusbounds[order][0])
        else:
            chromStart = leftbound + ((start + starts[order]) - locusbounds[order][0])
            if thickEdges[order][0] is not None:
                thickStart = leftbound + (thickEdges[order][0] - locusbounds[order][0])
            if thickEdges[order][1] is not None:
                thickEnd = leftbound + (thickEdges[order][1] - locusbounds[order][0])
        chromEnd = chromStart + totlen
        if hasThickEdges:
            thickStart = chromStart if thickStart is None else thickStart
            thickEnd = chromEnd if thickEnd is None else thickEnd

        genomicbounds.append((genomicchr, chromStart, chromEnd, locusdir))
        blocks = [BedBlock(chromStart + locusestarts[i], chromStart + locusestarts[i] + locusesizes[i])
                  for i in range(len(locusesizes))]
        bedinfo.append((genomicchr, chromStart, chromEnd, locusdir, thickStart, thickEnd, blocks))
    return genomicbounds, bedinfo

def check_too_close(numloci, genomicbounds, min_dist):
    for i in range(numloci - 1):
        if genomicbounds[i][0] == genomicbounds[i + 1][0] and genomicbounds[i][3] == genomicbounds[i + 1][3]:  # chrom and strand are same
            if genomicbounds[i][3] == '+' and 0 < genomicbounds[i + 1][1] - genomicbounds[i][2] < min_dist:
                return True
            # when in negative direction, bounds are flipped
            elif genomicbounds[i][3] == '-' and 0 < genomicbounds[i][2] - genomicbounds[i + 1][1] < min_dist:
                return True
    return False

def write_final_fusion_reads(readsfile, freadsfinal):
    last = False
    temp = readsfile.split('.prelim')
    freads = open(''.join(temp), 'w')

    if readsfile.split('.')[-1] == 'fasta' or readsfile.split('.')[-1] == 'fa':
        for line in open(readsfile):
            if line[0] == '>':
                readname = line[1:].rstrip().split()[0]
                if readname in freadsfinal:
                    last = True
                else:
                    last = False
            if last:
                freads.write(line)
    else:
        linecount = 0
        for line in open(readsfile):
            if linecount % 4 == 0:
                readname = line[1:].rstrip().split()[0]
                if readname in freadsfinal:
                    last = True
                else:
                    last = False
            if last:
                freads.write(line)
            linecount += 1
    freads.close()

def get_gene_id(gene_og_to_new, gene_og, curr_gene_count):
    if gene_og not in gene_og_to_new:
        gene_id = f'FLG{curr_gene_count:08d}'
        gene_info = tuple(gene_og.split('--'))
        gene_og_to_new[gene_og] = (gene_id, gene_info)
        curr_gene_count += 1
    gene_id, gene_info = gene_og_to_new[gene_og]
    return curr_gene_count, gene_id, gene_info

def convert_synthetic_isos(isoformsbed, readmapfile, readsfile, breakpointfile, outname, min_dist):
    isoreadsup = get_iso_to_reads(readmapfile)
    synthchrtoinfo = get_synth_info(breakpointfile)

    freadsfinal = set()
    out = open(outname, 'w')
    gene_og_to_new = {}
    curr_gene_count = 1
    for bed in BedReader(isoformsbed, fixScores=True, bedClass=FlairBed):
        iso, start = bed.name, bed.chromStart
        esizes = [len(blk) for blk in bed.blocks]
        estarts = [blk.start - start for blk in bed.blocks]
        synthinfo = [x.split('..') for x in synthchrtoinfo[bed.chrom].split('--')]
        synthinfo = [[y[0], y[1], int(y[2]), int(y[3])] for y in synthinfo]

        locusbounds = get_locus_bounds(synthinfo)
        if start < locusbounds[0][1] and locusbounds[-1][0] < bed.chromEnd:
            numloci = len(synthinfo)

            # FIXME add handling of existing thickStart and thickEnd values from synthetic aligned bed,
            hasThickEdges = bed.thickStart != bed.thickEnd
            starts, exonindexes, thickEdges = separate_exons_by_locus(esizes, estarts, numloci, locusbounds, start, bed.thickStart, bed.thickEnd)

            assigned_exons = sorted(i for indexes in exonindexes for i in indexes)
            if assigned_exons != list(range(len(esizes))):
                logging.warning('skipping %s: not all synthetic exons could be assigned to exactly one locus', iso)
                continue

            if None not in starts:
                genomicbounds, bedinfos = convert_to_genomic_coords(numloci, synthinfo, exonindexes, starts, locusbounds, esizes, estarts, start, iso, hasThickEdges, thickEdges)
                if not check_too_close(numloci, genomicbounds, min_dist):
                    for fusion_order, (genomicchr, chromStart, chromEnd, locusdir, thickStart, thickEnd, blocks) in enumerate(bedinfos):
                        curr_gene_count, fusion_id, ref_gene_info = get_gene_id(gene_og_to_new, bed.chrom, curr_gene_count)
                        flair_genes = []
                        for gene in ref_gene_info:
                            curr_gene_count, gene_id, ref_gene_ids = get_gene_id(gene_og_to_new, gene, curr_gene_count)
                            flair_genes.append(gene_id)
                        # convert reference genes that are chromosome regions to none in the ref_gene_mappings field
                        ref_gene_info = tuple([x if x[:3] != 'chr' else None for x in ref_gene_info])
                        outbed = FlairBed(genomicchr, chromStart, chromEnd, name=bed.name, score=bed.score, strand=locusdir, blocks=blocks,
                                          gene_id=fusion_id, ref_gene_mappings=ref_gene_info, thickStart=thickStart, thickEnd=thickEnd,
                                          read_support=bed.read_support, frac_support=bed.frac_support, productivity=bed.productivity, transcript_class='fusion',
                                          fused_genes=tuple(flair_genes), pos_in_fusion=fusion_order, samples=bed.samples)
                        outbed.write(out)
                    freadsfinal.update(isoreadsup[iso])
    out.close()
    write_final_fusion_reads(readsfile, freadsfinal)
