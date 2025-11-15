
import sys, os


def get_iso_to_reads(readmapfile):
    isoreadsup = {}
    for line in open(readmapfile):
        line = line.split('\t')
        reads = line[1].split(',')
        isoreadsup[line[0]] = set(reads)
    return isoreadsup

def get_synth_info(breakpointfile):
    synthchrtoinfo = {}
    for line in open(breakpointfile):
        line = line.rstrip().split('\t')
        synthchrtoinfo[line[0]] = '--'.join(line[-1].split('--')[1:])
    return synthchrtoinfo

def get_paralog_ref(paralogfile):
    grouptogenes = {}
    genetoparalogs = {}
    # 5  244     2       69321074        69338940        1       ENSG00000205572 SERF1B  small EDRK-rich factor 1B (centromeric) [Source:HGNC Symbol;Acc:10756]
    for line in open(paralogfile):
        line = line.split('\t')
        group = line[1]
        gname = line[6]
        if group not in grouptogenes: grouptogenes[group] = set()
        grouptogenes[group].add(gname)

    for group in grouptogenes:
        gset = grouptogenes[group]
        for g in gset:
            genetoparalogs[g] = group
    return genetoparalogs

def get_gene_name_conv(annotgtf):
    genetoname = {}
    for line in open(annotgtf):
        if line[0] != '#':
            line = line.split('\t', 3)
            if line[2] == 'gene':
                geneid = line[-1].split('gene_id "')[1].split('"')[0]
                if 'gene_name' in line[-1]:
                    genename = line[-1].split('gene_name "')[1].split('"')[0]
                else:
                    genename = geneid
                genetoname[geneid.split('.')[0]] = genename
    return genetoname

def identify_promiscuous_genes(isoformsbed, genetoparalogs):
    locustopartners = {}
    for line in open(isoformsbed):
        fnames = set(line.split('\t', 1)[0].split('--'))
        fnames = {x.split('.')[0] if x[:3] != 'chr' else x.split('-')[0] + '-' + str(round(int(x.split('-')[1]), -6))
                  for x in fnames}
        for i in fnames:
            other = fnames - {i, }
            newother = frozenset([genetoparalogs[g] if g in genetoparalogs else g for g in other])
            if i not in locustopartners: locustopartners[i] = set()
            locustopartners[i].add(newother)
    return locustopartners


def identify_fusion_problems(fusionchr, locustopartners, maxpromiscuity, genetoname, genetoparalogs, synthinfo, isoreadsup_iso):
    fgenes = set(
        [x.split('.')[0] if x[:3] != 'chr' else x.split('-')[0] + '-' + str(round(int(x.split('-')[1]), -6)) for x in
         fusionchr.split('--')])
    ispromiscuous, areparalogs, areig, allnames = False, [], [], set()
    for g in fgenes:
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
            if g[:3] == 'chr':
                areig.append(True)
            else:
                areig.append(False)
        if g in genetoparalogs:
            other = fgenes - {g, }
            other = {genetoparalogs[g2] if g2 in genetoparalogs else g2 for g2 in other}
            if genetoparalogs[g] in other or len(other) == 0:
                areparalogs.append(True)
            else:
                areparalogs.append(False)
        else:
            areparalogs.append(False)
    return len(isoreadsup_iso) >= 1 and not ispromiscuous and any(x==False for x in areparalogs) and any(x==False for x in areig) and len(allnames) > 1 and 'chrM' not in [x[1] for x in synthinfo]

def get_locus_bounds(synthinfo):
    locuslen = [abs(x[3] - x[2]) for x in synthinfo]
    locusbounds = []
    laststart = 0
    for i in range(len(locuslen)):
        locusbounds.append((laststart, locuslen[i] + laststart))
        laststart += locuslen[i]
    return locusbounds

def separate_exons_by_locus(esizes, estarts, numloci, locusbounds, start):
    starts = [None for x in range(numloci)]
    exonindexes = [[] for x in range(numloci)]
    for i in range(len(esizes)):
        thisstart, thisend = start + estarts[i], start + estarts[i] + esizes[i]

        for order in range(numloci):
            if locusbounds[order][0] < thisstart and thisend <= locusbounds[order][1]:
                if starts[order] == None:
                    starts[order] = estarts[i]  # thisstart #- locusbounds[order][0]
                exonindexes[order].append(i)
    return starts, exonindexes

def convert_to_genomic_coords(numloci, synthinfo, exonindexes, starts, locusbounds, esizes, estarts, start, iso):
    genomicbounds = []
    outlines = []
    for order in range(numloci):
        genename, genomicchr, leftbound, rightbound = synthinfo[order]
        locusesizes = [esizes[i] for i in exonindexes[order]]
        locusestarts = [estarts[i] for i in exonindexes[order]]
        locusestarts = [x - starts[order] for x in locusestarts]
        totlen = locusestarts[-1] + locusesizes[-1]
        locusdir = '+' if leftbound < rightbound else '-'
        if leftbound > rightbound:  # reverse direction
            temp = []
            for i in range(len(locusestarts) - 1, -1, -1):
                temp.append(totlen - (locusestarts[i] + locusesizes[i]))
            locusestarts = temp
            locusesizes = locusesizes[::-1]
            outline = [genomicchr, str(leftbound - (((start + starts[order]) - locusbounds[order][0]) + totlen)),
                       str(leftbound - ((start + starts[order]) - locusbounds[order][0])),
                       'gene' + str(order + 1) + '_' + iso,
                       '1000', locusdir, str(leftbound - (start + totlen)), str(leftbound - start), '0',
                       str(len(locusesizes)),
                       ','.join([str(x) for x in locusesizes]),
                       ','.join([str(x) for x in locusestarts])]
            genomicbounds.append((genomicchr, leftbound - ((start + starts[order]) - locusbounds[order][0]),
                                  leftbound - (((start + starts[order]) - locusbounds[order][0]) + totlen), locusdir))
        else:
            outline = [genomicchr, str(leftbound + ((start + starts[order]) - locusbounds[order][0])),
                       str(leftbound + ((start + starts[order]) - locusbounds[order][0]) + totlen),
                       'gene' + str(order + 1) + '_' + iso,
                       '1000', locusdir, str(leftbound + ((start + starts[order]) - locusbounds[order][0])),
                       str(leftbound + ((start + starts[order]) - locusbounds[order][0]) + totlen), '0',
                       str(len(locusesizes)),
                       ','.join([str(x) for x in locusesizes]),
                       ','.join([str(x) for x in locusestarts])]
            genomicbounds.append((genomicchr, leftbound + ((start + starts[order]) - locusbounds[order][0]),
                                  leftbound + ((start + starts[order]) - locusbounds[order][0]) + totlen, locusdir))
        outlines.append('\t'.join(outline) + '\n')
    return genomicbounds, outlines

def check_too_close(numloci, genomicbounds):
    istooclose = False
    for i in range(numloci - 1):
        if genomicbounds[i][0] == genomicbounds[i + 1][0] and genomicbounds[i][3] == genomicbounds[i + 1][3] and abs(
                genomicbounds[i][2] - genomicbounds[i + 1][1]) < 350000:
            istooclose = True
    return istooclose

def write_final_fusion_reads(readsfile, freadsfinal):
    last = False
    temp = readsfile.split('.')
    temp[-2] = 'fusionreads'
    freads = open('.'.join(temp), 'w')

    if readsfile.split('.')[-1] == 'fasta' or readsfile.split('.')[-1] == 'fa':
        for line in open(readsfile):
            if line[0] == '>':
                readname = line[1:].rstrip().split()[0]
                if readname in freadsfinal:
                    last = True
                else:
                    last = False
            if last: freads.write(line)
    else:
        linecount = 0
        for line in open(readsfile):
            if linecount % 4 == 0:
                readname = line[1:].rstrip().split()[0]
                if readname in freadsfinal:
                    last = True
                else:
                    last = False
            if last: freads.write(line)
            linecount += 1
    freads.close()


def convert_synthetic_isos(annotgtf, isoformsbed, readmapfile, readsfile, breakpointfile, outname, paralogfile, maxpromiscuity):

    isoreadsup = get_iso_to_reads(readmapfile)
    synthchrtoinfo = get_synth_info(breakpointfile)
    genetoparalogs = get_paralog_ref(paralogfile)
    genetoname = get_gene_name_conv(annotgtf)
    locustopartners = identify_promiscuous_genes(isoformsbed, genetoparalogs)

    freadsfinal = set()
    out = open(outname, 'w')
    for line in open(isoformsbed):
        line = line.split('\t')
        iso, start, esizes, estarts = line[3], int(line[1]), [int(x) for x in line[10].split(',')[:-1]], [int(x) for x in line[11].split(',')[:-1]]
        fusionchr = line[0]
        synthinfo = [x.split('..') for x in synthchrtoinfo[fusionchr].split('--')]
        synthinfo = [[y[0], y[1], int(y[2]), int(y[3])] for y in synthinfo]

        is_good_fusion = identify_fusion_problems(fusionchr, locustopartners, maxpromiscuity, genetoname, genetoparalogs, synthinfo, isoreadsup[iso])

        if is_good_fusion:

            locusbounds = get_locus_bounds(synthinfo)
            if start < locusbounds[0][1] and locusbounds[-1][0] < int(line[2]):

                numloci = len(synthinfo)
                starts, exonindexes = separate_exons_by_locus(esizes, estarts, numloci, locusbounds, start)

                if None not in starts:
                    genomicbounds, outlines = convert_to_genomic_coords(numloci, synthinfo, exonindexes, starts, locusbounds, esizes, estarts, start, iso)
                    if not check_too_close(numloci, genomicbounds):
                        for l in outlines:
                            out.write(l)
                        freadsfinal.update(isoreadsup[iso])
    out.close()
    write_final_fusion_reads(readsfile, freadsfinal)



