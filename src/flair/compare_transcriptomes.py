import sys

###Developed Jayni Ram, Alex Barnes, and the rest of their groups from the BME160 Winter 2025 class
###With support from Colette Felton

###usage: python3 compare_transcriptomes.py MANE MANE_isos.bed your_isos.bed output.bed
###The MANE mode assumes that your file 1 has exactly one isoform per gene
###OR
###usage: python3 compare_transcriptomes.py compfiles isofile1.bed isofile2.bed output.bed
###This mode does not assume anything about the transcriptomes, and tries to match each isoform in file 2 to an isoform in file 1
###I reccommend running this mode both with file1 vs file2 and file2 vs file1 to get the full picture
###The different modes also output different summary metrics

class Isoform:

    def __init__(self, fields):
        chrom = fields[0]
        chromStart = int(fields[1])
        chromEnd = int(fields[2])
        name = fields[3]  # isoform_gene = nomenclature
        gene = name.split('_')[-1]
        # can have multiple parts to the name separated by _
        isoform_name = '_'.join(name.split('_')[:-1])
        strand = fields[5]

        blockSizes = fields[10]
        blockSizes = blockSizes.rstrip(',').split(',')  # each number is stored as a String datatype
        blockSizes = [int(x) for x in blockSizes]  # converts to integer

        blockStarts = fields[11]
        blockStarts = blockStarts.rstrip(',').split(',')
        blockStarts = [int(x) for x in blockStarts]


        self.chrom=chrom
        self.chromStart=chromStart
        self.chromEnd=chromEnd
        self.name=name
        self.gene=gene
        self.isoform_name=isoform_name
        self.strand=strand
        self.blockStarts=blockStarts
        self.blockSizes=blockSizes
        self.exons = self.get_exons()
        self.junctions = []
        for i in range(len(self.exons)-1):
            self.junctions.append(self.exons[i][1])
            self.junctions.append(self.exons[i+1][0])
        self.fields = fields


    def get_exons(self):
        """Extracts exon start and end positions."""
        exons = []
        for i in range(len(self.blockStarts)):
            exon_start = self.chromStart + self.blockStarts[i]
            exon_end = exon_start + self.blockSizes[i]
            exons.append((exon_start, exon_end))
        if self.strand == '-':
            exons = exons[::-1]
            exons = [x[::-1] for x in exons]
        return exons


def process_file(file):
    isoinfo = {}
    for line in open(file):
        fields = line.rstrip().split('\t')  # list created
        isoform = Isoform(fields)
        if isoform.gene not in isoinfo:
            isoinfo[isoform.gene] = []
        isoinfo[isoform.gene].append(isoform)

    return isoinfo


def compareIndivIsos(isoform, mane_isoform):
    isoresults = []
    diffscore = 0
    shared_edges = set(isoform.junctions) & set(mane_isoform.junctions)
    diffjunccount = (max(len(isoform.junctions), len(mane_isoform.junctions)) - len(
        shared_edges))  # abs(len(isoform.junctions) - len(mane_isoform.junctions)) +
    if diffjunccount <= 5:
        refexons = mane_isoform.exons
        isoexons = isoform.exons
        refindex, isoindex = 0, 0
        while isoindex < len(isoexons):
            if refindex >= len(refexons):
                isoresults.append('additionalexonattend')
                diffscore += 2
            else:
                restart, reend = refexons[refindex]
                iestart, ieend = isoexons[isoindex]
                if restart != iestart and reend != ieend:
                    if refindex == len(refexons) - 1 and isoindex == len(isoexons) - 1:
                        isoresults.append('altlastexon')
                        diffscore += 3
                    elif isoindex == len(isoexons) - 1:
                        isoresults.append('altlastexon')
                    elif refexons[refindex] == isoexons[isoindex + 1] or refexons[refindex][0] == \
                            isoexons[isoindex + 1][0] or refexons[refindex][1] == \
                            isoexons[isoindex + 1][1]:
                        isoresults.append('insertedexonaftere' + str(refindex))
                        diffscore += 3
                        isoindex += 1
                        continue
                    elif refindex == len(refexons) - 1:
                        isoresults.append('altlastexon')
                        diffscore += 3
                    elif refexons[refindex + 1] == isoexons[isoindex] or refexons[refindex + 1][0] == \
                            isoexons[isoindex][0] or refexons[refindex + 1][1] == isoexons[isoindex][1]:
                        isoresults.append('skippedexon' + str(refindex + 1))
                        diffscore += 2
                        refindex += 1
                        continue
                    elif refexons[refindex + 1] == isoexons[isoindex + 1]:
                        if isoindex == 0:
                            isoresults.append('altfirstexon')
                            diffscore += 3
                        else:
                            isoresults.append('altexon' + str(refindex + 1))
                            diffscore += 3
                else:
                    if restart != iestart:
                        if isoindex == 0:
                            isoresults.append('altstart')
                            diffscore += 0.4
                        else:
                            isoresults.append('exon' + str(refindex + 1) + 'altacceptor')
                            diffscore += 1
                    elif reend != ieend:
                        ###check for intron retention
                        if refindex < len(refexons) - 1 and ieend == refexons[refindex + 1][1]:
                            isoresults.append('e' + str(refindex + 1) + 'toe' + str(refindex + 2) + 'intronretention')
                            diffscore += 2
                            isoindex += 1
                            refindex += 2
                            continue
                        else:
                            if isoindex == len(isoexons) - 1:
                                isoresults.append('altend')
                                diffscore += 0.4
                            else:
                                isoresults.append('exon' + str(refindex + 1) + 'altdonor')
                                diffscore += 1
            refindex += 1
            isoindex += 1
    else:
        isoresults.append('verydifferent')
        diffscore = 30

    return isoresults, diffscore

def getsimpleres(isoresults, removestartend=True):
    simpleisoresults = isoresults.copy()
    if removestartend:
        if 'altstart' in simpleisoresults: simpleisoresults.remove('altstart')
        if 'altend' in simpleisoresults: simpleisoresults.remove('altend')
    if removestartend and len(simpleisoresults) > 1: simpleisoresults = ['multiplesplicingevents']
    elif not removestartend and len(set(simpleisoresults)-{'altstart', 'altend'}) > 1: simpleisoresults = ['multiplesplicingevents']
    ###Add start and end change back maybe, not sure this option is important
    # elif not removestartend and len(set(simpleisoresults)-{'altstart', 'altend'}) == 1: simpleisoresults = list(set(simpleisoresults)-{'altstart', 'altend'})
    simpleisoresults = ','.join(sorted(simpleisoresults))
    if simpleisoresults == "": simpleisoresults = "refmatch"
    for i in range(10):
        simpleisoresults = simpleisoresults.replace(str(i), 'X')
    simpleisoresults = simpleisoresults.replace('XX', 'X')
    return simpleisoresults

def compare_isoforms_to_ref(file1, file2):
    """
    Compares each isoform against the corresponding reference isoform (file1).
    Records differences in start/end positions, exon structure, and splice sites.
    """
    results = []
    summary = {}
    for gene_id, isoforms in file2.items():
        if gene_id not in file1:
            for isoform in isoforms:
                isoresults = 'nogenematch'
                new_name = isoresults + '_' + isoform.name  # combines this information together in the name
                isoform.fields[3] = new_name
                results.append(isoform.fields)
                if isoresults not in summary: summary[isoresults] = 1
                else: summary[isoresults] += 1
        else:
            for isoform in isoforms:
                isoresults = []  # only initially contains splice site info
                for mane_isoform in file1[gene_id]:  # Reference isoforms for the gene
                    isoresults, diffscore = compareIndivIsos(isoform, mane_isoform)
                    break  # only gets first reference isoform and runs loop once
                simpleisoresults = getsimpleres(isoresults)
                isoresults = ','.join(isoresults)
                if isoresults == "": isoresults = "refmatch"
                if simpleisoresults not in summary: summary[simpleisoresults] = 1
                else: summary[simpleisoresults] += 1
                new_name = isoresults + '_' + isoform.name  # combines this information together in the name
                isoform.fields[3] = new_name
                results.append(isoform.fields)
    return results, summary

def compare_isoforms(file1, file2):
    results = []
    summary = {}
    file1isostofile2isos = {'nogenematch':[]}
    for gene_id, isoforms in file2.items():
        if gene_id not in file1:
            for isoform in isoforms:
                isoresults = 'nogenematch'
                new_name = isoresults + '_' + isoform.name  # combines this information together in the name
                isoform.fields[3] = new_name
                results.append(isoform.fields[3])
                if isoresults not in summary: summary[isoresults] = 1
                else: summary[isoresults] += 1
                file1isostofile2isos['nogenematch'].append(isoform.name)
        else:
            for isoform in isoforms:
                isocomps = []
                for mane_isoform in file1[gene_id]:  # Reference isoforms for the gene
                    isoresults, diffscore = compareIndivIsos(isoform, mane_isoform)
                    isocomps.append((diffscore, isoresults, mane_isoform.name))
                isocomps.sort()
                isoresults = isocomps[0][1]
                file1iso = isocomps[0][2]

                # simpleisoresults = getsimpleres(isoresults, removestartend=False)
                supersimpleres = getsimpleres(isoresults)
                isoresults = ','.join(isoresults)
                if isoresults == "": isoresults = "refmatch"
                if gene_id not in file1isostofile2isos: file1isostofile2isos[gene_id] = {'noisomatch':[]}

                if isoresults == 'verydifferent': file1isostofile2isos[gene_id]['noisomatch'].append(isoform.name)
                else:
                    if file1iso not in file1isostofile2isos[gene_id]: file1isostofile2isos[gene_id][file1iso] = []
                    if isoresults == 'refmatch': file1isostofile2isos[gene_id][file1iso].append(('refmatch', isoform.name))
                    elif supersimpleres == 'refmatch': file1isostofile2isos[gene_id][file1iso].append(('juststartendchange', isoform.name))
                    elif supersimpleres == 'multiplesplicingevents': file1isostofile2isos[gene_id][file1iso].append(('multiplesplicingevents', isoform.name))
                    else: file1isostofile2isos[gene_id][file1iso].append(('onesplicingevent', isoform.name))


                # if simpleisoresults not in summary: summary[simpleisoresults] = 1
                # else: summary[simpleisoresults] += 1
                new_name = isoresults + '_' + isoform.name  # combines this information together in the name
                isoform.fields[3] = new_name
                results.append(isoform.fields[3])
    allf1tof2comps = {'noisomatch':0}
    for gene in file1isostofile2isos:
        if gene == 'nogenematch': allf1tof2comps['nogenematch'] = len(file1isostofile2isos['nogenematch'])
        else:
            for iso in file1isostofile2isos[gene]:
                if iso == 'noisomatch': allf1tof2comps['noisomatch'] += len(file1isostofile2isos[gene][iso])
                else:
                    allisoevents = ','.join(sorted([x[0] for x in file1isostofile2isos[gene][iso]], reverse=True))
                    if allisoevents not in allf1tof2comps: allf1tof2comps[allisoevents] = 0
                    allf1tof2comps[allisoevents] += 1
    return results, allf1tof2comps


def main():

    # mode = 'compfiles'
    # file1 = 'ogflair.031325.isoforms.bed'
    # file2 = '031725.isoforms.bed'
    # output_file = 'test_032625_output.bed'

    mode = 'MANE'
    file1 = 'final/gencode.v38.annotation.MANEselectonly.bed'
    file2 = 'ogflair.031325.isoforms.bed'
    output_file = 'test_032625_output_MANE.bed'

    # mode = sys.argv[1]
    # file1 = sys.argv[2]
    # file2 = sys.argv[3]
    # output_file = sys.argv[4]


    refisos = process_file(file1)

    file2isos = process_file(file2)


    with open(output_file, "w") as file:  # opens file to write (w) in it
        if mode == 'MANE': results, summary = compare_isoforms_to_ref(refisos, file2isos)
        else: results, summary = compare_isoforms(refisos, file2isos)
        for i in results:  # writes out results of comparison_results
            file.write('\t'.join(i) + "\n")

    for key in summary:
        print(summary[key], key)  # prints both out with space between


if __name__ == "__main__":
    main()