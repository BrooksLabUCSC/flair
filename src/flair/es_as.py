#!/usr/bin/env python3

import sys
from flair import FlairInputDataError
from flair.es_events_tsv import EsEventsWriter
from flair.flair_bed import FlairBed
from flair.pycbio.hgdata.bed import BedReader


class Gene(object):
    '''
    '''
    def __init__(self, geneID=None, chromosome=None, strand=None):
        self.name = geneID
        self.isoforms = dict()
        self.strand = strand
        self.chrom = chromosome
        if self.strand == "+":
            self.acceptor, self.donor = 0, 1
        else:
            self.acceptor, self.donor = 1, 0

    def build_graph(self):  # noqa: C901 - FIXME: reduce complexity
        self.exonGraph = dict()
        self.knownJuncs = dict()
        self.spliceSites = dict()

        isos = sorted(list(self.isoforms.keys()))
        for i in isos:
            for num, j in enumerate(self.isoforms[i]):
                acceptor, donor = j[self.acceptor], j[self.donor]

                if num == 1 and len(self.isoforms[i]) == 2:
                    # two-exon isoform, only add the splice junction and not the exons
                    if acceptor not in self.spliceSites:
                        self.spliceSites[acceptor] = SpliceSite(acceptor, "acceptor")
                    previousDonor = self.isoforms[i][num - 1][self.donor]
                    if previousDonor not in self.spliceSites:
                        self.spliceSites[previousDonor] = SpliceSite(previousDonor, "donor")
                    j1 = (self.spliceSites[previousDonor], self.spliceSites[acceptor])
                    if j1 not in self.knownJuncs:
                        self.knownJuncs[j1] = set()
                    self.knownJuncs[j1].add(i)
                    continue
                elif num == 0 or num + 1 >= len(self.isoforms[i]):
                    # first or last exons cannot be skipped, continue
                    continue

                previousDonor = self.isoforms[i][num - 1][self.donor]
                nextAcceptor = self.isoforms[i][num + 1][self.acceptor]

                if acceptor not in self.spliceSites:
                    self.spliceSites[acceptor] = SpliceSite(acceptor, "acceptor")
                if donor not in self.spliceSites:
                    self.spliceSites[donor] = SpliceSite(donor, "donor")
                if nextAcceptor not in self.spliceSites:
                    self.spliceSites[nextAcceptor] = SpliceSite(nextAcceptor, "acceptor")
                if previousDonor not in self.spliceSites:
                    self.spliceSites[previousDonor] = SpliceSite(previousDonor, "donor")

                acceptor, donor = self.spliceSites[acceptor], self.spliceSites[donor]
                previousDonor = self.spliceSites[previousDonor]
                nextAcceptor = self.spliceSites[nextAcceptor]

                exon = tuple(sorted([acceptor.name, donor.name]))
                if exon not in self.exonGraph:
                    self.exonGraph[exon] = Exon(exon, acceptor, donor)

                self.exonGraph[exon].inclusionJuncs.add(((previousDonor, acceptor), (donor, nextAcceptor)))

                j1, j2 = (previousDonor, acceptor), (donor, nextAcceptor)
                if j1 not in self.knownJuncs:
                    self.knownJuncs[j1] = set()
                if j2 not in self.knownJuncs:
                    self.knownJuncs[j2] = set()
                self.knownJuncs[j1].add(i)
                self.knownJuncs[j2].add(i)

    def find_skipped_exons(self, writer):
        '''
        '''

        self.writer = writer
        for e, obj in self.exonGraph.items():
            donor, acceptor = obj.donor, obj.acceptor
            if donor is None or acceptor is None:
                # first or last exon
                continue
            inclusionIsos = set()
            exclusionIsos = set()
            for juncs in obj.inclusionJuncs:
                j1, j2 = juncs
                if j1 in self.knownJuncs and j2 in self.knownJuncs:
                    inclusionIsos = inclusionIsos.union(self.knownJuncs[j1].intersection(self.knownJuncs[j2]))
                if (j1[0], j2[-1]) in self.knownJuncs:
                    exclusionIsos = exclusionIsos.union(self.knownJuncs[j1[0], j2[-1]])
            inclusionIsos = sorted(list(inclusionIsos))
            exclusionIsos = sorted(list(exclusionIsos))
            self.writer.writeRow(("%s:%s-%s" % (self.chrom, acceptor.name, donor.name),
                                  self.strand, len(inclusionIsos), len(exclusionIsos),
                                  inclusionIsos, exclusionIsos))
            # print(self.chrom, "\t".join(str(x) for x in sorted([acceptor.name,donor.name])), "%s:%s-%s" % (self.chrom,acceptor.name,donor.name), self.name, self.strand, sep="\t")


class SpliceSite(object):
    '''
    '''
    def __init__(self, ssName=None, ssType=None):
        self.name = ssName
        self.ssType = ssType


class Exon(object):
    '''
    '''

    def __init__(self, exonID=None, acceptor=None, donor=None):
        self.name = exonID
        self.donor = donor
        self.acceptor = acceptor
        self.inclusionJuncs = set()


def main():

    flairIsoforms = sys.argv[1]
    genes = dict()
    for bed in BedReader(flairIsoforms, bedClass=FlairBed, fixScores=True):
        if bed.gene_id is None:
            raise FlairInputDataError(
                f"{flairIsoforms}: isoform {bed.name} has no gene_id column; exon "
                "skipping is called within a gene, so a FLAIR isoform BED is needed")
        exons = [(blk.start, blk.end) for blk in bed.blocks]

        if bed.gene_id not in genes:
            genes[bed.gene_id] = Gene(bed.gene_id, bed.chrom, bed.strand)
        genes[bed.gene_id].isoforms[bed.name] = exons

    with EsEventsWriter(None, outFh=sys.stdout) as writer:
        for gobj in genes.values():
            gobj.build_graph()
            gobj.find_skipped_exons(writer)


if __name__ == "__main__":
    main()
