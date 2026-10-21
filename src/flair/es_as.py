#!/usr/bin/env python3

import sys
from flair import FlairInputDataError
from flair.diffsplice_event_ids import junction_id
from flair.es_events_tsv import EsEventsWriter
from flair.flair_bed import FlairBed
from flair.pycbio.hgdata.bed import BedReader
from flair.pycbio.sys import cli


class Gene(object):
    '''
    '''
    def __init__(self, gene_id=None, chromosome=None, strand=None):
        self.name = gene_id
        self.isoforms = dict()
        self.strand = strand
        self.chrom = chromosome
        if self.strand == "+":
            self.acceptor, self.donor = 0, 1
        else:
            self.acceptor, self.donor = 1, 0

    def build_graph(self):  # noqa: C901 - FIXME: reduce complexity
        self.exon_graph = dict()
        self.known_juncs = dict()
        self.splice_sites = dict()

        isos = sorted(list(self.isoforms.keys()))
        for i in isos:
            for num, j in enumerate(self.isoforms[i]):
                acceptor, donor = j[self.acceptor], j[self.donor]

                if num == 1 and len(self.isoforms[i]) == 2:
                    # two-exon isoform, only add the splice junction and not the exons
                    if acceptor not in self.splice_sites:
                        self.splice_sites[acceptor] = SpliceSite(acceptor, "acceptor")
                    previous_donor = self.isoforms[i][num - 1][self.donor]
                    if previous_donor not in self.splice_sites:
                        self.splice_sites[previous_donor] = SpliceSite(previous_donor, "donor")
                    j1 = (self.splice_sites[previous_donor], self.splice_sites[acceptor])
                    if j1 not in self.known_juncs:
                        self.known_juncs[j1] = set()
                    self.known_juncs[j1].add(i)
                    continue
                elif num == 0 or num + 1 >= len(self.isoforms[i]):
                    # first or last exons cannot be skipped, continue
                    continue

                previous_donor = self.isoforms[i][num - 1][self.donor]
                next_acceptor = self.isoforms[i][num + 1][self.acceptor]

                if acceptor not in self.splice_sites:
                    self.splice_sites[acceptor] = SpliceSite(acceptor, "acceptor")
                if donor not in self.splice_sites:
                    self.splice_sites[donor] = SpliceSite(donor, "donor")
                if next_acceptor not in self.splice_sites:
                    self.splice_sites[next_acceptor] = SpliceSite(next_acceptor, "acceptor")
                if previous_donor not in self.splice_sites:
                    self.splice_sites[previous_donor] = SpliceSite(previous_donor, "donor")

                acceptor, donor = self.splice_sites[acceptor], self.splice_sites[donor]
                previous_donor = self.splice_sites[previous_donor]
                next_acceptor = self.splice_sites[next_acceptor]

                exon = tuple(sorted([acceptor.name, donor.name]))
                if exon not in self.exon_graph:
                    self.exon_graph[exon] = Exon(exon, acceptor, donor)

                self.exon_graph[exon].inclusion_juncs.add(((previous_donor, acceptor), (donor, next_acceptor)))

                j1, j2 = (previous_donor, acceptor), (donor, next_acceptor)
                if j1 not in self.known_juncs:
                    self.known_juncs[j1] = set()
                if j2 not in self.known_juncs:
                    self.known_juncs[j2] = set()
                self.known_juncs[j1].add(i)
                self.known_juncs[j2].add(i)

    def find_skipped_exons(self, writer):
        '''
        '''

        self.writer = writer
        for e, obj in self.exon_graph.items():
            donor, acceptor = obj.donor, obj.acceptor
            if donor is None or acceptor is None:
                # first or last exon
                continue
            inclusion_isos = set()
            exclusion_isos = set()
            for juncs in obj.inclusion_juncs:
                j1, j2 = juncs
                if j1 in self.known_juncs and j2 in self.known_juncs:
                    inclusion_isos = inclusion_isos.union(self.known_juncs[j1].intersection(self.known_juncs[j2]))
                if (j1[0], j2[-1]) in self.known_juncs:
                    exclusion_isos = exclusion_isos.union(self.known_juncs[j1[0], j2[-1]])
            inclusion_isos = sorted(list(inclusion_isos))
            exclusion_isos = sorted(list(exclusion_isos))
            self.writer.writeRow((junction_id(self.chrom, self.strand, acceptor.name, donor.name),
                                  self.strand, len(inclusion_isos), len(exclusion_isos),
                                  inclusion_isos, exclusion_isos))
            # print(self.chrom, "\t".join(str(x) for x in sorted([acceptor.name,donor.name])), "%s:%s-%s" % (self.chrom,acceptor.name,donor.name), self.name, self.strand, sep="\t")


class SpliceSite(object):
    '''
    '''
    def __init__(self, ss_name=None, ss_type=None):
        self.name = ss_name
        self.ss_type = ss_type


class Exon(object):
    '''
    '''

    def __init__(self, exon_id=None, acceptor=None, donor=None):
        self.name = exon_id
        self.donor = donor
        self.acceptor = acceptor
        self.inclusion_juncs = set()


def check_gene_placement(gene, iso, chrom, strand, flair_isoforms):
    """All isoforms of a gene must sit on the same chromosome and strand.  Skipped
    exons are found by graphing a gene's isoforms together, and isoforms of two
    different loci graphed as one gene produce events that are in neither."""
    if (gene.chrom, gene.strand) != (chrom, strand):
        raise FlairInputDataError(
            f"{flair_isoforms}: isoform {iso} places gene {gene.name} on "
            f"{chrom}{strand}, an earlier isoform placed it on "
            f"{gene.chrom}{gene.strand}; give each locus its own gene name")


def main():
    with cli.ErrorHandler():
        find_skipped_exons()


def find_skipped_exons():

    flair_isoforms = sys.argv[1]
    genes = dict()
    for bed in BedReader(flair_isoforms, bedClass=FlairBed, fixScores=True):
        if bed.gene_id is None:
            raise FlairInputDataError(
                f"{flair_isoforms}: isoform {bed.name} has no gene_id column; exon "
                "skipping is called within a gene, so a FLAIR isoform BED is needed")
        exons = [(blk.start, blk.end) for blk in bed.blocks]

        if bed.gene_id not in genes:
            genes[bed.gene_id] = Gene(bed.gene_id, bed.chrom, bed.strand)
        gene = genes[bed.gene_id]
        check_gene_placement(gene, bed.name, bed.chrom, bed.strand, flair_isoforms)
        gene.isoforms[bed.name] = exons

    with EsEventsWriter(None, outFh=sys.stdout) as writer:
        for gobj in genes.values():
            gobj.build_graph()
            gobj.find_skipped_exons(writer)


if __name__ == "__main__":
    main()
