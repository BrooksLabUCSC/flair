#!/usr/bin/env python3
"""Call skipped exon events from a FLAIR isoform BED, writing one row per
internal exon with the isoforms that include it and the isoforms that splice
over it.  See es_events_tsv for the columns."""

import argparse
import sys
from flair import FlairInputDataError
from flair.cli_help import FlairHelpFormatter
from flair.diffsplice_event_ids import junction_id
from flair.es_events_tsv import EsEventsWriter
from flair.flair_bed import FlairBedReader
from flair.pycbio.sys import cli


class SpliceSite(object):
    "one splice site of a gene, identified by its position"

    def __init__(self, ss_name=None, ss_type=None):
        self.name = ss_name
        self.ss_type = ss_type


class Exon(object):
    "an internal exon and the junction pairs that include it"

    def __init__(self, exon_id=None, acceptor=None, donor=None):
        self.name = exon_id
        self.donor = donor
        self.acceptor = acceptor
        self.inclusion_juncs = set()


class Gene(object):
    """The isoforms of one gene and the splice site graph built from them.
    Exons are held 5' to 3', so neighbors in an isoform's exon list are
    neighbors in transcription order on either strand."""

    def __init__(self, gene_id=None, chromosome=None, strand=None):
        self.name = gene_id
        self.isoforms = dict()
        self.strand = strand
        self.chrom = chromosome
        if self.strand == "+":
            self.acceptor, self.donor = 0, 1
        else:
            self.acceptor, self.donor = 1, 0
        self.exon_graph = dict()
        self.known_juncs = dict()
        self.splice_sites = dict()

    def add_isoform(self, name, exons):
        "exons are (start, end) pairs in genomic order"
        self.isoforms[name] = exons if self.strand == "+" else list(reversed(exons))

    def build_graph(self):
        for name in sorted(self.isoforms.keys()):
            self._add_isoform_to_graph(name)

    def _add_isoform_to_graph(self, name):
        exons = self.isoforms[name]
        if len(exons) == 2:
            # neither exon can be skipped, so record only the junction
            self._add_junction(name, exons[0], exons[1])
        else:
            for num in range(1, len(exons) - 1):
                self._add_internal_exon(name, exons[num - 1], exons[num], exons[num + 1])

    def _splice_site(self, pos, ss_type):
        if pos not in self.splice_sites:
            self.splice_sites[pos] = SpliceSite(pos, ss_type)
        return self.splice_sites[pos]

    def _add_junction(self, name, prev_exon, next_exon):
        "record the junction between two consecutive exons and return it"
        junc = (self._splice_site(prev_exon[self.donor], "donor"),
                self._splice_site(next_exon[self.acceptor], "acceptor"))
        self.known_juncs.setdefault(junc, set()).add(name)
        return junc

    def _add_internal_exon(self, name, prev_exon, exon, next_exon):
        juncs = (self._add_junction(name, prev_exon, exon),
                 self._add_junction(name, exon, next_exon))
        acceptor = self._splice_site(exon[self.acceptor], "acceptor")
        donor = self._splice_site(exon[self.donor], "donor")
        exon_id = tuple(sorted([acceptor.name, donor.name]))
        if exon_id not in self.exon_graph:
            self.exon_graph[exon_id] = Exon(exon_id, acceptor, donor)
        self.exon_graph[exon_id].inclusion_juncs.add(juncs)

    def find_skipped_exons(self, writer):
        for exon in self.exon_graph.values():
            self._write_exon_event(writer, exon)

    def _write_exon_event(self, writer, exon):
        inclusion_isos = set()
        exclusion_isos = set()
        for j1, j2 in exon.inclusion_juncs:
            inclusion_isos |= self._isos_with_both_juncs(j1, j2)
            exclusion_isos |= self.known_juncs.get((j1[0], j2[1]), set())
        writer.writeRow((junction_id(self.chrom, self.strand, exon.acceptor.name, exon.donor.name),
                         self.strand, len(inclusion_isos), len(exclusion_isos),
                         sorted(inclusion_isos), sorted(exclusion_isos)))

    def _isos_with_both_juncs(self, j1, j2):
        if (j1 in self.known_juncs) and (j2 in self.known_juncs):
            return self.known_juncs[j1] & self.known_juncs[j2]
        else:
            return set()


def check_gene_placement(gene, iso, chrom, strand, flair_isoforms):
    """All isoforms of a gene must sit on the same chromosome and strand.  Skipped
    exons are found by graphing a gene's isoforms together, and isoforms of two
    different loci graphed as one gene produce events that are in neither."""
    if (gene.chrom, gene.strand) != (chrom, strand):
        raise FlairInputDataError(
            f"{flair_isoforms}: isoform {iso} places gene {gene.name} on "
            f"{chrom}{strand}, an earlier isoform placed it on "
            f"{gene.chrom}{gene.strand}; give each locus its own gene name")


def read_genes(flair_isoforms):
    genes = dict()
    for bed in FlairBedReader(flair_isoforms, fixScores=True):
        if bed.gene_id is None:
            raise FlairInputDataError(
                f"{flair_isoforms}: isoform {bed.name} has no gene_id column; exon "
                "skipping is called within a gene, so a FLAIR isoform BED is needed")
        if bed.gene_id not in genes:
            genes[bed.gene_id] = Gene(bed.gene_id, bed.chrom, bed.strand)
        gene = genes[bed.gene_id]
        check_gene_placement(gene, bed.name, bed.chrom, bed.strand, flair_isoforms)
        gene.add_isoform(bed.name, [(blk.start, blk.end) for blk in bed.blocks])
    return genes


def call_skipped_exons(flair_isoforms, out_fh):
    genes = read_genes(flair_isoforms)
    with EsEventsWriter(None, outFh=out_fh) as writer:
        for gene in genes.values():
            gene.build_graph()
            gene.find_skipped_exons(writer)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=FlairHelpFormatter)
    parser.add_argument('flair_isoforms',
                        help="FLAIR isoform BED, with a gene_id column")
    return parser


def main():
    with cli.ErrorHandler():
        opts = build_parser().parse_args()
        call_skipped_exons(opts.flair_isoforms, sys.stdout)


if __name__ == "__main__":
    main()
