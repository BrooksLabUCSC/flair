"""The homopolymer filter of the indel caller, which measures the run of bases in the
genome on the 3' side of an indel"""
from flair.identify_vars import check_homopoly


class _Genome:
    "a genome of one sequence, read as pysam reads a FASTA"
    def __init__(self, seq):
        self.seq = seq

    def fetch(self, chrom, start, end):
        return self.seq[max(start, 0):end]


#            0-based  0         1
#                     012345678901
GENOME = _Genome('TTTTCAAAGTTT')
# the C is one-based 5, and the run of three A is one-based 6 to 8


def test_an_insertion_extending_a_run_of_three_is_a_homopolymer_at_any_vaf():
    # REF C, ALT CA: the A is inserted after the C, so the run to its right is AAA
    assert check_homopoly('INS', 'CA', 0.5, GENOME, 'chr1', 5)


def test_a_deletion_from_the_run_is_measured_from_after_the_deleted_base():
    # REF CA, ALT C: one A deleted, leaving two of them to the right
    assert not check_homopoly('DEL', 'CA', 0.5, GENOME, 'chr1', 5)
    assert check_homopoly('DEL', 'CA', 0.1, GENOME, 'chr1', 5)


def test_an_insertion_of_another_base_is_not_a_homopolymer():
    assert not check_homopoly('INS', 'CG', 0.5, GENOME, 'chr1', 5)


def test_an_insertion_of_two_different_bases_is_not_measured():
    assert not check_homopoly('INS', 'CAG', 0.5, GENOME, 'chr1', 5)


def test_an_indel_at_the_end_of_the_sequence_has_nothing_to_its_right():
    assert not check_homopoly('INS', 'TT', 0.5, GENOME, 'chr1', len(GENOME.seq))
