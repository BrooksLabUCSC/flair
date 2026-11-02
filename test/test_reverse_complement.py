"""
tests of the reverse complement of isoform_data
"""
import pytest
from flair import FlairInputDataError
from flair.isoform_data import COMPBASE, get_reverse_complement

# https://www.bioinformatics.org/sms/iupac.html, with U and the gap characters
IUPAC_CODES = 'ACGTURYSWKMBDHVN.-'


def test_plain_bases():
    assert get_reverse_complement('ACGTACGT') == 'ACGTACGT'
    assert get_reverse_complement('AAACCC') == 'GGGTTT'


def test_every_iupac_code():
    assert get_reverse_complement(IUPAC_CODES) == '-.NBDHVKMWSRYAACGT'


def test_lower_case_is_upper_cased():
    assert get_reverse_complement('acgtryn') == 'NRYACGT'


def test_uracil_complements_to_adenine():
    assert get_reverse_complement('U') == 'A'


def test_gaps_are_kept():
    assert get_reverse_complement('A-C.G') == 'C.G-T'


def test_complement_of_complement_is_the_sequence():
    assert get_reverse_complement(get_reverse_complement(IUPAC_CODES)) == IUPAC_CODES.replace('U', 'T')


def test_every_code_has_a_complement_that_is_itself_a_code():
    assert set(COMPBASE) == set(IUPAC_CODES)
    assert set(COMPBASE.values()) <= set(IUPAC_CODES)


def test_other_character_names_itself():
    with pytest.raises(FlairInputDataError, match="not a nucleotide: 'Z'"):
        get_reverse_complement('ACGZ')
