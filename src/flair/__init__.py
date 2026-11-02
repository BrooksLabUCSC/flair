import os
from flair.pycbio import NoStackError

VERSION = "3.0.1+master"

class FlairError(Exception):
    """General error condition in FLAIR"""
    pass

class FlairInputDataError(FlairError, NoStackError):
    """Error in FLAIR input data"""
    pass


# the IUPAC nucleotide codes to their complements, with U and the gap characters:
# https://www.bioinformatics.org/sms/iupac.html
COMPBASE = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
            'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N',
            '.': '.', '-': '-'}

def complement_base(base):
    try:
        return COMPBASE[base]
    except KeyError as ex:
        raise FlairInputDataError(f"not a nucleotide: '{base}'; sequences must be the IUPAC "
                                  f"codes ACGTURYSWKMBDHVN or a gap (. or -)") from ex

def reverse_complement(seq):
    "the reverse complement, upper case, of a sequence of IUPAC nucleotide codes"
    return ''.join([complement_base(base) for base in reversed(seq.upper())])

def set_unix_path():
    "add programs in package to PATH."

    # FIXME: 2025-04-02 markd
    # this should be replaced with converting the exec use an explicit
    # path to make code more obvious.
    os.environ["PATH"] = os.path.dirname(os.path.realpath(__file__)) + ':' + os.environ["PATH"]
