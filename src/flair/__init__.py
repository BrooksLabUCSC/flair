import os
from flair.pycbio import NoStackError

VERSION = "3.0.0"

class FlairError(Exception):
    """General error condition in FLAIR"""
    pass

class FlairInputDataError(FlairError, NoStackError):
    """Error in FLAIR input data"""
    pass

def set_unix_path():
    "add programs in package to PATH."

    # FIXME: 2025-04-02 markd
    # this should be replaced with converting the exec use an explicit
    # path to make code more obvious.
    os.environ["PATH"] = os.path.dirname(os.path.realpath(__file__)) + ':' + os.environ["PATH"]
