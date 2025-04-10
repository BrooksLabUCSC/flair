VERSION = "2.1.1"
import os

def set_unix_path():
    # get programs in package in PATH.

    # FIXME: 2025-04-02 markd
    # this should be replaced with converting the exec use an explict
    # path to make code more obvious.
    os.environ["PATH"] = os.path.dirname(__file__) + ':' + os.environ["PATH"]
