VERSION = "2.1.0"
import os

def set_unix_path():
    # get programs in package in PATH.
    #
    os.environ["PATH"] = os.path.dirname(__file__) + ':' + os.environ["PATH"]
