VERSION = "2.2.0"
import os

def set_unix_path():
    "add programs in package to PATH."

    # FIXME: 2025-04-02 markd
    # this should be replaced with converting the exec use an explicit
    # path to make code more obvious.
    os.environ["PATH"] = os.path.dirname(os.path.realpath(__file__)) + ':' + os.environ["PATH"]

def check_diffexp_dependencies():
    """Called by wrapper programs for doing the different expression analysis to
    check if the optional dependencies are install,"""
    try:
        import rpy2
        import numpy
    except ModuleNotFoundError as ex:
        raise Exception("A FLAIR dependency for differential expression analysis is not installed.\n"
                        "See FLAIR documentation for information on installing the needed packages") from ex
