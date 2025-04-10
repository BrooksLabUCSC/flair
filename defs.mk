#
# definitions for Makefiles
#
# include with
#   root = ..
#   include ${root}/defs.mk

.PRECIOUS:
.SECONDARY:

VERSION = 2.1.1
PACKAGE_NAME = flair-brookslab-${VERSION}
PACKAGE_FILE_PREFIX = flair_brookslab-${VERSION}

SHELL = /bin/bash
export BASHOPTS = -beEu -o pipefail
MAKEFLAGS += -rR

POETRY = poetry
PYTHON = python3
FLAKE8 = ${python} -m flake8
