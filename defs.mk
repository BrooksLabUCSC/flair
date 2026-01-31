#
# definitions for Makefiles
#
# include with
#   root = ..
#   include ${root}/defs.mk

.PRECIOUS:
.SECONDARY:

VERSION = 3.0.0+dev
PACKAGE_NAME = flair-brookslab-${VERSION}
PACKAGE_FILE_PREFIX = flair_brookslab-${VERSION}

SHELL = /bin/bash
export BASHOPTS = -beEu -o pipefail
MAKEFLAGS += -rRw

POETRY = poetry
PYTHON = python3
FLAKE8 = ${PYTHON} -m flake8
