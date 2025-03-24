#
# definitions for Makefiles
#
# include with
#   root = ..
#   include ${root}/defs.mk

.PRECIOUS:
.SECONDARY:

VERSION = 2.0.0

SHELL = /bin/bash
export BASHOPTS = -beEu -o pipefail
MAKEFLAGS += -rR

POETRY = poetry
PYTHON = python3
FLAKE8 = ${python} -m flake8
