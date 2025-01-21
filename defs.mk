#
# definitions for Makefiles
#
# include with
#   root = ..
#   include ${root}/defs.mk

.PRECIOUS:
.SECONDARY:

SHELL = /bin/bash
export BASHOPTS = -beEu -o pipefail
MAKEFLAGS += -rR

POETRY = poetry
PYTHON = python3
FLAKE8 = ${python} -m flake8
