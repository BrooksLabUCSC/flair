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

poetry = poetry
python = ${poetry} run python3
flake8 = ${python} -m flake8
