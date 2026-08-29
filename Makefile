root = .
include ${root}/defs.mk

# If you want to use the install flair rather than the one in this
# tree use:
#   make test use_installed_flair=yes
#
# or use
#   make test-installed
#

# flake8 walks the tree for *.py itself; PYPROGS adds executable scripts that
# have no .py extension.  Symbolic links are not reported by file(1), so the
# bin/ commands are not checked twice, only their targets in src/flair.
PYPROGS = $(shell file -F $$'\t' bin/* test/bin/* | awk -F'\t' '/Python script/{print $$1}')

FLAKE8_CHECK = . ${PYPROGS}

default:

doc:
	${MAKE} -C docs html

##
# test targets, the xx-installed test with the installed FLAIR rather than the
# tree.
##
.PHONY: test
test:
	${MAKE} -C test test use_installed_flair=${use_installed_flair}
test-base:
	${MAKE} -C test test-base use_installed_flair=${use_installed_flair}
test-installed:
	${MAKE} -C test test-installed
test-base-installed:
	${MAKE} -C test test-base-installed

# release testing: tests with only the conda environment and system directories
# in PATH, see test/Makefile
test-env-only:
	${MAKE} -C test test-env-only use_installed_flair=${use_installed_flair}


##
# lint check with flake8
#   see .flake8 for configuration and the list of excluded directories
##
lint: flake8
pycbio-lint: pycbio-flake8

flake8:
	${FLAKE8} --color=never ${FLAKE8_CHECK}

PYCBIO_DIR = src/flair/pycbio
pycbio-flake8:
	${FLAKE8} --color=never --config=${PYCBIO_DIR}/setup.cfg ${PYCBIO_DIR}


##
# test environment for pip install
##
pip_test_env = pip_test_env

PYPI_UPLOAD_URL = https://upload.pypi.org/legacy/
PYPI_INSTALL_URL = https://pypi.org/simple/
TESTPYPI_UPLOAD_URL = https://test.pypi.org/legacy/
TESTPYPI_INSTALL_URL = https://test.pypi.org/simple/

define pip_env_setup
	rm -rf ${pip_test_env}
	mkdir -p ${pip_test_env}
	${PYTHON} -m virtualenv --quiet ${pip_test_env}
endef
pip_env_act = source ${pip_test_env}/bin/activate


##
# release targets
##
build: clean
	poetry build

# test if pip install locally
test-pip:
	${pip_env_setup}
	${pip_env_act} && pip install --no-cache-dir ./dist/${PACKAGE_FILE_PREFIX}-py3-none-any.whl
	${pip_env_act} && ${MAKE} -C test test use_installed_flair=yes

# testpypy
publish-testpypi: build
	poetry publish -r testpypi

test-testpypi:
	${pip_env_setup}
	${pip_env_act} && pip install --no-cache-dir  --index-url=${TESTPYPI_INSTALL_URL} --extra-index-url=${PYPI_INSTALL_URL} flair-brookslab==${VERSION}
	${pip_env_act} && ${MAKE} -C test test use_installed_flair=yes


# pypy
publish-pypi: build
	poetry publish

test-pypi:
	${pip_env_setup}
	${pip_env_act} && pip install --no-cache-dir  --index-url=${PYPI_INSTALL_URL} flair-brookslab==${VERSION}
	${pip_env_act} && ${MAKE} -C test test use_installed_flair=yes


##
# clean targets
##
clean:
	rm -rf build/ dist/ ${pip_test_env}/ src/flair/__pycache__/
	cd test && ${MAKE} clean

real-clean: clean
	cd test && ${MAKE} real-clean
	${MAKE} -C docs clean


