root = .
include ${root}/defs.mk

# If you want to use the install flair rather than the one in this
# tree use:
#   make test use_installed_flair=yes
#
# or use
#   make test-installed

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
test-installed:
	${MAKE} -C test test use_installed_flair=yes

test-diffexpress:
	${MAKE} -C test test-diffexpress use_installed_flair=${use_installed_flair}
test-diffexpress-installed:
	${MAKE} -C test test-diffexpress use_installed_flair=yes


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


