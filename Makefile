root = .
include ${root}/defs.mk

# If you want to use the install flair rather than the one in this
# tree use:
#   make test use_installed_flair=yes


default:

doc:
	${MAKE} -C docs html

##
# test targets
##
.PHONY: test
test:
	${MAKE} -C test test use_installed_flair=${use_installed_flair}

test-diffexp:
	${MAKE} -C test test-diffexp use_installed_flair=${use_installed_flair}


##
# test environment for pip install
##
pip_test_env = pip_test_env

pypi_url = https://upload.pypi.org/legacy/
testpypi_url = https://test.pypi.org/legacy/

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
	${pip_env_act} && pip install --no-cache-dir ./dist/flair-${VERSION}-py3-none-any.whl
	${pip_env_act} && ${MAKE} -C test test use_installed_flair=yes

# testpypy
publish-testpypi: build
	poetry publish -r testpypi --build

test-testpypi:
	${pip_env_setup}
	${pip_env_act} && pip install --no-cache-dir  --index-url=${testpypi_url} flair-brookslab==${version}
	${pip_env_act} && ${MAKE} -C test test use_installed_flair=yes


# pypy
publish-pypi: build
	poetry publish --build

test-pypi:
	${pip_env_setup}
	${pip_env_act} && pip install --no-cache-dir  --index-url=${pypi_url} flair-brookslab==${version}
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


