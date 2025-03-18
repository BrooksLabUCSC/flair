root = .
include ${root}/defs.mk

# If you want to use the install flair rather than the one in this
# tree use:
#   make test use_installed_flair=yes


default:

.PHONY: test
test:
	cd test && ${MAKE} test use_installed_flair=${use_installed_flair}

test-diffexp:
	cd test && ${MAKE} test-diffexp use_installed_flair=${use_installed_flair}

clean:
	rm -rf build src/flair_brookslab.egg-info/ ls src/flair/__pycache__/
	cd test && ${MAKE} clean

real-clean: clean
	cd test && ${MAKE} real-clean
