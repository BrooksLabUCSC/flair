
# relative to src/flair
FLAKE8_SRC = \
	__init__.py

FLAKE8_CHECK = ${FLAKE8_SRC:%=src/flair/%}
