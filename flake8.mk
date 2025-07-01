
# relative to src/flair
FLAKE8_SRC = \
	__init__.py \
	flair_cli.py

FLAKE8_CHECK = ${FLAKE8_SRC:%=src/flair/%}
