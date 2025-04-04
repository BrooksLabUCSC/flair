# FLAIR Development

Python dependencies and package creation are managed with `poetry`.
A conda environment is use to provide the non-python environment.

## Creating the conda enviroment
Create a new conda environment for developing FLAIR from the top of
the tree and add the current packages with poetry.

```
conda env create --name flair-dev -f misc/flair_dev_conda_env.yaml
conda activate flair-dev
make clean
pip install -e .[dev]
```

To run the deprecated `diffExp` and `diffSplice` tests:
```
conda env update --name flair-dev --file misc/flair_diffexp_conda_env.yaml
pip install -e .[diffexp]
```

If you get warning like:
```
warning  libmamba Problem type not implemented SOLVER_RULE_STRICT_REPO_PRIORITY
```
and then it fails with conflicts, make sure your conda is current with:

```
conda update -n base -c conda-forge conda
conda install mamba -c conda-forge
```
and then change to flexible channel_priority:
```
 conda config --set channel_priority flexible
```

## Running FLAIR program during development.

FLAIR being a flair-dev conda environment or a virtual environment 
in pip editable mode using `pip install -e`.  However, the programs
it `bin/` will be copied rather than linked, so they are not
editable without a reinstall.  This will be fixed in the future.

Alternatively, one can add FLAIR programs to the PATH by
doing the following while in the top level directory of 
the FLAIR cloned repository:

```
export PATH=$(pwd)/bin:${PATH}
export PYTHONPATH=$(pwd)/src:${PYTHONPATH}
```

## Managing dependencies

Poetry is used to manage dependencies.

Poetry Cheat sheet:
* add dependency: `poetry add pysam`
* add a development dependency: `poetry add --dev flake8`
* show dependencies: `poetry show`
* check `pyproject.toml`: `poetry check`
* update dependencies to their latest version: `poetry update`
* verify update: `poetry show --latest`
* update `pyproject.toml` from `poetry.lock` file: `poetry sync`
* install dependencies in virtual : `poetry install`

Note: you need to commit after make changes to packages or updating.


## Testing:

To run the tests using the source in the tree:
```
make -O -j 64 test
```

To run the tests using the version of FLAIR found on PATH:
```
make -O -j 64 test use_installed_flair=yes
```

To run the deprecated diffExp and diffSplice tests:
```
make -O -j 64 test-diffexp
```

# Releasing flair

The following documents cover releasing FLAIR:

* [Releasing FLAIR](release.md)
