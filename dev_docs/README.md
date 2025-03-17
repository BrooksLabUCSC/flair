# FLAIR Development

Python dependencies and package creation are managed with `poetry`.
A conda environment is use to provide the non-python environment,
including R and dependencies. 


## Creating the conda enviroment
Create a new conda environment for developing FLAIR from the top of
the tree and add the current packages with poetry.
This is the recommended approach, as it can be difficult
to get the right R pieces installed and the `rpy2` python package built:
```
conda env create --name flair-dev -f misc/flair_dev_conda_env.yaml
conda activate flair-dev
```

To run the deprecated `diffExp` and `diffSplice` tests:
```
conda install --name flair-dev --file misc/flair_diffexp_conda_env.yaml
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


Note you need to commit after make changes to packages or updating.


## Testing:

To run the tests using the source in the tree:

```
cd test
make -O -j 64 test
```

To run the tests using the version of FLAIR found on PATH:

```
cd test
make -O -j 64 test use_installed_flair=yes
```

To run the deprecated diffExp and diffSplice tests:
```
cd test
make -O -j 64 test-diffexp
```
