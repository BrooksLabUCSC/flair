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
conda env create -f misc/flair_dev_conda_env.yaml
conda activate flair-dev
poetry install
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



Poetry Cheat sheet:
* install dependencies in virtual : `poetry install`
* add dependency: `poetry add pysam`
* add a development dependency: `poetry add --dev flake8`
* show dependencies: `poetry show`
* update dependencies to their latest version: `poetry update`
* verify update: `poetry show --latest`
* update `pyproject.toml` from `poetry.lock`: `poetry sync`


Note you need to commit after make changes to packages or updating.
