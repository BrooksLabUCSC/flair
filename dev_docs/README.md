# FLAIR Development

Python dependencies and package creation are managed with `poetry`.
A conda environment is use to provide the non-python environment.

## Creating the conda environment and installing FLAIR

**This is the recommend method FLAIR developers.**

Create a new conda environment for developing FLAIR from the top of
the tree and add the current packages with poetry.

```
conda env create --name flair-dev -f misc/flair_dev_conda_env.yaml --yes
conda activate flair-dev
make clean
pip install -e .[dev]
```

If you want to install only the dependencies without an FLAIR, use:
```
poetry install  --no-root
```
instead of `pip install`.

To run `diffexp` and `diffsplice`, install the dependencies with:
```
conda env update --name flair-dev --file misc/flair_diffexp_conda_env.yaml
pip install -e .[diffexp]
```
Note, the `pip` edit install is needed, as this environment install the
released version of FLAIR from PyPi.

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

## Running FLAIR program during development

The preferred way to run FLAIR is using flair-dev conda environment as
described above, which will have it in the PATH.  By using the `-e`
(`--editable`) option to `pip`, edits to the source tree will be reflected in
the conda environment without a need for another install.

## Managing dependencies

Poetry is used to manage dependencies and build PyPi packages.
The poetry lock facility to manage indirect dependencies is not
used due to incompatibility with conda.  Only top-level dependencies
are managed for FLAIR.

Protry
* add dependency to `pyproject.toml`: `poetry add pysam`
* add a development dependency to `pyproject.toml`: `poetry add --dev flake8`
* show dependencies: `poetry show`
* show current and latest version of dependencies: `poetry show --latest`
* show outdated dependencies: `poetry show --outdated`
* update a dependency: `poetry add pysam@<version>`
* update a dependency to the latest: `poetry add pysam@latest`
* check `pyproject.toml`: `poetry check`
* updated the lock file: `poetry update`

Normally `pip` is used to install rather than `poetry install`.  Updating the
lock file servers only serves as a reference due to Conda issues.

The `poetry sync` to install exact version of dependencies due to Conda has
incorrectly creating a `direct_url.json` in the `site-packages/* dist-info`
directories. Removing these files still results in failures for unknown
reasons, so  `poetry sync` is not used.

## Testing:

To run the tests using the source in the tree:
```
make -O -j 64 test
```

To run the tests using the version of FLAIR found on PATH:
```
make -O -j 64 test use_installed_flair=yes
```

To run the diffExp and diffSplice tests:
```
make -O -j 64 test-expdiff
```

To run all tests:
```
make -O -j 64 test-all
```

# Releasing flair

The following documents cover releasing FLAIR:

* [Releasing FLAIR](release.md)
