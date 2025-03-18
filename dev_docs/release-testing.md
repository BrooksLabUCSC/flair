# Pre-Release testing

These are tests to do before making the release.

## conda dev environment testing

1. delete and recreate the conda flair-dev environment
   ```
   conda deactivate  # if you are in a flair environment
   conda env remove --name flair-dev --yes
   conda env create --name flair-dev -f misc/flair_dev_conda_env.yaml --yes
   conda activate flair-dev
   make clean
   pip install -e .[dev]
   ```
2. run tests
   ```
   make -O -j 64 test
   ```
3. run deprecated diff expression tests 
   ```
   conda env update --name flair-dev --file misc/flair_diffexp_conda_env.yaml
   pip install -e .[diffexp]
   make -O -j 64 test-diffexp
   ```
