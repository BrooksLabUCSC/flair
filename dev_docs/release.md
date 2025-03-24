# Releasing FLAIR

For each Flair release we provide the following:

- pip package (basic flair; python dependencies; no bedtools, samtools, mapman)
- conda package (uses pip; full install)
- docker image (uses pip; full install)

It is important to do the below steps in order:

## Onetime account setup

### PyPi user/computer setup 
You will need a PyPi accounts for, both https://pypi.org/ and https://testpypi.pypi.org
and be added to the `flair-brookslab` project by a current owner.

-  Create an API token, see https://pypi.org/help/#apitoken
-  Store the token in your local `~/.pypirc`:
```
[distutils]
index-servers = 
    pypi
    testpypi
    
[pypi]
repository = https://upload.pypi.org/legacy/
username = __token__
password = <your-api-token>

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = <your-test-api-token>
```

Poetry also needs this information:
```
    poetry config pypi-token.pypi <your-api-token>>
    poetry config repositories.testpypi https://test.pypi.org/legacy/
    poetry config pypi-token.testpypi <your-test-api-token>>
```

### ReadTheDocs user setup

You must be have a ReadTheDocs account and registered as an admin on the flair
project check here: https://app.readthedocs.org/dashboard/

## 1. Update CHANGELOG.md

Edit CHANGELOG.md outline high-level changes.  This will contain new features,
incompatibilities, and high-visibility bug fixes.  It doesn't need to contain
minor changes.  Review the commit logs with

## 2. Update dependence
```
   poetry update
   poetry local
   git commit -am 'updated dependencies'
```

## 3. Update dependency versions in conda environments

Review dependency versions in misc/*_conda_env.yaml to see if they should be
updated.  Use `conda search <package>` to find versions.

The FLAIR version in `misc/flair_conda_env.yaml` is updated by bump-my-version

## 4. Create a clean conda flair-dev environment
```
   conda deactivate  # if you are in a flair environment
   conda env remove --name flair-dev --yes
   conda env create --name flair-dev -f misc/flair_dev_conda_env.yaml --yes
   conda activate flair-dev
   make clean
   pip install -e .[dev]
```
   
## 5. Run pre-release tests

```
   make -O -j 64 test
```
Repeat this on Apple ARM (M1, M2, ...) processor systems.
Push any pending changes to master HEAD and use git.

Include conda dependencies deprecated diff expression support and tests.
This does not work on Apple ARM systems.
```
   conda env update --name flair-dev --file misc/flair_diffexp_conda_env.yaml
   pip install -e .[diffexp]
   make -O -j 64 test-diffexp
```

# 1. Check documentation

```
    make doc
```
Review `docs/build/html/` in your web browser.
Push any changes to github master.

## 2.5 Verify and review ReadTheDocs

Readthedocs occasionally changes their requirements and when that happens the
build of the FLAIR documentation may start failing. To make sure that it's up
and running, log in at readthedocs.org.

   https://readthedocs.org/projects/flair/builds/

and check that the latest build. Builds are started immediately after every
push to the master branch.

Review the documentation on ReadTheDocs to make sure you are happy


## 7. Update version numbers in all relevant files

Make sure tree is committed at this point, it makes undoing problems
with setting version easier.

Version numbers should be increased following the major/minor/patch logic:
- patch: 1.6.3 to 1.6.4: if it's a bug fix
- minor: 1.6.3 to 1.7.0 if new functionality, minor incompatibles
- major: 1.6.3 to 2.0.0 if adding major new functionality and incompatible changes
   major revisions are also coordinated with papers

The following files contain version numbers:
```
  pyproject.toml
  defs.mk
  misc/Dockerfile
  misc/flair_conda_env.yaml
  src/flair/flair_cli.py
```

Use `bump-my-version` to increment the version numbers:
```
  bump-my-version bump major|minor|patch
```
Before you do this, make sure the current version is correctly listed in `.bumpversion.cfg`.

Check that version were updated:
```
   git diff
```

DO NOT COMMIT VERSION CHANGE YET.

## 3. Build distribution and test
```
   make build
   make -O -j 32 test-pip
```

## 4. Test PyPi with testpypi

To upload to testpypi and test

```
   make publish-testpypi
   make -O -j 32 test-testpypi
```

If an error occurs, you can not update testpypi, instead you must create a dev
version, such as `2.0.0.dev0` with:

```
   bump-my-version bump dev
```

Then repeat above publish testing.  Once all is working you can reset to the new
release version without `.dev0` the with a command in the form

```
   bump-my-version bump --current-version 2.0.0.dev0 --new-version 2.0.0
```

## 2. Git commit, tag push
```
   git status
   git commit -am "setting up release <current release version>"
   git tag -a v<version> -m "Release v<version>"
   git push
   git push --tags
```

## 4. Create the pip package and upload it

You can only do this once per release, so be sure to doe the testpipy test
first.  PyPi does not allow submission of the same release number twice.  If
something goes wrong, re-release with a new patch version.

The pypi package name is `flair-brookslab`.

```
   make publish-pypi
   make -O -j 32 test-pypi
```

## 1. Make the release on github.com

     https://github.com/BrooksLabUCSC/flair/releases

Select Draft a new release (top right) and follow instructions

Copy CHANGELOG.md entries to release description.


## 5. Update the conda recipe and submit

### THIS STEP MIGHT NOT BE NECESSARY ####
If you do not make any changes to flair's dependencies (scipy, pandas, etc) then
the biocondabot may detect the new release and update the conda package automatically. 
Simply wait a few days, then check the version at https://anaconda.org/bioconda/flair
#########################################

Full details are here: https://bioconda.github.io/contributor/index.html

1. Fork the bioconda recipes repo: https://github.com/bioconda/bioconda-recipes/fork
2. git clone that directory to your local computer
3 (optional, for when you make dependency changes). create a bioconda environment for testing:
      conda create -n bioconda -c conda-forge -c bioconda bioconda-utils pytorch
      conda activate bioconda
4. update the recipe in recipes/flair/meta.yaml with the new version number
   and the pypi url and md5, found at https://pypi.org/project/flair-brookslab/(current version)/#files
5. git commit, git push
6. submit a pull request via https://github.com/bioconda/bioconda-recipes/pulls
	This starts a testing process. Once all checks have passed and a green mark appears, 
	add this comment to the pull request:
	    @BiocondaBot please add label
	This should take care of the red 'Review Required' and 'Merging is Blocked' notifications
7. Delete your fork.

####### 6. Build the docker image locally using the updated Dockerfile and push it to dockerhub ######

Docker does allow you to resubmit the same version number, it will overwrite the image if you do.

Please submit the container both as brookslab/flair:<current version number> and brookslab/flair:latest
The use of 'latest' is heavily discouraged by Docker because it obscures the actual version. However, 
when people try pulling without a version number the docker pull command automatically looks for the 'latest' tag.
Dockerhub is smart enough to just add a second label to the same image, so submitting it twice does not
take a lot of time.

    ##### setup section #######################################
    
    Ask to be added to the dockergroup on your system if you aren't already.
    
    Create an account on hub.docker.com
    Ask (Jeltje or Cameron) to be added as a member to the brookslab organization
    
    ##### end of setup section ################################

From the ./misc directory, run
    docker build -t brookslab/flair:<current version number> .
    docker tag brookslab/flair:<current version number> brookslab/flair:latest
    docker push brookslab/flair:<current version number>
    docker push brookslab/flair:latest
The reason that the build takes long is that pysam doesn't have a fast installation method.
