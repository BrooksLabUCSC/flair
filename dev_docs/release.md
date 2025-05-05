# Releasing FLAIR

For each Flair release we provide the following:

- pip package (basic flair; python dependencies; no bedtools, samtools, minimap2)
- a bioconda package (uses pip; full install)
- docker image (uses pip; full install)

If this is your first time doing, be user and do the *Onetime account setup*.

If you need to make edits to this document after the release is tag and pushed,
don't worry, make a temporary copy to update and then commit after the release.
This file is instructions, not history.

It is important to do the below steps in order. It is suggested to copy this
markdown checklist to the GitHub release ticket and use it to track progress.

- [] 1. Update CHANGELOG.md
- [] 2. Update python dependencies
- [] 3. Update dependency versions in Conda and Docker
- [] 4. Create a clean conda flair-dev environment
- [] 5. Run pre-release tests
- [] 6. Check documentation
- [] 7. Verify and review ReadTheDocs
- [] 8. Update version numbers in all relevant files
- [] 9. Build distribution and test pipe install
- [] 10. Test Docker locally
- [] 11. Test PyPi with testpypi
- [] 12. Git commit, tag push
- [] 13. Create pip package and test
- [] 14. Make the release on github.com
- [] 15. Update the conda recipe and submit
- [] 16. Build the docker image and push it to dockerhub
- [] 17. Set the release as latest in the PyPi
- [] 18. Set the release latest in GitHub
- [] 19. Announce release


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

Create a ReadTheDocs account if you don't already have one and have an
existing FLAIR a
registered as an admin on the FLAIR
project check here: https://app.readthedocs.org/dashboard/

### DockerHub user setup

Create a DockerHub account if you don't have one and ask a current
lab admin to add you to the brookslab organization.

## 1. Update CHANGELOG.md and copyrights

Edit CHANGELOG.md outline high-level changes.  This will contain new features,
incompatibilities, and high-visibility bug fixes.  It doesn't need to contain
minor changes.  Review the commit logs with:

```
git log 2.0.0..HEAD --pretty=format:"%h %an: %s"
```

Copyrights are needed in:

- `LICENSE`
- `docs/source/conf.py`

Also update authors in `docs/source/conf.py`.

## 2. Update python dependencies
```
poetry update
poetry local
git commit -am 'updated dependencies'
```

## 3. Update dependency versions in Conda and Docker

Review dependency versions in misc/*_conda_env.yaml to see if they should be
updated.  Use `conda search <package>` to find versions.

The FLAIR version in `misc/flair_conda_env.yaml` is updated by bump-my-version

Edit `misc/Dockerfile` to have the same non-Python package dependencies as the conda
files.  Also ensure that the Docker Ubuntu versions is a current LTS release.
You can find the versions of packages matching the Ubuntu release at
https://packages.ubuntu.com/. 


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

Include conda dependencies deprecated expression diff support and tests.
This does not work on Apple ARM systems.
```
conda env update --name flair-dev --file misc/flair_diffexp_conda_env.yaml
pip install -e .[diffexp]
make -O -j 64 test-expdiff
```

## 6. Check documentation

```
make doc
```
Review `docs/build/html/index.html` in your web browser.
Push any changes to github master.

## 7. Verify and review ReadTheDocs

Commit and push the master branch so that it is rebuilt in ReadTheDocs as
the latest.  The version numbers will be changed at this point.

```
git commit -am 'Release work'
git push
```

Readthedocs occasionally changes their requirements and when that happens the
build of the FLAIR documentation may start failing. To make sure that it's up
and running, log in at readthedocs.org.

   https://readthedocs.org/projects/flair/builds/

and check that the latest build. Builds are started immediately after every
push to the master branch.

Review the documentation on ReadTheDocs to make sure you are happy.


## 8. Update version numbers in all relevant files

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
src/flair/__init__.py
```

Use `bump-my-version` to increment the version numbers:
```
bump-my-version bump <major|minor|patch>
```
Before you do this, make sure the current version is correctly listed in `.bumpversion.cfg`.

Check that version were updated:
```
git diff
```

**DO NOT COMMIT VERSION CHANGE YET**

## 9. Build distribution and test pipe install
```
make clean build
make -O -j 32 test-pip
```

## 10. Test Docker locally

Build without installing FLAIR:
```
docker build --build-arg FLAIR_INSTALL=no --network=host -t brookslab/flair:<version> misc
```
Note that is can be very slow or fail on network errors without `--network=host`.


Install FLAIR and test.
```
docker run --rm -it -v $(pwd):/mnt/flair --network=host brookslab/flair:<version> bash
% cd /mnt/flair
% pip install --break-system-packages .
% make clean
% make -O -j 64 test-installed
% make clean
% exit
```
## 11. Test PyPi with testpypi

To upload to testpypi and test

```
make publish-testpypi
make -O -j 32 test-testpypi
```

If an error occurs, delete can the testpypi release.


## 12. Git commit, tag push
```
   git status
   git commit -am "FLAIR <version> release"
   git tag -a v<version> -m "Release <version>"
   git push
   git push --tags
```

## 13. Create pip package, and test

You can only do this once per release, so be sure to doe the testpipy test
first.  PyPi does not allow submission of the same release number twice.
While you can delete a release, you can not reuse the file names. If something
goes wrong, re-release with a new patch version.

The pypi package name is `flair-brookslab`.

```
   make publish-pypi
   make -O -j 32 test-pypi
```

## 14. Make the release on github.com

     https://github.com/BrooksLabUCSC/flair/releases

Select Draft a new release (top right) and follow instructions

Copy CHANGELOG.md entry to release description.

Set these options:
- Set as the latest release 
- Create a discussion for this release 

## 15. Update the conda recipe and submit

Full details are here: https://bioconda.github.io/contributor/index.html

1. fork the bioconda recipes repo: https://github.com/bioconda/bioconda-recipes/fork
2. git clone that repo to your local computer
3. update the recipe `bioconda-recipes/recipes/recipes/flair/meta.yaml` with
   - new version number,
   - correct dependencies from pyproject.yaml, these must be explicitly listed in meta.yaml
     - NOTE: check numpy version against conda-forge
   - pypi URL and sha256 of the `.tar.gz` source file found at https://pypi.org/project/flair-brookslab/#files
4  create a local bioconda environment and test:
```
   cd ../bioconda-recipes/
   conda create -n bioconda-test -c conda-forge -c bioconda bioconda-utils pytorch --yes
   conda activate bioconda-test
   bioconda-utils build --mulled-test recipes config.yml --packages flair >&build.out
   conda deactivate
   conda env remove -n bioconda-test --yes
```

Inspect `build.out` for problems.  There maybe unimportant warning labels as
`ERR` that do not relate to the FLAIR and can be ignore.  If there are errors
in the build, add the option `--keep-old-work` to build command to save output
for inspection.

  
5. git commit; git push
6. submit a pull request via https://github.com/bioconda/bioconda-recipes/pulls
   This starts a testing process. Once all checks have passed and a green mark appears, 
   add this comment to the pull request:
	    @BiocondaBot please add label
   This should take care of the red 'Review Required' and 'Merging is Blocked' notifications

Note: THESE STEP MIGHT NOT BE NECESSARY

If you do not make any changes to flair's dependencies (scipy, pandas, etc) then
the biocondabot may detect the new release and update the conda package automatically. 
Simply wait a few days, then check the version at 

    https://bioconda.github.io/recipes/flair/README.html

## 16. Build the docker image and push it to dockerhub

Docker does allow you to resubmit the same version number, it will overwrite the image if you do.

Please submit the container both as brookslab/flair:<current version number> and brookslab/flair:latest
The use of 'latest' is heavily discouraged by Docker because it obscures the actual version. However, 
when people try pulling without a version number the docker pull command automatically looks for the 'latest' tag.
Dockerhub is smart enough to just add a second label to the same image, so submitting it twice does not
take a lot of time.

Dependencies in  `misc/DockerFile` should have already been updated.  The version
number should have been updated by `bump-my-version`.

From the ./misc directory, build the image:
```
docker build --network=host -t brookslab/flair:<version> misc
docker tag brookslab/flair:<version> brookslab/flair:latest
```

Test the Docker image

```
docker run --rm -it -v $(pwd):/mnt/flair --network=host brookslab/flair:<version> bash
% cd /mnt/flair
% make -O -j 64 test-installed
% exit
```

Push to Dockerhub:
```
docker login -u <username>
docker push brookslab/flair:<version>
docker push brookslab/flair:latest
```

The reason that the build takes long is that pysam doesn't have a fast installation method.

## 17. Set the release as the latest in PyPi

## 18. Set the release latest in GitHub
- Change from pre-release to latest release
- Create a discussion for this release 

## 19. Announce release

Mail an announcement, including `CHANGELOG.md` summary to
`flair-announce-group@ucsc.edu`.

🍷🍷
