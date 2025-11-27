# Releasing FLAIR

For each Flair release we provide the following:

- pip package (basic flair; python dependencies; no bedtools, samtools, minimap2)
- a bioconda package (uses pip; full install)
- docker image (uses pip; full install)

If this is your first time doing, be user and do the
[*Onetime account setup*](release-setup.md).

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
- [] 9. Build distribution and test pip install
- [] 10. Test Docker locally
- [] 11. Test PyPi with testpypi
- [] 12. Git commit, tag push
- [] 13. Create pip package and test
- [] 14. Make the release on github.com
- [] 15. Update the BioConda recipe and submit
- [] 16. Build the docker image and push it to dockerhub
- [] 17. Set the release as latest in the PyPi
- [] 18. Set the release latest in GitHub
- [] 19. Announce release
- [] 20. Set version in tree to have +master suffix

## 1. Update CHANGELOG.md and copyrights

Edit CHANGELOG.md outline high-level changes.  This will contain new features,
incompatibilities, and high-visibility bug fixes.  It doesn't need to contain
minor changes.  Review the commit logs with:

```
git log v2.0.0..HEAD --pretty=format:"%h %an: %s"
```

Copyrights are needed in:

- `LICENSE`
- `docs/source/conf.py`

Also update authors in `docs/source/conf.py`.

## 2. Update python dependencies
```
poetry update
poetry lock
git commit -am 'updated dependencies'
```

## 3. Update dependency versions in Conda and Docker

Review dependency versions in misc/*_conda_env.yaml to see if they should be
updated.  Use `conda search <package>` to find versions.

The FLAIR version in `misc/*.yaml` files are updated by bump-my-version.

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
make -k -O -j 64 test-installed
```
Repeat this on Apple ARM (M1, M2, ...) processor systems.

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

If the bump rules doesn't work right, just use 1--new-version` to force the
version.

General use of `bump-my-version`:

```
# show possible bumps
bump-my-version show-bump

# Force to any version, ignoring all rules
bump-my-version bump --new-version 3.0.0b1

# Increment alpha build number
bump-my-version bump build      # → 3.0.0a3
bump-my-version bump build      # → 3.0.0a4

# Jump to beta on demand
bump-my-version bump release    # → 3.0.0b1

# Increment beta build number
bump-my-version bump build      # → 3.0.0b2
bump-my-version bump build      # → 3.0.0b3

# Release final version
bump-my-version bump release    # → 3.0.0

# Start next version
bump-my-version bump patch      # → 3.0.1a1
bump-my-version bump minor      # → 3.1.0a1
bump-my-version bump major      # → 4.0.0a1
```

**DO NOT COMMIT VERSION CHANGE YET**

## 9. Build distribution and test pip install
```
make clean build
make -k -O -j 32 test-pip
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
% make -k -O -j 64 test-installed
% make clean
% exit
```

## 11. Test PyPi with testpypi

To upload to testpypi and test

```
make publish-testpypi
make -k -O -j 32 test-testpypi
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
make -k -O -j 32 test-pypi
```

## 14. Make the release on github.com

     https://github.com/BrooksLabUCSC/flair/releases

Select Draft a new release (top right) and follow instructions

Copy CHANGELOG.md entry to release description.

Add the three Conda `misc/*.yaml` files as assets of the release
to support direct creation of Conda environments via the URLs.

Set these options:
- Set as the latest release 
- Create a discussion for this release 

## 15. Update the BioConda recipe and submit

Full details are here: https://bioconda.github.io/contributor/index.html

1. fork the bioconda recipes repo: https://github.com/bioconda/bioconda-recipes/fork
   if you already have a fork, sync with master on GitHub
2. git clone or pull that repo to your local computer
3. update the recipe `bioconda-recipes/recipes/flair/meta.yaml` with
   - new version number,
   - correct dependencies from pyproject.yaml, these must be explicitly listed in meta.yaml
     - NOTE: check numpy version against conda-forge
   - pypi URL and sha256 of the `.tar.gz` source file found at https://pypi.org/project/flair-brookslab/#files
4  create a local bioconda environment and test:
```
   cd ../bioconda-recipes/
   conda env remove -n bioconda-test --yes
   conda create -n bioconda-test -c conda-forge -c bioconda bioconda-utils pytorch --yes
   conda activate bioconda-test
   bioconda-utils build --mulled-test recipes config.yml --packages flair >&build.out
   conda deactivate
```

Inspect `build.out` for problems.  There maybe unimportant warning labels as
`ERR` that do not relate to the FLAIR and can be ignore.  If there are errors
in the build, add the option `--keep-old-work` to build command to save output
for inspection.

NOTE: If there are NameResolutionErrors in the file, this is due to needing to 
`--network=host` with Docker, however there is no way to pass it in.
Just give up and push and let the bot check it/

  
5. git commit -am 'FLAIR v<version> release'
6. git push
7. submit a pull request via https://github.com/bioconda/bioconda-recipes/pulls
   Title of pull request should be "Update FLAIR <version> release"
   This starts a testing process. Once all checks have passed and a green mark appears, 
   Once the tests have passed, add this comment to the pull request:
	    @BiocondaBot please add label
   This should take care of the red 'Review Required' and 'Merging is Blocked' notifications

Note: THESE STEP MIGHT NOT BE NECESSARY

If you do not make any changes to flair's dependencies (scipy, etc) then
the biocondabot might detect the new release and update the conda package automatically. 
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
docker build --network=host -t brookslab/flair:<version> misc >& build.log

docker tag brookslab/flair:<version> brookslab/flair:latest
```

Test the Docker image

```
docker run --rm -it -v $(pwd):/mnt/flair --network=host brookslab/flair:<version> bash
% cd /mnt/flair
% make -k -O -j 64 test-installed
% make clean
% exit
```

Push to Dockerhub:
```
docker login -u <username>
docker push brookslab/flair:<version>
docker push brookslab/flair:latest
```

The reason that the build takes long is that pysam doesn't have a fast installation method.

## 17. Check the release is the latest in PyPi
This happens automatically bases on version numbers.
An alpha or beta release will not be made the latest.

## 18. Set the release latest in GitHub
- Change from pre-release to latest release
- Create a discussion for this release 

## 19. Announce release

Mail an announcement, including `CHANGELOG.md` summary to
`flair-announce-group@ucsc.edu`.

## 20. Set version in tree to have +master suffix

Add `+master` to the current version to distinguish tree checkouts 
```
bump-my-version replace --current-version=2.2.0 --new-version=2.2.0+master
```

IMPORTANT: then manually edit `.bumpversion.toml` to have this version,
as `replace` doesn't change this.


```
git commit -am 'set tree local version'
```


🍷🍷
