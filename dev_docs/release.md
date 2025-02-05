# Releasing FLAIR

For each Flair release we provide the following:

- pip package (basic flair; python dependencies; no bedtools, samtools, mapman)
- conda package (uses pip; full install)
- docker image (uses pip; full install)

It is important to do these in the correct order, namely:
1. Update version numbers in all relevant files (with bumpversion)
2. Git commit; git push
3. Make the release on github.com
4. Create the pip package and upload it
5. Update the conda recipe and submit
6. Build the docker image locally using the updated Dockerfile and push it to dockerhub

Once you have done this once and are properly set up, these six steps can be done quickly.

# 1. Update version numbers in all relevant files ###############

Version numbers should be increased following the major/minor/patch logic:
  - 1.6.3 to 1.6.4: if it's a bug fix
  - 1.6.3 to 1.7.0 if new functionality
  - 1.6.3 to 2.0.0 if adding major new functionality and incompatible changes
    major revisions are also coordinated with papers

As of this writing, the following files contain version numbers:
  ./misc/Dockerfile
  ./setup.cfg  (these are the pip instructions)
  ./src/flair/flair.py
You can update these manually, or use bump2version (which can be pip installed) like so:
  bump2version --allow-dirty --verbose major|minor|patch src/flair/flair.py misc/Dockerfile setup.cfg
Before you do this, make sure the current version is correctly listed in .bumpversion.cfg

####### 2. Git commit; git push #######################################

   # make sure that the version number did actually change
   git diff 
   git add setup.cfg src/flair/flair.py misc/Dockerfile
   git commit -m "setting up release <current release version>"
   git push


####### 2.5 Check readthedocs #########################################

Readthedocs occasionally changes their requirements and when that happens the build of the Flair
documentation may start failing. To make sure that it's up and running, 
log in at readthedocs.org (you must be registered there as an admin on the flair project)
and go to https://readthedocs.org/projects/flair/builds/
and check that the latest build Passed.
(Builds are started immediately after every push to the master branch)


####### 3. Make the release on github.com #############################

https://github.com/BrooksLabUCSC/flair/releases
Select Draft a new release (top right) and follow instructions
Please do describe what's new in this release, both for our users and for your own future sanity.
Hint: use git log to see commit messages.

####### 4. Create the pip package and upload it #######################

You can only do this once per release, so be careful. Pypi does not allow submission of the same release number twice.

    ##### the section below only needs to be done once ########
    
    Make a virtual environment with wheel and build
       python3 -m venv mytestvenv
       source mytestvenv/bin/activate
       pip install --upgrade pip
       pip install wheel build
    
    Create an account on pypi and ask the current owner (jeltje.van.baren@gmail.com as of this writing) 
    to be added to the flair-brookslab project
    
    Create an API token (needed for safe upload), see https://pypi.org/help/#apitoken
    Store the token in your local ~/.pypirc
    That file should look like this ($ cat ~/.pypirc):
       [distutils]
       index-servers =
           pypi
       
       [pypi]
       repository = https://upload.pypi.org/legacy/
       username = __token__
       password = pypi-wa(...extremely long string...)XgkH
    
    ##### end of setup section ################################

Create the package (from the main Flair directory)
   rm -rf dist/*
   python3 -m build
   # and upload to pypi
   python3 -m twine upload dist/*

####### 5. Update the conda recipe and submit #########################

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
