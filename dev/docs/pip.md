# Updating FLAIR'S pip package:

## initial setup

this section only needs to be done once

- make environment with wheel and build (only need to do this once)
    python3 -m venv mytestvenv
    source mytestvenv/bin/activate
    pip install --upgrade pip
    pip install wheel build

- make an account on pypi and ask the current owner (jeltje.van.baren@gmail.com as of this writing) 
  to be added to the flair-brookslab project

- create an API token (needed for safe upload), see https://pypi.org/help/#apitoken
- store the token in your local ~/.pypirc
- That file should look like this:
    [distutils]
    index-servers = pypi
    [pypi]
    repository = https://upload.pypi.org/legacy/
    username = __token__
    password = pypi-wa(...extremely long string...)XgkH

## FLAIR release

- update version where relevant
  - 1.6.3 to 1.6.4: if it's a bug fix
  - 1.6.3 to 1.7.0 if it's a bigger change
  - 1.6.3 to 2.0.0 if adding major new functionality and incompatible changes
    major revisions are also coordinated with papers
- files to update with bump2version
  - ./misc/Dockerfile
  - ./setup.cfg
  - ./src/flair/flair.py

   bump2version --allow-dirty --verbose patch flair.py misc/Dockerfile setup.cfg
   git add setup.cfg src/flair/flair.py misc/Dockerfile
   git commit -m "setting up release 1.6.4"
   # DO NOT PUSH YET

- create package (from main flair directory)
   rm -rf dist/*
   python3 -m build
   python3 -m twine upload dist/*


