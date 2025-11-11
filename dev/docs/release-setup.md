# Onetime account setup

Accounts need for a new person releasing FLAIR.

## PyPi user/computer setup 
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

## ReadTheDocs user setup

Create a ReadTheDocs account if you don't already have one and have an
existing FLAIR a
registered as an admin on the FLAIR
project check here: https://app.readthedocs.org/dashboard/

## DockerHub user setup

Create a DockerHub account if you don't have one and ask a current
lab admin to add you to the brookslab organization.

