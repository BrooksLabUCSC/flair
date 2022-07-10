Other environments
==================

Docker 
------

If the user wishes to run FLAIR using
`docker <https://docs.docker.com/>`__ instead of cloning this
repository, the following commands can be used:
``docker pull quay.io/brookslab/flair``
``docker run -w /usr/data -v [your_path_to_data]:/usr/data  -t -d [image_id]``
``docker exec [container_id] python3 /usr/local/flair/flair.py [your_command]``

Conda Environment 
-----------------

Users can run FLAIR within the conda environment provided in
``misc/flair_conda_env.yaml``. FLAIR should run smoothly in this
environment. Refer to the conda docs for how to `create an environment
from an environment.yml
file <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`__.

