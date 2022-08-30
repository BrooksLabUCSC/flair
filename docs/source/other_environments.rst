Other environments
==================

Docker 
------

Flair is available as a `docker <https://www.docker.com/>`_ image and can be run like so:

.. code:: sh

    docker pull brookslab/flair:latest
    docker run -w /usr/data -v [your_path_to_data]:/usr/data brookslab/flair:latest flair [align/correct/...]


Conda Environment 
-----------------

Basic Flair (align, correct, collapse) can be installed and run in a `conda <https://docs.conda.io/en/latest/>`_ environment using the ``.yaml`` file provided in ``misc/flair_basic_conda_env.yaml``:  

.. code:: sh

    conda env create -f misc/flair_basic_conda_env.yaml
    conda activate flair_basic_conda
    flair [align/correct/...]

**Important note**: While this method installs all of Flair, it does not install all dependencies for full Flair. Flair diffExp and diffSplice rely on R packages that we have not been able to install correctly using the conda method.

If R is installed on your system, add the necessary libraries like so:

.. code:: sh

    R -e 'update.packages(ask=FALSE, repos = "http://cran.us.r-project.org")'
    R -e 'install.packages(c("devtools", "BiocManager", "ggplot2", "qqman", "lazyeval"), repos = "http://cran.us.r-project.org")'
    R -e 'requireNamespace("BiocManager"); BiocManager::install(c("DRIMSeq", "stageR", "DESeq2", "apeglm"))'

You will also need to install one additional conda package:

.. code:: sh

    conda install -c bioconda pandas

And one pip package:

.. code:: sh
   pip install rpy2

