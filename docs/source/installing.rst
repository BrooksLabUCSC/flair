.. _installing-label:

Installing Flair
================

BioConda Install
~~~~~~~~~~~~~~~~

The easiest way to install Flair and all of its dependencies is via BioConda:

.. code:: sh

   conda create -n flair -c conda-forge -c bioconda flair
   conda activate flair
   flair [align/correct/...]

On Apple Silicon Mac systems (ARM64 CPUs: M1, M2, ...) you must use:

.. code:: sh

   CONDA_SUBDIR=osx-64 conda create -n flair
   conda activate flair
   conda config --env --set subdir osx-64
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda install flair

Note that mamba currently fails to install FLAIR on Mac ARM64.

The BioConda ``flair`` package does not include the ``R`` packages that ``diffexp``
and ``diffsplice`` need.  To use those modules, add them with:

.. code:: sh

   conda install -n flair -c conda-forge -c bioconda \
       r-argparse r-ggplot2 r-qqman r-lazyeval r-data.table \
       bioconductor-deseq2 bioconductor-drimseq bioconductor-stager bioconductor-apeglm

Conda Install from GitHub Release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BioConda release are often delayed by the manual review process.  To create
a Conda environment from a release in GitHub, use:

.. code:: sh

   conda env create -n flair --file https://github.com/BrooksLabUCSC/flair/releases/download/v3.0.0+dev/flair_conda_env.yaml
   conda activate flair

This environment includes the ``R`` packages that ``diffexp`` and ``diffsplice`` need.

Running FLAIR with Docker
~~~~~~~~~~~~~~~~~~~~~~~~~

FLAIR releases are also available as a Docker images, which includes the R
dependencies to run ``diffexp`` and ``diffsplice``.  Obtain and run the
FLAIR Docker with:

.. code:: sh

    docker pull brookslab/flair:3.0.0+dev
    docker run -v [your_path_to_data]:/data brookslab/flair:3.0.0+dev flair [align/correct/...]

.. _pip-install-label:

Pip install
~~~~~~~~~~~

To install with pip, you must first install these packages. This can be done
with Conda, yes system package manager or by downloading and compiling,

- `bedtools <https://github.com/arq5x/bedtools2/>`_
- `samtools <https://github.com/samtools/samtools/releases>`_
- `minimap2 <https://github.com/lh3/minimap2>`_

.. code:: sh

   pip install flair-brookslab

or to install the current code from github:

.. code:: sh

   pip install git+https://github.com/BrooksLabUCSC/flair.git

If you already have FLAIR install, you will need to uninstall it first
to update the code, since the version number will not have changed in the tree.

.. code:: sh

   pip uninstall flair-brookslab
   pip install git+https://github.com/BrooksLabUCSC/flair.git


The``flair diffexp`` and ``flair diffsplice`` modules require ``R`` , along
with these ``R`` packages. Some of these do not work on Apple Silicon.

- `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`__
- `ggplot2 <https://ggplot2.tidyverse.org>`__
- `qqman <https://cran.r-project.org/web/packages/qqman/index.html>`__
- `DRIMSeq <http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html>`__
- `stageR <http://bioconductor.org/packages/release/bioc/html/stageR.html>`__


Running from GitHub clone
~~~~~~~~~~~~~~~~~~~~~~~~~

The dependence must first be installed as for :ref:`pip-install-label`.

.. code:: sh

    git clone https://github.com/BrooksLabUCSC/flair.git
    cd flair
    ./bin/flair [align/correct/...]
