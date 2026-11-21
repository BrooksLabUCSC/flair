.. _installing-label:

Installing FLAIR
================

Installing with BioConda
~~~~~~~~~~~~~~~~~~~~~~~~

BioConda installs FLAIR and all of its dependencies:

.. code:: sh

   conda create -n flair -c conda-forge -c bioconda flair
   conda activate flair
   flair [transcriptome/quantify/...]

On Apple Silicon Mac systems (ARM64 CPUs: M1, M2, ...) you must use:

.. code:: sh

   CONDA_SUBDIR=osx-64 conda create -n flair
   conda activate flair
   conda config --env --set subdir osx-64
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda install flair

mamba fails to install FLAIR on Mac ARM64.

The BioConda ``flair`` package does not include the ``R`` packages that
``flair diffexp`` and ``flair diffsplice`` require.  Add them with:

.. code:: sh

   conda install -n flair -c conda-forge -c bioconda \
       r-argparse r-ggplot2 r-qqman r-lazyeval r-data.table \
       bioconductor-deseq2 bioconductor-drimseq bioconductor-stager bioconductor-apeglm

Installing a GitHub release with Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BioConda releases are often delayed by its manual review.  To create a Conda
environment from a release in GitHub, use:

.. code:: sh

   conda env create -n flair --file https://github.com/BrooksLabUCSC/flair/releases/latest/download/flair_conda_env.yaml
   conda activate flair

This environment includes the ``R`` packages that ``flair diffexp`` and
``flair diffsplice`` require.
For an older release, put its tag in place of ``latest/download``, as
``download/v3.0.1``.

Running FLAIR with Docker
~~~~~~~~~~~~~~~~~~~~~~~~~

FLAIR releases are also Docker images, which include the R dependencies of
``flair diffexp`` and ``flair diffsplice``.  Obtain and run one with:

.. code:: sh

    docker pull brookslab/flair:latest
    docker run -v [your_path_to_data]:/data brookslab/flair:latest flair [transcriptome/quantify/...]

``latest`` is the newest release; an older one is named by its version, as
``brookslab/flair:3.0.1``.

.. _pip-install-label:

Installing with pip
~~~~~~~~~~~~~~~~~~~

FLAIR runs under python 3.13 and 3.14, which is what ``requires-python`` in
``pyproject.toml`` spans; pip rejects the install under any other version.

These programs must be installed first, with Conda, with a system package manager
or by downloading and compiling them:

- `bedtools <https://github.com/arq5x/bedtools2/>`_
- `samtools <https://github.com/samtools/samtools/releases>`_
- `minimap2 <https://github.com/lh3/minimap2>`_
- bedToBigBed, from the `UCSC tools <https://hgdownload.soe.ucsc.edu/admin/exe/>`_,
  which writes the bigBed of ``flair transcriptome``, ``flair combine`` and
  ``flair fusion``
- `intron-prospector <https://github.com/diekhans/intron-prospector>`_, for flair
  fusion and for junctions from long reads

The conda environment also has `longshot <https://github.com/pjedge/longshot>`_,
which calls the variants of FLAIR's tests.  FLAIR does not run it itself.

.. code:: sh

   pip install flair-brookslab

or to install the current code from github:

.. code:: sh

   pip install git+https://github.com/BrooksLabUCSC/flair.git

An installed FLAIR must be uninstalled before the current code is installed over
it, since the version number in the tree does not change between releases.

.. code:: sh

   pip uninstall flair-brookslab
   pip install git+https://github.com/BrooksLabUCSC/flair.git


``flair diffexp`` and ``flair diffsplice`` require ``R`` and these ``R`` packages.
Not all of them have an Apple Silicon build; on Apple Silicon, install them into an
``osx-64`` environment, as above.

- `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`__
- `ggplot2 <https://ggplot2.tidyverse.org>`__
- `qqman <https://cran.r-project.org/web/packages/qqman/index.html>`__
- `DRIMSeq <https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html>`__
- `stageR <https://bioconductor.org/packages/release/bioc/html/stageR.html>`__


.. _whatshap-install-label:

WhatsHap
~~~~~~~~

`WhatsHap <https://whatshap.readthedocs.io>`_, or another phasing tool, prepares
the phased VCF and the haplotagged reads that ``flair isoalleles`` takes; see
:ref:`isoalleles-label`.  FLAIR does not run it itself.

It is in neither the conda environment files nor the docker image, since FLAIR does
not run it and bioconda builds it through python 3.13, which would pin the whole
environment to that version.  Install it into an environment of its own, so that it
does not constrain FLAIR's:

.. code:: sh

   conda create -n whatshap -c bioconda -c conda-forge 'whatshap>=2.8'
   conda activate whatshap

or with pip, which has wheels through python 3.13 and compiles from source above
that:

.. code:: sh

   pip install 'whatshap>=2.8'


Running from a GitHub clone
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dependencies must be installed first, as for :ref:`pip-install-label`.

.. code:: sh

    git clone https://github.com/BrooksLabUCSC/flair.git
    cd flair
    ./bin/flair [transcriptome/quantify/...]
