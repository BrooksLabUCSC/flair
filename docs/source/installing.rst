.. _installing-label:

Installing Flair
================

BioConda Install
~~~~~~~~~~~~~~~~

The easiest way to install Flair and all of its dependencies is via BioConda:

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

   conda env create -n flair --file https://github.com/BrooksLabUCSC/flair/releases/latest/download/flair_conda_env.yaml
   conda activate flair

This environment includes the ``R`` packages that ``diffexp`` and ``diffsplice`` need.
For an older release, put its tag in place of ``latest/download``, as
``download/v3.0.1``.

Running FLAIR with Docker
~~~~~~~~~~~~~~~~~~~~~~~~~

FLAIR releases are also available as a Docker images, which includes the R
dependencies to run ``diffexp`` and ``diffsplice``.  Obtain and run the
FLAIR Docker with:

.. code:: sh

    docker pull brookslab/flair:latest
    docker run -v [your_path_to_data]:/data brookslab/flair:latest flair [transcriptome/quantify/...]

``latest`` is the newest release; an older one is named by its version, as
``brookslab/flair:3.0.1``.

.. _pip-install-label:

Pip install
~~~~~~~~~~~

To install with pip, you must first install these packages. This can be done
with Conda, your system package manager or by downloading and compiling,

- `bedtools <https://github.com/arq5x/bedtools2/>`_
- `samtools <https://github.com/samtools/samtools/releases>`_
- `minimap2 <https://github.com/lh3/minimap2>`_
- bedToBigBed, from the `UCSC tools <https://hgdownload.soe.ucsc.edu/admin/exe/>`_,
  which writes the bigBed of flair transcriptome, combine and fusion
- `intron-prospector <https://github.com/diekhans/intron-prospector>`_, for flair
  fusion and for junctions from long reads

The conda environment also has `longshot <https://github.com/pjedge/longshot>`_,
which calls the variants of FLAIR's tests.  FLAIR doesn't run it itself.

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
- `DRIMSeq <https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html>`__
- `stageR <https://bioconductor.org/packages/release/bioc/html/stageR.html>`__


.. _whatshap-install-label:

WhatsHap
~~~~~~~~

`WhatsHap <https://whatshap.readthedocs.io>`_, or another phasing tool, prepares
the phased VCF and the haplotagged reads that ``flair isoalleles`` takes; see
:ref:`isoalleles-label`.  FLAIR doesn't run it itself.

It isn't in the conda environment, as bioconda builds it for python 3.12 alone,
which would hold the whole environment at that version.  The docker image has it.
Install it when you need it, into an environment of its own so that it doesn't
constrain FLAIR's:

.. code:: sh

   conda create -n whatshap -c bioconda -c conda-forge 'whatshap>=2.8'
   conda activate whatshap

or with pip, which has wheels through python 3.13 and compiles from source above
that:

.. code:: sh

   pip install 'whatshap>=2.8'


Running from GitHub clone
~~~~~~~~~~~~~~~~~~~~~~~~~

The dependence must first be installed as for :ref:`pip-install-label`.

.. code:: sh

    git clone https://github.com/BrooksLabUCSC/flair.git
    cd flair
    ./bin/flair [transcriptome/quantify/...]
