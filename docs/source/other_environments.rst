Other environments
==================

The easiest way to install Flair and all of its dependencies is via conda:

.. code:: sh

   conda create -n flair -c conda-forge -c bioconda flair
   conda activate flair
   flair [align/correct/...]

It is also possible to get the full Flair setup as a docker image:

.. code:: sh

    docker pull brookslab/flair:latest
    docker run -w /usr/data -v [your_path_to_data]:/usr/data brookslab/flair:latest flair [align/correct/...]

You can build a Conda environment from the source tree with
.. code:: sh

    git clone git@github.com:BrooksLabUCSC/flair.git
    cd flair
    git checkout v<version-tag>
    conda env create -n flair -f misc/flair_conda_env.yaml
    conda activate flair
    pip install .
    

Other methods (not recommended)
-------------------------------

Flair consists of six modules. The first three are ``align``,
``correct``, and ``collapse``. They are the most used, so we
will refer to them here as basic Flair.

The other three modules are ``quantify``, ``flair_diffExp``, and 
``flair_diffSplice``. Together with basic Flair these are called full Flair.
These three additional modules have more dependencies than basic Flair
so if you don't need them, you will not need the modules listed under 5.

There are other ways to install Flair:

* ``pip install flair-brookslab`` will install basic Flair and all necessary python modules (see below)
* Download `the latest release <https://github.com/BrooksLabUCSC/flair/releases>`_
* Use git to check out `the current flair repository <https://github.com/BrooksLabUCSC/flair.git>`_


Requirements
~~~~~~~~~~~~

1. `Bedtools <https://github.com/arq5x/bedtools2/>`_
2. `samtools <https://github.com/samtools/samtools/releases>`_
3. `minimap2 <https://github.com/lh3/minimap2>`_

If you do not use ``pip install`` or ``conda env create```, you will also need:

4. python v3.12+ and python modules: 

   * numpy=1.12.*
   * ncls
   * pybedtools
   * mappy
   * pysam=v0.8.4+

5. full Flair additional python modules:

  - Cython
  - pandas
  - rpy2=2.9.*
  - R
  - r-ggplot2=2.2.1
  - r-qqman
  - bioconductor-deseq2
  - bioconductor-drimseq
  - bioconductor-stager
  - matplotlib
  - seaborn


Pip install
~~~~~~~~~~~

``pip install flair-brookslab`` will put the latest Flair release in your ``$PATH``, as well
as the helper scripts discussed in this manual. It also installs all python modules
needed to run basic Flair. If you want to use full Flair, install the packages
listed under point 5 in the list above.


Download the latest release
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to the Flair `release page <https://github.com/BrooksLabUCSC/flair/releases>`_
and select one of the source code files under Assets. Exctract the file and navigate
to the resulting `flair` directory. Add Flair and the helper scripts to your ``$PATH``
for instance (in Linux) with ``export PATH=$(pwd)/bin:$PATH``. 

Make sure to (`pip`) install the python modules listed above. If you have conda, you can
create a basic Flair environment using

``conda env create -f misc/flair_basic_conda_env.yaml``


Download the latest code
~~~~~~~~~~~~~~~~~~~~~~~~

Check out `the current Flair repository <https://github.com/BrooksLabUCSC/flair.git>`_
from github. Please be aware that while this may have the latest bug fixes, it's quite
possible that new bugs were introduced. This method is only useful if you have 
`reported a problem <https://github.com/BrooksLabUCSC/flair/issues>`_ and a Flair developer
lets you know it has been fixed.

Once you have cloned the repository, navigate to the `/flair` directory. Follow the
steps as described under Download the latest release.

