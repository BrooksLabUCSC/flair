Installing Flair
================

Flair consists of six modules. The first three are ``align``,
``correct``, and ``collapse``. They are the most used, so we
will refer to them here as basic flair.

The other three modules are ``quantify``, ``diffExp``, and 
``diffSplice``. Together with basic flair these are called full flair.
These three additional modules have more dependencies than basic flair
so if you don't need them, you will not need the modules listed under 5.

There are three ways to install flair:

* ``pip install flair-brookslab`` will install basic flair and all necessary python modules (see below)
* Download `the latest release <https://github.com/BrooksLabUCSC/flair/releases>`_
* Use git to check out `the current flair repository <https://github.com/BrooksLabUCSC/flair.git>`_

Requirements
~~~~~~~~~~~~

1. `Bedtools <https://github.com/arq5x/bedtools2/>`_
2. `samtools <https://github.com/samtools/samtools/releases>`_
3. `minimap2 <https://github.com/lh3/minimap2>`_

If you do not use ``pip install``, you will also need:

4. python v3.6+ and python modules: 
   * numpy=1.9.*
   * tqdm
   * ncls
   * pybedtools
   * mappy
   * pysam=v0.8.4+
5. full flair additional python modules:
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

For ``conda`` users: ``yaml`` environment files can be found in flair's 
``/misc`` directory.

Pip install
~~~~~~~~~~~

``pip install flair-brookslab`` will put the latest flair release in your ``$PATH``, as well
as the helper scripts discussed in this manual. It also installs all python modules
needed to run basic flair. If you want to use full flair, install the packages
listed under point 5 in the list above.


Download the latest release
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to the Flair `release page <https://github.com/BrooksLabUCSC/flair/releases>`_
and select one of the source code files under Assets. Exctract the file and navigate
to the resulting `flair` directory. Add flair and the helper scripts to your ``$PATH``
for instance (in Linux) with ``export PATH=$(pwd)/bin:$PATH``. 

Make sure to (`pip`) install the python modules listed above. If you have conda, you can
create a (full or basic) flair environment using

``conda env create -f misc/flair_conda_env.yaml``

``conda env create -f misc/flair_basic_conda_env.yaml``


Download the latest code (not preferred)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check out `the current flair repository <https://github.com/BrooksLabUCSC/flair.git>`_
from github. Please be aware that while this may have the latest bug fixes, it's quite
possible that new bugs were introduced. This method is only useful if you have 
`reported a problem <https://github.com/BrooksLabUCSC/flair/issues>`_ and a Flair developer
lets you know it has been fixed.

Once you have cloned the repository, navigate to the `/flair` directory. Follow the
steps as described under Download the latest release.

