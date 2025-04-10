Other environments
==================

The easiest way to install FLAIR and all of its dependencies is via conda:

.. code:: sh

   conda create -n flair -c conda-forge -c bioconda flair
   conda activate flair
   flair [align/correct/...]

It is also possible to get the full FLAIR setup as a docker image:

.. code:: sh

    docker pull brookslab/flair:latest
    docker run -w /usr/data -v [your_path_to_data]:/usr/data brookslab/flair:latest flair [align/correct/...]

You can build a Conda environment from GitHub, for example for release v2.1.0:

.. code:: sh

    conda env create -n flair -f https://raw.githubusercontent.com/BrooksLabUCSC/flair/refs/tags/v2.1.0/misc/flair_conda_env.yaml
    conda activate flair
    
You can also build a Conda environment from the source tree with
.. code:: sh

    git clone git@github.com:BrooksLabUCSC/flair.git
    cd flair
    git checkout v<version-tag>
    conda env create -n flair -f misc/flair_conda_env.yaml
    conda activate flair
    pip install .
    

Manual installation methods
---------------------------

Requires not installable with pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. `Bedtools <https://github.com/arq5x/bedtools2/>`_
2. `samtools <https://github.com/samtools/samtools/releases>`_
3. `minimap2 <https://github.com/lh3/minimap2>`_

Pip install
~~~~~~~~~~~

``pip install flair-brookslab`` will put the latest FLAIR release in your ``$PATH``, as well
as the helper scripts discussed in this manual. It also installs all python modules
needed to run basic FLAIR. If you want to use full FLAIR, install the packages
listed under point 5 in the list above.


Download the latest release
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to the FLAIR `release page <https://github.com/BrooksLabUCSC/flair/releases>`_
and select one of the source code files under Assets. Exctract the file and navigate
to the resulting `flair` directory. Add FLAIR and the helper scripts to your ``$PATH``
for instance (in Linux) with ``export PATH=$(pwd)/bin:$PATH``. 

Make sure to (`pip`) install the python modules listed above. If you have conda, you can
create a basic FLAIR environment using

``conda env create -f misc/flair_conda_env.yaml``


Download the latest code
~~~~~~~~~~~~~~~~~~~~~~~~

Check out `the current FLAIR repository <https://github.com/BrooksLabUCSC/flair.git>`_
from github. Please be aware that while this may have the latest bug fixes, it's quite
possible that new bugs were introduced. This method is only useful if you have 
`reported a problem <https://github.com/BrooksLabUCSC/flair/issues>`_ and a FLAIR developer
lets you know it has been fixed.

Once you have cloned the repository, navigate to the `/flair` directory. Follow the
steps as described under Download the latest release.

flair_diffExp and flair_diffSplice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two other programs are available in the source tree, although they now deprecated:
``flair_diffExp``, and ``flair_diffSplice``.  The dependencies are not
include in the standard conda or pip installs.  The following packages are required

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

These maybe added to standard conda environment with:

```
conda env update --name flair --file misc/flair_diffexp_conda_env.yaml
pip install -e .[diffexp]
```
