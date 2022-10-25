Testing flair
=============

Prerequisites:

   * flair and flair scripts are in your ``$PATH`` (see below)
   * You have a copy of the ``flair/test`` directory (e.g. ``git clone git@github.com:BrooksLabUCSC/flair.git``)
   * GNU make

Flair is in your ``$PATH`` if you used ``conda install -c conda-forge -c bioconda flair``.

If you downloaded the latest release from github or cloned the flair repository:

.. code:: text

   export PATH=/path/to/flair/src/flair:/path/to/flair/bin:$PATH


Move to the ``flair/test`` directory, then run ``make test``.

If this is the first time, make will download some sequences from 
`the UCSC Genome Browser download page <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes>`__
and store them as test_input/genome.fa.

``make test`` tests all six flair modules and two helper programs. You can also test them individually using:

 - ``make test-align``
 - ``make test-correct``
 - ``make test-collapse``
 - ``make test-quantify``
 - ``make test-diffexp``
 - ``make test-diffsplice``
 - ``make test-predict-productivity``
 - ``make test-diff-iso-usage``

``make`` outputs a lot of information. If a test fails, it will stop with an error and not run any additional tests.
Errors look like this:

``make: *** [makefile:71: test-predict-productivity] Error 2``

You can usually find more information in the lines preceding the error. If you cannot figure out the problem, please 
`create a ticket <https://github.com/BrooksLabUCSC/flair/issues>`__.



