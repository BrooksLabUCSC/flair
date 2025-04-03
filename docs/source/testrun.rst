Testing flair
=============

Prerequisites:

   * You have all of the FLAIR dependencies install, normally through ``conda``
   * You have a copy of the ``flair/test`` directory (e.g. ``git clone git@github.com:BrooksLabUCSC/flair.git``)

Move to the ``flair/test`` directory, then run ``make test`` or run tests in parallel with
``make -O -j 32 test``.

By default, this uses the FLAIR code in the tree, to teat install FLAIR, use ``make test use_installed_flair=yes``.

If this is the first time, make will download some sequences from 
`the UCSC Genome Browser download page <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes>`__
and store them as test_input/genome.fa.

``make test`` tests all four FLAIR modules and two helper programs. You can also test them individually using:

 - ``make test-align``
 - ``make test-correct``
 - ``make test-collapse``
 - ``make test-quantify``
 - ``make test-predict-productivity``
 - ``make test-diff-iso-usage``

``make`` outputs a lot of information. If a test fails, it will stop with an error and not run any additional tests
unless you specific the ``-k`` option
Errors look like this:

``make: *** [makefile:71: test-predict-productivity] Error 2``

You can usually find more information in the lines preceding the error. If you cannot figure out the problem, please 
`create a ticket <https://github.com/BrooksLabUCSC/flair/issues>`__.



