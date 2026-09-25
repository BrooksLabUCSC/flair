Testing flair
=============

Prerequisites:

   * You have all of the FLAIR dependencies install, normally through ``conda``
   * You have a copy of the ``flair/test`` directory (e.g. ``git clone git@github.com:BrooksLabUCSC/flair.git``)

Move to the ``flair/test`` directory, then run ``make test`` or run tests in parallel with
``make -O -j 32 test``.

By default, this uses the FLAIR code in the tree. To test the installed FLAIR, use ``make test use_installed_flair=yes``.

If this is the first time, make will download some sequences from 
`the UCSC Genome Browser download page <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes>`__
and store them as test_input/genome.fa.

``make test`` tests every FLAIR subcommand and the helper programs. Tests that need ``R`` are
run by ``make test-with-R``; ``make test-base`` runs the rest. You can also run one group at a time:

 - ``make align-tests``
 - ``make transcriptome-tests``
 - ``make quantify-tests``
 - ``make combine-tests``
 - ``make fusion-tests``
 - ``make variants-tests``
 - ``make diffexp-tests``
 - ``make diffsplice-tests``
 - ``make diff-iso-usage-tests``
 - ``make plot-usage-tests``
 - ``make spliceevents-tests``
 - ``make partition-tests``
 - ``make lib-tests``

``make test-help`` checks that every ``--help`` output and the generated command line
documentation in ``docs/source/cli`` still match what the programs accept.

``make`` outputs a lot of information. If a test fails, it will stop with an error and not run any additional tests
unless you specify the ``-k`` option
Errors look like this:

``make: *** [makefile:71: test-predict-productivity] Error 2``

You can usually find more information in the lines preceding the error. If you cannot figure out the problem, please 
`create a ticket <https://github.com/BrooksLabUCSC/flair/issues>`__.



