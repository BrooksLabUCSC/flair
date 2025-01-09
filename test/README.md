# Testing your flair installation

Prerequisites:
 - GNU make

Normally these test use FLAIR code that is in the tree to test a conda install,
run with `make test use_installed_flair=yes`.

Flair is in your $PATH if you used `conda install -c conda-forge -c bioconda flair`.

Then run `make test`.

If this is the first time, make will download some sequences from 
[the UCSC Genome Browser download page](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/)
and store them as test_input/genome.fa.

`make test` tests all six flair modules and two helper programs. You can also test them individually using
 - `make test-align`
 - `make test-correct`
 - `make test-collapse`
 - `make test-quantify`
 - `make test-diffexp`
 - `make test-diffsplice`
 - `make test-predict-productivity`
 - `make test-diff-iso-usage`

`make` outputs a lot of information. If a test fails, it will stop with an error and not run any additional tests.
Errors look like this:

`make: *** [makefile:71: test-predict-productivity] Error 2`

You can usually find more information in the lines preceding the error. If you cannot figure out the problem, please 
[create a ticket](https://github.com/BrooksLabUCSC/flair/issues).



