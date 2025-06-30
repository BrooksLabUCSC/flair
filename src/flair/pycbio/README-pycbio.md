# pycbio select module import

This is an import of selective modules from
<https://github.com/diekhans/pycbio>.  

This is a large collection of personnel tools for computation biology. It is
not distributed in PyPi, and not conda installable, and not ready to generally
share as a dependable package.  Thus select packages and modules are
import into FLAIR.

Python modules are added by checking in the pycbio_import branch with minimal
modification, including changing package names to flair.pycbio.   The pycbio
tests are not imported.

Step to add or update modules.  The unmodified files are first commited,
the updated so that git has the edit history.


  git checkout pycbio_import
  git merge dev
  # copy in or update files and add
  git commit -am 'import of pycbio blah'
  # edit files to fixed imports, etc

  make pycbio-flake8
  make -j 32 test

  git commit -am 'migration of pycbio blah'
  git push
  git checkout dev
  git merge pycbio_import
