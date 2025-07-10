# pycbio select module import

This is an import of selective modules from
<https://github.com/diekhans/pycbio>.  

This is a large collection of personnel tools for computation biology. It is
not distributed in PyPi, and not conda installable, and not ready to generally
share as a dependable package.  Thus select packages and modules are
import into FLAIR.

Python modules are added by checking in the pycbio_import branch with minimal
modification for changing the names to flair.pycbio.  The pycbio tests are not
imported.

To add or update a pycbio module, copy it into appropriate place under src/flair/pycbio/
and run:

```
sed --in-place -r -e 's/^from pycbio/from flair.pycbio/'  src/flair/pycbio/thefile.py
```

Then test and commit.
