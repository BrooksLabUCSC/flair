# Major user-visible changes

## [2.1.2] 2025-04-17
* Address issue getting BioConda to work
* Bug fixes for collapse command line parsing

## [2.1.1] 2025-04-10
* converted all programs to use console scripts to allow BioConda to work

## [2.1.0] 2025-03-27
* Numerous bug fixes.
* Removed support for PSL format.
* Remove `flair 123` to run multiple modules at once.
* Compatibility with Python 3.12 
* Compatibility with Apple ARM64 systems.
* Deprecated `diffExp` and `diffSplice`, they will be removed in a future release.
  Lets us know if you use this functionality.  Their dependencies are no longer
  part of the conda package, they can be added the conda environment with
  `misc/flair_diffexp_conda_env.yaml`.
