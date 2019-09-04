Creating a stable release of the code
=====================================

Procedure is:

-  set the version number in phantom/build/Makefile
-  update the :doc:`release notes <releasenotes>`
-  select a gitsha that passes the nightly tests
-  merge this to the ‘staging’ branch
-  create pull request for merge to ‘stable’

The merge to stable can be approved only by developers with special
permission. Nobody has direct write permission to the stable branch, you
*must* issue a pull request from staging. The merge will only be allowed
if the pipeline test is passing
