Creating a stable release of the code
=====================================

Procedure is:

- set the version number in phantom/build/Makefile
- set the version number in phantom/docs/conf.py
- update the :doc:`release notes <releasenotes>`
- use git to tag the code version for the release

::

    git tag 'v2055.0.1'

- push the tag and let the github actions do the rest

::
     
    git push -v tags

