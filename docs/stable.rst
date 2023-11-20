How to use download and use stable code releases
================================================

The master branch of phantom develops rapidly, but if you want to only
:doc:`stable versions of phantom <releasenotes>`, the best way is to use
the \`stable’ code branch from the git repository. **Nobody has
permission** to push directly to the stable branch, it must be updated
via pull request, so it **cannot** be inadvertently broken

Obtaining a stable copy of the code
-----------------------------------

::

   $ git clone https://github.com/danieljprice/phantom
   $ cd phantom
   $ git checkout stable

Updating to the latest stable copy of the code
----------------------------------------------

Simply use \`git pull’ while on the stable branch:

::

   $ git checkout stable
   $ git pull

Switching back to the developer version
---------------------------------------

::

   $ cd phantom
   $ git checkout master
