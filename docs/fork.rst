Forking Phantom
===============

Why fork?
---------

You can work on a fork of Phantom if you wish to keep local changes
without having to push them into the public code. We discourage this in
general, because it is nearly impossible to merge changes in code which
diverged from the master code a long time ago. But, warning heeded, here
is how to do it…

How to fork
-----------

First, click “Fork this repository” on the main phantom repository in
bitbucket. You then have your own bitbucket copy of Phantom in
bitbucket.org/USERNAME/phantom

How to work on your fork
------------------------

Clone a copy of your fork onto a local machine

::

   git clone https://bitbucket.org/USERNAME/phantom

push and pull from your fork as you would with the regular phantom
repository

How to update your fork with the latest from the main phantom repo
------------------------------------------------------------------

First, make a remote branch that tracks the main repo:

::

   git remote add upstream https://bitbucket.org/danielprice/phantom
   git fetch upstream

Then every time you want to update, in your forked copy, type:

::

   git checkout master
   git fetch upstream
   git merge upstream/master
