Forking Phantom
===============

Why fork?
---------

Working on a fork allows you to work on your own copy of phantom, contributing code back to the main code via a "pull request"

How to fork
-----------

First, click “Fork this repository” on the main phantom repository in
github. You then have your own copy of Phantom in
github.com/USERNAME/phantom

How to work on your fork
------------------------

Clone a copy of your fork onto a local machine

::

   git clone https://github.com/USERNAME/phantom

push and pull from your fork as you would with the regular phantom
repository

How to update your fork with the latest from the main phantom repo
------------------------------------------------------------------

First, make a remote branch that tracks the main repo:

::

   git remote add upstream https://github.com/danieljprice/phantom
   git fetch upstream

Then every time you want to update, in your forked copy, type:

::

   git checkout master
   git fetch upstream
   git merge upstream/master
