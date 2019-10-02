Working with Phantom and git
============================

Make sure you have the git version control system installed. See also :doc:`how to work on a fork <fork>`.

Getting your first copy
-----------------------

Once you have a bitbucket account, use:

::

   git clone https://USERNAME@bitbucket.org/danielprice/phantom

this gets a copy of the entire phantom repository. Obviously replacing
USERNAME with your bitbucket username.

Setting your username and email address
---------------------------------------

Before you can push changes, you must ensure that your name and email
address are set, as follows:

::

   cd phantom
   git config --global user.name "Joe Bloggs"
   git config --global user.email "joe.bloggs@monash.edu"

Please use your full name in the format above, as this is what appears
in the commit logs (and in the AUTHORS file)

Install git-lfs
---------------

In order to read large datafiles stored in the phantom git repo, you
will need git-lfs installed:

::

   git lfs install

If your git installation is too old, you will get an error message:

::

   git: 'lfs' is not a git command. See 'git --help'.

In this case, please ask your system admin to update your copy of git.
Most phantom calculations do not rely on external data files, so not
having git-lfs installed is not a show-stopper. Keep calm and carry on.

Receiving updates
-----------------

Procedure is: stash your changes, pull the updates, replay your changes

::

   git stash
   git pull
   git stash pop

Committing changes
------------------

Submit changes to Phantom carefully! The first thing is to pull any
upstream changes as described above. Once you have done this, first
check what you will commit:

::

   git diff

then go through each subset of changes you have made and commit the
file(s) with a message:

::

   git commit -m 'changed units in dim file for problem x' src/main/dim_myprob.f90

and so on, for all the files that you want to commit. Then, when you’re
ready to push the changeset back to the global repository (i.e. :doc:`HAVING
MADE SURE YOUR CHANGES WILL NOT BREAK ANYONE ELSE’S STUFF <testing>`),
use

::

   git push

Note that you will only be allowed to push changes if you have already
updated your copy to the latest version.

At the moment, everyone in the developer group has permission to push
back changes to the Phantom repositories. This may change in the future
if I think it is a bad idea, in which case I’ll start pulling in changes
explicitly. Compilation failures are automatically checked by the
buildbot, which runs every night there have been commits made to the
code and checks for any compilation failures and/or new compiler
warnings and **emails those responsible**.
