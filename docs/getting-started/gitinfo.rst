Working with Phantom and git
============================

Make sure you have the git version control system installed.

Getting your first copy
-----------------------

Once you have a GitHub account, you must create your own :doc:`fork </developer-guide/fork>`.
This is done using the “fork” button (the big button on top right of the
repo page).  You can then clone your fork to your computer::

   git clone https://github.com/USERNAME/phantom.git

This gets a copy of the entire phantom repository. Obviously replace
USERNAME with your GitHub username.

Setting your username and email address
---------------------------------------

Before you can push changes, you must ensure that your name and email
address are set, as follows::

   cd phantom
   git config --global user.name "Joe Bloggs"
   git config --global user.email "joe.bloggs@monash.edu"

Please use your full name in the format above, as this is what appears
in the commit logs (and in the AUTHORS file)

Receiving updates from your fork
--------------------------------

Procedure is: stash your changes, pull the updates, reapply your changes::

   git stash
   git pull
   git stash pop

Receiving updates from the master branch
----------------------------------------

Before you can receive updates from the master branch, you must first link
your fork to the master branch::

   git remote add upstream https://github.com/danieljprice/phantom.git

This only needs to be done once.

To update, the procedure is: stash your changes, pull the updates,
reapply your changes::

   git stash
   git fetch upstream
   git merge upstream/master
   git stash pop

This will update your fork on your local machine only.

Committing changes to your fork
-------------------------------

Submit changes to Phantom carefully! The first thing is to pull any
upstream changes as described above. Once you have done this, first
check what you will commit::

   git diff

then go through each subset of changes you have made and commit the
file(s) with a message::

   git commit -m 'changed units in dim file for problem x' src/main/dim_myprob.f90

and so on, for all the files that you want to commit. Then, when you’re
ready to push the changeset back to your fork use::

   git push

Note that you will only be allowed to push changes if you have already
updated your copy to the latest version.

If you have just updated your code from the master repo, simply update
your fork via::

   git commit -m 'merge'
   git push

This will push all the remote changes to your forked version of Phantom.

Committing changes to the master branch
---------------------------------------

This is done through a “pull request”.  To do this,
you can click the “contribute” button on the GitHub page to request
that your changes be pulled into the master copy of Phantom. Please do
this frequently. Many small pull requests are much better than one giant
pull request!

Automated tests will be performed on all pull requests to ensure nothing gets broken. 
Once these pass and the code has been reviewed, the code can be merged.
