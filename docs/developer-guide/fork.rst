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

Clone a copy of your fork onto a local machine::

   git clone git@github.com:USERNAME/phantom

push and pull from your fork as you would with the regular phantom
repository

How to update your fork with the latest from the main phantom repo
------------------------------------------------------------------

First, make a remote branch that tracks the main repo::

   git remote add upstream https://github.com/danieljprice/phantom
   git fetch upstream

Then every time you want to update, in your forked copy, type::

   git checkout master
   git fetch upstream
   git merge upstream/master

How to push changes to your fork when you originally cloned the main phantom repo
---------------------------------------------------------------------------------

A common situation is to have checked out a copy from the original repository::

   git clone https://github.com/danieljprice/phantom

and you then make some changes to some files, which you commit to the local repo::

   cd phantom
   ...make some amazing code changes...
   git commit -m 'my amazing code change' file.f90

How should you contribute these back so everyone can benefit? First, you should
*create your fork* as described above. Then you can simply add your new
fork as a remote branch of the current repository::

   git remote add myfork git@github.com:USERNAME/phantom

Notice that in the above we used the ssh address for github, because you need WRITE
permission which is only possible over ssh. If you haven't already done it, you
will need to add your public ssh key to github. To do this go to your .ssh directory::

  cd ~/.ssh
  cat id_rsa.pub
  ... some long key is printed ...

copy everything that was printed above and paste it into the relevant box under
github->settings->SSH and GPG keys, with a name like "my-laptop" or whatever the 
machine you are currently working on is called. You will need to do this once
from every machine you want to push changes from.

If the key exchange was done successfully you should now be able to push your
local changes back to your fork::

   git push myfork

And finally, you can click the "contribute" button which will create a pull request
for your changes to go back to the main phantom repository. Please do this, we
are a community code and everyone benefits when you contribute even small things...

I only have a one line change, should I really issue a whole pull request?
---------------------------------------------------------------------------------

Yes, yes and yes. The easiest pull requests to merge are frequent small changes.
If you complete an entire PhD worth of work and THEN submit a giant pull request
built on changes to a copy of the code you checked out three years ago, it is
difficult (but not impossible) to merge. Frequent, small contributions 
are a much better strategy.

What if I break the code?
-------------------------

That's why we have a comprehensive test suite that runs on every pull request.
Nearly every pull request fails the first time, and requires some tweaking
and improvement to pass all of the integration requirements. But this
process can only begin once you open your pull request!






