running the bots on a pull request
==================================

Wondering how to "run the bots" as requested in the pull request template?
Look no further...

what are the bots and why should I run them?
--------------------------------------------

The bots perform some automated maintenance tasks in phantom, like
updating the file headers, fixing indentation errors and updating the
AUTHORS list. Running the bots yourself allows you to keep ownership over
the lines of code you wrote, e.g. when typing "git blame" on a subroutine
or when looking at who wrote particular lines in VSCode. Otherwise a small
adjustment in indentation means the lines will be attributed to someone else...

performing a dry run
---------------------

By default the bots perform a dry run, making no changes but showing you
what will be changed when you apply the changes.

You can run the bots yourself as follows::

   cd phantom/scripts
   ./bots.sh

doing it for real
---------------------

If you are happy with the dry run, and you have a clean repository (i.e. everything
is committed onto the current branch you are working on), proceed to ACTUALLY
run the bots as follows::

   cd phantom/scripts
   ./bots.sh --commit

Please note this will PULL, COMMIT and PUSH to whatever git repository your current
branch is tracking. If it is tracking the main branch of the main repo
of phantom your push will be disallowed. Instead you should push the changes
to your fork and issue a pull request.

apply but do not commit
-----------------------

In some circumstances you want to apply the changes to the files but leave them
as "modified" in the current directory, which means that the can be reversed
with "git restore" and/or committed manually::

   cd phantom/scripts
   ./bots.sh --apply

running only specific bots
--------------------------

You can run specific bots using the --only flag, e.g.::

   cd phantom/scripts
   ./bots.sh --only authors

or, to run multiple (or all)::

   cd phantom/scripts
   ./bots.sh --apply --only "tabs gt shout header whitespace authors endif"

fixing merge conflicts in the AUTHORS file
-------------------------------------------

you can easily recreate/update the AUTHORS file using the authors bot::

   cd phantom/scripts
   ./bots.sh --apply --only authors
   cd ..
   git add AUTHORS
