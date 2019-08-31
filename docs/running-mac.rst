Using phantom on macOS
======================

It is *not* recommended to perform 3D simulations on a laptop, but it
can be useful to install phantom in this manner for editing and
debugging code and for simple testing.

Software requirements
---------------------

You will need:

-  Mac Developer tools installed (i.e. get Xcode from the `App
   Store <https://itunes.apple.com/au/app/xcode/id497799835?mt=12>`__)

-  Xcode command line tools

-  `Xquartz <https://www.xquartz.org>`__

-  A Fortran compiler (either
   `gfortran <https://gcc.gnu.org/wiki/GFortranBinaries>`__ or
   `ifort <https://software.intel.com/en-us/fortran-compilers>`__)

-  A copy of `splash <http://users.monash.edu.au/~dprice/splash>`__

-  A reasonably up-to-date copy of git (Large data files in the
   phantom/data directory are handled using `git’s large file
   storage <https://git-lfs.github.com>`__. A recent version of git is
   needed to retrieve these correctly when you \`git clone’)

The easiest way to install the last three is to first install either
`Homebrew <https://brew.sh>`__ or `Macports <http://macports.org>`__.
Using homebrew you can simply type:

::

   brew tap danieljprice/all
   brew install splash

This will **automatically** install the gfortran compiler and the Xcode
command line tools if you haven’t already got them.

How to install phantom
----------------------

Once you have the above, you can just follow the regular :doc:`getting
started guide <running-first-calculation>`.

Known issues / troubleshooting
------------------------------

*Issue*: Default git on Mac does not allow git-lfs, resulting in the
error

::

   git-lfs filter-process: git-lfs: command not found
   fatal: the remote end hung up unexpectedly
   warning: Clone succeeded, but checkout failed.
   You can inspect what was checked out with ‘git status’.
   and retry the checkout with ‘git checkout -f HEAD’

*Solution:* Use homebrew to install the git-lfs package:

::

   brew install git-lfs
   git lfs install

*Issue*: Sometimes if you installed gfortran from binary (which installs
to /usr/local/bin) you can get a clash with the copy of \`as’ installed
by Macports (in /opt/local/bin), giving some problem like:

::

   FATAL:/opt/local/bin/../libexec/as/x86_64/as: I don’t understand ‘m’ flag

*Solution:* Delete or move the \`as’ binary installed by Macports:

::

   cd /opt/local/bin
   sudo mv as as-bak
