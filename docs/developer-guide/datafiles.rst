How to incorporate binary or ascii data files into the Phantom repo
===================================================================

If your setup or module relies on external data files (e.g. an ascii
table) to function, these files need to be distributed with Phantom in
order for your modules to be portable.

Small files
-----------

For *very small files* (under 100Kb), you can simply add these to the git
repository in a subdirectory of the phantom/data directory:

::

   $ cd phantom/data
   $ ls
   README          forcing/        galaxy_merger/      isolatedgalaxy/     star_data_files/
   eos/            gal_radii_cdfs/     galcen/         neutronstar/        velfield/

To use the datafile in your code, say from a file called ‘mydata.txt’
which is stored in the subdirectory ‘star_data_files/red_giant/’, you
should use the **find_phantom_datafile** routine from the **datafile**
module to retrieve the file:

::

   #!fortran
   use datafiles, only:find_phantom_datafile
   character(len=120) :: filename
   ...

   filename=find_phantom_datafile('mydata.txt','star_data_files/red_giant')

This routine searches for the file first in the current directory,
followed by the phantom/data directory once the **PHANTOM_DIR**
environment variable is specified, e.g.:

::

   $ export PHANTOM_DIR=~/phantom
   $ ./phantomsetup

Large files
-----------

For large data files, the procedure is as follows:

1. Do **NOT** add the file to the git repository. Instead, place a
   README file in the directory where the file belongs and add this to
   the git repository:

::

   git add data/star_data_files/red_giant/README
   git commit -m 'placeholder for mydata.txt' data/star_data_files/red_giant/README

2. Add the name of the file to the **.gitignore** file in the root-level
   phantom directory
3. Then, upload your file(s) to a repository on zenodo.org, obtain
   the DOI for the repository, and submit it in the phantom zenodo community

   https://zenodo.org/communities/phantom/

   This will allow the file to be downloaded
   by Phantom users at runtime
4. Edit the datafiles.f90 module to include the URL for the file in the
   **map_dir_to_web** routine. This routine maps the directory where the
   file should be located to the URL where the file can be downloaded
   from. For example, if the file is in the
   **star_data_files/red_giant** directory, you would add a case to the
   routine like so:

::

    case('data/star_data_files/red_giant')
         url = 'https://zenodo.org/records/12345678/files/'

    where the URL is the URL of the repository on zenodo.org

5. Implement the call in the code as previously using the
   find_phantom_datafile routine. This will automatically retrieve the
   file from the web into your phantom/data directory at runtime.
   Alternatively you can manually download the file to the
   appropriate folder
