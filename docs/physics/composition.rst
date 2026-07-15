Composition tracking in Phantom
=====================================

Tracking the chemical composition of the gas in phantom with fixed
abundances is straightforward, since phantom is a Lagrangian code and
the particle identifiers are preserved throughout the simulation.

In the non-MPI code without particle injection, the particles are 
always written to the dump files in the same order, so the particle
id is simply the particle index in the dump file. In the MPI code,
the particle id is stored in the 'iorig' array in the dump file.

Hence composition tracking can be done as a post-processing step.

Tracking chemical abundances with KROME
=====================================

Phantom contains routines to post-process simulation outputs with KROME, a chemical kinetics package that solves time-dependent chemical evolution according to a specified network of chemical reactions.
KROME has only been robustly tested in Phantom for post-processing, and not for on-the-fly chemistry.
The publicly available routines of Phantom+KROME are meant to be used to process wind simulations, and are not guaranteed to be accurate for other types of simulations.
Below we provide a guide on how to install KROME, and run phantomanalysis to calculate the evolution of chemical species in a simulation.

To calculate chemical abundances in Phantom dumps, you first need to install KROME and HDF5.

Downloading KROME
~~~~~~~~~~~~~~~~~~
You can download KROME from the publicly available repository (https://bitbucket.org/tgrassi/krome/src/master/). For an introduction on how KROME works, follow the tutorials provided in their documentation. 
Make sure to store it in a directory that will be accessible to your Phantom installation. You do not need to compile KROME right now, it will be compiled automatically when you compile Phantom with KROME support.
You should set the environment variable ``KROMEPATH`` to point to the installation directory of KROME: 
::
    export KROMEPATH=/path/to/krome

Installing HDF5
~~~~~~~~~~~~~~~~~~
HDF5 is a library and file format for storing large amounts of data. The KROME routines within Phantom make use of the HDF5 library to read and write output files, because chemical abundance files are typically too large
to be stored in ASCII format.
There are several ways to install HDF5, you can get it directly from your package manager (on Linux or macOS), you can download precompiled binaries from the official release on the GitHub page (https://github.com/HDFGroup/hdf5/releases/tag/2.1.1), 
or you can compile it from source (also using their GitHub repository). If you work on a cluster, it is likely to already be installed as a module, and you should be able to load the HDF5 module installed on the system (for questions, contact your cluster admin).

To decide how to install HDF5, **you must check which compiler you intend to use for Phantom**. When linking libraries, it is important that they are compiled with the same compiler and version, otherwise you may encounter linking errors when compiling your code.
Thus, if you intend to compile Phantom with gfortran, you should install HDF5 with gfortran, same goes for ifort or ifx. Pre-compiled binaries available on the HDF5 GitHub releases or your own package manager are compiled with gfortran or ifx, but not ifort.
If you are using ifort, you will need to compile HDF5 from source.
Another constraint to account for is that **KROME is only parallelized with ifort and ifx**, so if you intend to use KROME in parallel, you will need to compile HDF5 with ifort or ifx as well.
You can technically compile KROME with gfortran, but this will strongly limit the performance of the post-processing routines.

If you already have HDF5 installed, you can check which compiler was used to compile it by running the command ``h5fc -showconfig`` in the terminal.

Below we detail the three installation methods for gfortran, ifort, and ifx, and more specific instructions for macOS users.

gfortran
---------
If you're planning on compiling Phantom with gfortran, you can simply install HDF5 from your package manager, for example on Ubuntu/Debian, with the following commands:
::
    # Ubuntu/Debian
    sudo apt-get install libhdf5-serial-dev

You should then set the environment variable ``HDF5_DIR`` to point to the installation directory of HDF5, which is typically ``/usr/lib/x86_64-linux-gnu/hdf5/serial``:
::
    export HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial
You can check where it is installed by running the command ``h5fc -showconfig`` in the terminal, which will show details about the installation of HDF5.


ifx
---------
You can download ifx precompiled binaries from the official release on the GitHub page (https://github.com/HDFGroup/hdf5/releases/tag/2.1.1). 
Make sure to download the file containing "intel" in the filename, for example ``hdf5-2.1.1-ubuntu-2404-intel.tar.gz``.
The downloaded archive can be extracted:
::
    tar -xvf hdf5-2.1.1-ubuntu-2404-intel.tar.gz

which will create a directory containing another archive, which needs to be extracted as well:
::
    tar -xvf HDF5-2.1.1-Linux.tar.gz
The extracted directory contains a series of subdirectories, including a ``HDF5/``, which contains all the libraries and binaries needed to link HDF5 to Phantom. 
You can move this directory to a permanent location that is more accessible, and set the environment variable ``HDF5_DIR`` to point to it:
::
    export HDF5_DIR=/path/to/HDF5

ifort
-----
Since no precompiled binaries are available for ifort, you will need to compile HDF5 from source. You can download the source code from the official release on the GitHub page (https://github.com/HDFGroup/hdf5/releases/tag/2.1.1).
Extract the downloaded archive to a permanent location and store the location in the environment variable ``HDF5_DIR``:
::
    tar -xvf hdf5-2.1.1.tar.gz
    cd /path/to/hdf5-2.1.1
    export HDF5_DIR=/path/to/hdf5-2.1.1
Then, you can set up the compilation of HDF5 with ifort using the following commands:
::
    cmake --preset ci-StdShar-Intel -B build/ifort \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DCMAKE_INSTALL_PREFIX=$HDF5_DIR \
    -DHDF5_BUILD_FORTRAN=ON \
    -DHDF5_BUILD_JAVA=OFF \
    -DHDF5_ENABLE_JNI=OFF \
    -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_BUILD_TOOLS=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF \
    -DBUILD_TESTING=OFF \
    -DHDF5_ENABLE_PLUGIN_SUPPORT=OFF
You can then compile and install HDF5 with:
::
    cmake --build build/ifort -j
    cmake --install build/ifort
Finally, you need to specify two new environment variables to point to the HDF5 libraries and binaries:
::
    export HDF5INCLUDE=$HDF5_DIR/mod/shared
    export HDF5LIB=$HDF5_DIR/lib

macOS
-----
On macOS, you can only compile Phantom with gfortran, since Intel compilers are not supported anymore. 
This means that you can only compile Phantom and KROME with gfortran, and **you can only use KROME in serial**.
You can install HDF5 with Homebrew:
::
    brew install hdf5
Using the following command, you can store the path to HDF5 (likely something like ``/opt/homebrew/opt/hdf5``) and automatically store it to the environment variable ``HDF5_DIR``:
::
    export HDF5_DIR=$(brew --prefix hdf5)


Compiling Phantom with KROME
~~~~~~~~~~~~~~~~~~
Once both KROME and HDF5 are installed, you can compile Phantom with KROME libraries. 
You can use KROME post-processing on pre-existing models, or you can create one for the occasion.

Once you have selected the model to process, you need to create a setup file for KROME in the directory where your phantom makefile is stored.
This file, called krome.setup, contains compilation options for KROME. For wind models, it can be created using the following command:
::
    echo -e "-n=networks/react_umist\n-iRHS\n-noSinkCheck\n-noRecCheck\n-noTlimits\n-unsafe\n-skipODEthermo\n-skipJacobian" > krome.setup 

Each flag is explained in the KROME documentation, but the most important one is ``-n``, which specifies the chemical network to use. 
The chemical network files are stored in the ``networks/`` directory of the KROME directory, and you may choose any network, 
although **there is no guarantee of the chemical network being suitable for your simulation**.
The network we are using here, ``react_umist`` is a detailed network for winds of AGB stars based on observations and detailed 1D chemical simulations, 
but you may want to explore other networks depending on your needs.

Once the setup file is created, you can compile Phantom with KROME by running:
::
    make KROME=yes

The compilation will only work if the paths to HDF5 and KROME directories are set up properly, so remember to set the environment variables ``KROMEPATH`` and ``HDF5_DIR`` before compiling Phantom.
If you intend to use KROME with Phantom often, you may consider adding these variables to your ``~/.bashrc`` file. 

You can then compile the KROME routines for post-processing by running:
::
    make analysis KROME=yes ANALYSIS=analysis_krome.F90


.. note::
    - You need to have compiled Phantom with KROME support before attempting to compile or use the KROME analysis routines.
    - Forgetting to set KROME=yes will result in a regular compilation of Phantom **without** KROME support, and the analysis routines will not compile.
    - If you need to compile other phantom utilities (phantomsetup, phantommoddump, etc.) while you are using KROME analysis, you need to compile them **with** KROME support as well, otherwise you will get linking errors when trying to run the different routines.
    - Switching between KROME and non-KROME compilations of Phantom is possible, but you need to recompile Phantom with the appropriate flags each time you switch and **you need to clean the build directory** with ``make clean KROME=yes`` to avoid conflicts.
    - When using ``make KROME=yes``, the compilation of Phantom takes longer than usual, since it compiles KROME as well, thus you may want to compile Phantom with KROME only when you need to use the KROME routines. 
    - The compilation of KROME is quiet, so you will not see much output from the compilation process itself, and since it takes a few minutes to complete it may look like it is stuck.
    - If the compilation fails with missing routines or variables with "krome" in their name, KROME likely failed to compile properly. Run ``make clean KROME=yes``, check your environment variables then recompile Phantom with KROME.


Phantomanalysis with KROME
~~~~~~~~~~~~~~~~~~
Your simulation is now ready to be post-processed with KROME the same way you would post-process it with any other analysis routine:
::
    ./phantomanalysis dump_?????

The routines will read the first dump file, if it finds a chemistry file corresponding to the dump (dump_?????.h5), it will read the chemical abundances from that file and use them as initial conditions for the chemical evolution, otherwise it will use default initial abundances.
The default initial abundances are set in ``analysis_krome.F90`` (located in ``src/utils/``), and correspond to AGB stellar winds abundances.

In case you need to change the initial abundances, you can edit the subroutine ``set_initial_abundances`` in ``analysis_krome.F90``, recompile the analysis routines, and run ``phantomanalysis`` again.
Please copy the ``analysis_krome.F90`` file to a different name before editing it, so we can keep the original version intact as you make any changes you want to your copy. You then need to change how ``phantomanalysis`` is compiled:
::
    make analysis KROME=yes ANALYSIS=analysis_krome_mycopy.F90

If your analysis stopped before processing all the dumps, you can restart it by using a bit of shell scripting. You can store the following commands in a bash script and run it, or include them in the script you use to run ``phantomanalysis`` on your cluster:
::
    export FILE=dump_name
    export DIR=/path/to/dumps

    start=00010     #the last dump number that was processed
    files=$(printf '%s\n' $DIR/"$FILE"_* \
    | awk -F_ -v s="$start" '$NF ~ /^[0-9]{5}$/ && ($NF+0)>=s' \
    | sort)

    if [ -z "$files" ]; then
    echo "No files found in $DIR matching pattern $FILE_* with number >= $start"
    exit 1
    fi
    ./phantomanalysis $files
This script will find all the dump files in the directory that match the pattern and have a number greater than or equal to the last processed dump (that you must indicate by hand), and run ``phantomanalysis`` on the list of files.

.. WARNING:: 
    The KROME routines are not guaranteed to work for all types of simulations, and have only been robustly tested for post-processing wind simulations. 
    More importantly, no one but you can guarantee that the chemical network you are using is suitable for your simulation. 
    KROME developers provide a large number of chemical networks, but it is up to you to choose the one that is appropriate for your simulation, or build your own.
    Please make sure to check the KROME documentation and tutorials to understand how to choose a suitable chemical network, and how to set up initial abundances and physical conditions for your simulation.
