How to compute density on a set of particles with Phantom
=========================================================
The basic steps are:
 1. install and compile phantom
 2. create phantom-compatible initial conditions
 3. run phantom for zero timesteps to compute the density field
 4. perform visualisation and analysis (e.g. PDF calculation) on phantom data file

Installing phantom
-------------------
Grab a copy of phantom as usual::

    git clone https://github.com/danieljprice/phantom


Make a directory somewhere else (not inside the phantom directory) to do the work::

    mkdir $HOME/mycalc
    cd $HOME/mycalc

Compile phantom in this directory using the "writemake" script. We will also use the compile time configuration for the `turbdrive` setup, since we want a periodic box::

    ~/Codes/phantom/scripts/writemake.sh turbdrive > Makefile

We can compile the code as usual, but switching off the turbulent driving::

     export SYSTEM=gfortran
     make SRCTURB=""

Creating phantom-compatible initial conditions
----------------------------------------------

I copied the setup_empty.F90 from phantom/src/setup/ as a template for a new initial conditions file::

    cp ~/Codes/phantom/src/setup_empty.f90 ./setup_fromfile.f90

We then want to edit this file to read the supplied datafiles. I edited the file as follows::

    subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
     use units,     only:set_units
     use physcon,   only:pc,solarm
     use prompting, only:prompt
     use boundary,  only:set_boundary
     use  part, only:kill_particle,shuffle_part
     integer,           intent(in)    :: id
     integer,           intent(inout) :: npart
     integer,           intent(out)   :: npartoftype(:)
     real,              intent(out)   :: xyzh(:,:)
     real,              intent(out)   :: massoftype(:)
     real,              intent(out)   :: polyk,gamma,hfact
     real,              intent(inout) :: time
     character(len=20), intent(in)    :: fileprefix
     real,              intent(out)   :: vxyzu(:,:)
     integer :: ierr,iunit,idn
     character(len=60) :: filename

     !
     ! set code units
     !
     call set_units(dist=1.*pc,mass=solarm,G=1.)
     !
     ! set boundary for periodic boundaries
     !
     call set_boundary(0.,4.,0.,4.,0.,4.)
     !
     ! default time, sound speed squared and gamma
     !
     time = 0.
     polyk = 0.
     gamma = 1.

     npart = 0
     npartoftype(:) = 0
     massoftype = 1.
     filename = 'bin5_tracer.dat'

     call prompt('enter filename to read positions',filename)
     open(newunit=iunit,file=filename,status='old',iostat=ierr)
     do while (ierr == 0)
      npart = npart + 1
      read(iunit,*,iostat=ierr) idn,xyzh(1:3,npart)
      if (ierr /= 0) npart = npart - 1
     enddo
     xyzh(4,:) = 0.01
     print*,' read ',npart,' particles '
     polyk = 1.
     ! set particle number and mass (arbitrary)
     npartoftype(1) = npart
     massoftype(1) = 1.e6/npart

    end subroutine setpart

In the above we set the periodic boundaries to match those of the input file (in this example, :math:`x,y,z \in [0,4]`)

One problem we encountered was that some particles had identical positions. To avoid this I found one of the particles that matched and deleted the overlapping particles using the kill_particle and shuffle_part routines in phantom::

    print*,'ref=',xyzh(1:3,489543)
    do idn=1,npart
       if (norm2(xyzh(1:3,idn) - xyzh(1:3,489543)) < 1.e-6) then
         print*,idn,xyzh(1:3,idn)
         call kill_particle(idn)
      endif
    enddo
    call shuffle_part(npart)
    npartoftype(1) = npart
    massoftype(1) = 1.e6/npart

We then compile phantomsetup using the new setup file::

    make setup SRCTURB="" SETUPFILE=setup_fromfile.f90

Followed by running phantomsetup to create phantom-compatible initial conditions
::

    ./phantomsetup dens

This should result in creation of a phantom input file (dens.in) and an initial conditions dump (dens_0000.tmp)::

     -------->   TIME =    0.000    : full dump written to file dens_00000.tmp   <--------


     input file dens.in written successfully.
     To start the calculation, use:

     ./phantom dens.in

You should check that the files indeed exist
::

    $ ls
    Makefile phantom phantomsetup dens.in dens_00000.tmp

Running phantom to compute the density
--------------------------------------
To enable phantom to ONLY compute the density and force, we just need to set the parameter "nmax" to zero in the input file. Open dens.in in your text editor and change the corresponding line to read::

    nmax =           0    ! maximum number of timesteps (0=just get derivs and stop)

Then run the code as usual
::

    $ ./phantom dens.in

which should produce an output file WITH density calculated::

    $ ls dens*
    dens.in dens01.ev dens_00000

You can visualise the density field in this file using splash, e.g.::

    $ brew tap danieljprice/all
    $ brew install splash
    $ ssplash dens_00000 -r 6 -dev /xw

Computing the volume weighted PDF of the density
------------------------------------------------
To compute the volume-weighted density probability density function (PDF), we need to interpolate the density field to a mesh. One needs to be careful to resolve the smoothing length of each particle when doing this, otherwise the normalisation of the interpolation will not be correct. In Price & Federrath (2010) and Price, Federrath & Brunt (2011) we solved this by interpolating to an adaptive mesh and computing the PDF on this mesh. For this we can compile the phantom2pdf-amr utility::

    $ make phantom2pdf-amr

which complains that it needs some files from the SPLASH source code::

    ERROR: cannot find SPLASH directory needed for some source files

so we just need a clone of the splash source code::

    $ cd /tmp/
    $ git clone https://github.com/danieljprice/splash
    $ export SPLASH_DIR=/tmp/splash
    $ make phantom2pdf-amr

We can then run this utility as follows::

    $ ./phantom2pdf-amr dens_00000

The first run produces the error::

    ERROR opening analysis.in
    WRITING ANALYSIS OPTIONS TO analysis.in
    RERUN THE ANALYSIS AFTER EDITING THIS FILE

The file analysis.in contains the input parameters for phantom2pdf-amr, which you should edit to specify the box dimensions and the desired binning for the PDF::

    $ more analysis.in
    # box dimensions
                    xmin =       0.000    !  min boundary
                    xmax =       4.000    !  max boundary
                    ymin =       0.000    !  min boundary
                    ymax =       4.000    !  max boundary
                    zmin =       0.000    !  min boundary
                    zmax =       4.000    !  max boundary
    # analysis options
               rhologmin =        -15.    !  min ln(rho) for PDF
               rhologmax =         15.    !  max ln(rho) for PDF
              binspacing =       0.100    !  bin width in ln(rho) for PDF

The final problem we encountered was that the number of meshes exceeded the memory limits::

      ERROR: nmesh > maxmeshes (1e6): change parameter and recompile

To fix this, just edit relevant line in ~/phantom/src/utils/adaptivemesh.f90 as follows::

      integer, parameter :: maxmeshes = 5e6

and recompile phantom2pdf-amr::

     rm phantom2pdf-amr; make phantom2pdf-amr

A successful run should produce an output file::

     writing to dens_00000_pdf_lnrho.dat

With contents similar to::

     # volume weighted PDF, calculated using Phantomanalysis: PDF on AMR mesh, part of Phantom v1.3.0 (c) 2007-2019 The Authors
     #   300 bins evenly spaced in lnrho
        3.0590232050182579E-007   1.3448621777151119E-005
        3.3807434839047367E-007   3.1021988445730148E-004
        3.7362993798852602E-007   1.2105637898566894E-005

This can be plotted in Python, or with splash as follows::

     splash -ev dens_00000_pdf_lnrho.dat
