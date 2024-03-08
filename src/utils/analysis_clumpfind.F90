!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which runs the CLUMPFIND algorithm
!
!  This uses local potential minima (and neighbour lists)
!  to identify objects/clumps in the gas distribution
!  Clumps are matched across timesteps using standard merger tree algorithms
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, getneighbours, part, prompting, ptmass,
!   readwrite_dumps, sortutils, units
!
 use dim,             only:maxp
 use getneighbours,    only:generate_neighbour_lists, read_neighbours, write_neighbours, &
                           neighcount,neighb,neighmax
 implicit none
 character(len=20), parameter, public :: analysistype = 'clumpfind'

 real,    parameter :: sinkclumprad = 1.0    ! particles within sinkclumprad*hacc are automatically part of the clump
 real,    parameter :: rhomin_cgs   = 0.0    ! particles with rho < rhomin will not be considered for the lead particle in a clumps
 integer, parameter :: nclumpmax    = 1000   ! maximum number of clumps that can be analysed

 integer            :: runningclumpmax
 logical            :: checkbound, sinkpotential,skipsmalldumps, firstdump
 character(len=100) :: firstdumpfile, previousdumpfile, rawtag, proctag

 public :: do_analysis

 ! Derived type sphclump records clump data
 ! (Total particle counts, clump position/velocity, mass and extent)

 type sphclump

    integer           :: num,mbp,pointmass,ID,border
    real,dimension(3) :: r,v
    real              :: mass,size
    real              :: kinetic,potential,thermal,bound,virial

 end type sphclump

 type(sphclump), allocatable, dimension(:) :: clump,oldclump
 integer,        allocatable, dimension(:) :: member,oldmember
 real,           allocatable, dimension(:) :: dpoten
 integer                                   :: nclump

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,            only:iphase,maxphase,igas,get_partinfo,ihacc,poten,rhoh
 use readwrite_dumps, only:opened_full_dump
 use units,           only:unit_density
 use sortutils,       only:indexx

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter                :: neighcrit = 40
 integer, dimension(1)             :: maxclump
 integer, allocatable,dimension(:) :: neighclump, ipotensort
 integer            :: k,l,iclump,ipart,jpart,iamtypei,deletedclumps
 real               :: percent,percentcount,rhomin
 logical            :: existneigh,iactivei,iamdusti,iamgasi
 logical            :: write_raw_data,write_neighbour_list,fexists
 character(len=100) :: neighbourfile
 character(len=100) :: fmt

! Initialise variables
 rhomin         = rhomin_cgs/unit_density
 rawtag         = 'raw'
 proctag        = 'proc'
 checkbound     = .false.
 sinkpotential  = .false.
 skipsmalldumps = .false.
 firstdump      = .false.
 write_raw_data = .false.
 write_neighbour_list = .false.
! Read in input parameters from file (if it exists)
 call read_analysis_options(dumpfile)

! Print warning
 if (checkbound) then
    print*, '*************************************************************'
    print*, '*                                                           *'
    print*, '*                         WARNING!!!                        *'
    print*, '*                                                           *'
    print*, '*     Potential energy of a clump is the sum of poten,      *'
    print*, '*        Sum_i=1^Nclump  Sum_j=1^Npart Gm_jm_i/r_ij         *'
    print*, '*     rather than being the potential of just the clump,    *'
    print*, '*        Sum_i=1^Nclump  Sum_j=1^Nclump Gm_jm_i/r_ij        *'
    print*, '*     Therefore the potential energy is much too large!     *'
    print*, '*     Consider using analysis_clumpfindTD.F90               *'
    print*, '*                                                           *'
    print*, '*************************************************************'
 endif

! Skip small dumps (as they do not include velocity data)
 if (skipsmalldumps .and. .not.opened_full_dump) then
    print*, 'Skipping clumpfind analysis for small dumps (velocity data missing)'
    return
 endif

 !***************************************
 !1. Obtain Neighbour Lists
 !****************************************

 print '(a)', 'Obtaining Neighbour Lists'

 ! Check if a neighbour file is present

 neighbourfile = 'neigh_'//TRIM(dumpfile)
 inquire(file=neighbourfile,exist = existneigh)

 if (existneigh) then
    print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
    call read_neighbours(neighbourfile,npart)
 else
    ! If there is no neighbour file, generate the list
    print*, 'No neighbour file found: generating'

    call generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile,write_neighbour_list)

 endif

 print*, 'Neighbour lists acquired'

!****************************************************************************************
!-- 2. Sort particles by the modulus of their gravitational potential (descending order)
!-- Indexx requires real*8 input, and poten is real*4
!****************************************************************************************

 allocate(ipotensort(npart)) ! Array for sorting particles by potential
 allocate(dpoten(npart)) ! Holding array for potential (real*8)


 dpoten = DBLE(poten)

! Add potential contribution from all sinks first

 if (sinkpotential) call add_sink_potential(npart)

 dpoten = abs(dpoten)

! Now do sort
 call indexx(npart,dpoten,ipotensort)

!*************************************************************************
!-- 3. Begin CLUMPFIND procedure
!
! In brief, the algorithm travels down the sorted list of particles
! Any particle with a neighbour in clump i is also in clump i
! Any particle with no neighbours in clumps begins a new clump j

! Particles with neighbours in multiple clumps select the clump with most
! neighbours to join
!
!*************************************************************************

!
!--Allocate clump array data (uses 'sphclump' type)
!
 allocate(clump(nclumpmax))      ! Maximum number of allowed clumps = npart
 print*, 'Clump data allocated'
 allocate(member(npart))         ! Clump membership for each particle
 allocate(neighclump(nclumpmax)) ! For a set of particle neighbours, this tracks which clumps they are in
 print*, 'Catalogues allocated'

 nclump    = 0                   ! Records total number of clumps counted so far
 member(:) = 0                   ! Set membership for all particles to zero

! Set positions and velocities of clumps to zero
 do k=1,3
    clump(:)%r(k) = 0.0
    clump(:)%v(k) = 0.0
 enddo

! Set mass, maximum radius and total number of part`icles per clump to zero
 clump(:)%ID        = 0
 clump(:)%mass      = 0.0
 clump(:)%size      = 0.0
 clump(:)%num       = 0
 clump(:)%mbp       = 0
 clump(:)%kinetic   = 0.0
 clump(:)%potential = 0.0
 clump(:)%thermal   = 0.0
 clump(:)%pointmass = 0     ! Indicates if clump contains a pointmass (and which)
 clump(:)%border    = 0     ! Indicates if a clump has particles across the boundary
!
!--All sinks are automatically identified as clumps
!
 call create_sink_clumps(npart, xyzh)
!
!
!--Now begin CLUMPFIND analysis on SPH particles
!**************************************************************************************
! Algorithm:
! If particle not in a clump AND none of its neighbours in a clump, create a new clump
! If majority of particle's neighbours in a clump, particle also in that clump
!**************************************************************************************
!
!--Begin loop over particles
!
 percentcount = 10.0

 fmt = "(a,I4,a,I8,2X,1pg10.3,2X,1pg10.3,2x,1pg10.3)"

 over_parts: do l = 1,npart

    percent = 100.0*REAL(l)/REAL(npart)

    if (percent > percentcount) then
       write(*,'(I3," % complete")') int(percentcount)
       percentcount = percentcount +10.0
    ENDif

    ! Pick next particle in order of descending potential
    ipart = ipotensort(npart-l+1)

    ! skip particle if not dense enough
    if (rhoh(xyzh(4,ipart),particlemass) < rhomin) cycle over_parts

    ! Only calculate for gas particles
    if (maxphase==maxp) then
       call get_partinfo(iphase(ipart),iactivei,iamgasi,iamdusti,iamtypei)
       ! If particle isn't a gas particle, skip it
       if (.not.iamgasi) cycle over_parts
    endif

    ! If particle not already in a clump
    if (member(ipart)==0) then
       ! Reset array for counting neighbour particle clump membership
       neighclump(:) = 0

       ! Check if its neighbours are in a clump
       do jpart=1,neighcount(ipart)
          if (neighcount(ipart) > neighmax) then
             print*, 'finding clumps A.  neighcount(ipart) > neighmax.  aborting'
             stop
          endif
          iclump = member(neighb(ipart,jpart))
          if (iclump > 0) then
             neighclump(iclump) = neighclump(iclump)+1
          endif
       enddo

       ! Compare neighbour membership, identify clump with most neighbours as members
       ! This is now particle ipart's host clump
       if (sum(neighclump)/=0) then
          maxclump = maxloc(neighclump)
          call add_particle_to_clump(ipart,maxclump(1))
       endif
    endif

    ! If particle still does not belong to any other clump, then create a new clump for it
    if (member(ipart)==0 .and. neighcount(ipart) > neighcrit) call initialise_clump(ipart)

    ! If particle is already in a clump, then add any neighbours not in a clump
    if (member(ipart) > 0) then
       do jpart=1,neighcount(ipart)
          if (neighcount(ipart) > neighmax) then
             print*, 'finding clumps B.  neighcount(ipart) > neighmax.  aborting'
             stop
          endif
          if (member(neighb(ipart,jpart))==0) call add_particle_to_clump(neighb(ipart,jpart), member(ipart))
       enddo
    endif

 enddo over_parts

 deallocate(ipotensort)

! Calculate boundness of all clumps

 do iclump=1,nclump
    call calc_clump_boundness(iclump)
 enddo

 print*, 'Initial clump identification complete'

 print "(/,65('-')/)"

!
!--End initial clump identification process
!

 print*, 'Clumps identified without gravitational boundness checking:'

 call print_clump_data

!*****************************************************
!
!--4. Check clumps for gravitational boundness
!-- clumps have outermost particles deleted until they are bound (or totally deleted)
!
!*****************************************************

 deletedclumps = 0
 if (checkbound) then
    call test_clump_boundness(deletedclumps,npart,xyzh,particlemass)
 endif

!****************************************************
!
!--5. Write raw (unmatched) clump data to file
!
!****************************************************
 if (write_raw_data) then
    print*, 'Clumpfind complete for snapshot ', dumpfile
    print*, 'Writing unmatched clump data to raw files'
    call write_clump_data(nclump,deletedclumps,npart, time,dumpfile,rawtag)
 endif
!****************************************************
!
!--6. Link clumps in this dump with the previous dump
!
!****************************************************

! Check if this is the first file to be analysed
! If so, then we'll simply write the data to file
! as these are the original clumps

 inquire(file="proc_clumpcat_"//trim(previousdumpfile),exist=fexists)
 print*, "proc_clumpcat_"//trim(previousdumpfile),firstdump,.not.fexists
 if (firstdump .or. .not.fexists) then
    k = 0
    do iclump = 1,nclump
       if (clump(iclump)%ID > 0) then
          k = k + 1
          clump(iclump)%ID = k          ! Reset all the clump IDs for cleanlieness
       endif
    enddo
    runningclumpmax = k                 ! Update current maximum clump count and save it
    call amend_options_file(dumpfile)
 else
    ! Match clump IDs across dump files
    call merger_tree(npart,dumpfile)
 endif

! Now write to file
 call write_clump_data(nclump,deletedclumps,npart, time,dumpfile,proctag)

! Deallocate memory for next dump
 deallocate(neighcount,neighb,dpoten)
 deallocate(member,neighclump,clump)

end subroutine do_analysis
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Read in options for analysis from file
!+
!-----------------------------------------------------------------------
subroutine read_analysis_options(dumpfile)
 use prompting, only:prompt

 character(len=*), intent(in) :: dumpfile
 logical :: inputexist
 character(len=17) :: inputfile
 character(len=1) :: boundchoice, sinkchoice,skipchoice

! Check for existence of input file
 inputfile = 'clumpfind.options'
 inquire(file=inputfile, exist=inputexist)

 if (inputexist) then
    print '(a,a,a)', "File ",inputfile, " found: reading analysis options"
    open(10,file=inputfile, form='formatted')
    read(10,*) boundchoice
    read(10,*) sinkchoice
    read(10,*) skipchoice
    read(10,*) firstdumpfile

    ! If this is a file from a previous run, then we need to be careful
    ! about what the previous dump was
    if (trim(dumpfile)==trim(firstdumpfile)) then
       write(10,*) trim(dumpfile), "    Previous SPH dump analysed"
       runningclumpmax = 0
    else
       read(10,*) previousdumpfile
       read(10,*) runningclumpmax
    endif
    close(10)

    if (boundchoice=='y' .or. boundchoice=='Y') checkbound     = .true.
    if (sinkchoice =='y' .or. sinkchoice =='Y') sinkpotential  = .true.
    if (skipchoice =='y' .or. skipchoice =='Y') skipsmalldumps = .true.
    if (trim(firstdumpfile)==trim(dumpfile))    firstdump      = .true.
 else
    boundchoice = 'n'
    sinkchoice  = 'n'
    skipchoice  = 'n'
    firstdump   = .true.

    call prompt('Do you want to test for boundness?', checkbound, .true.)
    call prompt('Do you want to include sink contribution to the potential?', sinkpotential,.true.)
    call prompt('Do you want to skip small dumps (missing velocity data)? ',skipsmalldumps,.true.)

    if (checkbound)     boundchoice = 'y'
    if (sinkpotential)  sinkchoice  = 'y'
    if (skipsmalldumps) skipchoice  = 'y'

    ! Write choices to new inputfile
    open(10,file=inputfile, status='new', form='formatted')
    write(10,*) boundchoice, "      Test clumps for gravitational boundness?"
    write(10,*) sinkchoice,  "      Include sinks' contribution to potential?"
    write(10,*) skipchoice,  "      Skip small dumps (velocity data missing)?"
    write(10,*) trim(dumpfile), "    First SPH dump to be analysed"
    write(10,*) trim(dumpfile), "    Previous SPH dump analysed"
    close(10)
    previousdumpfile = ""
    runningclumpmax  = 0
 endif

 if (checkbound) then
    print '(a)', 'Clumps will be checked for gravitational boundness'
 else
    print '(a)', 'Clumps will NOT be checked for gravitational boundness'
 endif

 if (sinkpotential) then
    print '(a)', 'Contribution of sinks to the potential will be added'
 else
    print '(a)', 'Contribution of sinks to the potential will NOT be added'
 endif

 if (skipsmalldumps) then
    print '(a)', 'Skipping small dumps'
 else
    print '(a)', 'NOT skipping small dumps'
 endif

 if (firstdump) then
    print '(a)', 'This is the first dump to be analysed'
 else
    print '(a)', 'This is not the first dump - will attempt to link clump IDs to previous dump'
 endif

 if (checkbound .and. .not.(skipsmalldumps)) then
    print*, 'ERROR in read_analysis_options!'
    print*, 'Cannot check boundness using small dumps!'
    print*, 'Please rerun, either choosing to check boundness or to use small dumps, not both'
    stop
 endif

end subroutine read_analysis_options
!-----------------------------------------------------------------------
!+
! Amend the options file to update which dump was previously analysed
!+
!-----------------------------------------------------------------------
subroutine amend_options_file(dumpfile)
 character(len=*),intent(in) :: dumpfile

 ! Open the options file, and wind forward to the line of interest
 open(10,file='clumpfind.options', form='formatted')
 read(10,*)
 read(10,*)
 read(10,*)
 read(10,*)
 write(10,*) trim(dumpfile), "    Previous SPH dump analysed"
 write(10,*) runningclumpmax, " Current maximum clump number"
 close(10)

end subroutine amend_options_file
!-----------------------------------------------------------------------
!+
! Add contribution to potential from sinks
! This currently ignores periodicity
!+
!-----------------------------------------------------------------------
subroutine add_sink_potential(npart)
 use part,    only:xyzh,massoftype,igas,xyzmh_ptmass,nptmass
 use ptmass,  only:get_accel_sink_gas
 integer, intent(in) :: npart
 integer             :: ipart
 real                :: dumx,dumy,dumz, epoti

 print '(a)', 'Calculating sink contribution to the gravitational potential'
 do ipart=1,npart
    ! Dummy variables for gravitational forces (not used here)
    dumx  = 0.0
    dumy  = 0.0
    dumz  = 0.0
    epoti = 0.0

    call get_accel_sink_gas(nptmass,xyzh(1,ipart),xyzh(2,ipart),xyzh(3,ipart),xyzh(4,ipart),&
                            xyzmh_ptmass,dumx,dumy,dumz,epoti)
    dpoten(ipart) = dpoten(ipart) + massoftype(igas)*epoti
 enddo
end subroutine add_sink_potential

!-----------------------------------------------------------------------
!+
! Create clumps at the location of all sinks
!  Note: Sinks do not possess neighbours in Phantom
!        All particles within a given radius of the sink are added to the clump
!+
!-----------------------------------------------------------------------
subroutine create_sink_clumps(npart,xyzh)
 use part,     only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 integer             :: k,iclump,ipart,jpart,iptmass
 real                :: clumpsep,sep,dx,dy,dz
 logical             :: first_particle
 character(len=100)  :: fmt

 print "(/,65('-')/)"

 print*, 'Automatically identifying sink particles as progenitors of clumps'
 print'(a,F5.1,a)', 'All particles within ',sinkclumprad,' accretion radii  will be added to sink clumps'

 do iptmass=1,nptmass
    if (xyzmh_ptmass(4,iptmass) < 0.) cycle
    clumpsep = sinkclumprad*xyzmh_ptmass(ihacc,iptmass)
    nclump   = nclump + 1

    call check_clump_max

    clump(iptmass)%ID        = iptmass
    clump(iptmass)%pointmass = iptmass
    clump(iptmass)%size      = clumpsep
    clump(iptmass)%num       = 1            ! We are allowing sinks to be groups of one

    do k=1,3
       clump(iptmass)%r(k) = xyzmh_ptmass(k,iptmass)
       clump(iptmass)%v(k) = vxyz_ptmass(k,iptmass)
    enddo
    clump(iptmass)%mass = xyzmh_ptmass(4,iptmass)

    print '(a,I4,a)', 'Pointmass ',iptmass, ' begins a clump:'
    print '(1pg10.3,2X,1pg10.3,2X,1pg10.3,2X,1pg10.3)', clump(nclump)%r(1), &
        clump(nclump)%r(2), clump(nclump)%r(3), clump(nclump)%mass

    ! Test each particle for distance
    first_particle = .true.
    do ipart=1,npart

       dx = xyzh(1,ipart) - xyzmh_ptmass(1,iptmass)
       dy = xyzh(2,ipart) - xyzmh_ptmass(2,iptmass)
       dz = xyzh(3,ipart) - xyzmh_ptmass(3,iptmass)
#ifdef PERIODIC
       if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
       if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
       if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
       sep = sqrt( dx*dx + dy*dy + dz*dz )

       if (sep < clumpsep .and. member(ipart)==0) then

          call add_particle_to_clump(ipart,iptmass)

          if (first_particle) then
             clump(iptmass)%mbp = ipart
             first_particle = .false.
          endif

          ! Add all SPH particle's neighbours to this clump
          do jpart=1,neighcount(ipart)
             if (member(neighb(ipart,jpart))==0) call add_particle_to_clump(neighb(ipart,jpart), iptmass)
          enddo
       endif
    enddo ! End of loop over particles
 enddo ! End of loop over sinks

 print'(a,I2,a)', 'Clumps created for ',nptmass, ' pointmasses'
 print*, 'ID x y z mass particle number'
 fmt = '(i4,2x,1pg10.3,2x,1pg10.3,2x,1pg10.3,2x, 1pg10.3,2X, I7)'
 do iclump=1,nclump
    print fmt , iclump, clump(iclump)%r(1),clump(iclump)%r(2),clump(iclump)%r(3),clump(iclump)%mass, clump(iclump)%num
 enddo

 print "(/,65('-')/)"

end subroutine create_sink_clumps
!-----------------------------------------------------------------------
!+
! Create a new clump at the location of SPH particle ipart
!+
!-----------------------------------------------------------------------
subroutine initialise_clump(ipart)
 use part, only: xyzh, vxyzu, massoftype,igas
 integer, intent(in) :: ipart
 character(len=100)  :: fmt
 integer             :: k

 fmt    = "(a,I4,a,I8,2X,1pg10.3,2X,1pg10.3,2x,1pg10.3)"
 nclump = nclump + 1
 call check_clump_max

 clump(nclump)%ID  = nclump
 member(ipart)     = nclump
 clump(nclump)%num = clump(nclump)%num + 1

 print fmt, 'Clump ', nclump, ' created for particle ', ipart,&
     xyzh(1,ipart),xyzh(2,ipart),xyzh(3,ipart)

! Store most bound particle index here
 clump(nclump)%mbp  = ipart

! Initialise clump with the particle's properties
 clump(nclump)%mass = massoftype(igas)
 do k=1,3
    clump(nclump)%r(k) = xyzh(k,ipart)
    clump(nclump)%v(k) = vxyzu(k,ipart)
 enddo
 clump(nclump)%thermal   = vxyzu(4,ipart)*massoftype(igas)
 clump(nclump)%potential = dpoten(ipart)
 clump(nclump)%size      = 0.0

end subroutine initialise_clump
!-----------------------------------------------------------------------
!+
! Add a particle to the clump
!+
!-----------------------------------------------------------------------
subroutine add_particle_to_clump(ipart,iclump)
 use part,     only:xyzh,vxyzu,massoftype,igas
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 integer, intent(in) :: ipart,iclump
 real                :: sep,dx,dy,dz

 member(ipart) = iclump

! Increase membership and mass
 clump(iclump)%num  = clump(iclump)%num + 1
 clump(iclump)%mass = clump(iclump)%mass+ massoftype(igas)

! Update clump energy totals
 clump(iclump)%kinetic = clump(iclump)%kinetic + &
     0.5*massoftype(igas)*((vxyzu(1,ipart)-clump(iclump)%v(1))**2 + &
     (vxyzu(2,ipart)-clump(iclump)%v(2))**2 + &
     (vxyzu(3,ipart)-clump(iclump)%v(3))**2)

 clump(iclump)%potential = clump(iclump)%potential + dpoten(ipart)
 clump(iclump)%thermal   = clump(iclump)%thermal   + vxyzu(4,ipart)*massoftype(igas)

! Calculate distance from the local maximum of the potential
 dx = clump(iclump)%r(1) - xyzh(1,ipart)
 dy = clump(iclump)%r(2) - xyzh(2,ipart)
 dz = clump(iclump)%r(3) - xyzh(3,ipart)
#ifdef PERIODIC
 if (abs(dx) > 0.5*dxbound) then
    dx = dx - dxbound*SIGN(1.0,dx)
    clump(iclump)%border = 1
 endif
 if (abs(dy) > 0.5*dybound) then
    dy = dy - dybound*SIGN(1.0,dy)
    clump(iclump)%border = 1
 endif
 if (abs(dz) > 0.5*dzbound) then
    dz = dz - dzbound*SIGN(1.0,dz)
    clump(iclump)%border = 1
 endif
#endif
 sep = sqrt( dx*dx + dy*dy + dz*dz )

! Update maximum size of the clump if necessary
 if (sep > clump(iclump)%size) clump(iclump)%size = sep

end subroutine add_particle_to_clump
!-----------------------------------------------------------------------
!+
! Remove a particle from the clump
!+
!-----------------------------------------------------------------------
subroutine remove_particle_from_clump(ipart,iclump)
 use part, only:vxyzu,massoftype,igas
 integer, intent(in) :: ipart,iclump

 if (member(ipart)/=iclump) then
    print*, 'Particle ',ipart, ' not in clump ', iclump, ': not removing ', member(ipart)
    return
 endif

 member(ipart) = 0
 clump(iclump)%num  = clump(iclump)%num - 1
 clump(iclump)%mass = clump(iclump)%mass - massoftype(igas)

 ! Update energies
 clump(iclump)%kinetic = clump(iclump)%kinetic - &
       0.5*massoftype(igas)*((vxyzu(1,ipart)-clump(iclump)%v(1))**2 + &
       (vxyzu(2,ipart)-clump(iclump)%v(2))**2 + &
       (vxyzu(3,ipart)-clump(iclump)%v(3))**2)

 clump(iclump)%thermal   = clump(iclump)%thermal   - vxyzu(4,ipart)*massoftype(igas)
 clump(iclump)%potential = clump(iclump)%potential - dpoten(ipart)

end subroutine remove_particle_from_clump
!-----------------------------------------------------------------------
!+
! Test clumps for gravitational boundness
!+
!-----------------------------------------------------------------------
subroutine test_clump_boundness(deletedclumps,npart,xyzh,pmass)
 use part,      only: xyzmh_ptmass,ihacc
 use sortutils, only: indexx
#ifdef PERIODIC
 use boundary,  only:dxbound,dybound,dzbound
#endif
 integer, intent(inout) :: deletedclumps
 integer, intent(in)    :: npart
 real,    intent(in)    :: pmass,xyzh(:,:)
 integer, dimension(:), allocatable :: iseparr,iorig
 real,    dimension(:), allocatable :: separr
 integer                :: iclump,ipart,jpart,counter,originalnum
 real                   :: sep,dx,dy,dz

 print "(/,65('-')/)"
 print*, 'Checking clumps for gravitational boundness'

 ! Loop over clumps
 over_clumps: do iclump=1,nclump

    ! If a clump crosses the boundary, then recalculate the potentials of the member particles
    if (clump(iclump)%border==1) then
       call update_border_potential(iclump,npart,xyzh,xyzmh_ptmass,pmass)
       call calc_clump_boundness(iclump)
    endif

    ! Ignore clumps that are bound
    if (clump(iclump)%bound < 1) cycle over_clumps
    if (clump(iclump)%num == 1 .and. clump(iclump)%pointmass /= 0 ) cycle over_clumps

    originalnum = clump(iclump)%num

    ! If clump unbound, then calculate separation of all particles from centre
    allocate(separr(originalnum))
    allocate(iseparr(originalnum))
    allocate(iorig(originalnum))

    counter = 0
    separations: do ipart=1,npart
       if (member(ipart)/=iclump) cycle separations

       counter          = counter + 1
       iorig(counter)   = ipart
       iseparr(counter) = counter

       dx = clump(iclump)%r(1) - xyzh(1,ipart)
       dy = clump(iclump)%r(2) - xyzh(2,ipart)
       dz = clump(iclump)%r(3) - xyzh(3,ipart)
#ifdef PERIODIC
       if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
       if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
       if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
       sep = sqrt( dx*dx + dy*dy + dz*dz )
       separr(counter) = sep
    enddo separations

    if (counter /= originalnum) then
       print*, "Not all particles listed; likely a result of sink particles.  Resetting value."
       originalnum = counter
    endif

    ! Sort by separation
    call indexx(originalnum,separr,iseparr)

    ! Remove most distant particles until bound

    remove_particles: do jpart=1,originalnum
       ipart = iorig(iseparr(originalnum - jpart+1))

       call remove_particle_from_clump(ipart,iclump)
       call calc_clump_boundness(iclump)

       if (clump(iclump)%bound < 1.0) exit

    enddo remove_particles

    ! Re-add sink particles if they are all that is left
    if (clump(iclump)%num == 0 .and. clump(iclump)%pointmass /= 0 ) clump(iclump)%num = 1

    ! If all particles removed, clump not written to file - log number of deleted clumps
    if (clump(iclump)%num == 0) then
       print*, 'Clump ',iclump,' totally unbound: deleting'
       clump(iclump)%ID = -1
       deletedclumps    = deletedclumps + 1
    else
       ! Update the clump size since we deleted particles
       clump(iclump)%size = 0.
       if (clump(iclump)%pointmass/=0) then
          clump(iclump)%size = sinkclumprad*xyzmh_ptmass(ihacc,clump(iclump)%pointmass)
       endif

       separations_again: do ipart=1,npart
          if (member(ipart)/=iclump) cycle separations_again

          dx = clump(iclump)%r(1) - xyzh(1,ipart)
          dy = clump(iclump)%r(2) - xyzh(2,ipart)
          dz = clump(iclump)%r(3) - xyzh(3,ipart)
#ifdef PERIODIC
          if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
          if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
          if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
          sep = sqrt( dx*dx + dy*dy + dz*dz )
          clump(iclump)%size = max(clump(iclump)%size,sep)
       enddo separations_again

    endif

    ! Deallocate arrays for next clump
    if (allocated(separr))  deallocate(separr)
    if (allocated(iseparr)) deallocate(iseparr)
    if (allocated(iorig))   deallocate(iorig)

 enddo over_clumps

 print*, 'Boundness check complete: ', deletedclumps, ' clumps were deleted'
 print*, 'Remaining clumps: '

 call print_clump_data

end subroutine test_clump_boundness
!-----------------------------------------------------------------------
!+
! Calculate clump's boundness and virial parameter
!+
!-----------------------------------------------------------------------
subroutine calc_clump_boundness(iclump)
 integer,intent(in) :: iclump

 if ( clump(iclump)%potential > 0. ) then
    clump(iclump)%bound  = (clump(iclump)%thermal + clump(iclump)%kinetic) / clump(iclump)%potential
    clump(iclump)%virial = 2.0*clump(iclump)%kinetic/clump(iclump)%potential
 endif

end subroutine calc_clump_boundness
!-----------------------------------------------------------------------
!+
! Recalculate a clump's gravitational potential if it crosses the boundary
! This may not be the most accurate, but for speed and simplicity,
! re-calculate poten(i) via brute force using only the clump members
!+
!-----------------------------------------------------------------------
subroutine update_border_potential(iclump,npart,xyzh,xyzmh_ptmass,pmass)
#ifdef PERIODIC
 use boundary,  only:dxbound,dybound,dzbound
#endif
 integer, intent(in) :: iclump,npart
 real,    intent(in) :: pmass,xyzh(:,:),xyzmh_ptmass(:,:)
 integer             :: i,j,ii,jj,jpt,cID,nmember,imember(maxp)
 real                :: sep,dx,dy,dz,clumppoten,dpotenj_thread(maxp)

 cID        = clump(iclump)%ID
 jpt        = clump(iclump)%pointmass
 clumppoten = 0.
 nmember    = 0
 do i = 1,npart
    if (member(i)==cID) then
       nmember = nmember + 1
       imember(nmember) = i
    endif
 enddo
 print*, 'Recalculating potential for clump ',cID,' which has ',nmember,' members'
!$omp parallel default(none) &
!$omp shared(nmember,imember,dpoten,xyzh,xyzmh_ptmass,dxbound,dybound,dzbound,pmass,jpt) &
!$omp private(i,j,ii,jj,dx,dy,dz,sep,dpotenj_thread) &
!$omp reduction(+:clumppoten)
 dpotenj_thread = 0.
!$omp do schedule(runtime)
 do ii = 1,nmember
    i = imember(ii)
    dpoten(i) = 0.
    !--Add the potential from all members
    do jj = ii+1,nmember
       j = imember(jj)
       dx = xyzh(1,i) - xyzh(1,j)
       dy = xyzh(2,i) - xyzh(2,j)
       dz = xyzh(3,i) - xyzh(3,j)
#ifdef PERIODIC
       if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
       if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
       if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
       sep = sqrt( dx**2 + dy**2 + dz**2 )
       dpoten(i) = dpoten(i) + pmass**2/sep
       dpotenj_thread(j) = dpotenj_thread(j) + pmass**2/sep
    enddo
    !--Add the potential from the included sink
    if (jpt /= 0) then
       dx = xyzh(1,i) - xyzmh_ptmass(1,jpt)
       dy = xyzh(2,i) - xyzmh_ptmass(2,jpt)
       dz = xyzh(3,i) - xyzmh_ptmass(3,jpt)
#ifdef PERIODIC
       if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
       if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
       if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
       sep        = sqrt( dx**2 + dy**2 + dz**2 )
       dpoten(i)  = dpoten(i)  + pmass*xyzmh_ptmass(4,jpt)/sep
    endif
 enddo
!$omp enddo
!$omp critical(collatedata)
 do ii = 1,nmember
    i = imember(ii)
    dpoten(i) = dpoten(i) + dpotenj_thread(i)
 enddo
!$omp end critical(collatedata)
!$omp end parallel

!$omp parallel default(none) &
!$omp shared(nmember,imember,dpoten) &
!$omp private(i,ii) &
!$omp reduction(+:clumppoten)
!$omp do schedule(runtime)
 do ii = 1,nmember
    i = imember(ii)
    clumppoten = clumppoten + dpoten(i)
 enddo!$omp enddo
!$omp end parallel
 clump(iclump)%potential = clumppoten

end subroutine update_border_potential
!-----------------------------------------------------------------------
!+
! Print clump data to the screen
!+
!-----------------------------------------------------------------------
subroutine print_clump_data
 integer :: iclump
 character(len=100) :: fmt
 character(len=260) :: title

 fmt   = '(I4,X,I7,X, 10(1pg10.3,X))'
 title = 'ID    Npart     x        y        z       mass      size      '
 title = trim(title)//'   KE       PE      TE      Boundness      Virial Parameter'

 print '(a)', trim(title)

 do iclump=1,nclump
    if (clump(iclump)%num==0) cycle
    print fmt, clump(iclump)%ID, clump(iclump)%num, clump(iclump)%r(1), clump(iclump)%r(2), clump(iclump)%r(3), &
          clump(iclump)%mass, clump(iclump)%size, clump(iclump)%kinetic, clump(iclump)%potential, clump(iclump)%thermal, &
          clump(iclump)%bound, clump(iclump)%virial
 enddo

end subroutine print_clump_data
!-----------------------------------------------------------------------
!+
! Write the clump data to file
!+
!-----------------------------------------------------------------------
subroutine write_clump_data(nclump,deletedclumps,npart,time,dumpfile,tag)

 character(len=*), intent(in) :: dumpfile
 character(len=*), intent(in) :: tag
 real,             intent(in) :: time
 integer,          intent(in) :: nclump,npart,deletedclumps

 integer                      :: i,j,iclump
 character(len=100)           :: clumpfile

! Write catalogue of clump data

! These are tagged 'raw' if clump IDs not checked between dumps
! These are tagged 'proc' if files have been processed to check clump IDs

 clumpfile = trim(tag)//"_clumpcat_"//trim(dumpfile)

 open(10,file=clumpfile, status='unknown')
 write(10,*) nclump-deletedclumps, time

 do iclump=1,nclump
    if (clump(iclump)%ID > 0) then
       write(clumpfile,'(2a,I5.5)') trim(tag),"_clumpcat_clumpID",clump(iclump)%ID
       open(11,file=clumpfile, position='append')
       do j = 10,11
          write(j,*) clump(iclump)%ID, clump(iclump)%r(:),clump(iclump)%v(:),clump(iclump)%mass, &  !  1,2-4,5-7,8
                     clump(iclump)%size, clump(iclump)%num, &                                       !  9,10  (radius & number of particles)
                     clump(iclump)%mbp, clump(iclump)%pointmass, clump(iclump)%kinetic, &           ! 11,12,13
                     clump(iclump)%potential, clump(iclump)%thermal, &                              ! 14,15
                     clump(iclump)%bound, clump(iclump)%virial,time                                 ! 16,17,18
       enddo
       close(11)
    endif
 enddo
 close(10)

! Write complete list of particle clump membership to file
 clumpfile = trim(tag)//"_clumpmembers_"//trim(dumpfile)

 print*, 'Writing individual particle memberships to ',clumpfile
 do i = 1,npart
    if (member(i) > 0) member(i) = clump(member(i))%ID
 enddo

 open(10,file=clumpfile, form='unformatted')
 write(10) (member(i), i=1,npart)
 close(10)

end subroutine write_clump_data
!-----------------------------------------------------------------------
!+
! Read clump data from previous timestep
!+
!-----------------------------------------------------------------------
subroutine read_oldclump_data(noldclump,npart,oldtime,olddumpfile,tag)
 character(len=*), intent(in)    :: olddumpfile
 character(len=*), intent(in)    :: tag
 integer,          intent(in)    :: npart
 real,             intent(inout) :: oldtime
 integer,          intent(inout) :: noldclump
 integer                         :: i,iclump
 real                            :: time
 character(len=100)              :: clumpfile

! Read catalogue of clump data

! These are tagged 'raw' if clump IDs not checked between dumps
! These are tagged 'proc' if files have been processed to check clump IDs

 clumpfile = trim(tag)//"_clumpcat_"//trim(olddumpfile)

 open(10,file=clumpfile, status='unknown')
 read(10,*) noldclump, oldtime

 allocate(oldclump(noldclump))

 do iclump=1,noldclump
    read(10,*) oldclump(iclump)%ID,        oldclump(iclump)%r(:),      oldclump(iclump)%v(:),    &
               oldclump(iclump)%mass,      oldclump(iclump)%size,      oldclump(iclump)%num,     &
               oldclump(iclump)%mbp,       oldclump(iclump)%pointmass, oldclump(iclump)%kinetic, &
               oldclump(iclump)%potential, oldclump(iclump)%thermal,   oldclump(iclump)%bound,   &
               oldclump(iclump)%virial,    time
 enddo
 close(10)

! Read complete list of particle clump membership from file
 clumpfile = trim(tag)//"_clumpmembers_"//trim(olddumpfile)

 print*, 'Reading individual particle memberships from ',clumpfile

 allocate(oldmember(npart))

 open(10,file=clumpfile, form='unformatted')
 read(10) (oldmember(i), i=1,npart)
 close(10)

end subroutine read_oldclump_data
!-----------------------------------------------------------------------
!+
! Do merger tree calculation to link clumps in current dump to previous dump
!+
!-----------------------------------------------------------------------
subroutine merger_tree(npart,dumpfile)
 integer,         intent(in) :: npart
 character(len=*),intent(in) :: dumpfile
 integer                     :: IDmax,noldclump,ipart,iclump,jclump
 real                        :: memberfraction,oldtime

 print*, 'This is not the first dump - attempting merger tree calculation'
 print*, 'Reading data from previous dump ', previousdumpfile

 ! Read previous clump data
 call read_oldclump_data(noldclump,npart,oldtime,previousdumpfile,proctag)

 print '(I5,a,a)', noldclump, ' clumps found in previous dump ', previousdumpfile

 ! Clear IDs of clumps in this dump
 do jclump=1,nclump
    clump(jclump)%ID = 0
 enddo

 !*********************************************************
 ! Procedure:
 ! Find where previous clump's most bound particles (MBPs) are in current dump
 ! Check 50% of previous clump membership is also in the new clump containing MBP
 ! Assign new IDs
 !*********************************************************

 print*, 'Running merger tree on current and previous dumps'

 ! Loop over previous file's clumps
 do iclump = 1,noldclump
    if (oldclump(iclump)%pointmass /= 0) then
       ! Use sinks to match clumps
       do jclump = 1,nclump
          if ( oldclump(iclump)%pointmass == clump(jclump)%pointmass ) then
             clump(jclump)%ID = oldclump(iclump)%ID
          endif
       enddo
    else
       ! Use most bound particles to match clumps
       jclump = member(oldclump(iclump)%mbp)
       if (jclump==0) cycle ! What if iclump not found in new dump? TODO
       if (clump(jclump)%num==0 .or. clump(jclump)%pointmass > 0) cycle ! What if iclump not found in new dump? TODO

       ! Determine the amount of membership overlap between the two clumps
       memberfraction = 0.0
       do ipart=1,npart
          if (oldmember(ipart)==oldclump(iclump)%ID .and. member(ipart)==jclump) then
             memberfraction = memberfraction + 1.0
          endif
       enddo

       if (oldclump(iclump)%num > 0.0) then
          memberfraction = memberfraction/real(oldclump(iclump)%num)
       else
          memberfraction = 0.0
       endif

       ! If membership overlap > 0.5, then tag the two clumps as the same
       if (memberfraction > 0.5) then
          clump(jclump)%ID = oldclump(iclump)%ID
       endif
       print*, 'old clump id =',iclump,'. new/old clump membership overlap = ',memberfraction,'. N_old =',oldclump(iclump)%num
    endif
 enddo

 print*, 'Checking new clump IDs'

! Now check for clumps that have no progenitor in the previous timestep
 IDmax = maxval(clump(:)%ID)
 if (IDmax > runningclumpmax) runningclumpmax = IDmax

 do jclump=1,nclump
    if (clump(jclump)%ID==0 .and. clump(jclump)%num > 0) then
       runningclumpmax  = runningclumpmax + 1
       clump(jclump)%ID = runningclumpmax
    endif
 enddo

! Amend .options file in preparation for next dump
 call amend_options_file(dumpfile)

 print*, 'Merger Tree Done'
 call print_clump_data

! Deallocate arrays for old file
 deallocate(oldclump,oldmember)

end subroutine merger_tree
!-----------------------------------------------------------------------
!+
! Subroutine ensures that maximum clump number not exceeded
!+
!-----------------------------------------------------------------------
subroutine check_clump_max

 if (nclump > nclumpmax) then
    print '(a)', 'ERROR: total clump number exceeds maximum clump number:'
    print '(a,i6,1x,i6)', 'nclump, nclumpmax: ', nclump, nclumpmax
    print '(a)', 'Please recompile with higher nclumpmax'
    stop
 endif

end subroutine check_clump_max
!-----------------------------------------------------------------------
end module analysis
