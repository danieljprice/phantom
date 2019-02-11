!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine which runs the CLUMPFIND algorithm
!
!  This uses local potential minima (and neighbour lists)
!  to identify objects/clumps in the gas distribution
!  Clumps are matched across timesteps using standard merger tree algorithms
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, kernel, linklist, part, prompting, ptmass, sortutils
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'clumpfind'

 real, parameter :: sinkclumprad = 10.0

 integer, parameter :: neighmax = 1000
 integer, parameter :: maxcellcache = 50000
 integer, parameter :: nclumpmax = 1000

 integer, allocatable, dimension(:) :: neighcount
 integer, allocatable, dimension(:,:) :: neighb

 integer :: runningclumpmax
 real :: meanneigh, sdneigh, neighcrit
 logical :: neigh_overload, checkbound, sinkpotential,skipsmalldumps, firstdump
 character(100) :: firstdumpfile, previousdumpfile, rawtag, proctag

 public :: do_analysis

 ! Derived type sphclump records clump data
 ! (Total particle counts, clump position/velocity, mass and extent)

 type sphclump

    integer :: num, mbp, pointmass, ID
    real,dimension(3) :: r, v
    real :: mass, size
    real :: kinetic, potential, thermal, bound, virial

 end type sphclump

 type(sphclump), allocatable, dimension(:) :: clump, oldclump
 integer, allocatable, dimension(:) :: member, oldmember
 real, allocatable, dimension(:) :: dpoten
 integer :: nclump

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 use dim, only: maxp
 use part, only: iphase,maxphase, igas, get_partinfo,ihacc, poten
 use sortutils, only: indexx

 implicit none

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time

 integer :: iamtypei

 integer, parameter :: neighcrit = 40

 integer, dimension(1) :: maxclump
 integer, allocatable,dimension(:) :: neighclump, ipotensort

 integer :: iclump,ipart,jpart
 integer :: k,l, deletedclumps
 real :: percent,percentcount

 logical :: iactivei, iamdusti,iamgasi
 logical :: existneigh

 character(100) :: neighbourfile
 character(100) :: fmt

 ! Read in input parameters from file (if it exists)

 rawtag= 'raw'
 proctag = 'proc'

 checkbound = .false.
 sinkpotential = .false.
 skipsmalldumps = .false.
 firstdump = .false.

 call read_analysis_options(dumpfile)

 if(skipsmalldumps) then

    ! Skip small dumps (as they do not include velocity data)
    if(sum(vxyzu(1,:)) < tiny(vxyzu) .and. &
        sum(vxyzu(2,:))< tiny(vxyzu) .and. &
        sum(vxyzu(3,:)) < tiny(vxyzu)) then
       print*, 'Skipping clumpfind analysis for small dumps (velocity data missing)'
       return
    endif
 endif

 !***************************************
 !1. Obtain Neighbour Lists
 !****************************************

 print '(a)', 'Obtaining Neighbour Lists'

 ! Check if a neighbour file is present

 neighbourfile = 'neigh_'//TRIM(dumpfile)
 inquire(file=neighbourfile,exist = existneigh)

 if(existneigh.eqv..true.) then
    print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
    call read_neighbours(neighbourfile,npart)

 else

    ! If there is no neighbour file, generate the list
    print*, 'No neighbour file found: generating'

    call generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile)

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

 if(sinkpotential) call add_sink_potential(npart)

 dpoten = ABS(dpoten)

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

 allocate(clump(nclumpmax)) ! Maximum number of allowed clumps = npart

 print*, 'Clump data allocated'

 allocate(member(npart)) ! Clump membership for each particle


 allocate(neighclump(nclumpmax)) ! For a set of particle neighbours, this tracks which clumps they are in

 print*, 'Catalogues allocated'

 nclump = 0 ! Records total number of clumps counted so far

! Set membership for all particles to zero
 member(:) = 0

! Set positions and velocities of clumps to zero
 do k=1,3
    clump(:)%r(k) = 0.0
    clump(:)%v(k) = 0.0
 enddo

! Set mass, maximum radius and total number of particles per clump to zero
 clump(:)%ID = 0
 clump(:)%mass = 0.0
 clump(:)%size = 0.0
 clump(:)%num = 0
 clump(:)%kinetic = 0.0
 clump(:)%potential = 0.0
 clump(:)%thermal = 0.0

 clump(:)%pointmass=0 ! Indicates if clump contains a pointmass (and which)

!
!--All sinks are automatically identified as clumps

!

 call create_sink_clumps(npart, xyzh)

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

    if(percent>percentcount) THEN
       write(*,'(I3," % complete")') INT(percentcount)
       percentcount = percentcount +10.0
    ENDIF

    ! Pick next particle in order of descending potential

    ipart = ipotensort(npart-l+1)

    ! Only calculate for gas particles

    if(maxphase==maxp) then
       call get_partinfo(iphase(ipart), iactivei,iamdusti,iamtypei)
       iamgasi = (iamtypei ==igas)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi = .true.
    endif

    ! If particle isn't a gas particle, skip it
    if(.not.iamgasi) cycle over_parts

    ! If particle not already in a clump
    IF(member(ipart)==0) THEN
       !print*, 'Checking neighbours'

       ! Reset array for counting neighbour particle clump membership
       neighclump(:) = 0

       !        Check if its neighbours are in a clump
       do jpart=1,neighcount(ipart)
          if(neighcount(ipart) > neighmax) cycle
          iclump = member(neighb(ipart,jpart))

          IF(iclump>0) then
             neighclump(iclump) = neighclump(iclump)+1
          endif

       enddo

       ! Compare neighbour membership, identify clump with most neighbours as members
       ! This is now particle ipart's host clump

       IF(sum(neighclump)/=0) THEN
          maxclump = maxloc(neighclump)
          call add_particle_to_clump(ipart,maxclump(1))
       ENDIF

    ENDIF

    ! If particle still does not belong to any other clump, then create a new clump for it
    IF(member(ipart)==0.and.neighcount(ipart)>neighcrit) call initialise_clump(ipart)

! If particle is already in a clump, then add any neighbours not in a clump
    IF(member(ipart)>0) THEN

       do jpart=1,neighcount(ipart)
          if(neighcount(ipart) > neighmax) cycle
          IF(member(neighb(ipart,jpart))==0) call add_particle_to_clump(neighb(ipart,jpart), member(ipart))
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

    call test_clump_boundness(deletedclumps,npart,xyzh)


 endif

!****************************************************
!
!--5. Write raw (unmatched) clump data to file
!
!****************************************************
 print*, 'Clumpfind complete for snapshot ', dumpfile
 print*, 'Writing unmatched clump data to raw files'

 call write_clump_data(nclump,deletedclumps,npart, time,dumpfile,rawtag)

!****************************************************
!
!--6. Link clumps in this dump with the previous dump
!
!****************************************************

! Check if this is the first file to be analysed
! If so, then we'll simply write the data to file
! as these are the original clumps

 if (.not.firstdump) then

! Otherwise do the merger tree analysis
    call merger_tree(npart,dumpfile)

 else
    ! If this is the first dump, update current maximum clump count and save it
    runningclumpmax = nclump
    call amend_options_file(dumpfile)
 endif

! Now write to file

 call write_clump_data(nclump,deletedclumps,npart, time,dumpfile,proctag)

! Deallocate memory for next dump
 deallocate(neighcount,neighb,dpoten)
 deallocate(member,neighclump,clump)

end subroutine do_analysis

!---------------------------------
!+
! Read in options for analysis from file
!+
!---------------------------------

subroutine read_analysis_options(dumpfile)

 use prompting, only:prompt

 character(len=*), intent(in) :: dumpfile
 logical :: inputexist
 character(len=17) :: inputfile
 character(len=1) :: boundchoice, sinkchoice,skipchoice

! Check for existence of input file
 inputfile = 'clumpfind.options'
 inquire(file=inputfile, exist=inputexist)

 if(inputexist) then

    print '(a,a,a)', "File ",inputfile, " found: reading analysis options"

    open(10,file=inputfile, form='formatted')
    read(10,*) boundchoice
    read(10,*) sinkchoice
    read(10,*) skipchoice
    read(10,*) firstdumpfile



! If this is a file from a previous run, then we need to be careful
! about what the previous dump was

    if(trim(dumpfile)==trim(firstdumpfile)) then
       write(10,*) trim(dumpfile), "    Previous SPH dump analysed"
       runningclumpmax = 0
    else
       read(10,*) previousdumpfile
       read(10,*) runningclumpmax
    endif

    close(10)

    if(boundchoice=='y' .or. boundchoice=='Y') checkbound=.true.
    if(sinkchoice=='y' .or. sinkchoice=='Y') sinkpotential=.true.
    if(skipchoice=='y' .or. skipchoice=='Y') skipsmalldumps=.true.
    if(trim(firstdumpfile)==trim(dumpfile)) firstdump = .true.

 else

    boundchoice = 'n'
    sinkchoice = 'n'
    skipchoice = 'n'

    call prompt('Do you want to test for boundness?', checkbound, .true.)
    call prompt('Do you want to include sink contribution to the potential?', sinkpotential,.true.)
    call prompt('Do you want to skip small dumps (missing velocity data)? ',skipsmalldumps,.true.)

    if(checkbound) boundchoice = 'y'
    if(sinkpotential) sinkchoice = 'y'
    if(skipsmalldumps) skipchoice = 'y'

! Write choices to new inputfile

    open(10,file=inputfile, status='new', form='formatted')
    write(10,*) boundchoice, "      Test clumps for gravitational boundness?"
    write(10,*) sinkchoice,  "      Include sinks' contribution to potential?"
    write(10,*) skipchoice,  "      Skip small dumps (velocity data missing)?"
    write(10,*) trim(dumpfile), "    First SPH dump to be analysed"
    write(10,*) trim(dumpfile), "    Previous SPH dump analysed"
    close(10)
 endif

 if(checkbound) then
    print '(a)', 'Clumps will be checked for gravitational boundness'
 else
    print '(a)', 'Clumps will NOT be checked for gravitational boundness'
 endif

 if(sinkpotential) then
    print '(a)', 'Contribution of sinks to the potential will be added'
 else
    print '(a)', 'Contribution of sinks to the potential will NOT be added'
 endif

 if(skipsmalldumps) then
    print '(a)', 'Skipping small dumps'
 else
    print '(a)', 'NOT skipping small dumps'
 endif

 if(firstdump) then
    print '(a)', 'This is the first dump to be analysed'
 else
    print '(a)', 'This is not the first dump - will attempt to link clump IDs to previous dump'
 endif

 if(checkbound .and. .not.(skipsmalldumps)) then
    print*, 'ERROR in read_analysis_options!'
    print*, 'Cannot check boundness using small dumps!'
    print*, 'Please rerun, either choosing to check boundness or to use small dumps, not both'
    stop
 endif

end subroutine read_analysis_options

!------------------------------------------
!+
! Amend the options file to update which dump was previously analysed
!+
!------------------------------------------
subroutine amend_options_file(dumpfile)

 implicit none
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


!----------------------------------
!+
! Add contribution to potential from sinks
!+
!----------------------------------

subroutine add_sink_potential(npart)

 use part,    only: xyzh,massoftype,igas,xyzmh_ptmass, nptmass
 use ptmass,  only:get_accel_sink_gas

 integer, intent(in) :: npart
 integer :: ipart

 real :: dumx,dumy,dumz, epoti

 print '(a)', 'Calculating sink contribution to the gravitational potential'

 do ipart=1,npart

! Dummy variables for gravitational forces (not used here)
    dumx = 0.0
    dumy = 0.0
    dumz = 0.0

    epoti = 0.0

    call get_accel_sink_gas(nptmass,xyzh(1,ipart),xyzh(2,ipart),xyzh(3,ipart),xyzh(4,ipart),&
xyzmh_ptmass,dumx,dumy,dumz,epoti)

    dpoten(ipart) = dpoten(ipart) + massoftype(igas)*epoti

 enddo

 print '(a)', 'Done'

end subroutine add_sink_potential

!---------------------------------------------------------------
!+
! Create clumps at the location of all sinks
!--Sinks do not possess neighbours in Phantom
!--All particles within a given radius of the sink are added to the clump
!+
!----------------------------------------------------------------

subroutine create_sink_clumps(npart, xyzh)
 use part, only: xyzmh_ptmass, vxyz_ptmass, nptmass, ihacc
 implicit none

 integer, intent(in) :: npart
 real, intent(in) :: xyzh(:,:)

 integer :: iclump, ipart, jpart,iptmass,k
 real :: clumpsep,sep
 logical :: first
 character(100) :: fmt

 print "(/,65('-')/)"

 print*, 'Automatically identifying sink particles as progenitors of clumps'
 print'(a,F5.1,a)', 'All particles within ',sinkclumprad,' accretion radii  will be added to sink clumps'

 do iptmass=1,nptmass

    clumpsep = sinkclumprad*xyzmh_ptmass(ihacc,iptmass)

    nclump = nclump +1

    call check_clump_max

    clump(iptmass)%ID = iptmass
    clump(iptmass)%pointmass = iptmass

    do k=1,3
       clump(iptmass)%r(k) = xyzmh_ptmass(k,iptmass)
       clump(iptmass)%v(k) = vxyz_ptmass(k,iptmass)
    enddo
    clump(iptmass)%mass = xyzmh_ptmass(4,iptmass)

    print '(a,I4,a)', 'Pointmass ',iptmass, ' begins a clump:'
    print '(1pg10.3,2X,1pg10.3,2X,1pg10.3,2X,1pg10.3)', clump(nclump)%r(1), &
        clump(nclump)%r(2), clump(nclump)%r(3), clump(nclump)%mass

    ! Test each particle for distance

    do ipart=1,npart

       sep = sqrt((xyzh(1,ipart) - xyzmh_ptmass(1,iptmass))**2 + &
           (xyzh(2,ipart) - xyzmh_ptmass(2,iptmass))**2 + &
           (xyzh(3,ipart) - xyzmh_ptmass(3,iptmass))**2)

       if(sep < clumpsep .and.member(ipart)==0) then

          call add_particle_to_clump(ipart,iptmass)

          if(first) then
             clump(iptmass)%mbp = ipart
             first = .false.
          endif

          ! Add all SPH particle's neighbours to this clump

          do jpart=1,neighcount(ipart)
             IF(member(neighb(ipart,jpart))==0) call add_particle_to_clump(neighb(ipart,jpart), iptmass)
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

!----------------------------------
!+
! Create a new clump at the location of SPH particle ipart
!+
!----------------------------------
subroutine initialise_clump(ipart)

 use part, only: xyzh, vxyzu, massoftype,igas

 integer, intent(in) :: ipart
 character(len=100) :: fmt

 integer :: k

 fmt = "(a,I4,a,I8,2X,1pg10.3,2X,1pg10.3,2x,1pg10.3)"

 nclump = nclump +1
 call check_clump_max

 clump(nclump)%ID = nclump
 member(ipart) = nclump
 clump(nclump)%num = clump(nclump)%num + 1

 print fmt, 'Clump ', nclump, ' created for particle ', ipart,&
     xyzh(1,ipart),xyzh(2,ipart),xyzh(3,ipart)

! Store most bound particle index here
 clump(nclump)%mbp  = ipart

! Locate origin of clump (max of potential), and initialise its mass
 clump(nclump)%mass = massoftype(igas)

 do k=1,3
    clump(nclump)%r(k) = xyzh(k,ipart)
    clump(nclump)%v(k) = vxyzu(k,ipart)
 enddo

 clump(nclump)%thermal = vxyzu(4,ipart)*massoftype(igas)
 clump(nclump)%potential = dpoten(ipart)

 clump(nclump)%size = 0.0

end subroutine initialise_clump

!-------------------------------
!+
! Add a particle to the clump
!+
!------------------------------
subroutine add_particle_to_clump(ipart,iclump)
 use part, only: xyzh, vxyzu, massoftype,igas
 integer, intent(in) :: ipart,iclump

 integer :: k
 real :: sep

 member(ipart) = iclump

! Increase membership and mass

 clump(iclump)%num = clump(iclump)%num + 1
 clump(iclump)%mass = clump(iclump)%mass+ massoftype(igas)

! Update clump energy totals
 clump(iclump)%kinetic = clump(iclump)%kinetic + &
     0.5*massoftype(igas)*((vxyzu(1,ipart)-clump(iclump)%v(1))**2 + &
     (vxyzu(2,ipart)-clump(iclump)%v(2))**2 + &
     (vxyzu(3,ipart)-clump(iclump)%v(3))**2)

 clump(iclump)%potential = clump(iclump)%potential + dpoten(ipart)
 clump(iclump)%thermal = clump(iclump)%thermal + vxyzu(4,ipart)*massoftype(igas)

! Calculate distance from the local maximum of the potential
 sep = 0.0
 do k=1,3
    sep = sep + (clump(iclump)%r(k) - xyzh(k,ipart))**2
 enddo

 sep = sqrt(sep)

! Update maximum size of the clump if necessary
 IF(sep > clump(iclump)%size) clump(iclump)%size = sep

end subroutine add_particle_to_clump

subroutine remove_particle_from_clump(ipart,iclump)
 use part, only: vxyzu, massoftype,igas
 integer, intent(in) :: ipart,iclump

 if(member(ipart)/=iclump) then
    print*, 'Particle ',ipart, ' not in clump ', iclump, ': not removing ', member(ipart)
    return
 endif

 member(ipart) = 0
 clump(iclump)%num = clump(iclump)%num - 1
 clump(iclump)%mass = clump(iclump)%mass - massoftype(igas)

 ! Update energies
 clump(iclump)%kinetic = clump(iclump)%kinetic - &
       0.5*massoftype(igas)*((vxyzu(1,ipart)-clump(iclump)%v(1))**2 + &
       (vxyzu(2,ipart)-clump(iclump)%v(2))**2 + &
       (vxyzu(3,ipart)-clump(iclump)%v(3))**2)

 clump(iclump)%thermal = clump(iclump)%thermal - vxyzu(4,ipart)*massoftype(igas)
 clump(iclump)%potential = clump(iclump)%potential - dpoten(ipart)

end subroutine remove_particle_from_clump

!--------------------------------------------------------
!+
! Test clumps for gravitational boundness
!+
!---------------------------------------------------------

subroutine test_clump_boundness(deletedclumps,npart,xyzh)

 use sortutils, only: indexx

 implicit none

 integer, intent(inout) :: deletedclumps
 integer, intent(in) :: npart
 real, intent(in) :: xyzh(:,:)

 integer, dimension(:), allocatable :: iseparr,iorig
 real, dimension(:) , allocatable :: separr

 integer :: iclump, ipart, jpart, k,counter, originalnum
 real :: sep

 print "(/,65('-')/)"
 print*, 'Checking clumps for gravitational boundness'

 ! Loop over clumps
 over_clumps: do iclump=1,nclump

    ! Ignore clumps that are bound

    if(clump(iclump)%bound < 1) cycle over_clumps

    originalnum = clump(iclump)%num

    ! If clump unbound, then calculate separation of all particles from centre
    allocate(separr(originalnum))
    allocate(iseparr(originalnum))
    allocate(iorig(originalnum))

    counter = 0

    separations: do ipart=1,npart

       if(member(ipart)/=iclump) cycle separations

       counter = counter + 1
       iorig(counter) = ipart
       iseparr(counter) = counter

       sep = 0.0
       do k=1,3
          sep = sep + (clump(iclump)%r(k) - xyzh(k,ipart))**2
       enddo

       separr(counter) = sqrt(sep)
    enddo separations

    ! Sort by separation
    call indexx(originalnum,separr,iseparr)

    ! Remove most distant particles until bound

    remove_particles: do jpart=1,originalnum
       ipart = iorig(iseparr(originalnum - jpart+1))

       call remove_particle_from_clump(ipart,iclump)
       call calc_clump_boundness(iclump)

       if(clump(iclump)%bound < 1.0) exit

    enddo remove_particles

    ! If all particles removed, clump not written to file - log number of deleted clumps
    if(clump(iclump)%num ==0) then
       print*, 'Clump ',iclump,' totally unbound: deleting'
       clump(iclump)%ID = -1
       deletedclumps = deletedclumps+1
    endif

    ! Deallocate arrays for next clump
    if(allocated(separr)) deallocate(separr)
    if(allocated(iseparr)) deallocate(iseparr)
    if(allocated(iorig)) deallocate(iorig)

 enddo over_clumps


 print*, 'Boundness check complete: ', deletedclumps, ' clumps were deleted'
 print*, 'Remaining clumps: '

 call print_clump_data

end subroutine test_clump_boundness
!-------------------------------------------
!+
! Calculate clump's boundness and virial parameter
!+
!------------------------------------------

subroutine calc_clump_boundness(iclump)
 integer,intent(in) :: iclump

 clump(iclump)%bound = (clump(iclump)%thermal + clump(iclump)%kinetic)/ clump(iclump)%potential
 clump(iclump)%virial = 2.0*clump(iclump)%kinetic/clump(iclump)%potential

end subroutine calc_clump_boundness


!--------------------------------
!+
! Print clump data to the screen
!+
!--------------------------------

subroutine print_clump_data

 integer :: iclump
 character(len=100) :: fmt
 character(len=260) :: title

 fmt = '(I4,X,I7,X, 10(1pg10.3,X))'
 title = 'ID    Npart     x        y        z       mass      size      '
 title = trim(title)//'   KE       PE      TE      Boundness      Virial Parameter'

 print '(a)', trim(title)

 do iclump=1,nclump

    if(clump(iclump)%num==0) cycle

    print fmt, clump(iclump)%ID, clump(iclump)%num, clump(iclump)%r(1), clump(iclump)%r(2), clump(iclump)%r(3), &
clump(iclump)%mass, clump(iclump)%size, clump(iclump)%kinetic, clump(iclump)%potential, clump(iclump)%thermal, &
clump(iclump)%bound, clump(iclump)%virial

 enddo


end subroutine print_clump_data


!------------------------------
!+
! Write the clump data to file
!+
!------------------------------
subroutine write_clump_data(nclump,deletedclumps,npart, time,dumpfile,tag)

 character(len=*), intent(in) :: dumpfile
 character(len=*), intent(in) :: tag
 real, intent(in) :: time
 integer, intent(in) :: nclump,npart,deletedclumps

 integer :: i,iclump
 character(len=100) :: clumpfile

! Write catalogue of clump data

! These are tagged 'raw' if clump IDs not checked between dumps
! These are tagged 'proc' if files have been processed to check clump IDs

 clumpfile = trim(tag)//"_clumpcat_"//trim(dumpfile)

 OPEN(10,file=clumpfile, status='unknown')
 write(10,*) nclump-deletedclumps, time

 do iclump=1,nclump
    if(clump(iclump)%ID>0) then
       write(10,*) clump(iclump)%ID, clump(iclump)%r(:),clump(iclump)%v(:),clump(iclump)%mass, &
        clump(iclump)%size, clump(iclump)%num, &
        clump(iclump)%mbp, clump(iclump)%pointmass, clump(iclump)%kinetic, &
        clump(iclump)%potential, clump(iclump)%thermal, &
        clump(iclump)%bound, clump(iclump)%virial
    endif
 enddo
 close(10)

! Write complete list of particle clump membership to file
 clumpfile = trim(tag)//"_clumpmembers_"//trim(dumpfile)

 print*, 'Writing individual particle memberships to ',clumpfile

 OPEN(10,file=clumpfile, form='unformatted')
 write(10) (member(i), i=1,npart)
 close(10)

end subroutine write_clump_data


!------------------------------
!+
! Read clump data from previous timestep
!+
!------------------------------
subroutine read_oldclump_data(noldclump,npart, oldtime,olddumpfile, tag)

 character(len=*), intent(in) :: olddumpfile
 character(len=*), intent(in) :: tag

 integer, intent(in) :: npart
 real, intent(inout) :: oldtime
 integer, intent(inout) :: noldclump

 integer :: i,iclump
 character(len=100) :: clumpfile

! Read catalogue of clump data

! These are tagged 'raw' if clump IDs not checked between dumps
! These are tagged 'proc' if files have been processed to check clump IDs

 clumpfile = trim(tag)//"_clumpcat_"//trim(olddumpfile)

 OPEN(10,file=clumpfile, status='unknown')
 read(10,*) noldclump, oldtime

 allocate(oldclump(noldclump))

 do iclump=1,noldclump
    read(10,*) oldclump(iclump)%ID, oldclump(iclump)%r(:),oldclump(iclump)%v(:), &
   oldclump(iclump)%mass, oldclump(iclump)%size, oldclump(iclump)%num, &
   oldclump(iclump)%mbp, oldclump(iclump)%pointmass, oldclump(iclump)%kinetic, &
   oldclump(iclump)%potential, oldclump(iclump)%thermal, &
   oldclump(iclump)%bound, clump(iclump)%virial
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

!---------------------------------------------------------
!+
! Do merger tree calculation to link clumps in current dump to previous dump
!+
!---------------------------------------------------------

subroutine merger_tree(npart,dumpfile)

 implicit none
 integer, intent(in) :: npart
 character(len=*),intent(in) :: dumpfile

 integer :: IDmax, noldclump, ipart, iclump,jclump
 real :: memberfraction, oldtime

 print*, 'This is not the first dump - attempting merger tree calculation'
 print*, 'Reading data from previous dump ', previousdumpfile

 ! Read previous clump data
 call read_oldclump_data(noldclump,npart, oldtime,previousdumpfile, proctag)

 print '(I5,a,a)', noldclump, ' clumps found in previous dump ', previousdumpfile

 ! Clear IDs of clumps in this dump (except for sinks)
 do jclump=1,nclump
    if(clump(jclump)%pointmass==0) clump(jclump)%ID = 0
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

    ! If this clump begun by a sink, find it in the current dump --> automatically linked
    if(oldclump(iclump)%pointmass/=0) then

       ! Otherwise use the MBP
    else

       ! clump in current dump containing MBP from previous dump
       jclump = member(oldclump(iclump)%mbp)

       if(jclump==0 .or. clump(jclump)%num==0 .or.clump(jclump)%pointmass >0) cycle ! What if iclump not found in new dump? TODO

       ! Now check that 50% of the membership of iclump is in jclump
       memberfraction = 0.0

       do ipart=1,npart

          if(oldmember(ipart)/=iclump) cycle
          if(member(ipart)==jclump) then
             memberfraction = memberfraction + 1
          endif
       enddo

       if(oldclump(iclump)%num > 0.0) then
          memberfraction = memberfraction/real(oldclump(iclump)%num)
       else
          memberfraction = 0.0
       endif

       if(memberfraction>0.5) then
          clump(jclump)%ID = oldclump(iclump)%ID
       endif
    endif
 enddo

 print*, 'Checking new clump IDs'

! Now check for clumps that have no progenitor in the previous timestep
 IDmax = maxval(clump(:)%ID)
 if(IDmax > runningclumpmax) runningclumpmax = IDmax

 do jclump=1,nclump
    if(clump(jclump)%ID==0 .and. clump(jclump)%num>0) then
       runningclumpmax = runningclumpmax +1
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


!-------------------------------------------------------
!+
! Subroutine ensures that maximum clump number not exceeded
!+
!-------------------------------------------------------

subroutine check_clump_max

 if(nclump>nclumpmax) then
    print '(a)', 'ERROR: total clump number exceeds maximum clump number:'
    print '(a,i6,1x,i6)', 'nclump, nclumpmax: ', nclump, nclumpmax
    print '(a)', 'Please recompile with higher nclumpmax'
    stop
 endif

end subroutine check_clump_max

!--------------------------------------------------------
!+
! Generate neighbour lists for all particles
!+
!--------------------------------------------------------
subroutine generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile)

 use dim, only: maxneigh,maxp
 use kernel, only: radkern2
 use linklist, only: ncells, ifirstincell, set_linklist, get_neighbour_list
 use part, only: get_partinfo, igas, iboundary,maxphase, ll, iphase

 implicit none

 real, intent(in) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in) :: npart
 character(len=*), intent(in) :: dumpfile

 real,allocatable,dimension(:,:) :: dumxyzh

 integer :: i, j, iamtypei, icell, ineigh, nneigh
 real :: dx,dy,dz, rij2
 real :: hi1,hj1,hi21,hj21, q2i, q2j

 integer,save :: listneigh(maxneigh)
 real, save:: xyzcache(4,maxcellcache)
 !$omp threadprivate(xyzcache,listneigh)
 real :: fgrav(20)

 logical :: iactivei, iamdusti,iamgasi, ifilledcellcache

 character(100) :: neighbourfile

 !****************************************
 ! 1. Build kdtree and linklist
 ! --> global (shared) neighbour lists for all particles in tree cell
 !****************************************

 print*, 'Building kdtree and linklist: '

 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 print*, '- Done'

 print*, 'Allocating arrays for neighbour storage : '

 allocate(neighcount(npart))
 allocate(neighb(npart,neighmax))

 neighcount(:) = 0
 neighb(:,:) = 0

 print*, '- Done'
 print "(A,I5)", 'Maximum neighbour number allocated:  ', neighmax

 !***************************************
 ! 2. Assign neighbour lists to particles by searching shared list of host cell
 !***************************************

 print*, 'Creating neighbour lists for particles'

 ! Loop over cells

 !$omp parallel default(none) &
 !$omp shared(ncells,ll,ifirstincell,npart) &
 !$omp shared(xyzh,vxyzu,iphase) &
 !$omp shared(neighcount,neighb) &
 !$omp private(icell,i, j)&
 !$omp private(iamtypei,iamgasi,iamdusti,iactivei) &
 !$omp private(ifilledcellcache,nneigh) &
 !$omp private(hi1,hi21,hj1,hj21,rij2,q2i,q2j) &
 !$omp private(fgrav, dx,dy,dz)
 !$omp do schedule(runtime)
 over_cells: do icell=1,int(ncells)

    i = ifirstincell(icell)

    ! Skip empty/inactive cells
    if(i<=0) cycle over_cells

    ! Get neighbour list for the cell

#ifdef GRAVITY
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.,f=fgrav)
#else
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)
#endif
    ifilledcellcache = .true.

    ! Loop over particles in the cell

    over_parts: do while(i /=0)
       !print*, i, icell, ncells
       if(i<0) then ! i<0 indicates inactive particles
          i = ll(abs(i))
          cycle over_parts
       endif

       if(maxphase==maxp) then
          call get_partinfo(iphase(i), iactivei,iamdusti,iamtypei)
          iamgasi = (iamtypei ==igas)
       else
          iactivei = .true.
          iamtypei = igas
          iamdusti = .false.
          iamgasi = .true.
       endif

       ! Catches case where first particle is inactive
       if (.not.iactivei) then
          i = ll(i)
          cycle over_parts
       endif

       ! do not compute neighbours for boundary particles
       if(iamtypei ==iboundary) cycle over_parts


       ! Fill neighbour list for this particle

       neighcount(i) = 0

       over_neighbours: do ineigh = 1,nneigh
          !print*, i,ineigh, listneigh(ineigh)
          j = abs(listneigh(ineigh))

          ! Skip self-references
          if(i==j) cycle over_neighbours

          dx = xyzh(1,i) - xyzh(1,j)
          dy = xyzh(2,i) - xyzh(2,j)
          dz = xyzh(3,i) - xyzh(3,j)

          rij2 = dx*dx + dy*dy +dz*dz

          hi1 = 1.0/xyzh(4,i)
          hi21 = hi1*hi1

          q2i = rij2*hi21

          hj1 = 1.0/xyzh(4,j)
          hj21 = hj1*hj1
          q2j = rij2*hj21

          is_sph_neighbour: if(q2i < radkern2 .or. q2j < radkern2) then
             !$omp critical
             neighcount(i) = neighcount(i) + 1
             if(neighcount(i) <=neighmax) neighb(i,neighcount(i)) = j
             !$omp end critical
          endif is_sph_neighbour

       enddo over_neighbours
       ! End loop over neighbours

       i = ll(i)
    enddo over_parts
    ! End loop over particles in the cell

 enddo over_cells
 !$omp enddo
 !$omp end parallel

 ! End loop over cells in the kd-tree

 ! Do some simple stats on neighbour numbers

 meanneigh = 0.0
 sdneigh = 0.0
 neighcrit = 0.0

 call neighbours_stats(npart)

 !**************************************
 ! 3. Output neighbour lists to file
 !**************************************

 neighbourfile = 'neigh_'//TRIM(dumpfile)
 call write_neighbours(neighbourfile, npart)

 print*, 'Neighbour finding complete for file ', TRIM(dumpfile)
 deallocate(dumxyzh)
end subroutine generate_neighbour_lists



!--------------------------------------------------------------------
!+
! Calculates the mean and standard deviation of the neighbour number
! Also calculates a 5 sigma deviation from meanneigh
! (This is principally used as a diagnostic aid for structure finding
! algorithms that rely on the nearest neighbours, like CLUMPFIND)
!+
!--------------------------------------------------------------------
subroutine neighbours_stats(npart)

 implicit none
 integer, intent(in) :: npart
 integer :: ipart

 real :: minimum, maximum

 ! Calculate mean and standard deviation of neighbour counts

 maximum = maxval(neighcount)
 minimum = minval(neighcount)
 print*, 'The maximum neighbour count is ', maximum
 print*, 'The minimum neighbour count is ', minimum

 if(maximum > neighmax) then
    print*, 'WARNING! Neighbour count too large for allocated arrays'
 endif

 meanneigh = sum(neighcount)/REAL(npart)
 sdneigh = 0.0

!$omp parallel default(none) &
!$omp shared(neighcount,meanneigh,npart)&
!$omp private(ipart) &
!$omp reduction(+:sdneigh)
!$omp do schedule(runtime)
 do ipart=1,npart
    sdneigh = sdneigh+(neighcount(ipart)-meanneigh)**2
 enddo
 !$omp enddo
 !$omp end parallel

 sdneigh = sqrt(sdneigh/REAL(npart))

 print*, 'Mean neighbour number is ', meanneigh
 print*, 'Standard Deviation: ', sdneigh

 return
end subroutine neighbours_stats

!---------------------------------------
!+
! Reads in a pre-written neighbours file
!+
!---------------------------------------
subroutine read_neighbours(neighbourfile,npart)

 implicit none

 integer, intent(in) :: npart
 character(100), intent(in) ::neighbourfile

 integer :: i,j,neighcheck, tolcheck

 neigh_overload = .false.

 allocate(neighcount(npart))
 allocate(neighb(npart,neighmax))
 neighcount(:) = 0
 neighb(:,:) = 0

 print*, 'Reading neighbour file ', TRIM(neighbourfile)

 open(2, file= neighbourfile,  form = 'UNFORMATTED')

 read(2)  neighcheck, tolcheck, meanneigh,sdneigh,neighcrit

 if(neighcheck/=neighmax) print*, 'WARNING: mismatch in neighmax: ', neighmax, neighcheck

 read(2) (neighcount(i), i=1,npart)
 do i=1,npart

    if(neighcount(i) > neighmax) then
       neigh_overload = .true.
       read(2) (neighb(i,j), j=1,neighmax)
    else
       read(2) (neighb(i,j), j=1,neighcount(i))
    endif

 enddo
 close(2)

 call neighbours_stats(npart)

 if(neigh_overload) then
    print*, 'WARNING! File Read incomplete: neighbour count exceeds array size'
 else
    print*, 'File Read Complete'
 endif

end subroutine read_neighbours



!--------------------------------------------------------------------
!+
! Writes neighbour data to binary file
!+
!--------------------------------------------------------------------
subroutine write_neighbours(neighbourfile,npart)

 implicit none

 integer, intent(in) :: npart

 integer :: i,j
 character(100)::neighbourfile
 ! This is a dummy parameter, used to keep file format similar to other codes
 ! (Will probably delete this later)

 real, parameter :: tolerance = 2.0e0

 neigh_overload = .false.

 neighbourfile = TRIM(neighbourfile)

 print*, 'Writing neighbours to file ', neighbourfile

 OPEN (2, file=neighbourfile, form='unformatted')

 write(2)  neighmax, tolerance, meanneigh,sdneigh,neighcrit
 write(2) (neighcount(i), i=1,npart)
 do i=1,npart
    if(neighcount(i) > neighmax) then
       neigh_overload = .true.
       print*, 'neighbour overload: ', neighcount(i), neighmax
       write(2) (neighb(i,j), j=1,neighmax)
    else
       write(2) (neighb(i,j), j=1,neighcount(i))
    endif
 enddo

 close(2)


 if(neigh_overload) then
    print*, 'WARNING! File write incomplete: neighbour count exceeds array size'
 else
    print*, 'File Write Complete'
 endif

 return
end subroutine write_neighbours




end module
