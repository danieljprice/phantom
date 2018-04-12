!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup, producing a box/sphere with density distribution
! read from a grid file (no velocity data)
! The file is ASCII, with each cell's data written in format
!
! xface yface zface dx dy dz rho(x,y,z)
!
! (The faces indicate the 'leftward' or lower value face,
! and dx/dy/dz are the cell dimensions)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    Tinit         -- Initial temperature of system
!    gridfile      -- Filename for input densit grid
!    inputrotvel   -- Angular velocity of system (units of sound speed)
!    inputtransvel -- translational velocity of system (units of sound speed)
!    iseed         -- Random number seed
!    masstoflux    -- mass-to-magnetic flux ratio in units of critical value
!    np            -- total number of particles in Box
!    positiveBz    -- direction of B_z: T if B_z > 0
!    rcut          -- Radius of sphere to cut from grid (set <0 to avoid cutting)
!
!  DEPENDENCIES: boundary, dim, eos, externalforces, infile_utils, io,
!    kernel, options, part, physcon, prompting, ptmass, random,
!    setup_params, timestep, unifdis, units
!+
!--------------------------------------------------------------------------
module setup
 use part, only:mhd
 use dim,  only:use_dust
 implicit none
 public :: setpart

 private
 !--private module variables
 real :: xmini(3), xmaxi(3)
 real :: totmass,Tinit,rcut, inputtransvel,inputrotvel
 real :: box_velocity, masstoflux,dusttogas,pmass_dusttogas
 real(kind=8) :: udist,umass
 integer :: np,iseed,nx,ny,nz
 logical :: positiveBz
 character(len=1), parameter :: labelx(3) = (/'x','y','z'/)
 character(len=100) :: gridfile

 real,allocatable,dimension(:) :: xface,yface,zface,dx,dy,dz
 real,allocatable,dimension(:,:,:) :: rhogrid

contains

!----------------------------------------------------------------
!+
!  setup for a sphere-in-a-box
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,years,kboltz,mass_proton_cgs
 use setup_params, only:rhozero,npart_total,rmax,ihavesetupB
 use io,           only:master,fatal
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use units,        only:set_units,utime,unit_density,unit_Bfield
 use eos,          only:ieos,gmw
 use part,         only:Bxyz,Bextx,Bexty,Bextz,igas,idust,set_particle_type,kill_particle,shuffle_part
 use timestep,     only:dtmax,tmax
 use ptmass,       only:icreate_sinks,r_crit,h_acc
 use options,      only:twallmax, dtwallmax,nfulldump, iexternalforce
 use kernel,       only:hfact_default
 use random,       only:ran2
 use externalforces, only: iext_staticsine
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 real    :: vol_box,psep,cs
 real    :: przero, rpart,thetapart
 real    :: mpartmin,t_ff,area,Bzero,rmasstoflux_crit
 integer :: i,npmax, nprequired,ix,iy,iz
 logical :: iexist
 character(len=100) :: filename
 character(len=10)  :: string
 character(len=40)  :: fmt

 integer, parameter :: igridunit = 14

 npmax = size(xyzh(1,:))
 filename=trim(fileprefix)//'.setup'
 print "(/,1x,63('-'),1(/,a),/,1x,63('-'),/)",&
   '  Grid-to-box setup: Because why would you use a grid?'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename)
    call read_gridfile_header(igridunit,gridfile)
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'

    udist = 1.d16
    umass = solarm

    !
    !--prompt user for settings
    !

    ! NB - grid file contains unit data as well
    gridfile = 'grid.dat'
    call prompt('What is the grid filename? ',gridfile)

    call read_gridfile_header(igridunit,gridfile)

    !
    !--units
    !
    call set_units(dist=udist,mass=umass,G=1.d0)

    np    = int(size(xyzh(1,:))) ! Maximum particle number

    npmax = np
    call prompt('Enter the total number of particles in the box',np,0,npmax)

    Tinit = 10.0
    call prompt('What is the initial temperature of the system?',Tinit,0.0)

    rcut = -1.0
    call prompt("If you wish to cut out a sphere, what is its radius? (<0.0=no cut)",rcut, -1.0)
    inputtransvel = 0.0
    write(string,"(es10.3)") udist/utime
    call prompt("Enter the system's translational velocity in units of "//trim(adjustl(string))//" cm/s ",inputtransvel,0.0)

    inputrotvel = 0.0
    write(string,"(es10.3)") 1.0/utime
    call prompt("Enter the system's rotational velocity in units of "//trim(adjustl(string))//" s^-1 ",inputrotvel,0.0)

    iseed = -46854
    call prompt('Enter the random number seed ', iseed, -100000,0)

    !
    !--write default input file
    !
    call write_setupfile(filename) ! TODO rewrite setup file
    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
    stop
 else
    stop
 endif

 ! Read the grid file header (obtain unit data etc)
 !call read_gridfile_header(igridunit,gridfile,nx,ny,nz,udist,umass,udistchar,umasschar)
 ! Read the rest of the file

 call read_gridfile(igridunit,gridfile, mpartmin, totmass)

 !
 !--units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)

 !
 !--boundaries
 !

 xmini(1) = minval(xface)
 xmaxi(1) = maxval(xface)
 xmini(2) = minval(yface)
 xmaxi(2) = maxval(yface)
 xmini(3) = minval(zface)
 xmaxi(3) = maxval(zface)

 call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3))

 nprequired = int(totmass/mpartmin)

 print "(a, I9)", 'This requires a minimum particle number of ',nprequired

 if(nprequired > npmax) then
    print*, 'ERROR: Not enough particle memory to resolve grid'
    print*, 'Recompile and run with higher MAXP'
    stop
 endif

 if(nprequired > np) then
    print*, 'Input particle number insufficient to resolve grid'
    np = nprequired
    print*, 'New input particle number: ',np
 endif

 !
 !--general parameters
 !
 time        = 0.
 hfact       = hfact_default
 gamma       = 1.
 rmax        = rcut
 vol_box     = dxbound*dybound*dzbound
 rhozero     = totmass/vol_box
 t_ff        = sqrt(3.*pi/(32.*rhozero))
 npart_total = 0

!
!--Global system sound speed
!

 cs = sqrt(utime*gamma*kboltz*Tinit/(udist*gmw*mass_proton_cgs)) ! In code units

 polyk = cs*cs
 przero = cs*cs*rhozero

 !
 !--magnetic field -- (TODO - check the correct critical-mass-to-flux ratio for the entire box)
 !

 rmasstoflux_crit = 2./3.*0.53*sqrt(5./pi)
 if (mhd) then
    area = vol_box**0.666
    if (masstoflux > tiny(masstoflux)) then
       Bzero = totmass/(area*masstoflux*rmasstoflux_crit)
    else
       Bzero = 0.
    endif
    ihavesetupB = .true.
 else
    Bzero = 0.
 endif
 Bextx  = 0.
 Bexty  = 0.
 Bextz  = Bzero


 print "(a,i10)",' Input npart = ',np
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass,totmass*umass,' g'
 if(rcut > 0.0)  print fmt,' Cutting radius   : ',rcut,rcut*udist,' cm'
 print fmt,' Mean Density     : ',rhozero,rhozero*unit_density,' g/cm^3'
 print fmt,' Sound Speed      : ',cs,cs*udist/utime,' cm/s'
 print fmt,' Free fall time   : ',t_ff,t_ff*utime/years,' yrs'
 if (mhd) then
    print fmt,' B field (z)      : ',Bzero,Bzero*unit_Bfield*1.d6,' micro-G'
    print fmt,' Alfven speed     : ',Bzero/sqrt(rhozero),Bzero/sqrt(rhozero)*udist/utime,' cm/s'
    if (Bzero > 0.) then
       print fmt,' plasma beta      : ',przero/(0.5*Bzero*Bzero)
       print fmt,' mass-to-flux     : ',totmass/(area*Bzero)/rmasstoflux_crit
    endif
 endif
!!$ if (use_dust) then
!!$   print fmt,' dust-to-gas ratio: ',dusttogas,' '
!!$   print fmt,' dust-to-gas particle mass ratio: ',pmass_dusttogas,' '
!!$ endif
 print "(1x,50('-'))"

 !
 !--setup particles
 !

 !
 !--set particle properties
 !

 massoftype(igas)  = totmass/np

 !
 !--Loop over each grid cell, populating it with particles
 !

 do ix=1,nx
    do iy=1,ny
       do iz = 1,nz

          iseed = iseed - 1
          psep = (massoftype(igas)/rhogrid(ix,iy,iz))**(1.0/3.0)

          ! If rcut < 0, then don't cut out a sphere
          if(rcut<0.0) then
             call set_unifdis('random',id,master,xface(ix),xface(ix)+dx(ix),&
                  yface(iy),yface(iy)+dy(iy),zface(iz),zface(iz)+dz(iz),psep,&
                  hfact,npart,xyzh,nptot=npart_total,inputiseed=iseed,verbose=.false.)
          else
             call set_unifdis('random',id,master,xface(ix),xface(ix)+dx(ix),&
                  yface(iy),yface(iy)+dy(iy),zface(iz), zface(iz)+dz(iz),psep,&
                  hfact,npart,xyzh,rmax=rcut,nptot=npart_total,inputiseed=iseed,verbose=.false.)
          endif
       enddo
    enddo
 enddo

 print "(I8,a)", npart, ' particles located in grid'

 npartoftype(:)    = 0
 npartoftype(igas) = npart
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo

 ! Set up internal energy
 do i=1, npart
    if (size(vxyzu(:,i)) >= 4) vxyzu(4,i) = 1.5*polyk
 enddo

 !
 !--Now add velocities
 !

 do i=1,npart

    !
    !--Simple translational velocity in the x-direction
    !

    vxyzu(1,i) = inputtransvel

    !
    !-- Solid body rotation about the z-axis
    !

    rpart = sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i))

    thetapart = atan2(xyzh(2,i),xyzh(1,i))

    vxyzu(1,i) = vxyzu(1,i) - inputrotvel*xyzh(2,i)
    vxyzu(2,i) = inputrotvel*xyzh(1,i)

!
!-- Set up magnetic fields
!
    if (mhd) then
       Bxyz(:,i) = 0.
       if (positiveBz) then
          Bxyz(3,i) = real(Bzero,kind=kind(Bxyz))
       else
          Bxyz(3,i) = real(-Bzero,kind=kind(Bxyz))
       endif
    endif
 enddo
 !
 !--set default runtime parameters
 !
 tmax          = 10.75
 dtmax         = t_ff/100.
 ieos          = 2
 icreate_sinks = 1
 r_crit        = 5.e-2
 h_acc         = 1.e-2
 twallmax      = 604800
 dtwallmax     = 21600
 nfulldump     = 1
 iexternalforce = iext_staticsine

end subroutine setpart


!--------------------------------------------
!
! Read the header of the grid file
!
!--------------------------------------------
subroutine read_gridfile_header(igridunit, gridfile)
 implicit none

 integer, intent(in) :: igridunit
 character(len=100), intent(inout) :: gridfile
 character(len=5) :: masslabel,distlabel

 print "(a,a)", "Reading header of file ", trim(gridfile)

 open(igridunit,file=gridfile,form='formatted',status='old')
 read(igridunit,*) nx,ny,nz,umass,udist, masslabel,distlabel
 close(igridunit)

 print "(a,'(',I4,',',I4,',',I4,')')", 'Grid Cell Number (x,y,z): ', nx,ny,nz

 print "(a,a,a,1pg10.3,a)", 'Grid Distance Units: ', trim(distlabel), ' - ', udist, ' cm'
 print "(a,a,a,1pg10.3,a)", 'Grid Mass Units: ', trim(masslabel), ' - ', umass, ' g'

end subroutine read_gridfile_header

subroutine read_gridfile(igridunit,gridfile, mpartmin,totmass)

 implicit none
 integer, intent(in) :: igridunit
 character(len=100), intent(in) :: gridfile
 real, intent(out) :: mpartmin,totmass

 real :: mcell
 integer :: ix,iy,iz

 print "(a,a)", "Reading contents of file ", trim(gridfile)

 open(igridunit,file=gridfile,form='formatted',status='old')
 read(igridunit,*)

 allocate(xface(nx))
 allocate(yface(ny))
 allocate(zface(nz))

 allocate(dx(nx))
 allocate(dy(ny))
 allocate(dz(nz))

 allocate(rhogrid(nx,ny,nz))

 totmass = 0.0
 mpartmin = 1.0e30

! Read the rest of the file
 do ix=1,nx
    do iy=1,ny
       do iz=1,nz
          read(igridunit,*) xface(ix), yface(iy), zface(iz), &
              dx(ix), dy(iy), dz(iz), rhogrid(ix,iy,iz)

          mcell = rhogrid(ix,iy,iz)*dx(ix)*dy(iy)*dz(iz)
          if(mcell < mpartmin .and. mcell >0.0) mpartmin = mcell
          totmass = totmass + mcell

       enddo
    enddo
 enddo

 print "(a)", "Grid file read"
 print "(a, 1pg10.3)", 'Minimum particle mass to resolve the grid: ', mpartmin
 print "(a, 1pg10.3)", 'Total mass of the system: ', totmass

end subroutine read_gridfile

subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for grid-to-box setup routines (setup_fromgrid.f90)'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(gridfile, 'gridfile', 'Filename for input densit grid', iunit)
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','total number of particles in Box',iunit)

 call write_inopt(Tinit,'Tinit', 'Initial temperature of system', iunit)
 call write_inopt(rcut, 'rcut', 'Radius of sphere to cut from grid (set <0 to avoid cutting)',iunit)


 call write_inopt(inputtransvel, 'inputtransvel', 'translational velocity of system (units of sound speed)',iunit)

 call write_inopt(inputrotvel, 'inputrotvel', 'Angular velocity of system (units of sound speed)',iunit)

 if (mhd) then
    call write_inopt(masstoflux,'masstoflux','mass-to-magnetic flux ratio in units of critical value',iunit)
    call write_inopt(positiveBz,'positiveBz','direction of B_z: T if B_z > 0',iunit)
 endif
! if (use_dust) then
!    call write_inopt(dusttogas,'dusttogas','dust-to-gas ratio',iunit)
!    call write_inopt(pmass_dusttogas,'pmass_dusttogas','dust-to-gas particle mass ratio',iunit)
! endif
 call write_inopt(iseed,'iseed', 'Random number seed',iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 21
 integer :: ierr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(gridfile,'gridfile',db,ierr)
 call read_inopt(np,'np',db,ierr)
 call read_inopt(Tinit,'Tinit',db,ierr)
 call read_inopt(rcut,'rcut',db,ierr)
 call read_inopt(inputtransvel,'inputtransvel',db,ierr)
 call read_inopt(inputrotvel,'inputrotvel',db,ierr)

 if (mhd) then
    call read_inopt(masstoflux,'masstoflux',db,ierr)
    call read_inopt(positiveBz,'positiveBz',db,ierr)
 endif
! if (use_dust) then
!   call read_inopt(dusttogas,'dusttogas',db,ierr)
!   call read_inopt(pmass_dusttogas,'pmass_dusttogas',db,ierr)
! endif
 call read_inopt(iseed, 'iseed',db,ierr)
 call close_db(db)

end subroutine read_setupfile

end module setup

