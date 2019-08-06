!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup, producing a box filled with uniform density spheres in
! pressure equilibrium with a diffuse medium
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    box_velocity     -- translational velocity of box (units of sound speed)
!    clumpfrac        -- mass fraction of system in clumps
!    cs_clump         -- sound speed in clumps in code units
!    density_contrast -- density contrast in code units
!    iseed            -- Random number seed
!    masstoflux       -- mass-to-magnetic flux ratio in units of critical value
!    np               -- total number of particles in Box
!    positiveBz       -- direction of B_z: T if B_z > 0
!    r_clump          -- radius of clumps in code units
!    totmass          -- total mass of system in code units
!    udist            -- distance unit in cm
!    umass            -- mass unit in g
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
 real :: density_contrast,totmass,r_clump,cs_clump,clumpfrac
 real :: box_velocity, masstoflux,dusttogas,pmass_dusttogas
 real(kind=8) :: udist,umass
 integer :: np,iseed
 logical :: positiveBz
 character(len=1), parameter :: labelx(3) = (/'x','y','z'/)

contains

!----------------------------------------------------------------
!+
!  setup for a sphere-in-a-box
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,years
 use setup_params, only:rhozero,npart_total,rmax,ihavesetupB
 use io,           only:master,fatal
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use units,        only:set_units,utime,unit_density,unit_Bfield
 use eos,          only:polyk2,ieos
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

 real    :: vol_box,psep,psep_box, mclump
 real    :: vol_clump,dens_clump,dens_medium,cs_medium,przero
 real    :: totmass_clumps,t_ff,area,Bzero,rmasstoflux_crit,sep,icell
 integer, dimension(3) :: nx
 real, dimension(3) :: dr
 real,allocatable,dimension(:,:) :: xsphere
 integer :: i,iclump,ipart,jpart, nclumps, npart_clump, npart_medium
 integer :: npartclump,npmax,ninside
 logical :: iexist
 character(len=100) :: filename
 character(len=10)  :: string
 character(len=40)  :: fmt

 npmax = size(xyzh(1,:))
 filename=trim(fileprefix)//'.setup'
 print "(/,1x,63('-'),1(/,a),/,1x,63('-'),/)",&
   '  Multiple-spheres-in-box setup: Beyond Archimedes.'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename)
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    udist = 1.d16
    umass = solarm
    call prompt('Enter code units of distance in cm ',udist,0.)
    call prompt('Enter code units of mass in g ',umass,0.)
    !
    !--units
    !
    call set_units(dist=udist,mass=umass,G=1.d0)
    !
    !--prompt user for settings
    !

    np    = int(2.0/3.0*size(xyzh(1,:))) ! approx max number allowed in sphere given size(xyzh(1,:))
    npmax = np
    call prompt('Enter the total number of particles in the box',np,0,npmax)
    xmini = -2.
    xmaxi = 2.
    write(string,"(es10.3)") udist
    do i=1,3
       call prompt('Enter '//labelx(i)//'min of box in units of '//trim(adjustl(string))//' cm',xmini(i))
       call prompt('Enter '//labelx(i)//'max of box in units of '//trim(adjustl(string))//' cm',xmaxi(i),xmini(i))
    enddo

    r_clump = 0.1
    call prompt('Enter radius of clumps in units of '//trim(adjustl(string))//' cm',r_clump,0.)

    density_contrast = 30.0
    call prompt('Enter density contrast between clumps and box ',density_contrast,1.)

    totmass = 1.0
    write(string,"(es10.3)") umass
    call prompt('Enter total mass in box in units of '//trim(adjustl(string))//' g',totmass,0.)

    clumpfrac = 0.5
    call prompt('And what fraction of this mass is in clumps? ',clumpfrac,0.,1.)

    cs_clump = 0.19
    write(string,"(es10.3)") udist/utime
    call prompt('Enter sound speed in sphere in units of '//trim(adjustl(string))//' cm/s',cs_clump,0.)

    box_velocity = 50.0
    call prompt("Enter the box's translational velocity in units of the sound speed: ",box_velocity)

!!$    if (use_dust) then
!!$       dusttogas = 0.01
!!$       call prompt('Enter dust-to-gas ratio ',dusttogas,0.)
!!$       pmass_dusttogas = dusttogas*10.0
!!$       call prompt('Enter dust-to-gas particle mass ratio ',pmass_dusttogas,0.)
!!$    endif

    iseed = -46854
    call prompt('Enter the random number seed ', iseed, -100000,0)

    !
    !--write default input file
    !
    call write_setupfile(filename) ! TODO - rewrite setup file
    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
    stop
 else
    stop
 endif
 !
 !--units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !
 !--boundaries
 !
 call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3))
 !
 !--general parameters
 !
 time        = 0.
 hfact       = hfact_default
 gamma       = 1.
 rmax        = r_clump
 vol_box     = dxbound*dybound*dzbound
 vol_clump  = 4./3.*pi*r_clump**3
 rhozero     = totmass/vol_box
 t_ff        = sqrt(3.*pi/(32.*rhozero))
 npart_total = 0

! Calculate the density of the clumps, ensuring they are sub-Jeans
! Then obtain properties of the medium to ensure pressure equilibrium

 dens_clump = pi*cs_clump*cs_clump/(r_clump*r_clump)
 dens_medium = dens_clump/density_contrast
 cs_medium   = cs_clump*sqrt(density_contrast)

 ! Total mass of the system

 totmass = dens_medium*vol_box/(1.0-clumpfrac*(1.0-1.0/density_contrast))
 !totmass = 1.0
 ! Total mass of the system in clumps
 totmass_clumps = clumpfrac*totmass
 !totmass_box = (vol_box - vol_sphere)*dens_medium
 !totmass     = totmass_box + totmass_sphere

 przero = cs_clump**2*dens_clump

 !
 !--temperature set to give a pressure equilibrium
 !
 polyk  = cs_clump**2
 polyk2 = cs_medium**2


 print "(a,i10)",' Input npart = ',np
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass,totmass*umass,' g'
 print fmt,' Mass Fraction in clumps   : ', clumpfrac
 print fmt,' Radius of clumps : ',r_clump,r_clump*udist,' cm'
 print fmt,' Density sphere   : ',dens_clump,dens_clump*unit_density,' g/cm^3'
 print fmt,' Density medium   : ',dens_medium,dens_medium*unit_density,' g/cm^3'
 print fmt,' cs in sphere     : ',cs_clump,cs_clump*udist/utime,' cm/s'
 print fmt,' cs in medium     : ',cs_medium,cs_medium*udist/utime,' cm/s'
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

 print "(a, 1pg10.3)", 'The total mass of the system is ', totmass
 ! First, setup particles in clumps

 ! Calculate total number of clumps
 nclumps = int(clumpfrac*totmass/(vol_clump*dens_clump))
 mclump = clumpfrac*totmass/real(nclumps)

 allocate(xsphere(3,nclumps))

 print "(a, I3)", 'Total number of clumps to create: ', nclumps
 print "(a, 1pg10.3)", 'Clump Mass: ', mclump


 ! Given total number of particles and system mass, calculate number of particles per clump
 npart_clump = int(np*clumpfrac/nclumps)

 print "(a, I5, a)", 'Each clump will possess ', npart_clump, ' particles'
 psep = (vol_clump/npart_clump)**(1.0/3.0)


 ! Randomly select locations for all clumps
 ! To avoid clump overlap, we discretise the box with cell size = 1.5*r_clump
 ! We randomly select cells, we can avoid overlaps
 ! This is acceptable as long as the density contrast is >>1


 ! First calculate number of grid cells in each dimension
 do i=1,3
    nx(i) = int((xmaxi(i)-xmini(i))/(1.5*r_clump))
 enddo

 print*, nx(1)*nx(2)*nx(3)

 do iclump=1,nclumps
    do i=1,3
       ! Randomly select a cell index for dimension i
       icell = 1+ ran2(iseed)*(nx(i)-1) ! Add 0.5 as we want the centre of the cell

       xsphere(i,iclump) = xmini(i) + icell*1.5*r_clump

       ! Shift the sphere slightly from the centre

       xsphere(i,iclump) = xsphere(i,iclump) + (-1.0 + 2.0*ran2(iseed))*r_clump/10.0
       !if (abs(xsphere(i,iclump)-xmaxi(i))<r_clump) xsphere(i,iclump) = xmaxi(i)-1.1*r_clump
       !if (abs(xsphere(i,iclump)-xmini(i))<r_clump) xsphere(i,iclump) = xmini(i) + 1.1*r_clump
!       xsphere(i,iclump) = xmini(i) + ran2(iseed)*(xmaxi(i)-xmini(i)) Alternative case without discretising
    enddo
 enddo


 ! Loop over number of clumps

 jpart = 1
 do iclump = 1,nclumps

    ! Create a sphere from randomly distributed particles

    call set_unifdis('random',id,master,-r_clump,r_clump,-r_clump,r_clump,-r_clump,r_clump,psep,&
         hfact,npart,xyzh,rmax=r_clump, nptot = npart_total)

    print "(a,I5,a, I3)", 'Placed ',npart-jpart+1, ' particles in clump ', iclump

    ! Move the sphere to a random location inside the box
    ! (ensuring entire sphere is contained within box limits)

    print "(a)", "Moving clump to location:"
    print*, xsphere(:,iclump)

    ! Move particles, and give them internal energies along the way
    do ipart=jpart, npart
       do i=1,3
          xyzh(i,ipart) = xyzh(i,ipart)+xsphere(i,iclump)
       enddo
       if (size(vxyzu(:,ipart)) >= 4) vxyzu(4,ipart) = 1.5*polyk
    enddo

    jpart = jpart + (npart-jpart+1)

 enddo

 print "(a)", "All clumps placed in the box"

 npartclump = npart

 !--setup surrounding low density medium
 psep_box = (vol_box/(np-npartclump))**(1.0/3.0)  ! calculate psep in box

 !
 !-- Check that low density medium particle count does not exceed maxp (TODO)
 !

 npart_medium = int(vol_box/psep_box**(3.0))

 if (npart+ npart_medium > npmax) then
    print "(a)", 'Warning! Low density medium particles will exceed maxp'
    print*, npart, npart_medium, npmax, psep_box,vol_box
 endif

 print "(a,es10.3)",' Particle separation in clumps = ',psep
 print "(a,I6)",' Total Particles in Clumps = ', npart_clump
 print "(a,es10.3)",' Particle separation in box = ', psep_box
 print *,' Total Particles in Medium = ', npart_medium

 call set_unifdis('random',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep_box, &
                   hfact,npart,xyzh,nptot=npart_total)
 print "(a,i10,a)",' added ',npart-npartclump,' particles in low-density medium'
 print*, ""

 ! Set up internal energy for the medium
 do ipart=jpart, npart
    if (size(vxyzu(:,ipart)) >= 4) vxyzu(4,ipart) = 1.5*polyk2
 enddo

 !
 !--set particle properties
 !

 npartoftype(:)    = 0
 npartoftype(igas) = npart
 massoftype(igas)  = totmass/npart_total
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo

 print*, massoftype(igas)

 ! Check for medium particles that stray within the clumps
 ! If found, reset their internal energy to avoid spurious gradients

 ninside = 0
 do ipart = jpart, npart
    do iclump=1,nclumps
       do i=1,3
          dr(i) = xyzh(i,ipart) - xsphere(i,iclump)
       enddo

       sep = sqrt(dr(1)*dr(1) + dr(2)*dr(2)+dr(3)*dr(3))

       if (sep<r_clump) then
          ninside = ninside+1
          vxyzu(4,ipart) = 1.5*polyk

          exit
       endif
    enddo
 enddo

 np = npart

 print "(a, I6)", 'Medium particles found inside clumps: ', ninside
 print "(a)", 'These particles have had their sound speeds updated to match that of the clumps'


 npart_total = np

 !
 !--Set dust (TODO)
!!$ if (use_dust) then
!!$   !--particle separation in dust sphere & sdjust for close-packed lattice
!!$   psep = (vol_sphere/pmass_dusttogas)**(1./3.)/real(nx)
!!$   psep = psep*sqrt(2.)**(1./3.)
!!$   call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
!!$                     hfact,npart,xyzh,rmax=r_sphere,nptot=npart_total)
!!$   npartoftype(idust) = npart_total - npartoftype(igas)
!!$   massoftype(idust)  = totmass_sphere*dusttogas/npartoftype(idust)
!!$   !
!!$   do i = npartoftype(igas)+1,npart
!!$     call set_particle_type(i,idust)
!!$   enddo
!!$   !
!!$   print "(a,4(i10,1x))", ' particle numbers: (gas_total, gas_sphere, dust, total): ' &
!!$                        , npartoftype(igas),npartsphere,npartoftype(idust),npart
!!$   print "(a,2es10.3)"  , ' particle masses: (gas,dust): ',massoftype(igas),massoftype(idust)
!!$ else
!   print "(a,3(i10,1x))", ' particle numbers: (sphere, low-density medium, total): ' &
!                        , npartsphere, npart-npartsphere,npart
 print "(a,es10.3)",' particle mass = ',massoftype(igas)
!!$ endif



 !
 !--Simple translational velocity in the x-direction
 !

 do i=1,npart
    vxyzu(1,i) = box_velocity*cs_clump
 enddo

! Set B fields using standard set_Bfield routine outside of setpart
 if (mhd) ihavesetupB = .false.
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

subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for spheres-in-box setup routines (setup_clumpybox.f90)'
 write(iunit,"(/,a)") '# units'
 call write_inopt(udist,'udist','distance unit in cm',iunit)
 call write_inopt(umass,'umass','mass unit in g',iunit)
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','total number of particles in Box',iunit)
 write(iunit,"(/,a)") '# options for box'
 do i=1,3
    call write_inopt(xmini(i),labelx(i)//'min',labelx(i)//' min',iunit)
    call write_inopt(xmaxi(i),labelx(i)//'max',labelx(i)//' max',iunit)
 enddo
 write(iunit,"(/,a)") '# options for sphere'
 call write_inopt(r_clump,'r_clump','radius of clumps in code units',iunit)
 call write_inopt(density_contrast,'density_contrast','density contrast in code units',iunit)
 call write_inopt(totmass,'totmass','total mass of system in code units',iunit)
 call write_inopt(clumpfrac,'clumpfrac','mass fraction of system in clumps', iunit)
 call write_inopt(cs_clump,'cs_clump','sound speed in clumps in code units',iunit)
 call write_inopt(box_velocity, 'box_velocity', 'translational velocity of box (units of sound speed)',iunit)

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
 integer :: ierr,i
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(udist,'udist',db,ierr)
 call read_inopt(umass,'umass',db,ierr)
 call read_inopt(np,'np',db,ierr)
 do i=1,3
    call read_inopt(xmini(i),labelx(i)//'min',db,ierr)
    call read_inopt(xmaxi(i),labelx(i)//'max',db,ierr)
 enddo
 call read_inopt(r_clump,'r_clump',db,ierr)
 call read_inopt(density_contrast,'density_contrast',db,ierr)
 call read_inopt(totmass,'totmass',db,ierr)
 call read_inopt(clumpfrac,'clumpfrac',db,ierr)
 call read_inopt(cs_clump,'cs_clump',db,ierr)
 call read_inopt(box_velocity, 'box_velocity', db,ierr)

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

