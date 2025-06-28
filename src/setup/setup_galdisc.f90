!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! this module does setup for galactic discs with live or analytic stars
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - angvel            : *velocity for rotation curve [unit distance (cm/s)]*
!   - dust_to_gas_ratio : *dust to gas ratio*
!   - np                : *number of particles*
!   - partdist          : *particle distribution (r=random,o=observed)*
!   - rcyl              : *outer radius for particle setup [kpc]*
!   - rcylin            : *inner radius for particle setup [kpc]*
!   - thermal           : *initial temperature in Kelvin*
!   - totmass           : *total mass*
!   - use_bulge         : *include live bulge*
!   - use_gas           : *include gas*
!   - use_halo          : *include live halo*
!   - use_live_stars    : *use live star particles*
!   - use_star          : *include live stars*
!   - vset              : *velocity profile (r=rotation curve,f=flat,c=constant)*
!
! :Dependencies: datafiles, dim, extern_spiral, externalforces,
!   infile_utils, io, kernel, options, part, physcon, prompting, random,
!   set_dust, setup_params, units
!
 use options, only:use_dustfrac
 use dim,     only:maxp,maxvxyzu
 implicit none
 public :: setpart

 private

 ! Module variables for setup parameters
 logical :: use_live_stars,use_gas,use_star,use_bulge,use_halo
 integer :: np
 real    :: thermal,totmass,rcylin,rcyl,dust_to_gas_ratio
 character(len=1) :: partdist,vset
 real(kind=8) :: angvel

contains
!--------------------------------------------------------------------------
!
! This subroutine is a utility for setting up isolated galactic discs
!
! The code has two main options. The code can setup a gas disc alone,
! which should then be exposed to some galactic potentials.
! Alternatively it will assume a live stellar component (with optional
! live bulge and dark matter) where stars are simply N-body particles.
!
! If setting a live disc then the code reads asciifiles of the
! galsetup python routines included with phantom.
! Seperate files are needded for each component (gas, stars, bulge, halo).
! This routine can read in ANY IC files providing they are in the format
!     x y z m vx vy vz
! so long as the different components are stored in the different asciifiles.
! Important note: the setup needs to know the number of each type of
! galactic component which requires the small galsetic.txt file.
!
! Created by Alex Pettitt from a similar routine by Clare Dobbs modified by
! Daniel Price for use in phantom.
!
!--------------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,            only:use_dust,h2chemistry
 use setup_params,   only:rhozero
 use physcon,        only:Rg,pi,solarm,pc,kpc
 use units,          only:udist,unit_ergg,set_units
 use part,           only:abundance,iHI,dustfrac,ndusttypes
 use options,        only:iexternalforce,icooling,nfulldump
 use io,             only:fatal,master
 use infile_utils,   only:get_options
 use set_dust,       only:set_dustfrac
 use kernel,         only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: i,ierr,iseed
 real :: h2ratio,gmw

 !
 ! set code units
 !
 call set_units(dist=100.*pc,mass=1.d5*solarm,G=1.d0)
 !
 ! set input file options
 !
 iexternalforce = 8   !8=galactic disk potentials
 icooling       = 0   !1=cooling on, 0=off
 nfulldump      = 1
 hfact = hfact_default

 ! default parameters
 use_live_stars = .false.
 use_gas = .false.
 use_star = .false.
 use_bulge = .false.
 use_halo = .false.
 np = min(maxp,500000)
 thermal = 10000.
 totmass = 1.e4
 rcylin = 0.
 rcyl = 10.
 partdist = 'r'
 vset = 'r'
 angvel = 2.20e7
 dust_to_gas_ratio = 0.01

 ! Get setup parameters from file or interactive setup
 if (id==master) print "(/,a,/)",'  >>> Setting up galactic disc <<<'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 if (maxvxyzu >= 4) then
    gamma = 1.6666667 ! non-isothermal PV^gamma=K
 else
    gamma = 1.0 ! isothermal PV = K as T=const
 endif
 time = 0.

 ! thermal energy
 h2ratio = 0.
 gmw = (2.*h2ratio+(1.-2.*h2ratio)+0.4)/ &
       (0.1+h2ratio+(1.-2.*h2ratio))
 print*,' assuming mean molecular weight is ',gmw
 thermal = 3.0/2.0*thermal*Rg/gmw/unit_ergg
 polyk = 2./3.*thermal
 print*,'polyk (from thermal) = ',polyk

 ! convert radii to code units
 rcylin = rcylin*kpc/udist
 rcyl   = rcyl*kpc/udist

 if (use_live_stars) then
    call setup_live_stars(id,xyzh,vxyzu,massoftype,npartoftype,npart,totmass,thermal,hfact)
 else
    call setup_gas_only_disc(xyzh,massoftype,npartoftype,npart,totmass,rhozero,hfact,iseed)
    call setup_velocities(npart,xyzh,vxyzu,iexternalforce,iseed) ! setup velocities for gas-only disc
    ! use one fluid dust, and allow for adding dust to the gas
    if (use_dust) use_dustfrac = .true.
 endif

 ! set abundances
 if (h2chemistry) then
    do i=1,npart
       abundance(:,i)   = 0.
       abundance(iHI,i) = 1.  ! assume all atomic hydrogen initially
    enddo
 endif

 ! set dust fractions
 if (use_dustfrac) then
    ndusttypes = 1
    do i=1,npart
       call set_dustfrac(dust_to_gas_ratio,dustfrac(:,i))
    enddo
 endif

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  routine to read files containing tabulated density profile
!  requires that PHANTOM_DIR is set
!+
!-----------------------------------------------------------------------
subroutine readinr(file_in,filelen,cdf_ar,r_ar,ierr)
 use datafiles, only:find_phantom_datafile
 character(len=*), intent(in)  :: file_in
 integer,          intent(in)  :: filelen
 real,             intent(out) :: cdf_ar(filelen),r_ar(filelen)
 integer,          intent(out) :: ierr
 integer, parameter :: lu=11
 integer :: i
 real    :: cdf,r
 logical :: iexist
 character(len=120) :: filename

 filename = find_phantom_datafile(file_in,'gal_radii_cdfs')

 inquire(file=trim(filename),exist=iexist)
 if (.not.iexist) then
    print*,'ERROR: '//trim(filename)//' not found'
    ierr = 1
    return
 endif
 i=1
 open(unit=lu,file=trim(filename),status='old',iostat=ierr)
 if (ierr /= 0) then
    write(*,*) 'ERROR opening '//trim(filename)
    return
 endif

 ierr = 0
 do while (i<=filelen .and. ierr==0)
    read(lu,*,iostat=ierr) cdf,r
    cdf_ar(i) = cdf
    r_ar(i)   = r
    i=i+1
 enddo

 if (ierr > 0) then
    write(*,"(1x,a,i2,a)") 'ERROR reading '//trim(filename)//': aborting'
 elseif (ierr < 0) then
    write(*,*) 'ERROR reading '//trim(filename)//' (hit end of file/record): aborting'
 endif
 close(lu)

end subroutine readinr

!-----------------------------------------------------------------------
!+
!  Draw randomly from MW disc density profile tabulated in files
!+
!-----------------------------------------------------------------------
subroutine cdfplacement(noinbin,lenbin,cdf,r_ar,xyzh,itot,Npart,hi)
 use random, only:ran2
 real,    intent(inout) :: xyzh(:,:)
 integer, intent(inout) :: noinbin,itot
 real,    intent(in)    :: hi
 integer, intent(in)    :: lenbin,Npart
 real,    intent(in)    :: cdf(lenbin),r_ar(lenbin)
 real    :: thisran,phii,xi,yi,zi,zmax
 integer :: i,j,iseed

 zmax=10.*0.04        !-simple vertical dependance
!--Loop over particles assigned to this radii bin
 do i=1,noinbin,1
    thisran = ran2(iseed)
    !--Loop over the cdf & r arrays
    j=1
    do j=1,lenbin,1
       !--if we find the nearest random number partner, assign particle position
       if (thisran<=cdf(j)) then
          itot = itot+1
          phii = ran2(iseed)*2.0*3.14159    !Pluck a random azimuth.
          xi   = r_ar(j)*cos(phii)
          yi   = r_ar(j)*sin(phii)
          zi   = 2.0*(zmax*ran2(iseed) - 0.5*zmax)   !Dont care about z.
          xyzh(1,itot) = 10.*xi!+ran2(iseed)
          xyzh(2,itot) = 10.*yi!+ran2(iseed)
          xyzh(3,itot) = 10.*zi!+ran2(iseed)
          xyzh(4,itot) = hi
          exit
          !endif
          !The case where we reach end of the array (ie ran2~1.0)
       elseif (j==lenbin) then
          itot = itot+1
          phii = ran2(iseed)*2.0*3.14159    !Pluck a random azimuth.
          xi   = r_ar(lenbin)*cos(phii)
          yi   = r_ar(lenbin)*sin(phii)
          zi   = 2.0*(zmax*ran2(iseed) - 0.5*zmax)   !Dont care about z.
          xyzh(1,itot) = 10.*xi
          xyzh(2,itot) = 10.*yi
          xyzh(3,itot) = 10.*zi
          xyzh(4,itot) = hi
          exit
       endif
    enddo
 enddo

end subroutine cdfplacement

!-----------------------------------------------------------------------
!+
!  Calculate Gaussian random velocity dispersion
!+
!-----------------------------------------------------------------------
real function vdisp(iseed,disp)
 use random, only:ran2
 integer, intent(inout) :: iseed
 real,    intent(in)    :: disp
 real :: rand,prob

 rand = huge(rand)
 prob = 0.
 vdisp = 0.
 do while (prob < rand)
    vdisp = 40.*(ran2(iseed)-0.5)
    prob = exp(-(vdisp/disp)**2./2.)
    rand = ran2(iseed)
 enddo

end function vdisp

!-----------------------------------------------------------------------
!+
!  Add a component by reading particle data from file
!+
!-----------------------------------------------------------------------
subroutine add_component(file_in,itype,istart,iend,xyzh,vxyzu,massoftype,npartoftype,&
                         totmass,thermal,hfact,totvol,component_name)
 use part,      only:set_particle_type
 use units,     only:udist,umass,unit_velocity
 use datafiles, only:find_phantom_datafile
 character(len=*), intent(in)    :: file_in,component_name
 integer,          intent(in)    :: itype,istart,iend
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 integer,          intent(inout) :: npartoftype(:)
 real,             intent(inout) :: totmass
 real,             intent(in)    :: thermal,hfact,totvol
 character(len=120) :: filename
 integer :: i,lu
 real :: xis,yis,zis,mis,vxis,vyis,vzis,phaseis
 real :: rhozero,h

 filename = find_phantom_datafile(file_in,'isolatedgalaxy/arpic_lowrestest')

 open(newunit=lu,file=filename,form='formatted')
 i = istart
 do while(i <= iend)
    read(lu,*) xis,yis,zis,mis,vxis,vyis,vzis,phaseis
    !--Set from above converting from SI to cgs then code units
    xyzh(1,i) = xis*100./udist
    xyzh(2,i) = yis*100./udist
    xyzh(3,i) = zis*100./udist
    vxyzu(1,i) = vxis*100./unit_velocity
    vxyzu(2,i) = vyis*100./unit_velocity
    vxyzu(3,i) = vzis*100./unit_velocity
    mis = 1000.*mis/umass
    totmass = totmass + mis
    call set_particle_type(i,itype)
    if (maxvxyzu >= 4) then
       vxyzu(4,i) = thermal
    endif
    i = i + 1
 enddo
 close(lu)

 ! Calculate mass per particle, density and smoothing length
 massoftype(itype) = totmass/real(npartoftype(itype))
 rhozero = totmass/totvol
 h = hfact*(massoftype(itype)/rhozero)**(1./3.)

 ! Set smoothing lengths for all particles of this type
 do i=istart,iend
    xyzh(4,i) = h
 enddo

 print*,'Set '//trim(component_name)//' with mean density',rhozero,'smoothing length',h

end subroutine add_component

!-----------------------------------------------------------------------
!+
!  Setup velocities for galactic disc
!+
!-----------------------------------------------------------------------
subroutine setup_velocities(npart,xyzh,vxyzu,iexternalforce,iseed)
 use units,          only:unit_velocity
 use externalforces, only:externalforce,initialise_externalforces
 use extern_spiral,  only:LogDisc
 use dim,            only:maxvxyzu
 integer, intent(in)    :: npart,iexternalforce
 integer, intent(inout) :: iseed
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real :: radius2,disp
 real :: fextxi,fextyi,fextzi,poti,timei,dr,r
 real(kind=8) :: vcirc,phi
 integer :: i,ierr

 if (vset=='c') then
    !--Set parameters for rotation curve without radial dependance
    angvel = angvel / unit_velocity
 endif
 ierr=0
 timei=0.
 call initialise_externalforces(iexternalforce,ierr)
 do i=1,npart
    fextxi=0.
    fextyi=0.
    fextzi=0.
    poti  =0.
    radius2 = xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i)
    if (vset=='r') then
       !--Pull an initial velocity from the actual rotation curve defined by .in file
       call externalforce(iexternalforce,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i), &
                      timei,fextxi,fextyi,fextzi,poti)
       angvel=-sqrt( sqrt(fextxi**2+fextyi**2) *sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) ) )
    elseif (vset=='f') then
       !--Pull an initial velocity from a simple flat profile
       dr   = 1./sqrt(radius2+xyzh(3,i)*xyzh(3,i))     !1/r
       r    = 1./dr
       call LogDisc(xyzh(1,i),xyzh(2,i),xyzh(3,i),radius2,phi,fextxi,fextyi,fextzi)
       angvel=-sqrt( sqrt(fextxi**2+fextyi**2) *sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) ) )
    endif

    vcirc   = sqrt(radius2/(1.0+radius2))
    phi     = atan2(xyzh(2,i),xyzh(1,i))

    vxyzu(1,i) = -angvel*vcirc*sin(phi)
    vxyzu(2,i) = +angvel*vcirc*cos(phi)
    vxyzu(3,i) = 0.

    !--Add Gaussian random velocity dispersion
    !--For 5 km/s dispersion, use Gaussian with sigma=disp=5
    !--Set velocity dispersion parameter:
    disp=5.
    vxyzu(1,i) = vxyzu(1,i) + vdisp(iseed,disp)*100000./unit_velocity
    vxyzu(2,i) = vxyzu(2,i) + vdisp(iseed,disp)*100000./unit_velocity
    vxyzu(3,i) = vxyzu(3,i) + vdisp(iseed,disp)*100000./unit_velocity
    if (maxvxyzu >= 4) vxyzu(4,i) = thermal
 enddo

end subroutine setup_velocities

!-----------------------------------------------------------------------
!+
!  Setup gas-only disc
!+
!-----------------------------------------------------------------------
subroutine setup_gas_only_disc(xyzh,massoftype,npartoftype,npart,totmass,rhozero,hfact,iseed)
 use part,           only:igas
 use io,             only:fatal
 use physcon,        only:pi,kpc,solarm
 use units,          only:udist,umass
 use random,         only:ran2
 integer,          intent(out)   :: iseed
 real,             intent(inout) :: xyzh(:,:),massoftype(:)
 integer,          intent(inout) :: npartoftype(:)
 integer,          intent(out)   :: npart
 real,             intent(inout) :: totmass
 real,             intent(in)    :: hfact
 real,             intent(out)   :: rhozero
 real :: faclod,xmin,xmax,ymin,ymax,zmin,zmax
 real :: rmax,rcyl2,rcylin2,totvol
 real :: xmax5,ymax5,zmax5,xi,yi,zi,r2,h1
 integer :: i,itot,it
 integer, parameter, dimension(5) :: lenAr=(/1000,3000,2670,1830,4500/)
 real,    dimension(5) :: mratios=(/0.0008,0.0788,0.2293,0.1750,0.5161/)
 integer :: mratiosi(5)
 integer :: ierrf(5)
 real :: cdf1(lenAr(1)),rp1(lenAr(1))
 real :: cdf2(lenAr(2)),rp2(lenAr(2))
 real :: cdf3(lenAr(3)),rp3(lenAr(3))
 real :: cdf4(lenAr(4)),rp4(lenAr(4))
 real :: cdf5(lenAr(5)),rp5(lenAr(5))

 ! no live stars, gas-only disc
 rcylin = rcylin*kpc/udist
 rcyl   = rcyl*kpc/udist
 !--Height dist. comp to width of disk.
 faclod = 0.04
 zmax = rcyl*faclod
 zmin = -zmax
 xmax = rcyl
 ymax = rcyl
 xmin = -xmax
 ymin = -ymax
 rmax = sqrt(rcyl*rcyl + zmax*zmax)
 rcyl2 = rcyl*rcyl
 rcylin2 = rcylin*rcylin
 !--initialise random number generator
 iseed = -6485
 print "(3(a,f10.3),a)",' galactic disc setup... rmin = ',rcylin,' rmax = ',rcyl,' in units of ',udist/kpc,' kpc'
 xi = ran2(iseed)
 if (np> maxp) call fatal('setup','npart > maxp; use ./phantomsetup --maxp=10000000')
 npartoftype(1) = np
 npart = np

 print "(1x,a,es10.3,a,1pg10.3,a)",'Mass is in units of ',umass,' g (',umass/solarm,' solar masses)'
 massoftype(igas) = totmass/real(np)

 totvol = pi*(zmax-zmin)*(rcyl2 - rcylin2)
 rhozero = totmass/totvol
 print "(a,es10.3)",' mean density = ',rhozero
 h1 = hfact*(massoftype(1)/rhozero)**(1./3.)

 if (trim(partdist)=='o') then
    !--Loop for pseudo-random placement (from observed distribution)
    print "(a)",' Realistic gas distribution requires location of CDF(r) files:'
    itot=0
    call readinr('r0to1kpc.txt',lenAr(1),cdf1,rp1,ierrf(1))
    call readinr('r1to4kpc.txt',lenAr(2),cdf2,rp2,ierrf(2))
    call readinr('r4to6kpc.txt',lenAr(3),cdf3,rp3,ierrf(3))
    call readinr('r6to9kpc.txt',lenAr(4),cdf4,rp4,ierrf(4))
    call readinr('r9to13kpc.txt',lenAr(5),cdf5,rp5,ierrf(5))
    if (any(ierrf(1:5) /= 0)) call fatal('setup_galdisc','error reading cdf files')
    !And now place:
    do it=1,5,1
       mratiosi(it)=nint(mratios(it)*np)
    enddo
    print*,'Particles in each distribution bin:', mratiosi
    print*,'Total set:', mratiosi(1)+mratiosi(2)+mratiosi(3)+mratiosi(4)+mratiosi(5)
    call cdfplacement(mratiosi(1),lenAr(1),cdf1,rp1,xyzh,itot,npart,h1)
    call cdfplacement(mratiosi(2),lenAr(2),cdf2,rp2,xyzh,itot,npart,h1)
    call cdfplacement(mratiosi(3),lenAr(3),cdf3,rp3,xyzh,itot,npart,h1)
    call cdfplacement(mratiosi(4),lenAr(4),cdf4,rp4,xyzh,itot,npart,h1)
    call cdfplacement(mratiosi(5),lenAr(5),cdf5,rp5,xyzh,itot,npart,h1)
 else
    !--Loop for purely randomly placed particles
    xmax5 = 0.5*xmax
    ymax5 = 0.5*ymax
    zmax5 = 0.5*zmax
    i = 0
    over_npart: do while(i < npart)
       !--Initialise the particle placements randomly
       xi = 2.0*(xmax*ran2(iseed) - xmax5)
       yi = 2.0*(ymax*ran2(iseed) - ymax5)
       zi = 2.0*(zmax*ran2(iseed) - zmax5)
       r2 = xi*xi + yi*yi
       if (r2 > rcyl2 .or. r2 < rcylin2) cycle over_npart
       i = i + 1
       if (mod(i,1000000)==0) print*,i
       if (i > maxp) stop 'error, i > maxp; use ./phantomsetup --maxp=N'
       xyzh(1,i) = xi
       xyzh(2,i) = yi
       xyzh(3,i) = zi
       xyzh(4,i) = h1
    enddo over_npart
 endif
 print*,npart,' particles set, each of mass',massoftype(1)

end subroutine setup_gas_only_disc

!-----------------------------------------------------------------------
!+
!  Setup live stars from IC files
!+
!-----------------------------------------------------------------------
subroutine setup_live_stars(id,xyzh,vxyzu,massoftype,npartoftype,npart,totmass,thermal,hfact)
 use part,      only:igas,istar,ibulge,idarkmatter
 use io,        only:fatal,master
 use datafiles, only:find_phantom_datafile
 use physcon,   only:pi,kpc
 use units,     only:udist
 integer,          intent(in)    :: id
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 integer,          intent(inout) :: npartoftype(:)
 integer,          intent(out)   :: npart
 real,             intent(out)   :: totmass
 real,             intent(in)    :: thermal,hfact
 character(120) :: galsetupic
 real :: totmassD,totmassG,totmassB,totmassH
 real :: totvol,totvolB,totvolH
 real :: zmax,zmin

 ! live star setup
 npartoftype = 0
 massoftype = 0.0

 ! values for initial guesstimate of smoothing lentgths
 rcyl  = 13.*kpc/udist  ! Override for live stars
 zmax  = 0.5*kpc/udist
 zmin  = -zmax
 totvol  = pi*(zmax-zmin)*(rcyl*rcyl) ! disc
 totvolB = 4./3.*pi*(0.2*rcyl)**3     ! inner spheroid
 totvolH = 4./3.*pi*(5.0*rcyl)**3     ! halo

 galsetupic = find_phantom_datafile('galsetic.txt','isolatedgalaxy/arpic_lowrestest')
 print "(a)",' opening '//trim(galsetupic)
 call read_ic_parameters(galsetupic,npartoftype)
 print*,'Read in the IC parameter file: ',galsetupic
 print*,' ngas  =',npartoftype(igas)
 print*,' nstar =',npartoftype(istar)
 print*,' nbulge=',npartoftype(ibulge)
 print*,' nhalo =',npartoftype(idarkmatter)

 !--Checks:
 if (use_gas .and. (npartoftype(igas) <= 0)) then
    call fatal('setup_isogal','IC file and setup parameters inconsistent (gas)')
 elseif (use_star .and. (npartoftype(istar) <= 0)) then
    call fatal('setup_isogal','IC file and setup parameters inconsistent (star)')
 elseif (use_bulge .and. (npartoftype(ibulge) <= 0)) then
    call fatal('setup_isogal','IC file and setup parameters inconsistent (bulge)')
 endif

 !--gas loop
 totmassG=0.
 if (use_gas .and. npartoftype(igas) > 0) then
    call add_component('asciifile_G',igas,1,npartoftype(igas),xyzh,vxyzu,&
                       massoftype,npartoftype,totmassG,thermal,hfact,totvol,'gas')
 else
    npartoftype(igas)   = 0
    massoftype(igas)    = 0.
 endif

 !--disc loop
 totmassD=0.
 if (use_star .and. npartoftype(istar) > 0) then
    call add_component('asciifile_D',istar,npartoftype(igas)+1,&
                       npartoftype(igas)+npartoftype(istar),xyzh,vxyzu,&
                       massoftype,npartoftype,totmassD,thermal,hfact,totvol,'star')
 else
    npartoftype(istar)   = 0
    massoftype(istar)    = 0.
 endif

 !--bulge loop
 totmassB=0.
 if (use_bulge .and. npartoftype(ibulge) > 0) then
    call add_component('asciifile_B',ibulge,npartoftype(istar)+npartoftype(igas)+1,&
                       npartoftype(igas)+npartoftype(istar)+npartoftype(ibulge),xyzh,vxyzu,&
                       massoftype,npartoftype,totmassB,thermal,hfact,totvolB,'bulge')
 else
    npartoftype(ibulge)   = 0
    massoftype(ibulge)    = 0.
 endif

 !--halo loop
 totmassH=0.
 if (use_halo .and. npartoftype(idarkmatter) > 0) then
    call add_component('asciifile_H',idarkmatter,npartoftype(ibulge)+npartoftype(istar)+npartoftype(igas)+1,&
                       npartoftype(igas)+npartoftype(istar)+npartoftype(ibulge)+npartoftype(idarkmatter),xyzh,vxyzu,&
                       massoftype,npartoftype,totmassH,thermal,hfact,totvolH,'halo')
 else
    npartoftype(idarkmatter)   = 0
    massoftype(idarkmatter)    = 0.
 endif

 ! finalise
 npart    = npartoftype(igas)+npartoftype(istar)+npartoftype(ibulge)+npartoftype(idarkmatter)
 totmass  = totmassD + totmassB + totmassG + totmassH
 if (id==master) then
    print*,npartoftype(igas),' gas particles set, each of mass',massoftype(igas)
    print*,npartoftype(istar),' star particles set, each of mass',massoftype(istar)
    print*,npartoftype(ibulge),' bulge particles set, each of mass',massoftype(ibulge)
    print*,npartoftype(idarkmatter),' halo particles set, each of mass',massoftype(idarkmatter)
    print*,'Total parts:',npart
    print "(5(/,a))",'Important PSA:',&
          ' If you are using a static halo you need to ensure that the halo', &
          ' potential (mass, scale length etc.) matches those of the IC generator.',&
          ' If they dont then your galaxy rotation curve will be VERY wrong.', &
          '   <--By proceeding you acknowledge this at your own peril-->'
 endif

end subroutine setup_live_stars

!-----------------------------------------------------------------------
!+
!  Read IC parameter file
!+
!-----------------------------------------------------------------------
subroutine read_ic_parameters(filename,npartoftype)
 use part, only:igas,istar,ibulge,idarkmatter
 character(len=*), intent(in)    :: filename
 integer,          intent(out)   :: npartoftype(:)
 integer :: i,lu
 character(30) :: sometext

 open(newunit=lu,file=filename,form='formatted')
 do i=1,5
    if (i==1) then
       read(lu,*) sometext,npartoftype(igas)
    elseif (i==2) then
       read(lu,*) sometext,npartoftype(istar)
    elseif (i==3) then
       read(lu,*) sometext,npartoftype(ibulge)
    elseif (i==4) then
       read(lu,*) sometext,npartoftype(idarkmatter)
    endif
 enddo
 close(lu)

end subroutine read_ic_parameters

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for galactic disc setup'
 call write_inopt(use_live_stars,'use_live_stars','use live star particles',iunit)
 if (use_live_stars) then
    call write_inopt(use_gas,'use_gas','include gas',iunit)
    call write_inopt(use_star,'use_star','include live stars',iunit)
    call write_inopt(use_bulge,'use_bulge','include live bulge',iunit)
    call write_inopt(use_halo,'use_halo','include live halo',iunit)
 else
    call write_inopt(np,'np','number of particles',iunit)
    call write_inopt(partdist,'partdist','particle distribution (r=random,o=observed)',iunit)
    if (partdist=='r') then
       call write_inopt(rcylin,'rcylin','inner radius for particle setup [kpc]',iunit)
       call write_inopt(rcyl,'rcyl','outer radius for particle setup [kpc]',iunit)
    endif
    call write_inopt(totmass,'totmass','total mass',iunit)
    call write_inopt(vset,'vset','velocity profile (r=rotation curve,f=flat,c=constant)',iunit)
    if (vset=='c') then
       call write_inopt(angvel,'angvel','velocity for rotation curve [unit distance (cm/s)]',iunit)
    endif
    if (use_dustfrac) then
       call write_inopt(dust_to_gas_ratio,'dust_to_gas_ratio','dust to gas ratio',iunit)
    endif
 endif
 call write_inopt(thermal,'thermal','initial temperature in Kelvin',iunit)
 close(iunit)

end subroutine write_setupfile

!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(use_live_stars,'use_live_stars',db,errcount=nerr)
 if (use_live_stars) then
    call read_inopt(use_gas,'use_gas',db,errcount=nerr)
    call read_inopt(use_star,'use_star',db,errcount=nerr)
    call read_inopt(use_bulge,'use_bulge',db,errcount=nerr)
    call read_inopt(use_halo,'use_halo',db,errcount=nerr)
 else
    call read_inopt(np,'np',db,min=100,max=maxp,errcount=nerr)
    call read_inopt(partdist,'partdist',db,errcount=nerr)
    if (partdist=='r') then
       call read_inopt(rcylin,'rcylin',db,min=0.,errcount=nerr)
       call read_inopt(rcyl,'rcyl',db,min=0.,errcount=nerr)
    endif
    call read_inopt(totmass,'totmass',db,min=0.,errcount=nerr)
    call read_inopt(vset,'vset',db,errcount=nerr)
    if (vset=='c') then
       call read_inopt(angvel,'angvel',db,errcount=nerr)
    endif
    if (use_dustfrac) then
       call read_inopt(dust_to_gas_ratio,'dust_to_gas_ratio',db,min=0.,errcount=nerr)
    endif
 endif
 call read_inopt(thermal,'thermal',db,min=1.,max=1000000.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Interactive setup
!+
!-----------------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt

 ! live stars vs analytic potential
 call prompt('Using live star particles?',use_live_stars)

 if (use_live_stars) then
    print*,' Does the system contain the following live components:'
    print*,' [NOTE: this will overide the read in,'
    print*,'   e.g. if you want to turn off gas]'
    call prompt('Gas?',use_gas)
    call prompt('Stars?',use_star)
    call prompt('Bulge?',use_bulge)
    call prompt('Halo?',use_halo)
 else
    call prompt('Enter number of particles ',np,1,maxp)
    call prompt('random particle placement or observed distribution (r/o)?',partdist,list=(/'r','o'/))
    if (partdist=='r') then
       call prompt('Enter inner radius for particle setup [kpc]',rcylin,0.)
       call prompt('Enter outer radius for particle setup [kpc]',rcyl,0.)
    endif
    call prompt('Enter total mass in solar masses',totmass,0.)

    print "(3(a,/),a)",'Set velocities from:', &
               ' -[r] rotation curve', &
               ' -[f] simple log-disc profile', &
               ' -[c] enter constant value'
    call prompt('Enter selected velocity profile',vset,list=(/'r','f','c'/))
    if (vset=='c') then
       call prompt(' Enter velocity for rotation curve [unit distance (cm/s)]',angvel)
    endif
    if (use_dustfrac) then
       call prompt(' What is dust to gas ratio?',dust_to_gas_ratio,0.)
    endif
 endif
 call prompt('Enter initial temperature in Kelvin',thermal,1.,1000000.)

end subroutine setup_interactive

end module setup
