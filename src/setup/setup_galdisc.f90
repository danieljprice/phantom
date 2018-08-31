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
!  this module does setup for galactic discs with live or analytic stars
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: datafiles, dim, dust, extern_spiral, externalforces, io,
!    mpiutils, options, part, physcon, prompting, random, setup_params,
!    units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

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
 use dim,            only:maxp,maxvxyzu,use_dust
 use setup_params,   only:rhozero
 use physcon,        only:Rg,pi,solarm,pc,kpc
 use units,          only:umass,udist,utime,set_units
 use mpiutils,       only:bcast_mpi
 use random,         only:ran2
 use part,           only:h2chemistry,abundance,iHI,dustfrac,istar,igas,ibulge,&
                          idarkmatter,iunknown,set_particle_type,ndusttypes
 use options,        only:iexternalforce,icooling,nfulldump,use_dustfrac
 use externalforces, only:externalforce,initialise_externalforces
 use extern_spiral,  only:LogDisc
 use io,             only:fatal
 use prompting,      only:prompt
 use set_dust,       only:set_dustfrac
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: iseed,i,itot,it,ierr
 real :: rcyl,rcylin,faclod,xmin,xmax,ymin,ymax,zmin,zmax
 real :: rmax,rcyl2,rcylin2,totmass,totvol
 real :: xmax5,ymax5,zmax5,xi,yi,zi,r2,h1,gmw
 real :: radius2,vx1,vy1,vz1,prob,thermal,disp,h2ratio
 real :: fextxi,fextyi,fextzi,poti,timei,dr,r,rand
 real(kind=8) :: uergg,angvel,vcirc,phi
 real :: tempinput,dust_to_gas_ratio
 character(100)        :: partdist,vset,galsetupic
 integer, parameter, dimension(5) :: lenAr=(/1000,3000,2670,1830,4500/)
 real,    dimension(5) :: mratios=(/0.0008,0.0788,0.2293,0.1750,0.5161/)
 integer :: mratiosi(5)
 integer :: ierrf(5)
 real :: cdf1(lenAr(1)),rp1(lenAr(1))
 real :: cdf2(lenAr(2)),rp2(lenAr(2))
 real :: cdf3(lenAr(3)),rp3(lenAr(3))
 real :: cdf4(lenAr(4)),rp4(lenAr(4))
 real :: cdf5(lenAr(5)),rp5(lenAr(5))

 real:: totmassD,totmassG,totmassB,totmassH,totvolB,totvolH
 real:: rhozero1,rhozero2,rhozero3,rhozero4,h2,h3,h4
 real:: xis,yis,zis,mis,vxis,vyis,vzis,phaseis
 character(30) :: sometext
 integer       :: use_live_stars,yn_gas,yn_star,yn_bulge,yn_halo

!
!--initialising units and flags
!
! set code units
 call set_units(dist=100.*pc,mass=1.d05*solarm,G=1.)
!
!--set input file options
!--maxvxyzu(3-4) and therefore ieos(1-2) are set in dim_galdisc
 iexternalforce = 8   !8=galactic disk potentials
 icooling       = 0   !1=cooling on, 0=off
 nfulldump      = 1

 hfact = 1.2

!
!-------------------------Setting-energies-------------------------
!
 if (maxvxyzu >= 4) then
    !Non-Isothermal PV^gamma=K
    gamma = 1.6666667
 else
    !Isothermal PV = K as T=const.
    gamma = 1.0
 endif
 time  = 0.

 !thermal default at 10000K, 1K min and 1000000K max:
 thermal = 10000.
 call prompt('Enter initial temperature in Kelvin',thermal,1.,1000000.)
 call bcast_mpi(thermal)
 tempinput=thermal
 uergg = udist**2/utime**2
 h2ratio = 0.
 gmw = (2.*h2ratio+(1.-2.*h2ratio)+0.4)/ &
       (0.1+h2ratio+(1.-2.*h2ratio))
 print*,' assuming mean molecular weight is ',gmw
 thermal = 3.0/2.0*thermal*Rg/gmw/uergg
 polyk = 2./3.*thermal
 print*,'polyk (from thermal) = ',polyk

!
!--------------------------Determine setup method--------------------------
!
 use_live_stars = 0
 yn_gas  =0
 yn_star =0
 yn_bulge=0
 yn_halo =0
 call prompt('Using live star particles [1/0]?',use_live_stars)
 call bcast_mpi(use_live_stars)
!
!-----------------Setting-positions/velocities-for-live-disc---------------
!
 if (use_live_stars==1) then
    npartoftype = 0
    massoftype = 0.0
    print*,' Does the system contain the following live components:'
    print*,' [NOTE: this will overide the read in,'
    print*,'   e.g. if you want to turn off gas]'
    call prompt('Gas? [1/0]',yn_gas)
    call bcast_mpi(yn_gas)
    call prompt('Stars? [1/0]',yn_star)
    call bcast_mpi(yn_star)
    call prompt('Bulge? [1/0]',yn_bulge)
    call bcast_mpi(yn_bulge)
    call prompt('Halo? [1/0]',yn_halo)
    call bcast_mpi(yn_halo)

    !Remember : igas=1, istar=4, idm=5, ibulge=6
    !Values for initial guesstimate of smoothing lentgths
    rcyl  = 13000.*pc/udist
    zmax  = 0.5*kpc/udist
    zmin  = -zmax
    totvol   = pi*(zmax-zmin)*(rcyl*rcyl) !disc
    totvolB  = 4./3.*pi*(0.2*rcyl)**3     !inner spheroid
    totvolH  = 4./3.*pi*(5.0*rcyl)**3     !halo

    galsetupic = 'galsetic.txt'
    OPEN(21,file=galsetupic,form='formatted')
    do i=1,5
       if (i==1) then
          read(21,*)sometext,npartoftype(igas)
       elseif (i==2) then
          read(21,*)sometext,npartoftype(istar)
       elseif (i==3) then
          read(21,*)sometext,npartoftype(ibulge)
       elseif (i==4) then
          read(21,*)sometext,npartoftype(idarkmatter)
       endif
    enddo
    close(21)
    print*,'Read in the IC parameter file: ',galsetupic
    print*,' ngas  =',npartoftype(igas)
    print*,' nstar =',npartoftype(istar)
    print*,' nbulge=',npartoftype(ibulge)
    print*,' nhalo =',npartoftype(idarkmatter)


    !--Checks:
    if ((yn_gas==1).and.(npartoftype(igas) <= 0)) then
       call fatal('setup_isogal','IC file and setup parameters inconsistent (gas)')
    elseif ((yn_star==1).and.(npartoftype(istar) <= 0)) then
       call fatal('setup_isogal','IC file and setup parameters inconsistent (star)')
    elseif ((yn_bulge==1).and.(npartoftype(ibulge) <= 0)) then
       call fatal('setup_isogal','IC file and setup parameters inconsistent (bulge)')
    endif


    !--Gas loop
    totmassG=0.
    if (yn_gas==1) then
       OPEN(24,file='asciifile_G',form='formatted')
       i=1
       over_npartG: do while(i <= npartoftype(igas))
          read(24,*)xis,yis,zis,mis,vxis,vyis,vzis,phaseis
          !--Set from above converting from SI to cgs then code units
          xyzh(1,i) = xis*100./udist
          xyzh(2,i) = yis*100./udist
          xyzh(3,i) = zis*100./udist
          vxyzu(1,i) =vxis*100./(udist/utime)
          vxyzu(2,i) =vyis*100./(udist/utime)
          vxyzu(3,i) =vzis*100./(udist/utime)
          mis = 1000.*mis/umass
          totmassG = totmassG + mis
          call set_particle_type(i,igas)
          if (maxvxyzu >= 4) then
             vxyzu(4,i) = thermal
          endif
          i=i+1
       enddo over_npartG
       CLOSE(24)
       massoftype(igas)   = totmassG/real( npartoftype(igas))
       rhozero3 = totmassG/totvol
       h3 = hfact*(massoftype(igas) /rhozero3)**(1./3.)
       !--repeat loop to set h2 guess for gas.
       print*,'Set gas with mean density',rhozero3,'smoothing length',h3
       do i=1,npartoftype(igas)
          xyzh(4,i) = h3
       enddo
    else
       npartoftype(igas)   = 0
       massoftype(igas)    = 0.
    endif


    !--Disc loop
    totmassD=0.
    if (yn_star==1) then
       OPEN(22,file='asciifile_D',form='formatted')
       i= npartoftype(igas) + 1
       over_npartS: do while(i <= npartoftype(igas) + npartoftype(istar))
          read(22,*)xis,yis,zis,mis,vxis,vyis,vzis,phaseis
          !--Set from above converting from SI to cgs the code units
          xyzh(1,i) = xis*100./udist
          xyzh(2,i) = yis*100./udist
          xyzh(3,i) = zis*100./udist
          vxyzu(1,i) =vxis*100./(udist/utime)
          vxyzu(2,i) =vyis*100./(udist/utime)
          vxyzu(3,i) =vzis*100./(udist/utime)
          mis = 1000.*mis/umass
          totmassD = totmassD + mis
          call set_particle_type(i,istar)
          if (maxvxyzu >= 4) then
             vxyzu(4,i) = thermal
          endif
          i=i+1
       enddo over_npartS
       CLOSE(22)
       massoftype(istar)  = totmassD/real( npartoftype(istar))
       rhozero1 = totmassD/totvol
       h1 = hfact*(massoftype(istar)/rhozero1)**(1./3.)
       !--repeat loop to set h1 guess for stars.
       print*,'Set stars with mean density',rhozero1,'smoothing length',h1
       do i=npartoftype(igas)+ 1,npartoftype(istar)+npartoftype(igas)
          xyzh(4,i) = h1
       enddo
    else
       npartoftype(istar)   = 0
       massoftype(istar)    = 0.
    endif


    !--Bulge loop
    totmassB=0.
    if (yn_bulge==1) then
       OPEN(23,file='asciifile_B',form='formatted')
       i=npartoftype(istar)+npartoftype(igas) + 1
       over_npartB: do while(i <=npartoftype(igas)+npartoftype(istar)+npartoftype(ibulge))
          read(23,*)xis,yis,zis,mis,vxis,vyis,vzis,phaseis
          !--Set from above converting from SI to cgs the code units
          xyzh(1,i) = xis*100./udist
          xyzh(2,i) = yis*100./udist
          xyzh(3,i) = zis*100./udist
          vxyzu(1,i) =vxis*100./(udist/utime)
          vxyzu(2,i) =vyis*100./(udist/utime)
          vxyzu(3,i) =vzis*100./(udist/utime)
          mis = 1000.*mis/umass
          totmassB = totmassB + mis
          call set_particle_type(i,ibulge)
          if (maxvxyzu >= 4) then
             vxyzu(4,i) = thermal
          endif
          i=i+1
       enddo over_npartB
       CLOSE(23)
       massoftype(ibulge)  = totmassB/real( npartoftype(ibulge))
       rhozero2 = totmassB/totvolB
       h2 = hfact*(massoftype(ibulge)/rhozero2)**(1./3.)
       !--repeat loop to set h2 guess for bulge.
       print*,'Set bulge with mean density',rhozero2,'smoothing length',h2
       do i=npartoftype(istar)+npartoftype(igas)+ 1,npartoftype(ibulge)+npartoftype(istar)+npartoftype(igas)
          xyzh(4,i) = h2
       enddo
    else
       npartoftype(ibulge)   = 0
       massoftype(ibulge)    = 0.
    endif


    !--Halo loop
    totmassH=0.
    if (yn_halo==1) then
       OPEN(23,file='asciifile_H',form='formatted')
       i=npartoftype(ibulge)+npartoftype(istar)+npartoftype(igas) + 1
       over_npartH: do while(i <=npartoftype(igas)+npartoftype(istar)+npartoftype(ibulge)+npartoftype(idarkmatter))
          read(23,*)xis,yis,zis,mis,vxis,vyis,vzis,phaseis
          !--Set from above converting from SI to cgs the code units
          xyzh(1,i) = xis*100./udist
          xyzh(2,i) = yis*100./udist
          xyzh(3,i) = zis*100./udist
          vxyzu(1,i) =vxis*100./(udist/utime)
          vxyzu(2,i) =vyis*100./(udist/utime)
          vxyzu(3,i) =vzis*100./(udist/utime)
          mis = 1000.*mis/umass
          totmassH = totmassH + mis
          call set_particle_type(i,idarkmatter)
          if (maxvxyzu >= 4) then
             vxyzu(4,i) = thermal
          endif
          i=i+1
       enddo over_npartH
       CLOSE(23)
       massoftype(idarkmatter)  = totmassH/real( npartoftype(idarkmatter))
       rhozero4 = totmassH/totvolH
       h4 = hfact*(massoftype(idarkmatter)/rhozero4)**(1./3.)
       !--repeat loop to set h2 guess for bulge.
       print*,'Set halo with mean density',rhozero4,'smoothing length',h4
       do i=npartoftype(ibulge)+npartoftype(istar)+npartoftype(igas)+ 1,&
         npartoftype(idarkmatter)+npartoftype(ibulge)+npartoftype(istar)+npartoftype(igas)
          xyzh(4,i) = h4
       enddo
    else
       npartoftype(idarkmatter)   = 0
       massoftype(idarkmatter)    = 0.
    endif


    !--Finalise
    npart    = npartoftype(igas)+npartoftype(istar)+npartoftype(ibulge)+npartoftype(idarkmatter)
    totmass  = totmassD + totmassB + totmassG + totmassH
    print*,npartoftype(igas),' gas particles set, each of mass',massoftype(igas)
    print*,npartoftype(istar),' star particles set, each of mass',massoftype(istar)
    print*,npartoftype(ibulge),' bulge particles set, each of mass',massoftype(ibulge)
    print*,npartoftype(idarkmatter),' halo particles set, each of mass',massoftype(idarkmatter)
    print*,'Total parts:',npart
    print*,'Important PSA:'
    print*,' If you are using a static halo you need to ensure that the halo'
    print*,' potential (mass, scale length etc.) matches those of the IC generator.'
    print*,' If they dont then your galaxy rotation curve will be VERY wrong.'
    print*,'   <--By proceeding you acknowledge this at your own peril-->'


!
!--------------Setting-positions/velocities-for-gas-only-disc--------------
!
 else ! no live stars
    rcylin=0.
    rcyl  =10.
    partdist = 'r'
    call prompt('random particle placement or observed distribution (r/o)?',partdist,list=(/'r','o'/),noblank=.true.)
    !print*,' random particle placement or observed distribution (r/o)'
    !read*,partdist
    if (partdist=='r') then
       call prompt('Enter inner radius for particle setup [kpc]',rcylin,0.)
       call bcast_mpi(rcylin)
       call prompt('Enter outer radius for particle setup [kpc]',rcyl,0.)
       call bcast_mpi(rcyl)
    endif
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
    totmass = 1.e4            !Where the total mass of MW is 1E4 x umass= 1E9 Mo
    !--initialise random number generator
    iseed = -6485
    print "(a,i10)",' random seed = ',iseed
    print "(3(a,f10.3),a)",' galactic disc setup... rmin = ',rcylin,' rmax = ',rcyl,' in units of ',udist/kpc,' kpc'
    xi = ran2(iseed)
    npart = maxp
    call prompt('Enter number of particles ',npart,1,maxp)
    call bcast_mpi(npart)
    if (npart > maxp) call fatal('setup','npart > maxp')
    npartoftype(1) = npart

    print "(a,es10.3,a,1pg10.3,a)",'Mass is in units of ',umass,' g (',umass/solarm,' solar masses)'
    totmass = 1.e4
    call prompt('Enter total mass ',totmass,0.)
    call bcast_mpi(totmass)
    massoftype(igas) = totmass/real(npart)

    totvol = pi*(zmax-zmin)*(rcyl2 - rcylin2)
    rhozero = totmass/totvol
    print "(a,es10.3)",' mean density = ',rhozero
    h1 = hfact*(massoftype(1)/rhozero)**(1./3.)


    if (TRIM(partdist)=='o') then
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
          mratiosi(it)=nint(mratios(it)*npart)
       enddo
       print*,'Particles in each distribution bin:', mratiosi
       print*,'Total set:', mratiosi(1)+mratiosi(2)+mratiosi(3)+mratiosi(4)+mratiosi(5)
       call cdfplacement(mratiosi(1),lenAr(1),cdf1,rp1,xyzh,itot,npart,h1)
       call cdfplacement(mratiosi(2),lenAr(2),cdf2,rp2,xyzh,itot,npart,h1)
       call cdfplacement(mratiosi(3),lenAr(3),cdf3,rp3,xyzh,itot,npart,h1)
       call cdfplacement(mratiosi(4),lenAr(4),cdf4,rp4,xyzh,itot,npart,h1)
       call cdfplacement(mratiosi(5),lenAr(5),cdf5,rp5,xyzh,itot,npart,h1)
    else
       !--Loop for purely randomly placed particles:
       xmax5 = 0.5*xmax
       ymax5 = 0.5*ymax
       zmax5 = 0.5*zmax
       i = 0
       over_npart: do while(i < npart)
          !--Initialise the particle placements randomly.
          xi = 2.0*(xmax*ran2(iseed) - xmax5)
          yi = 2.0*(ymax*ran2(iseed) - ymax5)
          zi = 2.0*(zmax*ran2(iseed) - zmax5)
          r2 = xi*xi + yi*yi
          if (r2 > rcyl2 .or. r2 < rcylin2) cycle over_npart
          i = i + 1
          if (mod(i,1000000)==0) print*,i
          if (i > maxp) stop 'error, i > maxp'
          xyzh(1,i) = xi
          xyzh(2,i) = yi
          xyzh(3,i) = zi
          xyzh(4,i) = h1
       enddo over_npart
    endif
    print*,npart,' particles set, each of mass',massoftype(1)
    !
    !------------------------Setting-velocities------------------------
    !
    ! default should be idisc = 1 and iarms=1
    print "(3(a,/),a)",'Set velocities from:', &
               ' -[r] rotation curve', &
               ' -[f] simple log-disc profile', &
               ' -[c] enter constant value'
    vset = 'r'
    call prompt('Enter selected velocity profile',vset,list=(/'r','f','c'/),noblank=.true.)
    if (vset=='c') then
       !--Set parameters for rotation curve without radial dependance
       !  e.g. angvel = -2.20e+07 is a standard value
       angvel = 2.20e7
       call prompt(' Enter velocity for rotation curve [unit distance (cm/s)]',angvel)
       call bcast_mpi(angvel)
       angvel = angvel * utime/udist
    endif
    if (use_dust) use_dustfrac = .true.
    if (use_dustfrac) then
       dust_to_gas_ratio = 0.01
       call prompt(' What is dust to gas ratio?',dust_to_gas_ratio)
       call bcast_mpi(dust_to_gas_ratio)
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

       rand = huge(rand)
       prob = 0.
       do while (prob < rand)
          vx1=40.*(ran2(iseed)-0.5)
          prob=exp(-(vx1/disp)**2./2.)
          rand = ran2(iseed)
       enddo
       vxyzu(1,i)=vxyzu(1,i)+vx1*100000.*utime/udist

       rand = huge(rand)
       prob = 0.
       do while (prob < rand)
          vy1=40.*(ran2(iseed)-0.5)
          prob=exp(-(vy1/disp)**2./2.)
          rand = ran2(iseed)
       enddo
       vxyzu(2,i)=vxyzu(2,i)+vy1*100000.*utime/udist

       rand = huge(rand)
       prob = 0.
       do while (prob < rand)
          vz1=40.*(ran2(iseed)-0.5)
          prob=exp(-(vz1/disp)**2./2.)
          rand = ran2(iseed)
       enddo
       vxyzu(3,i)=vxyzu(3,i)+vz1*100000.*utime/udist

       if (maxvxyzu >= 4) then
          vxyzu(4,i) = thermal
       endif
    enddo

 endif

!
!------------------------Setting-abundances------------------------
!
 if (h2chemistry) then
    do i=1,npart
       abundance(:,i)   = 0.
       abundance(iHI,i) = 1.  ! assume all atomic hydrogen initially
    enddo
 endif
!
!------------------------Dust-fractions----------------------------
!
 if (use_dustfrac) then
    ndusttypes = 1
    do i=1,npart
       call set_dustfrac(dust_to_gas_ratio,dustfrac(:,i))
    enddo
 endif

 return
end subroutine setpart

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!--Routine to read files containing tabulated denisty profile
!--Requires that PHANTOM_DIR is set
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
    read(lu,*,iostat=ierr)cdf,r
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

 return
end subroutine


!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!--Draw radomly from MW disc density profile tabulated in files
subroutine cdfplacement(noinbin,lenbin,cdf,r_ar,xyzh,itot,Npart,hi)
 use random,  only:ran2
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

end module setup
