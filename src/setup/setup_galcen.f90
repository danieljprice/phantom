!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Alex Rimoldi
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, io, options, part, physcon, setup_params,
!    spherical, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for a galactic centre model
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master,fatal
 use spherical,    only:set_sphere
 use options,      only:ieos,iexternalforce
 use timestep,     only:tmax
 use centreofmass, only:reset_centreofmass
 use units,        only:udist,umass,utime,set_units
 use physcon,      only:pc,solarm,gg,years,pi
 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 ! real, allocatable,    dimension(:)   :: rnew
 real    :: psep,bgmass,snmass,totmass
 real    :: bgvol,bgrmax,bgrmin
 real    :: snvol,snrmax,snrmin
 real    :: ufac,uval,mbh,ro,accretion_radius,gcode
 real    :: enblast,spsoundsn
 integer :: i,np,nx,omega,ilattice,maxvxyzu,i1 !,iseed
 integer :: npinsn,npinbg
 integer, parameter :: ng  = 1024
 real(kind=8) :: uergg
 real :: r(ng) !,den
 ! logical :: dostretch

 call set_units(dist=pc,mass=solarm,G=1.)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 5./3.
 bgrmin  = 0. !0.001*pc/udist
 bgrmax  = 1.*pc/udist
 snrmin = 0.0
 snrmax = 1.0d-3*pc/udist
 bgmass = 10.*solarm/umass
 snmass = 1.*solarm/umass
 totmass = bgmass+snmass
 mbh = 4.0d6*solarm/umass
 iexternalforce = 0
 ilattice = 2
 ieos = 2
 tmax = 1.
!
!--setup particles
!
 np = 10000 !size(xyzh(1,:))
 npinsn = np/8. ! just to test
 npinbg = np-npinsn
 print *,'initial npinbg, npinsn = ',npinbg,npinsn

 maxvxyzu = size(vxyzu(:,1))
 bgvol = 4./3.*pi*bgrmax**3
 nx = int(npinbg**(1./3.))
 psep = bgvol**(1./3.)/real(nx)
 print*,'Setup for background gas: '
 print*,' background volume = ',bgvol,' particle separation = ',psep
 print*,' maximum background radius = ',bgrmax
 print*,' bgvol/psep**3 = ',bgvol/psep**3
 if (maxvxyzu < 4) then
    call fatal('setup','Setup requires thermal energy to be stored')
 endif
!
!--add stretched sphere for background
!
 npart = 0
 npart_total = 0
 call set_sphere('closepacked',id,master,bgrmin,bgrmax,psep,hfact,npart,xyzh,rhofunc=rhofunc,nptot=npart_total)
 npinbg = npart_total
!
!--set particle properties
!
 rhozero = bgmass/bgvol
 polyk   = 0.
 print*,' background mass = ',bgmass,' mean density = ',rhozero,' polyk = ',polyk
 print*,' free fall time = ',sqrt(3.*pi/(32.*rhozero))

 massoftype(igas) = bgmass/npinbg ! change this: scale to sn mass
 print*,' particle mass = ',massoftype(igas)

 gcode = gg*umass*utime**2/udist**3
 print *,'gcode = ',gcode
 ufac = gcode*mbh/((gamma-1.)*(omega+1.))/bgmass
 ! allocate(rnew(npart))
 do i=1,npart
    ! ro = sqrt(sum(xyzh(1:3,i)**2))
    vxyzu(4,i) = 0. !ufac*(1./ro - (ro**omega)/(rmax**(omega+1)))
 enddo
!
!--set sink properties
!
 nptmass = 1
 accretion_radius = 0.001*pc/udist
 xyzmh_ptmass(:,nptmass) = 0.
 xyzmh_ptmass(4,nptmass) = mbh
 xyzmh_ptmass(ihacc,nptmass) = accretion_radius
 xyzmh_ptmass(ihsoft,i1) = 0.0
 vxyz_ptmass(:,:) = 0.
!
!--reset centre of mass to the origin
!
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
!
!--now add an (offset) explosion
!

 snvol = 4./3.*pi*snrmax**3
 nx = int(npinsn**(1./3.))
 psep = snvol**(1./3.)/real(nx)
 print*,' Setup for blast wave: '
 print*,' initial volume = ',snvol,' particle separation = ',psep
 print*,' initial maximum radius = ',snrmax
 print*,' snvol/psep**3 = ',snvol/psep**3

 uergg = umass*udist**2/utime**2
 print *,'uergg = ',uergg
 enblast = 1.0e51/uergg

 call set_supernova('closepacked',id,master,snrmin,snrmax,psep,hfact,npart,xyzh,vxyzu,massoftype,enblast,npart_total)

 npinsn = npart_total-npinbg
 npartoftype(:) = 0
 npartoftype(igas) = npart_total

 ! sound speed = sqrt(gamma*P/rho) = sqrt(gamma*(gamma-1)*u)
 spsoundsn = sqrt(gamma*(gamma-1.)*enblast/(npinsn*massoftype(igas)))
 if (id==master) then
    print *,'enblast/(npinsn*massoftype(igas))',enblast/(npinsn*massoftype(igas))
    print *,'supernova sound crossing time = ',snrmax/spsoundsn
 endif

end subroutine setpart


real function rhofunc(r)
 real, intent(in) :: r

 rhofunc = 1./r

end function rhofunc


subroutine set_supernova(lattice,id,master,rmin,rmax,delta,hfact,np,xyzh,vxyzu,massoftype,enblast,nptot)

 use physcon,   only:pc
 use units,     only:udist
 use io,        only:error,fatal
 use spherical, only:set_sphere
 use part,      only:igas
 character(len=*), intent(in)    :: lattice
 integer,          intent(in)    :: id,master
 integer,          intent(inout) :: np
 real,             intent(in)    :: rmin,rmax,hfact,enblast
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(out)   :: vxyzu(:,:)
 real,             intent(in)    :: massoftype(:)
 real,             intent(inout) :: delta
 integer(kind=8),  intent(inout) :: nptot
 real :: xyz_origin(3)
 real               :: massinsn
 integer            :: i,npin
 integer(kind=8)    :: npintot

 npin = np
 npintot = nptot

 xyz_origin = (/0.25*pc/udist,0.,0./)

 call set_sphere(lattice,id,master,rmin,rmax,delta,hfact,np,xyzh,xyz_origin=xyz_origin,rhofunc=rhosnfunc,nptot=nptot)
 massinsn = massoftype(igas)*(nptot-npintot)
 if (id==master) then
    print *,'Number of particles in supernova: ',nptot-npintot,' mass = ',massinsn,' enblast/massinsn = ',enblast/massinsn
 endif
 do i = npin+1,np
    vxyzu(1:3,i) = 0.
    vxyzu(4,i) = enblast/massinsn
 enddo

end subroutine set_supernova


real function rhosnfunc(r)
 real, intent(in) :: r

 rhosnfunc = 1.

end function rhosnfunc


end module setup

