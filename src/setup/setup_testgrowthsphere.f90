!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Dust + gas setup:
! - gas particles on uniform lattice, v = 0
! - dust particles initially at same location (tiny offsets between them)
! - dust velocities from Gaussian distribution
! - explicit Gaussian RNG (Box–Muller)
!
 use setup_params, only:rhozero,npart_total
 use dim,          only:use_dustgrowth
 use units,        only:udist,unit_density,unit_velocity
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 integer, private :: ifrag,isnow
 integer, private :: iseed = -123456789
 real,    private :: deltax,polykset,dust_spread,sigma_v
 real,    private :: vfragSI,gsizemincgs
 real,    private :: grainsize(1),graindens(1)
! real,    private :: grainsizemin,vfrag,vref,grainsizecgs,graindenscgs
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform gas particle distribution + dust particles all stacked together
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu, &
                   polyk,gamma,hfact,time,fileprefix)

 use io,       only:master
 use boundary, only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound
 use unifdis,  only:set_unifdis
 use part,     only:set_particle_type,igas,idust,periodic,&
                    dustprop,dustgasprop,VrelVf,&
                    filfac,probastick,&
                    iamdust,&
                    kill_particle
 use units,    only:set_units
 use physcon,  only:pi,solarm,au,fourpi
 use random,   only:gauss_random
 use mpidomain,only:i_belong
 use dust,     only:idrag

 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:),vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 integer :: maxp,maxvxyzu,ndustx
 integer :: i,ngas,ndust,ipart,icompt
!integer :: ind_x,ind_y,ind_z
 real :: totmass
 !real :: x0,y0,z0,dx,dy,dz
 real    :: mprev(npart)
 real    :: filfacprev(npart)
 real    :: rmin,rmax
 real    :: x,y,z,r,rcompt

 !-----------------------
 ! Units and parameters
 !-----------------------
 call set_units(mass=solarm,dist=au,G=1.d0)

 time   = 0.
 gamma  = 1.
 hfact  = 1.2
 rhozero = 1.
 polykset = 1.

 npartx = 16
 ilattice = 1
 dust_spread = 1e-2   ! max offset from center
 sigma_v = 1.       ! velocity dispersion

 polyk  = polykset**2
 deltax = dxbound/npartx

 idrag = 0 ! no drag

 rmin = 0.1 !min sphere radius
 rmax = 0.40 !max sphere radius

!
!--dust growth - needed to get vrel
!
 if (use_dustgrowth) then
    ifrag = 0
    isnow = 0
    vfragSI = 15.
    gsizemincgs = 5.e-3

    grainsize = 0.1 !random value
    graindens = 3. !random value

    mprev(:) = 99.
    filfacprev(:) = 99.
 endif

! setup particles
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 npart = 0
 npartoftype(:) = 0
 npart_total = 0

 if (id==master) print *, ">>> Dust cluster + gas setup <<<"

 !--------------------------------------------------
 ! 1. GAS particles (uniform lattice)
 !--------------------------------------------------
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                  hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)

 ngas = int(npart_total,kind=4)
 do i=1,ngas
   call set_particle_type(i,igas)

   ! gas at rest
   vxyzu(1,i) = 0.
   vxyzu(2,i) = 0.
   vxyzu(3,i) = 0.

   if (xyzh(1,i)==0 .and. xyzh(2,i)==0 .and. xyzh(3,i)==0) then !assign random pos in case of stacking at origin
      xyzh(1,i) = gauss_random(iseed)*0.01
      xyzh(2,i) = gauss_random(iseed)*0.01
      xyzh(3,i) = gauss_random(iseed)*0.01
   endif
 enddo

 !--------------------------------------------------
 ! 2. DUST particles (clustered at one location)
 !--------------------------------------------------
 ndustx = 128  !used npartx*10 for analysis at timestep 0
 deltax = dxbound/ndustx
 ndust = 0
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                 hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong, rmin=rmin,rmax=rmax)
 ndust = int(npart_total,kind=4) - ngas

 icompt = 0
 rcompt = 1e70
 do i=1,ndust

   ipart = ngas + i
   call set_particle_type(ipart,idust)

   ! Velocity pointing towards origin with magnitude of 1
   x = xyzh(1,ipart)
   y = xyzh(2,ipart)
   z = xyzh(3,ipart)
   r = sqrt( x**2 + y**2 + z**2 )
   vxyzu(1,ipart) = -x/r
   vxyzu(2,ipart) = -y/r
   vxyzu(3,ipart) = -z/r

   ! Look for particle at origin
   if (r<rcompt) then
      rcompt = r
      icompt = ipart
   endif

   !--set dustprops
   if (use_dustgrowth) then
      !dustprop(:,i) = 0.
      dustprop(1,ipart) = fourpi/3.*graindens(1)*grainsize(1)**3
      dustprop(2,ipart) = graindens(1)
      filfac(ipart) = 0.
      probastick(ipart) = 1.
      dustgasprop(:,ipart) = 0.
      VrelVf(:,ipart)        = 0.
   endif
 enddo


 !npart = ngas + ndust

 ! Particle at origin with position and velocity of 0
 icompt = int(npart_total,kind=4)+1
 xyzh(1,icompt)  = 0.
 xyzh(2,icompt)  = 0.
 xyzh(3,icompt)  = 0.
 xyzh(4,icompt)  = 0.5  ! size of the box
 vxyzu(1,icompt) = 0.
 vxyzu(2,icompt) = 0.
 vxyzu(3,icompt) = 0.
 call set_particle_type(icompt,idust)
 dustprop(1,icompt) = fourpi/3.*graindens(1)*grainsize(1)**3
 dustprop(2,icompt) = graindens(1)
 filfac(icompt) = 0.
 probastick(icompt) = 1.
 dustgasprop(:,icompt) = 0.
 VrelVf(:,icompt)        = 0.
 npart_total = npart_total+1
 ndust = ndust+1
 npart = npart+1


 !--------------------------------------------------
 ! masses
 !--------------------------------------------------
 totmass = rhozero * (xmax-xmin)*(ymax-ymin)*(zmax-zmin)

 massoftype(igas)  = totmass/(2.*ngas)
 massoftype(idust) = totmass/(2.*ndust)

 npartoftype(igas)  = ngas
 npartoftype(idust) = ndust

end subroutine setpart

end module setup
