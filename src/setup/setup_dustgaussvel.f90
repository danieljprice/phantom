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
 integer, private :: ifrag,isnow,ivrelkin
 integer, private :: iseed = -123456789
 real,    private :: deltax,polykset,dust_spread,sigma_v
 !real,    private :: grainsizecgs,graindenscgs,vfragSI,gsizemincgs,grainsizemin,vfrag,vref
 real,    private :: vfragSI,gsizemincgs
 real,    private :: grainsize(1),graindens(1)
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

 integer         :: maxp,maxvxyzu,ndustx
 integer(kind=4) :: i,ipart!,ind_x,ind_y,ind_z
 integer(kind=4) :: ngas,ndust
 real            :: totmass
 !real            :: x0,y0,z0,dx,dy,dz
 real            :: mprev(npart)
 real            :: filfacprev(npart)

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

!
!--dust growth - needed to get vrel
!
 if (use_dustgrowth) then
    ifrag       = 0
    ivrelkin    = 0
    isnow       = 0
    vfragSI     = 15.
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
 ndustx = 16!64  !used npartx*10 for analysis at timestep 0
 deltax = dxbound/ndustx
 ndust = 0
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                 hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 ndust = int(npart_total,kind=4) - ngas
 !print*,npart_total,ngas,ndust
 ! cluster centre = box centre
! x0 = 0.9*xmin!0.5*(xmin+xmax)
! y0 = 0.9*xmin!0.5*(ymin+ymax)
! z0 = 0.9*xmin!0.5*(zmin+zmax)
!
! ! offsets
! dx = 0.9*(xmax-xmin) / ndust!**(1/3)
! dy = 0.9*(xmax-xmin) / ndust!**(1/3)
! dz = 0.9*(xmax-xmin) / ndust!**(1/3)
!
! ind_x = 0
! ind_y = 0
! ind_z = 0
 do i=1,ndust

   ipart = ngas + i
   call set_particle_type(ipart,idust)

   ! get inds
!   if (MOD(i,3)==0) then
!      ind_x = ind_x + 1
!   elseif (MOD(i,3)==1) then
!      ind_y = ind_y + 1
!   else
!      ind_z = ind_z + 1
!   endif

   !dx = dust_spread*gauss_random(iseed)
   !dy = dust_spread*gauss_random(iseed)
   !dz = dust_spread*gauss_random(iseed)
   !if (mod(ipart,2)==1) dx=-dx

   ! all dust at same location + tiny offsets
   !xyzh(1,ipart) = x0 + dx*ind_x
   !xyzh(2,ipart) = y0 + dy*ind_y
   !xyzh(3,ipart) = z0 + dz*ind_z

   ! smoothing length
   !xyzh(4,ipart) = 2. !xyzh(4,1) ! articficial to bypass warnings

   !call set_particle_type(ipart,idust)

   ! Gaussian velocity
   vxyzu(1,ipart) = sigma_v * gauss_random(iseed) !+ 3.0
   vxyzu(2,ipart) = sigma_v * gauss_random(iseed) !+ 3.0
   vxyzu(3,ipart) = sigma_v * gauss_random(iseed) !+ 3.0


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

! do i=1,npart
!    if (xyzh(4,i)==0.) then
!       call remove_particle_from_npartoftype(i,npoftype) !kill extra particles, to be coded properly
!    endif
! enddo

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
