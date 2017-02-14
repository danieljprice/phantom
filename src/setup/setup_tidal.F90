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
!  OWNER: Clément Bonnerot
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: physcon, setup_params, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for tidal debris
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:npart_total
 use units,          only:set_units
 use physcon,        only:solarr,solarm
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: totmass,mass
 integer :: i,maxp,maxvxyzu,num

 call set_units(dist=solarr,mass=solarm,G=1.)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 5./3.
!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))

 num = 0
 open (unit=12,form='formatted',file='input',status='old')
 do
    read (12,*, END=10)
    num = num + 1 ! count the number of particles
 enddo
 10 continue
 rewind(12)

 npart = num
 npart_total = num
 npartoftype(:) = 0
 npartoftype(1) = npart

 print*,' Setting up ',npart,' particles...'

 do i=1,npart
     read (12,*) xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),mass,xyzh(4,i),vxyzu(4,i)
 enddo

 massoftype(:) = 0
 massoftype(1) = mass
 totmass = npart*massoftype(1)

 print*,' mass of the particles = ',massoftype(1)
 print*,' total mass = ',totmass
 print*,' adiabatic exponent = ',gamma

end subroutine setpart

end module setup

