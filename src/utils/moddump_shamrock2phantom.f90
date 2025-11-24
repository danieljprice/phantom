!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! default moddump routine: does not make any modifications
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos,     only:ieos,polyk,gamma,qfacdisc
 use units,   only:set_units
 use physcon, only:solarm,au,gg,c
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: Mbh, Rg, R_in, HonR_in, omega_k, cs_in

 ! Shamrock units are Msun and au, and time units are seconds
 call set_units(mass=solarm,dist=au,time=1.d0)
 ieos = 3
 gamma = 1.
 polyk = 1.

 qfacdisc = 0.75

 Mbh = 1.e6*solarm
 Rg = 4.*gg*Mbh/c**2
 R_in = 4.*Rg
 HonR_in = 0.01

 omega_k = sqrt(gg*Mbh/R_in**3)
 cs_in = HonR_in*R_in*omega_k
 polyk = cs_in**2

 massoftype(1) = 1.e-10
 npartoftype(1) = npart
 print*,' GOT npart = ',npart
 print*,' GOT npartoftype = ',npartoftype

end subroutine modify_dump

end module moddump

