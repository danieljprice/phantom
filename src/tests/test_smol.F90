!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testsmol
!
! Unit tests of the equation of state module
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, growth_smol, io, part, physcon, set_dust, testutils,
!   units
!
 implicit none
 public :: test_smol

 private

contains
!----------------------------------------------------------
!+
!  Unit tests of smoluchowsky dust growth solver
!+
!----------------------------------------------------------
subroutine test_smol(ntests,npass)
 use io,        only:id,master,stdout
 use physcon,   only:solarm,au
 use units,     only:set_units
 use dim,       only:use_dust
 integer, intent(inout) :: ntests,npass

 if (.not.use_dust) then
    if (id==master) write(*,"(/,a,/)") '--> SKIPPING SMOLUCHOWSKY TEST (NEED DUST=yes)'
    return
 endif
#ifdef SMOL
 if (id==master) write(*,"(/,a,/)") '--> TESTING SMOLUCHOWSKY DUST GROWTH'

 call set_units(mass=solarm,dist=au,G=1.d0)
 call test_dustgrowth_smol(ntests, npass)

 if (id==master) write(*,"(/,a)") '<-- SMOLUCHOWSKY DUST GROWTH TEST COMPLETE'
#else
 if (id==master) write(*,"(/,a,/)") '--> SKIPPING SMOLUCHOWSKY TEST (NEED SMOL=yes)'
#endif

end subroutine test_smol

#ifdef SMOL
!----------------------------------------------------------
!+
!  test of dust growth via smoluchowsky module
!+
!----------------------------------------------------------
subroutine test_dustgrowth_smol(ntests, npass)
 use set_dust,  only:set_dustbinfrac
 use io,        only:id,master
 use testutils, only:checkval,update_test_scores
 use physcon,   only:micron,mm,years
 use units,     only:unit_density,utime
 use part,      only:dustfrac,grainsize,graindens,ndusttypes
 use growth_smol, only:grain_growth_smol
 integer, intent(inout) :: ntests,npass
 integer :: i,nfailed(1)
 real    :: rhoi,smincgs,smaxcgs,sindex,dt
 if (id==master) write(*,"(/,a)") '--> testing smoluchowsky solver'

 smincgs = 0.1*micron
 smaxcgs = 1.0*mm
 sindex  = -3.5
 i = 1
 ndusttypes = size(grainsize)
 graindens = 3./unit_density   ! 3 g/cm^3
 call set_dustbinfrac(smincgs,smaxcgs,sindex,dustfrac(:,i),grainsize(1:ndusttypes))

 print*,' dustfrac = ',dustfrac(:,i)
 rhoi = 1.e-18/unit_density
 dt = 1.*years/utime
 call grain_growth_smol(ndusttypes,dustfrac(:,i),rhoi,grainsize,graindens,dt)

 call update_test_scores(ntests,nfailed,npass)

end subroutine test_dustgrowth_smol
#endif

end module testsmol
