!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testcoala
!
! Unit tests of the COALA dust growth library
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
 public :: test_coala

 private

contains
!----------------------------------------------------------
!+
!  Unit tests of COALA Smoluchowsky dust growth solver
!+
!----------------------------------------------------------
subroutine test_coala(ntests,npass)
 use io,        only:id,master,stdout
 use physcon,   only:solarm,au
 use units,     only:set_units
 use dim,       only:use_dust
 integer, intent(inout) :: ntests,npass

 if (.not.use_dust) then
    if (id==master) write(*,"(/,a,/)") '--> SKIPPING COALA TEST (NEED DUST=yes)'
    return
 endif
#ifdef COALA
 if (id==master) write(*,"(/,a,/)") '--> TESTING COALA DUST GROWTH'

 call set_units(mass=solarm,dist=au,G=1.d0)
 call test_dustgrowth_coala(ntests, npass)

 if (id==master) write(*,"(/,a)") '<-- COALA DUST GROWTH TEST COMPLETE'
#else
 if (id==master) write(*,"(/,a,/)") '--> SKIPPING COALA TEST (NEED COALA=yes)'
#endif

end subroutine test_coala

#ifdef COALA
!----------------------------------------------------------
!+
!  test of dust growth via smoluchowsky module
!+
!----------------------------------------------------------
subroutine test_dustgrowth_coala(ntests, npass)
 use dim,          only:maxdustsmall
 use set_dust,     only:set_dustbinfrac
 use io,           only:id,master
 use testutils,    only:checkval,update_test_scores
 use physcon,      only:micron,mm,years
 use units,        only:unit_density,udist,unit_velocity
 use part,         only:dustfrac,grainsize,graindens,ndusttypes,ndustsmall,hrho,ics,eos_vars,&
                        xyzh,vxyzu,fxyzu,fext,dustevol,deltav,npart,igas,massoftype
 use growth_coala, only:init_growth_coala,get_growth_rate_coala
 integer, intent(inout) :: ntests,npass
 integer :: i,idust,ierr,nfailed(1)
 real    :: rhoi,smin,smax,sindex,dt
 real, allocatable :: dustevol_prev(:)
 if (id==master) write(*,"(/,a)") '--> testing COALA Smoluchowsky solver'

 smin = 0.1*micron/udist
 smax = 1.0*mm/udist
 sindex  = -3.5
 i = 1
 ndustsmall = maxdustsmall
 ndusttypes = ndustsmall
 graindens = 3./unit_density   ! 3 g/cm^3
 call set_dustbinfrac(smin,smax,sindex,dustfrac(:,i),grainsize(1:ndusttypes))

 do idust=1,ndusttypes
    if (id==master) print "(a,i2,a,1pg0.3)",'bin ',idust,': eps = ',dustfrac(idust,i)
 enddo
 rhoi = 1.e-18/unit_density
 massoftype(igas) = 1.e-8

 call init_growth_coala(ierr)
 call checkval(ierr,0,0,nfailed(1),'coala initialisation')  ! check ierr = 0
 call update_test_scores(ntests,nfailed,npass)

 npart = 1
 xyzh(1:3,i) = 0.
 xyzh(4,i) = hrho(rhoi,massoftype(igas))
 vxyzu(:,i) = 0.
 deltav(:,:,i) = 0.
 do idust=1,ndusttypes
    dustevol(idust,i) = sqrt(dustfrac(idust,i)/(1.0 - dustfrac(idust,i)))
 enddo
 dustevol_prev = dustevol(:,i)  ! automatically allocates memory
 fxyzu(:,i) = 0.
 fext(:,i) = 0.
 eos_vars(ics,i) = 2.e4/unit_velocity ! 200 m/s sound speed

 ! first check that dt = 0 does not change dustevol
 dt = 0.
 call get_growth_rate_coala(npart,xyzh,vxyzu,fxyzu,fext,grainsize,dustfrac,&
                            dustevol,deltav,dt,eos_vars)
 call checkval(ndusttypes,dustevol(:,i),dustevol_prev,epsilon(0.),nfailed(1),'dustevol = dustevol_prev with dt=0')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_dustgrowth_coala
#endif

end module testcoala
