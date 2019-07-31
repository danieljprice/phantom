!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testsetdisc
!
!  DESCRIPTION:
!  Unit tests of set_disc module
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: checksetup, deriv, dim, eos, io, options, part, setdisc,
!    testutils, timing, units
!+
!--------------------------------------------------------------------------
module testsetdisc
 implicit none
 public :: test_setdisc

 private

contains

subroutine test_setdisc(ntests,npass)
 use dim,        only:maxp
 use io,         only:id,master
 use part,       only:npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu,fext,Bevol,mhd, &
                      alphaind,maxalpha, &
                      divcurlv,divcurlB,dBevol,periodic,maxvxyzu,dustfrac,ddustevol,dustprop,ddustprop,temperature,&
                      pxyzu,dens,metrics
 use eos,        only:polyk,gamma
 use options,    only:ieos,alpha,alphau,alphaB
 use testutils,  only:checkval,checkvalf,checkvalbuf_start,checkvalbuf,checkvalbuf_end,update_test_scores
 use deriv,      only:derivs
 use timing,     only:getused,printused
 use setdisc,    only:set_disc
 use checksetup, only:check_setup
 use units,      only:set_units
 integer, intent(inout) :: ntests,npass
 integer :: nparttot
 integer :: nfailed(3),ncheck
 integer :: i,nerr,nwarn
 real :: time,dtext_dum
 real(kind=4) :: t1,t2
 logical :: testall
 real :: runit(3),xi(3)
 real :: rcyl,fr,vphi,sum,errmax

 if (periodic) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING DISC SETUP TESTS'
    return
 else
    if (id==master) write(*,"(/,a,/)") '--> TESTING DISC SETUP'
 endif

 testall  = .true.
 call set_units(mass=1.d0,dist=1.d0,G=1.d0)
!
!--test that centrifugal acceleration balances radial forces
!
 testvelocities: if (testall) then
    if (id==master) write(*,"(/,a)") '--> testing set_disc'
    !
    !  Set problem parameters
    !
    nparttot= min(size(xyzh(1,:)),1000000)
    gamma   = 1.0
    time    = 0.
    hfact = 1.2
    ieos = 3
    !
    ! set up the disc
    !
    call set_disc(id,master=master,&
                   nparttot = nparttot,&
                   npart   = npart,&
                   rmin    = 0.5, &
                   rmax    = 10.,&
                   p_index = 1.5,    &
                   q_index = 0.75,   &
                   HoverR  = 0.02, &
                   disc_mass = 1.e-4,   &
                   star_mass = 1.0,    &
                   gamma   = gamma,  &
                   particle_mass = massoftype(1), &
                   hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,polyk=polyk,&
                   ismooth=.true.,writefile=.false.)
    npartoftype(:) = 0
    npartoftype(1) = npart
    if (mhd) Bevol(:,:) = 0.
!
!--make sure AV is off
!
    alpha  = 0.
    alphau = 0.
    alphaB = 0.
    if (maxalpha==maxp)  alphaind = 0.
!
!--check that set_disc passes check_setup routine
!
    call check_setup(nerr,nwarn)
    call checkval(nerr,0,0,nfailed(1),'setup errors')
    call update_test_scores(ntests,nfailed(1:1),npass)

    call checkval(nwarn,0,0,nfailed(1),'setup warnings')
    call update_test_scores(ntests,nfailed(1:1),npass)
!
!--calculate derivatives
!
    call getused(t1)
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum,pxyzu,dens,metrics)
    call getused(t2)
    if (id==master) call printused(t1)
!
!--check force balance
!
    errmax = 0.
    nfailed = 0
    ncheck = 0
!    call checkvalbuf_start('centrifugal balance vphi**2/r = force_r')
    do i=1,npart
       xi    = xyzh(1:3,i)
       rcyl  = sqrt(dot_product(xi(1:2),xi(1:2)))
       runit = (/xi(1)/rcyl,xi(2)/rcyl,0./)
       fr    = dot_product(fxyzu(1:3,i),runit) - 1./rcyl**2
       vphi  = vxyzu(1,i)*(-xi(2)/rcyl) + vxyzu(2,i)*(xi(1)/rcyl)
       sum   = (vphi**2)/rcyl + fr
       call checkvalbuf(sum,0.,3.1e-1,'centrifugal balance vphi**2/r = force_r',nfailed(1),ncheck,errmax)
    enddo
    call checkvalbuf_end('vphi**2/r = force_r',ncheck,nfailed(1),errmax,3.1e-1)
    call update_test_scores(ntests,nfailed(1:1),npass)

 endif testvelocities

 if (id==master) write(*,"(/,a)") '<-- DISC SETUP TESTS COMPLETE'

end subroutine test_setdisc

end module testsetdisc
