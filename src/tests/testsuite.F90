!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: test
!
!  DESCRIPTION:
!   Instead of running a simulation this routine, when called,
!   initiates a series of internal tests on the code
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, io_summary, mpiutils, options, testcooling,
!    testcorotate, testderivs, testdust, testeos, testexternf,
!    testgeometry, testgnewton, testgravity, testgrowth, testindtstep,
!    testkdtree, testkernel, testlink, testmath, testnimhd, testptmass,
!    testrwdump, testsedov, testsetdisc, teststep, timing
!+
!--------------------------------------------------------------------------
module test
 implicit none
 public :: testsuite

 private

contains

subroutine testsuite(string,first,last,ntests,npass,nfail)
 use io,           only:iprint,id,master,iverbose
 use io_summary,   only:summary_initialise
 use testderivs,   only:test_derivs
 use teststep,     only:test_step
 use testlink,     only:test_link
 use testkdtree,   only:test_kdtree
 use testsedov,    only:test_sedov
 use testgravity,  only:test_gravity
 use testdust,     only:test_dust
 use testgrowth,   only:test_growth
 use testnimhd,    only:test_nonidealmhd
#ifdef FINVSQRT
 use testmath,     only:test_math
#endif
 use testkernel,   only:test_kernel
 use testptmass,   only:test_ptmass
 use testgnewton,  only:test_gnewton
 use testcorotate, only:test_corotate
 use testexternf,  only:test_externf
 use testindtstep, only:test_indtstep
 use testrwdump,   only:test_rwdump
 use testsetdisc,  only:test_setdisc
 use testeos,      only:test_eos
 use testcooling,  only:test_cooling
 use testgeometry, only:test_geometry
 use options,      only:set_default_options
 use timing,       only:get_timings,print_time
 use mpiutils,     only:barrier_mpi
 character(len=*), intent(in)    :: string
 logical,          intent(in)    :: first,last
 integer,          intent(inout) :: ntests,npass,nfail
 logical :: testall,dolink,dokdtree,doderivs,dokernel,dostep,dorwdump
 logical :: doptmass,dognewton,dosedov,doexternf,doindtstep,dogravity,dogeom
 logical :: dosetdisc,doeos,docooling,dodust,donimhd,docorotate,doany,dogrowth
#ifdef FINVSQRT
 logical :: usefsqrt,usefinvsqrt
#endif
 real(kind=4) :: twall1,tcpu1,twall2,tcpu2

 call summary_initialise

 iprint = 6
 iverbose = max(iverbose,2)
 if (first) then
    if (id==master) then
       write(*,"(/,a,/)") '--> RUNNING PHANTOM TEST SUITE'
       write(*,"(2x,a)") '"Nobody cares how fast you can calculate the wrong answer."'
       write(*,"(14x,a,/)") '-- Richard West (former UKAFF manager)'
    endif
    ntests = 0
    npass  = 0
    nfail  = 0
 endif
 call get_timings(twall1,tcpu1)
 testall    = .false.
 dokernel   = .false.
 dolink     = .false.
 dokdtree   = .false.
 doderivs   = .false.
 dostep     = .false.
 doptmass   = .false.
 dognewton  = .false.
 docorotate = .false.
 dosedov    = .false.
 doexternf  = .false.
 doindtstep = .false.
 dogravity  = .false.
 dorwdump   = .false.
 dosetdisc  = .false.
 doeos      = .false.
 dodust     = .false.
 dogrowth   = .false.
 donimhd    = .false.
 docooling  = .false.
 dogeom     = .false.
 if (index(string,'deriv')     /= 0) doderivs  = .true.
 if (index(string,'grav')      /= 0) dogravity = .true.
 if (index(string,'polytrope') /= 0) dogravity = .true.
 if (index(string,'directsum') /= 0) dogravity = .true.
 if (index(string,'dust')      /= 0) dodust    = .true.
 if (index(string,'growth')    /= 0) dogrowth  = .true.
 if (index(string,'nimhd')     /= 0) donimhd   = .true.
 if (index(string,'dump')      /= 0) dorwdump  = .true.
 if (index(string,'sink')      /= 0) doptmass  = .true.
 if (index(string,'cool')      /= 0) docooling = .true.
 if (index(string,'geom')      /= 0) dogeom    = .true.
 doany = any((/doderivs,dogravity,dodust,dogrowth,donimhd,dorwdump,doptmass,docooling,dogeom/))

 select case(trim(string))
 case('kernel','kern')
    dokernel = .true.
 case('link','tree')
    dolink = .true.
 case('kdtree','revtree')
    dokdtree = .true.
 case('step')
    dostep = .true.
 case('ptmass','sink')
    doptmass = .true.
 case('gnewton')
    dognewton = .true.
 case('corotate','corot','coriolis','centrifugal')
    docorotate = .true.
 case('externf','extern','extf')
    doexternf = .true.
 case('sedov','blast')
    dosedov = .true.
 case('indtstep','ind')
    doindtstep = .true.
 case('gravity','grav')
    dogravity = .true.
 case('dump','rwdump','dumprw')
    dorwdump = .true.
 case('setdisc','disc')
    dosetdisc = .true.
 case('eos')
    doeos = .true.
 case('dust')
    dodust = .true.
 case('growth')
    dogrowth = .true.
 case('nimhd')
    donimhd = .true.
 case default
    if (.not.doany) testall = .true.
 end select
#ifdef FINVSQRT
 call test_math(ntests,npass,usefsqrt,usefinvsqrt)
 call barrier_mpi()
#endif
!
!--test kernel module
!
 if (dokernel.or.testall) then
    call test_kernel(ntests,npass)
    call barrier_mpi()
 endif
!
!--test of linklist/neighbour finding module
!
 if (dolink.or.testall) then
    call test_link(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of derivs module
!
 if (doderivs.or.testall) then
    call test_derivs(ntests,npass,string)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of equation of state module
!
 if (doeos.or.testall) then
    call test_eos(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of equation of state module
!
 if (docooling.or.testall) then
    call test_cooling(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of kdtree module
!
 if (dokdtree.or.testall) then
    call test_kdtree(ntests,npass)
    call set_default_options
    call barrier_mpi()
 endif
!
!--test of self-gravity
!
 if (dogravity.or.testall) then
    call test_gravity(ntests,npass,string)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of dust module
!
 if (dodust.or.testall) then
    call test_dust(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of dust growth module
!
 if (dogrowth.or.testall) then
    call test_growth(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of non-ideal MHD
!
 if (donimhd.or.testall) then
    call test_nonidealmhd(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of readwrite_dumps module
!
 if (dorwdump.or.testall) then
    call test_rwdump(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of step module
!
 if (dostep.or.testall) then
    call test_step(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of ind tstep module
!
 if (doindtstep.or.testall) then
    call test_indtstep(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of external forces module
!
 if (doexternf.or.testall) then
    call test_externf(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of ptmass module
!
 if (doptmass.or.testall) then
    call test_ptmass(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of gnewton module
!
 if (dognewton.or.testall) then
    call test_gnewton(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of corotate module
!
 if (docorotate.or.testall) then
    call test_corotate(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of set_disc module
!
 if (dosetdisc.or.testall) then
    call test_setdisc(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--test of geometry module
!
 if (dogeom.or.testall) then
    call test_geometry(ntests,npass)
    call set_default_options ! restore defaults
    call barrier_mpi()
 endif
!
!--now do a "real" calculation, putting it all together (Sedov blast wave)
!
 if (dosedov.or.testall) then
    call test_sedov(ntests,npass)
    call barrier_mpi()
 endif

 if (last .and. id==master) then
    nfail = max(ntests - npass,0)
    write(*,"(/,a)") '<-- testing complete'
    call get_timings(twall2,tcpu2)
    call print_time(twall2-twall1,'total wall time = ')
    call print_time(tcpu2-tcpu1,  'total cpu time  = ')
    if (ntests > 0) then
       write(*,"(/,a,2(/,a,i3,a,i3,1x,2pf5.1,'%'),/)") &
        'SUMMARY OF ALL TESTS:', &
        'PASSED: ',npass,' of ',ntests,npass/real(ntests), &
        'FAILED: ',nfail,' of ',ntests,nfail/real(ntests)
    endif

    if (nfail==0) then
       write(*,"(5(a,/))") &
          " ____   _    ____ ____  ", &
          "|  _ \ / \  / ___/ ___| ", &
          "| |_) / _ \ \___ \___ \ ", &
          "|  __/ ___ \ ___) |__) |", &
          "|_| /_/   \_\____/____/ "

       write(*,"(a)") 'TEST SUITE PASSED'
       call system("say OK")
    else
       write(*,"(5(a,/))") &
          " _____ _    ___ _     ", &
          "|  ___/ \  |_ _| |    ", &
          "| |_ / _ \  | || |    ", &
          "|  _/ ___ \ | || |___ ", &
          "|_|/_/   \_\___|_____|"
       write(*,"(a)") 'TEST SUITE FAILED'
       call system("say fail")
    endif
 endif

end subroutine testsuite

end module test
