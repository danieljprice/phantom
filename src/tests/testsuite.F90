!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module test
!
! Instead of running a simulation this routine, when called,
!   initiates a series of internal tests on the code
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, io_summary, mpiutils, options, testapr,
!   testcooling, testcorotate, testdamping, testderivs, testdust, testeos,
!   testexternf, testgeometry, testgnewton, testgr, testgravity,
!   testgrowth, testindtstep, testiorig, testkdtree, testkernel, testlink,
!   testmath, testmpi, testnimhd, testpart, testpoly, testptmass,
!   testradiation, testrwdump, testsedov, testsetdisc, testsethier,
!   testsetstar, testsmol, teststep, testunits, testwind, timing
!
 implicit none
 public :: testsuite

 private

contains

subroutine testsuite(string,first,last,ntests,npass,nfail)
 use io,           only:iprint,id,master,iverbose,error
 use io_summary,   only:summary_initialise
 use testderivs,   only:test_derivs
 use teststep,     only:test_step
 use testlink,     only:test_link
 use testkdtree,   only:test_kdtree
 use testsedov,    only:test_sedov
 use testgravity,  only:test_gravity
 use testdust,     only:test_dust
 use testgrowth,   only:test_growth
 use testsmol,     only:test_smol
 use testpart,     only:test_part
 use testnimhd,    only:test_nonidealmhd
 use testapr,      only:test_apr
#ifdef FINVSQRT
 use testmath,     only:test_math
#endif
 use testkernel,   only:test_kernel
 use testptmass,   only:test_ptmass
#ifdef GR
 use testgr,       only:test_gr
#else
 use testgnewton,  only:test_gnewton
 use testcorotate, only:test_corotate
#endif
 use testexternf,  only:test_externf
 use testindtstep, only:test_indtstep
 use testrwdump,   only:test_rwdump
 use testsetdisc,  only:test_setdisc
 use testsetstar,  only:test_setstar
 use testsethier,  only:test_sethier
 use testeos,      only:test_eos
 use testcooling,  only:test_cooling
 use testgeometry, only:test_geometry
 use testwind,     only:test_wind
 use testiorig,    only:test_iorig
 use testpoly,     only:test_poly
 use testdamping,  only:test_damping
 use testradiation,only:test_radiation
 use testunits,    only:test_units
#ifdef MPI
 use testmpi,      only:test_mpi
#endif
 use timing,       only:get_timings,print_time
 use mpiutils,     only:barrier_mpi
 use dim,          only:do_radiation,use_apr
 character(len=*), intent(in)    :: string
 logical,          intent(in)    :: first,last
 integer,          intent(inout) :: ntests,npass,nfail
 logical :: testall,dolink,dokdtree,doderivs,dokernel,dostep,dorwdump,dosmol
 logical :: doptmass,dognewton,dosedov,doexternf,doindtstep,dogravity,dogeom
 logical :: dosetdisc,dosetstar,doeos,docooling,dodust,donimhd,docorotate,doany,dogrowth
 logical :: dogr,doradiation,dopart,dopoly,dompi,dohier,dodamp,dowind
 logical :: doiorig,doapr,dounits
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
       write(*,"(2x,a)") '"Trace, test and treat"'
       write(*,"(14x,a,/)") '-- South Korea'
       write(*,"(2x,a)") '"Testing a program demonstrates that it contains errors, never that it is correct"'
       write(*,"(14x,a,/)") '-- E. W. Dijkstra'
    endif
    ntests = 0
    npass  = 0
    nfail  = 0
 endif
 call get_timings(twall1,tcpu1)
 testall    = .false.
 dokernel   = .false.
 dolink     = .false.
 dopart     = .false.
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
 dosetstar  = .false.
 doeos      = .false.
 dodust     = .false.
 dogrowth   = .false.
 donimhd    = .false.
 docooling  = .false.
 dogeom     = .false.
 dogr       = .false.
 dosmol     = .false.
 doradiation = .false.
 dopoly     = .false.
 dompi      = .false.
 dohier     = .false.
 dodamp     = .false.
 dowind     = .false.
 doapr      = .false.
 doiorig    = .false.
 dounits    = .false.

 if (index(string,'deriv')     /= 0) doderivs  = .true.
 if (index(string,'grav')      /= 0) dogravity = .true.
 if (index(string,'part')      /= 0) dopart    = .true.
 if (index(string,'polytrope') /= 0) dogravity = .true.
 if (index(string,'directsum') /= 0) dogravity = .true.
 if (index(string,'dust')      /= 0) dodust    = .true.
 if (index(string,'growth')    /= 0) dogrowth  = .true.
 if (index(string,'nimhd')     /= 0) donimhd   = .true.
 if (index(string,'dump')      /= 0) dorwdump  = .true.
 if (index(string,'sink')      /= 0) doptmass  = .true.
 if (index(string,'cool')      /= 0) docooling = .true.
 if (index(string,'geom')      /= 0) dogeom    = .true.
 if (index(string,'gr')        /= 0) dogr      = .true.
 if (index(string,'smol')      /= 0) dosmol    = .true.
 if (index(string,'rad')       /= 0) doradiation = .true.
 if (index(string,'poly')      /= 0) dopoly    = .true.
 if (index(string,'mpi')       /= 0) dompi     = .true.
 if (index(string,'hier')      /= 0) dohier    = .true.
 if (index(string,'damp')      /= 0) dodamp    = .true.
 if (index(string,'wind')      /= 0) dowind    = .true.
 if (index(string,'iorig')     /= 0) doiorig   = .true.
 if (index(string,'ptmass')    /= 0) doptmass  = .true.
 if (index(string,'apr')       /= 0) doapr     = .true.
 if (index(string,'units')     /= 0) dounits   = .true.

 doany = any((/doderivs,dogravity,dodust,dogrowth,donimhd,dorwdump,&
               doptmass,docooling,dogeom,dogr,dosmol,doradiation,&
               dopart,dopoly,dohier,dodamp,dowind,doiorig,doapr,dounits/))

 select case(trim(string))
 case('kernel','kern')
    dokernel = .true.
 case('link','tree')
    dolink = .true.
 case('kdtree','revtree')
    dokdtree = .true.
 case('step')
    dostep = .true.
 case('ptmass','sink','fsi','chinchen','coin')
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
 case('setstar','star')
    dosetstar = .true.
 case('eos')
    doeos = .true.
 case('dust')
    dodust = .true.
 case('gr')
    dogr = .true.
 case('growth')
    dogrowth = .true.
 case('nimhd')
    donimhd = .true.
 case('wind')
    dowind = .true.
 case('iorig')
    doiorig = .true.
 case('mpi')
    dompi = .true.
 case('apr')
    doapr = .true.
 case('units')
    dounits = .true.
 case default
    if (.not.doany) testall = .true.
 end select
 call set_default_options_testsuite(iverbose) ! set defaults

#ifdef FINVSQRT
 call test_math(ntests,npass,usefsqrt,usefinvsqrt)
#endif

!
!--apr test
!
 if (use_apr.and.testall) then
    write(*,*) '-DAPR not currently compatible with test suite, recompile with APR=no'
    return
 elseif (use_apr.and.doapr) then
    call test_apr(ntests,npass)
 endif

!
!--test kernel module
!

 if (dokernel.or.testall) then
    call test_kernel(ntests,npass)
 endif
!
!--test of linklist/neighbour finding module
!
 if (dolink.or.testall) then
    call test_link(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of part module
!
 if (dopart .or. testall) then
    call test_part(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of derivs module
!
 if (doderivs.or.testall) then
    call test_derivs(ntests,npass,string)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of equation of state module
!
 if (doeos.or.testall) then
    call test_eos(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of equation of state module
!
 if (docooling.or.testall) then
    call test_cooling(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of kdtree module
!
 if (dokdtree.or.testall) then
    call test_kdtree(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of self-gravity
!
 if (dogravity.or.testall) then
    call test_gravity(ntests,npass,string)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of dust module
!
 if (dodust.or.testall) then
    call test_dust(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of dust growth module
!
 if (dogrowth.or.testall) then
    call test_growth(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of smoluchowsky growth solver
!
 if (dosmol.or.testall) then
    call test_smol(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of non-ideal MHD
!
 if (donimhd.or.testall) then
    call test_nonidealmhd(ntests,npass,string)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of readwrite_dumps module
!
 if (dorwdump.or.testall) then
    call test_rwdump(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of step module
!
 if (dostep.or.testall) then
    call test_step(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of ind tstep module
!
 if (doindtstep.or.testall) then
    call test_indtstep(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of ind tstep module
!
 if (dounits.or.testall) then
    call test_units(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of external forces module
!
 if (doexternf.or.testall) then
    call test_externf(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of ptmass module
!
 if (doptmass.or.testall) then
    call test_ptmass(ntests,npass,string)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif

#ifdef MPI
 if (dompi.or.testall) then
    call test_mpi(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
#endif

#ifdef GR
 if (dogr.or.testall) then
    call test_gr(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
#else
!
!--test of gnewton module
!
 if (dognewton.or.testall) then
    call test_gnewton(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of corotate module
!
 if (docorotate.or.testall) then
    call test_corotate(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
#endif
!
!--test of set_disc module
!
 if (dosetdisc.or.testall) then
    call test_setdisc(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of set_star module
!
 if (dosetstar.or.testall) then
    call test_setstar(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of set_hier module
!
 if (dohier.or.testall) then
    call test_sethier(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of geometry module
!
 if (dogeom.or.testall) then
    call test_geometry(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of polynomial solvers
!
 if (dopoly.or.testall) then
    call test_poly(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of damping module
!
 if (dodamp.or.testall) then
    call test_damping(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of radiation module
!
 if (doradiation.or.testall) then
    call test_radiation(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--test of wind module
!
 if (dowind.or.testall) then
    call test_wind(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif

!
!--test of particle id
!
 if (doiorig .or. testall) then
    call test_iorig(ntests,npass)
    call set_default_options_testsuite(iverbose) ! restore defaults
 endif
!
!--now do a "real" calculation, putting it all together (Sedov blast wave)
!
 if (dosedov.or.testall) then
    call test_sedov(ntests,npass)
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
       call system("say fantastic!")
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

!--Short subroutine to reset default options & iverbose
!  (note: iverbose is not set in set_default_options)
subroutine set_default_options_testsuite(iverbose)
 use options, only:set_default_options,iopacity_type
 integer, intent(inout) :: iverbose

 iverbose = max(iverbose,2)
 call set_default_options ! restore defaults
 iopacity_type = 0 ! do not use opacity tables in test suite

end subroutine set_default_options_testsuite

end module test
