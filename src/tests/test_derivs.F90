!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testderivs
!
!  DESCRIPTION:
!   Unit test of derivs module and densityforce routine
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, densityforce, deriv, dim, dust, eos, io, kernel,
!    linklist, mpiutils, nicil, options, part, physcon, testutils,
!    timestep_ind, timing, unifdis, units, viscosity
!+
!--------------------------------------------------------------------------
module testderivs
 use part, only:massoftype
 implicit none

 public :: test_derivs
 real, public :: grainsizek,graindensk
 real, parameter, private :: rhozero = 5.0

 private

contains

subroutine test_derivs(ntests,npass,string)
 use dim,          only:maxp,maxvxyzu,maxalpha,maxdvdx,ndivcurlv,nalpha,use_dust,&
                        maxdustsmall
 use boundary,     only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use eos,          only:polyk,gamma,use_entropy
 use io,           only:iprint,id,master,fatal,iverbose,nprocs
 use mpiutils,     only:reduceall_mpi
 use options,      only:tolh,alpha,alphau,alphaB,beta,ieos,psidecayfac,use_dustfrac
 use kernel,       only:radkern
 use part,         only:npart,npartoftype,igas,xyzh,hfact,vxyzu,fxyzu,fext,&
                        divcurlv,divcurlB,maxgradh,gradh,divBsymm,Bevol,dBevol,&
                        Bxyz,Bextx,Bexty,Bextz,alphaind,maxphase,rhoh,mhd,&
                        maxBevol,ndivcurlB,dvdx,dustfrac,ddustevol,temperature,&
                        idivv,icurlvx,icurlvy,icurlvz,idivB,icurlBx,icurlBy,icurlBz,deltav,dustprop,ddustprop,ndustsmall
#ifdef RADIATION
 use part,         only:radenevol,radenergy,radenflux,radthick
#endif
 use unifdis,      only:set_unifdis
 use physcon,      only:pi,au,solarm
 use deriv,        only:derivs
 use densityforce, only:get_neighbour_stats,densityiterate
 use linklist,     only:set_linklist
 use timing,       only:getused,printused
 use viscosity,    only:bulkvisc,shearparam,irealvisc
 use part,         only:iphase,isetphase,igas
 use nicil,        only:use_ambi
#ifdef IND_TIMESTEPS
 use timestep_ind, only:nactive
 use part,         only:ibin
#endif
#ifdef DUST
 use dust,         only:init_drag,idrag,K_code
 use part,         only:grainsize,graindens,ndustlarge,ndusttypes
#endif
 use units,        only:set_units
 use testutils,    only:checkval,checkvalf
 integer,          intent(inout) :: ntests,npass
 character(len=*), intent(in)    :: string
 real              :: psep,time,hzero,totmass
#ifdef IND_TIMESTEPS
 integer           :: itest,ierr,ierr2,nptest
 real              :: fracactive,speedup
 real(kind=4)      :: tallactive
 real, allocatable :: fxyzstore(:,:),dBdtstore(:,:)
#else
 integer           :: nactive
#endif
 real              :: psepblob,hblob,rhoblob,rblob,totvol,rtest
#ifdef PERIODIC
 integer           :: maxtrial,maxactual
 integer(kind=8)   :: nrhocalc,nactual,nexact
 real              :: trialmean,actualmean,realneigh
#endif
 real              :: rcut
 real              :: dtext_dum,rho1i,deint,demag,dekin,dedust,dmdust(maxdustsmall),dustfraci(maxdustsmall),tol
 real(kind=4)      :: t1,t2
 integer           :: nfailed(21),i,j,npartblob,nparttest
 integer           :: np,ieosprev,icurlvxi,icurlvyi,icurlvzi,ialphaloc,iu
 logical           :: testhydroderivs,testav,testviscderivs,testambipolar,testdustderivs
 logical           :: testmhdderivs,testdensitycontrast,testcullendehnen,testindtimesteps,testall
 real              :: vwavei,stressmax,rhoi,sonrhoi(maxdustsmall),drhodti,ddustevoli(maxdustsmall)
 integer(kind=8)   :: nptot
 real, allocatable :: dummy(:)
#ifdef IND_TIMESTEPS
 real              :: tolh_old
#endif
 logical           :: checkmask(maxp)

 if (id==master) write(*,"(a,/)") '--> TESTING DERIVS MODULE'

 testhydroderivs     = .false.
 testav              = .false.
 testviscderivs      = .false.
 testmhdderivs       = .false.
 testambipolar       = .false.
 testdustderivs      = .false.
 testdensitycontrast = .false.
 testindtimesteps    = .false.
 testcullendehnen    = .false.
 testall             = .false.
 select case(string)
 case('derivshydro','derivhydro','hydroderivs')
    testhydroderivs = .true.
 case('derivsav','derivav','avderivs')
    testav = .true.
 case('derivsvisc','derivvisc','viscderivs')
    testviscderivs = .true.
 case('derivsmhd','derivmhd','mhdderivs')
    testmhdderivs = .true.
 case('derivsdust','derivdust','dustderivs','dustderiv')
    testdustderivs = .true.
 case('derivsambi','derivambi','ambiderivs')
    testambipolar = .true.
 case('derivsdenscontrast','derivdens','derivscontrast')
    testdensitycontrast = .true.
 case('derivscd','derivcullendehnen','derivsswitch')
    testcullendehnen = .true.
 case('derivsindtimesteps','derivsinddts','derivsind')
    testindtimesteps = .true.
 case default
    testall = .true.
 end select

 iprint = 6
 iverbose = max(iverbose,2)
 psep = dxbound/100.
 npart = nint(dxbound/psep)*nint(dybound/psep)*nint(dzbound/psep)
#ifndef MPI
 if (npart > maxp) call fatal('testsuite','maxp too low to run derivs test suite')
#endif
 icurlvxi = icurlvx ! avoid compiler warnings
 icurlvyi = icurlvy
 icurlvzi = icurlvz
 iu = 4 ! avoid compiler warnings

 time = 0.
 npartoftype(:) = 0
 npart = 0
 totmass = rhozero*dxbound*dybound*dzbound

 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,hfact,npart,xyzh)
 np = npart

 if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
 nactive = npart
 npartoftype(1) = npart
 nptot = reduceall_mpi('+',npart)
 massoftype(1) = totmass/reduceall_mpi('+',npart)

#ifndef PERIODIC
 ! exclude particles near edge
 rcut = min(xmax,ymax,zmax) - 2.*radkern*hfact*psep
#else
 ! include all
 rcut = sqrt(huge(rcut))
#endif

 print*,'thread ',id,' npart = ',npart
 if (id==master) print "(a,g9.2)",' hfact = ',hfact

 hzero = hfact*(massoftype(1)/rhozero)**(1./3.)
!
!--make sure AV is off
!
 call reset_dissipation_to_zero
 tolh = 1.e-5 ! BEWARE: if you change this, some of the test results will be different

!
!--use isothermal or adiabatic equation of state for all tests
!
 if (maxvxyzu==4) then
    ieos = 2
    gamma = 5./3.
 else
    ieos = 1
    gamma = 1.0
 endif

 testhydro: if (testhydroderivs .or. testall) then
!
!--calculate pure hydro derivatives with velocity and
!  pressure distributions (no viscosity)
!
    if (id==master) write(*,"(/,a)") '--> testing Hydro derivatives '
    call set_velocity_and_energy
    call reset_mhd_to_zero
    !
    !--calculate derivatives
    !
    call getused(t1)
    call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
    call getused(t2)
    if (id==master) call printused(t1)
    call rcut_checkmask(rcut,xyzh,npart,checkmask)
    !
    !--check hydro quantities come out as they should do
    !
    nfailed(:) = 0
    call checkval(np,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)',checkmask)
    call checkvalf(np,xyzh,divcurlv(1,:),divvfunc,1.e-3,nfailed(2),'divv',checkmask)
    if (ndivcurlv >= 4) then
       call checkvalf(np,xyzh,divcurlv(icurlvxi,:),curlvfuncx,1.5e-3,nfailed(3),'curlv(x)',checkmask)
       call checkvalf(np,xyzh,divcurlv(icurlvyi,:),curlvfuncy,1.e-3,nfailed(4),'curlv(y)',checkmask)
       call checkvalf(np,xyzh,divcurlv(icurlvzi,:),curlvfuncz,1.e-3,nfailed(5),'curlv(z)',checkmask)
    endif
    if (maxgradh==maxp) then
       call checkval(np,gradh(1,:),1.01948,1.e-5,nfailed(6),'gradh',checkmask)
    endif
    if (maxvxyzu==4) then
       call checkvalf(np,xyzh,fxyzu(1,:),forcefuncx,1.e-3,nfailed(7),'force(x)',checkmask)
       call checkvalf(np,xyzh,fxyzu(2,:),forcefuncy,1.e-3,nfailed(8),'force(y)',checkmask)
       call checkvalf(np,xyzh,fxyzu(3,:),forcefuncz,1.e-3,nfailed(9),'force(z)',checkmask)
       if (use_entropy .or. ieos /= 2) then
          call checkval(np,fxyzu(iu,:),0.,epsilon(fxyzu),nfailed(10),'den/dt',checkmask)
       else
          call checkvalf(np,xyzh,fxyzu(iu,1:np)/((gamma-1.)*vxyzu(iu,1:np)),dudtfunc,1.e-3,nfailed(10),'du/dt',checkmask)
       endif
    endif
    !
    !--also check that the number of neighbours is correct
    !
#ifdef PERIODIC
    if (id==master) then
       call get_neighbour_stats(trialmean,actualmean,maxtrial,maxactual,nrhocalc,nactual)
       realneigh = 4./3.*pi*(hfact*radkern)**3
       call checkval(actualmean,real(int(realneigh)),tiny(0.),nfailed(11),'mean nneigh')
       call checkval(maxactual,int(realneigh),0,nfailed(12),'max nneigh')
       nexact = 2*nptot
       call checkval(nrhocalc,nexact,0,nfailed(13),'n density calcs')
       nexact = nptot*int(realneigh)
       call checkval(nactual,nexact,0,nfailed(14),'total nneigh')
    endif
#endif
    !
    !--check that the timestep bin has been set
    !
#ifdef IND_TIMESTEPS
    call checkval(all(ibin(1:npart) > 0),.true.,nfailed(15),'ibin > 0')
#endif

    ntests = ntests + 1
    if (all(nfailed==0)) npass = npass + 1

#ifdef IND_TIMESTEPS
    tallactive = t2-t1

    do itest=0,nint(log10(real(nptot)))-1
       nactive = 10**itest
       if (id==master) write(*,"(/,a,i10,a)") '--> testing Hydro derivatives (on ',nactive,' active particles)'
       call set_velocity_and_energy
       do i=1,npart
          if (i <= nactive/nprocs) then
             iphase(i) = isetphase(igas,iactive=.true.)
             xyzh(4,i) = hzero
          else
             iphase(i) = isetphase(igas,iactive=.false.)
          endif
       enddo
       call reset_mhd_to_zero
       !
       !--check timing for one active particle
       !
       call getused(t1)
       call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                   Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
       call getused(t2)
       if (id==master) then
          fracactive = nactive/real(npart)
          speedup = (t2-t1)/tallactive
          write(*,"(1x,'(',3(a,f9.5,'%'),')')") &
                'moved ',100.*fracactive,' of particles in ',100.*speedup, &
                ' of time, efficiency = ',100.*fracactive/speedup
       endif

       !
       ! Note that we check ALL values, including the inactives. That is we check
       ! that the inactives have preserved their values from last time they were
       ! calculated (finds bug of mistakenly setting inactives to zero)
       !
       nfailed(:) = 0
       call checkval(np,xyzh(4,1:np),hzero,3.e-4,nfailed(1),'h (density)',checkmask)
       call checkvalf(np,xyzh,divcurlv(1,1:np),divvfunc,1.e-3,nfailed(2),'divv',checkmask)
       if (ndivcurlv >= 4) then
          call checkvalf(np,xyzh,divcurlv(icurlvxi,1:np),curlvfuncx,1.5e-3,nfailed(3),'curlv(x)',checkmask)
          call checkvalf(np,xyzh,divcurlv(icurlvyi,1:np),curlvfuncy,1.e-3,nfailed(4),'curlv(y)',checkmask)
          call checkvalf(np,xyzh,divcurlv(icurlvzi,1:np),curlvfuncz,1.e-3,nfailed(5),'curlv(z)',checkmask)
       endif
       if (maxgradh==maxp) then
          call checkval(np,gradh(1,1:np),1.01948,1.e-5,nfailed(6),'gradh',checkmask)
       endif
       if (maxvxyzu==4) then
          call checkvalf(np,xyzh,fxyzu(1,1:np),forcefuncx,1.e-3,nfailed(7),'force(x)',checkmask)
          call checkvalf(np,xyzh,fxyzu(2,1:np),forcefuncy,1.e-3,nfailed(8),'force(y)',checkmask)
          call checkvalf(np,xyzh,fxyzu(3,1:np),forcefuncz,1.e-3,nfailed(9),'force(z)',checkmask)
          if (use_entropy .or. ieos /= 2) then
             call checkval(np,fxyzu(iu,1:np),0.,epsilon(fxyzu),nfailed(10),'den/dt',checkmask)
          else
             allocate(dummy(1:np))
             dummy(1:np) = fxyzu(iu,1:np)/((gamma-1.)*vxyzu(iu,1:np))
             call checkvalf(np,xyzh,dummy(1:np),dudtfunc,1.e-3,nfailed(11),'du/dt',checkmask)
             deallocate(dummy)
          endif
       endif

       ntests = ntests + 1
       if (all(nfailed==0)) npass = npass + 1
       !
       !--reset all particles to active for subsequent tests
       !
       call reset_allactive()
    enddo
#endif

 endif testhydro

 testavderivs: if (testav .or. testall) then
#ifdef IND_TIMESTEPS
    do itest=nint(log10(real(nptot))),0,-2
       nactive = 10**itest
#endif
!
!--check artificial viscosity terms (pressure + av)
!
       if (id==master) then
#ifdef DISC_VISCOSITY
          write(*,"(/,a)") '--> testing artificial viscosity terms (disc viscosity)'
#else
          if (maxalpha==maxp) then
             write(*,"(/,a)") '--> testing artificial viscosity terms (individual alpha)'
          else
             write(*,"(/,a)") '--> testing artificial viscosity terms (constant alpha)'
          endif
#endif
#ifdef IND_TIMESTEPS
          if (nactive /= npart) write(*,"(a,i10,a)") '    (on ',nactive,' active particles)'
#endif
       endif
       if (maxvxyzu < 4) polyk = 3.
       do i=1,npart
          vxyzu(1,i) = vx(xyzh(:,i))
          vxyzu(2,i) = vy(xyzh(:,i))
          vxyzu(3,i) = vz(xyzh(:,i))
          if (maxvxyzu==4) vxyzu(iu,i) = uthermconst(xyzh(:,i))
       enddo
       call set_active(npart,nactive,igas)
       call reset_mhd_to_zero
       call reset_dissipation_to_zero    ! turn off any other dissipation
       alpha  = 0.753 ! an arbitrary number that is not 1 or 0.
       if (maxalpha==maxp) alphaind(1,:) = real(alpha,kind=kind(alphaind))

       call getused(t1)
       call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
       if (id==master) call printused(t1)
       call rcut_checkmask(rcut,xyzh,npart,checkmask)
       nfailed(:) = 0
       call checkval(np,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)',checkmask)
       call checkvalf(np,xyzh,divcurlv(1,:),divvfunc,1.e-3,nfailed(2),'divv',checkmask)
       if (ndivcurlv >= 4) then
          call checkvalf(np,xyzh,divcurlv(icurlvxi,:),curlvfuncx,1.5e-3,nfailed(3),'curlv(x)',checkmask)
          call checkvalf(np,xyzh,divcurlv(icurlvyi,:),curlvfuncy,1.e-3,nfailed(4),'curlv(y)',checkmask)
          call checkvalf(np,xyzh,divcurlv(icurlvzi,:),curlvfuncz,1.e-3,nfailed(5),'curlv(z)',checkmask)
       endif
       call checkvalf(np,xyzh,fxyzu(1,:),forceavx,3.2e-2,nfailed(3),'art. visc force(x)',checkmask)
       call checkvalf(np,xyzh,fxyzu(2,:),forceavy,2.4e-2,nfailed(4),'art. visc force(y)',checkmask)
       call checkvalf(np,xyzh,fxyzu(3,:),forceavz,2.4e-2,nfailed(5),'art. visc force(z)',checkmask)

       ntests = ntests + 1
       if (all(nfailed==0)) npass = npass + 1
#ifdef IND_TIMESTEPS
       call reset_allactive()
    enddo
#endif

 endif testavderivs
!
!--check computation of d(divv)/dt) in Cullen & Dehnen switch
!
 testcdswitch: if (testcullendehnen .or. testall) then
    if (maxalpha==maxp .and. nalpha > 1) then
       if (id==master) write(*,"(/,a)") '--> testing ddivv/dt in Cullen & Dehnen switch'

       do i=1,npart
          vxyzu(1,i) = vx(xyzh(:,i))
          vxyzu(2,i) = vy(xyzh(:,i))
          vxyzu(3,i) = vz(xyzh(:,i))
          if (maxvxyzu==4) vxyzu(iu,i) = uthermconst(xyzh(:,i))
          ! set acceleration also
          fxyzu(1,i) = vx(xyzh(:,i))
          fxyzu(2,i) = vy(xyzh(:,i))
          fxyzu(3,i) = vz(xyzh(:,i))
          fext(:,i)  = 0.
       enddo
       nactive = np
       call set_active(npart,nactive,igas)
       call reset_mhd_to_zero
       call reset_dissipation_to_zero    ! turn off any other dissipation

       call getused(t1)
       ! ONLY call density, since we do not want accelerations being reset
       call set_linklist(npart,nactive,xyzh,vxyzu)
       call densityiterate(1,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,&
                           Bevol,stressmax,fxyzu,fext,alphaind,gradh&
#ifdef RADIATION
                           ,radenevol,radenflux,radenergy,radthick&
#endif
                           )
       if (id==master) call printused(t1)

       nfailed(:) = 0
       call checkval(np,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)')
       call checkvalf(np,xyzh,divcurlv(1,:),divvfunc,1.e-3,nfailed(2),'divv')
       if (ndivcurlv >= 4) then
          call checkvalf(np,xyzh,divcurlv(icurlvxi,:),curlvfuncx,1.5e-3,nfailed(3),'curlv(x)')
          call checkvalf(np,xyzh,divcurlv(icurlvyi,:),curlvfuncy,1.e-3,nfailed(4),'curlv(y)')
          call checkvalf(np,xyzh,divcurlv(icurlvzi,:),curlvfuncz,1.e-3,nfailed(5),'curlv(z)')
       endif
       if (nalpha >= 2) then
          ialphaloc = 2
          call checkvalf(np,xyzh,alphaind(ialphaloc,:),alphalocfunc,3.5e-4,nfailed(6),'alphaloc')
       endif

       ntests = ntests + 1
       if (all(nfailed==0)) npass = npass + 1
    else
       if (id==master) write(*,"(/,a)") '--> SKIPPING Cullen-Dehnen terms (need nalpha=2)'
    endif
 endif testcdswitch

 testvisc: if (testviscderivs .or. testall) then
!
!--check viscosity terms (no pressure)
!
    if (id==master) then
       if (maxdvdx==maxp) then
          write(*,"(/,a)") '--> testing physical viscosity terms (two first derivatives)'
       else
          write(*,"(/,a)") '--> testing physical viscosity terms (direct second derivatives)'
       endif
    endif
    polyk = 0.
    do i=1,npart
       vxyzu(1,i) = vx(xyzh(:,i))
       vxyzu(2,i) = vy(xyzh(:,i))
       vxyzu(3,i) = vz(xyzh(:,i))
       if (maxvxyzu==4) vxyzu(iu,i) = 0.
    enddo
    call reset_mhd_to_zero
    call reset_dissipation_to_zero
    irealvisc = 1
    shearparam = 6.66
    bulkvisc = 0.75

    call getused(t1)
    call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
    if (id==master) call printused(t1)
    call rcut_checkmask(rcut,xyzh,npart,checkmask)

    nfailed(:) = 0
    call checkval(np,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)',checkmask)
    call checkvalf(np,xyzh,divcurlv(1,:),divvfunc,1.e-3,nfailed(2),'divv',checkmask)
    if (ndivcurlv >= 4) then
       call checkvalf(np,xyzh,divcurlv(icurlvxi,:),curlvfuncx,1.5e-3,nfailed(3),'curlv(x)',checkmask)
       call checkvalf(np,xyzh,divcurlv(icurlvyi,:),curlvfuncy,1.e-3,nfailed(4),'curlv(y)',checkmask)
       call checkvalf(np,xyzh,divcurlv(icurlvzi,:),curlvfuncz,1.e-3,nfailed(5),'curlv(z)',checkmask)
    endif
    if (maxdvdx==maxp) then
       call checkvalf(np,xyzh,dvdx(1,:),dvxdx,1.7e-3,nfailed(6),  'dvxdx',checkmask)
       call checkvalf(np,xyzh,dvdx(2,:),dvxdy,2.5e-15,nfailed(7), 'dvxdy',checkmask)
       call checkvalf(np,xyzh,dvdx(3,:),dvxdz,2.5e-15,nfailed(8), 'dvxdz',checkmask)
       call checkvalf(np,xyzh,dvdx(4,:),dvydx,1.e-3,nfailed(9),   'dvydx',checkmask)
       call checkvalf(np,xyzh,dvdx(5,:),dvydy,2.5e-15,nfailed(10), 'dvydy',checkmask)
       call checkvalf(np,xyzh,dvdx(6,:),dvydz,1.e-3,nfailed(11),  'dvydz',checkmask)
       call checkvalf(np,xyzh,dvdx(7,:),dvzdx,2.5e-15,nfailed(9), 'dvzdx',checkmask)
       call checkvalf(np,xyzh,dvdx(8,:),dvzdy,1.5e-3,nfailed(10), 'dvzdy',checkmask)
       call checkvalf(np,xyzh,dvdx(9,:),dvzdz,2.5e-15,nfailed(11),'dvzdz',checkmask)
    endif

    call checkvalf(np,xyzh,fxyzu(1,:),forceviscx,4.e-2,nfailed(12),'viscous force(x)',checkmask)
    call checkvalf(np,xyzh,fxyzu(2,:),forceviscy,3.e-2,nfailed(13),'viscous force(y)',checkmask)
    call checkvalf(np,xyzh,fxyzu(3,:),forceviscz,3.1e-2,nfailed(14),'viscous force(z)',checkmask)
    !
    !--also check that the number of neighbours is correct
    !
#ifdef PERIODIC
    if (id==master) then
       call get_neighbour_stats(trialmean,actualmean,maxtrial,maxactual,nrhocalc,nactual)
       realneigh = 4./3.*pi*(hfact*radkern)**3
       call checkval(actualmean,real(int(realneigh)),tiny(0.),nfailed(15),'mean nneigh')
       call checkval(maxactual,int(realneigh),0,nfailed(16),'max nneigh')
       if (testall) then
          nexact = nptot  ! should be no iterations here
          call checkval(nrhocalc,nexact,0,nfailed(17),'n density calcs')
       endif
       nexact = nptot*int(realneigh)
       call checkval(nactual,nexact,0,nfailed(18),'total nneigh')
    endif
#endif
    !
    !--check that \sum m (du/dt + v.dv/dt) = 0.
    !  only applies if all particles active - with individual timesteps
    !
    if (maxvxyzu==4 .and. nactive==npart) then
       deint = 0.
       dekin = 0.
       do i=1,npart
          deint = deint + fxyzu(iu,i)
          dekin = dekin + dot_product(vxyzu(1:3,i),fxyzu(1:3,i))
       enddo
       deint = reduceall_mpi('+',deint)
       dekin = reduceall_mpi('+',dekin)
       nfailed(:) = 0
       if (maxdvdx==maxp) then
          tol = 1.52e-6
       else
          tol = 5.e-12
       endif
       call checkval(massoftype(1)*(deint + dekin),0.,tol,nfailed(19),'\sum v.dv/dt + du/dt = 0')

       ! also check that dissipation is positive definite
       call checkval(all(fxyzu(iu,1:np) >= 0.),.true.,nfailed(20),'du/dt >= 0 for all particles')
    endif

    ntests = ntests + 1
    if (all(nfailed==0)) npass = npass + 1
 endif testvisc

 testdust: if (testdustderivs .or. testall) then
!
!--check derivative terms for one-fluid dust
!
    if (use_dust) use_dustfrac=.true.
    if (use_dustfrac) then
       if (id==master) write(*,"(/,a)") '--> testing dust evolution terms'
#ifdef DUST
       idrag   = 2
       gamma   = 5./3.
       !--Warning, K_code is not well defined when using multiple dust grains
       !  and ONLY makes sense IFF all dust grains are identical (although
       !  potentially binned with unequal densities).
       !  K_code and K_k are related via: K_k = eps_k/eps*K_code)
       K_code = 10.
       grainsize = 0.01
       graindens = 3.
       ndustsmall = maxdustsmall
       ndustlarge = 0
       ndusttypes = ndustsmall + ndustlarge
       !need to set units if testing with physical drag
       !call set_units(dist=au,mass=solarm,G=1.d0)
       call init_drag(nfailed(1))
#endif
       polyk = 0.
       call reset_mhd_to_zero
       call reset_dissipation_to_zero
       call set_velocity_and_energy
       do i=1,npart
          do j=1,ndustsmall
             dustfrac(j,i) = real(dustfrac_func(xyzh(:,i)),kind=kind(dustfrac))
          enddo
       enddo

       call getused(t1)
       call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                   Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
       if (id==master) call printused(t1)

       nfailed(:) = 0
       call checkval(np,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)')
       call checkvalf(np,xyzh,divcurlv(1,:),divvfunc,1.e-3,nfailed(2),'divv')
       do j=1,1 !ndustsmall !--Only need one because all dust species are identical
#ifdef DUST
          grainsizek = grainsize(j)
          graindensk = graindens(j)
#endif
          call checkvalf(np,xyzh,ddustevol(j,:),ddustevol_func,4.e-5,nfailed(3),'deps/dt')
          if (maxvxyzu>=4) call checkvalf(np,xyzh,fxyzu(iu,:),dudtdust_func,5.e-4,nfailed(4),'du/dt')
          call checkvalf(np,xyzh,deltav(1,j,:),deltavx_func,2.3e-5,nfailed(5),'deltavx')
       enddo

       ntests = ntests + 1
       if (all(nfailed==0)) npass = npass + 1

       !
       !--check energy conservation, i.e. \sum m [v.dv/dt + (1 - epsilon)*du/dt - u deps/dt] = 0.
       !  this is equation (41) in Price & Laibe (2015)
       !
       if (maxvxyzu==4 .and. nactive==npart) then
          dekin  = 0.
          deint  = 0.
          dedust = 0.
          dmdust(:) = 0.
          do i=1,npart
             dustfraci(:)  = dustfrac(1:maxdustsmall,i)
             rhoi          = rhoh(xyzh(4,i),massoftype(igas))
             drhodti       = -rhoi*divcurlv(1,i)
!------------------------------------------------
!--sqrt(rho*epsilon) method
!             sonrhoi(:)    = sqrt(dustfrac(1:maxdustsmall,i)/rhoi)
!             ddustevoli(:) = 2.*sonrhoi(:)*ddustevol(:,i) - sonrhoi(:)**2*drhodti
!------------------------------------------------
!--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
             sonrhoi(:)    = sqrt(dustfraci(:)*(1.-dustfraci(:)))
             ddustevoli(:) = 2.*sonrhoi(:)*(1.-dustfraci(:))*ddustevol(:,i)
!------------------------------------------------
!--asin(sqrt(epsilon)) method
!             sonrhoi(:)    = asin(sqrt(dustfrac(1:maxdustsmall,i)))
!             ddustevoli(:) = 2.*cos(sonrhoi(:))*sin(sonrhoi(:))*ddustevol(:,i)
!------------------------------------------------
             dmdust(:)     = dmdust(:) + ddustevoli(:)
             dekin  = dekin  + dot_product(vxyzu(1:3,i),fxyzu(1:3,i))
             deint  = deint  + (1. - sum(dustfraci))*fxyzu(iu,i)
             dedust = dedust - vxyzu(iu,i)*sum(ddustevoli)
          enddo
          dmdust  = reduceall_mpi('+',dmdust)
          dekin   = reduceall_mpi('+',dekin)
          deint   = reduceall_mpi('+',deint)
          dedust  = reduceall_mpi('+',dedust)

          nfailed(:) = 0
          !print "(3(a,es17.10))",' dE_kin = ',dekin,' dE_therm = ',deint,' dE_dust = ',dedust
          call checkval(massoftype(1)*(dekin + deint + dedust),0.,6.5e-15,nfailed(1),'energy conservation (dE=0)')
          do i=1,ndustsmall
             call checkval(massoftype(1)*(dmdust(i)),0.,1.e-15,nfailed(2),'dust mass conservation')
          enddo
          ntests = ntests + 1
          if (nfailed(1)==0) npass = npass + 1
       endif
       !
       ! reset dustfrac to zero for subsequent tests
       !
       dustfrac(:,:) = 0.

    else
       if (id==master) write(*,"(/,a)") '--> SKIPPING dust evolution terms (need -DDUST)'
    endif
 endif testdust

!
!--calculate derivatives with MHD forces ON, zero pressure
!
 testmhd: if (testmhdderivs .or. testall) then
    ! obtain smoothing lengths
    call set_linklist(npart,nactive,xyzh,vxyzu)
    call densityiterate(2,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,&
                      Bevol,stressmax,fxyzu,fext,alphaind,gradh&
#ifdef RADIATION
                      ,radenevol,radenflux,radenergy,radthick&
#endif
                      )
#ifdef IND_TIMESTEPS
    do itest=nint(log10(real(nptot))),0,-2
       nactive = 10**itest
#endif
       polyk = 0.
       call reset_mhd_to_zero
       call reset_dissipation_to_zero
       if (mhd) then
          if (id==master) then
             write(*,"(/,a)") '--> testing MHD derivatives (using B/rho directly)'
             if (nactive /= np) write(*,"(a,i10,a)") '    (on ',nactive,' active particles)'
          endif
          Bextx = 2.0e-1
          Bexty = 3.0e-1
          Bextz = 0.5
          do i=1,npart
             vxyzu(1,i) = vx(xyzh(:,i))
             vxyzu(2,i) = vy(xyzh(:,i))
             vxyzu(3,i) = vz(xyzh(:,i))
             rho1i = 1.0/rhoh(xyzh(4,i),massoftype(igas))
             Bxyz(1,i) = Bx(xyzh(:,i))
             Bxyz(2,i) = By(xyzh(:,i))
             Bxyz(3,i) = Bz(xyzh(:,i))
             Bevol(1,i) = Bxyz(1,i) * rho1i
             Bevol(2,i) = Bxyz(2,i) * rho1i
             Bevol(3,i) = Bxyz(3,i) * rho1i
             if (maxvxyzu==4) vxyzu(iu,i) = 0.
          enddo
          call set_active(npart,nactive/nprocs,igas)
          call getused(t1)
          call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                   Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
          if (id==master) call printused(t1)
          !
          !--check that various quantities come out as they should do
          !
          nfailed(:) = 0
          call checkval(np,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)')

          call checkvalf(np,xyzh,divBsymm(:),divBfunc,2.e-3,nfailed(2),'divB (symm)')
          call checkvalf(np,xyzh,dBevol(1,:),dBxdt,2.e-3,nfailed(3),'dBx/dt')
          call checkvalf(np,xyzh,dBevol(2,:),dBydt,2.e-3,nfailed(4),'dBy/dt')
          call checkvalf(np,xyzh,dBevol(3,:),dBzdt,2.e-2,nfailed(5),'dBz/dt')

          call checkvalf(np,xyzh,fxyzu(1,:),forcemhdx,2.5e-2,nfailed(9),'mhd force(x)')
          call checkvalf(np,xyzh,fxyzu(2,:),forcemhdy,2.5e-2,nfailed(10),'mhd force(y)')
          call checkvalf(np,xyzh,fxyzu(3,:),forcemhdz,2.5e-2,nfailed(11),'mhd force(z)')
          if (ndivcurlB >= 1) then
             call checkvalf(np,xyzh,divcurlB(idivB,:),divBfunc,1.e-3,nfailed(12),'div B (diff)')
          endif
          if (ndivcurlB >= 4) then
             call checkvalf(np,xyzh,divcurlB(icurlBx,:),curlBfuncx,1.e-3,nfailed(13),'curlB(x)')
             call checkvalf(np,xyzh,divcurlB(icurlBy,:),curlBfuncy,1.e-3,nfailed(14),'curlB(y)')
             call checkvalf(np,xyzh,divcurlB(icurlBz,:),curlBfuncz,1.e-3,nfailed(15),'curlB(z)')
          endif
          ntests = ntests + 1
          if (all(nfailed==0)) npass = npass + 1

       endif
#ifdef IND_TIMESTEPS
       call reset_allactive()
    enddo
    do itest=nint(log10(real(nptot))),0,-2
       nactive = 10**itest
#endif
       if (mhd) then
          if (id==master) then
             write(*,"(/,a)") '--> testing artificial resistivity terms'
             if (nactive /= np) write(*,"(a,i10,a)") '    (on ',nactive,' active particles)'
          endif
          call reset_mhd_to_zero
          call reset_dissipation_to_zero
          alphaB = 0.214
          polyk = 0.
          ieosprev = ieos
          ieos  = 1  ! isothermal eos, so that the PdV term is zero
          do i=1,npart
             vxyzu(:,i) = 0.
             rho1i   = 1.0/rhoh(xyzh(4,i),massoftype(igas))
             Bevol(:,i) = 0.
             Bxyz(1,i) = Bx(xyzh(:,i))
             Bxyz(2,i) = By(xyzh(:,i))
             Bxyz(3,i) = Bz(xyzh(:,i))
             Bevol(1,i) = Bxyz(1,i) * rho1i
             Bevol(2,i) = Bxyz(2,i) * rho1i
             Bevol(3,i) = Bxyz(3,i) * rho1i
          enddo
          call set_active(npart,nactive,igas)
          call getused(t1)
          call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                   Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
          if (id==master) call printused(t1)
          call rcut_checkmask(rcut,xyzh,npart,checkmask)
          !
          !--check that various quantities come out as they should do
          !
          nfailed(:) = 0
          !
          !--resistivity test is very approximate
          !  To do a proper test, multiply by h/rij in densityforce
          !
          call checkvalf(np,xyzh,dBevol(1,:),dBxdtresist,3.7e-2,nfailed(1),'dBx/dt (resist)',checkmask)
          call checkvalf(np,xyzh,dBevol(2,:),dBydtresist,3.4e-2,nfailed(2),'dBy/dt (resist)',checkmask)
          call checkvalf(np,xyzh,dBevol(3,:),dBzdtresist,2.2e-1,nfailed(3),'dBz/dt (resist)',checkmask)

          ntests = ntests + 1
          if (all(nfailed==0)) npass = npass + 1
          !
          !--check that \sum m (du/dt + B/rho.dB/dt) = 0.
          !  only applies if all particles active - with individual timesteps
          !  we just hope that du/dt has not changed all that much on non-active particles
          !
          if (maxvxyzu==4 .and. nactive==npart) then
             deint = 0.
             demag = 0.
             do i=1,npart
                rho1i = 1./rhoh(xyzh(4,i),massoftype(1))
                deint = deint + fxyzu(iu,i)
                demag = demag + dot_product(Bevol(1:3,i),dBevol(1:3,i))*rho1i
             enddo
             nfailed(:) = 0
             call checkval(deint + demag,0.,2.7e-3,nfailed(1),'\sum du/dt + B.dB/dt = 0')
             ntests = ntests + 1
             if (nfailed(1)==0) npass = npass + 1
          endif

          !--restore ieos
          ieos = ieosprev

       endif
#ifdef IND_TIMESTEPS
       call reset_allactive()
    enddo
    tolh_old = tolh
    tolh = 1.e-7
    do itest=nint(log10(real(nptot))),0,-2
       nactive = 10**itest
#endif
       if (mhd .and. maxBevol==4) then
          if (id==master) then
             write(*,"(/,a)") '--> testing div B cleaning terms'
             if (nactive /= np) write(*,"(a,i10,a)") '    (on ',nactive,' active particles)'
          endif
          call reset_mhd_to_zero
          call reset_dissipation_to_zero
          psidecayfac = 0.8
          polyk = 2.
          ieosprev = ieos
          ieos  = 1  ! isothermal eos
          do i=1,npart
             vxyzu(1,i) = vx(xyzh(:,i))
             vxyzu(2,i) = vy(xyzh(:,i))
             vxyzu(3,i) = vz(xyzh(:,i))
             rho1i      = 1.0/rhoh(xyzh(4,i),massoftype(igas))
             Bxyz(1,i) = Bx(xyzh(:,i))
             Bxyz(2,i) = By(xyzh(:,i))
             Bxyz(3,i) = Bz(xyzh(:,i))
             Bevol(1,i) = Bxyz(1,i) * rho1i
             Bevol(2,i) = Bxyz(2,i) * rho1i
             Bevol(3,i) = Bxyz(3,i) * rho1i

             vwavei = sqrt(polyk + (Bxyz(1,i) * Bxyz(1,i) + Bxyz(2,i) * Bxyz(2,i) + Bxyz(3,i) * Bxyz(3,i)) &
                                   / rhoh(xyzh(4,i),massoftype(1)))
             Bevol(4,i) = psi(xyzh(:,i))/vwavei
             if (maxvxyzu==4) vxyzu(iu,i) = 0.
          enddo
          call set_active(npart,nactive,igas)
          call getused(t1)
          call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                   Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
          if (id==master) call printused(t1)
          !
          !--check that various quantities come out as they should do
          !
          nfailed(:) = 0
          call checkval(np,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)')
          call checkvalf(np,xyzh,divBsymm(:),divBfunc,1.e-3,nfailed(2),'divB')
          call checkvalf(np,xyzh,dBevol(1,:),dpsidx,8.5e-4,nfailed(3),'gradpsi_x')
          call checkvalf(np,xyzh,dBevol(2,:),dpsidy,8.5e-4,nfailed(4),'gradpsi_y')
          call checkvalf(np,xyzh,dBevol(3,:),dpsidz,2.e-3,nfailed(5),'gradpsi_z')
          !--can't do dpsi/dt check because we use vsigdtc = max over neighbours
          !call checkvalf(np,xyzh,dBevol(4,:),dpsidt,6.e-3,nfailed(6),'dpsi/dt')

          ntests = ntests + 1
          if (all(nfailed==0)) npass = npass + 1

          !--restore ieos
          ieos = ieosprev
       endif
#ifdef IND_TIMESTEPS
       call reset_allactive()
    enddo
    tolh = tolh_old
    do itest=nint(log10(real(nptot))),0,-2
       nactive = 10**itest
#endif
       if (mhd .and. use_ambi .and. testambipolar) then
          if (id==master) then
             write(*,"(/,a)") '--> testing Ambipolar diffusion terms'
             if (nactive /= np) write(*,"(a,i10,a)") '    (on ',nactive,' active particles)'
          endif
          call reset_mhd_to_zero
          call reset_dissipation_to_zero
          psidecayfac = 0.0
          polyk = 0.
          ieosprev = ieos
          ieos  = 1  ! isothermal eos
          do i=1,npart
             vxyzu(1,i) = vx(xyzh(:,i))
             vxyzu(2,i) = vy(xyzh(:,i))
             vxyzu(3,i) = vz(xyzh(:,i))
             rho1i      = 1.0/rhoh(xyzh(4,i),massoftype(igas))
             Bxyz(1,i) = Bx(xyzh(:,i))
             Bxyz(2,i) = By(xyzh(:,i))
             Bxyz(3,i) = Bz(xyzh(:,i))
             Bevol(1,i) = Bxyz(1,i) * rho1i
             Bevol(2,i) = Bxyz(2,i) * rho1i
             Bevol(3,i) = Bxyz(3,i) * rho1i
             if (maxBevol>=4) Bevol(4,i) = 0.
             if (maxvxyzu==4) vxyzu(iu,i) = 0.
          enddo
          call set_active(npart,nactive,igas)
          call getused(t1)
          call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                   Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
          if (id==master) call printused(t1)
          !
          !--check that various quantities come out as they should do
          !
          nfailed(:) = 0
          call checkval(np,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)')
          call checkvalf(np,xyzh,dBevol(1,:),dBambix,8.5e-4,nfailed(2),'dBambi_x')
          call checkvalf(np,xyzh,dBevol(2,:),dBambiy,8.5e-4,nfailed(3),'dBambi_y')
          call checkvalf(np,xyzh,dBevol(3,:),dBambiz,2.e-3,nfailed(4),'dBambi_z')

          ntests = ntests + 1
          if (all(nfailed==0)) npass = npass + 1

          !--restore ieos
          ieos = ieosprev
       endif

#ifdef IND_TIMESTEPS
       call reset_allactive()
    enddo
#endif

 endif testmhd

!
!--calculate hydro terms in setup with density contrast
!
!
!-- TODO: this test won't pass with MPI, because the particles are shuffled around,
!         and the 'test' particles cannot be identified using the current method
!
 testdenscontrast: if ((testdensitycontrast .or. testall) .and. (nprocs == 1)) then
    if (id==master) write(*,"(/,a)") '--> testing Hydro derivs in setup with density contrast '

    npart = 0
    psep = dxbound/50.
    rblob = 0.1
    rhoblob = 1000.*rhozero
    psepblob = psep*(rhozero/rhoblob)**(1./3.)
    rtest = rblob - 2.*hfact*psep
    !
    !--setup high density blob
    !
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psepblob,hfact,npart,xyzh,rmax=rtest)
    nparttest = npart

    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psepblob,hfact,npart,xyzh,rmin=rtest,rmax=rblob)
    npartblob = npart
    !
    !--setup surrounding medium
    !
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,hfact,npart,xyzh,rmin=rblob)
    npartoftype(1) = npart
    nptot = reduceall_mpi('+',npart)
    print*,' thread ',id,' npart = ',npart,' in blob = ',npartblob,' to test = ',nparttest

    totvol = dxbound*dybound*dzbound - 4./3.*pi*rblob**3
    totmass = rhozero*totvol

    massoftype(1) = totmass/reduceall_mpi('+',npart-npartblob)
    hzero = hfact*(massoftype(1)/rhozero)**(1./3.)
    hblob = hfact*(massoftype(1)/rhoblob)**(1./3.)
    call reset_dissipation_to_zero
    call set_velocity_and_energy
    call reset_mhd_to_zero
    !
    !--calculate derivatives
    !
    nactive = npart
    call set_active(npart,nactive,igas)
    call getused(t1)
    call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
    call getused(t2)
    if (id==master) call printused(t1)
    !
    !--check hydro quantities come out as they should do
    !
    nfailed(:) = 0
    call checkval(nparttest,xyzh(4,:),hblob,4.e-4,nfailed(1),'h (density)')
    call checkvalf(nparttest,xyzh,divcurlv(1,:),divvfunc,1.e-3,nfailed(2),'divv')
    if (ndivcurlv >= 4) then
       call checkvalf(nparttest,xyzh,divcurlv(icurlvxi,:),curlvfuncx,1.5e-3,nfailed(3),'curlv(x)')
       call checkvalf(nparttest,xyzh,divcurlv(icurlvyi,:),curlvfuncy,1.e-3,nfailed(4),'curlv(y)')
       call checkvalf(nparttest,xyzh,divcurlv(icurlvzi,:),curlvfuncz,1.e-3,nfailed(5),'curlv(z)')
    endif
    if (maxvxyzu==4) then
       call checkvalf(nparttest,xyzh,fxyzu(1,:),forcefuncx,1.e-3,nfailed(6),'force(x)')
       call checkvalf(nparttest,xyzh,fxyzu(2,:),forcefuncy,1.e-3,nfailed(7),'force(y)')
       call checkvalf(nparttest,xyzh,fxyzu(3,:),forcefuncz,1.e-3,nfailed(8),'force(z)')
       if (use_entropy .or. ieos /= 2) then
          call checkval(nparttest,fxyzu(iu,:),0.,epsilon(fxyzu),nfailed(9),'den/dt')
       else
          allocate(dummy(nparttest))
          dummy(1:nparttest) = fxyzu(iu,1:nparttest)/((gamma-1.)*vxyzu(iu,1:nparttest))
          call checkvalf(nparttest,xyzh,dummy(1:nparttest),dudtfunc,1.e-3,nfailed(9),'du/dt')
          deallocate(dummy)
       endif
    endif
    !
    !--also check that the number of neighbours is correct
    !
#ifdef PERIODIC
    if (id==master) then
       call get_neighbour_stats(trialmean,actualmean,maxtrial,maxactual,nrhocalc,nactual)
       realneigh = 57.466651861721814
       call checkval(actualmean,realneigh,1.e-17,nfailed(10),'mean nneigh')
       call checkval(maxactual,988,0,nfailed(11),'max nneigh')
       !
       !-- this test does not always give the same results: depends on how the tree is built
       !
       !  nexact = 1382952  ! got this from a reference calculation
       !  call checkval(nrhocalc,nexact,0,nfailed(12),'n density calcs')
       nexact = 37263216
       call checkval(nactual,nexact,0,nfailed(13),'total nneigh')
    endif
#endif

    ntests = ntests + 1
    if (all(nfailed==0)) npass = npass + 1

#ifdef IND_TIMESTEPS
    tallactive = t2-t1
    do itest=1,nint(log10(real(nparttest)))
       nactive = 10**itest
       if (nactive > nparttest) nactive = nparttest
       if (id==master) write(*,"(/,a,i6,a)") '--> testing Hydro derivs in setup with density contrast (nactive=',nactive,') '

       call set_active(npart,nactive,igas)
       call getused(t1)
       call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                    Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
       call getused(t2)
       if (id==master) then
          fracactive = nactive/real(npart)
          speedup = (t2-t1)/tallactive
          write(*,"(1x,'(',3(a,f9.5,'%'),')')") &
                 'moved ',100.*fracactive,' of particles in ',100.*speedup, &
                 ' of time, efficiency = ',100.*fracactive/speedup
       endif
       !
       !--check hydro quantities come out as they should do
       !
       nfailed(:) = 0
       call checkval(nparttest,xyzh(4,:),hblob,4.e-4,nfailed(1),'h (density)')
       call checkvalf(nparttest,xyzh,divcurlv(idivv,:),divvfunc,1.e-3,nfailed(2),'divv')
       if (ndivcurlv >= 4) then
          call checkvalf(nparttest,xyzh,divcurlv(icurlvxi,:),curlvfuncx,1.5e-3,nfailed(3),'curlv(x)')
          call checkvalf(nparttest,xyzh,divcurlv(icurlvyi,:),curlvfuncy,1.e-3,nfailed(4),'curlv(y)')
          call checkvalf(nparttest,xyzh,divcurlv(icurlvzi,:),curlvfuncz,1.e-3,nfailed(5),'curlv(z)')
       endif
       if (maxvxyzu==4) then
          call checkvalf(nparttest,xyzh,fxyzu(1,:),forcefuncx,1.e-3,nfailed(6),'force(x)')
          call checkvalf(nparttest,xyzh,fxyzu(2,:),forcefuncy,1.e-3,nfailed(7),'force(y)')
          call checkvalf(nparttest,xyzh,fxyzu(3,:),forcefuncz,1.e-3,nfailed(8),'force(z)')
          if (use_entropy .or. ieos /= 2) then
             call checkval(nparttest,fxyzu(iu,1:nparttest),0.,epsilon(fxyzu),nfailed(9),'den/dt')
          else
             allocate(dummy(nparttest))
             dummy(1:nparttest) = fxyzu(iu,1:nparttest)/((gamma-1.)*vxyzu(iu,1:nparttest))
             call checkvalf(nparttest,xyzh,dummy(1:nparttest),dudtfunc,1.e-3,nfailed(9),'du/dt')
             deallocate(dummy)
          endif
       endif
       ntests = ntests + 1
       if (all(nfailed==0)) npass = npass + 1
    enddo
#endif

 endif testdenscontrast
!
!--test force evaluation for individual timesteps when particles have very different smoothing lengths/ranges
!
 testinddts: if (testindtimesteps .or. testall) then
#ifdef IND_TIMESTEPS
    if (id==master) write(*,"(/,a,i6,a)") '--> testing force evaluation with ind_timesteps'
    polyk = 0.
    tolh  = 1.e-9
    call reset_mhd_to_zero
    call reset_dissipation_to_zero
    alpha  = 0.753 ! an arbitrary number that is not 1 or 0.
    if (maxalpha==maxp) alphaind(1,:) = real(alpha,kind=kind(alphaind))

    npart = 0
    call set_unifdis('random',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                      psep,hfact,npart,xyzh)

    !
    !--need to initialise dBevol to zero, otherwise if cleaning is not updated
    !  then test may give NaNs
    !
    if (mhd) dBevol(4,:) = 0.0

#ifdef MPI
    !
    !--call derivs once to ensure that particles are properly balanced
    !  before allocating arrays
    !
    nactive = npart
    call set_active(npart,nactive,igas)
    call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                 Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
#endif
    !
    !--first do the calculation with all particles active, then
    !  perform the force calculation with only a fraction of particles active
    !
    ierr2 = 0
    nptest = npart/10
    allocate(fxyzstore(maxvxyzu,nptest),stat=ierr)
    if (mhd) allocate(dBdtstore(maxBevol+1,nptest),stat=ierr2)
    if (ierr /= 0 .or. ierr2 /= 0) then
       write(*,*) 'ERROR allocating memory for force test, skipping...'
    else
       do itest=1,2
          if (itest==2) then
             nactive = nptest
             if (id==master) write(*,"(/,1x,a,i6,a)") 'evaluating derivs with ',nactive,' particles active...'
          else
             nactive = npart
             if (id==master) write(*,"(1x,a)") 'evaluating derivs with all particles active...'
          endif
          call set_velocity_and_energy
          do i=1,npart
             if (mhd) then
                Bevol(1,i) = Bx(xyzh(:,i))
                Bevol(2,i) = By(xyzh(:,i))
                Bevol(3,i) = Bz(xyzh(:,i))
                if (maxBevol >= 4) Bevol(4,i) = psi(xyzh(:,i))
             endif
          enddo
          call set_active(npart,nactive,igas)
          call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                       Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
          if (itest==1) then
             fxyzstore(:,1:nptest) = fxyzu(:,1:nptest)
             if (mhd) then
                dBdtstore(1:maxBevol,1:nptest) = dBevol(1:maxBevol,1:nptest)
                dBdtstore(maxBevol+1,1:nptest) = divBsymm(1:nptest)
             endif
          else
             nfailed(:) = 0
             call checkval(nptest,fxyzu(1,:),fxyzstore(1,1:nptest),1.e-4,nfailed(1),'force(x)')
             call checkval(nptest,fxyzu(2,:),fxyzstore(2,1:nptest),1.e-4,nfailed(2),'force(y)')
             call checkval(nptest,fxyzu(3,:),fxyzstore(3,1:nptest),1.e-4,nfailed(3),'force(z)')
             if (maxvxyzu >= 4) then
                call checkval(nptest,fxyzu(iu,:),fxyzstore(4,1:nptest),1.e-5,nfailed(4),'du/dt')
             endif
             if (mhd) then
                call checkval(nptest,dBevol(1,:),dBdtstore(1,1:nptest),1.e-5,nfailed(5),'dBx/dt')
                call checkval(nptest,dBevol(2,:),dBdtstore(2,1:nptest),1.e-5,nfailed(6),'dBy/dt')
                call checkval(nptest,dBevol(3,:),dBdtstore(3,1:nptest),1.e-5,nfailed(7),'dBz/dt')
                if (maxBevol >= 4) then
                   call checkval(nptest,dBevol(4,:),dBdtstore(4,1:nptest),1.e-5,nfailed(8),'dpsi/dt')
                endif
                call checkval(nptest,divBsymm,real(dBdtstore(maxBevol+1,1:nptest),kind=kind(divBsymm)),&
                              1.e-3,nfailed(9),'div B (symm)')
             endif
             ntests = ntests + 1
             if (all(nfailed==0)) npass = npass + 1
          endif
       enddo
    endif
    if (allocated(fxyzstore)) deallocate(fxyzstore)
    if (allocated(dBdtstore)) deallocate(dBdtstore)
#endif
 endif testinddts

 if (id==master) write(*,"(/,a)") '<-- DERIVS TEST COMPLETE'

contains

#ifdef IND_TIMESTEPS
subroutine reset_allactive
 !
 !--reset all particles to active for subsequent tests
 !
 do i=1,npart
    iphase(i) = isetphase(igas,iactive=.true.)
 enddo
 nactive = npart

end subroutine reset_allactive
#endif

subroutine set_active(npart,nactive,itype)
 integer, intent(in) :: npart, nactive, itype
 !
 !  set iphase for mixed active/inactive
 !
#ifdef IND_TIMESTEPS
 do i=1,npart
    if (i <= nactive) then
       iphase(i) = isetphase(itype,iactive=.true.)
    else
       iphase(i) = isetphase(itype,iactive=.false.)
    endif
 enddo
#endif
end subroutine set_active

!--------------------------------------
!+
!  reset all dissipation terms to zero
!+
!--------------------------------------
subroutine reset_dissipation_to_zero

 alpha = 0.
 alphau = 0.
 alphaB = 0.
 beta   = 0.
 if (maxalpha==maxp)  alphaind(:,:)  = 0.
 irealvisc = 0
 shearparam = 0.
 bulkvisc = 0.

end subroutine reset_dissipation_to_zero

!----------------------------------
!+
!  set vxyzu array using functions
!  ready for test suite
!+
!----------------------------------
subroutine set_velocity_and_energy
 integer :: iu
 iu = 4

 do i=1,npart
    vxyzu(1,i) = vx(xyzh(:,i))
    vxyzu(2,i) = vy(xyzh(:,i))
    vxyzu(3,i) = vz(xyzh(:,i))
    if (maxvxyzu >= 4) vxyzu(iu,i) = utherm(xyzh(:,i))
 enddo

end subroutine set_velocity_and_energy

!----------------------------------
!+
!  reset all MHD terms to zero
!+
!----------------------------------
subroutine reset_mhd_to_zero

 Bextx = 0.
 Bexty = 0.
 Bextz = 0.
 psidecayfac = 0.
 if (mhd) then
    Bevol(:,:) = 0.
    Bxyz(:,:)  = 0.
 endif
 if (use_dust) then
    dustfrac(:,:) = 0.
 endif

end subroutine reset_mhd_to_zero

end subroutine test_derivs

!----------------------------------------------------------------
!+
!  functional form for v and its derivatives
!+
!----------------------------------------------------------------
real function vx(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

 vx = 0.5/pi*dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function vx

real function vy(xyzhi)
 use boundary, only:xmin,dxbound,zmin,dzbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 vy = 0.5/pi*dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound) &
     - 0.5/pi*dzbound*sin(2.*pi*(xyzhi(3)-zmin)/dzbound)

end function vy

real function vz(xyzhi)
 use boundary, only:ymin,dybound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 vz = 0.05/pi*dybound*cos(4.*pi*(xyzhi(2)-ymin)/dybound)

end function vz

real function dvxdx(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 dvxdx = cos(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dvxdx

real function dvxdy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvxdy = 0.

end function dvxdy

real function dvxdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvxdz = 0.

end function dvxdz

real function dvydx(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 dvydx = cos(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dvydx

real function dvydy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvydy = 0.

end function dvydy

real function dvydz(xyzhi)
 use boundary, only:zmin,dzbound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

 dvydz = -cos(2.*pi*(xyzhi(3)-zmin)/dzbound)

end function dvydz

real function dvzdx(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvzdx = 0.

end function dvzdx

real function dvzdy(xyzhi)
 use boundary, only:ymin,dybound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

 dvzdy = -0.2*sin(4.*pi*(xyzhi(2)-ymin)/dybound)

end function dvzdy

real function dvzdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvzdz = 0.

end function dvzdz

!
!--second derivs
!
real function dvxdxdx(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

 dvxdxdx = -2.*pi/dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dvxdxdx

real function dvxdxdy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvxdxdy = 0.

end function dvxdxdy

real function dvxdxdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvxdxdz = 0.

end function dvxdxdz

real function dvxdydy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvxdydy = 0.

end function dvxdydy

real function dvxdydz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvxdydz = 0.

end function dvxdydz

real function dvxdzdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvxdzdz = 0.

end function dvxdzdz

real function dvydxdx(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

 dvydxdx = -2.*pi/dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dvydxdx

real function dvydxdy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvydxdy = 0.

end function dvydxdy

real function dvydxdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvydxdz = 0.

end function dvydxdz

real function dvydydy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvydydy = 0.

end function dvydydy

real function dvydydz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvydydz = 0.

end function dvydydz

real function dvydzdz(xyzhi)
 use boundary, only:zmin,dzbound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

 dvydzdz = 2.*pi/dzbound*sin(2.*pi*(xyzhi(3)-zmin)/dzbound)

end function dvydzdz

real function dvzdxdx(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvzdxdx = 0.

end function dvzdxdx

real function dvzdxdy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvzdxdy = 0.

end function dvzdxdy

real function dvzdxdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvzdxdz = 0.

end function dvzdxdz

real function dvzdydy(xyzhi)
 use boundary, only:ymin,dybound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

 dvzdydy = -0.8*pi/dybound*cos(4.*pi*(xyzhi(2)-ymin)/dybound)

end function dvzdydy

real function dvzdydz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvzdydz = 0.

end function dvzdydz

real function dvzdzdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dvzdzdz = 0.

end function dvzdzdz

!----------------------------------------------------------------
!+
!  functional form for divv
!+
!----------------------------------------------------------------
real function divvfunc(xyzhi)
 real, intent(in) :: xyzhi(4)

 divvfunc = dvxdx(xyzhi) + dvydy(xyzhi) + dvzdz(xyzhi)

end function divvfunc

!----------------------------------------------------------------
!+
!  functional form for curlv(x)
!+
!----------------------------------------------------------------
real function curlvfuncx(xyzhi)
 real, intent(in) :: xyzhi(4)

 curlvfuncx = dvzdy(xyzhi) - dvydz(xyzhi)

end function curlvfuncx

real function curlvfuncy(xyzhi)
 real, intent(in) :: xyzhi(4)

 curlvfuncy = dvxdz(xyzhi) - dvzdx(xyzhi)

end function curlvfuncy

real function curlvfuncz(xyzhi)
 real, intent(in) :: xyzhi(4)

 curlvfuncz = dvydx(xyzhi) - dvxdy(xyzhi)

end function curlvfuncz

real function graddivvfuncx(xyzhi)
 real, intent(in) :: xyzhi(4)

 graddivvfuncx = dvxdxdx(xyzhi) + dvydxdy(xyzhi) + dvzdxdz(xyzhi)

end function graddivvfuncx

real function graddivvfuncy(xyzhi)
 real, intent(in) :: xyzhi(4)

 graddivvfuncy = dvxdxdy(xyzhi) + dvydydy(xyzhi) + dvzdydz(xyzhi)

end function graddivvfuncy

real function graddivvfuncz(xyzhi)
 real, intent(in) :: xyzhi(4)

 graddivvfuncz = dvxdxdz(xyzhi) + dvydydz(xyzhi) + dvzdzdz(xyzhi)

end function graddivvfuncz

real function grad2vfuncx(xyzhi)
 real, intent(in) :: xyzhi(4)

 grad2vfuncx = dvxdxdx(xyzhi) + dvxdydy(xyzhi) + dvxdzdz(xyzhi)

end function grad2vfuncx

real function grad2vfuncy(xyzhi)
 real, intent(in) :: xyzhi(4)

 grad2vfuncy = dvydxdx(xyzhi) + dvydydy(xyzhi) + dvydzdz(xyzhi)

end function grad2vfuncy

real function grad2vfuncz(xyzhi)
 real, intent(in) :: xyzhi(4)

 grad2vfuncz = dvzdxdx(xyzhi) + dvzdydy(xyzhi) + dvzdzdz(xyzhi)

end function grad2vfuncz

!----------------------------------------------------------------
!+
!  functional form for ddivv/dt for use in Cullen/Dehnen switch
!+
!----------------------------------------------------------------
real function ddivvdtfunc(xyzhi)
 real, intent(in) :: xyzhi(4)

 ! for div a we assume a has been initialised using the
 ! same functions as v, so we can use divvfunc for this
 ddivvdtfunc = divvfunc(xyzhi) - (dvxdx(xyzhi)**2 + dvydy(xyzhi)**2 + dvzdz(xyzhi)**2 + &
    2.*(dvxdy(xyzhi)*dvydx(xyzhi) + dvxdz(xyzhi)*dvzdx(xyzhi) + dvydz(xyzhi)*dvzdy(xyzhi)))

end function ddivvdtfunc

real function alphalocfunc(xyzhi)
 use options,      only:alpha,alphamax
 use eos,          only:gamma,polyk
 use densityforce, only:get_alphaloc
 real, intent(in) :: xyzhi(4)
 real :: ddivvdti,spsoundi,xi_limiter,fac,curlv2

 ddivvdti = ddivvdtfunc(xyzhi)
 fac = -max(-divvfunc(xyzhi),0.)**2
 curlv2 = curlvfuncx(xyzhi)**2 + curlvfuncy(xyzhi)**2 + curlvfuncz(xyzhi)**2
 if (fac + curlv2 > 0.) then
    xi_limiter = fac/(fac + curlv2)
 else
    xi_limiter = 1.
 endif
 if (gamma < 1.0001) then
    spsoundi = sqrt(polyk)
 else
    spsoundi = sqrt(gamma*(gamma-1.)*uthermconst(xyzhi))
 endif
 alphalocfunc = get_alphaloc(ddivvdti,spsoundi,xyzhi(4),xi_limiter,alpha,alphamax)

end function alphalocfunc

!----------------------------------------------------------------
!+
!  functional form for dh/dt
!+
!----------------------------------------------------------------
real function dhdtfunc(xyzhi)
 use physcon,  only:pi
 use part,     only:rhoh,dhdrho
 real, intent(in) :: xyzhi(4)
 real :: drhodti

 drhodti = -rhoh(xyzhi(4),massoftype(1))*divvfunc(xyzhi)
 dhdtfunc = dhdrho(xyzhi(4),massoftype(1))*drhodti

end function dhdtfunc
!----------------------------------------------------------------
!+
!  functional form for du/dt = -P/rho (div v)
!  not including the term out the front = (gamma-1)*u
!+
!----------------------------------------------------------------
real function dudtfunc(xyzhi)
 real, intent(in) :: xyzhi(4)

 dudtfunc = -divvfunc(xyzhi)

end function dudtfunc

!----------------------------------------------------------------
!+
!  functional form for thermal energy
!+
!----------------------------------------------------------------
real function utherm(xyzhi)
 use boundary, only:xmin,dxbound,ymin,dybound,zmin,dzbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 utherm = 0.5/pi*(3. + sin(2.*pi*(xyzhi(1)-xmin)/dxbound) &
                     + cos(2.*pi*(xyzhi(2)-ymin)/dybound) &
                     + sin(2.*pi*(xyzhi(3)-zmin)/dzbound))

end function utherm

!----------------------------------------------------------------
!+
!  functional form for thermal energy
!  (constant => use for AV to get constant spsound)
!+
!----------------------------------------------------------------
real function uthermconst(xyzhi)
 real, intent(in) :: xyzhi(4)

 uthermconst = 4.0

end function uthermconst

!----------------------------------------------------------------
!+
!  functional form for hydrodynamic forces
!+
!----------------------------------------------------------------
real function dudx(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 dudx = 1./dxbound*cos(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dudx

real function dudy(xyzhi)
 use boundary, only:ymin,dybound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 dudy = -1./dybound*sin(2.*pi*(xyzhi(2)-ymin)/dybound)

end function dudy

real function dudz(xyzhi)
 use boundary, only:zmin,dzbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 dudz = 1./dzbound*cos(2.*pi*(xyzhi(3)-zmin)/dzbound)

end function dudz

real function del2u(xyzhi)
 use boundary, only:xmin,ymin,zmin,dxbound,dybound,dzbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)
 real :: dudxdx,dudydy,dudzdz

 dudxdx = -2.*pi/dxbound**2*sin(2.*pi*(xyzhi(1)-xmin)/dxbound)
 dudydy = -2.*pi/dybound**2*cos(2.*pi*(xyzhi(2)-ymin)/dybound)
 dudzdz = -2.*pi/dzbound**2*sin(2.*pi*(xyzhi(3)-zmin)/dzbound)
 del2u  = dudxdx + dudydy + dudzdz

end function del2u

real function forcefuncx(xyzhi)
 use eos,      only:gamma
 real, intent(in) :: xyzhi(4)

 ! fx = -(grad P)_x / rho,
 ! P = (gamma-1)*rho*u,
 ! grad P = (gamma-1)*rho*du/dx  (rho = const)

 forcefuncx = -(gamma-1.)*dudx(xyzhi)

end function forcefuncx

real function forcefuncy(xyzhi)
 use eos,      only:gamma
 real, intent(in) :: xyzhi(4)

 ! fy = -(grad P)_y / rho,
 ! P = (gamma-1)*rho*u,
 ! grad P = (gamma-1)*rho*du/dy  (rho = const)

 forcefuncy = -(gamma-1.)*dudy(xyzhi)

end function forcefuncy

real function forcefuncz(xyzhi)
 use eos,      only:gamma
 real, intent(in) :: xyzhi(4)

 ! fz = -(grad P)_z / rho,
 ! P = (gamma-1)*rho*u,
 ! grad P = (gamma-1)*rho*du/dz  (rho = const)

 forcefuncz = -(gamma-1.)*dudz(xyzhi)

end function forcefuncz

!----------------------------------------------------------------
!+
!  functional form for hydrodynamic + AV forces
!  (see e.g. Lodato & Price (2010) for translation
!   of SPH AV term into Navier-Stokes terms)
!+
!----------------------------------------------------------------
real function forceavx(xyzhi)
 use eos,     only:gamma,polyk
 use options, only:alpha
 real, intent(in) :: xyzhi(4)
 real :: spsoundi,fac,coeff1,coeff2

 if (gamma < 1.0001) then
    spsoundi = sqrt(polyk)
 else
    spsoundi = sqrt(gamma*(gamma-1.)*uthermconst(xyzhi))
 endif
#ifdef DISC_VISCOSITY
 fac = alpha*spsoundi*xyzhi(4)
#else
 fac = 0.5*alpha*spsoundi*xyzhi(4)
#endif
 coeff1 = fac*0.1
 coeff2 = fac*0.2
 forceavx = &
    coeff1*(dvxdxdx(xyzhi) + dvxdydy(xyzhi) + dvxdzdz(xyzhi)) &  ! del^2 v
   +coeff2*(dvxdxdx(xyzhi) + dvydxdy(xyzhi) + dvzdxdz(xyzhi))    ! grad (div v)

end function forceavx

real function forceavy(xyzhi)
 use eos,     only:gamma,polyk
 use options, only:alpha !,beta
 real, intent(in) :: xyzhi(4)
 real :: spsoundi,fac,coeff1,coeff2

 if (gamma < 1.0001) then
    spsoundi = sqrt(polyk)
 else
    spsoundi = sqrt(gamma*(gamma-1.)*uthermconst(xyzhi))
 endif
#ifdef DISC_VISCOSITY
 fac = alpha*spsoundi*xyzhi(4)
#else
 fac = 0.5*alpha*spsoundi*xyzhi(4)
#endif
 coeff1 = fac*0.1
 coeff2 = fac*0.2
 forceavy =  &
    coeff1*(dvydxdx(xyzhi) + dvydydy(xyzhi) + dvydzdz(xyzhi)) &  ! del^2 v
   +coeff2*(dvxdxdy(xyzhi) + dvydydy(xyzhi) + dvzdydz(xyzhi))    ! grad (div v)

end function forceavy

real function forceavz(xyzhi)
 use eos,     only:gamma,polyk
 use options, only:alpha
 real, intent(in) :: xyzhi(4)
 real :: spsoundi,fac,coeff1,coeff2

 if (gamma < 1.0001) then
    spsoundi = sqrt(polyk)
 else
    spsoundi = sqrt(gamma*(gamma-1.)*uthermconst(xyzhi))
 endif
#ifdef DISC_VISCOSITY
 fac = alpha*spsoundi*xyzhi(4)
#else
 fac = 0.5*alpha*spsoundi*xyzhi(4)
#endif
 coeff1 = fac*0.1
 coeff2 = fac*0.2
 forceavz =  &
    coeff1*(dvzdxdx(xyzhi) + dvzdydy(xyzhi) + dvzdzdz(xyzhi)) &  ! del^2 v
   +coeff2*(dvxdxdz(xyzhi) + dvydydz(xyzhi) + dvzdzdz(xyzhi))    ! grad (div v)

end function forceavz

!----------------------------------------------------------------
!+
!  functional form for viscous forces
!  (assuming constant viscosity parameters and constant density)
!+
!----------------------------------------------------------------
real function forceviscx(xyzhi)
 use viscosity, only:bulkvisc,shearfunc
 real, intent(in) :: xyzhi(4)
 real :: eta

 eta = shearfunc(xyzhi(1),xyzhi(2),xyzhi(3),0.)
 forceviscx = 1.*( &
                  eta*(dvxdxdx(xyzhi) + dvxdydy(xyzhi) + dvxdzdz(xyzhi)) &  ! del^2 v
 +(bulkvisc + eta/3.)*(dvxdxdx(xyzhi) + dvydxdy(xyzhi) + dvzdxdz(xyzhi)))   ! grad (div v)

end function forceviscx

real function forceviscy(xyzhi)
 use viscosity, only:bulkvisc,shearfunc
 real, intent(in) :: xyzhi(4)
 real :: eta

 eta = shearfunc(xyzhi(1),xyzhi(2),xyzhi(3),0.)
 forceviscy = 1.*( &
                  eta*(dvydxdx(xyzhi) + dvydydy(xyzhi) + dvydzdz(xyzhi)) &  ! del^2 v
 +(bulkvisc + eta/3.)*(dvxdxdy(xyzhi) + dvydydy(xyzhi) + dvzdydz(xyzhi)))   ! grad (div v)

end function forceviscy

real function forceviscz(xyzhi)
 use viscosity, only:bulkvisc,shearfunc
 real, intent(in) :: xyzhi(4)
 real :: eta

 eta = shearfunc(xyzhi(1),xyzhi(2),xyzhi(3),0.)
 forceviscz = 1.*( &
                  eta*(dvzdxdx(xyzhi) + dvzdydy(xyzhi) + dvzdzdz(xyzhi)) &  ! del^2 v
 +(bulkvisc + eta/3.)*(dvxdxdz(xyzhi) + dvydydz(xyzhi) + dvzdzdz(xyzhi)))   ! grad (div v)

end function forceviscz

!
!--spatial derivatives of shear viscosity parameter
!
real function detadx(xyzhi)
 real, intent(in) :: xyzhi(4)

 detadx = 0.

end function detadx

!--spatial derivative of shear viscosity parameter
real function detady(xyzhi)
 real, intent(in) :: xyzhi(4)

 detady = 0.

end function detady

!--spatial derivative of shear viscosity parameter
real function detadz(xyzhi)
 real, intent(in) :: xyzhi(4)

 detadz = 0.

end function detadz

!------------------
!+
!  Strain tensor
!+
!-----------------
real function sxx(xyzhi)
 real, intent(in) :: xyzhi(4)

 sxx = 2.*dvxdx(xyzhi)

end function sxx

real function sxy(xyzhi)
 real, intent(in) :: xyzhi(4)

 sxy = dvxdy(xyzhi) + dvydx(xyzhi)

end function sxy

real function sxz(xyzhi)
 real, intent(in) :: xyzhi(4)

 sxz = dvxdz(xyzhi) + dvzdx(xyzhi)

end function sxz

real function syy(xyzhi)
 real, intent(in) :: xyzhi(4)

 syy = 2.*dvydy(xyzhi)

end function syy

real function syz(xyzhi)
 real, intent(in) :: xyzhi(4)

 syz = dvydz(xyzhi) + dvzdy(xyzhi)

end function syz

real function szz(xyzhi)
 real, intent(in) :: xyzhi(4)

 szz = 2.*dvzdz(xyzhi)

end function szz

!----------------------------------------------------------------
!+
!  functional form for vector potential
!+
!----------------------------------------------------------------
real function Ax(xyzhi)
 use boundary, only:ymin,dybound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

! Bz = 1.5/pi*dybound*cos(2.*pi*(xyzhi(2)-ymin)/dybound) + 1.0
!
! Bz = dAy/dx - dAx/dy
!
 Ax = -7.5*dybound**2/(pi**2)*sin(2.*pi*(xyzhi(2)-ymin)/dybound)

end function Ax

real function Ay(xyzhi)
 use boundary, only:zmin,dzbound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

! NB this is non-zero div B
! Bx = 0.5/pi*dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound) &
!    - 0.5/pi*dzbound*cos(2.*pi*(xyzhi(3)-zmin)/dzbound) + 2.0
!
! Bx = dAz/dy - dAy/dz
!
 Ay = 2.5*dzbound**2/(pi**2)*sin(2.*pi*(xyzhi(3)-zmin)/dzbound)

end function Ay

real function Az(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon, only:pi
 real, intent(in) :: xyzhi(4)

! By = 0.5/pi*dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound) + 3.0
!
! By = dAx/dz - dAz/dx
!
 Az = 2.5*dxbound**2/(pi**2)*cos(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function Az

!----------------------------------------------------------------
!+
!  functional form for B and its derivatives
!+
!----------------------------------------------------------------
real function Bx(xyzhi)
 use boundary, only:xmin,dxbound,zmin,dzbound
 use physcon,  only:pi
 use part,     only:Bextx
 real, intent(in) :: xyzhi(4)

 Bx = -5./pi*dzbound*cos(2.*pi*(xyzhi(3)-zmin)/dzbound) + Bextx
! NB this is non-zero div B
 Bx = Bx + 0.5/pi*dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function Bx

real function By(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon,  only:pi
 use part,     only:Bexty
 real, intent(in) :: xyzhi(4)

 By = 5./pi*dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound) + Bexty

end function By

real function Bz(xyzhi)
 use boundary, only:ymin,dybound
 use physcon,  only:pi
 use part,     only:Bextz
 real, intent(in) :: xyzhi(4)

 Bz = 15./pi*dybound*cos(2.*pi*(xyzhi(2)-ymin)/dybound) + Bextz

end function Bz

real function dBxdx(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 ! Bx = 0.5/pi*dxbound*sin(2.*pi*(xyzh(1,i)-xmin)/dxbound)

 dBxdx = cos(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dBxdx

real function dBxdy(xyzhi)
 real, intent(in) :: xyzhi(4)

 ! Bx = 0.5/pi*dxbound*sin(2.*pi*(xyzh(1,i)-xmin)/dxbound)

 dBxdy = 0.

end function dBxdy

real function dBxdz(xyzhi)
 use boundary, only:zmin,dzbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 ! Bx = 0.5/pi*dxbound*sin(2.*pi*(xyzh(1,i)-xmin)/dxbound)
 !    - 0.5/pi*dzbound*cos(2.*pi*(xyzhi(3)-zmin)/dzbound)
 dBxdz = 10.*sin(2.*pi*(xyzhi(3)-zmin)/dzbound)

end function dBxdz

real function dBydx(xyzhi)
 use boundary, only:xmin,dxbound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 ! By = 0.5/pi*dxbound*sin(2.*pi*(xyzh(1,i)-xmin)/dxbound)

 dBydx = 10.*cos(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dBydx

real function dBydy(xyzhi)
 real, intent(in) :: xyzhi(4)

 ! By = 0.5/pi*dxbound*sin(2.*pi*(xyzh(1,i)-xmin)/dxbound)

 dBydy = 0.

end function dBydy

real function dBydz(xyzhi)
 real, intent(in) :: xyzhi(4)

 ! By = 0.5/pi*dxbound*sin(2.*pi*(xyzh(1,i)-xmin)/dxbound)

 dBydz = 0.

end function dBydz

real function dBzdx(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBzdx = 0.

end function dBzdx

real function dBzdy(xyzhi)
 use boundary, only:ymin,dybound
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 ! Bz = 15./pi*dybound*cos(2.*pi*(xyzhi(2)-ymin)/dybound) + Bextz
 dBzdy = -30.*sin(2.*pi*(xyzhi(2)-ymin)/dybound)

end function dBzdy

real function dBzdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBzdz = 0.

end function dBzdz

!
!--second derivatives of B (used to test resistivity)
!
real function dBxdxdx(xyzhi)
 use boundary, only:dxbound,xmin
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 dBxdxdx = -2.*pi/dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dBxdxdx

real function dBxdydy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBxdydy = 0.

end function dBxdydy

real function dBxdzdz(xyzhi)
 use boundary, only:dzbound,zmin
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

! dBxdz = 10.*sin(2.*pi*(xyzhi(3)-zmin)/dzbound)
 dBxdzdz = 20.*pi/dzbound*cos(2.*pi*(xyzhi(3)-zmin)/dzbound)

end function dBxdzdz

real function dBydxdx(xyzhi)
 use boundary, only:dxbound,xmin
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 !dBydx = 10.*cos(2.*pi*(xyzhi(1)-xmin)/dxbound)
 dBydxdx = -20.*pi/dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dBydxdx

real function dBydydy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBydydy = 0.

end function dBydydy

real function dBydzdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBydzdz = 0.

end function dBydzdz

real function dBzdxdx(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBzdxdx = 0.

end function dBzdxdx

real function dBzdydy(xyzhi)
 use boundary, only:dybound,ymin
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 !dBzdy = -30.*sin(2.*pi*(xyzhi(2)-ymin)/dybound)
 dBzdydy = -60.*pi/dybound*cos(2.*pi*(xyzhi(2)-ymin)/dybound)

end function dBzdydy

real function dBzdzdz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBzdzdz = 0.

end function dBzdzdz

real function dBxdtresist(xyzhi)
 use options, only:alphaB
 use part,    only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: vsig, rho1i

 vsig = 0. !valfven(xyzhi(1))
 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dBxdtresist = rho1i * (0.5*alphaB*xyzhi(4)*(vsig*(dBxdxdx(xyzhi) + dBxdydy(xyzhi) + dBxdzdz(xyzhi)) + &
   0.*(dvalfvendx(xyzhi)*dBxdx(xyzhi) + dvalfvendy(xyzhi)*dBxdy(xyzhi) + dvalfvendz(xyzhi)*dBxdz(xyzhi))))

end function dBxdtresist

real function dBydtresist(xyzhi)
 use options, only:alphaB
 use part,    only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: vsig, rho1i

 vsig = 0. !valfven(xyzhi)
 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dBydtresist = rho1i * (0.5*alphaB*xyzhi(4)*(vsig*(dBydxdx(xyzhi) + dBydydy(xyzhi) + dBydzdz(xyzhi)) + &
   0.*(dvalfvendx(xyzhi)*dBydx(xyzhi) + dvalfvendy(xyzhi)*dBydy(xyzhi) + dvalfvendz(xyzhi)*dBydz(xyzhi))))

end function dBydtresist

real function dBzdtresist(xyzhi)
 use options, only:alphaB
 use part,    only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: vsig, rho1i

 vsig = 0. !valfven(xyzhi)
 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dBzdtresist = 0.5*alphaB*xyzhi(4)*(vsig*(dBzdxdx(xyzhi) + dBzdydy(xyzhi) + dBzdzdz(xyzhi)) + &
   0.*(dvalfvendx(xyzhi)*dBzdx(xyzhi) + dvalfvendy(xyzhi)*dBzdy(xyzhi) + dvalfvendz(xyzhi)*dBzdz(xyzhi)))

end function dBzdtresist

!
!--functional form of div B
!
real function divBfunc(xyzhi)
 real, intent(in) :: xyzhi(4)

 divBfunc = dBxdx(xyzhi) + dBydy(xyzhi) + dBzdz(xyzhi)

end function divBfunc

!----------------------------------------------------------------
!+
!  functional form for curl B
!+
!----------------------------------------------------------------
real function curlBfuncx(xyzhi)
 real, intent(in) :: xyzhi(4)

 curlBfuncx = dBzdy(xyzhi) - dBydz(xyzhi)

end function curlBfuncx

real function curlBfuncy(xyzhi)
 real, intent(in) :: xyzhi(4)

 curlBfuncy = dBxdz(xyzhi) - dBzdx(xyzhi)

end function curlBfuncy

real function curlBfuncz(xyzhi)
 real, intent(in) :: xyzhi(4)

 curlBfuncz = dBydx(xyzhi) - dBxdy(xyzhi)

end function curlBfuncz

!----------------------------------------------------------------
!+
!  functional form for vector potential time derivatives
!+
!----------------------------------------------------------------
real function dAxdt(xyzhi)
 use part, only:Bexty,Bextz
 real, intent(in) :: xyzhi(4)

 dAxdt = -(Ax(xyzhi)*dvxdx(xyzhi) + Ay(xyzhi)*dvydx(xyzhi) + Az(xyzhi)*dvzdx(xyzhi)) &
         + vy(xyzhi)*Bextz - vz(xyzhi)*Bexty

end function dAxdt

real function dAydt(xyzhi)
 use part, only:Bextx,Bextz
 real, intent(in) :: xyzhi(4)

 dAydt = -(Ax(xyzhi)*dvxdy(xyzhi) + Ay(xyzhi)*dvydy(xyzhi) + Az(xyzhi)*dvzdy(xyzhi)) &
         + vz(xyzhi)*Bextx - vx(xyzhi)*Bextz

end function dAydt

real function dAzdt(xyzhi)
 use part, only:Bextx,Bexty
 real, intent(in) :: xyzhi(4)

 dAzdt = -(Ax(xyzhi)*dvxdz(xyzhi) + Ay(xyzhi)*dvydz(xyzhi) + Az(xyzhi)*dvzdz(xyzhi)) &
         + vx(xyzhi)*Bexty - vy(xyzhi)*Bextx

end function dAzdt

!----------------------------------------------------------------
!+
!  functional form for dB/dt
!+
!----------------------------------------------------------------
real function dBxdt(xyzhi)
 use part,    only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: rho1i

 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dBxdt = rho1i * (Bx(xyzhi)*dvxdx(xyzhi) + By(xyzhi)*dvxdy(xyzhi) &
             + Bz(xyzhi)*dvxdz(xyzhi))! - Bx(xyzhi)*divvfunc(xyzhi))

end function dBxdt

real function dBydt(xyzhi)
 use part,    only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: rho1i

 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dBydt = rho1i * (Bx(xyzhi)*dvydx(xyzhi) + By(xyzhi)*dvydy(xyzhi) &
             + Bz(xyzhi)*dvydz(xyzhi))! - By(xyzhi)*divvfunc(xyzhi))

end function dBydt

real function dBzdt(xyzhi)
 use part,    only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: rho1i

 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dBzdt = rho1i * (Bx(xyzhi)*dvzdx(xyzhi) + By(xyzhi)*dvzdy(xyzhi) &
             + Bz(xyzhi)*dvzdz(xyzhi))! - Bz(xyzhi)*divvfunc(xyzhi))

end function dBzdt

!----------------------------------------------------------------
!+
!  functional form for MHD forces
!+
!----------------------------------------------------------------
real function forcemhdx(xyzhi)
 real, intent(in) :: xyzhi(4)
 real :: fisox,fanisox

 fisox = Bx(xyzhi)*dBxdx(xyzhi) + By(xyzhi)*dBydx(xyzhi) + Bz(xyzhi)*dBzdx(xyzhi)
 fanisox = Bx(xyzhi)*dBxdx(xyzhi) + By(xyzhi)*dBxdy(xyzhi) + Bz(xyzhi)*dBxdz(xyzhi) !&
 !+ Bx(xyzhi)*divBfunc(xyzhi)

 forcemhdx = -1./rhozero*(fisox - fanisox)

end function forcemhdx

real function forcemhdy(xyzhi)
 real, intent(in) :: xyzhi(4)
 real :: fisoy,fanisoy

 fisoy = Bx(xyzhi)*dBxdy(xyzhi) + By(xyzhi)*dBydy(xyzhi) + Bz(xyzhi)*dBzdy(xyzhi)
 fanisoy = Bx(xyzhi)*dBydx(xyzhi) + By(xyzhi)*dBydy(xyzhi) + Bz(xyzhi)*dBydz(xyzhi) !&
 !+ By(xyzhi)*divBfunc(xyzhi)

 forcemhdy = -1./rhozero*(fisoy - fanisoy)

end function forcemhdy

real function forcemhdz(xyzhi)
 real, intent(in) :: xyzhi(4)
 real :: fisoz,fanisoz

 fisoz = Bx(xyzhi)*dBxdz(xyzhi) + By(xyzhi)*dBydz(xyzhi) + Bz(xyzhi)*dBzdz(xyzhi)
 fanisoz = Bx(xyzhi)*dBzdx(xyzhi) + By(xyzhi)*dBzdy(xyzhi) + Bz(xyzhi)*dBzdz(xyzhi) !&
 !+ Bz(xyzhi)*divBfunc(xyzhi)

 forcemhdz = -1./rhozero*(fisoz - fanisoz)

end function forcemhdz

!----------------------------------------------------------------
!+
!  Alfven speed and derivatives
!+
!----------------------------------------------------------------
real function valfven(xyzhi)
 use part, only:rhoh
 real, intent(in) :: xyzhi(4)

 valfven = sqrt((Bx(xyzhi)**2 + By(xyzhi)**2 + Bz(xyzhi)**2)/rhoh(xyzhi(4),massoftype(1)))

end function valfven

real function dvalfvendx(xyzhi)
 use part, only:rhoh
 real, intent(in) :: xyzhi(4)
 real :: rhoi, B2i

 rhoi = rhoh(xyzhi(4),massoftype(1))
 B2i  = Bx(xyzhi)**2 + By(xyzhi)**2 + Bz(xyzhi)**2
 dvalfvendx = 1./sqrt(rhoi*B2i)* &
             (Bx(xyzhi)*dBxdx(xyzhi) + By(xyzhi)*dBydx(xyzhi) + Bz(xyzhi)*dBzdx(xyzhi))

end function dvalfvendx

real function dvalfvendy(xyzhi)
 use part, only:rhoh
 real, intent(in) :: xyzhi(4)
 real :: rhoi, B2i

 rhoi = rhoh(xyzhi(4),massoftype(1))
 B2i  = Bx(xyzhi)**2 + By(xyzhi)**2 + Bz(xyzhi)**2
 dvalfvendy = 1./sqrt(rhoi*B2i)* &
             (Bx(xyzhi)*dBxdy(xyzhi) + By(xyzhi)*dBydy(xyzhi) + Bz(xyzhi)*dBzdy(xyzhi))

end function dvalfvendy

real function dvalfvendz(xyzhi)
 use part, only:rhoh
 real, intent(in) :: xyzhi(4)
 real :: rhoi, B2i

 rhoi = rhoh(xyzhi(4),massoftype(1))
 B2i  = Bx(xyzhi)**2 + By(xyzhi)**2 + Bz(xyzhi)**2
 dvalfvendz = 1./sqrt(rhoi*B2i)* &
             (Bx(xyzhi)*dBxdz(xyzhi) + By(xyzhi)*dBydz(xyzhi) + Bz(xyzhi)*dBzdz(xyzhi))

end function dvalfvendz

!----------------------------------------------------------------
!+
!  functional form for Psi and derivatives in test of div B cleaning
!+
!----------------------------------------------------------------
real function psi(xyzhi)
 use boundary, only:dxbound,dybound,dzbound,xmin,ymin,zmin
 use physcon,  only:pi
 real, intent(in) :: xyzhi(4)

 psi = 0.5/pi*dxbound*sin(2.*pi*(xyzhi(1)-xmin)/dxbound) &
     - 0.5/pi*dzbound*cos(2.*pi*(xyzhi(3)-zmin)/dzbound) &
     + 0.5/pi*dybound*sin(2.*pi*(xyzhi(2)-ymin)/dybound)

end function psi

real function dpsidt(xyzhi)
 use options, only:ieos,psidecayfac
 use eos,     only:equationofstate
 use part,    only:rhoh
 real, intent(in) :: xyzhi(4)
 real :: vsig,spsoundi,ponrhoi

 call equationofstate(ieos,ponrhoi,spsoundi,rhoh(xyzhi(4),massoftype(1)),xyzhi(1),xyzhi(2),xyzhi(3))
 vsig = sqrt(valfven(xyzhi)**2 + spsoundi**2)
 dpsidt = -vsig**2*divBfunc(xyzhi) - psi(xyzhi)*psidecayfac*vsig/xyzhi(4) &
          -0.5*psi(xyzhi)*divvfunc(xyzhi)

end function dpsidt

real function dpsidx(xyzhi)
 use boundary, only:dxbound,xmin
 use physcon,  only:pi
 use part,     only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: rho1i

 !--minus grad psi
 !  updated to be -1/rho grad psi (for B/rho evolution)
 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dpsidx = dBxdt(xyzhi) - rho1i * cos(2.*pi*(xyzhi(1)-xmin)/dxbound)

end function dpsidx

real function dpsidy(xyzhi)
 use boundary, only:dybound,ymin
 use physcon,  only:pi
 use part,     only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: rho1i

 !--minus grad psi
 !  updated to be -1/rho grad psi (for B/rho evolution)
 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dpsidy = dBydt(xyzhi) - rho1i * cos(2.*pi*(xyzhi(2)-ymin)/dybound)

end function dpsidy

real function dpsidz(xyzhi)
 use boundary, only:dzbound,zmin
 use physcon,  only:pi
 use part,     only:rhoh,massoftype,igas
 real, intent(in) :: xyzhi(4)
 real :: rho1i

 !--minus grad psi
 !  updated to be -1/rho grad psi (for B/rho evolution)
 rho1i = 1.0/rhoh(xyzhi(4),massoftype(igas))
 dpsidz = dBzdt(xyzhi) - rho1i * sin(2.*pi*(xyzhi(3)-zmin)/dzbound)

end function dpsidz

!----------------------------------------------------------------
!+
!  functional form for (dB/dt)_ambipolar
!+
!----------------------------------------------------------------
real function dBambix(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBambix = 0.

end function dBambix

real function dBambiy(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBambiy = 0.

end function dBambiy

real function dBambiz(xyzhi)
 real, intent(in) :: xyzhi(4)

 dBambiz = 0.

end function dBambiz

!----------------------------------------------------------------
!+
!  functional form for one-fluid dust derivatives
!+
!----------------------------------------------------------------
real function dustfrac_func(xyzhi)
 use physcon,  only:pi
 use part,     only:ndusttypes
 use boundary, only:dxbound,dybound,dzbound,xmin,ymin,zmin
 real, intent(in) :: xyzhi(4)

 dustfrac_func = 0.5/pi*(0.5 + 0.01*sin(4.*pi*(xyzhi(1)-xmin)/dxbound)  &
                             + 0.02*sin(2.*pi*(xyzhi(2)-ymin)/dybound)  &
                             + 0.05*cos(4.*pi*(xyzhi(3)-zmin)/dzbound)) &
                 *1./real(ndusttypes)

end function dustfrac_func

real function ddustevoldx(xyzhi)
 use physcon,  only:pi
 use part,     only:ndusttypes
 use boundary, only:dxbound,xmin
 real, intent(in) :: xyzhi(4)

 ddustevoldx = 2.*0.01/dxbound*cos(4.*pi*(xyzhi(1)-xmin)/dxbound) &
               *1./real(ndusttypes)

end function ddustevoldx

real function ddustevoldy(xyzhi)
 use physcon,  only:pi
 use part,     only:ndusttypes
 use boundary, only:dybound,ymin
 real, intent(in) :: xyzhi(4)

 ddustevoldy = 0.02/dybound*cos(2.*pi*(xyzhi(2)-ymin)/dybound) &
               *1./real(ndusttypes)

end function ddustevoldy

real function ddustevoldz(xyzhi)
 use physcon,  only:pi
 use part,     only:ndusttypes
 use boundary, only:dzbound,zmin
 real, intent(in) :: xyzhi(4)

 ddustevoldz = -2.*0.05/dzbound*sin(4.*pi*(xyzhi(3)-zmin)/dzbound) &
               *1./real(ndusttypes)

end function ddustevoldz

real function del2dustfrac(xyzhi)
 use physcon,  only:pi
 use part,     only:ndusttypes
 use boundary, only:dxbound,dybound,dzbound,xmin,ymin,zmin
 real, intent(in) :: xyzhi(4)

 del2dustfrac = (-8.*pi/dxbound**2*0.01*sin(4.*pi*(xyzhi(1)-xmin)/dxbound)  &
                 -2.*pi/dybound**2*0.02*sin(2.*pi*(xyzhi(2)-ymin)/dybound)  &
                 -8.*pi/dzbound**2*0.05*cos(4.*pi*(xyzhi(3)-zmin)/dzbound)) &
                *1./real(ndusttypes)

end function del2dustfrac

real function ddustevol_func(xyzhi)
 use eos,  only:gamma
 use part, only:ndusttypes
#ifdef DUST
 use dust, only:get_ts,idrag,K_code
#endif
 use part, only:rhoh
 real, intent(in) :: xyzhi(4)
 real :: dustfraci,uui,pri,tsi
 real :: gradu(3),gradeps(3),gradsumeps(3),gradp(3),gradts(3),gradepsts(3)
 real :: rhoi,rhogasi,rhodusti,spsoundi,del2P,du_dot_de,si
 real :: dustfracisum,del2dustfracsum
#ifdef DUST
 integer :: iregime
#endif

 rhoi       = rhoh(xyzhi(4),massoftype(1))
 dustfraci  = dustfrac_func(xyzhi)
 dustfracisum = real(ndusttypes)*dustfraci
 rhogasi    = (1. - dustfracisum)*rhoi
 rhodusti   = dustfracisum*rhoi
 uui        = utherm(xyzhi)
 pri        = (gamma-1.)*rhogasi*uui
 spsoundi   = gamma*pri/rhogasi

 gradu(1)   = dudx(xyzhi)
 gradu(2)   = dudy(xyzhi)
 gradu(3)   = dudz(xyzhi)
 gradeps(1) = ddustevoldx(xyzhi)
 gradeps(2) = ddustevoldy(xyzhi)
 gradeps(3) = ddustevoldz(xyzhi)
 gradsumeps = gradeps*real(ndusttypes)
 du_dot_de  = dot_product(gradu,gradsumeps)
 gradp(:)   = (gamma-1.)*(rhogasi*gradu - rhoi*uui*gradsumeps)
 del2dustfracsum = real(ndusttypes)*del2dustfrac(xyzhi)
 del2P = (gamma-1.)*rhoi*((1. - dustfracisum)*del2u(xyzhi) - 2.*du_dot_de - uui*del2dustfracsum)

 tsi   = 0.
#ifdef DUST
 call get_ts(idrag,grainsizek,graindensk,rhogasi,rhodusti,spsoundi,0.,tsi,iregime)
 !
 ! grad(ts) = grad((1-eps)*eps*rho/K_code)
 !          = rho/K_code*(1-2*eps)*grad(eps)          ! note the absence of eps_k
 !
 gradts(:) = rhoi/K_code*(1. - 2.*dustfracisum)*gradsumeps(:)
#else
 gradts(:) = 0.
#endif
 !
 ! deps_k/dt  = -1/rho \nabla.(eps_k ts (grad P))     ! note the presence of eps_k
 !            = -1/rho [eps_k ts \del^2 P + grad(eps_k ts).grad P]
 !            = -1/rho [eps_k ts \del^2 P + (eps_k*grad(ts) + ts*grad(eps_k)).grad P]
 !
 gradepsts(:) = dustfraci*gradts(:) + tsi*gradeps(:)

 !ddustevol_func = -1./rhoi*(dustfraci*tsi*del2P + dot_product(gradp,gradepsts))

!------------------------------------------------
!--sqrt(rho*epsilon) method
! si = sqrt(dustfraci*rhoi)
! ddustevol_func = -0.5/si*(dustfraci*tsi*del2P + dot_product(gradp,gradepsts)) - 0.5*si*divvfunc(xyzhi)
!------------------------------------------------
!--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
 si = sqrt(dustfraci/(1.-dustfraci))
 ddustevol_func = -0.5*((dustfraci*tsi*del2P + dot_product(gradp,gradepsts))/(rhoi*si*(1.-dustfraci)**2.))
!------------------------------------------------
!--asin(sqrt(epsilon)) method
! si = asin(sqrt(dustfraci))
! ddustevol_func = -0.5/(rhoi*sin(si)*cos(si))*(dustfraci*tsi*del2P + dot_product(gradp,gradepsts))
!------------------------------------------------

end function ddustevol_func

real function dudtdust_func(xyzhi)
 use eos,  only:gamma
 use part, only:ndusttypes
#ifdef DUST
 use dust, only:get_ts,idrag
#endif
 use part, only:rhoh
 real, intent(in) :: xyzhi(4)
 real :: dustfraci,uui,pri,tsi
 real :: gradp(3),gradu(3),gradeps(3),gradsumeps(3)
 real :: rhoi,rhogasi,rhodusti,spsoundi
 real :: dustfracisum
#ifdef DUST
 integer :: iregime
#endif

 rhoi       = rhoh(xyzhi(4),massoftype(1))
 dustfraci  = dustfrac_func(xyzhi)
 dustfracisum = real(ndusttypes)*dustfraci
 rhogasi    = (1. - dustfracisum)*rhoi
 rhodusti   = dustfracisum*rhoi
 uui        = utherm(xyzhi)
 gradu(1)   = dudx(xyzhi)
 gradu(2)   = dudy(xyzhi)
 gradu(3)   = dudz(xyzhi)
 gradeps(1) = ddustevoldx(xyzhi)
 gradeps(2) = ddustevoldy(xyzhi)
 gradeps(3) = ddustevoldz(xyzhi)
 gradsumeps = real(ndusttypes)*gradeps
 pri        = (gamma-1.)*rhogasi*uui
 spsoundi   = gamma*pri/rhogasi
 gradp(:)   = (gamma-1.)*(rhogasi*gradu - rhoi*uui*gradsumeps)
 tsi = 0.

#ifdef DUST
 call get_ts(idrag,grainsizek,graindensk,rhogasi,rhodusti,spsoundi,0.,tsi,iregime)
 if (iregime /= 0) stop 'iregime /= 0'
#endif
 ! this is equation (13) of Price & Laibe (2015) except
 ! that the sign on the second term is wrong in that paper
 ! (it is correct in Laibe & Price 2014a,b)
 dudtdust_func = -pri/rhogasi*divvfunc(xyzhi) &
                 +dustfracisum*tsi/rhogasi*dot_product(gradp,gradu)

end function dudtdust_func

real function deltavx_func(xyzhi)
 use eos,  only:gamma
 use part, only:ndusttypes
#ifdef DUST
 use dust, only:get_ts,idrag
#endif
 use part, only:rhoh
 real, intent(in) :: xyzhi(4)
 real :: rhoi,dustfraci,rhogasi,rhodusti,uui,pri,spsoundi,tsi,gradp
 real :: dustfracisum,gradsumeps,gradu
#ifdef DUST
 integer :: iregime
#endif

 rhoi       = rhoh(xyzhi(4),massoftype(1))
 dustfraci  = dustfrac_func(xyzhi)
 dustfracisum = real(ndusttypes)*dustfraci
 rhogasi    = (1.-dustfracisum)*rhoi
 rhodusti   = dustfracisum*rhoi
 gradsumeps = real(ndusttypes)*ddustevoldx(xyzhi)
 gradu      = dudx(xyzhi)
 uui        = utherm(xyzhi)
 pri        = (gamma-1.)*rhogasi*uui
 spsoundi   = gamma*pri/rhogasi
 tsi = 0.
#ifdef DUST
 call get_ts(idrag,grainsizek,graindensk,rhogasi,rhodusti,spsoundi,0.,tsi,iregime)
#endif
 gradp = (gamma-1.)*(rhogasi*gradu - rhoi*uui*gradsumeps)
 deltavx_func = tsi*gradp/rhogasi

end function deltavx_func

subroutine rcut_checkmask(rcut,xyzh,npart,checkmask)
 real,    intent(in)  :: rcut
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(in)  :: npart
 logical, intent(out) :: checkmask(:)
 real                 :: rcut2, xi,yi,zi,r2
 integer              :: i,ncheck

 ncheck = 0
 rcut2 = rcut*rcut
 checkmask(:) = .false.
 do i=1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    r2 = xi*xi + yi*yi + zi*zi
    if (r2 < rcut2) then
       checkmask(i) = .true.
       ncheck = ncheck + 1
    endif
 enddo

end subroutine rcut_checkmask

end module testderivs
