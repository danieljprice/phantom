!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module testradiation
!
! Unit tests for radiation hydro
!
! :References:
!    Whitehouse & Bate (2004), 353, 1078
!    Whitehouse, Bate & Monaghan (2005), 364, 1367
!    Biriukov (2019), PhD thesis, Monash Univ.
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, deriv, dim, domain, eos, io, kernel, mpiutils,
!   options, part, physcon, radiation_utils, readwrite_dumps,
!   step_lf_global, testutils, timestep, unifdis, units
!
 use part,      only:ithick,iradxi,ifluxx,ifluxy,ifluxz,ikappa
 use io,        only:id,master
 use testutils, only:checkval,update_test_scores,checkvalbuf,checkvalbuf_end
 implicit none

 public :: test_radiation
 private

contains
!-------------------------------------------------
!+
!  unit tests of radiation hydrodynamics
!+
!-------------------------------------------------
subroutine test_radiation(ntests,npass)
 use physcon, only:solarm,au
 use units,   only:set_units
 use dim,     only:do_radiation,periodic
 integer, intent(inout) :: ntests,npass

 if (.not.do_radiation) then
    if (id==master) write(*,"(/,a,/)") '--> SKIPPING RADIATION TEST (NEED RADIATION=yes)'
    return
 endif
 if (id==master) write(*,"(/,a,/)") '--> TESTING RADIATION MODULE'

 call set_units(dist=au,mass=solarm,G=1.d0)
 call test_exchange_terms(ntests,npass)

 if (.not.periodic) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF RADIATION DERIVS (need -DPERIODIC)'
 else
    call test_uniform_derivs(ntests,npass)
 endif

 if (id==master) write(*,"(/,a)") '<-- RADIATION TEST COMPLETE'

end subroutine test_radiation

!----------------------------------------------------
!+
!  unit tests of gas-radiation energy exchange terms
!+
!----------------------------------------------------
subroutine test_exchange_terms(ntests,npass)
 use radiation_utils, only:update_radenergy
 use units,      only:set_units,unit_ergg,unit_density,unit_opacity,utime
 use physcon,    only:au,solarm,seconds
 use dim,        only:maxp,periodic
 use options,    only:exchange_radiation_energy
 use io,         only:iverbose
 use part,       only:init_part,npart,rhoh,xyzh,fxyzu,vxyzu,massoftype,igas,&
                      iphase,maxphase,isetphase,rhoh,&
                      npartoftype,rad,radprop,maxvxyzu
 use kernel,     only:hfact_default
 use unifdis,    only:set_unifdis
 use eos,        only:gmw,gamma,polyk
 use boundary,   only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use mpiutils,   only:reduceall_mpi
 use domain,     only:i_belong
 real :: psep,hfact
 real :: pmassi,rhozero,totmass
 integer, intent(inout) :: ntests,npass
 real :: dt,t,physrho,rhoi,maxt,laste
 integer :: i,nerr(1)
 integer(kind=8) :: nptot
 logical, parameter :: write_output = .false.

 call init_part()
 iverbose = 0
 exchange_radiation_energy = .false.

 psep = 1./16.
 hfact = hfact_default
 npart = 0
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.5,0.5)
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  psep,hfact,npart,xyzh,periodic,mask=i_belong)
 rhozero = 1.e-7/unit_density  ! 1e-7 g/cm^3
 totmass = rhozero*(dxbound*dybound*dzbound)
 nptot = reduceall_mpi('+',npart)
 massoftype(igas) = totmass/nptot
 gamma = 5./3.
 gmw = 2.0
 polyk = 0.
 if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(1) = npart
 pmassi = massoftype(igas)

 do i=1,npart
    rhoi              = rhoh(xyzh(4,i),pmassi)
    rad(iradxi,i)     = 1e12/(unit_ergg*unit_density)/rhoi
    radprop(ikappa,i) = 0.4/unit_opacity
    vxyzu(4,i)        = 1e10/(unit_ergg*unit_density)
    vxyzu(4,i)        = vxyzu(4,i)/rhoi
    fxyzu(4,i)        = 0
 enddo

 maxt = 5e-7*seconds
 t = 0.
 rhoi    = rhoh(xyzh(4,1),pmassi)
 physrho = rhoi*unit_density
 i = 0
 do while(t < maxt/utime)
    dt = max(1d-18*seconds/utime,0.05d0*t)
    ! dt = maxt/utime
    call update_radenergy(1,xyzh,fxyzu,vxyzu,rad,radprop,dt)
    ! call solve_internal_energy_implicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    ! call solve_internal_energy_explicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    t = t + dt
    if (mod(i,10)==0) then
       laste = (vxyzu(4,1)*unit_ergg)*physrho
       if (write_output) write(24,*) t*utime, laste,(rad(iradxi,1)*unit_ergg)*physrho
    endif
    i = i + 1
 enddo
 call checkval(laste,21195027.055207778,1e-10,nerr(1),'energy exchange for gas cooling')
 call update_test_scores(ntests,nerr,npass)

 do i=1,npart
    rhoi              = rhoh(xyzh(4,i),pmassi)
    rad(iradxi,i)     = 1e12/(unit_ergg*unit_density)/rhoi
    radprop(ikappa,i) = 0.4/unit_opacity
    vxyzu(4,i)        = 1e2/(unit_ergg*unit_density)
    vxyzu(4,i)        = vxyzu(4,i)/rhoi
    fxyzu(4,i)        = 0
 enddo

 dt = 1e-11*seconds/utime
 t = 0.
 physrho = rhoi*unit_density
 i = 0
 do while(t < maxt/utime)
    dt = max(1d-18*seconds/utime,0.05d0*t)
    ! dt = maxt/utime
    call update_radenergy(1,xyzh,fxyzu,vxyzu,rad,radprop,dt)
    t = t + dt
    if (mod(i,10)==0) then
       laste = (vxyzu(4,1)*unit_ergg)*physrho
       if (write_output) write(25,*) t*utime, laste,(rad(iradxi,1)*unit_ergg)*physrho
    endif
    i = i + 1
 enddo
 call checkval(laste,21142367.365743987,1e-10,nerr(1),'energy exchange for gas heating')
 call update_test_scores(ntests,nerr,npass)

end subroutine test_exchange_terms

!---------------------------------------------------------
!+
!  unit tests of radiation derivatives: grad E and div F
!+
!---------------------------------------------------------
subroutine test_uniform_derivs(ntests,npass)
 use dim,             only:maxp
 use io,              only:id,master
 use part,            only:npart,xyzh,vxyzu,massoftype,igas,periodic,&
                           iphase,maxphase,isetphase,rhoh,npartoftype,&
                           rad,radprop,drad,ifluxx,maxvxyzu,init_part,fxyzu
 use kernel,          only:hfact_default
 use unifdis,         only:set_unifdis
 use units,           only:set_units,unit_opacity,get_c_code,get_steboltz_code,unit_velocity,unit_ergg
 use physcon,         only:Rg,pi,seconds
 use eos,             only:gamma,gmw
 use readwrite_dumps, only:write_fulldump
 use boundary,        only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax
 use deriv,           only:get_derivs_global
 use step_lf_global,  only:init_step,step
 use timestep,        only:dtmax
 use mpiutils,        only:reduceall_mpi
 use domain,          only:i_belong
 integer, intent(inout) :: ntests,npass
 real :: psep,hfact,a,c_code,cv1,rhoi,steboltz_code
 real :: dtext,pmassi, dt,t,kappa_code
 real :: Tref,xi0,D0,rho0,l0
 real :: dtnew,tmax
 real :: exact_grE,exact_DgrF,exact_xi
 real :: errmax_e,errmax_f,tol_e,tol_f,errmax_xi,tol_xi,de,dekin,degas,derad
 integer :: i,j
 integer :: nactive,nerr_e(1),ncheck_e,nerr_f(1),ncheck_f,nerr_xi(1),ncheck_xi
 integer(kind=8) :: nptot
 character(len=20) :: string !,filename

 psep = 1./32.
 hfact = hfact_default
 npart = 0
 call init_part()
 call set_boundary(-0.5,0.5,-0.1,0.1,-0.1,0.1)
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                  hfact,npart,xyzh,periodic,mask=i_belong)
 nptot = reduceall_mpi('+',npart)
 massoftype(igas) = 1./nptot*1e-25
 pmassi = massoftype(igas)
 if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(igas) = npart
 nactive = npart

 c_code = get_c_code()
 steboltz_code = get_steboltz_code()
 gamma = 5./3.
 gmw = 2.0
 cv1 = (gamma-1.)*gmw/Rg*unit_velocity**2
 a   = 4.*steboltz_code/c_code
 pmassi = massoftype(igas)
 kappa_code = 1.0/unit_opacity
 Tref = 100.

 rho0 = rhoh(xyzh(4,1),pmassi)
 xi0 = a*Tref**4.0/rho0
 do i=1,npart
    vxyzu(4,i) = (Tref/cv1)/(unit_ergg)
    radprop(ikappa,i) = kappa_code
    rad(iradxi,i)     = xi0*(1. + 1e-1*sin(xyzh(1,i)*2.*pi/(xmax-xmin)))
    radprop(ithick,i) = 1.
    ! etot = vxyzu(4,i) + radiation(iradxi,i)
    ! Tgas = vxyzu(4,i)*unit_ergg*cv1
    ! print*, vxyzu(4,i),radiation(iradxi,i),etot
    ! Trad = (rhoi*(etot-vxyzu(4,i))/a)**(1./4.)
    ! print*, Tref, Trad, Tgas
 enddo

 tmax = 5.e-22
 dtmax = tmax
 do i = 1,2
    call get_derivs_global(dt_new=dtnew)
 enddo

 nerr_e = 0
 ncheck_e = 0
 errmax_e = 0.
 tol_e = 1e-10

 nerr_f = 0
 ncheck_f = 0
 errmax_f = 0.
 tol_f = 2e-2

 l0 = 2*pi/(xmax-xmin)
 do i=1,npart
    rhoi = rhoh(xyzh(4,i),pmassi)
    D0  = c_code*(1./3)/kappa_code/rhoi
    exact_grE  =  xi0*rho0*0.1*l0   *cos(xyzh(1,i)*l0)
    exact_DgrF = -xi0*D0  *0.1*l0*l0*sin(xyzh(1,i)*l0)
    !print*,' got drad=',drad(iradxi,i), ' should be ',exact_DgrF
    call checkvalbuf(radprop(ifluxx,i),exact_grE,tol_e, '  grad{E}',nerr_e(1),ncheck_e,errmax_e)
    !  this test only works with fixed lambda = 1/3
    call checkvalbuf(drad(iradxi,i),exact_DgrF,tol_f,'D*grad{F}',nerr_f(1),ncheck_f,errmax_f)
 enddo
 call checkvalbuf_end('  grad{E}',ncheck_e,nerr_e(1),errmax_e,tol_e)
 call checkvalbuf_end('D*grad{F}',ncheck_f,nerr_f(1),errmax_f,tol_f)
 call update_test_scores(ntests,nerr_e,npass)
 call update_test_scores(ntests,nerr_f,npass)
 !
 ! check that energy is conserved (i.e. dEtot/dt = 0)
 !
 de = 0.; degas = 0.; derad = 0.; dekin = 0.
 do i=1,npart
    dekin = dekin + dot_product(vxyzu(1:3,i),fxyzu(1:3,i))  ! v.dv/dt
    degas = degas + fxyzu(4,i)       ! du/dt
    derad = derad + drad(iradxi,i)   ! dxi/dt
 enddo
 !print*,' GOT ',pmassi*dekin,pmassi*degas,pmassi*derad
 de = pmassi*(dekin + degas + derad)
 de = reduceall_mpi('+',de)
 call checkval(de,0.,2.6e-9,nerr_e(1),'dE/dt = 0')
 call update_test_scores(ntests,nerr_e,npass)
 !
 ! now solve diffusion as a function of time
 !
 t  = 0.
 dt = dtnew
 dtext = dt
 call init_step(npart,t,dtmax)
 i = 0
 do while(t < tmax)
    t = t + dt
    dtext = dt
    call step(npart,nactive,t,dt,dtext,dtnew)
    dt = dtnew
    i = i + 1

    if (mod(i,10) == 0) then
       nerr_xi = 0
       ncheck_xi = 0
       errmax_xi = 0.
       tol_xi = 3.5e-4
       do j = 1,npart
          rhoi = rhoh(xyzh(4,i),pmassi)
          D0  = c_code*(1./3)/kappa_code/rhoi
          exact_xi = xi0*(1.+0.1*sin(xyzh(1,i)*l0)*exp(-l0*l0*t*D0))
          write (string,"(a,i3.3,a)") 'xi(t_', i, ')'
          call checkvalbuf(rad(iradxi,i),exact_xi,tol_xi,trim(string),&
                           nerr_xi(1),ncheck_xi,errmax_xi)
       enddo
       call checkvalbuf_end(trim(string),ncheck_xi,nerr_xi(1),errmax_xi,tol_xi)
       call update_test_scores(ntests,nerr_xi,npass)
    endif
    !write (filename,'(A5,I3.3)') 'rad_test_', i
    !call write_fulldump(t,filename)
 enddo

 ! reset various things
 call init_part()

end subroutine test_uniform_derivs

end module testradiation
