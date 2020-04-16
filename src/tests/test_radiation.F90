!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testradiation
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Sergei Biriukov
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, densityforce, deriv, dim, eos, forces, io,
!    kernel, linklist, options, part, physcon, radiation_utils,
!    readwrite_dumps, step_lf_global, testutils, unifdis, units
!+
!--------------------------------------------------------------------------
module testradiation
 use part, only:ithick,iradxi,ifluxx,ifluxy,ifluxz,idflux,ikappa
 use io,   only:id,master
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
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING RADIATION MODULE'

 call set_units(dist=au,mass=solarm,G=1.d0)
 call test_exchange_terms(ntests,npass)

#ifndef PERIODIC
 if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF RADIATION DERIVS (need -DPERIODIC)'
#else

 call test_uniform_derivs(ntests,npass)
#endif

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
 use testutils,  only:checkval
 use physcon,    only:au,solarm,seconds
 use dim,        only:maxp
 use options,    only:exchange_radiation_energy
 use io,         only:iverbose
 use part,       only:init_part,npart,rhoh,xyzh,fxyzu,vxyzu,massoftype,igas,&
                      iphase,maxphase,isetphase,rhoh,&
                      npartoftype,radiation,maxvxyzu
 use kernel,     only:hfact_default
 use unifdis,    only:set_unifdis
 use eos,        only:gmw,gamma,polyk
 use boundary,   only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 real :: psep,hfact
 real :: pmassi,rhozero,totmass
 integer, intent(inout) :: ntests,npass
 real :: dt,t,physrho,rhoi,maxt,laste
 integer :: i,ierr
 logical, parameter :: write_output = .false.

 call init_part()
 fxyzu(:,:) = 0.
 iverbose = 1
 exchange_radiation_energy = .false.

 radiation(:,:) = 0.
 radiation(ithick,:) = 1.
 psep = 1./16.
 hfact = hfact_default
 npart = 0
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.5,0.5)
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,hfact,npart,xyzh)
 rhozero = 1.e-7/unit_density  ! 1e-7 g/cm^3
 totmass = rhozero*(dxbound*dybound*dzbound)
 massoftype(igas) = totmass/npart
 gamma = 5./3.
 gmw = 2.0
 polyk = 0.
 if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(1) = npart
 pmassi = massoftype(igas)

 do i=1,npart
    rhoi         = rhoh(xyzh(4,i),pmassi)
    radiation(iradxi,i) = 1e12/(unit_ergg*unit_density)/rhoi
    radiation(ikappa,i)  = 0.4/unit_opacity
    vxyzu(4,i)   = 1e10/(unit_ergg*unit_density)
    vxyzu(4,i)   = vxyzu(4,i)/rhoi
    fxyzu(4,i)  = 0
 enddo

 maxt = 5e-7*seconds
 t = 0.
 rhoi    = rhoh(xyzh(4,1),pmassi)
 physrho = rhoi*unit_density
 i = 0
 do while(t < maxt/utime)
    dt = max(1e-18*seconds/utime,0.05*t)
    ! dt = maxt/utime
    call update_radenergy(1,xyzh,fxyzu,vxyzu,radiation,dt)
    ! call solve_internal_energy_implicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    ! call solve_internal_energy_explicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    t = t + dt
    if (mod(i,10)==0) then
       laste = (vxyzu(4,1)*unit_ergg)*physrho
       if (write_output) write(24,*) t*utime, laste,(radiation(iradxi,1)*unit_ergg)*physrho
    endif
    i = i + 1
 enddo
 call checkval(laste,21195027.055207778,1e-10,ierr,'energy exchange for gas cooling')
 ntests = ntests + 1
 if (ierr == 0) npass = npass + 1

 do i=1,npart
    rhoi         = rhoh(xyzh(4,i),pmassi)
    radiation(iradxi,i) = 1e12/(unit_ergg*unit_density)/rhoi
    radiation(ikappa,i) = 0.4/unit_opacity
    vxyzu(4,i)   = 1e2/(unit_ergg*unit_density)
    vxyzu(4,i)   = vxyzu(4,i)/rhoi
    fxyzu(4,i)  = 0
 enddo

 dt = 1e-11*seconds/utime
 t = 0.
 physrho = rhoi*unit_density
 i = 0
 do while(t < maxt/utime)
    dt = max(1e-18*seconds/utime,0.05*t)
    ! dt = maxt/utime
    call update_radenergy(1,xyzh,fxyzu,vxyzu,radiation,dt)
    t = t + dt
    if (mod(i,10)==0) then
       laste = (vxyzu(4,1)*unit_ergg)*physrho
       if (write_output) write(25,*) t*utime, laste,(radiation(iradxi,1)*unit_ergg)*physrho
    endif
    i = i + 1
 enddo
 call checkval(laste,21142367.365743987,1e-10,ierr,'energy exchange for gas heating')

 ntests = ntests + 1
 if (ierr == 0) npass = npass + 1

end subroutine test_exchange_terms

!---------------------------------------------------------
!+
!  unit tests of radiation derivatives: grad E and div F
!+
!---------------------------------------------------------
subroutine test_uniform_derivs(ntests,npass)
 use dim,             only:maxp
 use io,              only:id,master
 use part,            only:npart,xyzh,vxyzu,massoftype,igas,&
                           iphase,maxphase,isetphase,rhoh,npartoftype,&
                           radiation,ifluxx,maxvxyzu,init_part
 use kernel,          only:hfact_default
 use unifdis,         only:set_unifdis
 use units,           only:set_units,unit_opacity,get_c_code,get_steboltz_code,unit_velocity,unit_ergg
 use physcon,         only:Rg,pi,seconds
 use eos,             only:gamma,gmw
 use readwrite_dumps, only:write_fulldump
 use boundary,        only:set_boundary
 use testutils,       only:checkvalbuf,checkvalbuf_end
 use deriv,           only:get_derivs_global
 use step_lf_global,  only:init_step,step
 use timestep,        only:dtmax
 integer, intent(inout) :: ntests,npass
 real :: psep,hfact,a,c_code,cv1,rhoi,steboltz_code
 real :: dtext,pmassi, dt,t,kappa_code
 real :: xmin,xmax,ymin,ymax,zmin,zmax,Tref,xi0,D0,rho0,l0
 real :: dtnew
 real :: exact_grE,exact_DgrF,exact_xi
 real :: errmax_e,errmax_f,tol_e,tol_f,errmax_xi,tol_xi
 integer :: i,j
 integer :: nactive,nerr_e,ncheck_e,nerr_f,ncheck_f,nerr_xi,ncheck_xi
 character(len=20) :: string

 psep = 1./32.
 hfact = hfact_default
 npart = 0
 xmin = -0.5
 xmax =  0.5
 ymin = -0.1
 ymax =  0.1
 zmin = -0.1
 zmax =  0.1
 call init_part()
 call set_boundary(xmin,xmax,ymin,ymax,zmin,zmax)
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,hfact,npart,xyzh)
 massoftype(igas) = 1./npart*1e-25
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
 radiation(ithick,:) = 1.
 kappa_code = 1.0/unit_opacity
 Tref = 100.

 rho0 = rhoh(xyzh(4,1),pmassi)
 xi0 = a*Tref**4.0/rho0
 do i=1,npart
    vxyzu(4,i) = (Tref/cv1)/(unit_ergg)
    radiation(ikappa,i) = kappa_code
    radiation(iradxi,i) = xi0*(1. + 1e-1*sin(xyzh(1,i)*2.*pi/(xmax-xmin)))
    ! etot = vxyzu(4,i) + radiation(iradxi,i)
    ! Tgas = vxyzu(4,i)*unit_ergg*cv1
    ! print*, vxyzu(4,i),radiation(iradxi,i),etot
    ! Trad = (rhoi*(etot-vxyzu(4,i))/a)**(1./4.)
    ! print*, Tref, Trad, Tgas
 enddo

 dt = 1e-23 !*seconds/utime
 t  = 0
 dtmax = dt
 dtext = dt
 do i = 1,2
    call get_derivs_global()
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

    call checkvalbuf(radiation(ifluxx,i),exact_grE,tol_e, '  grad{E}',nerr_e,ncheck_e,errmax_e)
    !  this test only works with fixed lambda = 1/3
    call checkvalbuf(radiation(idflux,i),exact_DgrF,tol_f,'D*grad{F}',nerr_f,ncheck_f,errmax_f)
 enddo
 call checkvalbuf_end('  grad{E}',ncheck_e,nerr_e,errmax_e,tol_e)
 call checkvalbuf_end('D*grad{F}',ncheck_f,nerr_f,errmax_f,tol_f)

 ntests = ntests + 1
 if (nerr_e == 0) npass = npass + 1
 ntests = ntests + 1
 if (nerr_f == 0) npass = npass + 1

 call init_step(npart,t,dtmax)
 do i = 1,50
    t = t + dt
    dtext = dt
    call step(npart,nactive,t,dt,dtext,dtnew)

    if (mod(i,10) == 0) then
       nerr_xi = 0
       ncheck_xi = 0
       errmax_xi = 0.
       tol_xi = 2e-2
       do j = 1,npart
          rhoi = rhoh(xyzh(4,i),pmassi)
          D0  = c_code*(1./3)/kappa_code/rhoi
          exact_xi = xi0*(1.+0.1*sin(xyzh(1,i)*l0)*exp(-l0*l0*t*D0))
          write (string,"(a,i2.2,a)") 'xi(t_', i, ')'
          call checkvalbuf(&
             radiation(iradxi,i),exact_xi,tol_xi,&
             trim(string),nerr_xi,ncheck_xi,errmax_xi)
       enddo
       call checkvalbuf_end(trim(string),ncheck_xi,nerr_xi,errmax_xi,tol_xi)
       ntests = ntests + 1
       if (nerr_e == 0) npass = npass + 1
    endif
    ! write (filename,'(A5,I2.2)') 'rad_test_', i
    ! call write_fulldump(t,filename)
 enddo

 ! reset various things
 call init_part()

end subroutine test_uniform_derivs

end module testradiation
