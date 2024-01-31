!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testradiation
!
! Unit tests for radiation hydro
!
! :References:
!   - Whitehouse & Bate (2004), 353, 1078
!   - Whitehouse, Bate & Monaghan (2005), 364, 1367
!   - Biriukov (2019), PhD thesis, Monash Univ.
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, deriv, dim, eos, io, kernel, linklist,
!   mpidomain, mpiutils, options, part, physcon, radiation_implicit,
!   radiation_utils, readwrite_dumps, step_lf_global, testutils, timestep,
!   unifdis, units
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
 use dim,     only:do_radiation,periodic,mpi
 use options, only:implicit_radiation
 use io,      only:iverbose
 use radiation_implicit, only:tol_rad
 integer, intent(inout) :: ntests,npass

 if (.not.do_radiation) then
    if (id==master) write(*,"(/,a,/)") '--> SKIPPING RADIATION TEST (NEED RADIATION=yes)'
    return
 endif
 if (id==master) write(*,"(/,a,/)") '--> TESTING RADIATION MODULE'

 call set_units(dist=au,mass=solarm,G=1.d0)
 call test_exchange_terms(ntests,npass,use_implicit=.false.)
 if (.not.mpi) call test_exchange_terms(ntests,npass,use_implicit=.true.)

 if (.not.periodic) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF RADIATION DERIVS (need -DPERIODIC)'
 else
    if (.not.mpi) call test_implicit_matches_explicit(ntests,npass)

    iverbose = -1
    implicit_radiation = .false.
    call test_radiation_diffusion(ntests,npass)

    if (.not.mpi) then
       implicit_radiation = .true.
       tol_rad = 1.e-6
       call test_radiation_diffusion(ntests,npass)
    endif
 endif

 if (id==master) write(*,"(/,a)") '<-- RADIATION TEST COMPLETE'

end subroutine test_radiation

!----------------------------------------------------
!+
!  unit tests of gas-radiation energy exchange terms
!+
!----------------------------------------------------
subroutine test_exchange_terms(ntests,npass,use_implicit)
 use radiation_utils, only:update_radenergy,kappa_cgs
 use units,      only:set_units,unit_ergg,unit_density,unit_opacity,utime
 use physcon,    only:au,solarm,seconds
 use dim,        only:maxp,periodic
 use options,    only:exchange_radiation_energy
 use io,         only:iverbose
 use part,       only:init_part,npart,rhoh,xyzh,fxyzu,vxyzu,massoftype,igas,&
                      iphase,maxphase,isetphase,rhoh,drad,&
                      npartoftype,rad,radprop,maxvxyzu
 use kernel,     only:hfact_default
 use unifdis,    only:set_unifdis
 use eos,        only:gmw,gamma,polyk,iopacity_type
 use boundary,   only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use mpiutils,   only:reduceall_mpi
 use mpidomain,  only:i_belong
 use radiation_implicit, only:do_radiation_implicit
 use linklist,   only:set_linklist
 real :: psep,hfact
 real :: pmassi,rhozero,totmass
 integer, intent(inout) :: ntests,npass
 logical, intent(in)    :: use_implicit
 real :: dt,t,physrho,rhoi,maxt,laste
 integer :: i,nerr(1),ndiff(1),ncheck,ierrmax,ierr,itest
 integer(kind=8) :: nptot
 logical, parameter :: write_output = .true.
 character(len=12) :: string,filestr

 call init_part()
 iverbose = 0
 exchange_radiation_energy = .false.
 iopacity_type = -1  ! preserve the opacity value
 kappa_cgs = 0.4

 if (use_implicit) then
    string = ' implicit'
    iopacity_type = 2
 else
    string = ' explicit'
 endif

 if (id==master) write(*,"(/,a)") '--> checking radiation exchange terms'//trim(string)

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

 if (use_implicit) call set_linklist(npart,npart,xyzh,vxyzu)

 !
 ! first version of the test: set gas temperature high and radiation temperature low
 ! so that gas cools towards radiation temperature (itest=1)
 !
 ! second version: gas heats towards radiation temperature (itest=2)
 !
 do itest = 1,2
    do i=1,npart
       rhoi              = rhoh(xyzh(4,i),pmassi)
       rad(iradxi,i)     = 1e12/(unit_ergg*unit_density)/rhoi
       radprop(ikappa,i) = kappa_cgs/unit_opacity
       if (itest==2) then
          vxyzu(4,i)     = 1e2/(unit_ergg*unit_density)
       else
          vxyzu(4,i)     = 1e10/(unit_ergg*unit_density)
       endif
       vxyzu(4,i)        = vxyzu(4,i)/rhoi
       fxyzu(4,i)        = 0.
    enddo

    if (write_output) then
       write(filestr,"(i2.2)") itest
       open(unit=23+itest,file='radiation-test'//trim(filestr)//trim(adjustl(string))//'.ev')
       write(23+itest,"(a)") '# time [s],Egas [erg/cm^3],Erad [erg/cm^3]'
    endif
    maxt = 5e-7*seconds
    t = 0.
    rhoi    = rhoh(xyzh(4,1),pmassi)
    physrho = rhoi*unit_density
    i = 0
    ndiff = 0
    ncheck = 0
    ierrmax = 0
    do while(t < maxt/utime)
       dt = max(1d-18*seconds/utime,0.05d0*t)
       !  dt = maxt/utime
       if (use_implicit) then
          if (i > 1) dt = 0.05*maxt/utime ! use large timesteps for implicit version
          call do_radiation_implicit(dt,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
          call checkvalbuf(ierr,0,0,'no errors from implicit solver',ndiff(1),ncheck,ierrmax)
       else
          call update_radenergy(1,xyzh,fxyzu,vxyzu,rad,radprop,dt)
       endif
       t = t + dt
       if (mod(i,10)==0 .or. use_implicit) then
          laste = (vxyzu(4,1)*unit_ergg)*physrho
          if (write_output) write(23+itest,*) t*utime, laste,(rad(iradxi,1)*unit_ergg)*physrho
       endif
       i = i + 1
    enddo
    if (ncheck > 0) then
       call checkvalbuf_end('no errors from implicit solver',ncheck,ndiff(1),ierrmax,0)
       call update_test_scores(ntests,ndiff,npass)
    endif
    if (itest==2) then
       call checkval(laste,21144463.0313597,3.3e-15,nerr(1),'energy exchange for gas heating'//trim(string))
    else
       call checkval(laste,21197127.9406196,2e-15,nerr(1),'energy exchange for gas cooling'//trim(string))
    endif
    call update_test_scores(ntests,nerr,npass)

    if (write_output) close(23+itest)
 enddo

end subroutine test_exchange_terms

!---------------------------------------------------------
!+
!  unit tests to check that derivatives computed
!  in the radiation_implicit routines match those
!  computed in the regular density/force routines
!  in an explicit calculation
!+
!---------------------------------------------------------
subroutine test_implicit_matches_explicit(ntests,npass)
 use part,     only:npart,xyzh,vxyzu,rad,radprop,drad,&
                    xyzh_label,init_part,radprop_label
 use boundary, only:xmin,xmax,ymin,ymax,zmin,zmax
 use options,  only:implicit_radiation,tolh
 use physcon,  only:pi
 use deriv,    only:get_derivs_global
 use io,       only:iverbose
 use timestep, only:dtmax
 use radiation_implicit, only:do_radiation_implicit
 integer, intent(inout) :: ntests,npass
 real(kind=kind(radprop)), allocatable :: flux_explicit(:,:)
 real :: kappa_code,c_code,xi0,rho0,errmax_e,tol_e,tolh_old,pmassi !,exact(9)
 integer :: i,j,itry,nerr_e(9),ierr

 if (id==master) write(*,"(/,a)") '--> checking implicit routine matches explicit for flux terms'

 implicit_radiation = .false.
 iverbose = 0
 call init_part()
 tolh_old = tolh
 tolh = 1.e-8
 dtmax = 5.e-22
 do itry = 1,2
    if (itry == 2) implicit_radiation = .true.
    call setup_radiation_diffusion_problem_sinusoid(kappa_code,c_code,xi0,rho0,pmassi)

    ! set some non-zero velocities
    do i=1,npart
       vxyzu(1,i) = 0.1*sin((xyzh(1,i)-xmin)*2.*pi/(xmax-xmin))
       vxyzu(2,i) = 0.1*cos((xyzh(2,i)-ymin)*2.*pi/(ymax-ymin))
       vxyzu(3,i) = 0.1*sin((xyzh(3,i)-zmin)*2.*pi/(zmax-zmin))
    enddo
    call get_derivs_global()
    if (itry==1) then
       call get_derivs_global()  ! twice to get density on neighbours correct

       !--allocate and copy flux
       allocate(flux_explicit(3,npart))
       flux_explicit = radprop(ifluxx:ifluxz,1:npart)
    else
       radprop = 0.
       call do_radiation_implicit(1.e-24,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
    endif
 enddo

 ! now check that things match
 tol_e = 1.e-15
 nerr_e = 0
 errmax_e = 0.
 do j=1,3
    call checkval(npart,radprop(ifluxx+j-1,:),flux_explicit(j,:),tol_e,nerr_e(j),radprop_label(ifluxx+j-1))
 enddo
 call update_test_scores(ntests,nerr_e,npass)

 vxyzu(1:3,:) = 0.  ! zero velocities again
 tolh = tolh_old

end subroutine test_implicit_matches_explicit

!---------------------------------------------------------
!+
!  unit tests of radiation derivatives: grad E and div F
!+
!---------------------------------------------------------
subroutine test_radiation_diffusion(ntests,npass)
 use part,            only:npart,xyzh,rhoh,vxyzu,&
                           rad,radprop,drad,ifluxx,maxvxyzu,fxyzu,init_part
 use boundary,        only:xmin,xmax
 use physcon,         only:pi
 use readwrite_dumps, only:write_fulldump
 use deriv,           only:get_derivs_global
 use step_lf_global,  only:init_step,step
 use timestep,        only:dtmax
 use options,         only:implicit_radiation,limit_radiation_flux,implicit_radiation_store_drad
 use mpiutils,        only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 real :: rhoi,dtext,pmassi,dt,t,kappa_code
 real :: xi0,D0,rho0,l0,c_code
 real :: dtnew,tmax
 real :: exact_grE,exact_DgrF,exact_xi
 real :: errmax_e,errmax_f,tol_e,tol_f,errmax_xi,tol_xi,de,dekin,degas,derad
 integer :: i,j
 integer :: nactive,nerr_e(1),ncheck_e,nerr_f(1),ncheck_f,nerr_xi(1),ncheck_xi
 character(len=20) :: string,filename
 logical, parameter :: write_output = .false.

 if (id==master) write(*,"(/,a)") '--> checking radiation diffusion on sine function'//get_tag(implicit_radiation)

 call init_part()
 call setup_radiation_diffusion_problem_sinusoid(kappa_code,c_code,xi0,rho0,pmassi)

 tmax = 5.e-22
 dtmax = tmax
 implicit_radiation_store_drad = .true.
 limit_radiation_flux = .false.
 do i = 1,2
    call get_derivs_global(dt_new=dtnew,dt=1.e-25)
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
    call checkvalbuf(radprop(ifluxx,i),exact_grE,tol_e, '  grad{E}',nerr_e(1),ncheck_e,errmax_e)
    !  this test only works with fixed lambda = 1/3
    if (.not.limit_radiation_flux) then
       call checkvalbuf(drad(iradxi,i),exact_DgrF,tol_f,'dxi/dt = -div F = div(D grad E)',nerr_f(1),ncheck_f,errmax_f)
    endif
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
 de = pmassi*(dekin + degas + derad)/xi0
 de = reduceall_mpi('+',de)
 call checkval(de,0.,4.e-9,nerr_e(1),'dE/dt = 0')
 call update_test_scores(ntests,nerr_e,npass)
 !
 ! now solve diffusion as a function of time
 !
 t  = 0.
 if (write_output) then
    write (filename,'(a5,i3.3)') 'rad_test_',0
    call write_fulldump(t,filename)
 endif
 dt = dtnew
 dtext = dt
 if (implicit_radiation) dt = 6.7e-24 ! force the explicit timestep
 implicit_radiation_store_drad = .false.

 call init_step(npart,t,dtmax)
 i = 0
 D0  = c_code*(1./3)/kappa_code/rho0
 if (write_output) print "(/,a,3(es10.3,a))", ' exact solution: ',xi0,'*(1 + 0.1*sin(x*',l0,')*exp(',-l0*l0*D0,'*t))'
 do while(t < tmax)
    t = t + dt
    dtext = dt
    call step(npart,nactive,t,dt,dtext,dtnew)
    dt = dtnew
    if (implicit_radiation) dt = 6.7e-24 ! force the explicit timestep

    i = i + 1

    if (mod(i,10) == 0) then
       nerr_xi = 0
       ncheck_xi = 0
       errmax_xi = 0.
       if (implicit_radiation) then
          tol_xi = 5.e-3
       else
          tol_xi = 3.5e-4
       endif
       do j = 1,npart
          rhoi = rhoh(xyzh(4,j),pmassi)
          D0  = c_code*(1./3)/kappa_code/rhoi
          exact_xi = xi0*(1.+0.1*sin(xyzh(1,j)*l0)*exp(-l0*l0*t*D0))
          write (string,"(a,i3.3,a)") 'xi(t_', i, ')'
          call checkvalbuf(rad(iradxi,j),exact_xi,tol_xi,trim(string),&
                           nerr_xi(1),ncheck_xi,errmax_xi)
       enddo
       call checkvalbuf_end(trim(string),ncheck_xi,nerr_xi(1),errmax_xi,tol_xi)
       call update_test_scores(ntests,nerr_xi,npass)
    endif
    if (write_output) then
       write (filename,'(a5,i3.3)') 'rad_test_', i
       call write_fulldump(t,filename)
    endif
 enddo

 ! reset various things
 call init_part()
 limit_radiation_flux = .true.
 drad = 0.

end subroutine test_radiation_diffusion

!---------------------------------------------------------
!+
!  subroutine to setup sinusoidal diffusion problem
!  (split so this can be called by different subtests)
!+
!---------------------------------------------------------
subroutine setup_radiation_diffusion_problem_sinusoid(kappa_code,c_code,xi0,rho0,pmassi)
 use dim,             only:maxp
 use io,              only:id,master
 use eos,             only:gamma,gmw,iopacity_type
 use part,            only:npart,hfact,xyzh,vxyzu,massoftype,igas,periodic,&
                           iphase,maxphase,isetphase,rhoh,npartoftype,&
                           rad,radprop,ifluxx,maxvxyzu
 use units,           only:unit_opacity,get_c_code,unit_ergg,get_radconst_code
 use physcon,         only:Rg,pi,seconds
 use mpidomain,       only:i_belong
 use mpiutils,        only:reduceall_mpi
 use boundary,        only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,          only:hfact_default
 use unifdis,         only:set_unifdis
 use radiation_utils, only:kappa_cgs
 real, intent(out) :: kappa_code,c_code,xi0,rho0,pmassi
 real :: psep,a,cv1,Tref
 integer :: i,nactive
 integer(kind=8) :: nptot

 psep = 1./32.
 hfact = hfact_default
 npart = 0
 call set_boundary(-0.5,0.5,-0.1,0.1,-0.1,0.1)
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                  hfact,npart,xyzh,periodic,mask=i_belong)
 nptot = reduceall_mpi('+',npart)

 rho0 = 2.5e-24
 massoftype(igas) = rho0*dxbound*dybound*dzbound/nptot
 pmassi = massoftype(igas)
 if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(igas) = npart
 nactive = npart

 iopacity_type = 2  ! constant opacity value
 c_code = get_c_code()
 gamma = 5./3.
 gmw = 2.0
 cv1 = (gamma-1.)*gmw/Rg*unit_ergg
 a   = get_radconst_code()
 pmassi = massoftype(igas)
 kappa_cgs = 1.0
 kappa_code = kappa_cgs/unit_opacity
 Tref = 100.

 xi0 = a*Tref**4/rho0
 do i=1,npart
    vxyzu(4,i) = (Tref/cv1)
    radprop(ikappa,i) = kappa_code
    rad(iradxi,i)     = xi0*(1. + 1e-1*sin(xyzh(1,i)*2.*pi/(xmax-xmin)))
    radprop(ithick,i) = 1.
    !etot = vxyzu(4,i) + rad(iradxi,i)
    !Tgas = vxyzu(4,i)*cv1
    !rhoi = rhoh(xyzh(4,i),pmassi)
    !Trad = (rho0*(etot-vxyzu(4,i))/a)**(1./4.)
    !print*, Tref, Trad, Tgas
 enddo

end subroutine setup_radiation_diffusion_problem_sinusoid
!---------------------------------------------------------
!+
!  function to return 'implicit' or 'explicit' string
!  based on logical flag
!+
!---------------------------------------------------------
function get_tag(implicit_flag) result(str)
 logical, intent(in) :: implicit_flag
 character(len=9) :: str

 if (implicit_flag) then
    str = ' implicit'
 else
    str = ' explicit'
 endif

end function get_tag

end module testradiation
