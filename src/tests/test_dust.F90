!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testdust
!
!  DESCRIPTION:
!   Unit tests of the dust module
!
!  REFERENCES:
!   Laibe & Price (2011),  MNRAS 418, 1491
!   Laibe & Price (2012a), MNRAS 420, 2345
!   Laibe & Price (2012b), MNRAS 420, 2365
!   Price & Laibe (2015),  MNRAS 451, 5332
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, deriv, dim, dust, energies, eos, io, kernel,
!    mpiutils, options, part, physcon, step_lf_global, testutils, timestep,
!    unifdis, units
!+
!--------------------------------------------------------------------------
module testdust
 use testutils, only:checkval
 use io,        only:id,master
 implicit none
 public :: test_dust

 private

contains

subroutine test_dust(ntests,npass)
#ifdef DUST
 use dust,      only:idrag,init_drag,get_ts,grainsize,graindens,&
                     set_dustfrac,smincgs,smaxcgs,sindex
 use physcon,   only:solarm,au
 use units,     only:set_units,unit_density
 use eos,       only:gamma
 use dim,       only:ndusttypes
 use mpiutils,  only:barrier_mpi
#endif
 integer, intent(inout) :: ntests,npass
#ifdef DUST
 integer :: i,nfailed(10),ierr,iregime
 real    :: dustfraci(ndusttypes),dustfracisum,rhoi,rhogasi,spsoundi,tsi(ndusttypes)
 real    :: dust_to_gas

 if (id==master) write(*,"(/,a)") '--> TESTING DUST MODULE'

 call set_units(mass=solarm,dist=au,G=1.d0)

 if (id==master) write(*,"(/,a)") '--> testing drag initialisation'

 nfailed = 0
 ntests = ntests + 1
 gamma = 5./3.
 do idrag=1,2
    call init_drag(ierr)
    call checkval(ierr,0,0,nfailed(idrag),'drag initialisation')
 enddo
 if (all(nfailed==0)) npass = npass + 1

 idrag = 1
 rhoi = 1.e-13/unit_density
 spsoundi = 1.
 if (ndusttypes>1) then
    dust_to_gas = 0.01
    call set_dustfrac(dust_to_gas,dustfraci,smincgs,smaxcgs,sindex)
 else
    dustfraci(:) = 0.5
 endif
 dustfracisum = sum(dustfraci)
 rhogasi = rhoi*(1. - dustfracisum)
 do i = 1,ndusttypes
    call get_ts(idrag,grainsize(i),graindens,rhogasi,rhoi*dustfracisum,spsoundi,0.,tsi(i),iregime)
 enddo
 call checkval(iregime,1,0,nfailed(1),'deltav=0 gives Epstein drag')
 ntests = ntests + 1
 if (all(nfailed==0)) npass = npass + 1

 !
 ! Test transition between Epstein/Stokes drag
 !
 call test_epsteinstokes(ntests,npass)
 call barrier_mpi()

 !
 ! DUSTYBOX test
 !
 call test_dustybox(ntests,npass)
 call barrier_mpi()

 !
 ! DUSTYDIFFUSE test
 !
 call test_dustydiffuse(ntests,npass)
 call barrier_mpi()

 if (id==master) write(*,"(/,a)") '<-- DUST TEST COMPLETE'
#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING DUST TEST (REQUIRES -DDUST)'
#endif

end subroutine test_dust

#ifdef DUST
!----------------------------------------------------
!+
!  Dustybox test from Laibe & Price (2011, 2012a,b)
!  This tests the *two fluid* dust algorithm
!+
!----------------------------------------------------
subroutine test_dustybox(ntests,npass)
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:igas,idust,npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,dustevol,ddustfrac,temperature,iphase,iamdust,maxtypes
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use energies,       only:compute_energies,ekin
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use eos,            only:ieos,polyk,gamma
 use dust,           only:K_code,idrag
 use options,        only:alpha,alphamax,use_dustfrac
 use unifdis,        only:set_unifdis
 use dim,            only:periodic,mhd,use_dust,ndusttypes
 use timestep,       only:dtmax
 use io,             only:iverbose
 use mpiutils,       only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: nx, itype, npart_previous, i, j, nsteps, ncheck(5), nerr(5)
 real :: deltax, dz, hfact, totmass, rhozero, errmax(5), dtext_dum
 real :: t, dt, dtext, dtnew
 real :: vg, vd, deltav, ekin_exact, fd
 real, parameter :: tol = 1.e-4, tolvg = 1.e-4, tolfg = 3.3e-3, tolfd = 3.3e-3

 if (periodic) then
    if (use_dustfrac .and. ndusttypes>1) then
       if (id==master) write(*,"(/,a)") '--> skipping DUSTYBOX because use_dustfrac = yes AND ndusttypes > 1'
       return
    else
       if (id==master) write(*,"(/,a)") '--> testing DUSTYBOX'
    endif
 else
    if (id==master) write(*,"(/,a)") '--> skipping DUSTYBOX (need -DPERIODIC)'
    return
 endif
 !
 ! setup for dustybox problem
 !
 nx = 32
 deltax = 1./nx
 dz = 2.*sqrt(6.)/nx
 call set_boundary(-0.5,0.5,-0.25,0.25,-dz,dz)
 hfact = hfact_default
 rhozero = 1.
 totmass = rhozero*dxbound*dybound*dzbound
 npart = 0

 do itype=1,2
    npart_previous = npart
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                     deltax,hfact,npart,xyzh,verbose=.false.)
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
       vxyzu(:,i) = 0.
       if (itype==idust) then
          vxyzu(1,i) = 1.
       endif
       fext(:,i) = 0.
       if (mhd) Bevol(:,i) = 0.
       if (use_dust) then
          dustevol(:,i) = 0.
          dustfrac(:,i) = 0.
       endif
    enddo
    npartoftype(itype) = npart - npart_previous
    npartoftypetot(itype) = reduceall_mpi('+',npartoftype(itype))
    massoftype(itype) = totmass/npartoftypetot(itype)
 enddo
 !
 ! runtime parameters
 !
 K_code = 0.35
 ieos = 1
 idrag = 2
 polyk = 1.
 gamma = 1.
 alpha = 0.
 alphamax = 0.
 iverbose = 0
 !
 ! call derivs the first time around
 !
 dt = 1.e-3
 nsteps = 100
 t = 0
 dtmax = nsteps*dt
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)
 !
 ! run dustybox problem
 !
 ncheck(:) = 0
 nerr(:) = 0
 errmax(:) = 0.
 call init_step(npart,t,dtmax)
 do i=1,nsteps
    t = t + dt
    dtext = dt
    !print*,'t = ',t
    call step(npart,npart,t,dt,dtext,dtnew)
    call compute_energies(t)

    deltav = exp(-2.*K_code*t)
    vg = 0.5*(1. - deltav)
    vd = 0.5*(1. + deltav)
    fd = K_code*(vg - vd)
    do j=1,npart
       if (iamdust(iphase(j))) then
          call checkvalbuf(vxyzu(1,j),vd,tol,'vd',nerr(1),ncheck(1),errmax(1))
          call checkvalbuf(fxyzu(1,j),fd,tolfd,'fd',nerr(2),ncheck(2),errmax(2))
       else
          call checkvalbuf(vxyzu(1,j),vg,tolvg,'vg',nerr(3),ncheck(3),errmax(3))
          call checkvalbuf(fxyzu(1,j),-fd,tolfg,'fg',nerr(4),ncheck(4),errmax(4))
       endif
    enddo
    !call checkval(npart/2-1,vxyzu(1,1:npart),vg,tolvg,nerr(2),'vg')
    ekin_exact = 0.5*totmass*(vd**2 + vg**2)
    !print*,' step ',i,'t = ',t,' ekin should be ',ekin_exact, ' got ',ekin,(ekin-ekin_exact)/ekin_exact
    call checkvalbuf(ekin,ekin_exact,tol,'ekin',nerr(5),ncheck(5),errmax(5))
 enddo
 call checkvalbuf_end('dust velocities match exact solution',ncheck(1),nerr(1),errmax(1),tol)
 call checkvalbuf_end('dust accel matches exact solution',   ncheck(2),nerr(2),errmax(2),tolfd)
 call checkvalbuf_end('gas velocities match exact solution',ncheck(3),nerr(3),errmax(3),tolvg)
 call checkvalbuf_end('gas accel matches exact solution',   ncheck(4),nerr(4),errmax(4),tolfg)
 call checkvalbuf_end('kinetic energy decay matches exact',ncheck(5),nerr(5),errmax(5),tol)

 ntests = ntests + 1
 if (all(nerr(1:5)==0)) npass = npass + 1

end subroutine test_dustybox

!----------------------------------------------------
!+
!  3D dust diffusion test from Price & Laibe (2015)
!+
!----------------------------------------------------
subroutine test_dustydiffuse(ntests,npass)
 use dim,       only:maxp,periodic,maxtypes,mhd,ndusttypes
 use part,      only:hfact,npart,npartoftype,massoftype,igas,dustfrac,ddustfrac,dustevol, &
                     xyzh,vxyzu,Bevol,dBevol,divcurlv,divcurlB,fext,fxyzu,set_particle_type,rhoh,temperature,&
                     dustprop,ddustprop
 use options,   only:use_dustfrac
 use kernel,    only:hfact_default
 use eos,       only:gamma,polyk,ieos
 use dust,      only:K_code,idrag
 use boundary,  only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax,set_boundary
 use io,        only:iverbose
 use unifdis,   only:set_unifdis
 use deriv,     only:derivs
 use testutils, only:checkvalbuf,checkvalbuf_end
 use mpiutils,  only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: nx,j,i,n,nsteps
 integer :: nerr(1),ncheck(1)
 integer :: eps_type
 real    :: errmax(1)
 real    :: deltax,rhozero,totmass,dt,dtnew,time,tmax
 real    :: epstot,epsi(ndusttypes),rc,rc2,r2,A,B,eta
 real    :: erri,exact,errl2,term,tol
 real,allocatable   :: ddustfrac_prev(:,:)
 logical, parameter :: do_output = .true.
 real,    parameter :: t_write(5) = (/0.1,0.3,1.0,3.0,10.0/)

 if (use_dustfrac .and. periodic) then
    if (id==master) write(*,"(/,a)") '--> testing DUSTYDIFFUSE'
 else
    if (id==master) write(*,"(/,a)") '--> skipping DUSTYDIFFUSE (need -DDUST and -DPERIODIC)'
    return
 endif
 !
 ! setup uniform box
 !
 nx = 32
 deltax = 1./nx
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.5,0.5)
 hfact = hfact_default
 rhozero = 3.
 totmass = rhozero*dxbound*dybound*dzbound
 time  = 0.
 npart = 0
 npartoftype(:) = 0
 iverbose = 2
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,verbose=.false.)
 npartoftype(igas) = npart
 npartoftypetot(igas) = reduceall_mpi('+',npartoftype(igas))
 massoftype(igas)  = totmass/npartoftypetot(igas)
 allocate(ddustfrac_prev(ndusttypes,npart))
 vxyzu = 0.
 if (mhd) Bevol = 0.
 !
 ! runtime options
 !
 K_code = 0.1
 ieos = 1
 idrag = 3
 polyk = 1.
 gamma = 1.
 iverbose = 0

 !
 ! setup dust fraction in the box
 !
 epstot = 0.1

 eps_type = 2
 select case(eps_type)
 case(1)
    !--Equal dust fractions
    epsi(:) = epstot/real(ndusttypes)
 case(2)
    !--Unequal dust fractions
    do i=1,ndusttypes
       epsi(i) = 1./real(i)
    enddo
    epsi = epstot/sum(epsi)*epsi
 case default
    stop 'eps_type not valid!'
 end select

 !--check that individual dust fractions add up to the total dust fraction
 if (abs(sum(epsi)-epstot)/epstot>1.e-14) then
    write(*,"(/,a)") 'ERROR! SUM(epsilon_k) /= epsilon'
    print*,'SUM(epsilon_k) = ',sum(epsi)
    print*,'       epsilon = ',epstot
 endif

 rc   = 0.25
 rc2  = rc**2
 do i=1,npart
    r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
    if (r2 < rc2) then
       dustfrac(:,i) = epsi(:)*(1. - r2/rc2)
    else
       dustfrac(:,i) = 0.
    endif
    call set_particle_type(i,igas)
 enddo

 ! factors in exact solution (Eq. 51 in PL15)
 B = rc2/epstot
 A = epstot*B**0.6
 eta = 0.1

 !
 ! evolve only the dust fraction with a simple predictor-corrector scheme
 !
 dt = 0.05
 tmax = 10.
 nsteps = nint(tmax/dt)
 dt = tmax/nsteps
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
             dustfrac,ddustfrac,temperature,time,dt,dtnew)

 if (do_output) call write_file(time,xyzh,dustfrac,npart)
 do i=1,npart
!------------------------------------------------
!--sqrt(rho*epsilon) method
!    dustevol(:,i) = sqrt(dustfrac(:,i)*rhoh(xyzh(4,i),massoftype(igas)))
!------------------------------------------------
!--asin(sqrt(epsilon)) method
    dustevol(:,i) = asin(sqrt(dustfrac(:,i)))
!------------------------------------------------
 enddo

 nerr = 0
 ncheck = 0
 errmax = 0.
 do j=1,nsteps
    time = j*dt
    !$omp parallel do private(i)
    do i=1,npart
       ddustfrac_prev(:,i) = ddustfrac(:,i)
       dustevol(:,i) = dustevol(:,i) + dt*ddustfrac(:,i)
!------------------------------------------------
!--sqrt(rho*epsilon) method
!       dustfrac(:,i) = dustevol(:,i)**2/rhoh(xyzh(4,i),massoftype(igas))
!------------------------------------------------
!--asin(sqrt(epsilon)) method
       dustfrac(:,i) = sin(dustevol(:,i))**2
!------------------------------------------------
    enddo
    !$omp end parallel do
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
                dustprop,ddustprop,dustfrac,ddustfrac,temperature,time,dt,dtnew)
    !$omp parallel do private(i)
    do i=1,npart
       dustevol(:,i) = dustevol(:,i) + 0.5*dt*(ddustfrac(:,i) - ddustfrac_prev(:,i))
    enddo
    !$omp end parallel do

    !
    ! check solution matches the analytic solution at each timestep
    !
    term = 10.*eta*time + B
    n = 0
    errl2 = 0.
    !$omp parallel do private(i,r2,exact,erri) reduction(+:errl2,n)
    do i=1,npart
       r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
       exact = A*abs(term)**(-0.6) - r2/term
       if (exact > 0.) then
          erri  = sum(dustfrac(:,i)) - exact
          errl2 = errl2 + erri*erri
          n = n + 1
       endif
    enddo
    !$omp end parallel do
    errl2 = sqrt(errl2/n)
    tol = 2.5e-3 !1.5e-3/(1. + time)  ! take tolerance down with time
    call checkvalbuf(errl2,0.,tol,'L2 err',nerr(1),ncheck(1),errmax(1))
    !
    ! write solution to file if necessary
    !
    if (do_output .and. any(abs(t_write-time) < 0.01*dt)) call write_file(time,xyzh,dustfrac,npart)
 enddo
 call checkvalbuf_end('dust diffusion matches exact solution',ncheck(1),nerr(1),errmax(1),tol)

 !
 ! clean up dog poo
 !
 dustevol  = 0.
 dustfrac  = 0.
 ddustfrac = 0.

end subroutine test_dustydiffuse

!---------------------------------------------------------
!+
!  check that the Epstein/Stokes transition is continuous
!+
!---------------------------------------------------------
subroutine test_epsteinstokes(ntests,npass)
 use dust,      only:idrag,init_drag,get_ts,grainsize,graindens,grainsizecgs
 use units,     only:unit_density,unit_velocity,utime
 use physcon,   only:years,kb_on_mh,pi
 use testutils, only:checkval,checkvalbuf,checkvalbuf_end
 integer, intent(inout) :: ntests,npass
 integer :: iregime,ierr,i,j,nfailed,ncheck
 integer, parameter :: npts=1001, nrhopts = 11
 real :: rhogas,spsoundi,tsi,ts1,deltav,tol
 real :: smin,smax,ds,rhomin,rhomax,drho,psi,exact,err,errmax
 logical :: write_output = .false.
 character(len=60) :: filename
 integer, parameter :: lu = 36

 if (id==master) write(*,"(/,a)") '--> testing Epstein/Stokes drag transition'

 idrag = 1
 spsoundi = 6.e4/unit_velocity
 deltav = 1.e-2*spsoundi
 !print*,' T = ',(6.e4)**2/kb_on_mh*2.,' delta v/cs = ',deltav/spsoundi
 smin = 1.e-6 ! cgs units
 smax = 1.e8
 rhomin = 1.e-19  ! cgs units
 rhomax = 1.e-09
 ds = (log10(smax) - log10(smin))/real(npts-1)
 drho = (log10(rhomax) - log10(rhomin))/real(nrhopts-1)

 do j=1,nrhopts
    rhogas = rhomin*10**((j-1)*drho)/unit_density
    if (write_output) then
       write(filename,"(a,1pe8.2,a)") 'ts-rho',rhogas*unit_density,'.out'
       open(unit=lu,file=filename,status='replace')
       print "(a)",' writing '//trim(filename)
    endif
    write(filename,"(a,1pe8.2)") 'rho=',rhogas*unit_density
    ncheck = 0
    nfailed = 0
    tol = 6.3e-2
    do i=1,npts
       grainsizecgs  = smin*10**((i-1)*ds)
       call init_drag(ierr)
       !--no need to test drag transition 'ndusttypes' times...once is enough
       call get_ts(idrag,grainsize(1),graindens,rhogas,0.,spsoundi,deltav**2,tsi,iregime)
       !print*,'s = ',grainsizecgs,' ts = ',tsi*utime/years,',yr ',iregime

       if (i > 1) call checkvalbuf((tsi-ts1)/abs(tsi),0.,tol,'ts is continuous into Stokes regime',nfailed,ncheck,errmax)
       ts1 = tsi
       if (write_output) write(lu,*) grainsizecgs,tsi*utime/years,iregime
    enddo
    if (write_output) close(lu)
    call checkvalbuf_end('ts is continuous into Stokes regime: '//trim(filename),ncheck,nfailed,errmax,tol)
    ntests = ntests + 1
    if (nfailed==0) npass = npass + 1
 enddo

 !
 ! graph of variation with delta v / cs
 !
 grainsizecgs = 0.1
 rhogas = 1.e-13/unit_density
 call init_drag(ierr)
 smin = 0.
 smax = 10.
 ds = (smax-smin)/real(npts-1)
 err = 0.
 if (write_output) open(unit=lu,file='ts-deltav.out',status='replace')
 do i=1,npts
    deltav = (smin + (i-1)*ds)*spsoundi
    call get_ts(idrag,grainsize(1),graindens,rhogas,0.,spsoundi,deltav**2,tsi,iregime)
    psi = sqrt(0.5)*deltav/spsoundi
    if (i==1) then
       ts1 = tsi
    else
       ! exact non-linear Epstein drag formula
       exact = 20./3./sqrt(pi)/(deltav/spsoundi* &
                ((1./psi + 1./(2.*psi**3))*exp(-psi**2) &
               + (1. + 1./psi**2 - 1./(4.*psi**4))*sqrt(pi)*erf(psi)))
       err = err + (tsi/ts1 - exact)**2
    endif
    if (write_output) write(lu,*) deltav/spsoundi,tsi/ts1,iregime
 enddo
 err = sqrt(err/npts)
 ntests = ntests + 1
 call checkval(err,0.,3.9e-3,nfailed,'Epstein drag formula matches non-linear solution')
 if (nfailed==0) npass = npass + 1

 if (write_output) close(lu)

end subroutine test_epsteinstokes

!---------------------------------------------------
!+
!  write an output file with r, dustfrac in table
!  this is what is shown in Phantom paper
!+
!---------------------------------------------------
subroutine write_file(time,xyzh,dustfrac,npart)
 use dim, only:ndusttypes,maxp
 real, intent(in)     :: time
 real, intent(in)    :: xyzh(:,:),dustfrac(:,:)
 integer, intent(in) :: npart
 character(len=30)   :: filename,str1,str2,fmt1
 integer :: i,lu
 real    :: r2,dustfracsum(maxp)

 write(str1,"(f5.1)") time
 if (ndusttypes>1) then
    write(str2,"(I5)") ndusttypes+2
    dustfracsum = sum(dustfrac,1)
 else
    write(str2,"(I5)") ndusttypes+1
 endif
 filename = 'dustfrac_t'//trim(adjustl(str1))//'.txt'
 fmt1 = '('//trim(adjustl(str2))//'F15.8)'
 open(newunit=lu,file=filename,status='replace')
 print*,' writing '//filename
 write(lu,*) time
 do i=1,npart
    r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
    if (ndusttypes>1) then
       write(lu,fmt1) sqrt(r2),dustfracsum(i),dustfrac(:,i)
    else
       write(lu,fmt1) sqrt(r2),dustfrac(:,i)
    endif
 enddo
 close(lu)

end subroutine write_file

#endif

end module testdust
