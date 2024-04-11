!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testgrowth
!
! Unit tests of the growth module
!
! :References:
!
! :Owner: Arnaud Vericel
!
! :Runtime parameters: None
!
! :Dependencies: boundary, checksetup, deriv, dim, dust, energies, eos,
!   growth, io, kernel, mpidomain, mpiutils, options, part, physcon,
!   step_lf_global, testdust, testutils, timestep, unifdis, units,
!   viscosity
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 use testdust,  only:test_dustybox
 implicit none
 public :: test_growth

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests for dust growth using Stepinksi & Valageas method
!+
!-----------------------------------------------------------------------
subroutine test_growth(ntests,npass)
 use dim,      only:use_dust,use_dustgrowth
 use growth,   only:init_growth,ifrag,isnow
 use physcon,  only:solarm,au
 use units,    only:set_units
 use mpiutils, only:barrier_mpi
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(5),ierr !don't forget the dimension of nfailed
 logical, dimension(2)  :: logic = (/.false., .true./)
 integer                :: i,j

 if (use_dust .and. use_dustgrowth) then
    if (id==master) write(*,"(/,a)") '--> TESTING DUSTGROWTH MODULE'
 else
    if (id==master) write(*,"(/,a)") '--> SKIPPING DUSTGROWTH TEST (REQUIRES -DDUST -DDUSTGROWTH)'
    return
 endif

 call set_units(mass=solarm,dist=au,G=1.d0)

 if (id==master) write(*,"(/,a)") '--> testing growth initialisation'

 nfailed = 0
 do ifrag=0,2
    do isnow=0,2
       call init_growth(ierr)
       call checkval(ierr,0,0,nfailed(ifrag+isnow+1),'growth initialisation')
    enddo
 enddo
 call update_test_scores(ntests,nfailed,npass)

 !
 ! The return of the dustybox test
 !
 call test_dustybox(ntests,npass)
 call barrier_mpi()

 !
 ! testing farmingbox with several config.
 !
 do i=1,2
    do j=1,2
       call test_farmingbox(ntests,npass,frag=logic(i),onefluid=logic(j))
       call barrier_mpi()
    enddo
 enddo

 if (id==master) write(*,"(/,a)") '<-- DUSTGROWTH TEST COMPLETE'

end subroutine test_growth

!-------------------
!-------------------
!-------------------

subroutine test_farmingbox(ntests,npass,frag,onefluid)
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:init_part,igas,idust,npart,xyzh,vxyzu,npartoftype,&
                          massoftype,set_particle_type,&
                          fxyzu,fext,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,dustevol,ddustevol,iphase,maxtypes,&
                          VrelVf,dustgasprop,Omega_k,alphaind,iamtype,&
                          ndustlarge,ndustsmall,rhoh,deltav,this_is_a_test,periodic, &
                          npartoftypetot,update_npartoftypetot
 use step_lf_global, only:step,init_step
 use deriv,          only:get_derivs_global
 use energies,       only:compute_energies
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use eos,            only:ieos,polyk,gamma,get_spsound
 use dust,           only:idrag,init_drag
 use growth,         only:ifrag,init_growth,isnow,vfrag,gsizemincgs,get_size
 use options,        only:alpha,alphamax,use_dustfrac
 use unifdis,        only:set_unifdis
 use dim,            only:periodic,mhd,use_dust,maxp,maxalpha
 use timestep,       only:dtmax
 use io,             only:iverbose
 use mpiutils,       only:reduceall_mpi
 use physcon,        only:au,solarm,Ro,pi,fourpi
 use viscosity,      only:shearparam
 use units,          only:set_units,udist,unit_density!,unit_velocity
 use mpidomain,      only:i_belong
 use checksetup,     only:check_setup

 integer, intent(inout) :: ntests,npass
 logical, intent(in)    :: frag,onefluid

 integer :: nx,nerror,nwarn
 integer :: itype,npart_previous,i,j,nsteps,modu,noutputs
 integer :: ncheck(4),nerr(4)
 real    :: errmax(4)
 integer :: ierr,iam
 integer, parameter :: ngrid = 20000

 logical :: do_output = .false.
 real    :: deltax,dz,hfact,totmass,rhozero
 real    :: Stcomp(ngrid),Stini(ngrid)
 real    :: cscomp(ngrid),tau(ngrid)
 real    :: s(ngrid),time,timelim(ngrid)
 real    :: sinit,dens,t,tmax,dt,dtext,dtnew,guillaume,dtgratio,rhog,rhod

 real, parameter :: tolst = 5.e-4
 real, parameter :: tolcs = 5.e-4
 real, parameter :: tols  = 5.e-4
 real, parameter :: tolrho = 5.e-4

 character(len=15) :: stringfrag
 character(len=15) :: stringmethod

 ! initialise particle arrays to zero
 call init_part()

 if (frag) then
    sinit       = 1./udist
    gsizemincgs = 1.e-3
    dtgratio    = 0.5
    stringfrag  = "fragmentation"
 else
    sinit       = 3.e-2/udist
    dtgratio    = 1.
    stringfrag  = "growth"
 endif

 if (onefluid) then
    use_dustfrac = .true.
    stringmethod = "one fluid"
    ndustsmall   = 1
    ndustlarge   = 0
    dtgratio     = 1.e-1
 else
    use_dustfrac = .false.
    stringmethod = "two fluid"
    ndustsmall   = 0
    ndustlarge   = 1
 endif
 dens  = 1./unit_density

 write(*,"(/,a)") '--> testing FARMINGBOX using: '//trim(stringfrag)//&
                ' and '//trim(stringmethod)//' dust method'
 !
 ! initialise
 !
 this_is_a_test = .true.

 !
 ! setup for dustybox problem
 !
 nx      = 16
 deltax  = 1./nx
 dz      = 2.*sqrt(6.)/nx
 call set_boundary(-0.5,0.5,-0.25,0.25,-dz,dz)
 hfact   = hfact_default
 rhozero = 1.e-11/unit_density
 totmass = rhozero*dxbound*dybound*dzbound
 if (onefluid) then
    rhog = rhozero * (1-dtgratio)
    rhod = dtgratio * rhozero
 else
    rhog = rhozero
    rhod = dtgratio * rhozero
 endif
 npart     = 0
 fxyzu     = 0.
 dustprop  = 0.
 ddustprop = 0.
 ddustevol = 0.
 dBevol    = 0.
 if (maxalpha==maxp) alphaind(:,:) = 0.

 !- setting gas particles
 itype = igas
 npart_previous = npart
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,periodic,verbose=.false.,mask=i_belong)
 do i=npart_previous+1,npart
    vxyzu(:,i)       = 0.
    fext(:,i)        = 0.
    if (mhd) Bevol(:,i) = 0.
    if (use_dust) then
       dustevol(:,i) = 0.
       dustfrac(:,i) = 0.
       deltav(:,:,i) = 0.
       dustgasprop(:,i) = 0.
       VrelVf(i)        = 0.
       if (use_dustfrac) then
          dustfrac(1,i) = dtgratio
          dustprop(1,i) = fourpi/3.*dens*sinit**3
          dustprop(2,i) = dens
       else
          dustprop(:,i) = 0.
          dustfrac(:,i) = 0.
       endif
    endif
    call set_particle_type(i,itype)
 enddo
 npartoftype(itype)    = npart - npart_previous
 call update_npartoftypetot
 massoftype(itype)     = totmass/npartoftypetot(itype)

 !- setting dust particles if not one fluid
 if (.not. use_dustfrac) then
    itype = idust
    npart_previous = npart
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                     deltax,hfact,npart,xyzh,periodic,verbose=.false.,mask=i_belong)
    do i=npart_previous+1,npart
       vxyzu(:,i)       = 0.
       fext(:,i)        = 0.
       if (mhd) Bevol(:,i) = 0.
       if (use_dust) then
          dustevol(:,i) = 0.
          dustfrac(:,i) = 0.
          dustprop(1,i) = fourpi/3.*dens*sinit**3
          dustprop(2,i) = dens
          dustgasprop(:,i) = 0.
          VrelVf(i)        = 0.
       endif
       call set_particle_type(i,itype)
    enddo
    npartoftype(itype)    = npart - npart_previous
    npartoftypetot(itype) = reduceall_mpi('+',npartoftype(itype))
    massoftype(itype)     = dtgratio*totmass/npartoftypetot(itype)
 endif

 !
 ! check that particle setup is sensible
 !
 call check_setup(nerror,nwarn)

 !
 ! runtime parameters
 !

 ieos          = 1
 idrag         = 1
 if (frag) then
    ifrag      = 1
    shearparam = 2.5e-2
 else
    ifrag      = 0
    shearparam = 1.e-2
 endif
 isnow        = 0
 vfrag        = 1.e-11
 gsizemincgs  = 1.e-2
 polyk        = 1.e-3
 gamma        = 1.
 alpha        = 0.
 alphamax     = 0.
 iverbose     = 0

 !- timestepping
 dt         = 1.e-3
 tmax       = 0.2
 nsteps     = int(tmax/dt)
 noutputs   = 150
 if (noutputs > nsteps) noutputs = nsteps
 modu       = int(nsteps/noutputs)
 dtmax      = nsteps*dt

 timelim(:) = 1.e3
 ncheck(:)  = 0
 nerr(:)    = 0
 errmax(:)  = 0.

 t          = 0.

 call init_drag(ierr)
 call init_growth(ierr)

 call get_derivs_global()

 call init_step(npart,t,dtmax)

 do j=1,npart
    iam = iamtype(iphase(j))
    if (iam == idust .or. (use_dustfrac .and. iam == igas)) then
       cscomp(j)        = get_spsound(ieos,xyzh(:,j),rhog,vxyzu(:,j))
       Stini(j)         = sqrt(pi*gamma/8)*dens*sinit/((rhog+rhod)*cscomp(j)) * Omega_k(j)
       Stcomp(j)        = Stini(j)
       tau(j)           = 1/(sqrt(2**1.5*Ro*shearparam)*Omega_k(j))*(rhog+rhod)/rhod/sqrt(pi*gamma/8.)
       s(j)             = sinit
       timelim(j)       = 2*sqrt(Stini(j))*(1.+Stini(j)/3.)*tau(j)
    endif
 enddo
 if (frag) write(*,"(a,f5.1,a)") "Analytical solution no longer valid after t = ", minval(timelim), " (size < 0)"

 !
 ! run farmingbox problem
 !
 do i=1,nsteps
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    t = t + dt
    if (do_output .and. mod(i,modu)==0) then
       call write_file_err(i,t,xyzh,dustprop(1,:)*udist,s*udist,&
                           dustgasprop(3,:),Stcomp,npart,"farmingbox_")
    endif
    do j=1,npart
       iam = iamtype(iphase(j))
       if (iam == idust .or. (iam == igas .and. use_dustfrac)) then
          if (frag) then
             time   = - t/tau(j) + 2.*sqrt(Stini(j))*(1.+Stini(j)/3.)
          else
             time   = 2.*sqrt(Stini(j))*(1.+Stini(j)/3.) + t/tau(j)
          endif
          guillaume = (8.+9.*time*time+3.*time*sqrt(16.+9.*time*time))**(1./3.)
          Stcomp(j) = guillaume/2. + 2./guillaume - 2
          s(j)      = Stcomp(j)/(sqrt(pi*gamma/8)*dens/((rhog+rhod)*cscomp(j))*Omega_k(j))
          if (onefluid) then
             call checkvalbuf(dustgasprop(3,j)/Stcomp(j),1.,tolst,'St',nerr(1),ncheck(1),errmax(1))
             call checkvalbuf(get_size(dustprop(1,j),dustprop(2,j))/s(j),1.,tols,'size',nerr(2),ncheck(2),errmax(2))
          else
             call checkvalbuf(dustgasprop(3,j)/Stcomp(j),1.,tolst,'St',nerr(1),ncheck(1),errmax(1))
             call checkvalbuf(get_size(dustprop(1,j),dustprop(2,j))/s(j),1.,tols,'size',nerr(2),ncheck(2),errmax(2))
             call checkvalbuf(dustgasprop(1,j)/cscomp(j),1.,tolcs,'csound',nerr(3),ncheck(3),errmax(3))
             call checkvalbuf(dustgasprop(2,j)/rhozero,1.,tolrho,'rhogas',nerr(4),ncheck(4),errmax(4))
          endif
       endif
    enddo
 enddo
 if (onefluid) then
    call checkvalbuf_end('Stokes number evaluation matches exact solution',ncheck(1),nerr(1),errmax(1),tolst)
    call checkvalbuf_end('size evaluation matches exact solution',ncheck(2),nerr(2),errmax(2),tols)
 else
    call checkvalbuf_end('Stokes number interpolation matches exact solution',ncheck(1),nerr(1),errmax(1),tolst)
    call checkvalbuf_end('size evaluation matches exact solution',ncheck(2),nerr(2),errmax(2),tols)
    call checkvalbuf_end('sound speed interpolation matches exact number',ncheck(3),nerr(3),errmax(3),tolcs)
    call checkvalbuf_end('rhogas interpolation matches exact number',ncheck(4),nerr(4),errmax(4),tolrho)
 endif

 call update_test_scores(ntests,nerr,npass)

end subroutine test_farmingbox

subroutine write_file_err(step,t,xyzh,gsize,gsize_exact,St,St_exact,npart,prefix)
 use part,                     only:iamdust,iphase,iamgas
 real, intent(in)              :: t
 real, intent(in)              :: xyzh(:,:)
 real, intent(in)              :: St(:),St_exact(:)
 real(kind=8), intent(in)      :: gsize(:),gsize_exact(:)
 character(len=*), intent(in)  :: prefix
 integer, intent(in)           :: npart,step
 character(len=30)             :: filename,str
 integer                       :: i,lu

 write(str,"(i000.4)") step
 filename = prefix//'dust_'//trim(adjustl(str))//'.txt'
 open(newunit=lu,file=filename,status='replace')
 write(lu,*) t
 do i=1,npart
    if (iamdust(iphase(i))) write(lu,*) xyzh(1,i),xyzh(2,i),xyzh(3,i),gsize(i),gsize_exact(i),&
        St(i),St_exact(i)
 enddo
 close(lu)

end subroutine write_file_err

end module testgrowth
