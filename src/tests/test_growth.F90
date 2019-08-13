!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testgrowth
!
!  DESCRIPTION:
!   Unit tests of the growth module
!
!  REFERENCES:
!
!  OWNER: Arnaud Vericel
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, deriv, dim, dust, energies, eos, growth, io,
!    kernel, mpiutils, options, part, physcon, step_lf_global, testutils,
!    timestep, unifdis, units, viscosity
!+
!--------------------------------------------------------------------------
module testgrowth
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_growth

 private

contains

subroutine test_growth(ntests,npass)
#ifdef DUST
#ifdef DUSTGROWTH
 use growth,      only:init_growth,get_growth_rate,ifrag,isnow
 use physcon,     only:solarm,au
 use units,       only:set_units
 use mpiutils,    only:barrier_mpi
#endif
#endif
 integer, intent(inout) :: ntests,npass

#ifdef DUST
#ifdef DUSTGROWTH
 integer :: nfailed(5),ierr !don't forget the dimension of nfailed

 if (id==master) write(*,"(/,a)") '--> TESTING DUSTGROWTH MODULE'

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
 call dustybox_thereturn(ntests,npass)
 call barrier_mpi()

 !
 ! the infamous Laibe 2008-inspired test
 !
 call The_big_laiboxi(ntests,npass) !- Nobody calls me Laiboxi. You got the wrong guy. I'm the dude, man.
 call barrier_mpi()

 if (id==master) write(*,"(/,a)") '<-- DUSTGROWTH TEST COMPLETE'
#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING DUSTGROWTH TEST (REQUIRES -DDUST -DDUSTGROWTH)'
#endif
#endif

end subroutine test_growth

#ifdef DUST
#ifdef DUSTGROWTH

subroutine dustybox_thereturn(ntests,npass)
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:igas,idust,npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,dustevol,ddustevol,temperature,iphase,iamdust,maxtypes,&
                          xyzmh_ptmass,dustproppred,VrelVf,dustgasprop,nptmass,Omega_k
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use energies,       only:compute_energies
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use eos,            only:ieos,polyk,gamma,qfacdisc,get_spsound
 use dust,           only:K_code,idrag,init_drag
 use growth,         only:ifrag
 use options,        only:alpha,alphamax
 use unifdis,        only:set_unifdis
 use dim,            only:periodic,mhd,use_dust
 use timestep,       only:dtmax
 use io,             only:iverbose
 use mpiutils,       only:reduceall_mpi
 use physcon,        only:au,solarm,Ro,pi
 use viscosity,      only:shearparam
 use units,          only:set_units,udist,unit_velocity,unit_density

 integer,intent(inout) :: ntests,npass

 integer(kind=8) :: npartoftypetot(maxtypes)

 integer         :: nx
 integer         :: itype
 integer         :: npart_previous
 integer         :: i
 integer         :: j
 integer         :: nsteps
 integer         :: modu
 integer         :: noutputs
 integer         :: ncheck(5)
 integer         :: nerr(5)
 integer         :: ierr

 logical         :: do_output = .false.

 real            :: deltax
 real            :: dz
 real            :: hfact
 real            :: totmass
 real            :: rhozero
 real            :: errmax(5)
 real            :: dtext_dum
 real            :: Stcomp(20000)
 real            :: cscomp(20000)
 real            :: s(20000)
 real            :: r
 real            :: sinit
 real            :: dens
 real            :: t
 real            :: tmax
 real            :: dt
 real            :: dtext
 real            :: dtnew
 real            :: dv
 real            :: Vt
 real            :: vini

 real, parameter :: tolst = 5.e-4
 real, parameter :: tolcs = 5.e-4
 real, parameter :: toldv = 5.e-4
 real, parameter :: tolrho = 5.e-4

 sinit = 1.
 dens  = 1.
 vini  = 1.

 write(*,"(/,a)")'--> testing DUSTYBOX 2, THE return'

 nptmass           = 1
 xyzmh_ptmass(4,1) = 1.

 !
 ! setup for dustybox problem
 !
 nx      = 32
 deltax  = 1./nx
 dz      = 2.*sqrt(6.)/nx
 call set_boundary(-0.5,0.5,-0.25,0.25,-dz,dz)
 call set_units(mass=solarm,dist=au,G=1.d0)
 hfact   = hfact_default
 rhozero = 1.
 totmass = rhozero*dxbound*dybound*dzbound
 npart   = 0

 do itype=igas,idust,(idust-igas)
    npart_previous = npart
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                     deltax,hfact,npart,xyzh,verbose=.false.)
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
       dustprop(:,i)    = 0.
       VrelVf(i)        = 0.
       dustgasprop(:,i) = 0.
       vxyzu(:,i)       = 0.
       if (itype == idust) then
          vxyzu(1,i)       = -1.*sqrt(vini/3.)
          vxyzu(2,i)       = -1.*sqrt(vini/2.)
          vxyzu(3,i)       = sqrt(vini/6.)
       endif
       fext(:,i)           = 0.
       if (mhd) Bevol(:,i) = 0.
       if (use_dust) then
          dustevol(:,i)    = 0.
          dustfrac(:,i)    = 0.
       endif
    enddo
    npartoftype(itype)     = npart - npart_previous
    npartoftypetot(itype)  = reduceall_mpi('+',npartoftype(itype))
    massoftype(itype)      = totmass/npartoftypetot(itype)
 enddo

 !
 ! runtime parameters
 !
 ieos       = 1
 idrag      = 2
 ifrag      = -1
 gamma      = 1.
 polyk      = 1.
 alpha      = 0.
 alphamax   = 0.
 iverbose   = 0
 shearparam = 1.
 dv         = 0.
 dt         = 1.e-3
 tmax       = 0.1
 nsteps     = int(tmax/dt)
 noutputs   = 100
 if (noutputs > nsteps) noutputs = nsteps
 modu       = int(nsteps/noutputs)
 ncheck(:)  = 0
 nerr(:)    = 0
 errmax(:)  = 0.

 t                   = 0.
 K_code              = 5.

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)
 !
 ! run dustybox problem
 !
 do i=1,nsteps
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    t = t + dt
    do j=1,npart
       if (iamdust(iphase(j))) then
          Stcomp(j) = 1./(2.*K_code) * Omega_k(j)
          cscomp(j) = get_spsound(ieos,xyzh(:,j),rhozero,vxyzu(:,j))
          dv        = vini*exp(-2.*K_code*t)
          call checkvalbuf(dustgasprop(1,j)/cscomp(j),1.,tolcs,'csound',nerr(1),ncheck(1),errmax(1))
          call checkvalbuf(dustgasprop(2,j)/rhozero,1.,tolrho,'rhogas',nerr(2),ncheck(2),errmax(2))
          call checkvalbuf(dustgasprop(3,j)/Stcomp(j),1.,tolst,'St',nerr(3),ncheck(3),errmax(3))
          call checkvalbuf(dustgasprop(4,j)/dv,1.,toldv,'dv',nerr(4),ncheck(4),errmax(4))
       endif
    enddo
    !if (do_output .and. mod(i,modu)==0) !call write_file_err(i,t,xyzh,dustprop,s,St,&
!Stcomp,csound,cscomp,dv,npart,"dustybox2_")
 enddo
 call checkvalbuf_end('sound speed interpolation match exact number',ncheck(1),nerr(1),errmax(1),tolcs)
 call checkvalbuf_end('rhogas interpolation match exact number',ncheck(2),nerr(2),errmax(2),tolrho)
 call checkvalbuf_end('Stokes number interpolation match exact number',ncheck(3),nerr(3),errmax(3),tolst)
 call checkvalbuf_end('dv interpolation match exact solution',ncheck(4),nerr(4),errmax(4),toldv)

 call update_test_scores(ntests,nerr(1:4),npass)

end subroutine dustybox_thereturn

!-------------------
!-------------------
!-------------------

subroutine The_big_laiboxi(ntests,npass)
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:igas,idust,npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,dustevol,ddustevol,temperature,iphase,iamdust,maxtypes,&
                          xyzmh_ptmass,VrelVf,dustgasprop,nptmass
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use energies,       only:compute_energies
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use eos,            only:ieos,polyk,gamma,get_spsound
 use dust,           only:K_code,idrag,init_drag
 use growth,         only:ifrag,init_growth
 use options,        only:alpha,alphamax
 use unifdis,        only:set_unifdis
 use dim,            only:periodic,mhd,use_dust
 use timestep,       only:dtmax
 use io,             only:iverbose
 use mpiutils,       only:reduceall_mpi
 use physcon,        only:au,solarm,Ro,pi
 use viscosity,      only:shearparam
 use units,          only:set_units,udist,unit_velocity,unit_density

 integer,intent(inout) :: ntests,npass

 integer(kind=8) :: npartoftypetot(maxtypes)

 integer         :: nx
 integer         :: itype
 integer         :: npart_previous
 integer         :: i
 integer         :: j
 integer         :: nsteps
 integer         :: modu
 integer         :: noutputs
 integer         :: ncheck(4)
 integer         :: nerr(4)
 integer         :: ierr

 logical         :: do_output = .false.

 real            :: deltax
 real            :: dz
 real            :: hfact
 real            :: totmass
 real            :: rhozero
 real            :: errmax(4)
 real            :: dtext_dum
 real            :: Stcomp(20000),Stini(20000)
 real            :: cscomp(20000),tau(20000)
 real            :: s(20000),time
 real            :: r
 real            :: sinit
 real            :: dens
 real            :: t
 real            :: tmax
 real            :: dt
 real            :: dtext
 real            :: dtnew
 real            :: dv
 real            :: Vt,guillaume
 real            :: vini

 real, parameter :: tolst = 5.e-4
 real, parameter :: tolcs = 5.e-4
 real, parameter :: tols  = 5.e-4
 real, parameter :: tolrho = 5.e-4

 sinit = 1.e-4/udist
 dens  = 1./unit_density
 vini  = 0.

 write(*,"(/,a)")'--> testing THE BIG LAIBOXI'

 !
 ! initialise
 !
 nptmass           = 1
 xyzmh_ptmass(4,1) = 1.
 !
 ! setup for dustybox problem
 !
 nx      = 32
 deltax  = 1./nx
 dz      = 2.*sqrt(6.)/nx
 call set_boundary(-0.5,0.5,-0.25,0.25,-dz,dz)
 call set_units(mass=solarm,dist=au,G=1.d0)
 hfact   = hfact_default
 rhozero = 1.e-12/unit_density
 totmass = rhozero*dxbound*dybound*dzbound
 npart   = 0

 do itype=igas,idust,(idust-igas)
    npart_previous = npart
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
    deltax,hfact,npart,xyzh,verbose=.false.)
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
       vxyzu(:,i)          = 0.
       dustgasprop(:,i)    = 0.
       VrelVf(i)           = 0.
       if (itype == idust) then
          dustprop(1,i)    = sinit
          dustprop(2,i)    = dens
       else
          dustprop(:,i)    = 0.
       endif
       fext(:,i)           = 0.
       if (mhd) Bevol(:,i) = 0.
       if (use_dust) then
          dustevol(:,i)    = 0.
          dustfrac(:,i)    = 0.
       endif
    enddo
    npartoftype(itype)     = npart - npart_previous
    npartoftypetot(itype)  = reduceall_mpi('+',npartoftype(itype))
    massoftype(itype)      = totmass/npartoftypetot(itype)
 enddo

 !
 ! runtime parameters
 !
 ieos       = 1
 idrag      = 1
 ifrag      = 0
 gamma      = 1.
 polyk      = (1000./unit_velocity)**2
 alpha      = 0.
 alphamax   = 0.
 iverbose   = 0
 shearparam = 5.e-2
 dv         = 0.
 dt         = 1.e-3
 tmax       = 0.2
 nsteps     = int(tmax/dt)
 noutputs   = 50
 if (noutputs > nsteps) noutputs = nsteps
 modu       = int(nsteps/noutputs)
 ncheck(:)  = 0
 nerr(:)    = 0
 errmax(:)  = 0.

 t          = 0.

 call init_drag(ierr)
 call init_growth(ierr)

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do j=1,npart
    r                = sqrt(xyzh(1,j)**2+xyzh(2,j)**2+xyzh(3,j)**2)
    cscomp(j)        = get_spsound(ieos,xyzh(:,j),rhozero,vxyzu(:,j))
    Stini(j)         = sqrt(pi*gamma/8)*dens*sinit/(2*rhozero*cscomp(j))/r**1.5
    Stcomp(j)        = Stini(j)
    dustgasprop(3,j) = Stini(j)
    tau(j)           = 1/(sqrt(2**1.5*Ro*shearparam)/r**1.5)*2/sqrt(pi*gamma/8.)
    s(j)             = Stini(j)/(sqrt(pi*gamma/8)*dens/(2*rhozero*cscomp(j))/r**1.5)
 enddo
 !
 ! run dustybox problem
 !
 do i=1,nsteps
    dtext = dt
 !if (do_output .and. mod(i,modu)==0) !call write_file_err(i,t,xyzh,dustprop*udist,s*udist,St,Stcomp,&
 !csound,cscomp,dv,npart,"Laiboxi_")
    call step(npart,npart,t,dt,dtext,dtnew)
    t = t + dt
    do j=1,npart
       if (iamdust(iphase(j))) then
          r         = sqrt(xyzh(1,j)**2+xyzh(2,j)**2+xyzh(3,j)**2)
          time      = t/tau(j) + 2.*sqrt(Stini(j))*(1.+Stini(j)/3.)
          guillaume = (8.+9.*time*time+3.*time*sqrt(16.+9.*time*time))**(1./3.)
          Stcomp(j) = guillaume/2. + 2./guillaume - 2.
          s(j)      = Stcomp(j)/(sqrt(pi*gamma/8)*dens/(2*rhozero*cscomp(j))/r**1.5)
          call checkvalbuf(dustgasprop(3,j)/Stcomp(j),1.,tolst,'St',nerr(1),ncheck(1),errmax(1))
          call checkvalbuf(dustprop(1,j)/s(j),1.,tols,'size',nerr(2),ncheck(2),errmax(2))
          call checkvalbuf(dustgasprop(1,j)/cscomp(j),1.,tolcs,'csound',nerr(3),ncheck(3),errmax(3))
          call checkvalbuf(dustgasprop(2,j)/rhozero,1.,tolrho,'rhogas',nerr(4),ncheck(4),errmax(4))
       endif
    enddo
 enddo
 call checkvalbuf_end('Stokes number interpolation match exact solution',ncheck(1),nerr(1),errmax(1),tolst)
 call checkvalbuf_end('size interpolation match exact solution',ncheck(2),nerr(2),errmax(2),tolcs)
 call checkvalbuf_end('sound speed interpolation match exact number',ncheck(3),nerr(3),errmax(3),tols)
 call checkvalbuf_end('rhogas interpolation match exact number',ncheck(4),nerr(4),errmax(4),tolrho)

 call update_test_scores(ntests,nerr(1:4),npass)

end subroutine The_big_laiboxi

!---------------------------------------------------
!+
!  write an output file with x, y, z ,
!  dustprop(1) (size) and dustprop(4) (vrel/vfrag)
!+
!---------------------------------------------------
subroutine write_file(step,t,xyzh,dustprop,cs,npart,prefix)
 real, intent(in)              :: t
 real, intent(in)              :: xyzh(:,:),dustprop(:,:),cs(:)
 character(len=*), intent(in)  :: prefix
 integer, intent(in)           :: npart,step
 character(len=30)             :: filename,str
 integer                       :: i,lu

 write(str,"(i000.4)") step
 filename = prefix//trim(adjustl(str))//'.txt'
 open(newunit=lu,file=filename,status='replace')
 write(lu,*) t
 do i=1,npart
    write(lu,*) xyzh(1,i),xyzh(2,i),xyzh(3,i),dustprop(1,i),dustprop(4,i),cs(i)
 enddo
 close(lu)

end subroutine write_file

subroutine write_file_err(step,t,xyzh,dustprop,scomp,St,Stcomp,cs,cscomp,dvcomp,npart,prefix)
 use part,                     only:iamdust,iphase,iamgas
 real, intent(in)              :: t,dvcomp
 real, intent(in)              :: xyzh(:,:),dustprop(:,:),cs(:),St(:),Stcomp(:),scomp(:),cscomp(:)
 character(len=*), intent(in)  :: prefix
 integer, intent(in)           :: npart,step
 character(len=30)             :: filename,str,filenamebis
 integer                       :: i,lu,lulu

 write(str,"(i000.4)") step
 filename = prefix//'dust'//trim(adjustl(str))//'.txt'
 !filenamebis = prefix//'gas'//trim(adjustl(str))//'.txt'
 open(newunit=lu,file=filename,status='replace')
 !open(newunit=lulu,file=filenamebis,status='replace')
 write(lu,*) t
 !write(lulu,*) t
 write(lu,*) "x   y   z   h   s   scomp   St   Stcomp   cs   cscomp   dv   dvcomp"
 !write(lulu,*) "x   y   z   h"
 do i=1,npart
    if (iamdust(iphase(i))) write(lu,*) xyzh(1,i),xyzh(2,i),xyzh(3,i),0.5,dustprop(1,i),&
    scomp(i),St(i),Stcomp(i),cs(i),cscomp(i),dustprop(4,i),dvcomp
 !if (iamgas(iphase(i))) write(lulu,*) xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i)
 enddo
 close(lu)
 !close(lulu)

end subroutine write_file_err

#endif
#endif

end module testgrowth
