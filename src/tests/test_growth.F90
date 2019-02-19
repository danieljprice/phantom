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
 use testutils, only:checkval
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
 ntests  = ntests + 1
 do ifrag=0,2
    do isnow=0,2
       call init_growth(ierr)
       call checkval(ierr,0,0,nfailed(ifrag+isnow+1),'growth initialisation')
    enddo
 enddo
 if (all(nfailed==0)) npass = npass + 1

 !
 ! GROWINGBOX test
 !
 call test_growingbox(ntests,npass)
 call barrier_mpi()
 !
 ! check stokes number interpolation
 !
 call check_stokes_number(ntests,npass)
 call barrier_mpi()

 if (id==master) write(*,"(/,a)") '<-- DUSTGROWTH TEST COMPLETE'
#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING DUSTGROWTH TEST (REQUIRES -DDUST -DDUSTGROWTH)'
#endif
#endif

end subroutine test_growth

#ifdef DUST
#ifdef DUSTGROWTH
!----------------------------------------------------
!+
!  Growingbox test
!  This tests the dust growth algorithm
!+
!----------------------------------------------------
subroutine test_growingbox(ntests,npass)
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:idust,igas,npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,ddustevol,temperature,iamdust,maxtypes,St,alphaind,csound
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use unifdis,        only:set_unifdis
 use eos,            only:ieos,polyk,gamma,get_spsound,init_eos,temperature_coef,gmw
 use viscosity,      only:shearparam
 use physcon,        only:au,solarm,Ro
 use dim,            only:periodic,maxp,maxalpha
 use timestep,       only:dtmax
 use io,             only:iverbose
 use units,          only:set_units
 use mpiutils,       only:reduceall_mpi
 use growth,         only:ifrag,get_vrelonvfrag,vfrag,isnow,vfragin,vfragout,rsnow,iinterpol,Tsnow

 integer, intent(inout) :: ntests,npass

 real                   :: deltax
 real                   :: dz
 real                   :: hfact
 real                   :: totmass
 real                   :: rhozero
 real                   :: errmax(9)
 real                   :: dtext_dum
 real                   :: t
 real                   :: dt
 real                   :: dtext
 real                   :: dtnew
 real                   :: csj
 real                   :: Stj
 real                   :: s
 real                   :: Vt
 real                   :: vrelonvfrag
 real                   :: csb(10000)
 real                   :: cs_snow
 real                   :: si
 real                   :: so
 real                   :: r
 real                   :: T08
 real                   :: tau
 real                   :: sinit
 real                   :: dens
 real                   :: s_in(10000)
 real                   :: s_out(10000)

 integer(kind=8)        :: npartoftypetot(maxtypes)
 integer                :: switch
 integer                :: nx
 integer                :: npart_previous
 integer                :: i
 integer                :: j
 integer                :: nsteps
 integer                :: ncheck(9)
 integer                :: nerr(9)
 integer                :: ierr

 real, parameter        :: tols = 5.e-4

 logical                :: do_output = .false.

 !- initialise
 csj         = 1.
 Stj         = 1.
 vrelonvfrag = 10.
 si          = 0.
 so          = 0.
 r           = 0.
 sinit       = 1.
 dens        = 1.
 ierr        = 0

 if (periodic) then
    if (id==master) write(*,"(/,a)") '--> testing GROWINGBOX'
 else
    if (id==master) write(*,"(/,a)") '--> skipping GROWINGBOX (need -DPERIODIC)'
    return
 endif
 !
 ! setup for growingbox problem
 !
 nx          = 32
 deltax      = 1./nx
 dz          = 2.*sqrt(6.)/nx

 call set_boundary(-0.5,0.5,-0.25,0.25,-dz,dz)

 hfact       = hfact_default
 rhozero     = 1.
 totmass     = rhozero*dxbound*dybound*dzbound
 npart       = 0
 vxyzu       = 0.
 fxyzu       = 0.
 fext        = 0.
 divcurlB    = 0.
 divcurlv    = 0.
 dustfrac    = 0.
 ddustprop   = 0.
 ddustevol   = 0.
 Bevol       = 0.
 dBevol      = 0.
 temperature = 0.
 if (maxalpha==maxp) alphaind(:,:) = 0.

 npart_previous = npart
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,verbose=.false.)
 call set_units(mass=solarm,dist=au,G=1.d0)
 do i=npart_previous+1,npart
    call set_particle_type(i,idust)
    dustprop(1,i)      = sinit
    dustprop(2,i)      = dens
    dustprop(3,i)      = 0.
    dustprop(4,i)      = 0.
    St(i)              = Stj
    csound(i)          = csj
 enddo
 npartoftype(idust)    = npart - npart_previous
 npartoftypetot(idust) = reduceall_mpi('+',npartoftype(idust))
 massoftype(idust)     = totmass/npartoftypetot(idust)
 massoftype(igas)      = 0.
 !
 ! runtime parameters
 !
 ifrag      = 0
 isnow      = 0
 shearparam = 1.
 iverbose   = 0
 ieos       = 1
 polyk      = 1.
 gamma      = 1.
 iinterpol  = .false.
 !
 !
 !
 dt         = 1.e-3
 nsteps     = 100
 t          = 0
 dtmax      = nsteps*dt
 !
 ! run growingbox problem
 !
 ncheck(:)  = 0
 nerr(:)    = 0
 errmax(:)  = 0.
 !
 ! usefull variables used along the tests
 !
 Vt         = sqrt(2**(0.5)*shearparam*Ro)*csj
 vfrag      = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vfrag < vrel : fragmentation
 tau        = 1/(sqrt(2**1.5*Ro*shearparam))
 switch     = abs(nsteps/2)
 !
 ! testing pure growth with St=cst & St=f(size)
 !
 write(*,"(/,a)")'------------------ pure growth (ifrag = 0, St = const) ------------------'
 !
 ! call deriv the first time around
 !
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
            Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)
 do i=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    s = sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    do j=1,npart
       call checkvalbuf(dustprop(1,j)/s,1.,tols,'size',nerr(1),ncheck(1),errmax(1))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(1),nerr(1),errmax(1),tols)
 write(*,"(/,a)")'------------------ ifrag = 0, St=f(size) inspired from Laibe et al. (2008) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 St(:)         = Stj
 csound(:)     = csj

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t     = t + dt
    dtext = dt
    T08   = t/tau+2*sqrt(sinit)*(1+sinit/3)
    call step(npart,npart,t,dt,dtext,dtnew)
    s     = (8+9*T08**2+3*T08*sqrt(16+9*T08**2))**(1./3.)/2 + 2/(8+9*T08**2+3*T08*sqrt(16+9*T08**2))**(1./3.) - 2
    St(:) = Stj*s
    do j=1,npart
       call checkvalbuf(dustprop(1,j)/s,1.,tols,'size',nerr(2),ncheck(2),errmax(2))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(2),nerr(2),errmax(2),tols)
 !
 ! testing pure fragmentation (ifrag = 1 & ifrag = 2)
 !
 write(*,"(/,a)")'------------------ pure fragmentation (ifrag = 1, St = const) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 St(:)         = Stj
 ifrag         = 1

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    s = sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    do j=1,npart
       call checkvalbuf(dustprop(1,j)/s,1.,tols,'size',nerr(3),ncheck(3),errmax(3))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(3),nerr(3),errmax(3),tols)

 write(*,"(/,a)")'------------------ pure fragmentation (ifrag = 2, St = const) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 ifrag         = 2

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t     = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    s     = sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*vrelonvfrag**2/(vrelonvfrag**2+1)*t
    do j=1,npart
       call checkvalbuf(dustprop(1,j)/s,1.,tols,'size',nerr(4),ncheck(4),errmax(4))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(4),nerr(4),errmax(4),tols)

 !
 ! testing growth inside the snow line and fragmentation outside of it
 !
 write(*,"(/,a)")'------------------ position based snow line ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 vfragin       = vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel < vfrag : growth
 vfragout      = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel > vfrag : fragmentation
 isnow         = 1
 rsnow         = 0.2

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t     = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    si    = sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    so    = sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    if (do_output) call write_file(i,dt,xyzh,dustprop/sinit,csb,npart,'snowline_pos_')
    do j=1,npart
       r  = sqrt(xyzh(1,j)**2+xyzh(2,j)**2+xyzh(3,j)**2)
       if (r < rsnow) call checkvalbuf(dustprop(1,j)/si,1.,10*tols,'size',nerr(6),ncheck(6),errmax(6))
       if (r > rsnow) call checkvalbuf(dustprop(1,j)/so,1.,10*tols,'size',nerr(7),ncheck(7),errmax(7))
    enddo
 enddo
 call checkvalbuf_end('size match exact solution (in)',ncheck(6),nerr(6),errmax(6),10*tols)
 call checkvalbuf_end('size match exact solution (out)',ncheck(7),nerr(7),errmax(7),10*tols)
 !
 ! testing growth inside the snow line and fragmentation outside of it
 !
 write(*,"(/,a)")'------------------ temperature based snow line ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 vfragin       = vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel < vfrag : growth
 vfragout      = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel > vfrag : fragmentation
 isnow         = 2
 Tsnow         = 300000.
 ieos          = 3
 polyk         = 0.1

 call init_eos(ieos,ierr)
 cs_snow       = sqrt(Tsnow/(temperature_coef*gmw))

 t = 0

 do j=1,npart
    csound(j) = get_spsound(ieos,xyzh(:,j),rhozero,vxyzu(:,j))
 enddo

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t           = t + dt
    dtext       = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    do j=1,npart
       csb(j)   = get_spsound(ieos,xyzh(:,j),rhozero,vxyzu(:,j))
       Vt       = sqrt(2**(0.5)*shearparam*Ro)*csb(j)
       s_in(j)  = sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
       s_out(j) = sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
       if (csb(j) > cs_snow) call checkvalbuf(dustprop(1,j)/s_in(j),1.,10*tols,'size',nerr(8),ncheck(8),errmax(8))
       if (csb(j) < cs_snow) call checkvalbuf(dustprop(1,j)/s_out(j),1.,10*tols,'size',nerr(9),ncheck(9),errmax(9))
    enddo
    if (do_output) call write_file(i,dt,xyzh,dustprop/sinit,csb/cs_snow,npart,'snowline_temp_')
 enddo
 call checkvalbuf_end('size match exact solution (in)',ncheck(8),nerr(8),errmax(8),10*tols)
 call checkvalbuf_end('size match exact solution (out)',ncheck(9),nerr(9),errmax(9),10*tols)

 if (all(nerr(1:9)==0)) npass = npass + 1
 ntests = ntests + 1

end subroutine test_growingbox

subroutine check_stokes_number(ntests,npass)
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:igas,idust,npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,dustevol,ddustevol,temperature,iphase,iamdust,maxtypes,St,xyzmh_ptmass,csound
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use energies,       only:compute_energies
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use eos,            only:ieos,polyk,gamma
 use dust,           only:K_code,idrag
 use growth,         only:ifrag,iinterpol
 use options,        only:alpha,alphamax
 use unifdis,        only:set_unifdis
 use dim,            only:periodic,mhd,use_dust
 use timestep,       only:dtmax
 use io,             only:iverbose
 use mpiutils,       only:reduceall_mpi
 use physcon,        only:au,solarm
 use units,          only:set_units,unit_velocity

 integer,intent(inout) :: ntests,npass

 integer(kind=8) :: npartoftypetot(maxtypes)

 integer         :: nx
 integer         :: itype
 integer         :: npart_previous
 integer         :: i
 integer         :: j
 integer         :: k
 integer         :: nsteps
 integer         :: ncheck(12)
 integer         :: nerr(12)

 real            :: deltax
 real            :: dz
 real            :: hfact
 real            :: totmass
 real            :: rhozero
 real            :: errmax(12)
 real            :: dtext_dum
 real            :: Stcomp
 real            :: cscomp
 real            :: r
 real            :: sinit
 real            :: dens
 real            :: t
 real            :: tmax
 real            :: dt
 real            :: dtext
 real            :: dtnew

 real, parameter :: tolst = 5.e-3
 real, parameter :: tolcs = 5.e-3

 sinit = 1.
 dens  = 1.

 write(*,"(/,a)")'--> testing INTERPOLATIONS'

 !
 ! initialise
 !
 dustprop(1,:)     = sinit
 dustprop(2,:)     = dens
 dustprop(3,:)     = 0.
 dustprop(4,:)     = 0.
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
       if (itype == idust) then
          vxyzu(:,i)       = 1.
       else
          vxyzu(:,i)       = 0.
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
 ieos      = 1
 idrag     = 2
 ifrag     = -1
 gamma     = 1.
 alpha     = 0.
 alphamax  = 0.
 iverbose  = 0
 iinterpol = .true.
 dt        = 1.e-3
 nsteps    = 100
 tmax      = nsteps*dt
 ncheck(:) = 0
 nerr(:)   = 0
 errmax(:) = 0.

 call init_step(npart,t,dtmax)

 !
 ! call derivs the first time around
 !
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
      Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)

 do k=1,6

    St(:)     = 0.
    csound(:) = 0.
    t         = 0.
    K_code    = 10.**(-k+3.5)
    cscomp    = (10.**(k/2.+1.))/unit_velocity
    polyk     = cscomp**2.

    write(*,"(/,a,es8.1,a,es8.2,a)")'-------- St ~ ',1./(2.*K_code),' , cs = ',cscomp*unit_velocity,' m/s --------'

    !
    ! run dustybox problem
    !
    do i=1,nsteps
       t     = t + dt
       dtext = dt
       call step(npart,npart,t,dt,dtext,dtnew)

       do j=1,npart
          if (iamdust(iphase(j))) then
             r      = sqrt(xyzh(1,j)**2+xyzh(2,j)**2+xyzh(3,j)**2)
             Stcomp = 1./(2.*K_code*r**(1.5))
             call checkvalbuf(St(j)/Stcomp,1.,tolst,'St',nerr(k),ncheck(k),errmax(k))
             call checkvalbuf(csound(j)/cscomp,1.,tolcs,'csound',nerr(k+6),ncheck(k+6),errmax(k+6))
          endif
       enddo
    enddo
    call checkvalbuf_end('Stokes number interpolation match exact solution',ncheck(k),nerr(k),errmax(k),tolst)
    call checkvalbuf_end('Sound speed interpolation match exact solution',ncheck(k+6),nerr(k+6),errmax(k+6),tolcs)
 enddo

 ntests = ntests + 1
 if (all(nerr(1:12)==0)) npass = npass + 1

end subroutine check_stokes_number
!---------------------------------------------------
!+
!  write an output file with x, y, z ,
!  dustprop(1) (size) and dustprop(4) (vrel/vfrag)
!+
!---------------------------------------------------
subroutine write_file(step,dt,xyzh,dustprop,cs,npart,prefix)
 real, intent(in)              :: dt
 real, intent(in)              :: xyzh(:,:),dustprop(:,:),cs(:)
 character(len=*), intent(in)  :: prefix
 integer, intent(in)           :: npart,step
 character(len=30)             :: filename,str
 integer                       :: i,lu

 write(str,"(i000.4)") step
 filename = prefix//trim(adjustl(str))//'.txt'
 open(newunit=lu,file=filename,status='replace')
 write(lu,*) step*dt
 do i=1,npart
    write(lu,*) xyzh(1,i),xyzh(2,i),xyzh(3,i),dustprop(1,i),dustprop(4,i),cs(i)
 enddo
 close(lu)

end subroutine write_file
#endif
#endif

end module testgrowth
