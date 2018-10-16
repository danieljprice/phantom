!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testnimhd
!
!  DESCRIPTION:
!   Unit tests of the non-ideal MHD algorithms; this requires the NICIL
!   module, but does not test its inner workings
!   Note that this module should also test super-timestepping, but that
!   is less straight forward
!
!  REFERENCES:
!   Wurster, Price & Ayliffe (2014),  MNRAS 444, 1104
!   Wurster, Price & Bate (2016), MNRAS 457, 1037
!   Wurster (2016), PASA 33, e041
!
!  OWNER: James Wurster
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, deriv, dim, eos, io, kernel, nicil, options,
!    part, physcon, step_lf_global, testutils, timestep, timestep_sts,
!    unifdis, units
!+
!--------------------------------------------------------------------------
module testnimhd
 use testutils, only:checkval
 use io,        only:id,master
#ifdef STS_TIMESTEPS
 use timestep_sts,   only:use_sts
#endif
 implicit none
 public :: test_nonidealmhd

 private
#ifndef STS_TIMESTEPS
 logical :: use_sts
#endif

contains
!--------------------------------------------------------------------------
subroutine test_nonidealmhd(ntests,npass)
 integer, intent(inout) :: ntests,npass
#ifndef MPI
#ifdef NONIDEALMHD

 if (id==master) write(*,"(/,a)") '--> TESTING NON-IDEAL MHD ALGORITHMS'
 !
 ! Test decay of a wave by ambipolar diffusion
 !
 call test_wavedamp(ntests,npass)
 !
 ! Test Standing Shock using the Hall effect
 !
 call test_standingshock(ntests,npass)
 !
 ! Test the passage of arrays required for non-constant eta
 !
 call test_narrays(ntests,npass)
 !
 if (id==master) write(*,"(/,a)") '<-- NON-IDEALMHD TEST COMPLETE'
#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING NON-IDEAL MHD TEST (REQUIRES -DNONIDEALMHD)'
#endif
#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING NON-IDEAL MHD TEST (MPI NOT SUPPORTED)'
#endif

end subroutine test_nonidealmhd

#ifndef MPI
!--------------------------------------------------------------------------
!+
!  Tests the decay of a wave using ambipolar diffusion
!  (e.g. Fig 3 of Wurster, Price & Ayliffe 2014, with the exact solution
!  given in Eqn. 49)
!+
!--------------------------------------------------------------------------
subroutine test_wavedamp(ntests,npass)
 use physcon,        only:pi
 use units,          only:set_units,utime,udist,umass,unit_Bfield
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:npart,xyzh,vxyzu,Bxyz,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustfrac,ddustevol,temperature,igas,alphaind,&
                          dustprop,ddustprop
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use testutils,      only:checkval
 use eos,            only:ieos,polyk,gamma
 use options,        only:alphaB,alpha,alphamax
 use unifdis,        only:set_unifdis
 use dim,            only:periodic
 use timestep,       only:dtmax,tmax
 use io,             only:iverbose
 use nicil,          only:nicil_initialise,eta_constant,eta_const_type,icnstsemi, &
                          use_ohm,use_hall,use_ambi,C_AD
#ifdef STS_TIMESTEPS
 use timestep_sts,   only:sts_initialise
 use timestep,       only:dtdiff,dtcourant
#endif
 integer, intent(inout) :: ntests,npass
 integer                :: i,j,nx,nsteps,ierr
 integer                :: nerr(4)
 real                   :: deltax,x_min,y_min,z_min,kx,rhozero,Bx0,vA,vcoef,totmass
 real                   :: t,dt,dtext_dum,dtext,dtnew
 real                   :: L2,h0,quada,quadb,quadc,omegaI,omegaR,Bzrms_num,Bzrms_ana
 real, parameter        :: tol     = 7.15d-5
 real, parameter        :: toltime = 4.00d-4
 logical                :: valid_dt
 logical                :: print_output = .false.

#ifndef ISOTHERMAL
 if (id==master) write(*,"(/,a)") '--> skipping wave dampening test (need -DISOTHERMAL)'
#endif
 if (periodic) then
    if (id==master) write(*,"(/,a)") '--> testing wave dampening test'
 else
    if (id==master) write(*,"(/,a)") '--> skipping wave dampening test (need -DPERIODIC)'
    return
 endif
 !
 ! initialise values for grid
 !
 nx      = 32
 deltax  = 1./nx
 x_min   = -0.5
 y_min   = x_min*sqrt(3.0)/2.0
 z_min   = x_min*sqrt(6.0)/3.0
 call set_boundary(x_min,-x_min,y_min,-y_min,z_min,-z_min)
 kx      = 2.0*pi/(xmax-xmin)
 rhozero = 1.0
 Bx0     = 1.0
 vA      = Bx0/sqrt(rhozero)
 vcoef   = 0.01*vA
 totmass = rhozero*dxbound*dybound*dzbound
 npart   = 0
#ifdef STS_TIMESTEPS
 call sts_initialise(ierr,dtdiff)
#else
 if (id==master) write(*,"(/,a)") '--> skipping super-timestepping portion of test (need -DSTS_TIMESTEPS)'
#endif
 !
 ! set particles
 !
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                   deltax,hfact_default,npart,xyzh,verbose=.false.)
 vxyzu = 0.0
 Bxyz  = 0.0
 Bevol = 0.0
 do i=1,npart
    call set_particle_type(i,igas)
    vxyzu(3,i) = vcoef*sin(kx*(xyzh(1,i)-xmin))
    Bxyz(1,i)  = Bx0
 enddo
 Bevol(1:3,:)      = Bxyz(1:3,:)/rhozero
 npartoftype(igas) = npart
 massoftype(igas)  = totmass/npartoftype(igas)
 !
 ! initialise runtime parameters
 !
 ieos           = 1
 polyk          = 1.0
 gamma          = 1.0
 alphaB         = 0.0
 alpha          = 0.0
 alphamax       = 0.0
 alphaind       = real(alpha,kind=kind(alphaind(1,1)))
 iverbose       = 0
 eta_const_type = icnstsemi
 eta_constant   = .true.
 use_ohm        = .false.
 use_hall       = .false.
 use_ambi       = .true.
 C_AD           = 0.010
 L2             = 0.0
 dt             = 4.0e-3
 nsteps         = 125
 t              = 0
 tmax           = nsteps*dt
 !
 ! initial values for the exact solution
 h0      = vcoef*Bx0/(vA*sqrt(2.0))
 quada   =  1.0
 quadb   =  (2.0*vA*pi)**2*C_AD
 quadc   = -(2.0*vA*pi)**2
 omegaI  = -0.5*quadb
 omegaR  =  0.5*sqrt( -quadb*quadb - 4.0*quada*quadc )
 !
 ! initialise the Nicil library
 call set_units(mass=1.d0,dist=1.d0)
 call nicil_initialise(real(utime),real(udist),real(umass),real(unit_Bfield),ierr)
 !
 ! call derivs the first time around
 use_sts = .true.
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)
 use_sts = .false.  ! Since we only want to run supertimestepping once to verify we get the correct dt
 !
 ! run wave damp problem
 !
 call init_step(npart,t,dtmax)
 nsteps   = 0
 valid_dt = .true.
 call init_step(npart,t,dtmax)
 do while (valid_dt .and. t <= tmax)
    t      = t + dt
    nsteps = nsteps + 1
    dtext  = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    ! calculate the averaged root-mean-square magnetic field in the z-direction
    Bzrms_num = 0.0
    do j = 1,npart
       Bzrms_num = Bzrms_num + Bxyz(3,j)**2
    enddo
    Bzrms_num = sqrt(Bzrms_num/npart)
    Bzrms_ana = h0*abs(sin(omegaR*t))*exp(omegaI*t)
    L2        = L2 + (Bzrms_ana - Bzrms_num)**2
    if (dtnew < dt) valid_dt = .false.
 enddo
 ! For printing outputs if further debugging is required.
 if (print_output) then
    open(unit=111,file='nimhd_wavedamp.dat')
    do i = 1,npart
       write(111,'(10Es18.6)') xyzh(:,i),Bxyz(1:3,i),vxyzu(1:3,i)
    enddo
    close(111)
 endif

 L2 = sqrt(L2/nsteps)
 call checkval(L2,0.0,tol,nerr(1),'L2 error on wave damp test')
 call checkval(valid_dt,.true.,nerr(2),'dt to ensure above valid default')
#ifdef STS_TIMESTEPS
 write(*,'(1x,a,3Es18.11)') 'dtcourant, dtdiff: ',dtcourant,dtdiff
 call checkval(dtcourant,4.51922587536d-3,toltime,nerr(3),'initial courant dt')
 call checkval(dtdiff,   2.88824049903d-2,toltime,nerr(4),'initial dissipation dt from sts')
#endif

 ntests = ntests + 1
 if (all(nerr==0)) npass = npass + 1

end subroutine test_wavedamp
!--------------------------------------------------------------------------
!+
!  Tests the Hall effect in a standing shock
!  (e.g. Fig C2 of Wurster, Price & Bate 2016)
!  This will also test that boundary particles are correctly implemented
!+
!--------------------------------------------------------------------------
subroutine test_standingshock(ntests,npass)
 use physcon,        only:pi
 use units,          only:set_units,utime,udist,umass,unit_Bfield
 use boundary,       only:set_boundary,ymin,ymax,zmin,zmax,dybound,dzbound
 use kernel,         only:hfact_default,radkern
 use part,           only:npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,hrho,rhoh,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustfrac,ddustevol,igas,iboundary,&
                          set_boundaries_to_active,alphaind,maxalpha,maxp,iphase,Bxyz,dustprop,ddustprop,temperature
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use testutils,      only:checkval
 use eos,            only:ieos,polyk,gamma
 use options,        only:alpha,alphamax,alphaB
 use unifdis,        only:set_unifdis,get_ny_nz_closepacked
 use dim,            only:periodic
 use timestep,       only:dtmax,tmax
 use io,             only:iverbose
 use nicil,          only:nicil_initialise,eta_constant,eta_const_type,icnstsemi, &
                          use_ohm,use_hall,use_ambi,C_OR,C_HE,C_AD
 integer, intent(inout) :: ntests,npass
 integer                :: i,nx,ny,nz,nsteps,ierr,idr,npts
 integer                :: nerr(5)
 real                   :: volume,totmass,fac,xleft,xright,dxleft,dxright,yleft,yright,zleft,zright,rhoi
 real                   :: t,dt,dtext_dum,dtext,dtnew
 real                   :: dexact,bexact,vexact,L2d,L2v,L2b,dx
 real                   :: leftstate(8),rightstate(8),exact_x(51),exact_d(51),exact_vx(51),exact_by(51)
 real, parameter        :: told = 2.1d-2, tolv=3.05d-2, tolb=1.1d-1
 logical                :: valid_dt
 logical                :: print_output = .false.
 logical                :: valid_bdy

#ifndef ISOTHERMAL
 if (id==master) write(*,"(/,a)") '--> skipping standing shock test (need -DISOTHERMAL)'
#endif
 if (periodic) then
    if (id==master) write(*,"(/,a)") '--> testing standing shock & boundary particles'
 else
    if (id==master) write(*,"(/,a)") '--> skipping standing shock test (cannot have -DPERIODIC)'
    return
 endif
 !
 ! initialise values for grid
 !
 nx             = 32
 tmax           = 1.0
 leftstate      = (/1.7942,0.017942,-0.9759,-0.6561,0.,1.,1.74885,0./)
 rightstate     = (/1.    ,0.01    ,-1.7510, 0.    ,0.,1.,0.6    ,0./)
 xleft          = -0.75
 xright         = -xleft - rightstate(3)*tmax
 dxleft         = -xleft/float(nx)
 dxright        = dxleft*(leftstate(1)/rightstate(1))**(1./3.)
 fac            = -6.*(int(1.99*radkern/6.) + 1)*max(dxleft,dxright)
 yleft          = fac*sqrt(0.75)
 zleft          = fac*sqrt(6.)/3.
 yright         = -yleft
 zright         = -zleft
 npart          = 0
 npartoftype(:) = 0
 use_sts        = .false.
 call set_boundary(xleft-1000.*dxleft,xright+1000.*dxright,yleft,yright,zleft,zright) ! false periodicity in x
 !
 ! set particles on the left half of the shock
 !
 call set_unifdis('closepacked',id,master,xleft,0.0,ymin,ymax,zmin,zmax,dxleft,hfact_default,npart,xyzh)
 volume                = -xleft*dybound*dzbound
 totmass               = volume*leftstate(1)
 massoftype(igas)      = totmass/npart
 massoftype(iboundary) = massoftype(igas)
 !
 ! set particles of the right half of the shock
 !
 call get_ny_nz_closepacked(dxright,ymin,ymax,zmin,zmax,ny,nz)
 volume  = xright*dybound*dzbound
 totmass = volume*rightstate(1)
 dxright = xright/((totmass/massoftype(igas))/(ny*nz))
 call set_unifdis('closepacked',id,master,0.0,xright,ymin,ymax,zmin,zmax,dxright,hfact_default,npart,xyzh,npy=ny,npz=nz)
 !
 ! set boundary particles, and set properties of the particles
 !
 fac   = nint(2.01*radkern*hfact_default)
 vxyzu = 0.0
 Bxyz  = 0.0
 Bevol = 0.0
 do i=1,npart
    ! set boundaries
    if ( (xyzh(1,i) < (xleft + fac*dxleft ) ) .or. (xyzh(1,i) > (xright - fac*dxright) ) ) then
       call set_particle_type(i,iboundary)
       npartoftype(iboundary) = npartoftype(iboundary) + 1
    else
       call set_particle_type(i,igas)
       npartoftype(igas) = npartoftype(igas) + 1
    endif
    ! set properties
    if (xyzh(1,i) > 0.) then
       xyzh(4,i)    = hrho(rightstate(1),massoftype(igas))
       vxyzu(1:3,i) = rightstate(3:5)
       Bxyz(1:3,i)  = rightstate(6:8)
       Bevol(1:3,i) = rightstate(6:8)/rightstate(1)
    else
       xyzh(4,i)    = hrho(leftstate(1),massoftype(igas))
       vxyzu(1:3,i) = leftstate(3:5)
       Bxyz(1:3,i)  = leftstate(6:8)
       Bevol(1:3,i) = leftstate(6:8)/leftstate(1)
    endif
 enddo
 !
 ! initialise runtime parameters
 !
 ieos     = 1
 use_ohm  = .true.
 use_hall = .true.
 use_ambi = .true.
 C_OR     =  1.12d-9
 C_HE     = -3.53d-2
 C_AD     =  7.83d-3
 alpha    =  0.0
 alphamax =  1.0
 alphaB   =  1.0
 fext     =  0.0
 polyk    =  0.01
 gamma    =  1.0
 dt       =  2.0d-3
 nsteps   =  500
 t        =  0
 iverbose =  0
 tmax     = nsteps*dt
 eta_const_type           = icnstsemi
 eta_constant             = .true.
 set_boundaries_to_active = .true.
 if (maxalpha==maxp) alphaind(1,:) = real(alpha,kind=4)
 !
 ! initialise the Nicil library
 call set_units(mass=1.0d0,dist=1.0d0)
 call nicil_initialise(real(utime),real(udist),real(umass),real(unit_Bfield),ierr)
 !
 ! call derivs the first time around
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)
 set_boundaries_to_active = .false.
 !
 ! run standing shock problem
 !
 nsteps   = 0
 valid_dt = .true.
 call init_step(npart,t,dtmax)
 do while (valid_dt .and. t <= tmax)
    t      = t + dt
    nsteps = nsteps + 1
    dtext  = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    if (dtnew < dt) valid_dt = .false.
 enddo
 !
 ! Compare to exact solution
 dx = 0.01
 exact_x =(/ 0.0000, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, &
             0.0900, 0.1000, 0.1100, 0.1200, 0.1300, 0.1400, 0.1500, 0.1600, 0.1700, &
             0.1800, 0.1900, 0.2000, 0.2100, 0.2200, 0.2300, 0.2400, 0.2500, 0.2600, &
             0.2700, 0.2800, 0.2900, 0.3000, 0.3100, 0.3200, 0.3300, 0.3400, 0.3500, &
             0.3600, 0.3700, 0.3800, 0.3900, 0.4000, 0.4100, 0.4200, 0.4300, 0.4400, &
             0.4500, 0.4600, 0.4700, 0.4800, 0.4900, 0.5000/)
 exact_d =(/ 1.4720, 1.4242, 1.3733, 1.3202, 1.2663, 1.2128, 1.1614, 1.1135, 1.0706, &
             1.0338, 1.0040, 0.9816, 0.9667, 0.9587, 0.9571, 0.9607, 0.9683, 0.9785, &
             0.9900, 1.0016, 1.0123, 1.0213, 1.0279, 1.0320, 1.0335, 1.0325, 1.0295, &
             1.0249, 1.0193, 1.0131, 1.0070, 1.0014, 0.9966, 0.9928, 0.9903, 0.9890, &
             0.9888, 0.9896, 0.9912, 0.9933, 0.9956, 0.9980, 1.0003, 1.0022, 1.0038, &
             1.0048, 1.0053, 1.0053, 1.0049, 1.0042, 1.0033/)
 exact_vx=(/-1.1895,-1.2294,-1.2750,-1.3262,-1.3827,-1.4437,-1.5077,-1.5725,-1.6355,&
            -1.6937,-1.7440,-1.7838,-1.8114,-1.8263,-1.8294,-1.8226,-1.8084,-1.7895,&
            -1.7686,-1.7481,-1.7296,-1.7145,-1.7034,-1.6966,-1.6942,-1.6958,-1.7007,&
            -1.7083,-1.7179,-1.7283,-1.7388,-1.7485,-1.7570,-1.7636,-1.7681,-1.7704,&
            -1.7708,-1.7693,-1.7666,-1.7628,-1.7586,-1.7544,-1.7504,-1.7470,-1.7444,&
            -1.7426,-1.7417,-1.7417,-1.7423,-1.7436,-1.7452/)
 exact_by=(/ 1.4984, 1.4432, 1.3769, 1.2982, 1.2062, 1.1009, 0.9834, 0.8568, 0.7262, &
             0.5988, 0.4835, 0.3893, 0.3240, 0.2920, 0.2934, 0.3239, 0.3762, 0.4413, &
             0.5107, 0.5772, 0.6357, 0.6828, 0.7170, 0.7380, 0.7464, 0.7435, 0.7309, &
             0.7109, 0.6855, 0.6573, 0.6285, 0.6014, 0.5779, 0.5593, 0.5466, 0.5401, &
             0.5395, 0.5439, 0.5524, 0.5636, 0.5762, 0.5889, 0.6007, 0.6107, 0.6185, &
             0.6237, 0.6263, 0.6265, 0.6247, 0.6211, 0.6164/)
 L2d  = 0.0
 L2v  = 0.0
 L2b  = 0.0
 npts = 0
 valid_bdy = .true.
 ! For printing outputs if further debugging is required.
 if (print_output) then
    open(unit=112,file='nimhd_shock_particles.dat')
    open(unit=113,file='nimhd_shock_solution.dat')
 endif
 do i = 1,npart
    rhoi = rhoh(xyzh(4,i),massoftype(igas))
    if (print_output) write(112,'(5Es18.6,I3)') xyzh(1:2,i),rhoi,vxyzu(1,i),Bxyz(2,i),iphase(i)
    if (exact_x(1) < xyzh(1,i) .and. xyzh(1,i) < exact_x(50) ) then
       npts   = npts + 1
       idr    = int(xyzh(1,i)/dx)+1
       dexact = exact_d (idr) + (exact_d (idr+1) - exact_d (idr))/dx*(exact_x(idr+1)-xyzh(1,i) )
       vexact = exact_vx(idr) + (exact_vx(idr+1) - exact_vx(idr))/dx*(exact_x(idr+1)-xyzh(1,i) )
       bexact = exact_by(idr) + (exact_by(idr+1) - exact_by(idr))/dx*(exact_x(idr+1)-xyzh(1,i) )
       L2d = L2d + (dexact - rhoi      )**2
       L2v = L2v + (vexact - vxyzu(1,i))**2
       L2b = L2b + (bexact - Bxyz(2,i) )**2
       if (print_output) then
          write(113,'(7Es18.6)') xyzh(1,i),rhoi,vxyzu(1,i),Bxyz(2,i),dexact,vexact,bexact
       endif
    endif
    if (xyzh(1,i) >  0.71 .and. rhoi > 0.99) valid_bdy = .false.
    if (xyzh(1,i) < -1.69 .and. rhoi > 1.78) valid_bdy = .false.
 enddo
 if (print_output) then
    close(112)
    close(113)
 endif

 if (npts > 0) then
    L2d = sqrt(L2d/npts)
    L2v = sqrt(L2v/npts)
    L2b = sqrt(L2b/npts)
 endif

 call checkval(L2d,0.0,  told,  nerr(1),'density error on standing shock, compared to analytics')
 call checkval(L2v,0.0,  tolv,  nerr(2),'v_x error on standing shock, compared to analytics')
 call checkval(L2b,0.0,  tolb,  nerr(3),'B_y error on standing shock, compared to analytics')
 call checkval(valid_dt, .true.,nerr(4),'dt to ensure above valid default')
 call checkval(valid_bdy,.true.,nerr(5),'Boundary particles are correctly initialised')

 ntests = ntests + 1
 if (all(nerr==0)) npass = npass + 1

end subroutine test_standingshock
!--------------------------------------------------------------------------
!+
!  Tests the correct passage of the number density arrays required for
!  non-ideal mhd
!+
!--------------------------------------------------------------------------
subroutine test_narrays(ntests,npass)
 use physcon,        only:pi,solarm
 use units,          only:set_units,utime,udist,umass,unit_Bfield,unit_density
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:npart,xyzh,vxyzu,Bxyz,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustfrac,ddustevol,igas,alphaind,&
                          n_R,n_electronT,rhoh,dustprop,ddustprop,eta_nimhd,iohm,ihall,iambi,temperature
 use deriv,          only:derivs
 use testutils,      only:checkval
 use eos,            only:ieos,init_eos,polyk,polyk2,gamma,get_temperature
 use options,        only:alphaB,alpha,alphamax
 use unifdis,        only:set_unifdis
 use dim,            only:periodic
 use io,             only:iverbose
 use nicil,          only:nicil_initialise,nicil_get_eta,eta_constant,use_ohm,use_hall,use_ambi,&
                          ion_rays,ion_thermal,unit_eta
 integer, intent(inout) :: ntests,npass
 integer, parameter     :: kmax = 2
 integer                :: i,k,nx,ierr
 integer                :: nerr(3*kmax)
 real                   :: deltax,x_min,y_min,z_min,totmass,cs_sphere,cs_medium
 real                   :: t,dtext_dum,Bi,rhoi,tempi
 real                   :: rho0(2),Bz0(2),eta_act(3,kmax)
 real, parameter        :: tol = 6.3e-5  ! 1.0e-7 (The higher tolerance is needed for some compilers during certain phases of the moon)
 !
 if (periodic) then
    if (id==master) write(*,"(/,a)") '--> testing calculation of non-constant eta'
 else
    if (id==master) write(*,"(/,a)") '--> skipping calculation of non-constant eta (need -DPERIODIC)'
    return
 endif
 !
 ! Initial density & magnetic field strength, along with expected eta coefficients (cgs)
 ! Note: these coefficients are calculated using nicil_get_one_value using the *output* rhoi & Bi
 !       since they shift slightly during derivs and step.
 ! Note: the first set of values is for the initial conditions of Wurster, Price & Bate (2016); the
 !       second set is high density/B-field to ensure thermal ionisation is working
 !
 call set_units(mass=solarm,dist=1.0d16,G=1.d0)
 rho0(1)      = 7.420d-18 /unit_density   ! [g/cm^3]
 Bz0(1)       = 8.130d-5  /unit_Bfield    ! [G]
 eta_act(1,1) = 1.14793940113d10          ! [cm^2/s] expected eta_ohm
 eta_act(2,1) = 3.40040077209d14          ! [cm^2/s] expected eta_hall
 eta_act(3,1) = 5.26247580402d17          ! [cm^2/s] expected eta_ambi
 rho0(2)      = 4.6d-3    /unit_density   ! [g/cm^3]
 Bz0(2)       = 1.92d2    /unit_Bfield    ! [G]
 eta_act(1,2) = 5.93454638765d8           ! [cm^2/s] expected eta_ohm
 eta_act(2,2) = 1.08059808926d4           ! [cm^2/s] expected eta_hall
 eta_act(3,2) = 4.17918319187d-3          ! [cm^2/s] expected eta_ambi
 !
 ! initialise values for grid
 !
 nx      = 8
 deltax  = 1./nx
 x_min   = -1.
 y_min   = x_min*sqrt(3.0)/2.0
 z_min   = x_min*sqrt(6.0)/3.0
 call set_boundary(x_min,-x_min,y_min,-y_min,z_min,-z_min)
 npart   = 0
 cs_sphere = 0.19       ! = 2.189E+04 cm/s
 cs_medium = 1.04       ! = 1.199E+05 cm/s
 !
 ! initialise runtime parameters
 !
 ieos         = 8
 gamma        = 1.0
 polyk        = cs_sphere**2
 polyk2       = cs_medium**2
 alphaB       = 0.0
 alpha        = 0.0
 alphamax     = 0.0
 iverbose     = 0
 eta_constant = .false.
 use_ohm      = .true.
 use_hall     = .true.
 use_ambi     = .true.
 ion_rays     = .true.
 ion_thermal  = .true.
 use_sts      = .false.
 ! initialise eos, & the Nicil library
 call init_eos(8,ierr)
 call nicil_initialise(real(utime),real(udist),real(umass),real(unit_Bfield),ierr)
 !
 !--Loop over both sets of calculations
 !
 do k = 1,kmax
    totmass      = rho0(k)*dxbound*dybound*dzbound
    n_R          = 0.0
    n_electronT  = 0.0
    !
    ! set particles
    !
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                      deltax,hfact_default,npart,xyzh,verbose=.false.)
    vxyzu = 0.0
    Bxyz  = 0.0
    Bevol = 0.0
    do i=1,npart
       call set_particle_type(i,igas)
       Bxyz(3,i) = Bz0(k)
    enddo
    Bevol(1:3,:)      = Bxyz(1:3,:)/rho0(k)
    npartoftype(igas) = npart
    massoftype(igas)  = totmass/npartoftype(igas)
    alphaind          = real(alpha,kind=kind(alphaind(1,1)))
    !
    ! call derivs, which will also calculate eta
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)
    !
    ! Calculate eta from NICIL
    rhoi  = rhoh(xyzh(4,1),massoftype(1))
    Bi    = sqrt( dot_product(Bevol(1:3,1),Bevol(1:3,1)) )*rhoi
    tempi = get_temperature(ieos,xyzh(1:3,1),rhoi,vxyzu(:,1))

    print*, ' '
    write(*,'(1x,a,3Es18.11)') 'Used   rho,B_z,temp (cgs): ',rhoi*unit_density,Bi*unit_Bfield,tempi
    write(*,'(1x,a,3Es18.11)') 'eta_ohm, eta_hall, eta_ambi (cgs): ', eta_nimhd(1:3,1)*unit_eta
    call checkval(eta_nimhd(iohm, 1)*unit_eta,eta_act(1,k),tol,nerr(3*(k-1)+1),'calculated non-constant eta_ohm')
    call checkval(eta_nimhd(ihall,1)*unit_eta,eta_act(2,k),tol,nerr(3*(k-1)+2),'calculated non-constant eta_hall')
    call checkval(eta_nimhd(iambi,1)*unit_eta,eta_act(3,k),tol,nerr(3*(k-1)+3),'calculated non-constant eta_ambi')

 enddo
 ntests = ntests + 1
 if (all(nerr==0)) npass = npass + 1

end subroutine test_narrays
!--------------------------------------------------------------------------
#endif

end module testnimhd
