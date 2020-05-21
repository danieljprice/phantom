!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Setup for 3D shock tube tests
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    dtg         -- Dust to gas ratio
!    dust_method -- 1=one fluid, 2=two fluid
!    gamma       -- Adiabatic index
!    nx          -- resolution (number of particles in x) for -xleft < x < xshock
!    polyk       -- square of the isothermal sound speed
!    xleft       -- x min boundary
!    xright      -- x max boundary
!
!  DEPENDENCIES: boundary, dim, dust, eos, infile_utils, io, kernel,
!    mpiutils, nicil, options, part, physcon, prompting, radiation_utils,
!    set_dust, setup_params, timestep, unifdis, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 integer :: nx, icase, dust_method
 real    :: xleft, xright, yleft, yright, zleft, zright
 real    :: dxleft, kappa
 character(len=100) :: shocktype
 character(len=100) :: latticetype = 'closepacked'
 integer :: nstates
 integer, parameter :: max_states = 8
 integer, parameter :: &
    idens = 1, &
    ipr   = 2, &
    ivx   = 3, &
    ivy   = 4, &
    ivz   = 5, &
    iBx   = 6, &
    iBy   = 7, &
    iBz   = 8

 character(len=4), parameter :: var_label(max_states) = &
   (/'dens','pr  ','vx  ','vy  ','vz  ','Bx  ','By  ','Bz  '/)


 real :: leftstate(max_states), rightstate(max_states)

 public  :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for shock tube problems in 3D
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:fatal,master,iprint,error
 use boundary,     only:ymin,zmin,ymax,zmax,set_boundary
 use mpiutils,     only:bcast_mpi
 use dim,          only:maxvxyzu,ndim,mhd,do_radiation,use_dust
 use options,      only:use_dustfrac
 use part,         only:labeltype,set_particle_type,igas,iboundary,hrho,Bxyz,mhd,&
                        periodic,dustfrac,gr,ndustsmall,ndustlarge,ndusttypes,ikappa
 use part,         only:rad,radprop,iradxi,ikappa
 use eos,          only:gmw
 use kernel,       only:radkern,hfact_default
 use timestep,     only:tmax
 use prompting,    only:prompt
 use set_dust,     only:set_dustfrac
 use units,        only:set_units,unit_opacity
 use dust,         only:idrag
 use unifdis,      only:is_closepacked,is_valid_lattice
#ifdef NONIDEALMHD
 use nicil,           only:rho_i_cnst
#endif
 use physcon,         only:au,solarm
 use radiation_utils, only:radiation_and_gas_temperature_equal
 use setshock,     only:set_shock,adjust_shock_boundaries
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npartoftype(:)
 integer,           intent(inout) :: npart
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real                             :: delta,gam1,xshock,fac,dtg
 real                             :: uuleft,uuright,xbdyleft,xbdyright,dxright,rholeft,rhoright
 integer                          :: i,ierr,nbpts,iverbose
 character(len=120)               :: shkfile, filename
 logical                          :: iexist,use_closepacked

 if (gr) call set_units(G=1.,c=1.)
 if (do_radiation) call set_units(dist=au,mass=solarm,G=1.d0)
 !
 ! quit if not periodic
 !
 if (.not.periodic) call fatal('setup','require PERIODIC=yes')
 !
 ! verify a legitimate lattice type has been chosen
 !
 if (.not.is_valid_lattice(latticetype)) then
    call fatal('setup','invalid lattice type')
 endif
 use_closepacked = is_closepacked(latticetype)
 !
 ! determine if an .in file exists
 !
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 !
 ! general parameters
 !
 time  = 0.0
 gamma = 5.0/3.0
 polyk = 0.1
 kappa = 1.e6
 if (.not.iexist) hfact = hfact_default
 nstates = max_states
 if (.not.mhd) nstates = max_states - 3
 if (use_dust) dust_method = 1 ! one fluid by default
 !
 ! setup particles
 !
 npart          = 0
 npart_total    = 0
 npartoftype(:) = 0

 !--zero dust to gas ratio in case dust is not being used
 dtg = 0.
 !
 ! read shock parameters from the .setup file.
 ! if file does not exist, then ask for user input
 !
 shkfile = trim(fileprefix)//'.setup'
 call read_setupfile(shkfile,iprint,nstates,gamma,polyk,dtg,ierr)
 if (ierr /= 0 .and. id==master) then
    call choose_shock(gamma,polyk,dtg,iexist) ! Choose shock
    call write_setupfile(shkfile,iprint,nstates,gamma,polyk,dtg)       ! write shock file with defaults
 endif

 rholeft  = get_conserved_density(leftstate)
 rhoright = get_conserved_density(rightstate)

 !
 ! choose dust method (from .setup file)
 !
 if (use_dust) then
    idrag = 2
    use_dustfrac = (dust_method == 1)
    if (dust_method==1) then
       ndustsmall = 1
       ndustlarge = 0
    elseif (dust_method==2) then
       ndustsmall = 0
       ndustlarge = 1
    endif
    ndusttypes = ndustsmall + ndustlarge
 endif
 !
 ! for one fluid dust, density is TOTAL density of gas + dust
 !
 if (use_dustfrac) then
    rholeft = rholeft*(1. + dtg)
    rhoright = rhoright*(1. + dtg)
 endif
 !
 ! print setup parameters
 !
 if (id==master) call print_shock_params(nstates)
 !
 ! adjust boundaries to allow space for boundary particles and inflow
 !
 dxleft = -xleft/float(nx)
 xshock = 0.5*(xleft + xright)
 call adjust_shock_boundaries(dxleft,dxright,radkern, &
      leftstate(ivx),rightstate(ivx),rholeft,rhoright,tmax,ndim,&
      xleft,xright,yleft,yright,zleft,zright,is_closepacked(latticetype))
 !
 ! set domain boundaries - use a longer box in x direction since
 ! it is not actually periodic in x
 !
 call set_boundary(xleft-1000.*dxleft,xright+1000.*dxright,yleft,yright,zleft,zright)
 !
 ! setup gas particles
 !
 iverbose = 1
 call set_shock(latticetype,id,master,igas,rholeft,rhoright,xleft,xright,ymin,ymax,zmin,zmax,&
                xshock,dxleft,hfact,npart,xyzh,massoftype,iverbose,ierr)
 if (ierr /= 0) call fatal('setup','errors in shock setup')

 ! define rhozero as density in left half; required for certain simulations (e.g. non-ideal MHD with constant resistivity)
 rhozero = rholeft
 !
 ! Fix the particles near x-boundary; else define as gas
 !
 nbpts     = 0
 fac       = nint(2.01*radkern*hfact)
 xbdyleft  = fac*dxleft
 xbdyright = fac*dxright
 do i=1,npart
    if ( (xyzh(1,i) < (xleft + xbdyleft ) ) .or. &
         (xyzh(1,i) > (xright - xbdyright) ) ) then
       call set_particle_type(i,iboundary)
       nbpts = nbpts + 1
    else
       call set_particle_type(i,igas)
    endif
 enddo
 massoftype(iboundary)  = massoftype(igas)
 npartoftype(iboundary) = nbpts
 npartoftype(igas)      = npart - nbpts
 write(iprint,'(1x,a,i8)') 'Setup_shock: npart     = ',npart
 write(iprint,'(1x,a,i8)') 'Setup_shock: ngas      = ',npartoftype(igas)
 write(iprint,'(1x,a,i8)') 'Setup_shock: nboundary = ',nbpts
 !
 ! now set particle properties
 !
 gam1 = gamma - 1.
 if (abs(gam1) > 1.e-3) then
    uuleft  = leftstate(ipr)/(gam1*leftstate(idens))
    uuright = rightstate(ipr)/(gam1*rightstate(idens))
 else
    uuleft  = 3.*leftstate(ipr)/(2.*leftstate(idens))
    uuright = 3.*rightstate(ipr)/(2.*rightstate(idens))
 endif

 Bxyz = 0.
 vxyzu = 0.
 dustfrac = 0.
 do i=1,npart
    delta = xyzh(1,i) - xshock
    if (delta > 0.) then
       xyzh(4,i)  = hrho(rhoright,massoftype(igas))
       vxyzu(1,i) = rightstate(ivx)
       vxyzu(2,i) = rightstate(ivy)
       vxyzu(3,i) = rightstate(ivz)
       if (maxvxyzu >= 4) vxyzu(4,i) = uuright
       if (mhd) Bxyz(1:3,i) = rightstate(iBx:iBz)
       if (do_radiation) then
          rad(iradxi,i) = radiation_and_gas_temperature_equal(rhoright,uuright,gamma,gmw)
       endif
    else
       xyzh(4,i)  = hrho(rholeft,massoftype(igas))
       vxyzu(1,i) = leftstate(ivx)
       vxyzu(2,i) = leftstate(ivy)
       vxyzu(3,i) = leftstate(ivz)
       if (maxvxyzu >= 4) vxyzu(4,i) = uuleft
       if (mhd) Bxyz(1:3,i) = leftstate(iBx:iBz)
       if (do_radiation) then
          rad(iradxi,i) = radiation_and_gas_temperature_equal(rholeft,uuleft,gamma,gmw)
       endif
    endif
    !
    ! one fluid dust: set dust fraction on gas particles
    !
    if (use_dustfrac) call set_dustfrac(dtg,dustfrac(:,i))
    !
    ! radiation: set opacity
    if (do_radiation) radprop(ikappa,i) = kappa/unit_opacity
 enddo
 if (mhd) ihavesetupB = .true.
 !
 ! set up dust (as separate set of particles)
 !
 if (use_dust .and. .not.use_dustfrac) then
    if (dtg > 0.) call set_dust_particles(dtg,npart,npartoftype,massoftype,xyzh,vxyzu,ierr)
    if (ierr /= 0) call error('setup','could not set up dust particles')
 endif
 write(iprint,'(1x,a,es16.8)') 'Setup_shock: mass of gas & boundary particles   = ', massoftype(igas)

#ifdef NONIDEALMHD
 !Modify ion density from fraction to physical (for ambipolar diffusion)
 if (.not.iexist) rho_i_cnst = rho_i_cnst * rhozero
#endif

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  setup dust particles (at the moment setup only allows for fixed
!  dust-to-gas ratio, so we just make copies of the gas particles)
!+
!-----------------------------------------------------------------------
subroutine set_dust_particles(dtg,npart,npartoftype,massoftype,xyzh,vxyzu,ierr)
 use part, only:iamtype,iphase,maxp,maxphase,igas,idust,iboundary,idustbound,set_particle_type
 use io,   only:iprint
 real,    intent(in)    :: dtg
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:),xyzh(:,:),vxyzu(:,:)
 integer, intent(out)   :: ierr
 integer :: i,j

 if (maxphase /= maxp) then
    print "(a)",' ERROR: cannot set dust particles (iphase not stored)'
    ierr = -1
    return
 endif
 ierr = 0
 j = npart
 do i=1,npart
    if (j+1 > maxp) then
       print*,' error: memory allocation too small for dust particles'
       npartoftype(idust) = 0
       ierr = 1
       return
    endif
    select case(iamtype(iphase(i)))
    case(igas)
       npartoftype(idust) = npartoftype(idust) + 1
       j = j + 1
       xyzh(:,j) = xyzh(:,i)
       vxyzu(1:3,j) = vxyzu(1:3,i)
       call set_particle_type(j,idust)
    case(iboundary)
       npartoftype(idustbound) = npartoftype(idustbound) + 1
       j = j + 1
       xyzh(:,j) = xyzh(:,i)
       vxyzu(1:3,j) = vxyzu(1:3,i)
       call set_particle_type(j,idustbound)
    end select
 enddo
 massoftype(idust) = dtg*massoftype(igas)
 massoftype(idustbound) = massoftype(idust)
 write(iprint,'(1x,a,i8)') 'Setup_shock: ndust     = ',npartoftype(idust)
 write(iprint,'(1x,a,i8)') 'Setup_shock: ndustbound= ',npartoftype(idustbound)
 write(iprint,'(1x,a,es16.8,/)') 'Setup_shock: mass of dust & dust boundary parts = ', massoftype(idust)

 npart = npart + npartoftype(idust) + npartoftype(idustbound)

end subroutine set_dust_particles


!-----------------------------------------------------------------------
!+
!  Choose which shock tube problem to set up
!+
!-----------------------------------------------------------------------
subroutine choose_shock (gamma,polyk,dtg,iexist)
 use io,        only:fatal,id,master
 use dim,       only:mhd,maxvxyzu,use_dust,do_radiation,mhd_nonideal,gr
 use eos,       only:equationofstate,ieos
 use physcon,   only:pi,Rg,au,solarm
 use options,   only:nfulldump,alpha,alphamax,alphaB,use_dustfrac
 use options,   only:alphau
 use timestep,  only:dtmax,tmax
 use prompting, only:prompt
 use dust,      only:K_code
#ifdef NONIDEALMHD
 use nicil,       only:use_ohm,use_hall,use_ambi,eta_constant,eta_const_type, &
                       C_OR,C_HE,C_AD,C_nimhd,icnstphys,icnstsemi,icnst
#endif
 use units,     only:udist,utime,unit_density,unit_pressure
 use eos,       only:gmw
 real,    intent(inout) :: gamma,polyk
 real,    intent(out)   :: dtg
 logical, intent(in)    :: iexist
 integer, parameter     :: nshocks = 11
 character(len=30)      :: shocks(nshocks)
 integer                :: i,choice
#ifdef NONIDEALMHD
 real                   :: gamma_AD,rho_i_cnst
#endif
 real                   :: const,uu,dens,pres,Tgas
 integer                :: relativistic_choice
 real                   :: uthermconst,densleft,densright,pondens,spsound,soundspeed
!
!--set default file output parameters
!
 if (.not. iexist) then
    tmax       = 0.20
    dtmax      = 0.01
    nfulldump  = 1
    alpha      = 1.0
    alphamax   = 1.0
    alphaB     = 1.0
#ifdef NONIDEALMHD
    use_ohm    = .false.
    use_hall   = .false.
    use_ambi   = .false.
    eta_constant = .true.
    eta_const_type = icnstsemi
#endif
 endif
 nx     = 256
 xleft  = -0.500
 yleft  = 0.0
 zleft  = 0.0
 xright = 0.0
 yright = 0.0
 zright = 0.0
 const  = sqrt(4.*pi)
!
!--list of shocks
!
 shocks(:) = 'none'
 shocks(1) = 'Sod shock'
 shocks(2) = 'Ryu 1a'
 shocks(3) = 'Ryu 1b'
 shocks(4) = 'Ryu 2a (w 7 discontinuities)'
 shocks(5) = 'Ryu 2b'
 shocks(6) = 'Brio-Wu (Ryu 5a)'
 shocks(7) = 'C-shock'
 shocks(8) = 'Steady shock'
 shocks(9) = 'Radiation shock'
 shocks(10) = 'Relativistic Sod shock'

 do i = 1, nshocks
    if (trim(shocks(i)) /= 'none') write(*,"(a5,i2,1x,a30)") 'Case ', i, shocks(i)
 enddo

 choice = 1
 if (mhd) then
    if (mhd_nonideal) then
       choice = 7
    else
       choice = 6
    endif
 endif
 if (gr) then
    choice = 10
 endif
 call prompt('Enter shock choice',choice,1,nshocks)
 icase = choice

 if (id==master) write(*,"('Setting up ',a)") trim(shocks(choice))
 select case (choice)
 case(1)
    !--Sod shock
    shocktype = 'Sod shock'
    gamma      = 5./3.
    leftstate  = (/1.000,1.0,0.,0.,0.,0.,0.,0./)
    rightstate = (/0.125,0.1,0.,0.,0.,0.,0.,0./)
    if (maxvxyzu < 4) call fatal('setup','Sod shock tube requires ISOTHERMAL=no')
 case(2)
    !--Ryu et al. shock 1a
    shocktype = 'Ryu et al. shock 1a'
    nx          = 128
    if (.not. iexist) then
       tmax      =   0.08
       dtmax     =   0.004
    endif
    gamma      =   5./3.
    leftstate  = (/1.,20.,10.,0.,0.,5./const,5./const,0./)
    rightstate = (/1.,1.,-10.,0.,0.,5./const,5./const,0./)
 case(3)
    !--Ryu et al. shock 1b
    shocktype = 'Ryu et al. shock 1b'
    if (.not. iexist) then
       tmax      =   0.03
       dtmax     =   0.0015
    endif
    gamma      =  5./3.
    leftstate  = (/1.0,1. ,0.,0.,0.,5./const,5./const,0./)
    rightstate = (/0.1,10.,0.,0.,0.,5./const,2./const,0./)
 case(4)
    !--Ryu et al. shock 2a
    shocktype  = "Ryu et al. shock 2a (with 7 discontinuities)"
    gamma      = 5./3.
    leftstate  = (/1.08,0.95,1.2,0.01,0.5,2./const,3.6/const,2./const/)
    rightstate = (/1.  ,1.  ,0. ,0.  ,0. ,2./const,4.0/const,2./const/)
 case(5)
    !--Ryu et al. shock 2b
    shocktype = 'Ryu et al. shock 2b'
    if (.not. iexist) then
       tmax     =   0.035
       dtmax    =   0.00175
    endif
    gamma      = 5./3.
    leftstate  = (/1.0,1. ,0.,0.,0.,3./const,6./const,0./)
    rightstate = (/0.1,10.,0.,2.,1.,3./const,1./const,0./)
 case(6)
    !--Brio-Wu shock
    shocktype = 'Brio/Wu (Ryu/Jones shock 5a)'
    if (.not. iexist) then
       tmax    = 0.1
       dtmax   = 0.005
    endif
    gamma      = 2.0
    leftstate  = (/1.000,1.0,0.,0.,0.,0.75, 1.,0./)
    rightstate = (/0.125,0.1,0.,0.,0.,0.75,-1.,0./)
 case(7)
    !--C-shock
    shocktype = 'C-shock'
    if (.not. iexist) then
       tmax       = 4.0e6
       dtmax      = 1.0e4
       alpha      = 0.0
#ifdef NONIDEALMHD
       use_ambi   = .true.
       gamma_AD   = 1.0
       rho_i_cnst = 1.0d-5
       C_AD       = 1.0/(gamma_AD*rho_i_cnst)
       C_nimhd    = 0.25/pi
#endif
    endif
    gamma      =  1.0
    polyk      =  0.01
    leftstate  = (/1.,0.006, 4.45,0.,0.,1./sqrt(2.),1./sqrt(2.),0./)
    rightstate = (/1.,0.006,-4.45,0.,0.,1./sqrt(2.),1./sqrt(2.),0./)
    nx         = 200
    xleft      = -4.00d6
 case(8)
    !--Steady shock (Falle 2003)
#ifdef NONIDEALMHD
    shocktype = 'Steady shock with large Hall Effect (Falle 2003; fig 3)'
    if (.not. iexist) then
       use_ohm      = .true.
       use_hall     = .true.
       use_ambi     = .true.
       C_OR         =  1.12d-9
       C_HE         = -3.53d-2
       C_AD         =  7.83d-3
    endif
#else
    shocktype = 'Steady shock (Falle 2003)'
#endif
    if (.not. iexist) then
       tmax     = 1.0
       alpha    = 0.0
       alphamax = 1.0
       alphaB   = 1.0
    endif
    nx         = 512
    polyk      = 0.01
    gamma      = 1.0
    leftstate  = (/1.7942,0.017942,-0.9759,-0.6561,0.,1.,1.74885,0./)
    rightstate = (/1.    ,0.01    ,-1.7510, 0.    ,0.,1.,0.6    ,0./)
    xleft      = -2.0
 case(9)
    ! if (.not.do_radiation) call fatal('setup','Radiation shock is only possible with "RADIATION=yes"')
    shocktype = 'Radiation shock'
    gamma = 5./3.
    gmw   = 2.38
    Tgas  = 1500.
    uu    = Tgas*Rg/(gamma - 1.0)/gmw
    dens  = 1.e-10
    pres  = (gamma-1.)*uu*dens/unit_pressure
    dens  = dens/unit_density
    ! (/'dens','pr  ','vx  ','vy  ','vz  ','Bx  ','By  ','Bz  '/)
    leftstate  = (/dens, pres,  3.2e5/(udist/utime), 0.,0.,0.,0.,0./)
    rightstate = (/dens, pres, -3.2e5/(udist/utime), 0.,0.,0.,0.,0./)
    tmax   = 1e9/utime
    dtmax  = 1e7/utime
    kappa  = 4e1
    xright = -3.2e5
    xleft  =  3.2e5
    call prompt('velocity right',xright, -1e30, 0.)
    call prompt('velocity left ',xleft, 0.,   1e30)
    leftstate(ivx)  = xleft/(udist/utime)
    rightstate(ivx) = xright/(udist/utime)
    xright =  1e15
    xleft  = -1e15
    call prompt('border right',xright,0.,  1e30)
    call prompt('border left',xleft,  -1e30, 0.)
    xright = xright/udist
    xleft  = xleft/udist
 case(10)
    !--Sod shock
    relativistic_choice = 1
    shocktype = "Mildly-Relativistic Sod shock"
    gamma      = 5./3.
    alphau     = 0.1
    leftstate(1:iBz)  = (/10.0,40./3.,0.,0.,0.,0.,0.,0./)
    rightstate(1:iBz) = (/1.00,1.e-6 ,0.,0.,0.,0.,0.,0./)
    write(*,"(a5,i2,1x,a20)") 'Case ', 1, 'Mildly relativistic'
    write(*,"(a5,i2,1x,a20)") 'Case ', 2, 'Ultra relativistic'
    write(*,"(a5,i2,1x,a20)") 'Case ', 3, 'Isothermal'
    call prompt('Enter relativistic shock choice',relativistic_choice,1,3)
    select case(relativistic_choice)
    case(2)
       shocktype = "Ultra-Relativistic Sod shock"
       leftstate(1:iBz)  = (/1.,1000.,0.,0.,0.,0.,0.,0./)
       rightstate(1:iBz) = (/1.,0.01 ,0.,0.,0.,0.,0.,0./)
    case(3)
       shocktype = "Isothermal relativistic shock"
       ieos        = 4
       soundspeed  = 0.1
       call prompt('Enter sound speed',soundspeed,0.,1.)
       uthermconst = soundspeed**2/(gamma-1.-soundspeed**2)
       polyk       = uthermconst
       densleft    = 10.
       densright   = 1.
       call equationofstate(ieos,pondens,spsound,densleft,0.,0.,0.)
       if (abs(spsound/soundspeed)-1.>1.e-10) call fatal('setup','eos soundspeed does not match chosen sound speed')
       leftstate(1:iBz)  = (/densleft,pondens*densleft,0.,0.,0.,0.,0.,0./)
       call equationofstate(ieos,pondens,spsound,densright,0.,0.,0.)
       rightstate(1:iBz) = (/densright,pondens*densright,0.,0.,0.,0.,0.,0./)
       if (abs(spsound/soundspeed)-1.>1.e-10) call fatal('setup','eos soundspeed does not match chosen sound speed')
    case default
    end select
    if (maxvxyzu < 4) call fatal('setup','Sod shock tube requires ISOTHERMAL=no')
 end select

 call prompt('Enter resolution (number of particles in x) for left half (x<0)',nx,8)

 if (abs(xright)  < epsilon(xright))  xright  = -xleft

 if (use_dust) then
    !--shock setup supports both one-fluid and two-fluid dust
    dtg = 1.
    K_code = 1000.
    call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
    use_dustfrac = (dust_method == 1)
    call prompt('Enter dust to gas ratio',dtg,0.)
    call prompt('Enter constant drag coefficient',K_code(1),0.)
 endif

 if (do_radiation) then
    call prompt('Enter kappa (total radiation opacity)',kappa,0.,1e6)
 endif

 return
end subroutine choose_shock

!------------------------------------------
!+
!  Pretty-print the shock parameters
!+
!------------------------------------------
subroutine print_shock_params(nstates)
 integer, intent(in) :: nstates
 integer             :: i

 write(*,"(/,1x,'Setup_shock: ',a,/,8(11x,a4,' L: ',f8.3,' R:',f8.3,/))") &
    trim(shocktype),(trim(var_label(i)),leftstate(i),rightstate(i),i=1,nstates)

end subroutine print_shock_params

!------------------------------------------
!+
!  Function to return conserved density
!+
!------------------------------------------
real function get_conserved_density(state) result(rho)
 use dim, only:gr
 real, intent(in) :: state(max_states)
 real :: lorentz,v2

 if (gr) then
    v2 = dot_product(state(ivx:ivz),state(ivx:ivz))
    lorentz = 1./sqrt(1.-v2)
    rho = lorentz*state(idens)
 else
    rho = state(idens)
 endif

end function get_conserved_density

!------------------------------------------
!+
!  Write setup parameters to input file
!+
!------------------------------------------
subroutine write_setupfile(filename,iprint,numstates,gamma,polyk,dtg)
 use infile_utils, only:write_inopt
 use dim,          only:tagline,maxvxyzu,use_dust,do_radiation
 integer,          intent(in) :: iprint,numstates
 real,             intent(in) :: gamma,polyk,dtg
 character(len=*), intent(in) :: filename
 integer, parameter           :: lu = 20
 integer                      :: i,ierr1,ierr2

 write(iprint,"(a)") ' Writing '//trim(filename)//' with initial left/right states'
 open(unit=lu,file=filename,status='replace',form='formatted')
 write(lu,"(a)") '# '//trim(tagline)
 write(lu,"(a)") '# input file for Phantom shock tube setup'

 write(lu,"(/,a)") '# shock tube name'
 call write_inopt(trim(shocktype),'name','',lu,ierr1)

 write(lu,"(/,a)") '# shock tube'
 do i=1,numstates
    call write_inopt(leftstate(i), trim(var_label(i))//'left', trim(var_label(i))//' (left)', lu,ierr1)
    call write_inopt(rightstate(i),trim(var_label(i))//'right',trim(var_label(i))//' (right)',lu,ierr2)
    if (ierr1 /= 0 .or. ierr2 /= 0) write(*,*) 'ERROR writing '//trim(var_label(i))
 enddo

 write(lu,"(/,a)") '# boundaries'
 call write_inopt(xleft,'xleft','x min boundary',lu,ierr1)
 call write_inopt(xright,'xright','x max boundary',lu,ierr2)
 if (ierr1 /= 0 .or. ierr2 /= 0) write(*,*) 'ERROR writing xmin, xmax'

 write(lu,"(/,a)") '# resolution'
 call write_inopt(nx,'nx','resolution (number of particles in x) for -xleft < x < xshock',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing nx'

 write(lu,"(/,a)") '# Equation-of-state properties'
 call write_inopt(gamma,'gamma','Adiabatic index',lu,ierr1)
 if (maxvxyzu==3) then
    call write_inopt(polyk,'polyk','square of the isothermal sound speed',lu,ierr1)
 endif
 if (ierr1 /= 0 .or. ierr2 /= 0) write(*,*) 'ERROR writing gamma, polyk'

 if (use_dust) then
    write(lu,"(/,a)") '# dust properties'
    call write_inopt(dust_method,'dust_method','1=one fluid, 2=two fluid',lu,ierr1)
    call write_inopt(dtg,'dtg','Dust to gas ratio',lu,ierr2)
    if (ierr1 /= 0 .or. ierr2 /= 0) write(*,*) 'ERROR writing dust options'
 endif

 if (do_radiation) then
    write(lu,"(/,a)") '# radiation properties'
    call write_inopt(kappa,'kappa','opacity in cm^2/g',lu,ierr1)
 endif

 close(unit=lu)

end subroutine write_setupfile

!------------------------------------------
!+
!  Read setup parameters from input file
!+
!------------------------------------------
subroutine read_setupfile(filename,iprint,numstates,gamma,polyk,dtg,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db,read_inopt
 use dim,          only:maxvxyzu,use_dust,do_radiation
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(in)  :: iprint,numstates
 integer,          intent(out) :: ierr
 real,             intent(out) :: gamma,polyk,dtg
 integer                       :: i,nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(iprint, '(1x,2a)') 'Setup_shock: Reading setup options from ',trim(filename)

 nerr = 0
 shocktype = "Name not read from .setup"
 call read_inopt(shocktype,'name',db)
 do i=1,numstates
    call read_inopt(leftstate(i), trim(var_label(i))//'left',db,errcount=nerr)
    call read_inopt(rightstate(i),trim(var_label(i))//'right',db,errcount=nerr)
 enddo
 call read_inopt(xleft,'xleft',db,errcount=nerr)
 call read_inopt(xright,'xright',db,min=xleft,errcount=nerr)
 call read_inopt(nx,'nx',db,min=1,errcount=nerr)

 call read_inopt(gamma,'gamma',db,min=1.,errcount=nerr)
 if (maxvxyzu==3) call read_inopt(polyk,'polyk',db,min=0.,errcount=nerr)

 if (use_dust) then
    call read_inopt(dust_method,'dust_method',db,min=1,errcount=nerr)
    call read_inopt(dtg,'dtg',db,min=0.,errcount=nerr)
 endif

 if (do_radiation) then
    call read_inopt(kappa,'kappa',db,min=0.,errcount=nerr)
 endif

 if (nerr > 0) then
    print "(1x,a,i2,a)",'Setup_shock: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile

end module setup
