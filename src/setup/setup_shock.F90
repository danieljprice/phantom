!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
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
!    gamma  -- Adiabatic index
!    nx     -- resolution (number of particles in x) for -xleft < x < xshock
!    polyk  -- square of the isothermal sound speed
!    xleft  -- x min boundary
!    xright -- x max boundary
!
!  DEPENDENCIES: boundary, dim, infile_utils, io, kernel, mpiutils, nicil,
!    options, part, physcon, prompting, set_dust, setup_params, timestep,
!    unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 !
 integer :: nx, icase
 real    :: xleft, xright, yleft, yright, zleft, zright
 real    :: dxleft
 character(len=100) :: shocktype
 character(len=100) :: latticetype = 'closepacked'
 logical :: use_closepacked

 integer, parameter :: max_states = 8
 integer            :: nstates
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
 !
 public  :: setpart
 !
 private
 !
contains

!----------------------------------------------------------------
!+
!  setup for shock tube problems in 3D
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:fatal,master,iprint
 use unifdis,      only:set_unifdis,get_ny_nz_closepacked
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,set_boundary
 use mpiutils,     only:bcast_mpi
 use dim,          only:maxp,maxvxyzu,ndim,mhd
 use options,      only:use_dustfrac
 use part,         only:labeltype,set_particle_type,igas,iboundary,hrho,Bxyz,mhd,periodic,dustfrac
 use kernel,       only:radkern,hfact_default
 use timestep,     only:tmax
 use prompting,    only:prompt
 use set_dust,     only:set_dustfrac_from_inopts
#ifdef NONIDEALMHD
 use nicil,          only:rho_i_cnst
#endif
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npartoftype(:)
 integer,           intent(inout) :: npart
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real                             :: totmass
 real                             :: xminleft(ndim),xmaxleft(ndim),xminright(ndim),xmaxright(ndim)
 real                             :: delta,gam1,xshock,fac,dtg
 real                             :: uuleft,uuright,volume,xbdyleft,xbdyright,dxright
 integer                          :: i,ierr,nbpts,ny,nz
 character(len=120)               :: shkfile, filename
 logical                          :: iexist
 !
 ! quit if not periodic
 !
 if (.not.periodic) call fatal('setup','require PERIODIC=yes')
 !
 ! verify a legitimate lattice type has been chosen
 !
 if (trim(latticetype)/='random' .and. trim(latticetype)/='cubic' .and. &
     trim(latticetype)/='closepacked' .and. trim(latticetype)/='hcp') then
    call fatal('setup','invalid lattice type')
 endif
 if (trim(latticetype)=='closepacked') then
    use_closepacked = .true.
 else
    use_closepacked = .false.
 endif
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
 if (.not.iexist) hfact = hfact_default
 nstates = max_states
 if (.not.mhd) nstates = max_states - 3
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
 dxleft = -xleft/float(nx)
 xshock = 0.5*(xleft + xright)
 !
 ! adjust boundaries to allow space for boundary particles and inflow
 !
 call adjust_shock_boundaries(dxleft,dxright,radkern, &
      leftstate(ivx),rightstate(ivx),leftstate(idens),rightstate(idens),tmax,ndim)
 !
 ! print setup parameters
 !
 if (id==master) call print_shock_params(nstates)
 !
 ! set domain boundaries - use a longer box in x direction since
 ! it is not actually periodic in x
 !
 call set_boundary(xleft-1000.*dxleft,xright+1000.*dxright,yleft,yright,zleft,zright)
 !
 ! set limits of the different domains
 !
 xminleft(:)  = (/xleft,ymin,zmin/)
 xmaxleft(:)  = (/xright,ymax,zmax/)
 xminright(:) = (/xleft,ymin,zmin/)
 xmaxright(:) = (/xright,ymax,zmax/)
 if ( (ymax-ymin) > (xmax-xmin) ) call fatal('setup','(ymax-ymin) > (xmax-xmin), but shock is in x-direction')
 if ( (zmax-zmin) > (xmax-xmin) ) call fatal('setup','(zmax-zmin) > (xmax-xmin), but shock is in x-direction')
 !
 ! setup the particles
 !
 if (abs(leftstate(idens)-rightstate(idens)) > epsilon(0.)) then
    ! then divide the x axis into two halves at xshock
    xmaxleft(1)  = xshock
    xminright(1) = xshock
    write(iprint,'(1x,3(a,es16.8))') 'Setup_shock: left half  ',xminleft(1), ' to ',xmaxleft(1), ' with dx_left  = ',dxleft

    ! set up a uniform lattice
    call set_unifdis(latticetype,id,master,xminleft(1),xmaxleft(1),xminleft(2), &
                     xmaxleft(2),xminleft(3),xmaxleft(3),dxleft,hfact,npart,xyzh)  ! set left half

    ! set particle mass
    volume           = product(xmaxleft-xminleft)
    totmass          = volume*leftstate(idens)*(1. + dtg)
    massoftype(igas) = totmass/npart
    if (id==master) print*,' particle mass = ',massoftype(igas)

    if (use_closepacked) then
       ! now adjust spacing on right hand side to get correct density given the particle mass
       volume  = product(xmaxright - xminright)
       totmass = volume*rightstate(idens)*(1. + dtg)
       call get_ny_nz_closepacked(dxright,xminright(2),xmaxright(2),xminright(3),xmaxright(3),ny,nz)
       dxright = (xmaxright(1) - xminright(1))/((totmass/massoftype(igas))/(ny*nz))
    endif

    ! now set up box for right half
    write(iprint,'(1x,3(a,es16.8))') 'Setup_shock: right half ',xminright(1),' to ',xmaxright(1),' with dx_right = ',dxright

    call set_unifdis(latticetype,id,master,xminright(1),xmaxright(1), &
         xminright(2),xmaxright(2),xminright(3),xmaxright(3),dxright,hfact,npart,xyzh,npy=ny,npz=nz) ! set right half

    ! define rhozero as average density; required for certain simulations (e.g. non-ideal MHD with constant resistivity)
    rhozero = (leftstate(idens)*product(xmaxleft-xminleft) + rightstate(idens)*product(xmaxright - xminright)) &
               / product(xmaxright - xminleft)
 else  ! set all of volume if densities are equal
    write(iprint,'(3(a,es16.8))') 'Setup_shock: one density  ',xminleft(1), ' to ',xmaxright(1), ' with dx  = ',dxleft
    call set_unifdis(latticetype,id,master,xminleft(1),xmaxleft(1),xminleft(2), &
                     xmaxleft(2),xminleft(3),xmaxleft(3),dxleft,hfact,npart,xyzh)
    volume           = product(xmaxleft-xminleft)
    rhozero          = leftstate(idens)
    dxright          = dxleft
    massoftype(igas) = leftstate(idens)*volume/real(npart)
 endif
 !
 ! Fix the particles near x-boundary; else define as gas
 !
 nbpts     = 0
 fac       = nint(2.01*radkern*hfact)
 xbdyleft  = fac*dxleft
 xbdyright = fac*dxright
 dustfrac  = 0.
 do i=1,npart
    if ( (xyzh(1,i) < (xminleft (1) + xbdyleft ) ) .or. &
         (xyzh(1,i) > (xmaxright(1) - xbdyright) ) ) then
       call set_particle_type(i,iboundary)
       nbpts = nbpts + 1
    else
       call set_particle_type(i,igas)
    endif
    !
    ! one fluid dust: set dust fraction on gas particles
    if (use_dustfrac) call set_dustfrac_from_inopts(dtg,ipart=i)
 enddo
 massoftype(iboundary)  = massoftype(igas)
 npartoftype(iboundary) = nbpts
 npartoftype(igas)      = npart - nbpts
 write(iprint,'(1x,a,i8)') 'Setup_shock: npart     = ',npart
 write(iprint,'(1x,a,i8)') 'Setup_shock: ngas      = ',npartoftype(igas)
 write(iprint,'(1x,a,i8)') 'Setup_shock: nboundary = ',nbpts
 write(iprint,'(1x,a,es16.8,/)') 'Setup_shock: mass of gas & boundary particles = ', massoftype(igas)
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
 do i=1,npart
    delta = xyzh(1,i) - xshock
    if (delta > 0.) then
       xyzh(4,i)  = hrho(rightstate(idens),massoftype(igas))
       vxyzu(1,i) = rightstate(ivx)
       vxyzu(2,i) = rightstate(ivy)
       vxyzu(3,i) = rightstate(ivz)
       if (maxvxyzu >= 4) vxyzu(4,i) = uuright
       if (mhd) Bxyz(1:3,i) = rightstate(iBx:iBz)
    else
       xyzh(4,i)  = hrho(leftstate(idens),massoftype(igas))
       vxyzu(1,i) = leftstate(ivx)
       vxyzu(2,i) = leftstate(ivy)
       vxyzu(3,i) = leftstate(ivz)
       if (maxvxyzu >= 4) vxyzu(4,i) = uuleft
       if (mhd) Bxyz(1:3,i) = leftstate(iBx:iBz)
    endif
 enddo
 if (mhd) ihavesetupB = .true.

#ifdef NONIDEALMHD
 !Modifiy ion density from fraction to physical (for ambipolar diffusion)
 if (.not.iexist) rho_i_cnst = rho_i_cnst * rhozero
#endif

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Adjust the shock boundaries to allow for inflow/outflow
!+
!-----------------------------------------------------------------------
subroutine adjust_shock_boundaries(dxleft,dxright,radkern,vxleft,vxright, &
                                   densleft,densright,tmax,ndim)
 real,    intent(in)    :: dxleft,radkern,vxleft,vxright,densleft,densright,tmax
 real,    intent(out)   :: dxright
 integer, intent(in)    :: ndim
 real :: fac

 if (vxleft > tiny(vxleft)) then
    xleft = xleft - vxleft*tmax
 endif
 dxright = dxleft*(densleft/densright)**(1./ndim) ! NB: dxright here is only approximate
 if (vxright < -tiny(vxright)) then
    xright = xright - vxright*tmax
 endif
 ! try to give y boundary that is a multiple of 6 particle spacings in the low density part
 fac = -6.*(int(1.99*radkern/6.) + 1)*max(dxleft,dxright)
 if (use_closepacked) then
    yleft   = fac*sqrt(0.75)
    zleft   = fac*sqrt(6.)/3.
 else
    yleft   = fac
    zleft   = fac
 endif
 yright  = -yleft
 zright  = -zleft

end subroutine adjust_shock_boundaries

!-----------------------------------------------------------------------
!+
!  Choose which shock tube problem to set up
!+
!-----------------------------------------------------------------------
subroutine choose_shock (gamma,polyk,dtg,iexist)
 use io,        only:fatal,id,master
 use dim,       only:mhd,maxvxyzu,use_dust
 use physcon,   only:pi
 use options,   only:nfulldump,alpha,alphamax,alphaB
 use timestep,  only:dtmax,tmax
 use prompting, only:prompt
 use set_dust,  only:interactively_set_dust
#ifdef NONIDEALMHD
 use nicil,       only:use_ohm,use_hall,use_ambi,eta_constant,eta_const_type, &
                       C_OR,C_HE,C_AD,C_nimhd,icnstphys,icnstsemi,icnst
#endif
 real,    intent(inout) :: gamma,polyk
 real,    intent(out)   :: dtg
 logical, intent(in)    :: iexist
 integer, parameter     :: nshocks = 10
 character(len=30)      :: shocks(nshocks)
 integer                :: i, choice,dust_method
 real                   :: const !, dxright
#ifdef NONIDEALMHD
 real                   :: gamma_AD,rho_i_cnst
#endif
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

 do i = 1, nshocks
    if (trim(shocks(i)) /= 'none') write(*,"(a5,i2,1x,a30)") 'Case ', i, shocks(i)
 enddo

 choice = 1
#ifdef MHD
#ifdef NONIDEALMHD
 choice = 7
#else
 choice = 6
#endif
#endif
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
    gamma      =  1.0
    leftstate  = (/1.7942,0.017942,-0.9759,-0.6561,0.,1.,1.74885,0./)
    rightstate = (/1.    ,0.01    ,-1.7510, 0.    ,0.,1.,0.6    ,0./)
    xleft      = -2.0
 end select

 call prompt('Enter resolution (number of particles in x) for left half (x<0)',nx,8)

 if (abs(xright)  < epsilon(xright))  xright  = -xleft

 if (use_dust) then
    dust_method  = 1
    call interactively_set_dust(dtg,imethod=dust_method,Kdrag=.true.)
    if (dust_method == 2) call fatal('setup','shock setup currently does not support two-fluid dust')
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
!  Write setup parameters to input file
!+
!------------------------------------------
subroutine write_setupfile(filename,iprint,numstates,gamma,polyk,dtg)
 use infile_utils, only:write_inopt
 use dim,          only:tagline,maxvxyzu
 use options,      only:use_dustfrac
 use set_dust,     only:write_dust_setup_options
 integer,          intent(in) :: iprint,numstates
 real,             intent(in) :: gamma,polyk,dtg
 character(len=*), intent(in) :: filename
 integer, parameter           :: lu = 20
 integer                      :: i,ierr1,ierr2,dust_method

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

 if (use_dustfrac) then
    !write(lu,"(/,a)") '# dust properties'
    !call write_inopt(dtg,'dtg','Dust to gas ratio',lu,ierr1)
    !if (ierr1 /= 0) write(*,*) 'ERROR writing dtg'

    dust_method = 1  !--one-fluid is currently forced for dustyshock
    call write_dust_setup_options(lu,dtg,imethod=dust_method,isimple=.true.)
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
 use dim,          only:maxvxyzu,use_dust
 use set_dust,     only:read_dust_setup_options
 use io,           only:fatal
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(in)  :: iprint,numstates
 integer,          intent(out) :: ierr
 real,             intent(out) :: gamma,polyk,dtg
 integer                       :: i,nerr,dust_method
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
 call read_inopt(nx,'nx',db,min=8,errcount=nerr)

 call read_inopt(gamma,'gamma',db,min=1.,errcount=nerr)
 if (maxvxyzu==3) call read_inopt(polyk,'polyk',db,min=0.,errcount=nerr)

 if (use_dust) then
    call read_dust_setup_options(db,nerr,dtg,imethod=dust_method,isimple=.true.)
    if (dust_method == 2) call fatal('setup','shock setup currently does not support two-fluid dust')
 endif

 if (nerr > 0) then
    print "(1x,a,i2,a)",'Setup_shock: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile
!-----------------------------------------------------------------------
end module setup
