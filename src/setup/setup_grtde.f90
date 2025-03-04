!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for general relativistic tidal disruption event
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters:
!   - beta           : *penetration factor*
!   - dumpsperorbit  : *number of dumps per orbit*
!   - ecc_bh         : *eccentricity (1 for parabolic)*
!   - mhole          : *mass of black hole (solar mass)*
!   - norbits        : *number of orbits*
!   - provide_params : *manually specify the position and velocity of the star(s)*
!   - sep_initial    : *initial separation from BH in tidal radii*
!   - theta_bh       : *inclination of orbit (degrees)*
!   - vx1            : *vel x star 1*
!   - vx2            : *vel x star 2*
!   - vy1            : *vel y star 1*
!   - vy2            : *vel y star 2*
!   - vz1            : *vel z star 1*
!   - vz2            : *vel z star 2*
!   - x1             : *pos x star 1*
!   - x2             : *pos x star 2*
!   - y1             : *pos y star 1*
!   - y2             : *pos y star 2*
!   - z1             : *pos z star 1*
!   - z2             : *pos z star 2*
!
! :Dependencies: eos, externalforces, gravwaveutils, infile_utils, io,
!   kernel, metric, mpidomain, options, part, physcon, relaxstar,
!   setbinary, setorbit, setstar, setup_params, systemutils, timestep,
!   units, vectorutils
!

 use setstar,  only:star_t
 use setorbit, only:orbit_t
 use metric,   only:mass1,a
 implicit none
 public :: setpart

 real    :: mhole,beta,ecc_bh,norbits,theta_bh,sep_initial
 real    :: x1,y1,z1,x2,y2,z2
 real    :: vx1,vy1,vz1,vx2,vy2,vz2
 integer :: dumpsperorbit,nstar
 logical :: relax,write_profile
 logical :: provide_params
 logical :: use_gr_ic
 integer, parameter :: max_stars = 2
 type(star_t)  :: star(max_stars)
 type(orbit_t) :: orbit

 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,&
                     gravity,eos_vars,rad,gr,nsinkproperties
 use setbinary, only:set_binary
 use setorbit,   only:set_defaults_orbit,set_orbit
 use setstar,   only:shift_star,set_defaults_stars,set_stars,shift_stars
 use units,     only:set_units
 use physcon,   only:solarm,pi,solarr
 use io,        only:master,fatal,warning
 use timestep,  only:tmax,dtmax
 use eos,       only:ieos,X_in,Z_in
 use kernel,    only:hfact_default
 use mpidomain, only:i_belong
 use externalforces, only:accradius1,accradius1_hard
 use vectorutils,    only:rotatevec
 use gravwaveutils,  only:theta_gw,calc_gravitwaves
 use setup_params,   only:rhozero,npart_total
 use systemutils,    only:get_command_option
 use options,        only:iexternalforce
 use units,          only:in_code_units
 use, intrinsic                   :: ieee_arithmetic
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename
 integer :: ierr,np_default
 integer :: nptmass_in
 integer :: i
 logical :: iexist,use_var_comp
 real    :: rtidal,rp,semia,period,hacc1,hacc2
 real    :: vxyzstar(3),xyzstar(3)
 real    :: r0,vel,lorentz
 real    :: vhat(3),x0,y0
 real    :: semi_maj_val
 real    :: mstars(max_stars),rstars(max_stars),haccs(max_stars)
 real    :: xyzmh_ptmass_in(nsinkproperties,2),vxyz_ptmass_in(3,2),angle
 real    :: alpha, delta_v, epsilon_target, tol
 integer :: max_iters
!
!-- general parameters
!
 hfact = hfact_default
 time  = 0.
 polyk = 1.e-10    ! <== uconst
 gamma = 5./3.
 angle = 0.
 xyzmh_ptmass_in(:,:) = 0.
 vxyz_ptmass_in(:,:)  = 0.
 if (.not.gravity) call fatal('setup','recompile with GRAVITY=yes')
!
!-- space available for injected gas particles
!
 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.
 nptmass        = 0
 nstar          = 1
!
!-- Default runtime parameters
!
 mhole           = 1.e6  ! (solar masses)
 call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
 call set_defaults_stars(star)
 call set_defaults_orbit(orbit)
 star(:)%m       = '1.*msun'
 star(:)%r       = '1.*solarr'
 np_default      = 1e6
 star%np         = int(get_command_option('np',default=np_default)) ! can set default value with --np=1e5 flag (mainly for testsuite)
!  star%iprofile   = 2
 beta            = 5.
 ecc_bh          = 0.8
 norbits         = 5.
 dumpsperorbit   = 100
 theta_bh        = 0.
 sep_initial     = 10.  ! initial separation in units of tidal radius
 write_profile   = .true.
 use_var_comp    = .false.
 relax           = .true.
 use_gr_ic       = .true. ! Whether initial velocity condition computed in GR is used for parabolic orbits.

!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Tidal disruption in GR'
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
 endif
 !
 !--set up and relax the stellar profiles for one or both stars
 !
 call set_stars(id,master,nstar,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,&
                X_in,Z_in,relax,use_var_comp,write_profile,&
                rhozero,npart_total,i_belong,ierr)

 do i=1,nstar
    rstars(i) = in_code_units(star(i)%r,ierr,unit_type='length')
    if (ierr /= 0) call fatal('setup','could not convert rstar to code units',i=i)
    mstars(i) = in_code_units(star(i)%m,ierr,unit_type='mass')
    if (ierr /= 0) call fatal('setup','could not convert mstar to code units',i=i)
    haccs(i) = in_code_units(star(i)%hacc,ierr,unit_type='mass')
    if (ierr /= 0) call fatal('setup','could not convert hacc to code units',i=i)
 enddo

 if (star(1)%iprofile == 0 .and. nstar == 1) then
    xyzmh_ptmass_in(4,1) = mstars(1)
    xyzmh_ptmass_in(5,1) = haccs(1)
 endif

 !
 !--set the stars around each other first if nstar > 1 (Assuming binary system)
 !
 if (nstar > 1 .and. (.not. provide_params)) then
    nptmass_in = 0
    call set_orbit(orbit,mstars(1),mstars(2),haccs(1),haccs(2),&
                   xyzmh_ptmass_in,vxyz_ptmass_in,nptmass_in,(id==master),ierr)

    if (ierr /= 0) call fatal ('setup_binary','error in call to set_orbit')
    if (ierr /= 0) call fatal('setup','errors in set_star')
 endif
 !
 !--place star / stars into orbit around the black hole
 !
 ! Calculate tidal radius
 if (nstar == 1) then
    ! for single star around the BH, the tidal radius is given by
    ! RT = rr * (MM / mm)**(1/3) where rr is rstar, MM is mass of BH and mm is mass of star
    rtidal          = rstars(1) * (mass1/mstars(1))**(1./3.)
    rp              = rtidal/beta
 else
    semi_maj_val = in_code_units(orbit%elems%semi_major_axis,ierr,unit_type='length')
    ! for a binary, tidal radius is given by
    ! orbit.an * (3 * MM / mm)**(1/3) where mm is mass of binary and orbit.an is semi-major axis of binary
    rtidal          = semi_maj_val * (3.*mass1 / (mstars(1) + mstars(2)))**(1./3.)
    rp              = rtidal/beta
 endif

 if (gr) then
    accradius1_hard = 5.*mass1
    accradius1      = accradius1_hard
 else
    if (mass1 > 0.) then
       accradius1_hard = 6.
       accradius1      = accradius1_hard
    endif
 endif

 a        = 0.
 theta_bh = theta_bh*pi/180.

 if (id==master) then
    print "(4(/,1x,a,1f10.3))",' black hole mass: ',mass1, &
                               '    tidal radius: ',rtidal,&
                               '            beta: ',beta,  &
                               'accretion radius: ',accradius1_hard
 endif

 if (.not. provide_params) then
    if (id==master) then
       do i=1,nstar
          print "(a,i1,a,es10.3)",'   mass of star ',i,': ',mstars(i)
          print "(a,i1,a,es10.3)",' radius of star ',i,': ',rstars(i)
       enddo
    endif

    xyzstar  = 0.
    vxyzstar = 0.
    period   = 0.

    if (ecc_bh<1.) then
       !
       !-- Set a binary orbit given the desired orbital parameters to get the position and velocity of the star
       !
       semia    = rp/(1.-ecc_bh)
       period   = 2.*pi*sqrt(semia**3/mass1)
       hacc1    = rstars(1)/1.e8    ! Something small so that set_binary doesnt warn about Roche lobe
       hacc2    = hacc1
       ! apocentre = rp*(1.+ecc_bh)/(1.-ecc_bh)
       ! trueanom = acos((rp*(1.+ecc_bh)/r0 - 1.)/ecc_bh)*180./pi
       call set_binary(mass1,mstars(1),semia,ecc_bh,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                       posang_ascnode=0.,arg_peri=90.,incl=0.,f=-180.)
       vxyzstar(:) = vxyz_ptmass(1:3,2)
       xyzstar(:)  = xyzmh_ptmass(1:3,2)
       nptmass  = 0

       call rotatevec(xyzstar,(/0.,1.,0./),-theta_bh)
       call rotatevec(vxyzstar,(/0.,1.,0./),-theta_bh)

    elseif (abs(ecc_bh-1.) < tiny(0.)) then
       !
       !-- Setup a parabolic orbit
       !
       r0       = sep_initial*rtidal      ! A default starting distance from the black hole.
       period   = 2.*pi*sqrt(r0**3/mass1) !period not defined for parabolic orbit, so just need some number
       y0       = -2.*rp + r0
       x0       = sqrt(r0**2 - y0**2)
       xyzstar(:)  = (/-x0,y0,0./)
       vel      = sqrt(2.*mass1/r0)
       vhat     = (/2.*rp,-x0,0./)/sqrt(4.*rp**2 + x0**2)
       vxyzstar(:) = vel*vhat
       if (rtidal <= 0.) then
          vxyzstar(:) = (/0.,0.,0./)
       elseif (use_gr_ic) then
          ! TO-DO: Obtaining initial conditions in GR for non-parabolic orbits.
          ! Parameters for gradient descent
          delta_v = 1.0d-5         ! Small velocity change for gradient computation
          alpha = 1.0d-3           ! Learning rate for gradient descent
          epsilon_target = 1.0d0   ! Target specific energy
          tol = 1.0d-9             ! Convergence tolerance
          max_iters = 1e5          ! Maximum number of iterations

          ! Perform gradient descent to obtain initial velocity in GR (epsilon=1)
          ! Assumes non-spinning SMBH (a = 0.0d0). TO-DO: Obtain velocity for SMBH spin >0.
          call refine_velocity(-x0, y0, 0.0d0, vxyzstar(1), vxyzstar(2), vxyzstar(3), &
                               mass1, 0.0d0, r0, epsilon_target, alpha, delta_v, tol, max_iters)
       endif

       call rotatevec(xyzstar,(/0.,1.,0./),theta_bh)
       call rotatevec(vxyzstar,(/0.,1.,0./),theta_bh)

    else
       call fatal('setup','please choose a valid eccentricity (0<ecc_bh<=1)',var='ecc_bh',val=ecc_bh)
    endif

    lorentz = 1./sqrt(1.-dot_product(vxyzstar(:),vxyzstar(:)))
    if (lorentz>1.1) call warning('setup','Lorentz factor of star greater than 1.1, density may not be correct')

    tmax      = norbits*period
    dtmax     = period/dumpsperorbit
 endif

 if (id==master) then
    print "(/,a)",       ' STAR SETUP:'
    print "(a,3f10.3)"  ,'         Position = ',xyzstar
    print "(a,3f10.3)"  ,'         Velocity = ',vxyzstar
    print "(a,1f10.3)"  ,' Lorentz factor   = ',lorentz
    print "(a,1f10.3)"  ,' Polytropic gamma = ',gamma
    print "(a,3f10.3,/)",'       Pericentre = ',rp
 endif
 !
 !--shift stars / sink particles
 !
 if (provide_params) then
    xyzmh_ptmass_in(1:3,1)  = (/x1,y1,z1/)
    xyzmh_ptmass_in(4,1) = mstars(1)
    xyzmh_ptmass_in(5,1) = haccs(1)
    vxyz_ptmass_in(:,1) = (/vx1,vy1,vz1/)

    if (nstar > 1) then
       xyzmh_ptmass_in(1:3,2)  = (/x2,y2,z2/)
       xyzmh_ptmass_in(4,2) = mstars(2)
       xyzmh_ptmass_in(5,2) = haccs(2)
       vxyz_ptmass_in(:,2) = (/vx2,vy2,vz2/)
    endif
 else
    do i=1,nstar
       xyzmh_ptmass_in(1:3,i) = xyzmh_ptmass_in(1:3,i) + xyzstar(:)
       vxyz_ptmass_in(1:3,i)  = vxyz_ptmass_in(1:3,i) + vxyzstar(:)
    enddo
 endif

 call shift_stars(nstar,star,xyzmh_ptmass_in(1:3,1:nstar),vxyz_ptmass_in(1:3,1:nstar),&
                  xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,&
                  npartoftype,nptmass)

 if (id==master) print "(/,a,i10,/)",' Number of gas particles = ',npart

 !
 !--set a few options for the input file
 !
 calc_gravitwaves = .true.
 if (abs(ecc_bh-1.) > epsilon(0.)) then
    theta_gw = theta_bh*180./pi
 else
    theta_gw = -theta_bh*180./pi
 endif

 if (.not.gr .and. iexternalforce == 0) iexternalforce = 1
 if (npart == 0 .and. nptmass == 0) call fatal('setup','no particles setup')
 if (ierr /= 0) call fatal('setup','ERROR during setup')

end subroutine setpart

!----------------------------------------------------------------
!+
!  write .setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setstar,      only:write_options_star,write_options_stars
 use relaxstar,    only:write_options_relax
 use setorbit,     only:write_options_orbit
 use eos,          only:ieos

 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')

 write(iunit,"(a)") '# input file for tidal disruption setup'
 call write_inopt(mhole,  'mhole', 'mass of black hole (solar mass)',  iunit)
 call write_options_stars(star,relax,write_profile,ieos,iunit,nstar)

 write(iunit,"(/,a)") '# options for orbit around black hole'
 call write_inopt(provide_params,'provide_params','manually specify the position and velocity of the star(s)',iunit)
 if (.not. provide_params) then
    call write_inopt(beta,         'beta',         'penetration factor',             iunit)
    call write_inopt(ecc_bh,       'ecc_bh',       'eccentricity (1 for parabolic)', iunit)
    call write_inopt(theta_bh,     'theta_bh',     'inclination of orbit (degrees)', iunit)
    if (ecc_bh >= 1.) then
       call write_inopt(sep_initial,'sep_initial', 'initial separation from BH in tidal radii',iunit)
    endif
    if (abs(ecc_bh-1.) < tiny(0.)) then
       call write_inopt(use_gr_ic,'use_gr_ic', 'whether initial velocity condition computed in GR is used',iunit)
    endif
    call write_inopt(norbits,      'norbits',      'number of orbits',               iunit)
    call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',      iunit)
    if (nstar > 1) then
       call write_options_orbit(orbit,iunit)
    endif
 else
    write(iunit,"(/,a)") '# manual entry of position and velocity for star(s)'
    call write_params(iunit,nstar)
 endif

 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  read .setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 use setstar,      only:read_options_star,read_options_stars
 use relaxstar,    only:read_options_relax
 use physcon,      only:solarm,solarr
 use units,        only:set_units,umass
 use setorbit,     only:read_options_orbit
 use eos,          only:ieos
 character(len=*), intent(in)    :: filename
 integer,          intent(out)   :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 !
 !--read black hole mass and use it to define code units
 !
 call read_inopt(mhole,'mhole',db,min=0.,errcount=nerr)
!  call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
 ! This ensures that we can run simulations with BH's as massive as 1e9 msun.
 ! A BH of mass 1e9 msun would be 1e3 in code units when umass is 1e6*solar masses.
 mass1 = mhole*solarm/umass
 !
 !--read star options and convert to code units
 !
 call read_options_stars(star,ieos,relax,write_profile,db,nerr,nstar)
 !
 !--read options for orbit around black hole
 !
 call read_inopt(provide_params,'provide_params',db,errcount=nerr)
 if (.not. provide_params) then
    call read_inopt(beta,           'beta',           db,min=0.,errcount=nerr)
    call read_inopt(ecc_bh,         'ecc_bh',         db,min=0.,max=1.,errcount=nerr)
    call read_inopt(theta_bh,       'theta_bh',       db,       errcount=nerr)
    if (ecc_bh >= 1.) then
       call read_inopt(sep_initial,'sep_initial',db,errcount=nerr)
    endif
    if (abs(ecc_bh-1.) < tiny(0.)) then
       call read_inopt(use_gr_ic,'use_gr_ic',db,errcount=nerr)
    endif
    call read_inopt(norbits,        'norbits',        db,min=0.,errcount=nerr)
    call read_inopt(dumpsperorbit,  'dumpsperorbit',  db,min=0 ,errcount=nerr)
    if (nstar > 1) then
       call read_options_orbit(orbit,db,nerr)
    endif
 else
    call read_params(db,nerr,nstar)
 endif
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!--------------------------------------------------------------------------
!+
!  write parameters for specifying star positions and velocities manually
!+
!--------------------------------------------------------------------------
subroutine write_params(iunit,nstar)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit,nstar

 call write_inopt(x1, 'x1', 'pos x star 1',iunit)
 call write_inopt(y1, 'y1', 'pos y star 1',iunit)
 call write_inopt(z1, 'z1', 'pos z star 1',iunit)
 call write_inopt(vx1,'vx1','vel x star 1',iunit)
 call write_inopt(vy1,'vy1','vel y star 1',iunit)
 call write_inopt(vz1,'vz1','vel z star 1',iunit)
 if (nstar > 1) then
    call write_inopt(x2, 'x2', 'pos x star 2',iunit)
    call write_inopt(y2, 'y2', 'pos y star 2',iunit)
    call write_inopt(z2, 'z2', 'pos z star 2',iunit)
    call write_inopt(vx2,'vx2','vel x star 2',iunit)
    call write_inopt(vy2,'vy2','vel y star 2',iunit)
    call write_inopt(vz2,'vz2','vel z star 2',iunit)
 endif

end subroutine write_params

!--------------------------------------------------------------------------
!+
!  read parameters for specifying star positions and velocities manually
!+
!--------------------------------------------------------------------------
subroutine read_params(db,nerr,nstar)
 use infile_utils, only:inopts,read_inopt
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 integer,      intent(in)    :: nstar

 call read_inopt(x1, 'x1', db,errcount=nerr)
 call read_inopt(y1, 'y1', db,errcount=nerr)
 call read_inopt(z1, 'z1', db,errcount=nerr)
 call read_inopt(vx1,'vx1',db,errcount=nerr)
 call read_inopt(vy1,'vy1',db,errcount=nerr)
 call read_inopt(vz1,'vz1',db,errcount=nerr)
 if (nstar > 1) then
    call read_inopt(x2, 'x2', db,errcount=nerr)
    call read_inopt(y2, 'y2', db,errcount=nerr)
    call read_inopt(z2, 'z2', db,errcount=nerr)
    call read_inopt(vx2,'vx2',db,errcount=nerr)
    call read_inopt(vy2,'vy2',db,errcount=nerr)
    call read_inopt(vz2,'vz2',db,errcount=nerr)
 endif

end subroutine read_params

!--------------------------------------------------------------------------
!+
!  Obtain velocity in Boyer–Lindquist coordinates using gradient descent method. Iterations are
!  performed until the target specific energy epsilon is reached
!  References: Tejeda et al. (2017), MNRAS 469, 4483–4503
!+
!--------------------------------------------------------------------------
subroutine refine_velocity(x, y, z, vx, vy, vz, M_h, a, r, epsilon_target, alpha, delta_v, tol, max_iters)
 real, intent(in) :: x, y, z, M_h, a, r, epsilon_target, alpha, delta_v, tol
 integer, intent(in) :: max_iters
 real, intent(inout) :: vx, vy, vz
 real :: epsilon_0, epsilon_x, epsilon_y
 real :: d_eps_dx, d_eps_dy
 real :: temp_vx, temp_vy
 real :: sign_epsilon
 integer :: iter

 iter = 0
 do while (iter < max_iters)
  ! Compute epsilon at current velocity
  call compute_epsilon(x, y, z, vx, vy, vz, M_h, a, r, epsilon_0)

  ! Check convergence
  if (abs(epsilon_0 - epsilon_target) < tol) exit

  ! Determine sign of epsilon - epsilon_target
  sign_epsilon = sign(1.0d0, epsilon_0 - epsilon_target)

  ! We only iterate in the x-y plane.
  ! Compute epsilon_x by changing vx by a small delta
  temp_vx = vx + delta_v
  call compute_epsilon(x, y, z, temp_vx, vy, vz, M_h, a, r, epsilon_x)
  d_eps_dx = (epsilon_x - epsilon_0) / delta_v

  ! Compute epsilon_y by changing vy by a small delta
  temp_vy = vy + delta_v
  call compute_epsilon(x, y, z, vx, temp_vy, vz, M_h, a, r, epsilon_y)
  d_eps_dy = (epsilon_y - epsilon_0) / delta_v

  ! Update velocities using gradient descent with sign adjustment
  vx = vx - alpha * sign_epsilon * d_eps_dx
  vy = vy - alpha * sign_epsilon * d_eps_dy

  iter = iter + 1
 end do
 call compute_epsilon(x, y, z, vx, vy, vz, M_h, a, r, epsilon_0)
 if (iter == max_iters) then
  print *, 'Warning: Gradient descent did not converge after ', max_iters, ' iterations.'
 end if
end subroutine refine_velocity

!--------------------------------------------------------------------------
!+
!  Compute GR total specific energy
!  References: Tejeda et al. (2017), MNRAS 469, 4483–4503
!+
!--------------------------------------------------------------------------
subroutine compute_epsilon(x, y, z, vx, vy, vz, M_h, a, r, epsilon_0)
 real, intent(in) :: x, y, z, vx, vy, vz, M_h, a, r
 real, intent(out) :: epsilon_0
 real :: xy_term, r_dot_term, Gamma_flat, Gamma_curved, Gamma_0
 real :: rho_squared, Delta

 rho_squared = r**2 + (a**2 * z**2) / r**2
 Delta = r**2 - 2 * M_h * r + a**2
 xy_term = x * vy - y * vx
 r_dot_term = r**2 * (x * vx + y * vy) + (r**2 + a**2) * z * vz
 Gamma_flat = 1.0d0 - vx**2 - vy**2 - vz**2
 Gamma_curved = (2.0d0 * M_h * r / rho_squared) * ((1.0d0 - a * xy_term / (r**2 + a**2))**2 + &
                        r_dot_term**2 / (r**2 * Delta * (r**2 + a**2)))
 Gamma_0 = 1.0d0 / sqrt(Gamma_flat - Gamma_curved)

 epsilon_0 = Gamma_0 * (1.0d0 - (2.0d0 * M_h * r / rho_squared) * (1.0d0 - a * xy_term / (r**2 + a**2)))
end subroutine compute_epsilon

end module setup
