!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! setup
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters:
!   - beta          : *beta*
!   - ecc           : *ecc*
!   - ellipse_start : *ellipse start at rp?*
!   - impact_param  : *impact parameter*
!   - mhole         : *mass of black hole (solar mass)*
!   - semi_maj_val  : *semi major axis value of binary*
!   - start_at_rp   : *start at rp or before?*
!   - start_sep     : *how far from rp?*
!   - use_gr_ic     : *whether initial velocity condition computed in GR is used*
!
! :Dependencies: eos, infile_utils, io, kernel, metric, mpidomain, orbits,
!   part, physcon, relaxstar, setorbit, setstar, setup_params, systemutils,
!   units
!

 use setstar, only:star_t
 use metric,  only:mass1,a
 implicit none
 public :: setpart

 real               :: beta,impact_param,mhole,ecc,start_sep,semi_maj_val
 integer            :: nstar
 integer, parameter :: max_stars = 2
 type(star_t)       :: star(max_stars)
 logical            :: relax,write_profile,use_gr_ic,start_at_rp,ellipse_start

 private

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use kernel,         only:hfact_default
 use eos,            only:ieos,X_in,Z_in
 use units,          only:set_units,in_code_units
 use setstar,        only:set_defaults_stars,set_stars,shift_stars
 use orbits,         only:refine_velocity
 use part,           only:gravity,eos_vars,nptmass,vxyz_ptmass,xyzmh_ptmass,rad,&
                          nsinkproperties,ihacc
 use io,             only:master,fatal,warning
 use physcon,        only:solarm,pi,solarr
 use systemutils,    only:get_command_option
 use mpidomain,      only:i_belong
 use setup_params,   only:rhozero,npart_total
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
 logical :: iexist,use_var_comp
 real    :: mstars(max_stars),rstars(max_stars),haccs(max_stars)
 real    :: r_binary,rtidal,rp,vhat(3),x0,r0,y0,vel,m_binary
 real    :: alpha,delta_v,epsilon_target,tol,a_val
 real    :: vxyzstar(3),xyzstar(3)
 real    :: rx,ry,rz,rel_r(3),rel_v_mag,rel_v(3)
 real    :: x1(3),x2(3),v1(3),v2(3),x0bin(3,2),v0bin(3,2)
 integer :: max_iters
 integer :: ierr,np_default,i
!
!-- general parameters
!
 hfact = hfact_default
 time  = 0.d0
 polyk = 1.e-10    ! <== uconst
 gamma = 5./3.
 ieos  = 2
 if (.not.gravity) call fatal('setup','recompile with GRAVITY=yes')
!
!-- space available for injected gas particles
!
 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.d0
 vxyzu(:,:)     = 0.d0
 nstar          = 2
!
!-- Default runtime parameters
!
 mhole           = 1.e6  ! (solar masses)
 call set_units(mass=mhole*solarm,c=1.d0,G=1.d0)
 call set_defaults_stars(star)
 star(:)%m       = '1.*msun'
 star(:)%r       = '1.*solarr'
 np_default      = 1e6
 ecc             = 1.0
 start_sep       = 5.0
 star%np         = int(get_command_option('np',default=np_default)) ! can set default value with --np=1e5 flag (mainly for testsuite)
 beta            = 1.0
 impact_param    = 0.5 ! units of rsun
 semi_maj_val    = 10.130607675019033
 relax           = .true.
 use_var_comp    = .false.
 write_profile   = .true.
 ellipse_start   = .false.
 start_at_rp     = .true.
 use_gr_ic       = .true. ! Whether initial velocity condition computed in GR is used for parabolic orbits.
!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Collision of stars around BH'
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
!
!--save the radius and mass of the stars
!
 do i=1,nstar
    rstars(i) = in_code_units(star(i)%r,ierr,unit_type='length')
    if (ierr /= 0) call fatal('setup','could not convert rstar to code units',i=i)
    mstars(i) = in_code_units(star(i)%m,ierr,unit_type='mass')
    if (ierr /= 0) call fatal('setup','could not convert mstar to code units',i=i)
    haccs(i) = in_code_units(star(i)%hacc,ierr,unit_type='mass')
    if (ierr /= 0) call fatal('setup','could not convert hacc to code units',i=i)
 enddo
 print*, rstars(1:2),'rstars',star(1)%r,star(2)%r
!
!--determine tidal radius of the binary
!
 r_binary = rstars(1) + rstars(2) ! similar to a TDE of a single star we set rstar = r1 + r2, this is also the relative distance

 rtidal  = semi_maj_val * (mass1 / (mstars(1) + mstars(2)))**(1./3.)
 rp = rtidal / beta
 r0 = rp
 !
 !--determine position and velocity of a binary system at the point of collision
 !
 ! for parabolic orbit, set the binary at the point of collision at rp or before tidal radius

 if (ecc < 1.d0) then ! elliptical orbits
    a_val = rp / (1.d0 - ecc)
    if (ellipse_start) then
       r0 = rp  ! place on pericentre
    else
       r0 = a_val * (1.d0 + ecc)  ! place on apocentre
    endif
    vel = sqrt(mass1 * ((2.d0 / r0 - 1.d0 / a_val)))
    x0  = 0.d0
    y0  = r0
    xyzstar(:)  = (/x0,-y0,0./)
    vhat        = (/1.d0,0.d0,0.d0/)
    vxyzstar(:) = vel*vhat
 else
    if (start_at_rp) then
       r0  = rp
       vel = sqrt(2.d0 * mass1 / r0)
       x0  = 0.d0
       y0  = r0
       xyzstar(:)  = (/x0,-y0,0./)
       vhat        = (/1.d0,0.d0,0.d0/)
       vxyzstar(:) = vel*vhat
    else
       r0  = start_sep * rtidal ! set the binary at the point of collision at 5*rp
       vel = sqrt(2.d0 * mass1 / r0)
       y0  = -2.d0 * rp + r0
       x0  = sqrt(r0**2 - y0**2)
       xyzstar(:)  = (/-x0,y0,0./)
       vel         = sqrt(2.*mass1/r0)
       vhat        = (/2.*rp,-x0,0./)/sqrt(4.*rp**2 + x0**2)
       vxyzstar(:) = vel*vhat
    endif
 endif

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
    if (ecc == 1.d0) then
       call refine_velocity(-x0, y0, 0.0d0, vxyzstar(1), vxyzstar(2), vxyzstar(3), &
                            mass1, 0.0d0, r0, epsilon_target, alpha, delta_v, tol, max_iters)
    endif
 endif
 vel = norm2(vxyzstar)
!
!--calculate position and velocity of the stars around the BH
!
 ry = sqrt(r_binary**2 - impact_param**2) ! using pythagoras theorem to determine x position
 rx = impact_param ! y position is same as impact parameter at point of collision
 rz = 0.d0
 rel_r = (/rx,ry,rz/)

 m_binary = mstars(1) + mstars(2)
 rel_v_mag = sqrt(2 * m_binary / r_binary) ! using sqrt(2*G*(m1+m2)/(r1+r2))
 rel_v   = (/0.d0,rel_v_mag,0.d0/)

 x1 =  (mstars(2) / m_binary) * rel_r + xyzstar
 x2 = -(mstars(1) / m_binary) * rel_r + xyzstar

 v1 = -(mstars(2) / m_binary) * rel_v + vxyzstar
 v2 =  (mstars(1) / m_binary) * rel_v + vxyzstar

 x0bin(:,1) = x1
 x0bin(:,2) = x2
 v0bin(:,1) = v1
 v0bin(:,2) = v2

 call shift_stars(nstar,star,x0bin,v0bin,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                 npart,npartoftype,nptmass)

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

 write(iunit,"(a)") '# input file for binary collision around BH setup'
 call write_inopt(mhole,          'mhole', 'mass of black hole (solar mass)',iunit)
 call write_inopt(beta,           'beta',           'beta',                  iunit)
 call write_inopt(ecc,            'ecc',            'ecc',                   iunit)
 call write_inopt(semi_maj_val,   'semi_maj_val',   'semi major axis value of binary', iunit)
 call write_inopt(impact_param,   'impact_param',   'impact parameter',      iunit)
 call write_inopt(use_gr_ic,      'use_gr_ic',      'whether initial velocity condition computed in GR is used',iunit)
 if (ecc < 1.d0) then
    call write_inopt(ellipse_start,   'ellipse_start',   'ellipse start at rp?',      iunit)
    call write_options_stars(star,relax,write_profile,ieos,iunit,nstar)
 else
    call write_inopt(start_at_rp,     'start_at_rp',     'start at rp or before?',    iunit)
    if (.not. start_at_rp) then
       call write_inopt(start_sep,   'start_sep',       'how far from rp?',          iunit)
    endif
    call write_options_stars(star,relax,write_profile,ieos,iunit,nstar)
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
 mass1 = mhole*solarm/umass
 !
 !--read star options and convert to code units
 !
 call read_inopt(beta,           'beta',           db,min=0.,errcount=nerr)
 call read_inopt(ecc,            'ecc',            db,errcount=nerr)
 call read_inopt(semi_maj_val,   'semi_maj_val',   db,errcount=nerr)
 call read_inopt(impact_param,   'impact_param',   db,min=0.,max=1.,errcount=nerr)
 call read_inopt(use_gr_ic,      'use_gr_ic',      db,errcount=nerr)
 if (ecc < 1.d0) then
    call read_inopt(ellipse_start,   'ellipse_start', db,errcount=nerr)
    call read_options_stars(star,ieos,relax,write_profile,db,nerr,nstar)
 else
    call read_inopt(start_at_rp,    'start_at_rp',    db,errcount=nerr)
    if (.not. start_at_rp) then
       call read_inopt(start_sep,  'start_sep',      db,min=0.,errcount=nerr)
    endif
    call read_options_stars(star,ieos,relax,write_profile,db,nerr,nstar)

 endif
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile


end module setup
