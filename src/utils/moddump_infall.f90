!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Adds either a sphere or cylinder of material to infall
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - add_turbulence     : *add turbulence (0=no, 1=yes)*
!   - b                  : *impact parameter [code units]*
!   - b_frac             : *impact parameter b as fraction of b_crit*
!   - cloud_control_mode : *cloud control mode (0=manual, 1=N sets size, 2=size sets mass)*
!   - eccentricity       : *eccentricity*
!   - ieos_infall        : *eos to set after infall (6=isothermal about sink, 14=binary)*
!   - in_mass            : *infall mass per cloud (sym_infall > 0 doubles total)*
!   - in_orbit           : *orbit type (0=parabolic, 1=hyperbolic)*
!   - in_shape           : *infall material shape (0=sphere, 1=ellipse)*
!   - incx               : *rotation on x axis (deg)*
!   - incy               : *rotation on y axis (deg)*
!   - incz               : *rotation on z axis (deg)*
!   - isink              : *index of the sink the eos is centred on*
!   - m_gas              : *gas particle mass in Msun for empty simulations (0 = use existing)*
!   - n_add              : *number of particles per cloud (sym_infall > 0 doubles total)*
!   - r_a                : *semi-major axis of ellipse*
!   - r_close            : *closest approach [code units]*
!   - r_in               : *radius of shape (or semi-minor axis)*
!   - r_init             : *initial radial distance*
!   - r_slope            : *density power law index*
!   - r_soft             : *softening radius*
!   - rho_mode           : *density mode (0=current, 1=Dullemond Eq4/Eq5)*
!   - rms_mach           : *rms Mach number*
!   - tfact              : *tfact*
!   - v_inf              : *velocity at infinity [km/s]*
!
! :Dependencies: datafiles, dim, eos, infile_utils, io, kernel,
!   moddump_utils, options, part, partinject, physcon, prompting,
!   spherical, stretchmap, units, vectorutils, velfield
!
 use moddump_utils, only:prompt_for_params,write_moddump_header
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer, parameter :: nr = 200
 ! runtime parameters, written to/read from the prefix.moddump file
 integer :: in_shape = 1, in_orbit = 1, add_turbulence = 0, sym_infall = 0
 integer :: rho_mode = 0, cloud_control_mode = 0, n_add = 0
 real    :: in_mass = 0.001, r_in = 250.0, r_a = 3500.0, r_init = 4000.0, r_close = 100.0
 real    :: v_inf = 1.0, b = 0.0, b_frac = 1.0, ecc = 0.0
 real    :: incx = 0.0, incy = 0.0, incz = 0.0, rms_mach = 1.0, tfact = 0.0
 real    :: r_slope = 0.0
 real    :: r_soft = 100.0
 real    :: m_gas = 0.0
 integer :: ieos_infall = 6

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:use_dust
 use partinject,     only:add_or_update_particle
 use options,        only:use_dustfrac
 use part,           only:igas,xyzmh_ptmass,nptmass,gravity, &
                          dustfrac
 use units,          only:udist,utime,get_G_code
 use io,             only:id,master,fatal
 use spherical,      only:set_sphere,set_ellipse
 use stretchmap,     only:rho_func
 use kernel,         only:hfact_default
 use physcon,        only:pi,au
 use vectorutils,    only:cross_product3D,rotatevec
 use eos,            only:ieos,isink,get_spsound
 use velfield,       only:set_velfield_from_cubes
 use datafiles,      only:find_phantom_datafile

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, allocatable :: xyzh_add(:,:),vxyzu_add(:,:)
 integer :: ipart,i,np,ierr
 integer(kind=8) :: nptot
 real    :: pmass,delta,big_omega,b_crit
 real    :: vp(3), xp(3), rellipsoid(3)
 real    :: dma,n0,pf,m0,x0,y0,z0,r0,vx0,vy0,vz0,mtot,tiny_number,n1
 real    :: y1,x1,dx,x_prime,y_prime
 real    :: unit_velocity,G,rms_in,vol_obj,rhoi,spsound,factor,my_vrms,vxi,vyi,vzi
 real    :: v_inf_code,mu,hyp_h,vr0,vt0,v0
 real    :: rho_cloud_cgs, rho_cloud, mu_cloud, r_equiv
 real    :: rhat(3),that(3)
 real    :: pos_sym(3),vel_sym(3),L_part(3),L_cloud(3),L_sym(3),cloud_centre(3)
 real    :: sym_axis(3),sym_angle
 real    :: dustfrac_tmp
 logical :: lrhofunc,empty_sim
 character(len=20), parameter :: filevx = 'cube_v1.dat'
 character(len=20), parameter :: filevy = 'cube_v2.dat'
 character(len=20), parameter :: filevz = 'cube_v3.dat'
 character(len=120)           :: filex,filey,filez
 procedure(rho_func), pointer :: prhofunc

 ! Local working variables; retain parameters read by the driver.
 big_omega = 0.
 tiny_number = 1e-4
 lrhofunc = .false.
 empty_sim = .false.
 ierr = 0
 my_vrms = 0.
 pf = 0.
 rho_cloud_cgs = 0.
 rho_cloud = 0.
 mu_cloud = 2.3
 r_equiv = 0.
 v_inf_code = 0.
 mu = 0.
 pmass = massoftype(igas)
 x0 = 0.
 y0 = 0.
 z0 = 0.
 vol_obj = 0.

 ! The driver has read the .moddump file; prompt only on the first run.
 if (prompt_for_params) call read_interactive_moddumpfile()

 ! udist default is cm
 unit_velocity = udist/utime ! cm/s
 G = get_G_code()

 if (npartoftype(igas) <= 0) then
    empty_sim = .true.
    pmass = 0.0
 endif

 if (gravity) then
    write(*,*) "Disc self-gravity is on. Including disc mass in cloud orbit calculation."
    mtot=sum(xyzmh_ptmass(4,:)) + npartoftype(igas)*massoftype(igas)
 else
    mtot=sum(xyzmh_ptmass(4,:))
 endif

 ! gas particle mass for empty simulations
 if ((cloud_control_mode == 1 .or. cloud_control_mode == 2) .and. empty_sim) then
    write(*,*) "WARNING: Simulation has no mass, massoftype(igas) is not set."
    if (m_gas <= 0.) call fatal('moddump_infall',&
       'empty simulation: set m_gas (gas particle mass in Msun) in the .moddump file')
    pmass = m_gas
    massoftype(igas) = pmass
 elseif (cloud_control_mode == 0 .and. empty_sim) then
    if (n_add <= 0) call fatal('moddump_infall','n_add must be > 0 for empty simulations')
    pmass = in_mass/real(n_add)
 endif

 if (cloud_control_mode == 0) then
    n_add = int(in_mass/pmass)
 endif

 if (cloud_control_mode == 1 .and. .not. empty_sim) then
    if (n_add <= 0 .and. pmass > 0.) n_add = int(in_mass/pmass)
 endif

 if (cloud_control_mode == 1) then
    ! n_add sets size (basically in_mass), this sets r_equiv
    in_mass = real(n_add)*pmass
    r_equiv = (in_mass/0.01)**(1./2.3) * 5000.
 elseif (in_shape == 0) then
    r_equiv = r_in
 else
    r_equiv = (r_in*r_in*r_a)**(1./3.)
 endif

 if (cloud_control_mode == 2) then
    ! size sets mass, we already set r_in, r_a, so get the mass
    in_mass = 0.01*(r_equiv/5000.)**2.3
    if (pmass > 0.) n_add = int(in_mass/pmass)
 endif

 if (n_add <= 0) call fatal('moddump_infall','number of particles to add must be > 0')

 vol_obj = (4.0/3.0)*pi*r_equiv**3

 if (cloud_control_mode /= 0) then
    r_slope = 0.
    lrhofunc = .false.
 endif
 if (cloud_control_mode == 0 .and. r_slope > tiny_number) then
    prhofunc => rhofunc
    lrhofunc = .true.
 endif

 if (in_orbit == 1) then
    write(*,*) "Hyperbolic orbit, see Dullemond+2019 for parameter definitions."
    if (v_inf <= tiny_number) call fatal('moddump_infall','v_inf must be > 0 for a hyperbolic orbit')
    v_inf_code = v_inf*1.e5/unit_velocity
    mu = G*mtot
    b_crit = mu/v_inf_code**2
    write(*,*) "Critical impact parameter, b_crit, is ", b_crit, " code units, ", b_crit*udist/au, " au"
    b = b_frac * b_crit
    ecc = sqrt(1 + b**2/b_crit**2)
    r_close = b * sqrt((ecc-1)/(ecc+1))
    write(*,*) "Eccentricity of the cloud is ", ecc
    write(*,*) "Closest approach of cloud center will be ", r_close, " code units, ", r_close*udist/au, " au."
 endif

 select case (in_orbit)
 case (0)
    ! Parabolic orbit, taken from set_flyby
    dma = r_close
    n0  = r_init/r_close
    !--focal parameter dma = pf/2
    pf = 2*dma

    !--define m0 = -x0/dma such that r0 = n0*dma
    !  companion starts at negative x and y
    !  positive root of 1/8*m**4 + m**2 + 2(1-n0**2) = 0
    !  for n0 > 1
    m0 = 2*sqrt(n0-1.0)

    !--perturber initial position
    x0 = -m0*dma
    y0 = dma*(1.0-(x0/pf)**2)
    z0 = 0.0
    xp = (/x0,y0,z0/)
 case (1)
    ! Dullemond+2019
    ! Initial position is x=r_init and y=b (impact parameter)
    if (in_shape==1) then
       x0 = -(r_init+r_in)
    else
       x0 = -r_init
    endif
    y0 = -b
    z0 = 0.0
    xp = (/x0, y0, z0/)
 end select
 write(*,*) "Initial centre is: ", xp

 allocate(xyzh_add(4,n_add+int(0.1*n_add)),vxyzu_add(4,n_add+int(0.1*n_add)))
 delta = 1.0 ! no idea what this is
 nptot = n_add + npartoftype(igas)
 np = 0

 if (in_shape == 0) then
    if (lrhofunc) then
       call set_sphere('random',id,master,0.,r_in,delta,hfact_default,np,xyzh_add,xyz_origin=xp,&
                     np_requested=n_add, nptot=nptot, rhofunc=prhofunc)
    else
       call set_sphere('random',id,master,0.,r_in,delta,hfact_default,np,xyzh_add,xyz_origin=xp,&
                     np_requested=n_add, nptot=nptot)
    endif
    write(*,*) "The sphere has been succesfully initialised."
 elseif (in_shape == 1) then
    rellipsoid(1) = r_in
    rellipsoid(2) = r_a
    rellipsoid(3) = r_in
    call set_ellipse('random',id,master,rellipsoid,delta,hfact_default,xyzh_add,np,&
                    np_requested=n_add, nptot=nptot)

    ! Need to correct the ellipse
    do i = 1,n_add
       xyzh_add(1, i) = xyzh_add(1, i) + xp(1)
       xyzh_add(2, i) = xyzh_add(2, i) + xp(2)
       xyzh_add(3, i) = xyzh_add(3, i) + xp(3)
    enddo
    if (in_orbit == 0) then
       do i = 1,n_add
          x1 = xyzh_add(1, i)
          y1 = xyzh_add(2, i)
          dma = r_close
          n0  = (sqrt(xp(1)**2 + xp(2)**2))/dma
          pf = 2*dma
          n1 = (xp(2)-y1)/dma
          m0 = 2*sqrt(n0-n1-1.0)
          x0 = -m0*dma
          y0 = dma*(1.0-(x0/pf)**2)
          dx = xyzh_add(1, i) - xp(1)
          y_prime = 4*dma/x0 *dx
          x_prime = dx

          xyzh_add(1, i) = x0 + x_prime
          xyzh_add(2, i) = y0 + y_prime
       enddo
    endif
    write(*,*) "The ellipse has been succesfully initialised."
 endif

 !--Set velocities (from pre-made velocity cubes)
 vxyzu_add(:,:) = 0.

 if (add_turbulence==1) then

    write(*,"(1x,a)") 'Setting up velocity field on the particles...'

    filex = find_phantom_datafile(filevx,'velfield')
    filey = find_phantom_datafile(filevy,'velfield')
    filez = find_phantom_datafile(filevz,'velfield')

    call set_velfield_from_cubes(xyzh_add,vxyzu_add,n_add,filex,filey,filez,1.,tfact*r_in,.false.,ierr)

    if (ierr /= 0) call fatal('setup','error setting up velocity field')

    if (in_shape == 0) then
       vol_obj = (4.0/3.0)*pi*r_in**3
    elseif (in_shape == 1) then
       vol_obj = (4.0/3.0)*pi*r_in*r_in*r_a
    endif

    rhoi = in_mass/vol_obj
    spsound = get_spsound(ieos,xp,rhoi,vxyzu_add(:,1)) ! eos_type,xyzi,rhoi,vxyzui
    rms_in = spsound*rms_mach

    do i=1,n_add
       vxi  = vxyzu_add(1,i)
       vyi  = vxyzu_add(2,i)
       vzi  = vxyzu_add(3,i)
       my_vrms = my_vrms + vxi*vxi + vyi*vyi + vzi*vzi
    enddo

    ! Normalise velocity field
    my_vrms = sqrt(1/real(n_add) * my_vrms)
    factor = rms_in/my_vrms
    do i=1,n_add
       vxyzu_add(1:3,i) = vxyzu_add(1:3,i)*factor
    enddo
 endif

 ! Set up velocities
 if (in_orbit == 0) then

    !--perturber initial velocity
    r0  = sqrt(x0**2+y0**2+z0**2)
    vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
    vy0 = -(x0/r0)*sqrt(mtot/pf)
    vz0 = 0.0
    vp  = (/vx0,vy0,vz0/)
    if (in_shape == 0) then
       ! Initiate initial velocity of the particles in the shape
       vxyzu_add(1, :) = vxyzu_add(1, :) + vx0
       vxyzu_add(2, :) = vxyzu_add(2, :) + vy0
       vxyzu_add(3, :) = vxyzu_add(3, :) + vz0

    elseif (in_shape == 1) then
       do i=1,n_add
          x0 = xyzh_add(1, i)
          y0 = xyzh_add(2, i)
          z0 = xyzh_add(3, i)

          r0  = sqrt(x0**2+y0**2+z0**2)
          vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
          vy0 = -(x0/r0)*sqrt(mtot/pf)
          vz0 = 0.0

          vxyzu_add(1, i) = vxyzu_add(1, i) + vx0
          vxyzu_add(2, i) = vxyzu_add(2, i) + vy0
          vxyzu_add(3, i) = vxyzu_add(3, i) + vz0

       enddo
       vxyzu_add(4, :) = vxyzu(4, 1)
    endif
 elseif (in_orbit == 1) then
    ! Dullemond+2019
    ! Initial centre velocity from the finite-radius hyperbolic orbit.
    r0 = sqrt(dot_product(xp,xp))
    if (r0 <= r_close) call fatal('moddump_infall','initial hyperbolic position must be outside r_close')
    hyp_h = b*v_inf_code
    v0 = sqrt(v_inf_code**2 + 2.0*mu/r0)
    vt0 = hyp_h/r0
    if (vt0 > v0) call fatal('moddump_infall','invalid hyperbolic orbit: tangential speed exceeds total speed')
    vr0 = -sqrt(max(0.0,v0**2 - vt0**2))
    rhat = xp/r0
    that = (/ -rhat(2), rhat(1), 0.0 /)
    vx0 = vr0*rhat(1) + vt0*that(1)
    vy0 = vr0*rhat(2) + vt0*that(2)
    vz0 = 0.0
    vp = (/vx0, vy0, vz0/)
    vxyzu_add(1, :) = vxyzu_add(1, :) + vx0
    vxyzu_add(2, :) = vxyzu_add(2, :) + vy0
    vxyzu_add(3, :) = vxyzu_add(3, :) + vz0
 endif

 write(*,*) "Initial velocity of object centre is ", vp

 if (use_dust) then
    if (use_dustfrac) then
       write(*,*) "Detected one-fluid dust in the simulation, adding smallest dust to infall."
       ! Set the dustfrac to the global dust to gas ratio

       dustfrac_tmp = sum(dustfrac)/npartoftype(igas)

       write(*,*) "The total dustfrac is ", dustfrac_tmp

       ! Set the dustfrac to the dustfrac of the smallest bin
       dustfrac_tmp = sum(dustfrac(1,:))/npartoftype(igas)
       write(*,*) "The single bin dustfrac is ", dustfrac_tmp

    endif
 endif

 ! Rotate the new cloud into its final orientation before constructing an
 ! angular-momentum-balancing companion.  Applying the same proper rotation
 ! to position and velocity preserves the cloud structure and orbital energy.
 do i = 1,n_add
    ! Rotate particle to correct position and velocity
    ! First rotate to get the right initial position
    ! Need to do this due to the parabolic orbit notation
    ! xyzh_add(4,i) = 1.0
    if (in_orbit == 0) then
       call rotatevec(xyzh_add(1:3,i),(/0.,-1.,0./),pi)
       call rotatevec(vxyzu_add(1:3,i),(/0.,-1.,0./),pi)
    endif

    ! Now rotate around x axis
    call rotatevec(xyzh_add(1:3,i),(/1.,0.,0./),incx*pi/180.)
    call rotatevec(vxyzu_add(1:3,i),(/1.,0.,0./),incx*pi/180.)

    call rotatevec(xyzh_add(1:3,i),(/0.,1.,0./),incy*pi/180.)
    call rotatevec(vxyzu_add(1:3,i),(/0.,1.,0./),incy*pi/180.)

    call rotatevec(xyzh_add(1:3,i),(/0.,0.,1./),incz*pi/180.)
    call rotatevec(vxyzu_add(1:3,i),(/0.,0.,1./),incz*pi/180.)
 enddo

 if (sym_infall > 0) then
    L_cloud = 0.
    cloud_centre = 0.
    do i = 1,n_add
       call cross_product3D(xyzh_add(1:3,i),vxyzu_add(1:3,i),L_part)
       L_cloud = L_cloud + L_part
       cloud_centre = cloud_centre + xyzh_add(1:3,i)
    enddo
    cloud_centre = cloud_centre/real(n_add)
    call get_symmetry_rotation(sym_infall,L_cloud,cloud_centre,sym_axis,sym_angle,L_sym)
    write(*,*) 'Original cloud angular momentum (common particle-mass factor omitted): ',L_cloud
    write(*,*) 'Companion cloud angular momentum: ',L_sym
    write(*,*) 'The companion doubles n_add and in_mass for the total infall.'
 endif

 ! Add the original cloud and, when requested, its rotated companion.
 ! n_add and in_mass describe each individual cloud, matching the convention
 ! used by the Pineda streamer injector.
 ipart = npart ! The initial particle number (post shuffle)
 do i = 1,n_add
    ! Add the particle
    ipart = ipart + 1
    call  add_or_update_particle(igas, xyzh_add(1:3,i), vxyzu_add(1:3,i), xyzh_add(4,i), &
                                vxyzu_add(4,i), ipart, npart, npartoftype, xyzh, vxyzu)
    if (use_dust) then
       if (use_dustfrac) then
          dustfrac(1, ipart) = dustfrac_tmp
       endif
    endif

    if (sym_infall > 0) then
       pos_sym = xyzh_add(1:3,i)
       vel_sym = vxyzu_add(1:3,i)
       call rotatevec(pos_sym,sym_axis,sym_angle)
       call rotatevec(vel_sym,sym_axis,sym_angle)
       ipart = ipart + 1
       call add_or_update_particle(igas,pos_sym,vel_sym,xyzh_add(4,i),vxyzu_add(4,i), &
                                  ipart,npart,npartoftype,xyzh,vxyzu)
       if (use_dust .and. use_dustfrac) dustfrac(1,ipart) = dustfrac_tmp
    endif

 enddo

 ! The input disc typically uses ieos=3 (locally isothermal about the
 ! origin), which no longer makes sense once infall has been added. Swap
 ! it for the user-requested eos: 6 (locally isothermal about sink isink)
 ! or 14 (binary).
 if (ieos==3) then
    select case (ieos_infall)
    case (6)
       if (isink > nptmass) call fatal('moddump_infall',&
          'ieos_infall=6 requires isink <= nptmass')
       write(*,*) "WARNING: Changing ieos from 3 to 6, centred on sink ", isink
       ieos = 6
    case (14)
       write(*,*) "WARNING: Changing ieos from 3 to 14."
       ieos = 14
    end select
 endif
 write(*,*)  " ###### Added infall successfully ###### "
 deallocate(xyzh_add,vxyzu_add)

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  Return the pi rotation that gives a companion cloud with the
!  requested angular-momentum symmetry.  A proper rotation R applied
!  to both r and v gives L' = R L exactly, including any turbulent
!  contribution to the cloud angular momentum.
!+
!----------------------------------------------------------------
subroutine get_symmetry_rotation(mode,L,cloud_centre,axis,angle,L_target)
 use physcon,     only:pi
 use vectorutils, only:cross_product3D,make_perp_frame
 integer, intent(in) :: mode
 real,    intent(in) :: L(3),cloud_centre(3)
 real,    intent(out) :: axis(3),angle,L_target(3)
 real :: axis_norm,L_norm,L_xy,scale,eperp(3)

 L_norm = sqrt(dot_product(L,L))
 angle = pi

 select case(mode)
 case(1)
    ! A pi rotation about z maps (Lx,Ly,Lz) to (-Lx,-Ly,Lz).
    axis = (/0.,0.,1./)
    L_target = (/-L(1),-L(2),L(3)/)
 case(2)
    ! Choose an axis perpendicular to L, so a pi rotation maps L to -L.
    call cross_product3D(cloud_centre,L,axis)
    L_target = -L
 case(3)
    ! Rotating by pi about the projection of L into the xy plane maps
    ! (Lx,Ly,Lz) to (Lx,Ly,-Lz).
    axis = (/L(1),L(2),0./)
    L_target = (/L(1),L(2),-L(3)/)
 end select

 axis_norm = sqrt(dot_product(axis,axis))
 select case(mode)
 case(2)
    scale = sqrt(dot_product(cloud_centre,cloud_centre))*L_norm
 case(3)
    scale = L_norm
 case default
    scale = 1.
 end select
 if (axis_norm <= 100.*epsilon(1.)*max(scale,tiny(1.))) then
    if (L_norm > tiny(1.)) then
       call make_perp_frame(L/L_norm,axis,eperp)
       axis_norm = sqrt(dot_product(axis,axis))
    else
       ! A zero-angular-momentum cloud is already balanced; use a harmless
       ! fixed rotation to keep the companion spatially distinct when possible.
       axis = (/0.,0.,1./)
       axis_norm = 1.
       L_target = 0.
    endif
 endif
 axis = axis/axis_norm

 ! In mode 3 an exactly z-directed L has no xy projection.  The fallback
 ! perpendicular axis above correctly maps Lz to -Lz.
 L_xy = sqrt(L(1)*L(1) + L(2)*L(2))
 if (mode == 3 .and. L_xy <= 100.*epsilon(1.)*max(L_norm,tiny(1.))) L_target = -L

end subroutine get_symmetry_rotation

real function rhofunc(r)
 real, intent(in) :: r

 rhofunc = 1./(abs(r) + r_soft)**(r_slope)

end function rhofunc

!----------------------------------------------------------------
!+
!  set parameters interactively (when no .moddump file is found)
!+
!----------------------------------------------------------------
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt
 use eos,       only:isink
 use part,      only:npartoftype,igas

 call prompt('Enter the infall material shape (0=sphere, 1=ellipse)',in_shape,0,1)
 call prompt('Enter cloud control mode (0=manual mass+size (rho not fixed), 1=mass/n_add set radius (fixed rho), '&
            //'2=radius sets mass/n_add (fixed rho))',cloud_control_mode,0,2)

 if (cloud_control_mode == 0 .or. cloud_control_mode == 2) then
    if (in_shape == 0) then
       call prompt('Enter radius of shape:', r_in, 0.1)
    else
       call prompt('Enter semi-minor axis of ellipse:', r_in, 0.1)
       call prompt('Enter semi-major axis of ellipse:', r_a, 0.1)
    endif
 endif

 if (cloud_control_mode == 0) then
    call prompt('Enter infall mass in Msun:', in_mass, 0.0)
    if (npartoftype(igas) <= 0) call prompt('Enter number of particles to add:',n_add,1)
    call prompt('Enter value of power-law density along radius:', r_slope, 0.0)
    if (r_slope > 0.) call prompt('Enter softening radius:', r_soft, 0.1)
 elseif (cloud_control_mode == 1) then
    call prompt('Enter number of particles to add:', n_add, 0)
 endif

 if (cloud_control_mode == 1 .or. cloud_control_mode == 2) then
    call prompt('Gas particle mass in Msun for empty simulations (0 = use existing):', m_gas, 0.0)
 endif

 call prompt('Enter orbit type (0=parabolic, 1=hyperbolic)', in_orbit, 0, 1)
 if (in_orbit == 0) then
    call prompt('Enter closest approach in au:', r_close, 0.)
 else
    call prompt('Enter cloud velocity at infinity, v_inf, in km/s:', v_inf, 0.0)
    call prompt('Enter impact parameter b as a ratio of b_crit:', b_frac, 0.0)
 endif

 call prompt('Enter initial radial distance in au:', r_init, 0.)

 call prompt('Add turbulence to the gas? (0=no, 1=yes)', add_turbulence, 0, 1)
 if (add_turbulence == 1) then
    call prompt('Enter rms Mach number:', rms_mach, 0., 20.)
    call prompt('Enter tfact:', tfact, 0.0)
 endif

 call prompt('Enter rotation on x axis (deg):', incx, -360., 360.)
 call prompt('Enter rotation on y axis (deg):', incy, -360., 360.)
 call prompt('Enter rotation on z axis (deg):', incz, -360., 360.)

 call prompt('Balance infall angular momentum (0=no, 1=Lx/Ly, 2=all, 3=Lz)',sym_infall,0,3)

 call prompt('Enter eos to set after infall (6=isothermal about a sink, 14=binary):', ieos_infall)
 do while (ieos_infall /= 6 .and. ieos_infall /= 14)
    call prompt('Please enter either 6 or 14:', ieos_infall)
 enddo
 if (ieos_infall == 6) call prompt('Enter index of the sink the eos is centred on:', isink, 1)

end subroutine read_interactive_moddumpfile

!----------------------------------------------------------------
!+
!  write the moddump parameters to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 use eos,          only:isink
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# Infall parameters'
 call write_inopt(in_shape,'in_shape','infall material shape (0=sphere, 1=ellipse)',iunit)
 call write_inopt(in_orbit,'in_orbit','orbit type (0=parabolic, 1=hyperbolic)',iunit)
 call write_inopt(rho_mode,'rho_mode','density mode (0=current, 1=Dullemond Eq4/Eq5)',iunit)
 call write_inopt(cloud_control_mode,'cloud_control_mode','cloud control mode (0=manual, 1=N sets size, 2=size sets mass)',iunit)
 call write_inopt(n_add,'n_add','number of particles per cloud (sym_infall > 0 doubles total)',iunit)
 call write_inopt(in_mass,'in_mass','infall mass per cloud (sym_infall > 0 doubles total)',iunit)
 call write_inopt(m_gas,'m_gas','gas particle mass in Msun for empty simulations (0 = use existing)',iunit)
 call write_inopt(r_in,'r_in','radius of shape (or semi-minor axis)',iunit)
 if (in_shape==1) call write_inopt(r_a,'r_a','semi-major axis of ellipse',iunit)

 call write_inopt(r_slope,'r_slope','density power law index',iunit)
 if (abs(r_slope) > 0.) call write_inopt(r_soft,'r_soft','softening radius',iunit)

 call write_inopt(r_init,'r_init','initial radial distance',iunit)

 if (in_orbit==0) call write_inopt(r_close,'r_close','closest approach',iunit)
 if (in_orbit==1) then
    call write_inopt(v_inf,'v_inf','velocity at infinity [km/s]',iunit)
    call write_inopt(b_frac,'b_frac','impact parameter b as fraction of b_crit',iunit)
    call write_inopt(b,'b','impact parameter [code units]',iunit)
    call write_inopt(ecc,'ecc','eccentricity',iunit)
    call write_inopt(r_close,'r_close','closest approach [code units]',iunit)
 endif

 call write_inopt(incx,'incx','rotation on x axis (deg)',iunit)
 call write_inopt(incy,'incy','rotation on y axis (deg)',iunit)
 call write_inopt(incz,'incz','rotation on z axis (deg)',iunit)
 call write_inopt(sym_infall,'sym_infall', &
                  'companion cloud angular momentum (0=no, 1=Lx/Ly, 2=all, 3=Lz)',iunit)

 call write_inopt(add_turbulence,'add_turbulence','add turbulence (0=no, 1=yes)',iunit)
 if (add_turbulence==1) then
    call write_inopt(rms_mach,'rms_mach','rms Mach number',iunit)
    call write_inopt(tfact,'tfact','tfact',iunit)
 endif

 call write_inopt(ieos_infall,'ieos_infall','eos to set after infall (6=isothermal about sink, 14=binary)',iunit)
 if (ieos_infall==6) call write_inopt(isink,'isink','index of the sink the eos is centred on',iunit)

 close(iunit)

end subroutine write_moddump

!----------------------------------------------------------------
!+
!  read the moddump parameters from the .moddump file; ierr counts
!  missing or invalid options (so get_moddump_options can top up the file
!  and ask the user to edit it)
!+
!----------------------------------------------------------------
subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use eos,          only:isink
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 23
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return

 ! control options: read these first as they determine which
 ! dependent options below are required
 call read_inopt(in_shape,'in_shape',db,errcount=nerr,min=0,max=1)
 call read_inopt(in_orbit,'in_orbit',db,errcount=nerr,min=0,max=1)
 call read_inopt(rho_mode,'rho_mode',db,errcount=nerr,min=0)
 call read_inopt(cloud_control_mode,'cloud_control_mode',db,errcount=nerr,min=0,max=2)
 call read_inopt(add_turbulence,'add_turbulence',db,errcount=nerr,min=0,max=1)
 call read_inopt(sym_infall,'sym_infall',db,errcount=nerr,min=0,max=3)

 call read_inopt(n_add,'n_add',db,errcount=nerr,min=0)
 call read_inopt(in_mass,'in_mass',db,errcount=nerr,min=0.)
 call read_inopt(m_gas,'m_gas',db,errcount=nerr,min=0.)
 call read_inopt(r_in,'r_in',db,errcount=nerr,min=0.)
 if (in_shape==1) call read_inopt(r_a,'r_a',db,errcount=nerr,min=0.)

 call read_inopt(r_slope,'r_slope',db,errcount=nerr,min=0.)
 if (abs(r_slope) > 0.) call read_inopt(r_soft,'r_soft',db,errcount=nerr,min=0.)

 call read_inopt(r_init,'r_init',db,errcount=nerr,min=0.)

 if (in_orbit==0) call read_inopt(r_close,'r_close',db,errcount=nerr,min=0.)
 if (in_orbit==1) then
    call read_inopt(v_inf,'v_inf',db,errcount=nerr,min=0.)
    call read_inopt(b_frac,'b_frac',db,errcount=nerr,min=0.)
    call read_inopt(b,'b',db,errcount=nerr,min=0.)
    call read_inopt(ecc,'ecc',db,errcount=nerr,min=0.)
    call read_inopt(r_close,'r_close',db,errcount=nerr,min=0.)
 endif

 call read_inopt(incx,'incx',db,errcount=nerr)
 call read_inopt(incy,'incy',db,errcount=nerr)
 call read_inopt(incz,'incz',db,errcount=nerr)

 if (add_turbulence==1) then
    call read_inopt(rms_mach,'rms_mach',db,errcount=nerr,min=0.)
    call read_inopt(tfact,'tfact',db,errcount=nerr,min=0.)
 endif

 call read_inopt(ieos_infall,'ieos_infall',db,errcount=nerr)
 if (ieos_infall /= 6 .and. ieos_infall /= 14) then
    print*,'ERROR: ieos_infall must be 6 or 14'
    nerr = nerr + 1
 endif
 if (ieos_infall==6) call read_inopt(isink,'isink',db,errcount=nerr,min=1)

 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine init_moddump()
 use eos, only:isink

 isink = 1
end subroutine init_moddump

end module moddump
