!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Adds either a sphere or cylinder of material to infall
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, io, part, partinject,
!   physcon, prompting, set_sphere, vectorutils, units
!
 implicit none

 integer,parameter :: nr = 200
 real              :: r_slope = 0.0
 real              :: r_soft = 100.0

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use partinject, only:add_or_update_particle
 use part,       only:igas,isdead_or_accreted,xyzmh_ptmass
 use units,      only:udist,utime,get_G_code
 use io,         only:id,master,fatal
 use spherical,  only:set_sphere,set_ellipse
 use stretchmap, only:rho_func
 use physcon,    only:pi
 use vectorutils,only:rotatevec
 use prompting,      only:prompt
 use centreofmass,   only:reset_centreofmass,get_total_angular_momentum
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, dimension(:,:), allocatable :: xyzh_add,vxyzu_add(:,:)
 integer :: in_shape,in_orbit,ipart,i,n_add,np
 integer(kind=8) :: nptot
 integer, parameter :: iunit = 23
 real    :: r_close,in_mass,hfact,pmass,delta,r_init,r_in,r_a,inc,big_omega
 real    :: v_inf,b,b_frac,theta_def,b_crit,a,ecc
 real    :: vp(3), xp(3), rot_axis(3), rellipsoid(3)
 real    :: dma,n0,pf,m0,x0,y0,z0,r0,vx0,vy0,vz0,mtot,tiny_number
 real    :: unit_velocity, G
 logical :: lrhofunc
 procedure(rho_func), pointer :: prhofunc

 r_close = 100.
 in_mass = 0.01
 r_in = 100.0
 r_init = 1000.0
 in_orbit = 1
 in_shape = 0
 r_slope = 0.0
 inc = 0.
 big_omega = 0.
 tiny_number = 1e-4
 lrhofunc = .false.
 v_inf = 1.0
 theta_def = 90.0
 b_frac = 1.0

 ! udist default is cm
 unit_velocity = udist/utime ! cm/s
 G = get_G_code()

 ! Gas particle properties
 pmass = massoftype(igas)

 #ifdef GRAVITY
   write(*,*), "Disc self-gravity is on. Including disc mass in cloud orbit calculation."
   mtot=sum(xyzmh_ptmass(4,:)) + npartoftype(igas)*massoftype(igas)
 #else
   mtot=sum(xyzmh_ptmass(4,:))
 #endif

 ! Prompt user for infall material shape
 call prompt('Enter the infall material shape (0=sphere, 1=ellipse, 2=cylinder)',in_shape,0,2)

 if (in_shape == 2) then
   write(*,*), "Cylinder not yet implemented."
   stop
 endif

 if (in_shape == 0) then
 call prompt('Enter radius of shape:', r_in, 0.1)
 elseif (in_shape == 1) then
   call prompt('Enter semi-minor axis of ellipse:', r_in, 0.1)
   call prompt('Enter semi-major axis of ellipse:', r_a, 0.1)
   rellipsoid(1) = r_in
   rellipsoid(2) = r_a
   rellipsoid(3) = r_in
 endif
 call prompt('Enter infall mass in Msun:', in_mass, 0.01)
 call prompt('Enter value of power-law density along radius:', r_slope, 0.0)
 call prompt('Enter initial radial distance in au:', r_init, 0.0)


 if (r_slope > tiny_number) then
    prhofunc => rhofunc
    lrhofunc = .true.
    call prompt('Enter softening radius:', r_soft, 0.1)
 endif

 ! Prompt user for the infall material orbit
 call prompt('Enter orbit type (0=bound, 1=parabolic, 2=hyperbolic)', in_orbit,0)

 if (in_orbit == 0) then
   write(*,*), "Bound orbit not yet implemented."
   stop
 endif

 if (in_orbit == 1) then
   print*, "Parabolic orbit"
   call prompt('Enter closest approach in au:', r_close, 0.)
endif

if (in_orbit == 2) then
  write(*,*), "Hyperbolic orbit, see Dullemond+2019 for parameter definitions."
  call prompt('Enter cloud velocity at infinity, v_inf, in km/s:', v_inf, 0.0)

  v_inf = v_inf * (100 * 1000) ! to cm/s
  v_inf =  v_inf / unit_velocity ! Change to code units
  b_crit = mtot * G / v_inf**2
  write(*,*), "Critical impact parameter, b_crit, is ", b_crit, " au"

  call prompt('Enter impact parameter b as a ratio of b_crit:', b_frac, 0.0)
  b = b_frac * b_crit
  ecc = sqrt(1 + b**2/b_crit**2)
  r_close = b * sqrt((ecc-1)/(ecc+1))
  write(*,*), "Eccentricity of the cloud is ", ecc
  write(*,*), "Closest approach of cloud center will be ", r_close, " au."
endif

! Incline the infall
write(*,*), "Notation: 0 degrees is prograde on xy-plane. 180 degrees is retrograde."
call prompt('Enter inclination of infall:', inc, 0., 360.)
call prompt('Enter position angle of ascending node:', big_omega, 0., 360.)

 if (in_orbit == 1) then
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
 elseif (in_orbit == 2) then
   ! Dullemond+2019
   ! Initial position is x=r_init and y=b (impact parameter)
   x0 = r_init
   y0 = b
   z0 = 0.0
   xp = (/x0, y0, z0/)
endif
  write(*,*), "Initial cloud centre is: ", xp

 ! Number of injected particles is given by existing particle mass and total
 ! added disc mass

 n_add = int(in_mass/pmass)
 write(*,*), "Number of particles that will be added ", n_add
 allocate(xyzh_add(4,n_add),vxyzu_add(4,n_add))
 hfact = 1.2
 delta = 1.0 ! no idea what this is
 nptot = n_add + npartoftype(igas)
 np = 0
 if (in_shape == 0) then
    if (lrhofunc) then
      call set_sphere('random',id,master,0.,r_in,delta,hfact,np,xyzh_add,xyz_origin=xp,&
        np_requested=n_add, nptot=nptot,rhofunc=prhofunc)
    else
       call set_sphere('random',id,master,0.,r_in,delta,hfact,np,xyzh_add,xyz_origin=xp,&
          np_requested=n_add, nptot=nptot)
    endif
 write(*,*), "The sphere has been succesfully initialised."
 elseif (in_shape == 1) then
   !set_ellipse(lattice,id,master,r_ellipsoid,delta,hfact,xyzh,np,nptot,np_requested,mask)
   call set_ellipse('random',id,master,rellipsoid,delta,hfact,xyzh_add,np,xyz_origin=xp,&
     np_requested=n_add, nptot=nptot)
   write(*,*), "The ellipse has been succesfully initialised."
 endif


 ! Set up velocities
 if (in_orbit == 1) then

    !--perturber initial velocity
    r0  = sqrt(x0**2+y0**2+z0**2)
    vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
    vy0 = -(x0/r0)*sqrt(mtot/pf)
    vz0 = 0.0
    vp  = (/vx0,vy0,vz0/)

    if (in_shape == 0) then
      ! Initiate initial velocity of the particles in the shape
      vxyzu_add(1, :) = vx0
      vxyzu_add(2, :) = vy0
      vxyzu_add(3, :) = vz0

    elseif (in_shape == 1) then
      do i=1,n_add
        x0 = xyzh_add(1, i)
        y0 = xyzh_add(2, i)
        z0 = xyzh_add(3, i)

        r0  = sqrt(x0**2+y0**2+z0**2)
        vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
        vy0 = -(x0/r0)*sqrt(mtot/pf)
        vz0 = 0.0

        vxyzu_add(1, i) = vx0
        vxyzu_add(2, i) = vy0
        vxyzu_add(3, i) = vz0
      enddo
    endif

    vxyzu_add(4, :) = vxyzu(4, 1)

  elseif (in_orbit == 2) then
    ! Dullemond+2019
    ! Initial velocity, all initially in x direction
    a = -mtot/v_inf**2
    vx0 = sqrt(mtot*(2/r_init - 1/a))
    vy0 = 0.0
    vz0 = 0.0
    vp = (/vx0, vy0, vz0/)
    vxyzu_add(1, :) = vx0
    vxyzu_add(2, :) = vy0
    vxyzu_add(3, :) = vz0
    vxyzu_add(4, :) = vxyzu(4, 1)
  endif

  write(*,*), "Initial velocity of cloud centre is ", vp

 ! Now rotate and add those new particles to existing disc
 ! Default is prograde orbit
 inc = pi-inc*pi/180.
 big_omega = big_omega*pi/180.

 ! Convention: clock-wise rotation in the zx-plane
 rot_axis = (/sin(big_omega),-cos(big_omega),0./)
 ipart = npart ! The initial particle number (post shuffle)
 do i = 1,n_add
    ! Rotate particle to correct position and velocity
    call rotatevec(xyzh_add(1:3,i),rot_axis,inc)
    call rotatevec(vxyzu_add(1:3,i), rot_axis,inc)

    ! Add the particle
    ipart = ipart + 1
    call  add_or_update_particle(igas, xyzh_add(1:3,i), vxyzu_add(1:3,i), xyzh_add(4,i), &
                         vxyzu_add(4,i), ipart, npart, npartoftype, xyzh, vxyzu)
 enddo

 write(*,*),  " ###### Added infall successfully ###### "

 deallocate(xyzh_add,vxyzu_add)

 return
end subroutine modify_dump

real function rhofunc(r)
 real, intent(in) :: r

 rhofunc = 1./(r + r_soft)**(r_slope)

end function rhofunc

end module moddump
