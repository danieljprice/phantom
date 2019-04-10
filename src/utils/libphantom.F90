! $Id$
!----------------------------------------------------------------
! These subroutines are part of the Phantom SPH code
! Use is by specific, written permission of the author
! (c) 2013 Daniel Price and Steven Toupin
!----------------------------------------------------------------
!+
!  Allows phantom to be driven from an external application
!+
!----------------------------------------------------------------

!
! Add gas or boundary particle
!
subroutine inject_or_update_particle(particle_number, mass, position, velocity, h, u, boundary)
 use partinject, only:add_or_update_particle
 use part,       only:igas,iboundary,npart,npartoftype,xyzh,vxyzu
 use units, only: umass, udist, utime
 implicit none
 integer, intent(in) :: particle_number
 double precision, intent(in) :: mass, position(3), velocity(3), h, u
 logical, intent(in) :: boundary
 !logical :: nodisabled

 integer :: i, itype, ierr

 if (boundary) then
    itype = iboundary
 else
    itype = igas
 endif

 !nodisabled = .false.

 !call get_npart(npart, nodisabled)
 call set_part_mass(mass, ierr)
 ! npart = npart + 1
 ! If ierr = 0, it means mass is set successfully.
 ! If it is 1, it means the mass could not be set, because there were
 ! already particles in the simulation.

 call add_or_update_particle(itype,position(:)/udist,velocity(:)/(udist/utime),&
     h/udist,u/(udist**2/utime**2),particle_number,npart,npartoftype,xyzh,vxyzu)

end subroutine

!
! Inject gas or boundary particles
!
subroutine inject_or_update_particles(ifirst, n, position, velocity, h, u, boundary)
 use partinject, only:add_or_update_particle
 use part,       only:igas,iboundary,npart,npartoftype,xyzh,vxyzu
 use units, only: umass, udist, utime
 implicit none
 integer, intent(in) :: ifirst, n
 double precision, intent(in) :: position(3,n), velocity(3,n), h(n), u(n)
 logical, intent(in) :: boundary

 integer :: i, itype

 if (boundary) then
    itype = iboundary
 else
    itype = igas
 endif
 do i=1,n
    call add_or_update_particle(itype,position(:,i)/udist,velocity(:,i)/(udist/utime),&
     h(i)/udist,u(i)/(udist**2/utime**2),ifirst+i-1,npart,npartoftype,xyzh,vxyzu)
 enddo
end subroutine

!
! Inject sink particle
!
subroutine inject_or_update_sink_particle(sink_number, position, velocity, mass, radius)
 use partinject, only:add_or_update_sink
 use units, only: umass, udist, utime
 implicit none
 integer, intent(in) :: sink_number
 double precision, intent(in) :: position(3), velocity(3), mass, radius

 call add_or_update_sink(position/udist,velocity/(udist/utime),radius/udist,mass/umass,sink_number)
end subroutine

!
! Inject sphere
!
subroutine inject_or_update_sphere(ifirst, resolution, center, radius, center_velocity, expansion_velocity, h, u, angles, boundary)
 use icosahedron, only: compute_matrices, compute_corners, pixel2vector
 implicit none
 integer, intent(in) :: ifirst, resolution
 double precision, intent(in) :: center(3), radius
 double precision, intent(in) :: center_velocity(3), expansion_velocity
 double precision, intent(in) :: h, u, angles(3)
 logical, intent(in) :: boundary

 real :: R(0:19,3,3), v(0:11,3)
 integer :: i, N
 real :: base_sphere(3,40*resolution*(resolution-1)+12), sphere(3,40*resolution*(resolution-1)+12)
 real :: particles(3,40*resolution*(resolution-1)+12), velocities(3,40*resolution*(resolution-1)+12)
 real :: h_part(40*resolution*(resolution-1)+12), u_part(40*resolution*(resolution-1)+12)
 real :: rot_m(3,3)
 real :: c_x, s_x, c_y, s_y, c_z, s_z

 ! Construct base uniform sphere
 N = 40*resolution*(resolution-1)+12
 call compute_matrices(R)
 call compute_corners(v)
 do i=0,N-1
    call pixel2vector(i, resolution, R, v, base_sphere(:,i+1))
 enddo

 ! Build rotation matrix
 c_x = cos(angles(1))
 s_x = sin(angles(1))
 c_y = cos(angles(2))
 s_y = sin(angles(2))
 c_z = cos(angles(3))
 s_z = sin(angles(3))
 rot_m(1,1) = c_y*c_z
 rot_m(1,2) = -c_y*s_z
 rot_m(1,3) = -s_y
 rot_m(2,1) = -s_x*s_y*c_z + c_x*s_z
 rot_m(2,2) = s_x*s_y*s_z + c_x*c_z
 rot_m(2,3) = -s_x*c_y
 rot_m(3,1) = c_x*s_y*c_z + s_x*s_z
 rot_m(3,2) = -c_x*s_y*s_z + s_x*c_z
 rot_m(3,3) = c_x*c_y

 ! Rotate base sphere
 sphere(1,:) = base_sphere(1,:)*rot_m(1,1) + base_sphere(2,:)*rot_m(1,2) + base_sphere(3,:)*rot_m(1,3)
 sphere(2,:) = base_sphere(1,:)*rot_m(2,1) + base_sphere(2,:)*rot_m(2,2) + base_sphere(3,:)*rot_m(2,3)
 sphere(3,:) = base_sphere(1,:)*rot_m(3,1) + base_sphere(2,:)*rot_m(3,2) + base_sphere(3,:)*rot_m(3,3)

 ! Build arrays
 do i=1,N
    particles(:,i) = center + sphere(:,i)*radius
 enddo
 do i=1,N
    velocities(:,i) = center_velocity + sphere(:,i)*expansion_velocity
 enddo
 h_part = h
 u_part = u

 ! Inject particles
 call inject_or_update_particles(ifirst, N, particles, velocities, h_part, u_part, boundary)
end subroutine

!
! Plot density cross section
!
subroutine plot_rho_xysec(z, xmin, ymin, npx, npy, dx, dy, npart, massofgas, hfact, xyzh, pixmap)
 use libphantomsplash, only:interpolate3D_fastxsec
 implicit none
 double precision, intent(in) :: z, xmin, ymin
 integer, intent(in) :: npx, npy
 double precision, intent(in) :: dx, dy
 integer, intent(in) :: npart
 double precision, intent(in) :: massofgas, hfact, xyzh(4,npart)
 double precision, intent(out) :: pixmap(npx, npy)

 real :: w(npart), datsmooth(npx, npy)

 w(:) = massofgas/xyzh(4,:)**3

 call interpolate3D_fastxsec(xyzh,1.,w,npart,&
     xmin,ymin,z,datsmooth,npx,npy,dx,dy,.false.)
 pixmap = dble(datsmooth)
end subroutine

!
! Plot log density cross section
!
subroutine plot_log_rho_xysec(z, xmin, ymin, npx, npy, dx, dy, npart, massofgas, hfact, xyzh, pixmap)
 use libphantomsplash, only:interpolate3D_fastxsec
 implicit none
 double precision, intent(in) :: z, xmin, ymin
 integer, intent(in) :: npx, npy
 double precision, intent(in) :: dx, dy
 integer, intent(in) :: npart
 double precision, intent(in) :: massofgas, hfact, xyzh(4,npart)
 double precision, intent(out) :: pixmap(npx, npy)

 double precision :: v, minpositive, logminpositive
 integer :: ix, iy

 call plot_rho_xysec(z, xmin, ymin, npx, npy, dx, dy, npart, massofgas, hfact, xyzh, pixmap)

 ! 1st pass: look for minimum positive value
 minpositive = huge(minpositive)
 do iy=1,npy
    do ix=1,npx
       v = pixmap(ix,iy)
       if (v  >  0.) then
          minpositive = min(minpositive, v)
       endif
    enddo
 enddo
 logminpositive = log10(minpositive)

 ! 2nd pass: set minimum value where pixmap is negative and log
 do iy=1,npy
    do ix=1,npx
       v = pixmap(ix,iy)
       if (v  >  0.) then
          pixmap(ix,iy) = log10(v)
       else
          pixmap(ix,iy) = logminpositive
       endif
    enddo
 enddo
end subroutine


!
! Get the boundaries
!
subroutine get_boundaries(xmin, xmax, ymin, ymax, zmin, zmax)
 use part, only:xyzh,npart
 use units, only:udist
 implicit none
 double precision, intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax

 integer :: i
 real :: x, y, z, h

 xmin = xyzh(1,1)
 xmax = xyzh(1,1)
 ymin = xyzh(2,1)
 ymax = xyzh(2,1)
 zmin = xyzh(3,1)
 zmax = xyzh(3,1)
 do i = 2, npart
    h = xyzh(4,i)
    if (h  >  0.) then
       x = xyzh(1,i)
       y = xyzh(2,i)
       z = xyzh(3,i)
       xmin = min(x,xmin)
       xmax = max(x,xmax)
       ymin = min(y,ymin)
       ymax = max(y,ymax)
       zmin = min(z,zmin)
       zmax = max(z,zmax)
    endif
 enddo
 xmin = xmin*udist
 xmax = xmax*udist
 ymin = ymin*udist
 ymax = ymax*udist
 zmin = zmin*udist
 zmax = zmax*udist
end subroutine

!
! Initialize Phantom
!
subroutine code_init()
 use initial, only:initialise
 implicit none

 call initialise()
end subroutine

!
! Set default parameters
!
subroutine set_defaults()
 use options, only: set_default_options
 implicit none

 call set_default_options()
end subroutine


!
! Initialize a simulation
!
subroutine init(len_infile, infile, logfile, evfile, dumpfile, tmax_in, dtmax_in)
 use initial, only:startrun
 use io, only:set_io_unit_numbers
 use evolvesplit, only:evol_init
 use timestep, only:tmax,dt,dtmax
 use units, only:utime
 implicit none
 integer,            intent(in)  :: len_infile
 character(len=len_infile) :: infile
 character(len=120), intent(out) :: logfile,evfile,dumpfile
 double precision, intent(in) :: tmax_in, dtmax_in

 call set_io_unit_numbers
#ifndef AMUSE
 if (index(infile,'.in')==0) then
    infile = trim(infile)//'.in'
 endif
#endif
 call startrun(infile,logfile,evfile,dumpfile)
 ! Overrides tmax and dtmax
 tmax = tmax_in/utime
 dtmax = dtmax_in/utime
 dt = dtmax
 call evol_init()
end subroutine

!
! Overrides tmax and dtmax
!
subroutine override_tmax_dtmax(tmax_in, dtmax_in)
 use timestep,only:tmax,dtmax
 use units,only:utime
 use evolvesplit, only:evol_init
 implicit none
 double precision, intent(in) :: tmax_in, dtmax_in

 tmax = tmax_in/utime
 dtmax = dtmax_in/utime
 call evol_init()
end subroutine

!
! Set dt
!
subroutine set_dt(dt_in)
 use timestep,only:dt
 use units,only:utime
 implicit none
 double precision, intent(in) :: dt_in

 dt = dt_in/utime
end subroutine

!
! Finalize a simulation
!
subroutine finalize()
 use initial, only:endrun
 implicit none

 call endrun()
end subroutine

!
! Initialize a timestep
!
subroutine init_step_wrapper()
 use evolvesplit, only:init_step
 implicit none
 call init_step()
end subroutine

!
! Finalize a timestep
!
subroutine finalize_step_wrapper(len_infile, infile, len_logfile, logfile, &
                                 len_evfile, evfile, len_dumpfile, dumpfile)
 use evolvesplit, only:finalize_step
 implicit none
 integer,                     intent(in)    :: len_infile, len_logfile, len_evfile, len_dumpfile
 character(len=len_infile),   intent(in)    :: infile
 character(len=len_logfile),  intent(inout) :: logfile
 character(len=len_evfile),   intent(inout) :: evfile
 character(len=len_dumpfile), intent(inout) :: dumpfile

 call finalize_step(infile, logfile, evfile, dumpfile)
end subroutine

!
! Calculate new timestep
!
subroutine calculate_timestep()
 use timestep, only:time,tmax,dtmax,dt,dtforce,dtcourant,dterr
 implicit none
 real :: dtexact
 dtexact = tmax - time + epsilon(dtmax)
 if (dtexact <= epsilon(dtmax) .or. dtexact >= (1.0-1e-8)*dtmax ) then
     dtexact = dtmax + epsilon(dtmax)
 endif
 dt = min(dtforce,dtcourant,dterr,dtmax+epsilon(dtmax),dtexact)
end subroutine

!
! Get stepping and timing information
!
subroutine get_time_info(time_out, tmax_out, nsteps_out, nmax_out, dt_out)
 use timestep,         only:time,tmax,nmax,nsteps,dt
 use units,            only:utime
 implicit none
 double precision, intent(out) :: time_out, tmax_out, dt_out
 integer, intent(out) :: nsteps_out, nmax_out

 time_out = dble(time*utime)
 tmax_out = dble(tmax*utime)
 nsteps_out = nsteps
 nmax_out = nmax
 dt_out = dble(dt*utime)
end subroutine

!
! Advance the simulation in time
!
subroutine step_wrapper()
 use timestep,         only:time,dt,dtextforce
 use step_lf_global, only:step
 use part, only:npart
#ifdef IND_TIMESTEPS
 use timestep_ind, only:nactive
#endif
 implicit none
 real :: dtnew
#ifndef IND_TIMESTEPS
 integer :: nactive

 nactive= npart
#endif

 call step(npart,nactive,time,dt,dtextforce,dtnew)
end subroutine

!
! Call inject particles routine
!
subroutine inject_particles_wrapper()
#ifdef INJECT_PARTICLES
 use timestep,    only:time
 use evolvesplit, only:dtlast
 use inject,      only:inject_particles
 implicit none

 call inject_particles(time,dtlast)
#endif
end subroutine

!
! Read a dump file
! Inspired by subroutine from phantomanalysis.f90
!
subroutine read_dump_wrapper(len_dumpfile, dumpfile, headeronly, time_dp, hfact_dp, massofgas_dp, ierr)
 use part, only:hfact,massoftype,igas
 use io, only:iprint,idisk1,set_io_unit_numbers
 use readwrite_dumps, only:read_dump
 use units, only:umass
 implicit none
 integer,                     intent(in) :: len_dumpfile
 character(len=len_dumpfile), intent(in) :: dumpfile
 logical,                     intent(in) :: headeronly
 double precision, intent(out) :: time_dp, hfact_dp, massofgas_dp
 integer, intent(out) :: ierr

 real :: time

 call set_io_unit_numbers
 iprint = 6
 call read_dump(dumpfile,time,hfact,idisk1,iprint,0,1,ierr,headeronly)
 time_dp = dble(time)
 hfact_dp = dble(hfact)
 massofgas_dp = dble(massoftype(igas)*umass)
end subroutine

!
! Get the mass of gas particles
!
subroutine get_massofgas(massofgas)
 use part, only:massoftype,igas
 use units, only:umass
 implicit none
 double precision, intent(out) :: massofgas

 massofgas = dble(massoftype(igas)*umass)
end subroutine

!
! Get hfact
!
!subroutine get_hfact(hfact_out)
! use part, only:hfact
! implicit none
! double precision, intent(out) :: hfact_out
!
! hfact_out = dble(hfact)
!end subroutine

!
! Get npart
!
subroutine get_npart(npart_out, nodisabled)
 use part, only:npart,xyzh
 implicit none
 integer, intent(out) :: npart_out
 logical, intent(in)  :: nodisabled

 integer :: i

 if (nodisabled) then
    npart_out = 0
    do i=1,npart
       if (xyzh(4,i)  >  0.) then
          npart_out = npart_out + 1
       endif
    enddo
 else
    npart_out = npart
 endif
end subroutine

!
! Get xyzh
!
subroutine get_part_xyzh(npart_in, part_xyzh, nodisabled, ierr)
 use part, only:npart,xyzh
 use units, only:udist
 implicit none
 integer, intent(in) :: npart_in
 double precision, dimension(4,npart_in), intent(out) :: part_xyzh
 logical, intent(in)  :: nodisabled
 integer, intent(out) :: ierr

 integer :: i, j, n

 if (nodisabled) then
    n = 0
    do i = 1,npart
       if (xyzh(4,i)  >  0.) then
          n = n + 1
          if (n  >  npart_in) then
             ierr = 1
             exit
          endif
          do j=1,4
             part_xyzh(j,n) = dble(xyzh(j,i)*udist)
          enddo
       endif
    enddo
 else
    if (npart_in == npart) then
       do i = 1,npart
          do j = 1,4
             part_xyzh(j,i) = dble(xyzh(j,i)*udist)
          enddo
       enddo
       ierr = 0
    else
       ierr = 1
    endif
 endif
end subroutine


!
! Get vxyz
!
subroutine get_part_vxyz(npart_in, part_vxyz, nodisabled, ierr)
 use part, only:npart,xyzh,vxyzu
 use units, only:udist,utime
 implicit none
 integer, intent(in) :: npart_in
 double precision, dimension(3,npart_in), intent(out) :: part_vxyz
 logical, intent(in)  :: nodisabled
 integer, intent(out) :: ierr

 integer :: i, n

 if (nodisabled) then
    n = 0
    do i=1,npart
       if (xyzh(4,i)  >  0.) then
          n = n + 1
          if (n  >  npart_in) then
             ierr = 1
             exit
          endif
          part_vxyz(1:3,n) = dble(vxyzu(1:3,i)*udist/utime)
       endif
    enddo
 else
    if (npart_in == npart) then
       part_vxyz(1:3,1:npart) = dble(vxyzu(1:3,1:npart)*udist/utime)
       ierr = 0
    else
       ierr = 1
    endif
 endif
end subroutine

!
! Get magnetic field, if possible
!
subroutine get_part_bxyz(npart_in, part_bxyz, nodisabled, ierr)
 use part,  only:npart,xyzh,Bxyz,mhd
 use units, only:unit_Bfield
 implicit none
 integer, intent(in) :: npart_in
 double precision, dimension(3,npart_in), intent(out) :: part_bxyz
 logical, intent(in)  :: nodisabled
 integer, intent(out) :: ierr
 integer :: i, n

 if (mhd) then
    if (nodisabled) then
       n = 0
       do i=1,npart
          if (xyzh(4,i)  >  0.) then
             n = n + 1
             if (n  >  npart_in) then
                ierr = 1
                exit
             endif
             part_bxyz(1:3,n) = dble(Bxyz(1:3,i)*unit_Bfield)
          endif
       enddo
    else
       if (npart_in == npart) then
          part_bxyz(1:3,1:npart) = dble(Bxyz(1:3,1:npart)*unit_Bfield)
          ierr = 0
       else
          ierr = 1
       endif
    endif
 else
    ierr = 2
 endif
end subroutine

!
! Get u, if possible
!
subroutine get_part_u(npart_in, part_u, nodisabled, ierr)
 use part, only:npart,xyzh,vxyzu,maxvxyzu
 use units, only:udist,utime
 implicit none
 integer, intent(in) :: npart_in
 double precision, dimension(npart_in), intent(out) :: part_u
 logical, intent(in)  :: nodisabled
 integer, intent(out) :: ierr
 integer :: i, n

 if (maxvxyzu == 4) then
    if (nodisabled) then
       n = 0
       do i=1,npart
          if (xyzh(4,i)  >  0.) then
             n = n + 1
             if (n  >  npart_in) then
                ierr = 1
                exit
             endif
             part_u(n) = vxyzu(4,i)
          endif
       enddo
    else
       if (npart_in == npart) then
          part_u(1:npart) = dble(vxyzu(4,1:npart)*udist**2/utime**2)
       else
          ierr = 1
       endif
    endif
 else
    ierr = 2
 endif
end subroutine

!
! Get temperature, if stored
!
subroutine get_part_temp(npart_in, part_temp, nodisabled, ierr)
 use part,  only:npart,xyzh,temperature,store_temperature
 implicit none
 integer, intent(in) :: npart_in
 double precision, dimension(npart_in), intent(out) :: part_temp
 logical, intent(in)  :: nodisabled
 integer, intent(out) :: ierr
 integer :: i, n

 if (store_temperature) then
    if (nodisabled) then
       n = 0
       do i = 1, npart
          if (xyzh(4,i) > 0.) then
             n = n + 1
             if (n > npart_in) then
                ierr = 1
                exit
             endif
             part_temp(n) = temperature(i)
          endif
       enddo
    else
       if (npart_in == npart) then
          part_temp(1:npart) = dble(temperature(1:npart))
       else
          ierr = 1
       endif
    endif
 else
    ierr = 2
 endif
end subroutine

!
! Get nptmass
!
subroutine get_nptmass(nptmass_out)
 use part, only:nptmass
 implicit none
 integer, intent(out) :: nptmass_out

 nptmass_out = nptmass
end subroutine

!
! Get ptmass properties: location, accretion radius, mass
!
subroutine get_ptmass_xyzmh(nptmass_in, ptmass_xyzmh_out, ierr)
 use part, only:nptmass,xyzmh_ptmass
 use units, only:udist,umass
 implicit none
 integer, intent(in) :: nptmass_in
 double precision, dimension(5,nptmass_in) :: ptmass_xyzmh_out
 integer, intent(out) :: ierr

 integer :: i

 if (nptmass_in == nptmass) then
    do i=1,nptmass
       ptmass_xyzmh_out(:,i) = xyzmh_ptmass(1:5,i)*(/ udist, udist, udist, umass, udist /)
    enddo
 else
    ierr = 1
 endif
end subroutine

!
! Get ptmass velocities
!
subroutine get_ptmass_vxyz(nptmass_in, ptmass_vxyz_out, ierr)
 use part, only:nptmass,vxyz_ptmass
 use units, only:udist,utime
 implicit none
 integer, intent(in) :: nptmass_in
 double precision, dimension(3,nptmass_in) :: ptmass_vxyz_out
 integer, intent(out) :: ierr

 if (nptmass_in == nptmass) then
    ptmass_vxyz_out(1:3,1:nptmass) = vxyz_ptmass(1:3,1:nptmass)*(udist/utime)
 else
    ierr = 1
 endif
end subroutine

!
! Get ptmass spin
!
subroutine get_ptmass_spinxyz(nptmass_in, ptmass_spinxyz, ierr)
 use part, only:nptmass,xyzmh_ptmass,ispinx,ispiny,ispinz
 use units, only:udist,umass,utime
 implicit none
 integer, intent(in) :: nptmass_in
 double precision, dimension(3,nptmass_in) :: ptmass_spinxyz
 integer, intent(out) :: ierr

 integer :: i

 if (nptmass_in == nptmass) then
    do i=1,nptmass
       ptmass_spinxyz(1,i) = xyzmh_ptmass(ispinx,i)*(udist**2*umass/utime)
       ptmass_spinxyz(2,i) = xyzmh_ptmass(ispiny,i)*(udist**2*umass/utime)
       ptmass_spinxyz(3,i) = xyzmh_ptmass(ispinz,i)*(udist**2*umass/utime)
    enddo
 else
    ierr = 1
 endif
end subroutine

!
! Set the mass of particles
!
subroutine set_part_mass(newmass, ierr)
 use part, only:npart,massoftype
 use units, only:umass
 implicit none
 double precision, intent(in) :: newmass
 integer, intent(out) :: ierr

 if (npart == 0) then
    massoftype(:) = newmass/umass
    ierr = 0
 else
    ierr = 1 ! Can only set mass of particles if there is no particle in the simulation
 endif
end subroutine

!
! Get the units of the simulation
!
subroutine get_units(udist_out, umass_out, utime_out, udens_out, umagfd_out)
 use units, only:udist,umass,utime,unit_density,unit_Bfield
 implicit none
 double precision, intent(out) :: udist_out, umass_out, utime_out, udens_out, umagfd_out

 udist_out = udist
 umass_out = umass
 utime_out = utime
 udens_out = unit_density
 umagfd_out = unit_Bfield
end subroutine

!
! Disable particles that lie outside a box
!
subroutine delete_particles_outside_box_wrapper(xmin_in, xmax_in, ymin_in, ymax_in, zmin_in, zmax_in)
 use part, only: delete_particles_outside_box
 use units, only: udist
 implicit none
 double precision, intent(in) :: xmin_in, xmax_in, ymin_in, ymax_in, zmin_in, zmax_in

 double precision :: xmin, xmax, ymin, ymax, zmin, zmax
 xmin = xmin_in/udist
 xmax = xmax_in/udist
 ymin = ymin_in/udist
 ymax = ymax_in/udist
 zmin = zmin_in/udist
 zmax = zmax_in/udist
 call delete_particles_outside_box(xmin, xmax, ymin, ymax, zmin, zmax)
end subroutine

!
! Disable particles that lie outside a sphere
!
subroutine delete_particles_outside_sphere_wrapper(center, radius)
 use part, only: delete_particles_outside_sphere
 use units, only: udist
 implicit none
 double precision, intent(in) :: center(3), radius

 call delete_particles_outside_sphere(center/udist, radius/udist)
end subroutine

!!
!! Below this are AMUSE helper subroutines
!!

!
! Initialize Phantom and set default parameters
!
subroutine amuse_initialize_code()
    use dim, only:maxp,maxp_hard
    use memory, only:allocate_memory
    use physcon, only:pc,solarm
    use units, only:set_units
    implicit none
    call allocate_memory(maxp_hard)
    call code_init()
    call set_defaults()
    call set_units(dist=0.1d0*pc,mass=solarm,G=1.)
end subroutine

subroutine amuse_cleanup_code()
    implicit none
    call finalize()
end subroutine

! New particles
subroutine amuse_new_sph_particle(i, mass, x, y, z, vx, vy, vz, u, h)
    use part, only:igas,npart,npartoftype,xyzh,vxyzu,massoftype
    use partinject, only:add_or_update_particle
    implicit none
    integer :: n, i, itype
    double precision :: mass, x, y, z, vx, vy, vz, u, h
    double precision :: position(3), velocity(3)
  
    itype = igas
    i = npart + 1
    position(1) = x
    position(2) = y
    position(3) = z
    velocity(1) = vx
    velocity(2) = vy
    velocity(3) = vz
    if (npartoftype(itype) == 0) then
        massoftype(itype) = mass
    endif
    call add_or_update_particle(itype,position,velocity,h, &
        u,i,npart,npartoftype,xyzh,vxyzu)
end subroutine

subroutine amuse_new_dm_particle(i, mass, x, y, z, vx, vy, vz, radius)
    use part, only:idarkmatter,npart,npartoftype,xyzh,vxyzu,massoftype
    use partinject, only:add_or_update_particle
    implicit none
    integer :: n, i, itype
    double precision :: mass, x, y, z, vx, vy, vz, radius, u
    double precision :: position(3), velocity(3)
  
    u=0
    itype = idarkmatter
    i = npart + 1
    position(1) = x
    position(2) = y
    position(3) = z
    velocity(1) = vx
    velocity(2) = vy
    velocity(3) = vz
    if (npartoftype(itype) == 0) then
        massoftype(itype) = mass
    endif

    call add_or_update_particle(itype,position,velocity,radius, &
        u,i,npart,npartoftype,xyzh,vxyzu)
end subroutine

subroutine amuse_new_sink_particle(i, mass, x, y, z, vx, vy, vz, radius)
    use part, only:npart
    use partinject, only:add_or_update_sink
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz, radius
    double precision :: position(3), velocity(3)
  
    i = npart + 1
    position(1) = x
    position(2) = y
    position(3) = z
    velocity(1) = vx
    velocity(2) = vy
    velocity(3) = vz
    call add_or_update_sink(position,velocity,radius,mass,i)
end subroutine

subroutine amuse_delete_particle(i)
    use part, only:kill_particle
    call kill_particle(i)
end subroutine

subroutine amuse_get_potential_energy(epot_out)
    use energies, only:epot
    implicit none
    double precision, intent(out) :: epot_out
    epot_out = epot
end subroutine

subroutine amuse_get_kinetic_energy(ekin_out)
    use energies, only:ekin
    implicit none
    double precision, intent(out) :: ekin_out
    ekin_out = ekin
end subroutine

subroutine amuse_get_thermal_energy(etherm_out)
    use energies, only:etherm
    implicit none
    double precision, intent(out) :: etherm_out
    etherm_out = etherm
end subroutine

subroutine amuse_get_time_step(dt_out)
    use timestep, only:dt
    implicit none
    double precision, intent(out) :: dt_out
    dt_out = dt
end subroutine

subroutine amuse_get_number_of_sph_particles(n)
    use part, only:npartoftype,igas
    implicit none
    integer, intent(out) :: n
    logical :: nodisabled
    nodisabled = .false.
    n = npartoftype(igas)
end subroutine

subroutine amuse_get_number_of_particles(n)
    use part, only:npart
    implicit none
    integer, intent(out) :: n
    logical :: nodisabled
    nodisabled = .false.
    ! should look at npartoftype(itype) in part.f90 module instead?, as npart can only increase
    call get_npart(n, nodisabled)
end subroutine

subroutine amuse_get_time(time_out)
    use timestep, only:time
    implicit none
    double precision, intent(out) :: time_out
    time_out = time
end subroutine

subroutine amuse_get_density(i, rho)
    use part, only:rhoh,iphase,massoftype,xyzh
    implicit none
    integer :: i
    double precision :: pmassi
    double precision, intent(out) :: rho
    pmassi = massoftype(abs(iphase(i)))
    rho = rhoh(xyzh(4,i), pmassi)
end subroutine

subroutine amuse_get_pressure(i, p)
    use part, only:rhoh,iphase,massoftype,xyzh
    use eos, only:ieos,equationofstate
    implicit none
    integer :: i, eos_type
    double precision :: pmassi, ponrho, rho, spsound, x, y, z
    double precision, intent(out) :: p
    eos_type = ieos
    pmassi = massoftype(abs(iphase(i)))
    call amuse_get_density(i, rho)
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    call equationofstate(eos_type,ponrho,spsound,rho,x,y,z)
    p = ponrho * rho
end subroutine

subroutine amuse_get_mass(i, part_mass)
    use part, only:iphase,massoftype
    implicit none
    double precision, intent(out) :: part_mass
    integer :: i
    !TODO: Need something different for sinks ("ptmass")
    part_mass = massoftype(abs(iphase(i)))
end subroutine

subroutine amuse_get_state_gas(i, mass, x, y, z, vx, vy, vz, u, h)
    implicit none
    integer :: i
    double precision, intent(inout) :: mass, x, y, z, vx, vy, vz, u, h
    call amuse_get_mass(i, mass)
    call amuse_get_position(i, x, y, z)
    call amuse_get_velocity(i, vx, vy, vz)
    call amuse_get_internal_energy(i, u)
    call amuse_get_smoothing_length(i, h)
end subroutine

subroutine amuse_get_state_dm(i, mass, x, y, z, vx, vy, vz, radius)
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz, radius
    call amuse_get_mass(i, mass)
    call amuse_get_position(i, x, y, z)
    call amuse_get_velocity(i, vx, vy, vz)
    call amuse_get_radius(i, radius)
end subroutine

subroutine amuse_get_position(i, x, y, z)
    use part, only:xyzh
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: x, y, z
    x = xyzh(1, i)
    y = xyzh(2, i)
    z = xyzh(3, i)
end subroutine

subroutine amuse_get_velocity(i, vx, vy, vz)
    use part, only:vxyzu
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: vx, vy, vz
    vx = vxyzu(1, i)
    vy = vxyzu(2, i)
    vz = vxyzu(3, i)
end subroutine

subroutine amuse_get_smoothing_length(i, h)
    use part, only:xyzh
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: h
    h = xyzh(4, i)
end subroutine

subroutine amuse_get_radius(i, radius)
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: radius
    call amuse_get_smoothing_length(i, radius)
end subroutine

subroutine amuse_get_internal_energy(i, u)
    use dim, only:maxvxyzu
    use part, only:vxyzu
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: u

    if (maxvxyzu >= 4) then
        u = vxyzu(4, i)
    else
        u = 0
    endif
end subroutine

subroutine amuse_get_dtmax(dtmax_out)
    use timestep, only:dtmax
    implicit none
    double precision, intent(out) :: dtmax_out
    dtmax_out = dtmax
end subroutine

subroutine amuse_set_time_step(dt_in)
    use timestep, only:dt
    implicit none
    double precision, intent(in) :: dt_in
    dt = dt_in
end subroutine

subroutine amuse_set_dtmax(dtmax_in)
    use timestep, only:dtmax
    implicit none
    double precision, intent(in) :: dtmax_in
    dtmax = dtmax_in
end subroutine

subroutine amuse_set_mass(i, part_mass)
    use part, only:iphase,massoftype
    implicit none
    double precision, intent(in) :: part_mass
    integer :: i
    ! Need to do something different for sinks ("ptmass")
    massoftype(abs(iphase(i))) = part_mass
end subroutine

subroutine amuse_set_state_gas(i, mass, x, y, z, vx, vy, vz, u, h)
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz, u, h
    call amuse_set_mass(i, mass)
    call amuse_set_position(i, x, y, z)
    call amuse_set_velocity(i, vx, vy, vz)
    call amuse_set_internal_energy(i, u)
    call amuse_set_smoothing_length(i, h)
end subroutine

subroutine amuse_set_state_dm(i, mass, x, y, z, vx, vy, vz, radius)
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz, radius
    call amuse_set_mass(i, mass)
    call amuse_set_position(i, x, y, z)
    call amuse_set_velocity(i, vx, vy, vz)
    call amuse_set_radius(i, radius)
end subroutine

subroutine amuse_set_position(i, x, y, z)
    use part, only:xyzh
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: x, y, z
    xyzh(1, i) = x
    xyzh(2, i) = y
    xyzh(3, i) = z
end subroutine

subroutine amuse_set_velocity(i, vx, vy, vz)
    use part, only:vxyzu
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: vx, vy, vz
    vxyzu(1, i) = vx
    vxyzu(2, i) = vy
    vxyzu(3, i) = vz
end subroutine

subroutine amuse_set_smoothing_length(i, h)
    use part, only:xyzh
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: h
    xyzh(4, i) = h
end subroutine

subroutine amuse_set_radius(i, radius)
    implicit none
    integer, intent(in) :: i
    double precision :: radius
    call amuse_set_smoothing_length(i, radius)
end subroutine

subroutine amuse_set_internal_energy(i, u)
    use dim, only:maxvxyzu
    use part, only:vxyzu
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: u
    if (maxvxyzu >= 4) then
        vxyzu(4, i) = u
    endif
end subroutine

subroutine amuse_evolve_model(tmax_in)
    use timestep, only:tmax, time, dt, dtmax, rhomaxnow
    use evolvesplit, only:init_step, finalize_step
    use units, only:utime, udist, umass
    use options, only:rhofinal1
    implicit none
    double precision, intent(in) :: tmax_in
    logical :: maximum_density_reached
    integer :: len_infile
    character(len=120) :: infile, logfile, evfile, dumpfile
    integer :: steps_this_loop

    steps_this_loop = 0
  
    infile = '/dev/null'
    logfile = '/dev/null'
    evfile = '/dev/null'
    dumpfile = '/dev/null'
    
    tmax = tmax_in
    
    timestepping: do while ( &
        (time < tmax) .and. &
        ((rhomaxnow*rhofinal1 < 1.0) .or. (steps_this_loop < 1)) &
    )
        call init_step()
        call calculate_timestep()
        call step_wrapper()
        call finalize_step(infile, logfile, evfile, dumpfile)
        steps_this_loop = steps_this_loop + 1
    enddo timestepping

end subroutine


!
! Setters and getters for parameters
!
subroutine amuse_set_c_courant(C_cour_in)
    use timestep, only:C_cour
    implicit none
    double precision, intent(in) :: C_cour_in
    C_cour = C_cour_in
end subroutine

subroutine amuse_set_c_force(C_force_in)
    use timestep, only:C_force
    implicit none
    double precision, intent(in) :: C_force_in
    C_force = C_force_in
end subroutine

subroutine amuse_set_tolv(tolv_in)
    use timestep, only:tolv
    implicit none
    double precision, intent(in) :: tolv_in
    tolv = tolv_in
end subroutine

subroutine amuse_set_hfact(hfact_in)
    use part, only:hfact
    implicit none
    double precision, intent(in) :: hfact_in
    hfact = hfact_in
end subroutine

subroutine amuse_set_tolh(tolh_in)
    use options, only:tolh
    implicit none
    double precision, intent(in) :: tolh_in
    tolh = tolh_in
end subroutine

subroutine amuse_set_tree_accuracy(tree_accuracy_in)
    use kdtree, only:tree_accuracy
    implicit none
    double precision, intent(in) :: tree_accuracy_in
    tree_accuracy = tree_accuracy_in
end subroutine

subroutine amuse_set_alpha(alpha_in)
    use options, only:alpha
    implicit none
    double precision, intent(in) :: alpha_in
    alpha = alpha_in
end subroutine

subroutine amuse_set_alphamax(alphamax_in)
    use options, only:alphamax
    implicit none
    double precision, intent(in) :: alphamax_in
    alphamax = alphamax_in
end subroutine

subroutine amuse_set_beta(beta_in)
    use options, only:beta
    implicit none
    double precision, intent(in) :: beta_in
    beta = beta_in
end subroutine

subroutine amuse_set_avdecayconst(avdecayconst_in)
    use options, only:avdecayconst
    implicit none
    double precision, intent(in) :: avdecayconst_in
    avdecayconst = avdecayconst_in
end subroutine

subroutine amuse_set_idamp(idamp_in)
    use options, only:idamp
    implicit none
    integer, intent(in) :: idamp_in
    idamp = idamp_in
end subroutine

subroutine amuse_set_ieos(ieos_in)
    use eos, only:ieos
    implicit none
    integer, intent(in) :: ieos_in
    ieos = ieos_in
end subroutine

subroutine amuse_set_mu(mu_in)
    use eos, only:gmw
    implicit none
    double precision, intent(in) :: mu_in
    gmw = mu_in
end subroutine

subroutine amuse_set_rhofinal(rhofinal_in)
    use options, only:rhofinal_cgs, rhofinal1
    use units, only:unit_density
    implicit none
    double precision, intent(in) :: rhofinal_in
    rhofinal_cgs = rhofinal_in * unit_density
    if (rhofinal_cgs > 0.) then
        rhofinal1 = unit_density/rhofinal_cgs
    else
        rhofinal1 = 0.0
    endif
end subroutine

subroutine amuse_set_rho_crit(rho_crit_in)
    use ptmass, only:rho_crit_cgs
    use units, only:unit_density
    implicit none
    double precision, intent(in) :: rho_crit_in
    rho_crit_cgs = rho_crit_in * unit_density
end subroutine

subroutine amuse_set_r_crit(r_crit_in)
    use ptmass, only:r_crit
    implicit none
    double precision, intent(in) :: r_crit_in
    r_crit = r_crit_in
end subroutine

subroutine amuse_set_h_acc(h_acc_in)
    use ptmass, only:h_acc
    implicit none
    double precision, intent(in) :: h_acc_in
    h_acc = h_acc_in
end subroutine

subroutine amuse_set_h_soft_sinkgas(h_soft_sinkgas_in)
    use ptmass, only:h_soft_sinkgas
    implicit none
    double precision, intent(in) :: h_soft_sinkgas_in
    h_soft_sinkgas = h_soft_sinkgas_in
end subroutine

subroutine amuse_set_h_soft_sinksink(h_soft_sinksink_in)
    use ptmass, only:h_soft_sinksink
    implicit none
    double precision, intent(in) :: h_soft_sinksink_in
    h_soft_sinksink = h_soft_sinksink_in
end subroutine

subroutine amuse_set_f_acc(f_acc_in)
    use ptmass, only:f_acc
    implicit none
    double precision, intent(in) :: f_acc_in
    f_acc = f_acc_in
end subroutine

subroutine amuse_set_iexternalforce(iexternalforce_in)
    use options, only:iexternalforce
    implicit none
    integer, intent(in) :: iexternalforce_in
    iexternalforce = iexternalforce_in
end subroutine

subroutine amuse_set_irealvisc(irealvisc_in)
    use viscosity, only:irealvisc
    implicit none
    integer, intent(in) :: irealvisc_in
    irealvisc = irealvisc_in
end subroutine

subroutine amuse_set_shearparam(shearparam_in)
    use viscosity, only:shearparam
    implicit none
    double precision, intent(in) :: shearparam_in
    shearparam = shearparam_in
end subroutine

subroutine amuse_set_bulkvisc(bulkvisc_in)
    use viscosity, only:bulkvisc
    implicit none
    double precision, intent(in) :: bulkvisc_in
    bulkvisc = bulkvisc_in
end subroutine

subroutine amuse_set_gamma(gamma_in)
    use eos, only:gamma
    implicit none
    double precision, intent(in) :: gamma_in
    gamma = gamma_in
end subroutine

subroutine amuse_get_c_courant(C_cour_out)
    use timestep, only:C_cour
    implicit none
    double precision, intent(out) :: C_cour_out
    C_cour_out = C_cour
end subroutine

subroutine amuse_get_c_force(C_force_out)
    use timestep, only:C_force
    implicit none
    double precision, intent(out) :: C_force_out
    C_force_out = C_force
end subroutine

subroutine amuse_get_tolv(tolv_out)
    use timestep, only:tolv
    implicit none
    double precision, intent(out) :: tolv_out
    tolv_out = tolv
end subroutine

subroutine amuse_get_hfact(hfact_out)
    use part, only:hfact
    implicit none
    double precision, intent(out) :: hfact_out
    hfact_out = hfact
end subroutine

subroutine amuse_get_tolh(tolh_out)
    use options, only:tolh
    implicit none
    double precision, intent(out) :: tolh_out
    tolh_out = tolh
end subroutine

subroutine amuse_get_tree_accuracy(tree_accuracy_out)
    use kdtree, only:tree_accuracy
    implicit none
    double precision, intent(out) :: tree_accuracy_out
    tree_accuracy_out = tree_accuracy
end subroutine

subroutine amuse_get_alpha(alpha_out)
    use options, only:alpha
    implicit none
    double precision, intent(out) :: alpha_out
    alpha_out = alpha
end subroutine

subroutine amuse_get_alphamax(alphamax_out)
    use options, only:alphamax
    implicit none
    double precision, intent(out) :: alphamax_out
    alphamax_out = alphamax
end subroutine

subroutine amuse_get_beta(beta_out)
    use options, only:beta
    implicit none
    double precision, intent(out) :: beta_out
    beta_out = beta
end subroutine

subroutine amuse_get_avdecayconst(avdecayconst_out)
    use options, only:avdecayconst
    implicit none
    double precision, intent(out) :: avdecayconst_out
    avdecayconst_out = avdecayconst
end subroutine

subroutine amuse_get_idamp(idamp_out)
    use options, only:idamp
    implicit none
    integer, intent(out) :: idamp_out
    idamp_out = idamp
end subroutine

subroutine amuse_get_ieos(ieos_out)
    use eos, only:ieos
    implicit none
    integer, intent(out) :: ieos_out
    ieos_out = ieos
end subroutine

subroutine amuse_get_mu(mu_out)
    use eos, only:gmw
    implicit none
    double precision, intent(out) :: mu_out
    mu_out = gmw
end subroutine

subroutine amuse_get_rhofinal(rhofinal_out)
    use options, only:rhofinal_cgs
    use units, only:unit_density
    implicit none
    double precision, intent(out) :: rhofinal_out
    rhofinal_out = rhofinal_cgs / unit_density
end subroutine

subroutine amuse_get_rho_crit(rho_crit_out)
    use ptmass, only:rho_crit_cgs
    use units, only:unit_density
    implicit none
    double precision, intent(out) :: rho_crit_out
    rho_crit_out = rho_crit_cgs / unit_density
end subroutine

subroutine amuse_get_r_crit(r_crit_out)
    use ptmass, only:r_crit
    implicit none
    double precision, intent(out) :: r_crit_out
    r_crit_out = r_crit
end subroutine

subroutine amuse_get_h_acc(h_acc_out)
    use ptmass, only:h_acc
    implicit none
    double precision, intent(out) :: h_acc_out
    h_acc_out = h_acc
end subroutine

subroutine amuse_get_h_soft_sinkgas(h_soft_sinkgas_out)
    use ptmass, only:h_soft_sinkgas
    implicit none
    double precision, intent(out) :: h_soft_sinkgas_out
    h_soft_sinkgas_out = h_soft_sinkgas
end subroutine

subroutine amuse_get_h_soft_sinksink(h_soft_sinksink_out)
    use ptmass, only:h_soft_sinksink
    implicit none
    double precision, intent(out) :: h_soft_sinksink_out
    h_soft_sinksink_out = h_soft_sinksink
end subroutine

subroutine amuse_get_f_acc(f_acc_out)
    use ptmass, only:f_acc
    implicit none
    double precision, intent(out) :: f_acc_out
    f_acc_out = f_acc
end subroutine

subroutine amuse_get_iexternalforce(iexternalforce_out)
    use options, only:iexternalforce
    implicit none
    integer, intent(out) :: iexternalforce_out
    iexternalforce_out = iexternalforce
end subroutine

subroutine amuse_get_irealvisc(irealvisc_out)
    use viscosity, only:irealvisc
    implicit none
    integer, intent(out) :: irealvisc_out
    irealvisc_out = irealvisc
end subroutine

subroutine amuse_get_shearparam(shearparam_out)
    use viscosity, only:shearparam
    implicit none
    double precision, intent(out) :: shearparam_out
    shearparam_out = shearparam
end subroutine

subroutine amuse_get_bulkvisc(bulkvisc_out)
    use viscosity, only:bulkvisc
    implicit none
    double precision, intent(out) :: bulkvisc_out
    bulkvisc_out = bulkvisc
end subroutine

subroutine amuse_get_gamma(gamma_out)
    use eos, only:gamma
    implicit none
    double precision, intent(out) :: gamma_out
    gamma_out = gamma
end subroutine
