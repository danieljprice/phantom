!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module get_apr_level
!
! Module that holds the routines to return get_apr. This is where you set
! the shape of your APR region.
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: apr_region, dim, io, utils_apr
!
 use dim, only:use_apr,gr
 use apr_region
 use utils_apr

 implicit none

 public :: get_apr
 public :: create_or_update_apr_clump

 ! Declare the procedure pointer without initial assignment
 procedure(get_apr_sphere), pointer :: get_apr
 procedure(split_dir_one), pointer :: split_dir_func

contains

!-----------------------------------------------------------------------
!+
!  routine to set the get_apr and split_dir functions correctly
!+
!-----------------------------------------------------------------------
subroutine set_get_apr()

 ! Initialize the procedure pointer
 if (.not.associated(get_apr)) get_apr => get_apr_sphere
 if (.not.associated(split_dir_func)) split_dir_func => split_dir_one

 ! For the apr type, chose the region shape
 if (apr_type == 6) then
    ref_dir = -1 ! need to enforce this for this one
 else
    get_apr => get_apr_sphere
 endif

 ! Here set the requirements for the apr_type to read in the right values
! for apr_types that read in a particle number
 if (apr_type == 2) then
    track_part(1) = track_part_in
 endif

 ! for apr_types that read in the centre values from the *.in file
 if (apr_type == 1) apr_centre(1:3,1) = apr_centre_in(1:3) ! from the .in file

 ! set the split direction function
 if (gr) then
   split_dir_func => split_dir_gr
 elseif (split_dir == 1) then
   split_dir_func => split_dir_one
 elseif (split_dir == 2) then
   split_dir_func => split_dir_two
 else
   split_dir_func => split_dir_three
 endif

end subroutine set_get_apr

!-----------------------------------------------------------------------
!+
!  routine to return the adaptive particle refinement level based on position
!  and the boundaries set by the apr_* arrays for a spherical region
!+
!-----------------------------------------------------------------------
pure subroutine get_apr_sphere(pos,icentre,apri)
 use io, only:fatal
 use apr_region, only:apr_region_is_circle
 real, intent(in)     :: pos(3)
 integer, intent(in)  :: icentre
 integer, intent(out) :: apri
 integer :: jj, kk
 real :: dx,dy,dz,r

 apri = -1 ! to prevent compiler warnings

 do jj = 1,apr_max
    if (ref_dir == 1) then
       kk = apr_max - jj + 1       ! going from apr_max -> 1
    else
       kk = jj                    ! going from 1 -> apr_max
    endif
    dx = pos(1) - apr_centre(1,icentre)
    dy = pos(2) - apr_centre(2,icentre)
    dz = pos(3) - apr_centre(3,icentre)

    if (apr_region_is_circle) dz = 0.

    r = sqrt(dx**2 + dy**2 + dz**2)
    if (r < apr_regions(kk)) then
       apri = kk
       return
    endif
 enddo

end subroutine get_apr_sphere

!-----------------------------------------------------------------------
!+
!  small routine to set up the new particles happily
!+
!-----------------------------------------------------------------------
subroutine set_new_splitpart(i,i_new,v,sep)
 use part, only:xyzh, vxyzu, apr_level, copy_particle_all
 use dim,  only:ind_timesteps
 real, intent(in) :: v(3), sep
 integer, intent(in) :: i,i_new

 real :: x_add, y_add, z_add
 integer(kind=1) :: aprnew

 x_add = sep*v(1)*xyzh(4,i)
 y_add = sep*v(2)*xyzh(4,i)
 z_add = sep*v(3)*xyzh(4,i)

 aprnew = apr_level(i) + int(1,kind=1) ! to prevent compiler warnings

 !--create the new particle
 call copy_particle_all(i,i_new,new_part=.true.)
 xyzh(1,i_new) = xyzh(1,i) + x_add
 xyzh(2,i_new) = xyzh(2,i) + y_add
 xyzh(3,i_new) = xyzh(3,i) + z_add
 vxyzu(:,i_new) = vxyzu(:,i)
 xyzh(4,i_new) = xyzh(4,i)*(0.5**(1./3.))
 apr_level(i_new) = aprnew
 if (ind_timesteps) call put_in_smallest_binn(i_new)

 ! Edit the old particle that was sent in and kept
 xyzh(1,i) = xyzh(1,i) - x_add
 xyzh(2,i) = xyzh(2,i) - y_add
 xyzh(3,i) = xyzh(3,i) - z_add
 apr_level(i) = aprnew
 xyzh(4,i) = xyzh(4,i)*(0.5**(1./3.))
 if (ind_timesteps) call put_in_smallest_binn(i)


end subroutine set_new_splitpart

!-----------------------------------------------------------------------
!+
!  core splitpart routine for split_dir = 1
!+
!-----------------------------------------------------------------------
subroutine split_dir_one(i,i_new,sep)
 use part, only:xyzh
 use apr_region, only:apr_region_is_circle
 use vectorutils, only:cross_product3D,rotatevec
 use physcon, only:pi
 integer, intent(in) :: i,i_new
 real, intent(inout) :: sep
 real :: dx, dy, dz, u(3), w(3), theta, mag_v, v(3), angle
 integer, save :: nangle = 1

 angle = nangle*(1./sqrt(2.)) - nint(nangle*(1./sqrt(2.)))
 !nangle = nangle + 1 ! for next round

 ! Calculate the plane that the particle must be split along
 ! to be tangential to the splitting region. Particles are split
 ! on this plane but rotated randomly on it.

 dx = xyzh(1,i) - apr_centre(1,icentre)
 dy = xyzh(2,i) - apr_centre(2,icentre)
 if (.not.apr_region_is_circle) then
    dz = xyzh(3,i) - apr_centre(3,icentre)       ! for now, split about centre

    ! Calculate a vector, v, that lies on the plane
    u = (/1.0,0.5,1.0/)
    w = (/dx,dy,dz/)
    call cross_product3D(u,w,v)

    ! rotate it around the normal to the plane by a random amount
    theta = angle*2.*pi
    call rotatevec(v,w,theta)

    mag_v = sqrt(dot_product(v,v))
    if (mag_v > tiny(mag_v)) then
       v = v/mag_v
    else
       v = 0.
    endif
 else
    dz = 0.
    u = 0.
    w = 0.
    v = 0.
    theta = atan2(dy,dx) + 0.5*pi
    v(1) = cos(theta)
    v(2) = sin(theta)
 endif

 ! Now apply it
 call set_new_splitpart(i,i_new,v,sep)

end subroutine split_dir_one

!-----------------------------------------------------------------------
!+
!  core splitpart routine for split_dir = 3
!+
!-----------------------------------------------------------------------
subroutine split_dir_three(i,i_new,sep)
 use part, only:xyzh
 use apr_region, only:apr_region_is_circle
 use physcon, only:pi
 integer, intent(in) :: i,i_new
 real, intent(inout) :: sep
 real :: dx, dy, dz, u(3), w(3), theta, mag_v, a, b, c, v(3), angle(3)
 integer, save :: nangle = 1

 ! set the "random" vector directions using some irrational numbers
 ! this just increments a different irrational number for each
 ! and ensures it's between -0.5 -> 0.5
 angle(1) = nangle*(1./sqrt(2.)) - nint(nangle*(1./sqrt(2.)))
 angle(2) = nangle*(sqrt(2.) - 1.) - nint(nangle*(sqrt(2.) - 1.))
 angle(3) = nangle*(pi - 3.) - nint(nangle*(pi - 3.))
!nangle = nangle + 1 ! for next round

 if (.not.apr_region_is_circle) then
    ! No directional splitting, so just create a unit vector in a random direction
    a = angle(1)
    b = angle(2)
    c = angle(3)
    v = (/a, b, c/)

    mag_v = sqrt(dot_product(v,v))
    if (mag_v > tiny(mag_v)) then
       v = v/mag_v
    else
       v = 0.
    endif
 else
    dx = xyzh(1,i) - apr_centre(1,icentre)
    dy = xyzh(2,i) - apr_centre(2,icentre)
    dz = 0.
    u = 0.
    w = 0.
    v = 0.
    theta = atan2(dy,dx) + 0.5*pi
    v(1) = cos(theta)
    v(2) = sin(theta)
 endif

 ! Now apply it
 call set_new_splitpart(i,i_new,v,sep)

end subroutine split_dir_three

!-----------------------------------------------------------------------
!+
!  core splitpart routine for split_dir = 3
!+
!-----------------------------------------------------------------------
subroutine split_dir_two(i,i_new,sepin)
 use part, only:xyzh,apr_level,copy_particle_all,aprmassoftype,igas, &
                vxyzu
 use vectorutils, only:cross_product3D,rotatevec
 use physcon, only:pi
 use dim,  only:ind_timesteps
 integer, intent(in) :: i,i_new
 real, intent(inout) :: sepin
 real :: pmass, uold, hnew, sep

 sep = sepin*xyzh(4,i)

 apr_level(i) = apr_level(i) + int(1,kind=1) ! to prevent compiler warnings
 call copy_particle_all(i,i_new,new_part=.true.)
 pmass = aprmassoftype(igas,apr_level(i))

 uold = vxyzu(4,i)
 hnew = xyzh(4,i)*(0.5**(1./3.))

 ! new part forward
 xyzh(4,i_new) = hnew ! set new smoothing length
 call integrate_geodesic(pmass,xyzh(:,i_new),vxyzu(:,i_new),sep,1.)
 if (ind_timesteps) call put_in_smallest_binn(i_new)

 ! old part backward
 ! switch direction
 vxyzu(1:3,i) = -vxyzu(1:3,i)
 xyzh(4,i) = hnew
 call integrate_geodesic(pmass,xyzh(:,i),vxyzu(:,i),sep,1.)
 ! switch direction back
 vxyzu(1:3,i) = -vxyzu(1:3,i)
 vxyzu(4,i) = uold
 if (ind_timesteps) call put_in_smallest_binn(i)

end subroutine split_dir_two

!-----------------------------------------------------------------------
!+
!  core splitpart routine for split_dir = gr only
!+
!-----------------------------------------------------------------------
subroutine split_dir_gr(i,i_new,sepin)
 use part, only:xyzh,apr_level,copy_particle_all,aprmassoftype,igas,itemp, &
                vxyzu,pxyzu,dens,eos_vars,igasP,igamma,metrics,metricderivs,fext
 use vectorutils, only:cross_product3D,rotatevec
 use physcon, only:pi
 use dim,  only:ind_timesteps
 use metric_tools, only:pack_metric,pack_metricderivs
 use extern_gr, only:get_grforce
 integer, intent(in) :: i,i_new
 real, intent(inout) :: sepin
 real :: pmass, uold, hnew, sep

 sep = sepin*xyzh(4,i)

 apr_level(i) = apr_level(i) + int(1,kind=1) ! to prevent compiler warnings
 call copy_particle_all(i,i_new,new_part=.true.)
 pmass = aprmassoftype(igas,apr_level(i))

 uold = vxyzu(4,i)
 hnew = xyzh(4,i)*(0.5**(1./3.))

 ! new part forward
 xyzh(4,i_new) = hnew ! set new smoothing length
 call integrate_geodesic_gr(pmass,xyzh(:,i_new),vxyzu(:,i_new),dens(i_new),eos_vars(igasP,i_new), &
                         eos_vars(igamma,i_new),eos_vars(itemp,i_new),pxyzu(:,i_new),sep)
 call pack_metric(xyzh(1:3,i_new),metrics(:,:,:,i_new))
 call pack_metricderivs(xyzh(1:3,i_new),metricderivs(:,:,:,i_new))
 call get_grforce(xyzh(:,i_new),metrics(:,:,:,i_new),metricderivs(:,:,:,i_new), &
                  vxyzu(1:3,i_new),dens(i_new),vxyzu(4,i_new),eos_vars(igasP,i_new),fext(1:3,i_new))
 if (ind_timesteps) call put_in_smallest_binn(i_new)

 ! old part backward
 ! switch direction
 vxyzu(1:3,i) = -vxyzu(1:3,i)
 pxyzu(1:3,i) = -pxyzu(1:3,i)
 xyzh(4,i) = hnew
 call integrate_geodesic_gr(pmass,xyzh(:,i),vxyzu(:,i),dens(i),eos_vars(igasP,i),eos_vars(igamma,i),eos_vars(itemp,i), &
                        pxyzu(:,i),sep)
 ! switch direction back
 vxyzu(1:3,i) = -vxyzu(1:3,i)
 pxyzu(1:3,i) = -pxyzu(1:3,i)
 call pack_metric(xyzh(1:3,i),metrics(:,:,:,i))
 call pack_metricderivs(xyzh(1:3,i),metricderivs(:,:,:,i))
 call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i), &
                  vxyzu(1:3,i),dens(i),vxyzu(4,i),eos_vars(igasP,i),fext(1:3,i))
 if (ind_timesteps) call put_in_smallest_binn(i)

end subroutine split_dir_gr

!-----------------------------------------------------------------------
!+
!  routine to put a particle on the shortest timestep
!+
!-----------------------------------------------------------------------
subroutine put_in_smallest_binn(i)
 use timestep_ind, only:nbinmax
 use part,         only:ibin
 integer, intent(in) :: i

 ibin(i) = nbinmax

end subroutine put_in_smallest_binn

!-----------------------------------------------------------------------
!+
!  Integrate particle along the geodesic
!  Update vel and metric for best energy conservation
!+
!-----------------------------------------------------------------------
subroutine integrate_geodesic(pmass,xyzh,vxyzu,dist,timei)
 use options,        only:iexternalforce
 use externalforces, only:externalforce,externalforce_vdependent
 use part,           only:rhoh
 real, intent(inout) :: xyzh(:),vxyzu(:)
 real, intent(in)    :: dist,timei,pmass
 real :: fext(3),fextv(3)
 real :: t,tend,v,dt,dens
 real :: xyz(3),vxyz(1:3),poti,uui

 xyz       = xyzh(1:3)
 vxyz      = vxyzu(1:3)
 fext      = 0.
 uui       = vxyzu(4)
 dens      = rhoh(xyzh(4),pmass)

 v = sqrt(dot_product(vxyz,vxyz))
 tend = dist/v

 if (iexternalforce > 0) then
    call externalforce(iexternalforce,xyz(1),xyz(2),xyz(3),xyzh(4), &
                                timei,fext(1),fext(2),fext(3),poti,dt)
    call externalforce_vdependent(iexternalforce,xyz,vxyz,fextv,poti,dens,uui)
    fext = fext + fextv
 endif

 t = 0.

 do while (t <= tend)
    dt = min(0.1,tend*0.1,dt,0.1*dt)
    t    = t + dt

    vxyz = vxyz + dt*fext
    xyz = xyz + dt*vxyz

    if (iexternalforce > 0) then
       call externalforce(iexternalforce,xyz(1),xyz(2),xyz(3),xyzh(4), &
                                timei,fext(1),fext(2),fext(3),poti,dt)
       call externalforce_vdependent(iexternalforce,xyz,vxyz,fextv,poti,dens,uui)
       fext = fext + fextv
    endif
 enddo

 xyzh(1:3) = xyz(1:3)
 vxyzu(1:3) = vxyz(1:3)
end subroutine integrate_geodesic

!-----------------------------------------------------------------------
!+
!  Integrate particle along the geodesic
!  Update vel and metric for best energy conservation
!+
!-----------------------------------------------------------------------
subroutine integrate_geodesic_gr(pmass,xyzh,vxyzu,dens,pr,gamma,temp,pxyzu,dist)
 use extern_gr,      only:get_grforce
 use metric_tools,   only:pack_metric,pack_metricderivs
 use eos,            only:ieos,equationofstate
 use cons2primsolver,only:conservative2primitive
 use io,             only:warning
 use part,           only:rhoh,ien_type
 real, intent(inout) :: xyzh(:),vxyzu(:),pxyzu(:)
 real, intent(inout) :: dens,pr,gamma,temp,pmass
 real, intent(in)    :: dist
 real :: metrics(0:3,0:3,2),metricderivs(0:3,0:3,3),fext(3)
 real :: t,tend,v,dt
 real :: xyz(3),pxyz(3),eni,vxyz(1:3),uui,rho,spsoundi,pondensi
 integer :: ierr

 xyz       = xyzh(1:3)
 pxyz      = pxyzu(1:3)
 eni       = pxyzu(4)
 vxyz      = vxyzu(1:3)
 uui       = vxyzu(4)
 rho       = rhoh(xyzh(4),pmass)

 v = sqrt(dot_product(vxyz,vxyz))
 tend = dist/v

 call pack_metric(xyz,metrics(:,:,:))
 call pack_metricderivs(xyz,metricderivs(:,:,:))
 call get_grforce(xyzh(:),metrics(:,:,:),metricderivs(:,:,:),vxyz,dens,uui,pr,fext(1:3),dt)

 t = 0.

 do while (t <= tend)
    dt = min(0.1,tend*0.1,dt)
    t    = t + dt
    pxyz = pxyz + dt*fext

    call conservative2primitive(xyz,metrics(:,:,:),vxyz,dens,uui,pr,&
                                       temp,gamma,rho,pxyz,eni,ierr,ien_type)
    if (ierr > 0) call warning('cons2primsolver [in integrate_geodesic (a)]','did not converge')

    xyz = xyz + dt*vxyz
    call pack_metric(xyz,metrics(:,:,:))
    call pack_metricderivs(xyz,metricderivs(:,:,:))

    call equationofstate(ieos,pondensi,spsoundi,dens,xyzh(1),xyzh(2),xyzh(3),temp,uui)
    pr = pondensi*dens

    call get_grforce(xyzh(:),metrics(:,:,:),metricderivs(:,:,:),vxyzu(1:3),dens,vxyzu(4),pr,fext(1:3),dt)
 enddo

 xyzh(1:3) = xyz(1:3)
 vxyzu(1:3) = vxyz(1:3)
 pxyzu(1:3) = pxyz(1:3)

end subroutine integrate_geodesic_gr


end module get_apr_level
