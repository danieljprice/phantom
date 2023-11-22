!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module raytracer
!
! This module contains all routines required to:
!   - perform radial ray tracing starting from the primary star
!   - calculate optical depts along the rays given the opacity distribution
!   - interpolate optical depths to all SPH particles
! Applicable both for single star as well as binary models
!
! WARNING: This module has only been tested on phantom wind setups
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: healpix, kernel, linklist, part, units
!
 use healpix

 implicit none
 public :: get_all_tau

 private

contains

 !--------------------------------------------------------------------------
 !+
 !  MAIN ROUTINE
 !  Returns the optical depth to each SPH particle, using the uniform outwards
 !  ray-tracing scheme.
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: nptmass:         The number of sink particles
 !  IN: xyzm_ptmass:     The array containing the properties of the sink particle
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa_cgs:       The array containing the opacities of all SPH particles
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !+
 !  OUT: tau:            The array of optical depths for each SPH particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, kappa_cgs, order, tau)
 use part, only: iReff
 integer, intent(in) :: npart, order, nptmass
 real, intent(in)    :: kappa_cgs(:), xyzh(:,:), xyzmh_ptmass(:,:)
 real, intent(out)   :: tau(:)

 if (nptmass == 2 ) then
    call get_all_tau_companion(npart, xyzmh_ptmass(1:3,1), xyzh, kappa_cgs, &
              xyzmh_ptmass(iReff,1), xyzmh_ptmass(1:3,2), xyzmh_ptmass(iReff,2), order, tau)
 else
    call get_all_tau_single(npart, xyzmh_ptmass(1:3,1), xyzh, kappa_cgs, xyzmh_ptmass(iReff,1), order, tau)
 endif
end subroutine get_all_tau

 !--------------------------------------------------------------------------
 !+
 !  Calculates the optical depth to each SPH particle, using the uniform outwards
 !  ray-tracing scheme for models containing a single star
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the kappa of all SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !+
 !  OUT: taus:           The array of optical depths to each SPH particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_single(npart, primary, xyzh, kappa, Rstar, order, tau)
 use part, only : isdead_or_accreted
 integer, intent(in) :: npart,order
 real, intent(in)    :: primary(3), kappa(:), Rstar, xyzh(:,:)
 real, intent(out)   :: tau(:)

 integer  :: i, nrays, nsides
 real     :: ray_dir(3),part_dir(3)
 real, dimension(:,:), allocatable  :: rays_dist, rays_tau
 integer, dimension(:), allocatable :: rays_dim
 integer, parameter :: ndim = 200

 nrays = 12*4**order ! The number of rays traced given the healpix order
 nsides = 2**order   ! The healpix nsides given the healpix order

 allocate(rays_dist(ndim, nrays))
 allocate(rays_tau(ndim, nrays))
 allocate(rays_dim(nrays))

 !-------------------------------------------
 ! CONSTRUCT the RAYS given the ORDER
 ! and determine the optical depth along them
 !-------------------------------------------

!$omp parallel default(none) &
!$omp private(ray_dir) &
!$omp shared(nrays,nsides,primary,kappa,xyzh,Rstar,rays_dist,rays_tau,rays_dim)
!$omp do
 do i = 1, nrays
    !returns ray_dir, the unit vector identifying ray (index i-1 becase healpix starts counting at index 0)
    call pix2vec_nest(nsides, i-1, ray_dir)
    !calculate the properties along the ray
    call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,rays_tau(:,i),rays_dist(:,i),rays_dim(i))
 enddo
!$omp enddo
!$omp end parallel


 !_----------------------------------------------
 ! DETERMINE the optical depth for each particle
 ! using the values available on the rays
 !-----------------------------------------------

!$omp parallel default(none) &
!$omp private(part_dir) &
!$omp shared(npart,primary,nsides,xyzh,ray_dir,rays_dist,rays_tau,rays_dim,tau)
!$omp do
 do i = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       part_dir = xyzh(1:3,i)-primary
       call interpolate_tau(nsides, part_dir, rays_tau, rays_dist, rays_dim, tau(i))
    else
       tau(i) = -99.
    endif
 enddo
!$omp enddo
!$omp end parallel

end subroutine get_all_tau_single

 !--------------------------------------------------------------------------
 !+
 !  Calculates the optical depth to each SPH particle, using the uniform outwards
 !  ray-tracing scheme for models containing primary star and a companion
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the opacity of all the SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: companion:       The xyz coordinates of the companion
 !  IN: Rcomp:           The radius of the companion
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !+
 !  OUT: tau:            The array of optical depths for each SPH particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_companion(npart, primary, xyzh, kappa, Rstar, companion, Rcomp, order, tau)
 use part, only : isdead_or_accreted
 integer, intent(in) :: npart, order
 real, intent(in)    :: primary(3), companion(3), kappa(:), Rstar, xyzh(:,:), Rcomp
 real, intent(out)   :: tau(:)

 integer  :: i, nrays, nsides
 real     :: normCompanion,theta0,phi,cosphi,sinphi,theta,sep,root
 real     :: ray_dir(3),part_dir(3),uvecCompanion(3)
 real, dimension(:,:), allocatable  :: dirs
 real, dimension(:,:), allocatable  :: rays_dist, rays_tau
 integer, dimension(:), allocatable :: rays_dim
 integer, parameter :: ndim = 200

 nrays = 12*4**order ! The number of rays traced given the healpix order
 nsides = 2**order   ! The healpix nsides given the healpix order

 allocate(dirs(3, nrays))
 allocate(rays_dist(ndim, nrays))
 allocate(rays_tau(ndim, nrays))
 allocate(rays_dim(nrays))

 uvecCompanion = companion-primary
 normCompanion = norm2(uvecCompanion)
 uvecCompanion = uvecCompanion/normCompanion
 theta0        = asin(Rcomp/normCompanion)
 phi           = atan2(uvecCompanion(2),uvecCompanion(1))
 cosphi        = cos(phi)
 sinphi        = sin(phi)

 !-------------------------------------------
 ! CONSTRUCT the RAYS given the ORDER
 ! and determine the optical depth along them
 !-------------------------------------------

!$omp parallel default(none) &
!$omp private(ray_dir,theta,root,sep) &
!$omp shared(nrays,nsides,primary,kappa,xyzh,Rstar,Rcomp,rays_dist,rays_tau,rays_dim) &
!$omp shared(uvecCompanion,normCompanion,cosphi,sinphi,theta0)
!$omp do
 do i = 1, nrays
    !returns ray_dir, the unit vector identifying ray (index i-1 becase healpix starts counting at index 0)
    call pix2vec_nest(nsides, i-1, ray_dir)
    !rotate ray vectors by an angle = phi so the main axis points to the companion (This is because along the
    !main axis (1,0,0) rays are distributed more uniformally
    ray_dir = (/cosphi*ray_dir(1) - sinphi*ray_dir(2),sinphi*ray_dir(1) + cosphi*ray_dir(2), ray_dir(3)/)
    theta   = acos(dot_product(uvecCompanion, ray_dir))
    !the ray intersects the companion: only calculate tau up to the companion
    if (theta < theta0) then
       root  = sqrt(Rcomp**2-normCompanion**2*sin(theta)**2)
       sep   = normCompanion*cos(theta)-root
       call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,rays_tau(:,i),rays_dist(:,i),rays_dim(i), sep)
    else
       call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,rays_tau(:,i),rays_dist(:,i),rays_dim(i))
    endif
 enddo
!$omp enddo
!$omp end parallel

 !-----------------------------------------------
 ! DETERMINE the optical depth for each particle
 ! using the values available on the rays
 !-----------------------------------------------

!$omp parallel default(none) &
!$omp private(part_dir) &
!$omp shared(npart,primary,cosphi,sinphi,nsides,xyzh,ray_dir,rays_dist,rays_tau,rays_dim,tau)
!$omp do
 do i = 1, npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       !vector joining the source to the particle
       part_dir = xyzh(1:3,i)-primary
       part_dir = (/cosphi*part_dir(1) + sinphi*part_dir(2),-sinphi*part_dir(1) + cosphi*part_dir(2), part_dir(3)/)
       call interpolate_tau(nsides, part_dir, rays_tau, rays_dist, rays_dim, tau(i))
    else
       tau(i) = -99.
    endif
 enddo
!$omp enddo
!$omp end parallel
end subroutine get_all_tau_companion

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth of a particle.
 !  Search for the four closest rays to a particle, perform four-point
 !  interpolation of the optical depts from these rays. Weighted by the
 !  inverse square of the perpendicular distance to the rays.
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: nsides:          The healpix nsides of the simulation
 !  IN: vec:             The vector from the primary to a point
 !  IN: rays_tau:        2-dimensional array containing the cumulative optical
 !                       depts along each ray
 !  IN: rays_dist:       2-dimensional array containing the distances from the
 !                       primary along each ray
 !  IN: rays_dim:        The vector containing the number of points defined along each ray
 !+
 !  OUT: tau:            The interpolated optical depth at the particle's location
 !+
 !--------------------------------------------------------------------------
subroutine interpolate_tau(nsides, vec, rays_tau, rays_dist, rays_dim, tau)
 integer, intent(in) :: nsides, rays_dim(:)
 real, intent(in)    :: vec(:), rays_dist(:,:), rays_tau(:,:)
 real, intent(out)   :: tau

 integer :: rayIndex, neighbours(8), nneigh, i, k
 real    :: tautemp, ray(3), vectemp(3), weight, tempdist(8), distRay_sq, vec_norm2
 logical :: mask(8)

 vec_norm2 = norm2(vec)
 !returns rayIndex, the index of the ray vector that points to the particle (direction vec)
 call vec2pix_nest(nsides, vec, rayIndex)
 !returns ray(3), the unit vector identifying the ray with index number rayIndex
 call pix2vec_nest(nsides, rayIndex, ray)
 vectemp       = vec - vec_norm2*ray
 distRay_sq    = dot_product(vectemp,vectemp)
 call get_tau_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), tautemp)
 if (distRay_sq > 0.) then
    tau    = tautemp/distRay_sq
    weight = 1./distRay_sq
 else
    ! the particle sits exactly on the ray, no need to get the neighbours
    tau    = tautemp
    return
 endif

 !returns the number nneigh and list of vectors (n) neighbouring the ray number index
 call neighbours_nest(nsides, rayIndex, neighbours, nneigh)
 !for each neighbouring ray calculate its distance to the particle
 do i=1,nneigh
    call pix2vec_nest(nsides, neighbours(i), ray)
    vectemp     = vec - vec_norm2*ray
    tempdist(i) = dot_product(vectemp,vectemp)
 enddo
 neighbours       = neighbours+1
 mask             = .true.
 if (nneigh <8) mask(nneigh+1:8) = .false.
 !take tau contribution from the 3 closest rays
 do i=1,3
    k       = minloc(tempdist,1,mask)
    mask(k) = .false.
    call get_tau_on_ray(vec_norm2, rays_tau(:,neighbours(k)), &
                rays_dist(:,neighbours(k)), rays_dim(neighbours(k)), tautemp)
    tau    = tau + tautemp/tempdist(k)
    weight = weight + 1./tempdist(k)
 enddo
 tau = tau / weight
end subroutine interpolate_tau


 !--------------------------------------------------------------------------
 !+
 !  Interpolation of the optical depth for an arbitrary point on the ray,
 !  with a given distance to the starting point of the ray.
 !+
 !  IN: distance:        The distance from the staring point of the ray to a
 !                       point on the ray
 !  IN: tau_along_ray:   The vector of cumulative optical depths along the ray
 !  IN: dist_along_ray:  The vector of distances from the primary along the ray
 !  IN: len:             The length of listOfTau and listOfDist
 !+
 !  OUT: tau:            The optical depth to the given distance along the ray
 !+
 !--------------------------------------------------------------------------
subroutine get_tau_on_ray(distance, tau_along_ray, dist_along_ray, len, tau)
 real, intent(in)    :: distance, tau_along_ray(:), dist_along_ray(:)
 integer, intent(in) :: len
 real, intent(out)   :: tau

 integer :: L, R, m ! left, right and middle index for binary search

 if (distance  <  dist_along_ray(1)) then
    tau = 0.
 elseif (distance  >  dist_along_ray(len)) then
    tau = 99.
 else
    L = 2
    R = len-1
    !bysection search for the index of the closest ray location to the particle
    do while (L < R)
       m = (L + R)/2
       if (dist_along_ray(m) > distance) then
          R = m
       else
          L = m + 1
       endif
    enddo
    !interpolate linearly ray properties to get the particle's optical depth
    tau = tau_along_ray(L-1)+(tau_along_ray(L)-tau_along_ray(L-1))/ &
                  (dist_along_ray(L)-dist_along_ray(L-1))*(distance-dist_along_ray(L-1))
 endif
end subroutine get_tau_on_ray

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth along a given ray
 !+
 !  IN: primary:         The location of the primary star
 !  IN: ray:             The unit vector of the direction in which the
 !                       optical depts will be calculated
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the particles opacity
 !  IN: Rstar:           The radius of the primary star
 !+
 !  OUT: tau_along_ray:  The vector of cumulative optical depts along the ray
 !  OUT: dist_along_ray: The vector of distances from the primary along the ray
 !  OUT: len:            The length of tau_along_ray and dist_along_ray
 !+
 !  OPT: maxDistance:    The maximal distance the ray needs to be traced
 !+
 !--------------------------------------------------------------------------
subroutine ray_tracer(primary, ray, xyzh, kappa, Rstar, tau_along_ray, dist_along_ray, len, maxDistance)
 use units, only:umass,udist
 real, intent(in)     :: primary(3), ray(3), Rstar, xyzh(:,:), kappa(:)
 real, optional       :: maxDistance
 real, intent(out)    :: dist_along_ray(:), tau_along_ray(:)
 integer, intent(out) :: len

 real    :: dr, next_dr, h, dtaudr, previousdtaudr, nextdtaudr, distance
 integer :: inext, i

 h = Rstar/100.
 inext=0
 do while (inext==0)
    h = h*2.
    !find the next point along the ray : index next
    call find_next(primary+Rstar*ray, h, ray, xyzh, kappa, previousdtaudr, dr, inext)
 enddo

 i = 1
 tau_along_ray(i)  = 0.
 distance          = Rstar
 dist_along_ray(i) = distance
 do while (hasNext(inext,tau_along_ray(i),distance,maxDistance))
    distance = distance+dr
    call find_next(primary + distance*ray, xyzh(4,inext), ray, xyzh, kappa, nextdtaudr, next_dr, inext)
    i = i + 1
    dtaudr            = (nextdtaudr+previousdtaudr)/2.
    previousdtaudr    = nextdtaudr
    !fix units for tau (kappa is in cgs while rho & r are in code units)
    tau_along_ray(i)  = tau_along_ray(i-1)+dr*dtaudr*umass/(udist**2)
    dist_along_ray(i) = distance
    dr                = next_dr
 enddo
 len = i
end subroutine ray_tracer

logical function hasNext(inext, tau, distance, maxDistance)
 integer, intent(in) :: inext
 real, intent(in)    :: distance, tau
 real, optional      :: maxDistance
 real, parameter :: tau_max = 99.
 if (present(maxDistance)) then
    hasNext = inext /= 0 .and. distance < maxDistance .and. tau < tau_max
 else
    hasNext = inext /= 0 .and. tau < tau_max
 endif
end function hasNext

 !--------------------------------------------------------------------------
 !+
 !  First finds the local optical depth derivative at the starting point, then finds the next
 !                       point on a ray and the distance to this point
 !+
 !  IN: inpoint:         The coordinate of the initial point projected on the
 !                       ray for which the opacity and the next point will be
 !                       calculated
 !  IN: h:               The smoothing length at the initial point
 !  IN: ray:             The unit vector of the direction in which the next
 !                       point will be calculated
 !  IN: xyzh:            The array containing the particles position+smoothing length
 !  IN: kappa:           The array containing the particles opacity
 !  IN: inext:           The index of the initial point
 !                       (this point will not be considered as possible next point)
 !+
 !  OUT: dtaudr:         The local optical depth derivative at the given location (inpoint)
 !  OUT: distance:       The distance to the next point
 !  OUT: inext:          The index of the next point on the ray
 !+
 !--------------------------------------------------------------------------
subroutine find_next(inpoint, h, ray, xyzh, kappa, dtaudr, distance, inext)
 use linklist, only:getneigh_pos,ifirstincell,listneigh
 use kernel,   only:radkern,cnormk,wkern
 use part,     only:hfact,rhoh,massoftype,igas
 real,    intent(in)    :: xyzh(:,:), kappa(:), inpoint(:), ray(:), h
 integer, intent(inout) :: inext
 real,    intent(out)   :: distance, dtaudr

 integer, parameter :: nmaxcache = 0
 real  :: xyzcache(0,nmaxcache)

 integer  :: nneigh, i, prev
 real     :: dmin, vec(3), dr, raydistance, q, norm_sq

 prev     = inext
 inext    = 0
 distance = 0.

 !for a given point (inpoint), returns the list of neighbouring particles (listneigh) within a radius h*radkern
 call getneigh_pos(inpoint,0.,h*radkern,3,listneigh,nneigh,xyzh,xyzcache,nmaxcache,ifirstincell)

 dtaudr = 0
 dmin = huge(0.)
 !loop over all neighbours
 do i=1,nneigh
    vec     = xyzh(1:3,listneigh(i)) - inpoint
    norm_sq = dot_product(vec,vec)
    q       = sqrt(norm_sq)/xyzh(4,listneigh(i))
    !add optical depth contribution from each particle
    dtaudr = dtaudr+wkern(q*q,q)*kappa(listneigh(i))*rhoh(xyzh(4,listneigh(i)), massoftype(igas))

    ! find the next particle : among the neighbours find the particle located the closest to the ray
    if (listneigh(i)  /=  prev) then
       dr = dot_product(vec,ray) !projected distance along the ray
       if (dr>0.) then
          !distance perpendicular to the ray direction
          raydistance = norm_sq - dr**2
          if (raydistance < dmin) then
             dmin     = raydistance
             inext    = listneigh(i)
             distance = dr
          endif
       endif
    endif
 enddo
 dtaudr = dtaudr*cnormk/hfact**3
end subroutine find_next
end module raytracer
