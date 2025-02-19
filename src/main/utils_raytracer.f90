!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module raytracer
!
! This module contains all routines required to:
!   - perform radial ray tracing starting from the primary star only
!   - calculate optical depth along the rays given the opacity distribution
!   - interpolate optical depths to all SPH particles
! Applicable both for single and binary star wind simulations
!
! WARNING: This module has only been tested on phantom wind setup
!
! :References: Esseldeurs M., Siess L. et al, 2023, A&A, 674, A122
!
! :Owner: Mats Esseldeurs
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

 !------------------------------------------------------------------------------------
 !+
 !  MAIN ROUTINE
 !  Returns the optical depth at each particle's location using an outward ray-tracing scheme
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
 !------------------------------------------------------------------------------------
subroutine get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, kappa_cgs, order, tau)
 use part,   only: iReff
 integer, intent(in) :: npart, order, nptmass
 real, intent(in)    :: kappa_cgs(:), xyzh(:,:), xyzmh_ptmass(:,:)
 real, intent(out)   :: tau(:)
 real :: Rinject

 Rinject = xyzmh_ptmass(iReff,1)
 if (nptmass == 2 ) then
    call get_all_tau_companion(npart, xyzmh_ptmass(1:3,1), xyzmh_ptmass(iReff,1), xyzh, kappa_cgs, &
            Rinject, xyzmh_ptmass(1:3,2), xyzmh_ptmass(iReff,2), order, tau)
 else
    call get_all_tau_single(npart, xyzmh_ptmass(1:3,1), xyzmh_ptmass(iReff,1), xyzh,&
         kappa_cgs, Rinject, order, tau)
 endif
end subroutine get_all_tau

 !---------------------------------------------------------------------------------
 !+
 !  Calculates the optical depth at each particle's location, using the uniform outward
 !  ray-tracing scheme for models containing a single star
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the kappa of all SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: Rinject:         The particles injection radius
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !+
 !  OUT: taus:           The array of optical depths to each SPH particle
 !+
 !---------------------------------------------------------------------------------
subroutine get_all_tau_single(npart, primary, Rstar, xyzh, kappa, Rinject, order, tau)
 use part, only : isdead_or_accreted
 integer, intent(in) :: npart,order
 real, intent(in)    :: primary(3), kappa(:), Rstar, Rinject, xyzh(:,:)
 real, intent(out)   :: tau(:)

 integer  :: i, nrays, nsides
 real     :: ray_dir(3),part_dir(3)
 real, dimension(:,:), allocatable  :: rays_dist, rays_tau
 integer, dimension(:), allocatable :: rays_dim
 integer, parameter :: ndim = 200 ! maximum number of points along the ray where tau is calculated

 nrays = 12*4**order ! The number of rays traced given the healpix order
 nsides = 2**order   ! The healpix nsides given the healpix order

 allocate(rays_dist(ndim, nrays)) ! distance from the central star of the points on the rays
 allocate(rays_tau(ndim, nrays))  ! value of tau at each point along each ray
 allocate(rays_dim(nrays))        ! effective number of points on the ray (< ndim)

 !-------------------------------------------
 ! CONSTRUCT the RAYS given the HEALPix ORDER
 ! and determine the optical depth along them
 !-------------------------------------------

 !$omp parallel default(none) &
 !$omp private(ray_dir) &
 !$omp shared(nrays,nsides,primary,kappa,xyzh,Rstar,Rinject,rays_dist,rays_tau,rays_dim)
 !$omp do
 do i = 1, nrays
    !returns ray_dir, the unit vector identifying a ray (index i-1 because healpix starts counting from index 0)
    call pix2vec_nest(nsides, i-1, ray_dir)
    !calculate the properties along the ray (tau, distance, number of points)
    call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,Rinject,rays_tau(:,i),rays_dist(:,i),rays_dim(i))
 enddo
 !$omp enddo
 !$omp end parallel


 !_----------------------------------------------
 ! DETERMINE the optical depth for each particle
 ! using the values available on the HEALPix rays
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
 !  Calculate the optical depth at each particle's location, using the uniform outward
 !  ray-tracing scheme for models containing a primary star and a companion
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the opacity of all the SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: Rinject:         The particles injection radius
 !  IN: companion:       The xyz coordinates of the companion
 !  IN: Rcomp:           The radius of the companion
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !+
 !  OUT: tau:            The array of optical depths for each SPH particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_companion(npart, primary, Rstar, xyzh, kappa, Rinject, companion, Rcomp, order, tau)
 use part, only : isdead_or_accreted
 integer, intent(in) :: npart, order
 real, intent(in)    :: primary(3), companion(3), kappa(:), Rstar, Rinject, xyzh(:,:), Rcomp
 real, intent(out)   :: tau(:)

 integer  :: i, nrays, nsides
 real     :: normCompanion,theta0,phi,cosphi,sinphi,theta,sep,root
 real     :: ray_dir(3),part_dir(3),uvecCompanion(3)
 real, dimension(:,:), allocatable  :: rays_dist, rays_tau
 integer, dimension(:), allocatable :: rays_dim
 integer, parameter :: ndim = 200 ! maximum number of points along the ray where tau is calculated

 nrays = 12*4**order ! The number of rays traced given the healpix order
 nsides = 2**order   ! The healpix nsides given the healpix order

 allocate(rays_dist(ndim, nrays)) ! distance from the central star of the points on the rays
 allocate(rays_tau(ndim, nrays))  ! value of tau at each point along each ray
 allocate(rays_dim(nrays))        ! effective number of points on the ray (< ndim)

 uvecCompanion = companion-primary
 normCompanion = norm2(uvecCompanion)
 uvecCompanion = uvecCompanion/normCompanion
 theta0        = asin(Rcomp/normCompanion)
 phi           = atan2(uvecCompanion(2),uvecCompanion(1))
 cosphi        = cos(phi)
 sinphi        = sin(phi)

 !-------------------------------------------
 ! CONSTRUCT the RAYS given the HEALPix ORDER
 ! and determine the optical depth along them
 !-------------------------------------------

 !$omp parallel default(none) &
 !$omp private(ray_dir,theta,root,sep) &
 !$omp shared(nrays,nsides,primary,kappa,xyzh,Rstar,Rinject,Rcomp,rays_dist,rays_tau,rays_dim) &
 !$omp shared(uvecCompanion,normCompanion,cosphi,sinphi,theta0)
 !$omp do
 do i = 1, nrays
    !returns ray_dir, the unit vector identifying a ray (index i-1 because healpix starts counting from index 0)
    call pix2vec_nest(nsides, i-1, ray_dir)
    !rotate ray vectors by an angle = phi so the main axis points to the companion (This is because along the
    !main axis (1,0,0) rays are distributed more uniformally
    ray_dir = (/cosphi*ray_dir(1) - sinphi*ray_dir(2),sinphi*ray_dir(1) + cosphi*ray_dir(2), ray_dir(3)/)
    theta   = acos(dot_product(uvecCompanion, ray_dir))
    !the ray intersects the companion: only calculate tau up to the companion
    if (theta < theta0) then
       root  = sqrt(Rcomp**2-normCompanion**2*sin(theta)**2)
       sep   = normCompanion*cos(theta)-root
       call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,Rinject,rays_tau(:,i),rays_dist(:,i),rays_dim(i), sep)
    else
       call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,Rinject,rays_tau(:,i),rays_dist(:,i),rays_dim(i))
    endif
 enddo
 !$omp enddo
 !$omp end parallel

 !-----------------------------------------------
 ! DETERMINE the optical depth for each particle
 ! using the values available on the HEALPix rays
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
 !  Calculate the optical depth at the SPH particle's location.
 !  Search for the four closest rays to a particle, perform four-point
 !  interpolation of the optical depth from these rays. Weighted by the
 !  inverse square of the perpendicular distance to the rays.
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: nsides:          The healpix nsides of the simulation
 !  IN: vec:             The vector from the primary to the particle
 !  IN: rays_tau:        2-dimensional array containing the cumulative optical
 !                       depth along each ray
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
 !returns rayIndex, the index of the ray vector of the HEALPix cell that points to the particle (direction vec)
 call vec2pix_nest(nsides, vec, rayIndex)
 !returns ray(3), the unit vector identifying the ray with index number rayIndex
 call pix2vec_nest(nsides, rayIndex, ray)
 !compute optical depth along ray rayIndex(+1)
 call get_tau_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), tautemp)
 !determine distance of the particle to the HEALPix ray
 vectemp       = vec - vec_norm2*ray
 distRay_sq    = dot_product(vectemp,vectemp)
 if (distRay_sq > 0.) then
    tau    = tautemp/distRay_sq
    weight = 1./distRay_sq
 else
    ! the particle sits exactly on the ray, no need to interpolate with the neighbours
    tau    = tautemp
    return
 endif

 !returns the number nneigh and list of vectors (n) neighbouring the ray number rayIndex
 call neighbours_nest(nsides, rayIndex, neighbours, nneigh)
 !for each neighbouring ray calculate its distance to the particle
 do i=1,nneigh
    call pix2vec_nest(nsides, neighbours(i), ray)
    vectemp     = vec - vec_norm2*ray
    tempdist(i) = dot_product(vectemp,vectemp)
 enddo
 neighbours     = neighbours+1
 mask           = .true.
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
 !  at a given distance to the starting point of the ray (primary star).
 !+
 !  IN: distance:        The distance from the staring point of the ray to a
 !                       point on the ray
 !  IN: tau_along_ray:   The vector of cumulative optical depths along the ray
 !  IN: dist_along_ray:  The vector of distances from the primary along the ray
 !  IN: len:             The length of tau_along_ray and dist_along_ray
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
    tau = tau_along_ray(1)
 elseif (distance  >  dist_along_ray(len)) then
    tau = tau_along_ray(len)
 else
    L = 2
    R = len
    !bysection search for the index of the closest points on the ray to the specified location
    do while (L < R)
       m = (L + R)/2
       if (dist_along_ray(m) > distance) then
          R = m
       else
          L = m + 1
       endif
    enddo
    !linear interpolation of the optical depth at the the point's location
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
 !                       optical depth will be calculated
 !  IN: xyzh:            The array containing the particles position+smoothing lenght
 !  IN: kappa:           The array containing the particles opacity
 !  IN: Rstar:           The radius of the primary star
 !  IN: Rinject:         The particles injection radius
 !+
 !  OUT: tau_along_ray:  The vector of cumulative optical depth along the ray
 !  OUT: dist_along_ray: The vector of distances from the primary along the ray
 !  OUT: len:            The length of tau_along_ray and dist_along_ray
 !+
 !  OPT: maxDistance:    The maximal distance the ray needs to be traced
 !+
 !--------------------------------------------------------------------------
subroutine ray_tracer(primary, ray, xyzh, kappa, Rstar, Rinject, tau_along_ray, dist_along_ray, len, maxDistance)
 use units, only:unit_opacity
 use part,  only:itauL_alloc
 real, intent(in)     :: primary(3), ray(3), Rstar, Rinject, xyzh(:,:), kappa(:)
 real, optional       :: maxDistance
 real, intent(out)    :: dist_along_ray(:), tau_along_ray(:)
 integer, intent(out) :: len
 real, parameter      :: tau_max = 99.

 real    :: dr, next_dr, h, dtaudr, previousdtaudr, nextdtaudr, distance
 integer :: inext, i, L, R, m ! left, right and middle index for binary search

 h = Rinject/100.
 inext=0
 do while (inext==0)
    h = h*2.
    !find the next point along the ray : index inext
    call find_next(primary+Rinject*ray, h, ray, xyzh, kappa, previousdtaudr, dr, inext)
 enddo

 i = 1
 tau_along_ray(i)  = 0.
 distance          = Rinject
 dist_along_ray(i) = distance
 do while (hasNext(inext,tau_along_ray(i),distance,maxDistance))
    distance = distance+dr
    call find_next(primary + distance*ray, xyzh(4,inext), ray, xyzh, kappa, nextdtaudr, next_dr, inext)
    i = i + 1
    if (itauL_alloc > 0) nextdtaudr = nextdtaudr*(Rstar/distance)**2
    dtaudr            = (nextdtaudr+previousdtaudr)/2.
    previousdtaudr    = nextdtaudr
    !fix units for tau (kappa is in cgs while rho & r are in code units)
    tau_along_ray(i)  = tau_along_ray(i-1) + real(dr*dtaudr/unit_opacity)
    dist_along_ray(i) = distance
    dr                = next_dr
 enddo

 if (itauL_alloc == 0 .and. present(maxDistance)) then
    i = i + 1
    tau_along_ray(i)  = tau_max
    dist_along_ray(i) = maxDistance
 endif
 len = i

 if (itauL_alloc > 0) then
    !reverse integration start from zero inward
    tau_along_ray(1:len) = tau_along_ray(len) - tau_along_ray(1:len)
    !find the first point where tau_lucy < 2/3
    if (tau_along_ray(1) > 2./3.) then
       L = 1
       R = len
       !bysection search for the index of the closest point to tau = 2/3
       do while (L < R)
          m = (L + R)/2
          if (tau_along_ray(m) < 2./3.) then
             R = m
          else
             L = m + 1
          endif
       enddo
       tau_along_ray(1:L) = 2./3.
       !The photosphere is located between ray grid point L and L+1, may be useful information!
    endif
 endif
end subroutine ray_tracer

logical function hasNext(inext, tau, distance, maxDistance)
 integer, intent(in) :: inext
 real, intent(in)    :: tau, distance
 real, optional      :: maxDistance
 real                :: tau_max = 99.
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
 !  IN: inpoint:         The coordinate of the initial point projected on the ray
 !                       for which the opacity and the next point will be calculated
 !  IN: h:               The smoothing length at the initial point
 !  IN: ray:             The unit vector of the direction in which the next
 !                       point will be calculated
 !  IN: xyzh:            The array containing the particles position+smoothing length
 !  IN: kappa:           The array containing the particles opacity
 !  IN: inext:           The index of the initial point
 !                       (this point will not be considered as possible next point)
 !+
 !  OUT: dtaudr:         The radial optical depth derivative at the given location (inpoint)
 !  OUT: distance:       The distance to the next point
 !  OUT: inext:          The index of the next point on the ray
 !+
 !--------------------------------------------------------------------------
subroutine find_next(inpoint, h, ray, xyzh, kappa, dtaudr, distance, inext)
 use linklist, only:getneigh_pos,ifirstincell,listneigh
 use kernel,   only:radkern,cnormk,wkern
 use part,     only:hfact,rhoh,massoftype,igas
 use dim,      only:maxpsph
 real,    intent(in)    :: xyzh(:,:), kappa(:), inpoint(:), ray(:), h
 integer, intent(inout) :: inext
 real,    intent(out)   :: distance, dtaudr

 integer, parameter :: nmaxcache = 0
 real  :: xyzcache(0,nmaxcache)

 integer  :: nneigh, i, prev,j
 real     :: dmin, vec(3), dr, raydistance, q, norm_sq

 prev     = inext
 inext    = 0
 distance = 0.

 !for a given point (inpoint), returns the list of neighbouring particles (listneigh) within a radius h*radkern
 call getneigh_pos(inpoint,0.,h*radkern,3,listneigh,nneigh,xyzh,xyzcache,nmaxcache,ifirstincell)

 dtaudr = 0.
 dmin = huge(0.)
 !loop over all neighbours
 do i=1,nneigh
    j = listneigh(i)
    if (j > maxpsph) cycle
    vec     = xyzh(1:3,j) - inpoint
    norm_sq = dot_product(vec,vec)
    q       = sqrt(norm_sq)/xyzh(4,j)
    !add optical depth contribution from each particle
    dtaudr = dtaudr+wkern(q*q,q)*kappa(j)*rhoh(xyzh(4,j), massoftype(igas))

    ! find the next particle : among the neighbours find the particle located the closest to the ray
    if (j  /=  prev) then
       dr = dot_product(vec,ray) !projected distance along the ray
       if (dr>0.) then
          !distance perpendicular to the ray direction
          raydistance = norm_sq - dr**2
          if (raydistance < dmin) then
             dmin     = raydistance
             inext    = j
             distance = dr
          endif
       endif
    endif
 enddo
 dtaudr = dtaudr*cnormk/hfact**3
end subroutine find_next
end module raytracer
