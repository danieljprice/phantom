!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
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
! :Owner: Esseldeurs Mats
!
! :Runtime parameters:
!   - TODO
!
! :Dependencies: linklist, kernel, part, healpix
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
   !  IN: primary:         The xyz coordinates of the primary star
   !  IN: xyzh:            The xyzh of all the SPH particles
   !  IN: opacities:       The array containing the opacities of all the SPH particles
   !  IN: Rstar:           The radius of the primary star
   !  IN: order:           The healpix order which is used for the uniform ray sampling
   !+
   !  OUT: taus:           The array of optical depths to each SPH particle
   !+
   !  OPT: companion:      The xyz coordinates of the companion
   !  OPT: Rcomp:          The radius of the companion
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau(primary, xyzh, opacities, Rstar, order, taus, companion, Rcomp)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), opacities(:), Rstar, xyzh(:,:)
      real, optional      :: Rcomp, companion(3)
      real, intent(out)   :: taus(:)

      if (present(companion) .and. present(Rcomp)) then
         call get_all_tau_companion(primary, xyzh, opacities, Rstar, companion, Rcomp, order, taus)
      else
         call get_all_tau_single(primary, xyzh, opacities, Rstar, order, taus)
      endif
   end subroutine get_all_tau
   
   !--------------------------------------------------------------------------
   !+
   !  Calculates the optical depth to each SPH particle, using the uniform outwards
   !  ray-tracing scheme for models containing a single star
   !
   !  Relies on healpix, for more information: https://healpix.sourceforge.io/
   !+
   !  IN: primary:         The xyz coordinates of the primary star
   !  IN: xyzh:            The xyzh of all the SPH particles
   !  IN: opacities:       The array containing the opacities of all the SPH particles
   !  IN: Rstar:           The radius of the primary star
   !  IN: order:           The healpix order which is used for the uniform ray sampling
   !+
   !  OUT: taus:           The array of optical depths to each SPH particle
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_single(primary, xyzh, opacities, Rstar, order, taus)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), opacities(:), Rstar, xyzh(:,:)
      real, intent(out)   :: taus(:)
      
      integer  :: i, nrays, nsides, len, npart
      real     :: ray_dir(3)
      real, dimension(:,:), allocatable  :: AllDistances, AllTaus
      real, dimension(:), allocatable    :: tau, dist
      integer, dimension(:), allocatable :: listOfLens
      integer, parameter :: ndim = 200
      nrays = 12*4**order ! The number of rays traced given the healpix order
      nsides = 2**order   ! The healpix nsides given the healpix order
      
      allocate(AllDistances(ndim, nrays))
      allocate(AllTaus(ndim, nrays))
      allocate(listOfLens(nrays))
      allocate(tau(ndim))
      allocate(dist(ndim))
     
      !-------------------------------------------
      ! CONSTRUCT the RAYS given the ORDER
      ! and determine the optical depth along them
      !-------------------------------------------
      
      !$omp parallel do private(ray_dir,tau,dist,len)
      do i = 1, nrays
         tau   = 0.
         dist  = 0.
         !returns ray_dir(3), the unit vector identifying ray number i-1; i-1 as healpix starts counting at 0 and fortran at 1
         call pix2vec_nest(nsides, i-1, ray_dir)
         !calculate the optical depth along a given ray, len is the number of points on the ray
         call ray_tracer(primary, ray_dir, xyzh, opacities, Rstar, tau, dist, len)
         AllTaus(:,i)      = tau
         AllDistances(:,i) = dist
         listOfLens        = len
      enddo
      !$omp end parallel do

      !_----------------------------------------------
      ! DETERMINE the optical depth for each particle
      ! using the values available on the rays
      !-----------------------------------------------

      taus  = 0.
      npart = size(taus)
      !$omp parallel do private(dir)
      do i = 1, npart
         ray_dir = xyzh(1:3,i)-primary
         call interpolate_tau(nsides, ray_dir, AllTaus, AllDistances, listOfLens, taus(i))
      enddo
      !$omp end parallel do
   end subroutine get_all_tau_single

   !--------------------------------------------------------------------------
   !+
   !  Calculates the optical depth to each SPH particle, using the uniform outwards
   !  ray-tracing scheme for models containing primary star and a companion
   !
   !  Relies on healpix, for more information: https://healpix.sourceforge.io/
   !+
   !  IN: primary:         The xyz coordinates of the primary star
   !  IN: xyzh:            The xyzh of all the SPH particles
   !  IN: opacities:       The array containing the opacities of all the SPH particles
   !  IN: Rstar:           The radius of the primary star
   !  IN: order:           The healpix order which is used for the uniform ray sampling
   !  IN: companion:       The xyz coordinates of the companion
   !  IN: Rcomp:           The radius of the companion
   !+
   !  OUT: taus:           The array of optical depths to each SPH particle
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_companion(primary, xyzh, opacities, Rstar, companion, Rcomp, order, taus)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), companion(3), opacities(:), Rstar, xyzh(:,:), Rcomp
      real, intent(out)   :: taus(:)
      
      integer  :: i, nrays, nsides, len, npart
      real     :: normCompanion = 0., theta0 = 0., uvecCompanion(3) = 0.
      real     :: theta, sep, root, ray_dir(3), phi = 0., cosphi = 0., sinphi = 0.
      real, dimension(:,:), allocatable  :: dirs
      real, dimension(:,:), allocatable  :: AllDistances, AllTaus
      real, dimension(:), allocatable    :: tau, dist
      integer, dimension(:), allocatable :: listOfLens
      integer, parameter :: ndim = 200
      nrays = 12*4**order ! The number of rays traced given the healpix order
      nsides = 2**order   ! The healpix nsides given the healpix order
      
      allocate(dirs(3, nrays))
      allocate(AllDistances(ndim, nrays))
      allocate(AllTaus(ndim, nrays))
      allocate(listOfLens(nrays))
      allocate(tau(ndim))
      allocate(dist(ndim))

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

      !$omp parallel do private(ray_dir,tau,dist,sep,root,theta,len)
      do i = 1, nrays
         tau   = 0.
         dist  = 0.
         !returns ray_dir, the unit vector identifying ray number i-1; i-1 as healpix starts counting at 0 and fortran at 1
         call pix2vec_nest(nsides, i-1, ray_dir)
         ray_dir = (/cosphi*ray_dir(1) - sinphi*ray_dir(2),sinphi*ray_dir(1) + cosphi*ray_dir(2), ray_dir(3)/)
         theta   = acos(dot_product(uvecCompanion, ray_dir))
         !the ray intersects the companion: only calculate tau up to the companion
         if (theta < theta0) then
            root  = sqrt(Rcomp**2-normCompanion**2*sin(theta)**2)
            sep   = normCompanion*cos(theta)-root
            call ray_tracer(primary, ray_dir, xyzh, opacities, Rstar, tau, dist, len, sep)
         else
            call ray_tracer(primary, ray_dir, xyzh, opacities, Rstar, tau, dist, len)
         endif
         AllTaus(:,i)      = tau
         AllDistances(:,i) = dist
         listOfLens(i)     = len
      enddo
      !$omp end parallel do

      !-----------------------------------------------
      ! DETERMINE the optical depth for each particle
      ! using the values available on the rays
      !-----------------------------------------------

      taus  = 0.
      npart = size(taus)
      !$omp parallel do private(ray_dir)
      do i = 1, npart
         ray_dir = xyzh(1:3,i)-primary
         ray_dir = (/cosphi*ray_dir(1) + sinphi*ray_dir(2),-sinphi*ray_dir(1) + cosphi*ray_dir(2), ray_dir(3)/)
         call interpolate_tau(nsides, ray_dir, AllTaus, AllDistances, listOfLens, taus(i))
      enddo
      !$omp end parallel do
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
   !  IN: AllTaus:         2-dimensional array containing the cumulative optical
   !                       depts along each ray
   !  IN: AllDistances:    2-dimensional array containing the distances from the
   !                       primary along each ray
   !  IN: listOflens:      1-dimensional array where the element on each index corresponds to
   !                       the length of the array in AllDistances on this index 
   !+
   !  OUT: tau:            The interpolated optical depth
   !+
   !--------------------------------------------------------------------------
   subroutine interpolate_tau(nsides, vec, AllTaus, AllDistances, listOfLens, tau)
      integer, intent(in) :: nsides, listOfLens(:)
      real, intent(in)    :: vec(:), AllDistances(:,:), AllTaus(:,:)
      real, intent(out)   :: tau
   
      integer :: rayIndex, neighIndex, neighbours(8), nneigh, i
      real    :: tautemp, ray(3), vectemp(3), weight, tempdist(8), distRay, vec_norm2
      logical :: mask(8)

      vec_norm2 = norm2(vec)
      !returns index, the index of the ray vector that points to the particle (direction vec)
      call vec2pix_nest(nsides, vec, rayIndex)
      !returns ray(3), the unit vector identifying the ray number index
      call pix2vec_nest(nsides, rayIndex, ray)
      vectemp       = vec - vec_norm2*ray
      distRay = dot_product(vectemp,vectemp)
      call get_tau_on_ray(vec_norm2, AllTaus(:,rayIndex+1), AllDistances(:,rayIndex+1), listOfLens(rayIndex+1), tautemp)
      tau           = tautemp/distRay
      weight        = 1./distRay

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
         neighIndex       = minloc(tempdist,1,mask)
         mask(neighIndex) = .false.
         call get_tau_on_ray(norm2(vec), AllTaus(:,neighbours(neighIndex)), &
                AllDistances(:,neighbours(neighIndex)), listOfLens(neighbours(neighIndex)), tautemp)
         tau    = tau + tautemp/tempdist(neighIndex)
         weight = weight + 1./tempdist(neighIndex)
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
   !  IN: listOfTaus:      The array of cumulative optical depts along the ray
   !  IN: listOfDistances: The array of distances from the primary along the ray
   !  IN: len:             The length of listOfTaus and listOfDists
   !+
   !  OUT: tau:            The optical depth to the given distance along the ray
   !+
   !--------------------------------------------------------------------------
   subroutine get_tau_on_ray(distance, listOfTaus, listOfDistances, len, tau)
      real, intent(in)    :: distance, listOfTaus(:), listOfDistances(:)
      integer, intent(in) :: len
      real, intent(out)   :: tau
   
      integer :: L, R, m ! left, right and middle index for binary search
   
      if (distance .lt. listOfDistances(1)) then
         tau = 0.
      else if (distance .gt. listOfDistances(len)) then
         tau = 1e10
      else
         L = 2
         R = len-1
         !bysection search for the index of the closest ray location to the particle
         do while (L < R)
            m = (L + R)/2
            if (listOfDistances(m) > distance) then
                  R = m
            else
                  L = m + 1
            end if
         enddo
         !interpolate linearly ray properties to get the particle's optical depth
         tau = listOfTaus(L-1)+(listOfTaus(L)-listOfTaus(L-1))/ &
                  (listOfDistances(L)-listOfDistances(L-1))*(distance-listOfDistances(L-1))
      endif
   end subroutine get_tau_on_ray
    
   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depts along a given ray
   !+
   !  IN: primary:         The location of the primary star
   !  IN: ray:             The unit vector of the direction in which the 
   !                       optical depts will be calculated
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: opacities:       The list of the opacities of the particles
   !  IN: Rstar:           The radius of the star
   !+
   !  OUT: listOfTaus:     The array of cumulative optical depts along the ray
   !  OUT: listOfDistances:The array of distances from the primary along the ray
   !  OUT: len:            The length of listOfTaus and listOfDistances
   !+
   !  OPT: maxDistance:    The maximal distance the ray needs to be traced
   !+
   !--------------------------------------------------------------------------
   subroutine ray_tracer(primary, ray, xyzh, opacities, Rstar, listOfTaus, listOfDistances, len, maxDistance)
      real, intent(in)     :: primary(3), ray(3), Rstar, xyzh(:,:), opacities(:)
      real, optional       :: maxDistance
      real, intent(out)    :: listOfDistances(:), listOfTaus(:)
      integer, intent(out) :: len
      
      real    :: distance, nextDistance, h, dtaudr, previousdtaudr, nextdtaudr, totalDistance
      integer :: Inext, i

      h = Rstar/100.
      Inext=0
      do while (Inext==0)
         h = h*2.
         !find the next point along the ray : index next
         call find_next(primary+Rstar*ray, h, ray, xyzh, opacities, previousdtaudr, distance, Inext)
      enddo

      i = 1
      listOfTaus(1)      = 0.
      totalDistance      = Rstar
      listOfDistances(i) = TotalDistance
      do while (hasNext(Inext,totalDistance,maxDistance))
         totalDistance = totalDistance+distance
         call find_next(primary + totalDistance*ray, xyzh(4,Inext), ray, xyzh, opacities, nextdtaudr, nextDistance, Inext)
         i = i + 1
         dtaudr             = (nextdtaudr+previousdtaudr)/2.
         previousdtaudr     = nextdtaudr
         listOfTaus(i)      = listOfTaus(i-1)+(distance)*dtaudr
         listOfDistances(i) = TotalDistance
         distance           = nextDistance
      enddo
      len = i
   end subroutine ray_tracer

   logical function hasNext(Inext, distance, maxDistance)
      integer, intent(in) :: Inext
      real, intent(in)    :: distance
      real, optional      :: maxDistance
      if (present(maxDistance)) then
         hasNext = Inext /= 0 .and. distance .lt. maxDistance
      else
         hasNext = Inext /= 0
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
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: opacities:       The list of the opacities of the particles
   !  IN: Inext:           The index of the initial point
   !                       (this point will not be considered as possible next point)
   !+
   !  OUT: dtaudr:         The local optical depth derivative at the given location (inpoint)
   !  OUT: distance:       The distance to the next point
   !  OUT: Inext:          The next point on the ray
   !+
   !--------------------------------------------------------------------------
   subroutine find_next(inpoint, h, ray, xyzh, opacities, dtaudr, distance, Inext)
      use linklist, only:getneigh_pos,ifirstincell,listneigh
      use kernel,   only:radkern,cnormk,wkern
      use part,     only:hfact,rhoh,massoftype,igas
      real, intent(in)       :: xyzh(:,:), opacities(:), inpoint(:), ray(:), h
      integer, intent(inout) :: Inext
      real, intent(out)      :: distance, dtaudr

      integer, parameter :: maxcache = 0
      real, allocatable  :: xyzcache(:,:)
      
      integer  :: nneigh, i, prev
      real     :: dmin, vec(3), dr, raydistance, q, norm_sq

      prev     = Inext
      Inext    = 0
      distance = 0.

      !for a given point (inpoint), returns the list of neighbouring particles (listneigh) within a radius h*radkern
      call getneigh_pos(inpoint,0.,h*radkern,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)

      dtaudr = 0
      dmin = huge(0.)
      !loop over all neighbours
      do i=1,nneigh
         vec     = xyzh(1:3,listneigh(i)) - inpoint
         norm_sq   = dot_product(vec,vec)
         q       = sqrt(norm_sq)/xyzh(4,listneigh(i))
         !add optical depth contribution from each particle
         dtaudr = dtaudr+wkern(q*q,q)*opacities(listneigh(i))*rhoh(xyzh(4,listneigh(i)), massoftype(igas))

                  ! find the next particle : among the neighbours find the particle located the closest to the ray
if (listneigh(i) .ne. prev) then
            dr = dot_product(vec,ray) !projected distance along the ray
            if (dr>0.) then
               !distance perpendicular to the ray direction
               raydistance = norm_sq - dr**2
               if (raydistance < dmin) then
                  dmin     = raydistance
                  Inext    = listneigh(i)
                  distance = dr
               end if
            end if
         endif
      enddo
      dtaudr = dtaudr*cnormk/hfact**3
   end subroutine find_next
end module raytracer