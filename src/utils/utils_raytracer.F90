module raytracer
   use healpix
   implicit none
   public :: get_all_tau_outwards
   private
 contains

   !*********************************************************************!
   !***************************   OUTWARDS   ****************************!
   !*********************************************************************!

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the uniform outwards
   !  ray-tracing scheme
   !+
   !  IN: primary:         coordinates of the primary star
   !  IN: xyzh:            array containing the particles position+smooting lenght
   !  IN: opacities:       array containing the particles opacity
   !  IN: Rstar:           The radius of the star
   !  IN: order:           The order in which the rays are sampled
   !+
   !  OUT: taus:           array containing the particles optical depth
   !+
   !  OPT: companion:      companion coordinates
   !  OPT: Rcomp:          companion radius
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards(primary, xyzh, opacities, Rstar, order, taus, companion, Rcomp)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), opacities(:), Rstar, xyzh(:,:)
      real, optional      :: Rcomp, companion(3)
      real, intent(out)   :: taus(:)

      if (present(companion) .and. present(Rcomp)) then
         call get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, companion, Rcomp, order, taus)
      else
         call get_all_tau_outwards_single(primary, xyzh, opacities, Rstar, order, taus)
      endif
   end subroutine get_all_tau_outwards

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the uniform outwards
   !  ray-tracing scheme in case only a single star is present
   !+
   !  IN: primary:         coordinates of the primary star
   !  IN: xyzh:            array containing the particles position+smoothing length
   !  IN: opacities:       array containing the particles opacity
   !  IN: Rstar:           primary star radius
   !  IN: order:           The order of the HEALPIX method
   !+
   !  OUT: taus:           array containing the particles optical depth
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards_single(primary, xyzh, opacities, Rstar, order, taus)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), opacities(:), Rstar, xyzh(:,:)
      real, intent(out)   :: taus(:)

      integer  :: i, nrays, nsides, len, npart
      real     :: ray_dir(3)
      real, dimension(:,:),  allocatable :: ray_dists, ray_taus
      real, dimension(:),    allocatable :: tau, dist
      integer, dimension(:), allocatable :: ray_dims

      integer, parameter :: ndim = 200
      nrays = 12*4**order
      nsides = 2**order

      allocate(ray_dists(ndim, nrays))
      allocate(ray_taus(ndim, nrays))
      allocate(ray_dims(nrays))
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
         ray_taus(:,i)  = tau
         ray_dists(:,i) = dist
         ray_dims(i)    = len
      enddo
      !$omp end parallel do


      !_----------------------------------------------
      ! DETERMINE the optical depth for each particle
      ! using the values available on the rays
      !-----------------------------------------------

      taus  = 0.
      npart = size(taus)
      !$omp parallel do private(ray_dir)
      do i = 1,npart
         ray_dir = xyzh(1:3,i)-primary
         call ray_polation(nsides, ray_dir, ray_taus, ray_dists, ray_dims, taus(i))
      enddo
      !$omp end parallel do
   end subroutine get_all_tau_outwards_single


   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the uniform outwards
   !  ray-tracing scheme concerning a binary system
   !+
   !  IN: primary:         coordinates of the primary star
   !  IN: xyzh:            array containing the particles position+smoothing length
   !  IN: opacities:       array containing the particles opacity
   !  IN: Rstar:           primary star radius
   !  IN: order:           The order of the HEALPIX method
   !  IN: companion:       coordinates of the companion star
   !  IN: Rcomp:           companion star radius
   !+
   !  OUT: taus:           array containing the particles optical depth
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, companion, Rcomp, order, taus)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), companion(3), opacities(:), Rstar, xyzh(:,:), Rcomp
      real, intent(out)   :: taus(:)

      integer  :: i, nrays, nsides, len, npart
      real     :: normCompanion = 0., theta0 = 0., unitCompanion(3) = 0.
      real     :: theta, sep, root, ray_dir(3), phi = 0., cosphi = 0., sinphi = 0.
      real, dimension(:,:), allocatable  :: ray_dists, ray_taus
      real, dimension(:), allocatable    :: tau, dist
      integer, dimension(:), allocatable :: ray_dims

      integer, parameter :: ndim = 200
      nrays = 12*4**order
      nsides = 2**order

      allocate(ray_dists(ndim, nrays))
      allocate(ray_taus(ndim, nrays))
      allocate(ray_dims(nrays))
      allocate(tau(ndim))
      allocate(dist(ndim))

      unitCompanion = companion-primary
      normCompanion = norm2(unitCompanion)
      theta0        = asin(Rcomp/normCompanion)
      unitCompanion = unitCompanion/normCompanion
      phi    = atan2(unitCompanion(2),unitCompanion(1))
      cosphi = cos(phi)
      sinphi = sin(phi)

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
         theta = acos(dot_product(unitCompanion, ray_dir))
         !the ray intersects the companion: only calculate tau up to the companion
         if (theta < theta0) then
            !root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+Rcomp**2)
            root = sqrt(Rcomp**2-normCompanion**2*sin(theta)**2)
            sep   = normCompanion*cos(theta)-root
            call ray_tracer(primary, ray_dir, xyzh, opacities, Rstar, tau, dist, len, sep)
         else
            call ray_tracer(primary, ray_dir, xyzh, opacities, Rstar, tau, dist, len)
         endif
         ray_taus(:,i)  = tau
         ray_dists(:,i) = dist
         ray_dims(i)    = len
      enddo
      !$omp end parallel do


      !-----------------------------------------------
      ! DETERMINE the optical depth for each particle
      ! using the values available on the rays
      !-----------------------------------------------

      taus = 0.
      npart = size(taus)
      !$omp parallel do private(ray_dir)
      do i = 1, npart
         ray_dir = xyzh(1:3,i)-primary
         ray_dir = (/cosphi*ray_dir(1) + sinphi*ray_dir(2),-sinphi*ray_dir(1) + cosphi*ray_dir(2), ray_dir(3)/)
         call ray_polation(nsides, ray_dir, ray_taus, ray_dists, ray_dims, taus(i))
      enddo
      !$omp end parallel do
    end subroutine get_all_tau_outwards_companion

   !--------------------------------------------------------------------------
   !+
   !  Calculate the particle's optical using the rays
   !+
   !  IN: nsides:          the healpix nsides of the simulation
   !  IN: vec:             the vector originating from the primary star and pointing to the particle
   !  IN: ray_taus:        array containing the optical depths at each locations an the rays
   !  IN: ray_dists:       array containing the locations along the ray where tau has been calculated
   !  IN: ray_dims:        array containing the number of point along the ray where tau has been calculated
   !+
   !  OUT: tau:            interpolated particle's optical depth
   !+
   !--------------------------------------------------------------------------
   subroutine ray_polation(nsides, vec, ray_taus, ray_dists, ray_dims, tau)
      integer, intent(in) :: nsides, ray_dims(:)
      real, intent(in)    :: vec(:), ray_dists(:,:), ray_taus(:,:)
      real, intent(out)   :: tau

      integer :: index, n(8), nneigh, i
      real    :: tau_ray, ray(3), vectemp(3), suminvdist2, tempdist(8), sq_dist2ray, vec_norm2
      logical :: mk(8)

      vec_norm2 = norm2(vec)

      !returns index, the index of the ray vector that points to the particle (direction vec)
      call vec2pix_nest(nsides, vec, index)
      !returns ray(3), the unit vector identifying the ray number index
      call pix2vec_nest(nsides, index, ray)
      vectemp     = vec - vec_norm2*ray
      sq_dist2ray = dot_product(vectemp,vectemp) !LS it is an approximation. It should be the length of the arc :
                                                 !vec_norm2*angle_between_ray_and_vec but I guess it doesn't mak any difference
      call get_tau_on_ray(vec_norm2, ray_taus(:,index+1), ray_dists(:,index+1), ray_dims(index+1), tau_ray)
      tau         = tau_ray/sq_dist2ray
      suminvdist2 = 1./sq_dist2ray

      !returns the number nneigh and list of vectors (n) neighbouring the ray number index
      call neighbours_nest(nsides, index, n, nneigh)
      !for each neighbouring ray calculate its distance to the particle
      do i=1,nneigh
         call pix2vec_nest(nsides, n(i), ray)
         vectemp     = vec - vec_norm2*ray
         tempdist(i) = dot_product(vectemp,vectemp)
      enddo
      n=n+1
      mk = .true.
      if (nneigh <8) mk(nneigh+1:8) = .false.
      !take tau contribution from the 3 closest rays
      do i=1,3
         index       = minloc(tempdist,1,mk)
         mk(index)   = .false.
         call get_tau_on_ray(vec_norm2, ray_taus(:,n(index)), ray_dists(:,n(index)), ray_dims(n(index)), tau_ray)
         tau         = tau + tau_ray/tempdist(index)
         suminvdist2 = suminvdist2 + 1./tempdist(index)
      enddo
      tau = tau / suminvdist2
   end subroutine ray_polation


   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth at a given distance on the ray
   !+
   !  IN: dist:            The radial distance of the particle
   !  IN: ray_tau:         array containing the optical depths at each locations along the rays
   !  IN: ray_dist:        array containing the locations along the ray where tau has been calculated
   !  IN: len:             length of ray_tau and ray_dist
   !+
   !  OUT: tau:            the optical depth at distance dist
   !+
   !--------------------------------------------------------------------------
   subroutine get_tau_on_ray(dist, ray_tau, ray_dist, len, tau)
      real, intent(in)    :: dist, ray_tau(:), ray_dist(:)
      integer, intent(in) :: len
      real, intent(out)   :: tau

      integer :: L, R, m

      if (dist .lt. ray_dist(1)) then
         tau = 0.
      else if (dist .gt. ray_dist(len)) then
         tau = 1e10
      else
         L = 2
         R = len-1
         !bysection search for the index of the closest ray location to the particle
         do while (L < R)
            m = int((L + R)/2) !LS missing int
            if (ray_dist(m) > dist) then
                  R = m
            else
                  L = m + 1
            end if
         enddo
         !interpolate linearly ray properties to get the particle's optical depth
         tau = ray_tau(L-1)+(ray_tau(L)-ray_tau(L-1))/(ray_dist(L)-ray_dist(L-1))*(dist-ray_dist(L-1))
      endif
    end subroutine get_tau_on_ray

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth along a given ray
   !+
   !  IN: primary:         The location of the primary star
   !  IN: ray:             The direction of the ray that needs to be traced
   !  IN: xyzh:            array containing the particles position+smoothing length
   !  IN: opacities:       array containing the particles opacity
   !  IN: Rstar:           stellar radius
   !+
   !  OUT: taus:           array containing the particles optical depth
   !  OUT: ray_dist:       array containing the locations along the ray where tau has been calculated
   !  OUT: len:            The number of point along the ray where tau has been calculated = length(ray_dist)
   !+
   !  OPT: maxDist:        The maximal distance the ray needs to be traced
   !+
   !--------------------------------------------------------------------------
   subroutine ray_tracer(primary, ray, xyzh, opacities, Rstar, taus, ray_dist, len, maxDist)
      real, intent(in)     :: primary(3), ray(3), Rstar, xyzh(:,:), opacities(:)
      real, optional       :: maxDist
      real, intent(out)    :: ray_dist(:), taus(:)
      integer, intent(out) :: len

      real    :: dist, nextDist, h, dtau_dr, prev_dtau_dr, next_dtau_dr, totalDist
      integer :: next, i

      h = Rstar/100.
      next=0
      do while (next==0)
         h = h*2.
         !find the next point along the ray : index next
         call find_next(primary+Rstar*ray, h, ray, xyzh, opacities, prev_dtau_dr, dist, next)
      enddo

      i = 1
      taus (1)    = 0.
      totalDist   = Rstar
      ray_dist(i) = TotalDist
      do while (hasNext(next,totalDist,maxDist))
         totalDist = totalDist+dist
         call find_next(primary + totalDist*ray, xyzh(4,next), ray, xyzh, opacities, next_dtau_dr, nextDist, next)
         i = i + 1
         dtau_dr         = (next_dtau_dr+prev_dtau_dr)/2.
         prev_dtau_dr    = next_dtau_dr
         taus(i)         = taus(i-1)+(dist)*dtau_dr
         ray_dist(i)     = TotalDist
         dist            = nextDist
      enddo
      len = i
   end subroutine ray_tracer

   logical function hasNext(next, dist, maxDist)
      integer, intent(in) :: next
      real, intent(in)    :: dist
      real, optional      :: maxDist
      if (present(maxDist)) then
         hasNext = next /= 0 .and. dist .lt. maxDist
      else
         hasNext = next /= 0
      endif
   end function hasNext


   !*********************************************************************!
   !****************************   COMMON   *****************************!
   !*********************************************************************!

   !--------------------------------------------------------------------------
   !+
   !  First finds the opacity at the starting point, then finds the next
   !                       point on a ray and the distance to this point
   !+
   !  IN: inpoint:         The initial point for which the kappa_rho and the
   !                       next point will be calculated
   !  IN: h:               The smoothing length at the initial point
   !  IN: ray:             The unit vector of the direction in which the next
   !                       point will be calculated
   !  IN: xyzh:            array containing the particles position+smoothing length
   !  IN: opacities:       array containing the particles opacity
   !  IN: next:            The index of the initial point
   !                       (this point will not be considered as possible next point)
   !+
   !  OUT: dtau_dr:        optical depth increment (kappa*rho) at the particle's location
   !  OUT: dist:           The distance along the ray to the next point
   !  OUT: next:           The index of the next point on the ray
   !+
   !--------------------------------------------------------------------------
   subroutine find_next(inpoint, h, ray, xyzh, opacities, dtau_dr, dist, next)
      use linklist, only:getneigh_pos,ifirstincell,listneigh
      use kernel,   only:radkern,cnormk,wkern
      use part,     only:hfact, rhoh, massoftype, igas
      real, intent(in)       :: xyzh(:,:), opacities(:), inpoint(:), ray(:), h
      integer, intent(inout) :: next
      real, intent(out)      :: dist, dtau_dr

      integer, parameter :: maxcache = 0
      real, allocatable  :: xyzcache(:,:)

      integer  :: nneigh, i, prev
      real     :: dmin, vec(3), dr_ray, dh_ray, q, norm_sq

      prev=next
      next=0
      dist=0.

      !for a given point (inpoint), returns the list of neighbouring particles (listneigh) within a radius h*radkern
      call getneigh_pos(inpoint,0.,h*radkern, 3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)

      dtau_dr = 0.
      dmin = huge(0.)

      !loop over all neighbours
      do i=1,nneigh
         vec     = xyzh(1:3,listneigh(i)) - inpoint
         norm_sq = dot_product(vec,vec)
         q       = sqrt(norm_sq)/xyzh(4,listneigh(i))
         !add optical depth contribution from each particle
         dtau_dr = dtau_dr+wkern(q*q,q)*opacities(listneigh(i))*rhoh(xyzh(4,listneigh(i)), massoftype(igas))

         ! find the next particle : among the neighbours find the particle located the closest to the ray
         if (listneigh(i) .ne. prev) then
            dr_ray = dot_product(vec,ray) !projected distance along the ray
            if (dr_ray>0.) then
               !distance perpendicular to the ray direction
               dh_ray = norm_sq - dr_ray**2
               if (dh_ray < dmin) then
                  dmin = dh_ray
                  next = listneigh(i)
                  dist = dr_ray
               end if
            end if
         endif
      enddo
      dtau_dr = dtau_dr*cnormk/hfact**3
   end subroutine find_next
end module raytracer
