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
   !  IN: primary:         The location of the primary star
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: opacities:       The list of the opacities of the particles
   !  IN: Rstar:           The radius of the star
   !  IN: order:           The order in which the rays are sampled
   !+
   !  OUT: taus:           The list of optical depths for each particle
   !+
   !  OPT: companion:      The location of the companion
   !  OPT: R:              The radius of the companion
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards(primary, xyzh, opacities, Rstar, order, taus, companion, R)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), opacities(:), Rstar, xyzh(:,:)
      real, optional      :: R, companion(3)
      real, intent(out)   :: taus(:)

      if (present(companion) .and. present(R)) then
         call get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, order, taus, companion, R)
      else
         call get_all_tau_outwards_single(primary, xyzh, opacities, Rstar, order, taus)
      endif
   end subroutine get_all_tau_outwards
   
   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the uniform outwards
   !  ray-tracing scheme concerning only a single star
   !+
   !  IN: primary:         The location of the primary star
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: opacities:       The list of the opacities of the particles
   !  IN: Rstar:           The radius of the star
   !  IN: order:           The order in which the rays are sampled
   !+
   !  OUT: taus:           The list of optical depths for each particle
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards_single(primary, xyzh, opacities, Rstar, order, taus)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), opacities(:), Rstar, xyzh(:,:)
      real, intent(out)   :: taus(:)
      
      integer  :: i, nrays, nsides, len
      real     :: dir(3)
      real, dimension(:,:), allocatable  :: dirs, listsOfDists, listsOfTaus
      real, dimension(:), allocatable    :: tau, dists
      integer, dimension(:), allocatable :: listOfLens
      nrays = 12*4**order
      nsides = 2**order
      
      allocate(dirs(3, nrays))
      allocate(listsOfDists(200, nrays))
      allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
      allocate(listOfLens(nrays))
      allocate(tau(size(listsOfDists(:,1))))
      allocate(dists(size(listsOfDists(:,1))))
     
      !$omp parallel do private(dir,tau,dists,len)
      do i = 1, nrays
         tau=0.
         dists=0.
         call pix2vec_nest(nsides, i-1, dir)
         dirs(:,i) = dir
         call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, len)
         listsOfTaus(:,i) = tau
         listsOfDists(:,i) = dists
         listOfLens = len
      enddo
      !$omp end parallel do

      taus = 0.
      !$omp parallel do private(dir)
      do i = 1, size(taus)
         dir = xyzh(1:3,i)-primary
         call ray_polation(nsides, dir, listsOfTaus, listsOfDists, listOfLens, taus(i))
      enddo
      !$omp end parallel do
   end subroutine get_all_tau_outwards_single

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the uniform outwards
   !  ray-tracing scheme concerning a binary system
   !+
   !  IN: primary:         The location of the primary star
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: opacities:       The list of the opacities of the particles
   !  IN: Rstar:           The radius of the star
   !  IN: order:           The order in which the rays are sampled
   !  IN: companion:       The location of the companion
   !  IN: R:               The radius of the companion
   !+
   !  OUT: taus:           The list of optical depths for each particle
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, order, taus, companion, R)
      integer, intent(in) :: order
      real, intent(in)    :: primary(3), companion(3), opacities(:), Rstar, xyzh(:,:), R
      real, intent(out)   :: taus(:)
      
      integer  :: i, nrays, nsides, len
      real     :: normCompanion = 0., theta0 = 0., unitCompanion(3) = 0.
      real     :: theta, root, dist, dir(3), phi = 0., cosphi = 0., sinphi = 0.
      real, dimension(:,:), allocatable  :: dirs
      real, dimension(:,:), allocatable  :: listsOfDists, listsOfTaus
      real, dimension(:), allocatable    :: tau, dists
      integer, dimension(:), allocatable :: listOfLens
      nrays = 12*4**order
      nsides = 2**order
      
      allocate(dirs(3, nrays))
      allocate(listsOfDists(200, nrays))
      allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
      allocate(listOfLens(nrays))
      allocate(tau(size(listsOfDists(:,1))))
      allocate(dists(size(listsOfDists(:,1))))

      unitCompanion = companion-primary
      normCompanion = norm2(unitCompanion)
      theta0 = asin(R/normCompanion)
      unitCompanion = unitCompanion/normCompanion
      phi = atan2(unitCompanion(2),unitCompanion(1))
      cosphi = cos(phi)
      sinphi = sin(phi)
     
      !$omp parallel do private(dir,tau,dists,dist,root,theta,len)
      do i = 1, nrays
         tau=0.
         dists=0.
         call pix2vec_nest(nsides, i-1, dir)
         dirs(:,i) = dir
         dir = (/cosphi*dir(1) - sinphi*dir(2),sinphi*dir(1) + cosphi*dir(2), dir(3)/)
         theta = acos(dot_product(unitCompanion, dir))
         if (theta < theta0) then
            root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
            dist = normCompanion*cos(theta)-root
            call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, len, dist)
         else
            call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, len)
         endif
         listsOfTaus(:,i) = tau
         listsOfDists(:,i) = dists
         listOfLens(i) = len
      enddo
      !$omp end parallel do

      taus = 0.
      !$omp parallel do private(dir)
      do i = 1, size(taus)
         dir = xyzh(1:3,i)-primary
         dir = (/cosphi*dir(1) + sinphi*dir(2),-sinphi*dir(1) + cosphi*dir(2), dir(3)/)
         call ray_polation(nsides, dir, listsOfTaus, listsOfDists, listOfLens, taus(i))
      enddo
      !$omp end parallel do
   end subroutine get_all_tau_outwards_companion

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of a particle with a distance on the ray
   !+
   !  IN: dist:            The distance of a point on the ray
   !  IN: listOfTaus:      The distribution of optical depths throughout the ray
   !  IN: listOfDists:     The distribution of distances throughout the ray
   !+
   !  OUT: tau:            The list of optical depths for the particle
   !+
   !--------------------------------------------------------------------------
   subroutine ray_polation(nsides, vec, listsOfTaus, listsOfDists, listOfLens, tau)
      integer, intent(in) :: nsides, listOfLens(:)
      real, intent(in)    :: vec(:), listsOfDists(:,:), listsOfTaus(:,:)
      real, intent(out)   :: tau
   
      integer :: index, n(8), nneigh, i
      real    :: tautemp, vectemp(3), suminvdist2, tempdist(8), tempdistIndex
      logical :: mk(8)
      

      ! call vec2pix_nest(nsides, vec, index)
      ! index = index + 1
      ! call get_tau_outwards(norm2(vec), listsOfTaus(:,index), listsOfDists(:,index), listOfLens(index), tau)
      
      call vec2pix_nest(nsides, vec, index)
      call pix2vec_nest(nsides, index, vectemp)
      vectemp = vec - norm2(vec)*vectemp
      tempdistIndex = dot_product(vectemp,vectemp)
      call get_tau_outwards(norm2(vec), listsOfTaus(:,index+1), listsOfDists(:,index+1), listOfLens(index+1), tautemp)
      tau = tautemp/tempdistIndex
      suminvdist2 = 1./tempdistIndex

      call neighbours_nest(nsides, index, n, nneigh)
      do i=1,nneigh
         call pix2vec_nest(nsides, n(i), vectemp)
         vectemp = vec - norm2(vec)*vectemp
         tempdist(i) = dot_product(vectemp,vectemp)
      enddo
      n=n+1
      mk = .true.
      mk(nneigh+1:8) = .false.
      do i=1,3
         index = minloc(tempdist,1,mk)
         mk(index) = .false.
         call get_tau_outwards(norm2(vec), listsOfTaus(:,n(index)), listsOfDists(:,n(index)), listOfLens(n(index)), tautemp)
         tau = tau + tautemp/tempdist(index)
         suminvdist2 = suminvdist2 + 1./tempdist(index)
      enddo
      tau = tau / suminvdist2
   end subroutine ray_polation
   
    
   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of a particle with a distance on the ray
   !+
   !  IN: dist:            The distance of a point on the ray
   !  IN: listOfTaus:      The distribution of optical depths throughout the ray
   !  IN: listOfDists:     The distribution of distances throughout the ray
   !+
   !  OUT: tau:            The list of optical depths for the particle
   !+
   !--------------------------------------------------------------------------
   subroutine get_tau_outwards(dist, listOfTaus, listOfDists, len, tau)
      real, intent(in)    :: dist, listOfTaus(:), listOfDists(:)
      integer, intent(in) :: len
      real, intent(out)   :: tau
   
      integer :: L, R, m
   
      if (dist .lt. listOfDists(1)) then
         tau = 0.
      else if (dist .gt. listOfDists(len)) then
         tau = 1e10
      else
         L = 2
         R = len-1
         do while (L < R)
            m = (L + R)/2
            if (listOfDists(m) > dist) then
                  R = m
            else
                  L = m + 1
            end if
         enddo
         tau = listOfTaus(L-1)+(listOfTaus(L)-listOfTaus(L-1))/(listOfDists(L)-listOfDists(L-1))*(dist-listOfDists(L-1))
      endif
   end subroutine get_tau_outwards
    
   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth along a given ray
   !+
   !  IN: primary:         The location of the primary star
   !  IN: ray:             The direction of the ray that needs to be traced
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: opacities:       The list of the opacities of the particles
   !  IN: Rstar:           The radius of the star
   !+
   !  OUT: taus:           The distribution of optical depths throughout the ray
   !  OUT: listOfDists:    The distribution of distances throughout the ray
   !+
   !  OPT: maxDist:        The maximal distance the ray needs to be traced
   !+
   !--------------------------------------------------------------------------
   subroutine ray_tracer(primary, ray, xyzh, opacities, Rstar, taus, listOfDist, len, maxDist)
      use linklist, only:getneigh_pos,ifirstincell,listneigh
      use kernel,   only:radkern
      real, intent(in)     :: primary(3), ray(3), Rstar, xyzh(:,:), opacities(:)
      real, optional       :: maxDist
      real, intent(out)    :: listOfDist(:), taus(:)
      integer, intent(out) :: len
      
      integer, parameter :: maxcache = 0
      real, allocatable  :: xyzcache(:,:)
      real :: dist, h, opacity, previousOpacity, nextOpacity
      integer :: nneigh, next, i

      dist = Rstar
      listOfDist(1)=dist
      h = Rstar/100.

      next=0
      do while (next==0)
         h = h*2.
         call getneigh_pos(primary+Rstar*ray,0.,h,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
         call find_next(primary, ray, dist, xyzh, listneigh, next, nneigh)
      enddo
      call calc_opacity(primary+Rstar*ray, xyzh, opacities, listneigh, nneigh, previousOpacity)

      i = 1
      do while (hasNext(next,dist,maxDist))
         i = i + 1
         call getneigh_pos(primary + dist*ray,0.,xyzh(4,next)*radkern, &
                           3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
         call calc_opacity(primary + dist*ray, xyzh, opacities, listneigh, nneigh, nextOpacity)
         opacity = (nextOpacity+previousOpacity)/2
         previousOpacity=nextOpacity
         taus(i) = taus(i-1)+(dist-listOfDist(i-1))*opacity
         listOfDist(i)=dist
         call find_next(primary, ray, dist, xyzh, listneigh, next,nneigh)
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
   !  Find the next point on a ray
   !+
   !  IN: inpoint:         The initial point for which the next point will 
   !                       be calculated
   !  IN: ray:             The unit vector of the direction in which the next 
   !                       point will be calculated
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: neighbors:       A list containing the indices of the neighbors of 
   !                       the initial point
   !  IN: next:            The index of the initial point
   !                       (this point will not be considered as possible next point)
   !  IN: nneighin:        The amount of neighbors
   !+
   !  OUT: next:           The next point on the ray
   !+
   !--------------------------------------------------------------------------
   subroutine find_next(inpoint, ray, dist, xyzh, neighbors, next, nneighin)
      integer, intent(in)    :: neighbors(:)
      real, intent(in)       :: xyzh(:,:), inpoint(:), ray(:)
      integer, intent(inout) :: next
      real, intent(inout)    :: dist
      integer, optional      :: nneighin
      
      real     :: trace_point(3), dmin, vec(3), tempdist, raydist
      real     :: nextdist
      integer  :: i, nneigh, prev
      
      dmin = huge(0.)
      if (present(nneighin)) then
         nneigh = nneighin
      else
         nneigh = size(neighbors)
      endif

      prev=next
      next=0
      nextDist=dist
      trace_point = inpoint + dist*ray

      i = 1
      do while (i <= nneigh .and. neighbors(i) /= 0)
         if (neighbors(i) .ne. prev) then
            vec=xyzh(1:3,neighbors(i)) - trace_point
            tempdist = dot_product(vec,ray)
            if (tempdist>0.) then
               raydist = dot_product(vec,vec) - tempdist**2
               if (raydist < dmin) then
                  dmin = raydist
                  next = neighbors(i)
                  nextdist = dist+tempdist
               end if
            end if
         endif
         i = i+1
      enddo
      dist=nextdist
   end subroutine find_next

   !--------------------------------------------------------------------------
   !+
   !  Calculate the opacity in a given location
   !+
   !  IN: r0:              The location where the opacity will be calculated
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: opacities:       The list of the opacities of the particles
   !  IN: neighbors:       A list containing the indices of the neighbors of 
   !                       the initial point
   !  IN: nneigh:          The amount of neighbors
   !+
   !  OUT: opacity:        The opacity at the given location
   !+
   !--------------------------------------------------------------------------
   subroutine calc_opacity(r0, xyzh, opacities, neighbors, nneigh, opacity)
      use kernel,   only:cnormk,wkern
      use part,     only:hfact
      real, intent(in)    :: r0(:), xyzh(:,:), opacities(:)
      integer, intent(in) :: neighbors(:), nneigh
      real, intent(out)   :: opacity

      integer :: i
      real    :: fact, q

      fact = cnormk/hfact**3
      opacity=0
      do i=1,nneigh
         q = norm2(r0 - xyzh(1:3,neighbors(i)))/xyzh(4,neighbors(i))
         opacity=opacity+fact*wkern(q*q,q)*opacities(neighbors(i))
      enddo
   end subroutine calc_opacity
end module raytracer