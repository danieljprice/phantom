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
         call get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, companion, R, order, taus)
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
   subroutine get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, companion, R, order, taus)
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
   !  Calculate the optical depth of a particle
   !+
   !  IN: nsides:          the healpix nsides of the simulation
   !  IN: vec:             the vector from the primary to the point
   !  IN: listOfTaus:      The distribution of optical depths throughout the ray
   !  IN: listOfDists:     The distribution of distances throughout the ray
   !  IN: listOfLens:      The list of lengths of listOfTaus and listOfDists
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
   !  IN: len:             The length of listOfTaus and listOfDists
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
   !  Calculate the optical depts along a given ray
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
      real, intent(in)     :: primary(3), ray(3), Rstar, xyzh(:,:), opacities(:)
      real, optional       :: maxDist
      real, intent(out)    :: listOfDist(:), taus(:)
      integer, intent(out) :: len
      
      real    :: dist, nextDist, h, opacity, previousOpacity, nextOpacity, totalDist
      integer :: next, i

      h = Rstar/100.
      next=0
      do while (next==0)
         h = h*2.
         call find_next(primary+Rstar*ray, h, ray, xyzh, opacities, previousOpacity, dist, next)
      enddo

      i = 1
      totalDist = Rstar
      listOfDist(i)=TotalDist
      do while (hasNext(next,totalDist,maxDist))
         totalDist = totalDist+dist
         call find_next(primary + totalDist*ray, xyzh(4,next), ray, xyzh, opacities, nextOpacity, nextDist, next)
         i = i + 1
         opacity = (nextOpacity+previousOpacity)/2
         previousOpacity=nextOpacity
         taus(i) = taus(i-1)+(dist)*opacity
         listOfDist(i)=TotalDist
         dist=nextDist
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
   !  IN: inpoint:         The initial point for which the opacity and the 
   !                       next point will be calculated
   !  IN: h:               The smoothing point at the initial point
   !  IN: ray:             The unit vector of the direction in which the next 
   !                       point will be calculated
   !  IN: xyzh:            The xyzh of all the particles
   !  IN: opacities:       The list of the opacities of the particles
   !  IN: next:            The index of the initial point
   !                       (this point will not be considered as possible next point)
   !+
   !  OUT: opacity:        The opacity at the given location
   !  OUT: dist:           The distance to the next point
   !  OUT: next:           The next point on the ray
   !+
   !--------------------------------------------------------------------------
   subroutine find_next(inpoint, h, ray, xyzh, opacities, opacity, dist, next)
      use linklist, only:getneigh_pos,ifirstincell,listneigh
      use kernel,   only:radkern,cnormk,wkern
      use part,     only:hfact, rhoh, massoftype, igas
      real, intent(in)       :: xyzh(:,:), opacities(:), inpoint(:), ray(:), h
      integer, intent(inout) :: next
      real, intent(out)      :: dist, opacity

      integer, parameter :: maxcache = 0
      real, allocatable  :: xyzcache(:,:)
      
      integer  :: nneigh, i, prev
      real     :: dmin, vec(3), tempdist, raydist, q, norm2

      prev=next
      next=0
      dist=0.

      call getneigh_pos(inpoint,0.,h*radkern, 3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)

      opacity=0
      dmin = huge(0.)
      do i=1,nneigh
         vec=xyzh(1:3,listneigh(i)) - inpoint
         norm2 = dot_product(vec,vec)
         q = sqrt(norm2)/xyzh(4,listneigh(i))
         opacity=opacity+wkern(q*q,q)*opacities(listneigh(i))*rhoh(xyzh(4,listneigh(i)), massoftype(igas))

         if (listneigh(i) .ne. prev) then
            tempdist = dot_product(vec,ray)
            if (tempdist>0.) then
               raydist = norm2 - tempdist**2
               if (raydist < dmin) then
                  dmin = raydist
                  next = listneigh(i)
                  dist = tempdist
               end if
            end if
         endif
      enddo
      opacity = opacity*cnormk/hfact**3
   end subroutine find_next
end module raytracer