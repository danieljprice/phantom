module raytracer
   use healpix
   implicit none
   public :: get_all_tau_outwards, get_all_tau_inwards, get_all_tau_optimised
   private
   contains

#define REALTIME
#define SMOOTHING
   
   !*********************************************************************!
   !***************************   OPTIMISED   ***************************!
   !*********************************************************************!
  
#ifdef REALTIME
   subroutine get_all_tau_optimised(primary, xyzh, opacities, &
                                   Rstar, minOrder, refineLevel, taus, companion, R, maxDist)
      integer, intent(in) :: primary, minOrder, refineLevel
#else
   subroutine get_all_tau_optimised(primary, xyzh, neighbors, opacities, &
                                   Rstar, minOrder, refineLevel, taus, companion, R, maxDist)
      integer, intent(in) :: primary, neighbors(:,:), minOrder, refineLevel
#endif
      integer, optional   :: companion
      real, intent(in)    :: opacities(:), xyzh(:,:), Rstar
      real, optional      :: R, maxDist
      real, intent(out)   :: taus(:)
   
      integer     :: i, nrays, nsides, index
      real        :: normCompanion, theta0, unitCompanion(3), theta, root, dist, vec(3), dir(3), phi, cosphi, sinphi
      real, dimension(:,:), allocatable :: dirs
      real, dimension(:,:), allocatable :: listsOfDists, listsOfTaus
      integer, dimension(:), allocatable :: indices
      real, dimension(:), allocatable    :: tau, dists

      if (present(companion) .and. present(R)) then
         unitCompanion = xyzh(1:3,companion)-xyzh(1:3,primary)
         normCompanion = norm2(unitCompanion)
         theta0 = asin(R/normCompanion)
         unitCompanion = unitCompanion/normCompanion
         phi = atan(unitCompanion(2)/unitCompanion(1))
         cosphi = cos(phi)
         sinphi = sin(phi)

         nrays = 12*4**minOrder + 12*(refineLevel+minOrder) + int(24*2**(refineLevel+minOrder)*theta0)
         nsides = 2**(minOrder+refineLevel)
         allocate(dirs(3, nrays))
         allocate(indices(12*4**(minOrder+refineLevel)))
         allocate(listsOfDists(200, nrays))
         allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
         allocate(tau(size(listsOfDists(:,1))))
         allocate(dists(size(listsOfDists(:,1))))

         call get_rays(normCompanion,R,minOrder,refineLevel,dirs, indices, nrays)
         !$omp parallel do private(tau,dist,dir,dists,root,theta)
         do i = 1, nrays
            tau=0.
            dists=0.
            dir = dirs(:,i)
            dir = (/cosphi*dir(1) - sinphi*dir(2),sinphi*dir(1) + cosphi*dir(2), dir(3)/)
            theta = acos(dot_product(unitCompanion, dir))
            if (theta < theta0) then
               root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
               dist = normCompanion*cos(theta)-root
#ifdef REALTIME
               call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, dist)
#else
               call ray_tracer(primary, dir, xyzh, opacities, neighbors, Rstar, tau, dists, dist)
#endif
            else
#ifdef REALTIME
               call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, maxDist)
#else
               call ray_tracer(primary, dir, xyzh, opacities, neighbors, Rstar, tau, dists, maxDist)
#endif
            endif
            listsOfTaus(:,i) = tau
            listsOfDists(:,i) = dists
         enddo
         !$omp end parallel do
         
         taus = 0.
         !$omp parallel do private(index,vec)
         do i = 1, size(taus(:))
            vec = xyzh(1:3,i)-xyzh(1:3,primary)
            vec = (/cosphi*vec(1) + sinphi*vec(2),-sinphi*vec(1) + cosphi*vec(2), vec(3)/)
            call vec2pix_nest(nsides, vec, index)
            index = indices(index + 1)
            call get_tau_outwards(dot_product(vec, dirs(:,index)), listsOfTaus(:,index), listsOfDists(:,index), taus(i))
         enddo
         !$omp end parallel do

      else
#ifdef REALTIME
         call get_all_tau_outwards_single(primary, xyzh, opacities, &
         Rstar, minOrder+refineLevel, taus, maxDist)
#else
         call get_all_tau_outwards_single(primary, xyzh, neighbors, opacities, &
         Rstar, minOrder+refineLevel, taus, maxDist)
#endif
      endif
   end subroutine get_all_tau_optimised

   subroutine get_rays(d,R,minOrder,refineLevel, vecs, indices, nrays)
      real, intent(in)     :: d, R
      integer, intent(in)  :: minOrder, refineLevel
      real, intent(out)    :: vecs(:,:)
      integer, intent(out) :: indices(:),nrays

      real    :: theta, dist
      real, dimension(:,:), allocatable  :: circ
      integer :: i, j, k, minNsides, minNrays, ind, ind2,n
      integer, dimension(:), allocatable  :: refine,refine2

      minNsides = 2**minOrder
      minNrays = 12*4**minOrder
      if (refineLevel == 0) then
         do i=1, minNrays
            call pix2vec_nest(minNsides,i-1, vecs(:,i))
            indices(i) = i
         enddo
         nrays = minNrays
         return
      endif
      theta = asin(R/d)
      dist = d*cos(theta)
      n = int(theta*6*2**(minOrder+refineLevel))+4
      allocate(circ(n,3))
      allocate(refine(n))
      allocate(refine2(n*2))
      refine2 = 0
      do i=1, n
         circ(i,1) = dist*cos(theta)
         circ(i,2) = dist*sin(theta)*cos(2*PI*i/n)
         circ(i,3) = dist*sin(theta)*sin(2*PI*i/n)
      enddo

      do i=1, n
         call vec2pix_nest(minNsides,circ(i,:),refine(i))
      enddo
      
      ind = 1
      ind2 = 1
      do i = 1, minNrays
         if (any(i-1 == refine)) then
            refine2(ind2) = i-1
            ind2 = ind2+1
         else
            call pix2vec_nest(minNsides,i-1,vecs(:,ind))
            indices(4**(refineLevel)*(i-1)+1:4**(refineLevel)*(i)) = ind
            ind = ind+1
         endif
      enddo
      
      i=1
      do k=1,refineLevel-1
         do j=1, n
            call vec2pix_nest(2**(minOrder+k),circ(j,:),refine(j))
         enddo
         do i = i, ind2-1
            do j = 1,4
               if (any(4*refine2(mod(i-1,n*2)+1)+j-1 == refine)) then
                  refine2(mod(ind2-1,n*2)+1) = 4*refine2(mod(i-1,n*2)+1)+j-1
                  ind2 = ind2+1
               else
                  call pix2vec_nest(2**(minOrder+k),4*refine2(mod(i-1,n*2)+1)+j-1,vecs(:,ind))
                  indices(4**(refineLevel-k)*(4*refine2(mod(i-1,n*2)+1)+j-1)+1: &
                        4**(refineLevel-k)*(4*(refine2(mod(i-1,n*2)+1))+j)) = ind
                  ind = ind+1
               endif
            enddo
         enddo
      enddo
      
      do i = i, ind2-1
         do j = 1,4
            call pix2vec_nest(2**(minOrder+refineLevel),4*refine2(mod(i-1,n*2)+1)+j-1,vecs(:,ind))
            indices(4*(refine2(mod(i-1,n*2)+1))+j) = ind
            ind = ind+1
         enddo
      enddo
      nrays = ind-1
   end subroutine get_rays

   !*********************************************************************!
   !***************************   OUTWARDS   ****************************!
   !*********************************************************************!

#ifdef REALTIME
   subroutine get_all_tau_outwards(primary, xyzh, opacities, Rstar, order, taus, companion, R, maxDist)
      integer, intent(in) :: primary, order
#else
   subroutine get_all_tau_outwards(primary, xyzh, neighbors, opacities, Rstar, order, taus, companion, R, maxDist)
      integer, intent(in) :: primary, neighbors(:,:), order
#endif
      integer, optional   :: companion
      real, intent(in)    :: opacities(:), Rstar, xyzh(:,:)
      real, optional      :: R, maxDist
      real, intent(out)   :: taus(:)

      if (present(companion) .and. present(R)) then
#ifdef REALTIME
         call get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, order, taus, companion, R, maxDist)
#else
         call get_all_tau_outwards_companion(primary, xyzh, neighbors, opacities, Rstar, order, taus, companion, R, maxDist)
#endif
      else
#ifdef REALTIME
         call get_all_tau_outwards_single(primary, xyzh, opacities, Rstar, order, taus, maxDist)
#else
         call get_all_tau_outwards_single(primary, xyzh, neighbors, opacities, Rstar, order, taus, maxDist)
#endif
      endif
   end subroutine get_all_tau_outwards
   
#ifdef REALTIME
   subroutine get_all_tau_outwards_single(primary, xyzh, opacities, Rstar, order, taus, maxDist)
      integer, intent(in) :: primary, order
#else
   subroutine get_all_tau_outwards_single(primary, xyzh, neighbors, opacities, Rstar, order, taus, maxDist)
      integer, intent(in) :: primary, neighbors(:,:), order
#endif
      real, intent(in)    :: opacities(:), Rstar, xyzh(:,:)
      real, optional      :: maxDist
      real, intent(out)   :: taus(:)
      
      integer  :: i, nrays, nsides, index
      real     :: vec(3)
      real     :: dist, dir(3)
      real, dimension(:,:), allocatable  :: dirs, listsOfDists, listsOfTaus
      real, dimension(:), allocatable    :: tau, dists
      nrays = 12*4**order
      nsides = 2**order
      
      allocate(dirs(3, nrays))
      allocate(listsOfDists(200, nrays))
      allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
      allocate(tau(size(listsOfDists(:,1))))
      allocate(dists(size(listsOfDists(:,1))))
     
      !$omp parallel do private(dir,tau,dists,dist)
      do i = 1, nrays
         tau=0.
         dists=0.
         call pix2vec_nest(nsides, i-1, dir)
         dirs(:,i) = dir
#ifdef REALTIME
         call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, maxDist)
#else
         call ray_tracer(primary, dir, xyzh, opacities, neighbors, Rstar, tau, dists, maxDist)
#endif
         listsOfTaus(:,i) = tau
         listsOfDists(:,i) = dists
      enddo
      !$omp end parallel do

      taus = 0.
      !$omp parallel do private(index,vec)
      do i = 1, size(taus)
         vec = xyzh(1:3,i)-xyzh(1:3,primary)
         call vec2pix_nest(nsides, vec, index)
         index = index + 1
         call get_tau_outwards(dot_product(vec, dirs(:,index)), listsOfTaus(:,index), listsOfDists(:,index), taus(i))
      enddo
      !$omp end parallel do
   end subroutine get_all_tau_outwards_single

#ifdef REALTIME
   subroutine get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, order, taus, companion, R, maxDist)
      integer, intent(in) :: primary, order, companion
#else
   subroutine get_all_tau_outwards_companion(primary, xyzh, neighbors, opacities, Rstar, order, taus, companion, R, maxDist)
      integer, intent(in) :: primary, neighbors(:,:), order, companion
#endif
      real, intent(in)    :: opacities(:), Rstar, xyzh(:,:), R, maxDist
      real, intent(out)   :: taus(:)
      
      integer  :: i, nrays, nsides, index
      real     :: normCompanion = 0., theta0 = 0., unitCompanion(3) = 0., vec(3)
      real     :: theta, root, dist, dir(3), phi = 0., cosphi = 0., sinphi = 0.
      real, dimension(:,:), allocatable  :: dirs
      real, dimension(:,:), allocatable  :: listsOfDists, listsOfTaus
      real, dimension(:), allocatable    :: tau, dists
      nrays = 12*4**order
      nsides = 2**order
      
      allocate(dirs(3, nrays))
      allocate(listsOfDists(200, nrays))
      allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
      allocate(tau(size(listsOfDists(:,1))))
      allocate(dists(size(listsOfDists(:,1))))

      unitCompanion = xyzh(1:3,companion)-xyzh(1:3,primary)
      normCompanion = norm2(unitCompanion)
      theta0 = asin(R/normCompanion)
      unitCompanion = unitCompanion/normCompanion
      phi = atan(unitCompanion(2)/unitCompanion(1))
      cosphi = cos(phi)
      sinphi = sin(phi)
     
      !$omp parallel do private(dir,tau,dists,dist,root,theta)
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
#ifdef REALTIME
            call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, dist)
#else
            call ray_tracer(primary, dir, xyzh, opacities, neighbors, Rstar, tau, dists, dist)
#endif
         else
#ifdef REALTIME
            call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, maxDist)
#else
            call ray_tracer(primary, dir, xyzh, opacities, neighbors, Rstar, tau, dists, maxDist)
#endif
         endif
         listsOfTaus(:,i) = tau
         listsOfDists(:,i) = dists
      enddo
      !$omp end parallel do

      taus = 0.
      !$omp parallel do private(index,vec,dir)
      do i = 1, size(taus)
         vec = xyzh(1:3,i)-xyzh(1:3,primary)
         vec = (/cosphi*vec(1) + sinphi*vec(2),-sinphi*vec(1) + cosphi*vec(2), vec(3)/)
         call vec2pix_nest(nsides, vec, index)
         index = index + 1
         call get_tau_outwards(dot_product(vec, dirs(:,index)), listsOfTaus(:,index), listsOfDists(:,index), taus(i))
      enddo
      !$omp end parallel do
   end subroutine get_all_tau_outwards_companion
    
   subroutine get_tau_outwards(dist, listOfTaus, listOfDists, tau)
      real, intent(in)    :: dist, listOfTaus(:), listOfDists(:)
      real, intent(out)   :: tau
   
      integer :: i
   
      if (dist .gt. listOfDists(1)) then
         i = 2
         do while (dist .gt. listOfDists(i) .and. listOfDists(i) .gt. 0.)
            i = i+1
         enddo
         if (dist .gt. listOfDists(i)) then
            tau = 1e10
         else
            tau = listOfTaus(i-1)+(listOfTaus(i)-listOfTaus(i-1))/(listOfDists(i)-listOfDists(i-1))*(dist-listOfDists(i-1))
         endif
      else
         tau = 0.
      endif
   end subroutine get_tau_outwards
    
#ifdef REALTIME
   subroutine ray_tracer(point, ray, xyzh, opacities, Rstar, taus, listOfDist, maxDist)
#else
   subroutine ray_tracer(point, ray, xyzh, opacities, neighbors, Rstar, taus, listOfDist, maxDist)
#endif
      use linklist, only:getneigh_pos,ifirstincell,listneigh
#ifdef REALTIME
      use kernel,   only:radkern
#endif
      real, intent(in)     :: ray(3), Rstar, xyzh(:,:), opacities(:)
      real, optional       :: maxDist
#ifdef REALTIME
      integer, intent(in)  :: point
#else
      integer, intent(in)  :: point, neighbors(:,:)
#endif
      real, intent(out)    :: listOfDist(:), taus(:)
      
      integer, parameter :: maxcache = 0
      real, allocatable  :: xyzcache(:,:)
#ifdef SMOOTHING
      real :: dist, h, opacity, previousOpacity, nextOpacity
      integer :: nneigh, next, i
#else
      real :: dist, h, opacity
      integer :: nneigh, previous, next, i

      previous = point
#endif
      dist = Rstar
      listOfDist(1)=dist
      h = Rstar/100.

      next=0
      do while (next==0)
         h = h*2.
         call getneigh_pos(xyzh(1:3,point)+Rstar*ray,0.,h,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
         call find_next(xyzh(1:3,point), ray, dist, xyzh, listneigh, next, nneigh)
      enddo
#ifdef SMOOTHING
      call calc_opacity(xyzh(1:3,point)+Rstar*ray, xyzh, opacities, listneigh, nneigh, previousOpacity)
#endif

      i = 1
      do while (hasNext(next,dist,maxDist))
         i = i + 1
#ifdef REALTIME
         call getneigh_pos(xyzh(1:3,point) + dist*ray,0.,xyzh(4,next)*radkern, &
                           3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
#else
         listneigh = neighbors(next,:)
         nneigh = nneigh+1
         listneigh(nneigh) = next
#endif
#ifdef SMOOTHING
         call calc_opacity(xyzh(1:3,point) + dist*ray, xyzh, opacities, listneigh, nneigh, nextOpacity)
         opacity = (nextOpacity+previousOpacity)/2
         previousOpacity=nextOpacity
#else
         opacity = (opacities(next) + opacities(previous))/2
         previous = next
#endif
         taus(i) = taus(i-1)+(dist-listOfDist(i-1))*opacity
         listOfDist(i)=dist
#ifdef REALTIME
         call find_next(xyzh(1:3,point), ray, dist, xyzh, listneigh, next,nneigh)
#else
         call find_next(xyzh(1:3,point), ray, dist, xyzh, neighbors(next,:), next)
#endif
      enddo
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
   !****************************   INWARDS   ****************************!
   !*********************************************************************!
   
   subroutine get_all_tau_inwards(primary, xyzh, neighbors, opacities, Rstar, taus, companion, R)
      real, intent(in)    :: opacities(:), Rstar, xyzh(:,:)
      real, optional      :: R
      integer, intent(in) :: primary, neighbors(:,:)
      integer, optional   :: companion
      real, intent(out)   :: taus(:)
      
      if (present(companion) .and. present(R)) then
         call get_all_tau_inwards_companion(primary, xyzh, neighbors, opacities, Rstar, taus, companion, R)
      else
         call get_all_tau_inwards_single(primary, xyzh, neighbors, opacities, Rstar, taus)
      endif
   end subroutine get_all_tau_inwards

   subroutine get_all_tau_inwards_single(primary, xyzh, neighbors, opacities, Rstar, taus)
      real, intent(in)    :: opacities(:), Rstar, xyzh(:,:)
      integer, intent(in) :: primary, neighbors(:,:)
      real, intent(out)   :: taus(:)
      
      integer :: i
      real    :: tau
      
      !$omp parallel do private(tau)
      do i = 1, size(taus)
         call get_tau_inwards(i, primary, xyzh, neighbors, opacities, Rstar, tau)
         taus(i) = tau
      enddo
      !$omp end parallel do
   end subroutine get_all_tau_inwards_single

   subroutine get_all_tau_inwards_companion(primary, xyzh, neighbors, opacities, Rstar, taus, companion, R)
      real, intent(in)    :: opacities(:), Rstar, xyzh(:,:), R
      integer, intent(in) :: primary, neighbors(:,:), companion
      real, intent(out)   :: taus(:)
      
      integer :: i
      real    :: normCompanion, theta0, unitCompanion(3), norm, theta, root, norm0, tau
      
      normCompanion = 0.
      theta0 = 0.
      unitCompanion = 0.
      normCompanion = norm2(xyzh(1:3,companion)-xyzh(1:3,primary))
      theta0 = asin(R/normCompanion)
      unitCompanion = (xyzh(1:3,companion)-xyzh(1:3,primary))/normCompanion
      
      !$omp parallel do private(tau,norm,theta,root,norm0)
      do i = 1, size(taus)
         norm = norm2(xyzh(1:3,i)-xyzh(1:3,primary))
         theta = acos(dot_product(unitCompanion, xyzh(1:3,i)-xyzh(1:3,primary))/norm)
         if (theta < theta0) then
            root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
            norm0 = normCompanion*cos(theta)-root
            if (norm > norm0) then
                  tau = 1e10
            else
               call get_tau_inwards(i, primary, xyzh, neighbors, opacities, Rstar, tau)
            endif
         else
            call get_tau_inwards(i, primary, xyzh, neighbors, opacities, Rstar, tau)
         endif
         taus(i) = tau
      enddo
      !$omp end parallel do
   end subroutine get_all_tau_inwards_companion
    
   subroutine get_tau_inwards(secondary, primary, xyzh, neighbors, opacities, Rstar, tau)
#ifdef SMOOTHING
      use linklist, only:getneigh_pos,ifirstincell,listneigh
      use kernel,   only:radkern
#endif
      real, intent(in)    :: xyzh(:,:), opacities(:), Rstar
      integer, intent(in) :: primary, secondary, neighbors(:,:)
      real, intent(out)   :: tau

#ifdef SMOOTHING
      integer :: i, next, previous, nneigh
      integer, parameter :: maxcache = 0
      real, allocatable  :: xyzcache(:,:)
      real    :: ray(3), nextDist, previousDist, maxDist, opacity, previousOpacity, nextOpacity
#else
      integer :: i, next, previous
      real    :: ray(3), nextDist, previousDist, maxDist, opacity
#endif

      ray = xyzh(1:3,primary) - xyzh(1:3,secondary)
      maxDist = norm2(ray)
      ray = ray / maxDist
      maxDist=max(maxDist-Rstar,0.)
      next=secondary
      nextOpacity = opacities(next)
      nextDist=0.
      
      tau = 0.
      i=1
      do while (next /= primary .and. next /=0)
         i = i + 1
         previous = next
         previousDist = nextDist
         call find_next(xyzh(1:3,secondary), ray, nextDist, xyzh, neighbors(next,:), next)
         if (nextDist .gt. maxDist) then
               next = primary
               nextDist = maxDist
#ifdef SMOOTHING
         endif
         call getneigh_pos(xyzh(1:3,secondary) + nextDist*ray,0.,xyzh(4,previous)*radkern, &
                           3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
         previousOpacity=nextOpacity
         call calc_opacity(xyzh(1:3,secondary) + nextDist*ray, xyzh, opacities, listneigh, nneigh, nextOpacity)
         opacity = (nextOpacity+previousOpacity)/2
#else
            opacity = opacities(previous)
         else
         opacity = (opacities(next)+opacities(previous))/2
         endif
#endif
         tau = tau + (nextDist-previousDist)*opacity
      enddo
   end subroutine get_tau_inwards

   !*********************************************************************!
   !****************************   COMMON   *****************************!
   !*********************************************************************!
   
   subroutine find_next(inpoint, ray, dist, xyzh, neighbors, next, nneighin)
      integer, intent(in)  :: neighbors(:)
      real, intent(in)     :: xyzh(:,:), inpoint(:), ray(:)
      integer, intent(inout) :: next
      real, intent(inout)  :: dist
      integer, optional    :: nneighin
      
      real                 :: trace_point(3), dmin, vec(3), tempdist, raydist
      real                 :: nextdist
      integer              :: i, nneigh, prev
      
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

#ifdef SMOOTHING
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
#endif
end module raytracer