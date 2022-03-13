module raytracer
   use healpix
   implicit none
   public :: get_all_tau_outwards, get_all_tau_inwards, get_all_tau_adaptive
   private
   contains
   
   !*********************************************************************!
   !***************************   ADAPTIVE   ****************************!
   !*********************************************************************!
  

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the adaptive ray-
   !  tracing scheme
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_adaptive(primary, xyzh, opacities, &
                                   Rstar, minOrder, refineLevel, taus, companion, R, maxDist)
      integer, intent(in) :: primary, minOrder, refineLevel
      integer, optional   :: companion
      real, intent(in)    :: opacities(:), xyzh(:,:), Rstar
      real, optional      :: R, maxDist
      real, intent(out)   :: taus(:)
   
      integer     :: i, nrays, nsides, index
      real        :: normCompanion, theta0, unitCompanion(3), theta, root, dist, vec(3), dir(3)
      real, dimension(:,:), allocatable :: dirs
      real, dimension(:,:), allocatable :: listsOfDists, listsOfTaus
      integer, dimension(:), allocatable :: indices
      real, dimension(:), allocatable    :: tau, dists

      if (present(companion) .and. present(R)) then
         unitCompanion = xyzh(1:3,companion)-xyzh(1:3,primary)
         normCompanion = norm2(unitCompanion)
         theta0 = asin(R/normCompanion)
         unitCompanion = unitCompanion/normCompanion

         nrays = int(12*4**minOrder*(-1+2**(refineLevel) + 2**(refineLevel+1))/2)
         nsides = 2**(minOrder+refineLevel)
         allocate(dirs(3, nrays))
         allocate(indices(12*4**(minOrder+refineLevel)))
         allocate(listsOfDists(200, nrays))
         allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
         allocate(tau(size(listsOfDists(:,1))))
         allocate(dists(size(listsOfDists(:,1))))

         call get_rays(xyzh(1:3,primary), xyzh(1:3, companion), xyzh, size(taus), minOrder, refineLevel, R, dirs, indices)
         !$omp parallel do private(tau,dist,dir,dists,root,theta)
         do i = 1, nrays
            tau=0.
            dists=0.
            dir = dirs(:,i)
            theta = acos(dot_product(unitCompanion, dir))
            if (theta < theta0) then
               root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
               dist = normCompanion*cos(theta)-root
               call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, dist)
            else
               call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, maxDist)
            endif
            listsOfTaus(:,i) = tau
            listsOfDists(:,i) = dists
         enddo
         !$omp end parallel do
         
         taus = 0.
         !$omp parallel do private(index,vec)
         do i = 1, size(taus(:))
            vec = xyzh(1:3,i)-xyzh(1:3,primary)
            call vec2pix_nest(nsides, vec, index)
            index = indices(index + 1)
            call get_tau_outwards(dot_product(vec, dirs(:,index)), listsOfTaus(:,index), listsOfDists(:,index), taus(i))
         enddo
         !$omp end parallel do

      else
         call get_all_tau_outwards_single(primary, xyzh, opacities, &
         Rstar, minOrder+refineLevel, taus, maxDist)
      endif
   end subroutine get_all_tau_adaptive

   !--------------------------------------------------------------------------
   !+
   !  Return all the directions of the rays that need to be traced for the
   !  adaptive ray-tracing scheme
   !+
   !--------------------------------------------------------------------------
   subroutine get_rays(primary, secondary, xyzh, npart, minOrder, refineLevel, R, vecs, indices)
      real, intent(in)     :: primary(:), secondary(:), xyzh(:,:), R
      integer, intent(in)  :: minOrder, refineLevel, npart
      real, intent(out)    :: vecs(:,:)
      integer, intent(out) :: indices(:)

      real    :: theta, dist, phi, cosphi, sinphi
      real, dimension(:,:), allocatable  :: circ
      integer :: i, j, minNsides, minNrays, ind,n, maxOrder, max, distr(12*4**(minOrder+refineLevel))
      integer, dimension(:,:), allocatable  :: distrs

      vecs = 0.
      indices = 0

      minNsides = 2**minOrder
      minNrays = 12*4**minOrder
      if (refineLevel == 0) then
         do i=1, minNrays
            call pix2vec_nest(minNsides,i-1, vecs(:,i))
            indices(i) = i
         enddo
      endif
      maxOrder = minOrder+refineLevel
      distr = 0
      !$omp parallel do private(ind)
      do i = 1, npart
         call vec2pix_nest(2**maxOrder, xyzh(1:3, i)-primary, ind)
         distr(ind+1) = distr(ind+1)+1
      enddo
      max = maxval(distr)

      dist = norm2(primary-secondary)
      theta = asin(R/dist)
      phi = atan2(secondary(2)-primary(2),secondary(1)-primary(1))
      phi = 0.76! TODO: FIX
      cosphi = cos(phi)
      sinphi = sin(phi)
      dist = dist*cos(theta)
      n = int(theta*6*2**(minOrder+refineLevel))+4
      allocate(circ(n,3))
      do i=1, n
         circ(i,1) = dist*cos(theta)
         circ(i,2) = dist*sin(theta)*cos(twopi*i/n)
         circ(i,3) = dist*sin(theta)*sin(twopi*i/n)
         circ(i,:) = (/cosphi*circ(i,1) - sinphi*circ(i,2),sinphi*circ(i,1) + cosphi*circ(i,2), circ(i,3)/)
      enddo

      do i=1, n
         call vec2pix_nest(minNsides,circ(i,:),ind)
         distr(ind) = max
      enddo

      allocate(distrs(12*4**(minOrder+refineLevel),refineLevel+1))
      distrs = 0
      distrs(:,1) = distr
      do i = 1, refineLevel
         do j = 1, 12*4**(maxOrder-i)
            distrs(j,i+1) = distrs(4*j,i)+distrs(4*j+1,i)+distrs(4*j+2,i)+distrs(4*j+3,i)
         enddo
      enddo
      max = maxval(distrs(:,refineLevel+1))+1

      ind=1
      do i=0, refineLevel-1
         call merge_argsort(distrs(1:12*4**(minOrder+i),refineLevel-i+1), distr)
         do j=1, 6*4**minOrder*2**(i)
            call pix2vec_nest(2**(minOrder+i), distr(j)-1, vecs(:,ind))
            indices(4**(refineLevel-i)*(distr(j)-1)+1:4**(refineLevel-i)*distr(j)) = ind
            ind=ind+1
            distrs(4*(distr(j)-1)+1:4*(distr(j)), refineLevel-i) = max
         enddo
         do j = j+1, 12*4**(minOrder+i)
            if (distrs(distr(j),refineLevel-i+1) == max) then
               distrs(4*(distr(j)-1)+1:4*(distr(j)), refineLevel-i) = max
            endif
         enddo
      enddo
      do i=1, 12*4**maxOrder
         if (distrs(i,1) .ne. max) then
            call pix2vec_nest(2**maxOrder, i-1, vecs(:,ind))
            indices(i) = ind
            ind=ind+1
         endif
      enddo
   end subroutine get_rays

   !--------------------------------------------------------------------------
   !+
   !  Routine that returns the arguments of the sorted array
   !  Source: https://github.com/Astrokiwi/simple_fortran_argsort/blob/master/sort_test.f90
   !+
   !--------------------------------------------------------------------------
   subroutine merge_argsort(r,d)
      integer, intent(in), dimension(:) :: r
      integer, intent(out), dimension(size(r)) :: d
    
      integer, dimension(size(r)) :: il

      integer :: stepsize
      integer :: i,j,n,left,k,ksize
    
      n = size(r)
    
      do i=1,n
          d(i)=i
      end do
    
      if ( n==1 ) return
    
      stepsize = 1
      do while (stepsize<n)
          do left=1,n-stepsize,stepsize*2
              i = left
              j = left+stepsize
              ksize = min(stepsize*2,n-left+1)
              k=1
        
              do while ( i<left+stepsize .and. j<left+ksize )
                  if ( r(d(i))<r(d(j)) ) then
                      il(k)=d(i)
                      i=i+1
                      k=k+1
                  else
                      il(k)=d(j)
                      j=j+1
                      k=k+1
                  endif
              enddo
        
              if ( i<left+stepsize ) then
                  ! fill up remaining from left
                  il(k:ksize) = d(i:left+stepsize-1)
              else
                  ! fill up remaining from right
                  il(k:ksize) = d(j:left+ksize-1)
              endif
              d(left:left+ksize-1) = il(1:ksize)
          end do
          stepsize=stepsize*2
      end do

      return
  end subroutine 

   !*********************************************************************!
   !***************************   OUTWARDS   ****************************!
   !*********************************************************************!

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the uniform outwards
   !  ray-tracing scheme
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards(primary, xyzh, opacities, Rstar, order, taus, companion, R, maxDist)
      integer, intent(in) :: primary, order
      integer, optional   :: companion
      real, intent(in)    :: opacities(:), Rstar, xyzh(:,:)
      real, optional      :: R, maxDist
      real, intent(out)   :: taus(:)

      if (present(companion) .and. present(R)) then
         call get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, order, taus, companion, R, maxDist)
      else
         call get_all_tau_outwards_single(primary, xyzh, opacities, Rstar, order, taus, maxDist)
      endif
   end subroutine get_all_tau_outwards
   
   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the uniform outwards
   !  ray-tracing scheme concerning only a single star
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards_single(primary, xyzh, opacities, Rstar, order, taus, maxDist)
      integer, intent(in) :: primary, order
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
         call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, maxDist)
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

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the uniform outwards
   !  ray-tracing scheme concerning a binary system
   !+
   !--------------------------------------------------------------------------
   subroutine get_all_tau_outwards_companion(primary, xyzh, opacities, Rstar, order, taus, companion, R, maxDist)
      integer, intent(in) :: primary, order, companion
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
      phi = atan2(unitCompanion(2),unitCompanion(1))
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
            call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, dist)
         else
            call ray_tracer(primary, dir, xyzh, opacities, Rstar, tau, dists, maxDist)
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
    
   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of a particle with a distance on the ray
   !+
   !--------------------------------------------------------------------------
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
    
   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth along a given ray
   !+
   !--------------------------------------------------------------------------
   subroutine ray_tracer(point, ray, xyzh, opacities, Rstar, taus, listOfDist, maxDist)
      use linklist, only:getneigh_pos,ifirstincell,listneigh
      use kernel,   only:radkern
      real, intent(in)     :: ray(3), Rstar, xyzh(:,:), opacities(:)
      real, optional       :: maxDist
      integer, intent(in)  :: point
      real, intent(out)    :: listOfDist(:), taus(:)
      
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
         call getneigh_pos(xyzh(1:3,point)+Rstar*ray,0.,h,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
         call find_next(xyzh(1:3,point), ray, dist, xyzh, listneigh, next, nneigh)
      enddo
      call calc_opacity(xyzh(1:3,point)+Rstar*ray, xyzh, opacities, listneigh, nneigh, previousOpacity)

      i = 1
      do while (hasNext(next,dist,maxDist))
         i = i + 1
         call getneigh_pos(xyzh(1:3,point) + dist*ray,0.,xyzh(4,next)*radkern, &
                           3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
         call calc_opacity(xyzh(1:3,point) + dist*ray, xyzh, opacities, listneigh, nneigh, nextOpacity)
         opacity = (nextOpacity+previousOpacity)/2
         previousOpacity=nextOpacity
         taus(i) = taus(i-1)+(dist-listOfDist(i-1))*opacity
         listOfDist(i)=dist
         call find_next(xyzh(1:3,point), ray, dist, xyzh, listneigh, next,nneigh)
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

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the inwards ray-
   !  tracing scheme
   !+
   !--------------------------------------------------------------------------
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

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the inwards ray-
   !  tracing scheme concerning only a single star
   !+
   !--------------------------------------------------------------------------
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

   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth of each particle, using the inwards ray-
   !  tracing scheme concerning a binary system
   !+
   !--------------------------------------------------------------------------
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
    
   !--------------------------------------------------------------------------
   !+
   !  Calculate the optical depth for a given particle, using the inwards ray-
   !  tracing scheme
   !+
   !--------------------------------------------------------------------------
   subroutine get_tau_inwards(secondary, primary, xyzh, neighbors, opacities, Rstar, tau)
      use linklist, only:getneigh_pos,ifirstincell,listneigh
      use kernel,   only:radkern
      real, intent(in)    :: xyzh(:,:), opacities(:), Rstar
      integer, intent(in) :: primary, secondary, neighbors(:,:)
      real, intent(out)   :: tau

      integer :: i, next, previous, nneigh
      integer, parameter :: maxcache = 0
      real, allocatable  :: xyzcache(:,:)
      real    :: ray(3), nextDist, previousDist, maxDist, opacity, previousOpacity, nextOpacity

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
         endif
         call getneigh_pos(xyzh(1:3,secondary) + nextDist*ray,0.,xyzh(4,previous)*radkern, &
                           3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
         previousOpacity=nextOpacity
         call calc_opacity(xyzh(1:3,secondary) + nextDist*ray, xyzh, opacities, listneigh, nneigh, nextOpacity)
         opacity = (nextOpacity+previousOpacity)/2
         tau = tau + (nextDist-previousDist)*opacity
      enddo
   end subroutine get_tau_inwards

   !*********************************************************************!
   !****************************   COMMON   *****************************!
   !*********************************************************************!
   
   !--------------------------------------------------------------------------
   !+
   !  Find the next point on a ray
   !+
   !--------------------------------------------------------------------------
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

   !--------------------------------------------------------------------------
   !+
   !  Calculate the opacity in a given location
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