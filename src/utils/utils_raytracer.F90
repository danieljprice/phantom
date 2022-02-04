module raytracer
   use healpix
   implicit none
   public :: get_all_tau_outwards, get_all_tau_inwards, get_all_tau_optimised
   private
   
!#define INLINE
!#define INLINE2
  contains
  
  !*********************************************************************!
  !***************************   OPTIMISED   ***************************!
  !*********************************************************************!
  
  subroutine get_all_tau_optimised(primary, points, xyzh, neighbors, opacities, &
                                   Rstar, minOrder, refineLevel, taus, companion, R, maxDist)
   integer, intent(in) :: primary, neighbors(:,:), minOrder, refineLevel
   integer, optional   :: companion
   real, intent(in)    :: points(:,:), opacities(:), xyzh(:,:), Rstar
   real, optional      :: R, maxDist
   real, intent(out)   :: taus(:)
  
   integer     :: i, nrays, nsides, index
   real        :: normCompanion, theta0, unitCompanion(size(points(:,1))), theta, root, dist
   real, dimension(:,:), allocatable :: dirs
   real, dimension(:,:), allocatable :: listsOfDists, listsOfTaus
   integer, dimension(:), allocatable :: indices
   real, dimension(:), allocatable    :: tau, dists
   integer, dimension(:), allocatable :: listOfPoints

   if (present(companion) .and. present(R)) then
    normCompanion = norm2(points(:,companion)-points(:,primary))
    theta0 = asin(R/normCompanion)
    unitCompanion = (points(:,companion)-points(:,primary))/normCompanion
    nrays = 12*4**minOrder + 12*refineLevel-9 + int(24*2**(refineLevel+log(theta0)/log(2.)))
    nsides = 2**(minOrder+refineLevel)
    allocate(dirs(size(points(:,1)), nrays))
    allocate(indices(12*4**(minOrder+refineLevel)))
    allocate(listsOfDists(200, nrays))
    allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
    allocate(tau(size(listsOfDists(:,1))))
    allocate(dists(size(listsOfDists(:,1))))
    allocate(listOfPoints(size(listsOfDists(:,1))))
    call get_rays(normCompanion,R,minOrder,refineLevel,dirs, indices, nrays)

    !$omp parallel do private(tau,dists,dist,root,theta,listOfPoints)
    do i = 1, nrays
      listOfPoints=0
      tau=0.
      dists=0.
      if (present(companion) .and. present(R)) then
         theta = acos(dot_product(unitCompanion, dirs(:,i)))
         if (theta < theta0) then
            root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
            dist = normCompanion*cos(theta)-root
            call ray_tracer(primary, dirs(:,i), points, xyzh, neighbors, Rstar, listOfPoints, dists, dist)
         else
           call ray_tracer(primary, dirs(:,i), points, xyzh, neighbors, Rstar, listOfPoints, dists, maxDist)
         endif
      else
         call ray_tracer(primary, dirs(:,i), points, xyzh, neighbors, Rstar, listOfPoints, dists, maxDist)
      endif
      call get_tau_list(listOfPoints, dists, opacities, tau)
      !$omp critical
      listsOfTaus(:,i) = tau
      listsOfDists(:,i) = dists
      !$omp end critical
   enddo
   
   taus = 0.
   !$omp parallel do private(index)
   do i = 1, size(taus(:))
     call vec2pix_nest(nsides, points(:,i)-points(:,primary), index)
     index = indices(index + 1)
     if (index /= 0) then
     call get_tau_outwards(points(:,i), points(:,primary), listsOfTaus(:,index), listsOfDists(:,index), dirs(:,index), taus(i))
     endif
   enddo

   else
      call get_all_tau_outwards(primary, points, xyzh, neighbors, opacities, &
      Rstar, minOrder+refineLevel, taus, companion, R, maxDist)
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
 allocate(refine2(n))
 do i=1, n
    circ(i,1) = dist*cos(theta)
    circ(i,2) = dist*sin(theta)*cos(2*PI*i/size(circ(:,i)))
    circ(i,3) = dist*sin(theta)*sin(2*PI*i/size(circ(:,i)))
 enddo
 
 do i=1, size(circ(:,i))
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
    do j=1, size(circ(:,1))
       call vec2pix_nest(2**(minOrder+k),circ(j,:),refine(j))
    enddo
    do i = i, ind2-1
       do j = 1,4
          if (any(4*refine2(mod(i-1,n)+1)+j-1 == refine)) then
             refine2(mod(ind2-1,n)+1) = 4*refine2(mod(i-1,n)+1)+j-1
             ind2 = ind2+1
          else
             call pix2vec_nest(2**(minOrder+k),4*refine2(mod(i-1,n)+1)+j-1,vecs(:,ind))
             indices(4**(refineLevel-k)*(4*refine2(mod(i-1,n)+1)+j-1)+1:4**(refineLevel-k)*(4*(refine2(mod(i-1,n)+1))+j)) = ind
             ind = ind+1
          endif
       enddo
    enddo
 enddo
 do i = i, ind2-1
    do j = 1,4
       call pix2vec_nest(2**(minOrder+refineLevel),4*refine2(mod(i-1,n)+1)+j-1,vecs(:,ind))
       indices(4*(refine2(mod(i-1,n)+1))+j) = ind
       ind = ind+1
    enddo
 enddo
 nrays = ind-1
  end subroutine get_rays

  !*********************************************************************!
  !***************************   OUTWARDS   ****************************!
  !*********************************************************************!
  
  subroutine get_all_tau_outwards(primary, points, xyzh, neighbors, opacities, Rstar, order, taus, companion, R, maxDist)
     integer, intent(in) :: primary, neighbors(:,:), order
     integer, optional   :: companion
     real, intent(in)    :: points(:,:), opacities(:), Rstar, xyzh(:,:)
     real, optional      :: R, maxDist
     real, intent(out)   :: taus(:)
    
     integer  :: i, nrays, nsides, index
     real     :: normCompanion = 0., theta0 = 0., unitCompanion(3) = 0., vec(3)
     real     :: theta, root, dist, dir(size(points(:,1))), phi = 0., cosphi = 0., sinphi = 0.
     real, dimension(:,:), allocatable  :: dirs
     real, dimension(:,:), allocatable  :: listsOfDists, listsOfTaus
     real, dimension(:), allocatable    :: tau, dists
     integer, dimension(:), allocatable :: listOfPoints
     nrays = 12*4**order
     nsides = 2**order
     
     allocate(dirs(size(points(:,1)), nrays))
     allocate(listsOfDists(200, nrays))
     allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
     allocate(tau(size(listsOfDists(:,1))))
     allocate(dists(size(listsOfDists(:,1))))
     allocate(listOfPoints(size(listsOfDists(:,1))))

     if (present(companion) .and. present(R)) then
      unitCompanion = points(:,companion)-points(:,primary)
      normCompanion = norm2(unitCompanion)
      theta0 = asin(R/normCompanion)
      unitCompanion = unitCompanion/normCompanion
      phi = atan(unitCompanion(2)/unitCompanion(1))
      cosphi = cos(phi)
      sinphi = sin(phi)
     endif
     
     !$omp parallel do private(dir,tau,dists,dist,root,theta,listOfPoints)
     do i = 1, nrays
       call pix2vec_nest(nsides, i-1, dir)
      !  dir = (/cosphi*dir(1) - sinphi*dir(2),sinphi*dir(1) + cosphi*dir(2), dir(3)/)
       listOfPoints=0
       tau=0.
       dists=0.
       if (present(companion) .and. present(R)) then
          theta = acos(dot_product(unitCompanion, dir))
          if (theta < theta0) then
             root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
             dist = normCompanion*cos(theta)-root
             call ray_tracer(primary, dir, points, xyzh, neighbors, Rstar, listOfPoints, dists, dist)
          else
            call ray_tracer(primary, dir, points, xyzh, neighbors, Rstar, listOfPoints, dists, maxDist)
          endif
       else
          call ray_tracer(primary, dir, points, xyzh, neighbors, Rstar, listOfPoints, dists, maxDist)
       endif
       call get_tau_list(listOfPoints, dists, opacities, tau)
       !$omp critical
       dirs(:,i) = dir
       listsOfTaus(:,i) = tau
       listsOfDists(:,i) = dists
       !$omp end critical
      enddo

     taus = 0.
     !$omp parallel do private(index,vec)
     do i = 1, size(taus)
       vec = points(:,i)-points(:,primary)
      !  vec = (/cosphi*vec(1) + sinphi*vec(2),-sinphi*vec(1) + cosphi*vec(2), vec(3)/)
       call vec2pix_nest(nsides, vec, index)
       index = index + 1
       call get_tau_outwards(points(:,i), points(:,primary), listsOfTaus(:,index), listsOfDists(:,index), dirs(:,index), taus(i))
     enddo
    end subroutine get_all_tau_outwards
    
    subroutine get_tau_outwards(point, primary, listOfTaus, listOfDists, ray, tau)
       real, intent(in)    :: point(:), primary(:), listOfTaus(:), listOfDists(:), ray(:)
       real, intent(out)   :: tau
    
       integer :: i
       real    :: dist
       
       dist = dot_product(point-primary, ray)
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
    
    subroutine get_tau_list(listOfPoints, listOfDists, opacities, listOfTaus)
       integer, intent(in) :: listOfPoints(:)
       real, intent(in)    :: listOfDists(:), opacities(:)
       real, intent(out)   :: listOfTaus(:)
    
       integer     :: i
       real        :: opacity
    
       listOfTaus = 0.
       i = 2
       do while (listOfPoints(i) /= 0)
          opacity = (opacities(listOfPoints(i)) + opacities(listOfPoints(i-1)))/2
          listOfTaus(i) = listOfTaus(i-1)+(listOfDists(i)-listOfDists(i-1))*opacity
          i = i+1
       enddo
    end subroutine get_tau_list
    
    subroutine ray_tracer(point, ray, points, xyzh, neighbors, Rstar, listOfPoints, listOfDist, maxDist)
     use linklist, only:getneigh_pos,ifirstincell,listneigh
     real, intent(in)     :: points(:,:), ray(size(points(:,1))), Rstar, xyzh(:,:)
     real, optional       :: maxDist
     integer, intent(in)  :: point, neighbors(:,:)
     real, intent(out)    :: listOfDist(:)
     integer, intent(out) :: listOfPoints(:)
    
     integer, parameter :: maxcache = 0
     real, allocatable  :: xyzcache(:,:)
     integer :: nneigh
     real :: dist, h
     integer :: next, i

     dist = Rstar
     listOfPoints(1)=point
     listOfDist(1)=dist
     h = Rstar/100.

     next=0
     do while (next==0)
      h = h*2.
      call getneigh_pos(points(:,point)+Rstar*ray,0.,h,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
      call find_next(points(:,point), ray, dist, points, listneigh, next, nneigh)
     enddo

     i = 1
     do while (hasNext(next,dist,maxDist))
        i = i + 1
        listOfPoints(i) = next
        listOfDist(i)=dist
        call find_next(points(:,point), ray, dist, points, neighbors(next,:), next)
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
    
    subroutine get_all_tau_inwards(primary, points, xyzh, neighbors, opacities, Rstar, taus, companion, R)
     real, intent(in)    :: points(:,:), opacities(:), Rstar, xyzh(:,:)
     real, optional      :: R
     integer, intent(in) :: primary, neighbors(:,:)
     integer, optional   :: companion
     real, intent(out)   :: taus(:)
    
     integer :: i
     real    :: normCompanion, theta0, unitCompanion(3), norm, theta, root, norm0, tau
     integer :: k, next, previous
     real    :: ray(3), nextDist, previousDist, maxDist
#ifdef INLINE2
     real :: vec(3),trace_point(3),tempdist,raydist,dmin,dist
     integer :: j,nneigh,nextj
#endif
     
     normCompanion = 0.
     theta0 = 0.
     unitCompanion = 0.
     if (present(companion) .and. present(R)) then
      normCompanion = norm2(points(:,companion)-points(:,primary))
      theta0 = asin(R/normCompanion)
      unitCompanion = (points(:,companion)-points(:,primary))/normCompanion
     endif
    
    !$omp parallel do &
    !$omp default(none)&
    !$omp private(tau,norm,theta,root,norm0)&
#ifdef INLINE2
    !$omp private(j,vec,trace_point,nneigh,tempdist,raydist,dmin,dist,nextj) &
#endif
    !$omp shared(taus,companion,R,theta0,primary,opacities,Rstar,&
#ifdef INLINE
    !$omp        neighbors,points,unitCompanion,normCompanion) &
    !$omp private(ray,maxDist,next,nextDist,k,previous,previousDist)
#else
     !$omp        neighbors,points,unitCompanion,normCompanion)
#endif   

     do i = 1, size(taus)
      if (present(companion) .and. present(R)) then
         norm = norm2(points(:,i)-points(:,primary))
         theta = acos(dot_product(unitCompanion, points(:,i)-points(:,primary))/norm)
         if (theta < theta0) then
            root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
            norm0 = normCompanion*cos(theta)-root
            if (norm > norm0) then
                  tau = 1e10
            else
               call get_tau_inwards(i, primary, points, neighbors, opacities, Rstar, tau)
            endif
         else
#ifdef INLINE
            ray = points(:,primary) - points(:,i)
            maxDist = norm2(ray)
            ray = ray / maxDist
            maxDist=max(maxDist-Rstar,0.)
            next=i
            nextDist=0.
            
            tau = 0.
            k=1
            do while (next /= primary .and. next /=0)
                k = k + 1
                previous = next
                previousDist = nextDist

#ifdef INLINE2               
!    subroutine find_next(inpoint, ray, dist, points, neighbors, next, nneighin)
                dmin = huge(0.)
                nneigh = size(neighbors(next,:))
                Dist=nextdist
                nextj = next
                next=0
                trace_point = points(:,i) + nextDist*ray

                j = 1
                dist = nextdist
                do while (j <= nneigh .and. neighbors(nextj,j) /= 0)
                    vec=points(:,neighbors(nextj,j)) - trace_point
                    tempdist = dot_product(vec,ray)
                    if (tempdist>0.) then
                        raydist = dot_product(vec,vec) - tempdist**2
                        if (raydist < dmin) then
                            dmin = raydist
                            next = neighbors(nextj,j)
                            nextdist = dist+tempdist
                        end if
                    end if
                    j = j+1
                enddo
#else                
                
                call find_next(points(:,i), ray, nextDist, points, neighbors(next,:), next)
#endif
                if (nextDist .gt. maxDist) then
                    next = primary
                    nextDist = maxDist
                end if
                tau = tau + (nextDist-previousDist)*(opacities(next)+opacities(previous))/2
            enddo
#else
            call get_tau_inwards(i, primary, points, neighbors, opacities, Rstar, tau)
#endif
         endif
      else
         call get_tau_inwards(i, primary, points, neighbors, opacities, Rstar, tau)
      endif
      taus(i) = tau
    enddo
    !$omp  end parallel do
    end subroutine get_all_tau_inwards
    
    subroutine get_tau_inwards(secondary, primary, points, neighbors, opacities, Rstar, tau)
     real, intent(in)    :: points(:,:), opacities(:), Rstar
     integer, intent(in) :: primary, secondary, neighbors(:,:)
     real, intent(out)   :: tau
    
     integer :: i, next, previous
     real    :: ray(size(points(:,1))), nextDist, previousDist, maxDist

     ray = points(:,primary) - points(:,secondary)
     maxDist = norm2(ray)
     ray = ray / maxDist
     maxDist=max(maxDist-Rstar,0.)
     next=secondary
     nextDist=0.
     
     tau = 0.
     i=1
     do while (next /= primary .and. next /=0)
        i = i + 1
        previous = next
        previousDist = nextDist
        call find_next(points(:,secondary), ray, nextDist, points, neighbors(next,:), next)
        if (nextDist .gt. maxDist) then
            next = primary
            nextDist = maxDist
        end if
        tau = tau + (nextDist-previousDist)*(opacities(next)+opacities(previous))/2
     enddo
    end subroutine get_tau_inwards

    !*********************************************************************!
    !****************************   COMMON   *****************************!
    !*********************************************************************!
    
    subroutine find_next(inpoint, ray, dist, points, neighbors, next, nneighin)
     integer, intent(in)  :: neighbors(:)
     real, intent(in)     :: points(:,:), inpoint(:), ray(:)
     integer, intent(out) :: next
     real, intent(inout)  :: dist
     integer, optional    :: nneighin
    
     real                 :: trace_point(3), dmin, vec(3), tempdist, raydist
     real                 :: nextdist
     integer              :: i, nneigh
     
     dmin = huge(0.)
     if (present(nneighin)) then
      nneigh = nneighin
     else
      nneigh = size(neighbors)
     endif

     next=0
     nextDist=dist
     trace_point = inpoint + dist*ray

     i = 1
     do while (i <= nneigh .and. neighbors(i) /= 0)
        vec=points(:,neighbors(i)) - trace_point
        tempdist = dot_product(vec,ray)
        if (tempdist>0.) then
            raydist = dot_product(vec,vec) - tempdist**2
            if (raydist < dmin) then
                dmin = raydist
                next = neighbors(i)
                nextdist = dist+tempdist
            end if
        end if
        i = i+1
     enddo
     dist=nextdist
    end subroutine find_next
    
    end module raytracer
