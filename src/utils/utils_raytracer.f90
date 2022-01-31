module raytracer
   use healpix
   implicit none
   public :: get_all_tau_outwards, get_all_tau_inwards
   private
  contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!   OUTWARDS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine get_all_tau_optimised(primary, points, xyzh, neighbors, opacities, Rstar, &
                                   minOrder, refineLevel, taus, companion, R, maxDist)
     integer, intent(in) :: primary, neighbors(:,:), minOrder, refineLevel
     integer, optional   :: companion
     real, intent(in)    :: points(:,:), opacities(:), Rstar, xyzh(:,:)
     real, optional      :: R, maxDist
     real, intent(out)   :: taus(:)
    
     integer     :: i, nrays, nsides, index
     real        :: normCompanion, theta0, unitCompanion(size(points(:,1))), theta, root, dist
     real, dimension(:,:), allocatable :: dirs
     integer, dimension(:,:), allocatable :: listsOfPoints
     real, dimension(:,:), allocatable :: listsOfDists, listsOfTaus
     logical     :: threeD
     
     threeD = size(points(:,1)) .eq. 3
    
     if (threeD) then
       nrays = 12*4**minOrder
       nsides = 2**minOrder
     else
       nrays = 12*4**minOrder
       nsides = 2**minOrder
     endif
     
     allocate(dirs(size(points(:,1)), nrays))
     allocate(listsOfPoints(10000, nrays))
     allocate(listsOfDists(size(listsOfPoints(:,1)), nrays))
     allocate(listsOfTaus(size(listsOfPoints(:,1)), nrays))

     if (present(companion) .and. present(R)) then
      normCompanion = norm2(points(:,companion)-points(:,primary))
      theta0 = asin(R/normCompanion)
      unitCompanion = (points(:,companion)-points(:,primary))/normCompanion
     endif
     
     do i = 1, nrays
       if (threeD) then
          call pix2vec_nest(nsides, i-1, dirs(:,i))
       else
          call pix2vec_nest_2d(nsides, i-1, dirs(:,i))
       endif
       if (present(companion) .and. present(R)) then
          theta = acos(dot_product(unitCompanion, dirs(:,i)))
          if (theta < theta0) then
             root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
             dist = normCompanion*cos(theta)-root
             call ray_tracer(primary, dirs(:,i), points, xyzh, neighbors, Rstar, listsOfPoints(:,i), listsOfDists(:,i), dist)
          else
            call ray_tracer(primary, dirs(:,i), points, xyzh, neighbors, Rstar, listsOfPoints(:,i), listsOfDists(:,i), maxDist)
          endif
       else
          call ray_tracer(primary, dirs(:,i), points, xyzh, neighbors, Rstar, listsOfPoints(:,i), listsOfDists(:,i), maxDist)
       endif
       call get_tau_list(listsOfPoints(:,i), listsOfDists(:,i), opacities, listsOfTaus(:,i))
     enddo
     taus = 0.
     do i = 2, size(points(1,:))
       if (threeD) then
          call vec2pix_nest(nsides, points(:,i)-points(:,primary), index)
       else
          call vec2pix_nest_2d(nsides, points(:,i)-points(:,primary), index)
       endif
       index = index + 1
       call get_tau_outwards(points(:,i), points(:,primary), listsOfTaus(:,index), listsOfDists(:,index), dirs(:,index), taus(i))
     enddo
    end subroutine get_all_tau_optimised

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!   OUTWARDS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine get_all_tau_outwards(primary, points, xyzh, neighbors, opacities, Rstar, order, taus, companion, R, maxDist)
     integer, intent(in) :: primary, neighbors(:,:), order
     integer, optional   :: companion
     real, intent(in)    :: points(:,:), opacities(:), Rstar, xyzh(:,:)
     real, optional      :: R, maxDist
     real, intent(out)   :: taus(:)
    
     integer  :: i, nrays, nsides, index
     real     :: normCompanion, theta0, unitCompanion(size(points(:,1))), theta, root, dist, dir(size(points(:,1)))
     real, dimension(:,:), allocatable  :: dirs
     real, dimension(:,:), allocatable  :: listsOfDists, listsOfTaus
     real, dimension(:), allocatable    :: tau, dists
     integer, dimension(:), allocatable :: listOfPoints
     logical     :: threeD

     threeD = size(points(:,1)) .eq. 3
    
     if (threeD) then
       nrays = 12*4**order
       nsides = 2**order
     else
       nrays = 12*4**order
       nsides = 2**order
     endif
     
     allocate(dirs(size(points(:,1)), nrays))
     allocate(listsOfDists(200, nrays))
     allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
     allocate(tau(size(listsOfDists(:,1))))
     allocate(dists(size(listsOfDists(:,1))))
     allocate(listOfPoints(size(listsOfDists(:,1))))

     if (present(companion) .and. present(R)) then
      normCompanion = norm2(points(:,companion)-points(:,primary))
      theta0 = asin(R/normCompanion)
      unitCompanion = (points(:,companion)-points(:,primary))/normCompanion
     endif
     
     !$omp parallel do default(none) &
     !$omp private(dir,tau,dists,dist,root,theta,listOfPoints) &
     !$omp shared(threeD,nsides,companion,R,unitCompanion,theta0, &
     !$omp normCompanion,points,xyzh,neighbors,Rstar,primary,maxDist,opacities,dirs,listsOfTaus,listsOfDists)
     do i = 1, nrays
       if (threeD) then
          call pix2vec_nest(nsides, i-1, dir)
       else
          call pix2vec_nest_2d(nsides, i-1, dir)
       endif
       listOfPoints=0
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
     !$omp parallel do private(index)
     do i = 1, size(taus)
       if (threeD) then
          call vec2pix_nest(nsides, points(:,i)-points(:,primary), index)
       else
          call vec2pix_nest_2d(nsides, points(:,i)-points(:,primary), index)
       endif
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
       i = 1
       do while (dist .gt. listOfDists(i) .and. listOfDists(i) .gt. 0.)
          i = i+1
       enddo
       
       if (dist .gt. listOfDists(i)) then
          tau = 1e10
       else if (i == 1) then
          tau = 0.
       else
          tau = listOfTaus(i-1)+(listOfTaus(i)-listOfTaus(i-1))/(listOfDists(i)-listOfDists(i-1))*(dist-listOfDists(i-1))
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
     real :: dist, r
     integer :: next, i

     r = Rstar/100.
     call getneigh_pos(points(:,point)+Rstar*ray,0.,r,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
     do while (nneigh==0)
      r = r*2.
      call getneigh_pos(points(:,point)+Rstar*ray,0.,r,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
     enddo
     next = point
     dist = Rstar

     i = 1
     listOfPoints(i) = next
     listOfDist(i)=dist
     call find_next(points(:,point), ray, dist, points, listneigh, next)
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!   INWARDS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine get_all_tau_inwards(primary, points, xyzh, neighbors, opacities, Rstar, taus, companion, R)
     real, intent(in)    :: points(:,:), opacities(:), Rstar, xyzh(:,:)
     real, optional      :: R
     integer, intent(in) :: primary, neighbors(:,:)
     integer, optional   :: companion
     real, intent(out)   :: taus(:)
    
     integer :: i
     real    :: normCompanion, theta0, unitCompanion(size(points(:,1))), norm, theta, root, norm0, tau
     if (present(companion) .and. present(R)) then
      normCompanion = norm2(points(:,companion)-points(:,primary))
      theta0 = asin(R/normCompanion)
      unitCompanion = (points(:,companion)-points(:,primary))/normCompanion
     endif
    
     !$omp parallel do private(tau)
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
            call get_tau_inwards(i, primary, points, neighbors, opacities, Rstar, tau)
         endif
      else
         call get_tau_inwards(i, primary, points, neighbors, opacities, Rstar, tau)
      endif
      !$omp critical
      taus(i) = tau 
      !$omp end critical
     enddo
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
     
     tau = 0
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
        if (tau<0) then
        endif
     enddo
    end subroutine get_tau_inwards

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!   COMMON   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine find_next(inpoint, ray, dist, points, neighbors, next)
     integer, intent(in)  :: neighbors(:)
     real, intent(in)     :: points(:,:), inpoint(:), ray(:)
     integer, intent(out) :: next
     real, intent(inout)  :: dist
    
     real                 :: trace_point(size(inpoint)), min, vec(size(inpoint)), tempdist, raydist
     real                 :: nextdist
     integer              :: i
     character(len=3)     :: inf="INF"
     read(inf,*) min
     
     next=0
     trace_point = inpoint + dist*ray
     i = 1
     do while (i <= size(neighbors) .and. neighbors(i) /= 0)
        vec=points(:,neighbors(i)) - trace_point
        tempdist = dot_product(vec,ray)
        if (tempdist>0) then
            raydist = dot_product(vec,vec) - tempdist**2
            if (raydist < min) then
                min = raydist
                next = neighbors(i)
                nextdist = dist+tempdist
            end if
        end if
        i = i+1
     enddo
     dist=nextdist
    end subroutine find_next
    
    end module raytracer