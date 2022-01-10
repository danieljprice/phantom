module raytracer
   use healpix
   implicit none
   public :: get_all_tau_outwards, get_all_tau_inwards
   private
  contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!   OUTWARDS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine get_all_tau_optimised(primary, points, neighbors, opacities, minOrder, refineLevel, taus, companion, R, maxDist)
     integer, intent(in) :: primary, neighbors(:,:), minOrder, refineLevel
     integer, optional   :: companion
     real, intent(in)    :: points(:,:), opacities(:)
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
             call ray_tracer(primary, dirs(:,i), points, neighbors, listsOfPoints(:,i), listsOfDists(:,i), dist)
          else
            call ray_tracer(primary, dirs(:,i), points, neighbors, listsOfPoints(:,i), listsOfDists(:,i), maxDist)
          endif
       else
          call ray_tracer(primary, dirs(:,i), points, neighbors, listsOfPoints(:,i), listsOfDists(:,i), maxDist)
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
  
  subroutine get_all_tau_outwards(primary, npart, points, neighbors, opacities, order, taus, companion, R, maxDist)
     integer, intent(in) :: primary, npart, neighbors(:,:), order
     integer, optional   :: companion
     real, intent(in)    :: points(:,:), opacities(:)
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
       nrays = 12*4**order
       nsides = 2**order
     else
       nrays = 12*4**order
       nsides = 2**order
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
             call ray_tracer(primary, dirs(:,i), points, neighbors, listsOfPoints(:,i), listsOfDists(:,i), dist)
          else
            call ray_tracer(primary, dirs(:,i), points, neighbors, listsOfPoints(:,i), listsOfDists(:,i), maxDist)
          endif
       else
          call ray_tracer(primary, dirs(:,i), points, neighbors, listsOfPoints(:,i), listsOfDists(:,i), maxDist)
       endif
       call get_tau_list(listsOfPoints(:,i), listsOfDists(:,i), opacities, listsOfTaus(:,i))
     enddo

     taus = 0.
     do i = 1, size(taus)-2
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
       i = 2
       do while (dist .gt. listOfDists(i) .and. listOfDists(i) .gt. 0.)
          i = i+1
       enddo
       
       if (dist .gt. listOfDists(i)) then
          tau = listOfTaus(i-1)!1000.
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
    
    subroutine ray_tracer(point, ray, points, neighbors, listOfPoints, listOfDist, maxDist)
     real, intent(in)     :: points(:,:), ray(size(points(:,1)))
     real, optional       :: maxDist
     integer, intent(in)  :: point, neighbors(:,:)
     real, intent(out)    :: listOfDist(:)
     integer, intent(out) :: listOfPoints(:)
    
     real :: dist
     integer :: next, i
     dist=0.

     listOfPoints(1) = point
     listOfDist(1) = 0.
     call find_next(points(:,point), ray, dist, points, neighbors(point,:), next)
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!   INWARDS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine get_all_tau_inwards(primary, points, neighbors, opacities, taus, companion, R)
     real, intent(in)    :: points(:,:), opacities(:)
     real, optional      :: R
     integer, intent(in) :: primary, neighbors(:,:)
     integer, optional   :: companion
     real, intent(out)   :: taus(:)
    
     integer :: i, j, n = 10, mod, index
     real    :: normCompanion, theta0, unitCompanion(size(points(:,1))), norm, theta, root, norm0
     mod = int(size(taus)/n)+1
     
     if (present(companion) .and. present(R)) then
      normCompanion = norm2(points(:,companion)-points(:,primary))
      theta0 = asin(R/normCompanion)
      unitCompanion = (points(:,companion)-points(:,primary))/normCompanion
     endif
    
     print*,0,'%'
     do i = 1, n
       do j = 1, mod
            index = mod*(i-1)+j
            if (index <= size(taus)) then
               if (present(companion) .and. present(R)) then
                  norm = norm2(points(:,index)-points(:,primary))
                  theta = acos(dot_product(unitCompanion, points(:,index)-points(:,primary))/norm)
                  if (theta < theta0) then
                     root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+R**2)
                     norm0 = normCompanion*cos(theta)-root
                     if (norm > norm0) then
                           taus(index) = 1000.
                     else
                        call get_tau_inwards(index, primary, points, neighbors, opacities, taus(index))
                     endif
                  else
                     call get_tau_inwards(index, primary, points, neighbors, opacities, taus(index))
                  endif
               else
                  call get_tau_inwards(index, primary, points, neighbors, opacities, taus(index))
               endif
            end if
       enddo
       print*,i*100/n,'%'
     enddo
    end subroutine get_all_tau_inwards
    
    subroutine get_tau_inwards(secondary, primary, points, neighbors, opacities, tau)
     real, intent(in)    :: points(:,:), opacities(:)
     integer, intent(in) :: primary, secondary, neighbors(:,:)
     real, intent(out)   :: tau
    
     integer :: listOfPoints(10000), i
     real    :: listOfDist(size(listOfPoints)), ds, opacity
     listOfPoints = 0
     listOfDist = 0.
     
     call ray_between_points(secondary, primary, points, neighbors, listOfPoints, listOfDist)
     tau = 0
     i=2
     do while (listOfPoints(i) /= 0)
        ds = listOfDist(i)-listOfDist(i-1)
        opacity = (opacities(listOfPoints(i))+opacities(listOfPoints(i-1)))/2
        tau = tau + opacity * ds
        i = i+1
     enddo
    end subroutine get_tau_inwards
    
    subroutine ray_between_points(secondary, primary, points, neighbors, listOfPoints, listOfDist)
     real, intent(in)     :: points(:,:)
     integer, intent(in)  :: primary, secondary, neighbors(:,:)
     real, intent(out)    :: listOfDist(:)
     integer, intent(out) :: listOfPoints(:)
    
     real    :: ray(size(points(:,1)))
     real    :: dist, maxdist
     integer :: next, i
     dist=0.
     
     ray = points(:,primary) - points(:,secondary)
     maxdist = norm2(ray)
     ray = ray / maxdist
     next=secondary
     i=1
     do while (next /= primary)
        listOfPoints(i) = next
        listOfDist(i) = dist
        i = i + 1
        call find_next(points(:,secondary), ray, dist, points, neighbors(next,:), next)
        if (dist .gt. maxdist) then
            next = primary
            dist = maxdist
        end if
     enddo
     listOfPoints(i) = next
     listOfDist(i) = dist
    end subroutine ray_between_points

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