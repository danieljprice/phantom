!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module raytracer_all
!
! raytracer_all
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
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the kappa of all SPH particles
 !  IN: Rstar:           The radius of the star
 !  IN: minOrder:        The minimal order in which the rays are sampled
 !  IN: refineLevel:     The amount of orders in which the rays can be
 !                       sampled deeper
 !  IN: refineScheme:    The refinement scheme used for adaptive ray selection
 !+
 !  OUT: taus:           The list of optical depths for each particle
 !+
 !  OPT: companion:      The xyz coordinates of the companion
 !  OPT: Rcomp:          The radius of the companion
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_adaptive(npart, primary, xyzh, kappa, Rstar, minOrder,&
                                 refineLevel, refineScheme, taus, companion, Rcomp)
 integer, intent(in) :: npart, minOrder, refineLevel, refineScheme
 real, intent(in)    :: primary(3), kappa(:), xyzh(:,:), Rstar
 real, optional      :: Rcomp, companion(3)
 real, intent(out)   :: taus(:)

 integer     :: i, nrays, nsides, index
 real        :: normCompanion, theta0, unitCompanion(3), theta, root, dist, vec(3), dir(3)
 real, dimension(:,:), allocatable :: dirs
 real, dimension(:,:), allocatable :: listsOfDists, listsOfTaus
 integer, dimension(:), allocatable :: indices, rays_dim
 real, dimension(:), allocatable    :: tau, dists

 if (present(companion) .and. present(Rcomp)) then
    unitCompanion = companion-primary
    normCompanion = norm2(unitCompanion)
    theta0 = asin(Rcomp/normCompanion)
    unitCompanion = unitCompanion/normCompanion

    call get_rays(npart, primary, companion, Rcomp, xyzh, minOrder, refineLevel, refineScheme, dirs, indices, nrays)
    allocate(listsOfDists(200, nrays))
    allocate(listsOfTaus(size(listsOfDists(:,1)), nrays))
    allocate(tau(size(listsOfDists(:,1))))
    allocate(dists(size(listsOfDists(:,1))))
    allocate(rays_dim(nrays))

    !$omp parallel do private(tau,dist,dir,dists,root,theta)
    do i = 1, nrays
       tau=0.
       dists=0.
       dir = dirs(:,i)
       theta = acos(dot_product(unitCompanion, dir))
       if (theta < theta0) then
          root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+Rcomp**2)
          dist = normCompanion*cos(theta)-root
          call ray_tracer(primary, dir, xyzh, kappa, Rstar, tau, dists, rays_dim(i), dist)
       else
          call ray_tracer(primary, dir, xyzh, kappa, Rstar, tau, dists, rays_dim(i))
       endif
       listsOfTaus(:,i) = tau
       listsOfDists(:,i) = dists
    enddo
    !$omp end parallel do

    nsides = 2**(minOrder+refineLevel)
    taus = 0.
    !$omp parallel do private(index,vec)
    do i = 1, npart
       vec = xyzh(1:3,i)-primary
       call vec2pix_nest(nsides, vec, index)
       index = indices(index + 1)
       call get_tau_on_ray(norm2(vec), listsOfTaus(:,index), listsOfDists(:,index), rays_dim(index), taus(i))
    enddo
    !$omp end parallel do

 else
    call get_all_tau_outwards_single(npart, primary, xyzh, kappa, &
      Rstar, minOrder+refineLevel, 0, taus)
 endif
end subroutine get_all_tau_adaptive

 !--------------------------------------------------------------------------
 !+
 !  Return all the directions of the rays that need to be traced for the
 !  adaptive ray-tracing scheme
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: companion:       The xyz coordinates of the companion
 !  IN: Rcomp:           The radius of the companion
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: minOrder:        The minimal order in which the rays are sampled
 !  IN: refineLevel:     The amount of orders in which the rays can be
 !                       sampled deeper
 !  IN: refineScheme:    The refinement scheme used for adaptive ray selection
 !+
 !  OUT: rays:           A list containing the rays that need to be traced
 !                       in the adaptive ray-tracing scheme
 !  OUT: indices:        A list containing a link between the index in the
 !                       deepest order and the rays in the adaptive ray-tracing scheme
 !  OUT: nrays:          The number of rays after the ray selection
 !+
 !--------------------------------------------------------------------------
subroutine get_rays(npart, primary, companion, Rcomp, xyzh, minOrder, refineLevel, refineScheme, rays, indices, nrays)
 integer, intent(in)  :: npart, minOrder, refineLevel, refineScheme
 real, intent(in)     :: primary(3), companion(3), xyzh(:,:), Rcomp
 real, allocatable, intent(out)    :: rays(:,:)
 integer, allocatable, intent(out) :: indices(:)
 integer, intent(out) :: nrays

 real    :: theta, dist, phi, cosphi, sinphi
 real, dimension(:,:), allocatable  :: circ
 integer :: i, j, minNsides, minNrays, ind,n, maxOrder, max, distr(12*4**(minOrder+refineLevel))
 integer, dimension(:,:), allocatable  :: distrs

 maxOrder = minOrder+refineLevel
 nrays = 12*4**(maxOrder)
 allocate(rays(3, nrays))
 allocate(indices(12*4**(maxOrder)))
 rays = 0.
 indices = 0

 !If there is no refinement, just return the uniform ray distribution
 minNsides = 2**minOrder
 minNrays = 12*4**minOrder
 if (refineLevel == 0) then
    do i=1, minNrays
       call pix2vec_nest(minNsides,i-1, rays(:,i))
       indices(i) = i
    enddo
    return
 endif

 !Fill a list to have the number distribution in angular space
 distr = 0
 !$omp parallel do private(ind)
 do i = 1, npart
    call vec2pix_nest(2**maxOrder, xyzh(1:3, i)-primary, ind)
    distr(ind+1) = distr(ind+1)+1
 enddo
 max = maxval(distr)

 !Make sure the companion is described using the highest refinement
 dist = norm2(primary-companion)
 theta = asin(Rcomp/dist)
 phi = atan2(companion(2)-primary(2),companion(1)-primary(1))
 cosphi = cos(phi)
 sinphi = sin(phi)
 dist = dist*cos(theta)
 n = int(theta*6*2**(minOrder+refineLevel))+4
 allocate(circ(n,3))
 do i=1, n !Define boundary of the companion
    circ(i,1) = dist*cos(theta)
    circ(i,2) = dist*sin(theta)*cos(twopi*i/n)
    circ(i,3) = dist*sin(theta)*sin(twopi*i/n)
    circ(i,:) = (/cosphi*circ(i,1) - sinphi*circ(i,2),sinphi*circ(i,1) + cosphi*circ(i,2), circ(i,3)/)
 enddo
 do i=1, n !Make sure the boundary is maximally refined
    call vec2pix_nest(2**maxOrder,circ(i,:),ind)
    distr(ind+1) = max
 enddo

 !Calculate the number distribution in all the orders needed
 allocate(distrs(12*4**(minOrder+refineLevel),refineLevel+1))
 distrs = 0
 distrs(:,1) = distr
 do i = 1, refineLevel
    do j = 1, 12*4**(maxOrder-i)
       distrs(j,i+1) = distrs(4*j,i)+distrs(4*j+1,i)+distrs(4*j+2,i)+distrs(4*j+3,i)
    enddo
 enddo
 max = maxval(distrs(:,refineLevel+1))+1

 !Fill the rays array walking through the orders
 ind=1

 ! refine half in each order
 if (refineScheme == 1) then
    do i=0, refineLevel-1
       call merge_argsort(distrs(1:12*4**(minOrder+i),refineLevel-i+1), distr)
       do j=1, 6*4**minOrder*2**(i)
          call pix2vec_nest(2**(minOrder+i), distr(j)-1, rays(:,ind))
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

    ! refine overdens regions in each order
 elseif (refineScheme == 2) then
    do i=0, refineLevel-1
       call merge_argsort(distrs(1:12*4**(minOrder+i),refineLevel-i+1), distr)
       j=1
       do while (distrs(distr(j),refineLevel-i+1)<npart/(12*4**(minOrder+i)))
          call pix2vec_nest(2**(minOrder+i), distr(j)-1, rays(:,ind))
          indices(4**(refineLevel-i)*(distr(j)-1)+1:4**(refineLevel-i)*distr(j)) = ind
          ind=ind+1
          distrs(4*(distr(j)-1)+1:4*(distr(j)), refineLevel-i) = max
          j=j+1
       enddo
       do j = j, 12*4**(minOrder+i)
          if (distrs(distr(j),refineLevel-i+1) == max) then
             distrs(4*(distr(j)-1)+1:4*(distr(j)), refineLevel-i) = max
          endif
       enddo
    enddo
 endif

 do i=1, 12*4**maxOrder
    if (distrs(i,1)  /=  max) then
       call pix2vec_nest(2**maxOrder, i-1, rays(:,ind))
       indices(i) = ind
       ind=ind+1
    endif
 enddo
 nrays = ind-1
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
 enddo

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
    enddo
    stepsize=stepsize*2
 enddo
end subroutine merge_argsort

 !*********************************************************************!
 !***************************   OUTWARD    ****************************!
 !*********************************************************************!

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth of each particle, using the uniform outwards
 !  ray-tracing scheme
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the kappa of all SPH particles
 !  IN: Rstar:           The radius of the star
 !  IN: order:           The order in which the rays are sampled
 !  IN: raypolation:     The interpolation scheme used for the ray interpolation
 !+
 !  OUT: taus:           The list of optical depths for each particle
 !+
 !  OPT: companion:      The location of the companion
 !  opt: Rcomp:          The radius of the primary star
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_outwards(npart, primary, xyzh, kappa, Rstar, order, raypolation, taus, companion, Rcomp)
 integer, intent(in) :: npart, order, raypolation
 real, intent(in)    :: primary(3), kappa(:), Rstar, xyzh(:,:)
 real, optional      :: Rcomp, companion(3)
 real, intent(out)   :: taus(:)

 if (present(companion) .and. present(Rcomp)) then
    call get_all_tau_outwards_companion(npart, primary, xyzh, kappa, Rstar, companion, Rcomp, order, raypolation, taus)
 else
    call get_all_tau_outwards_single(npart, primary, xyzh, kappa, Rstar, order, raypolation, taus)
 endif
end subroutine get_all_tau_outwards

 !--------------------------------------------------------------------------
 !+
 !  Calculates the optical depth to each SPH particle, using the uniform outwards
 !  ray-tracing scheme for models containing a single star
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the kappa of all SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !  IN: raypolation:     The interpolation scheme used for the ray interpolation
 !+
 !  OUT: taus:           The array of optical depths to each SPH particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_outwards_single(npart, primary, xyzh, kappa, Rstar, order, raypolation, tau)
 use part, only : isdead_or_accreted
 integer, intent(in) :: npart,order, raypolation
 real, intent(in)    :: primary(3), kappa(:), Rstar, xyzh(:,:)
 real, intent(out)   :: tau(:)

 integer  :: i, nrays, nsides
 real     :: ray_dir(3),part_dir(3)
 real, dimension(:,:), allocatable  :: rays_dist, rays_tau
 integer, dimension(:), allocatable :: rays_dim
 integer, parameter :: ndim = 200

 nrays = 12*4**order ! The number of rays traced given the healpix order
 nsides = 2**order   ! The healpix nsides given the healpix order

 allocate(rays_dist(ndim, nrays))
 allocate(rays_tau(ndim, nrays))
 allocate(rays_dim(nrays))

 !-------------------------------------------
 ! CONSTRUCT the RAYS given the ORDER
 ! and determine the optical depth along them
 !-------------------------------------------

 !$omp parallel default(none) &
 !$omp private(ray_dir) &
 !$omp shared(nrays,nsides,primary,kappa,xyzh,Rstar,rays_dist,rays_tau,rays_dim)
 !$omp do
 do i = 1, nrays
    !returns ray_dir, the unit vector identifying ray (index i-1 becase healpix starts counting at index 0)
    call pix2vec_nest(nsides, i-1, ray_dir)
    !calculate the properties along the ray
    call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,rays_tau(:,i),rays_dist(:,i),rays_dim(i))
 enddo
 !$omp enddo
 !$omp end parallel


 !_----------------------------------------------
 ! DETERMINE the optical depth for each particle
 ! using the values available on the rays
 !-----------------------------------------------

 !$omp parallel default(none) &
 !$omp private(part_dir) &
 !$omp shared(npart,primary,nsides,xyzh,ray_dir,rays_dist,rays_tau,rays_dim,raypolation,tau)
 !$omp do
 do i = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       part_dir = xyzh(1:3,i)-primary
       call interpolate_tau(nsides, part_dir, rays_tau, rays_dist, rays_dim, raypolation, tau(i))
    else
       tau(i) = -99.
    endif
 enddo
 !$omp enddo
 !$omp end parallel
end subroutine get_all_tau_outwards_single

 !--------------------------------------------------------------------------
 !+
 !  Calculates the optical depth to each SPH particle, using the uniform outwards
 !  ray-tracing scheme for models containing primary star and a companion
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the opacity of all the SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: companion:       The xyz coordinates of the companion
 !  IN: Rcomp:           The radius of the companion
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !  IN: raypolation:     The interpolation scheme used for the ray interpolation
 !+
 !  OUT: tau:            The array of optical depths for each SPH particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_outwards_companion(npart, primary, xyzh, kappa, Rstar, companion, Rcomp, order, raypolation, tau)
 use part, only : isdead_or_accreted
 integer, intent(in) :: npart, order, raypolation
 real, intent(in)    :: primary(3), companion(3), kappa(:), Rstar, xyzh(:,:), Rcomp
 real, intent(out)   :: tau(:)

 integer  :: i, nrays, nsides
 real     :: normCompanion,theta0,phi,cosphi,sinphi,theta,sep,root
 real     :: ray_dir(3),part_dir(3),uvecCompanion(3)
 real, dimension(:,:), allocatable  :: dirs
 real, dimension(:,:), allocatable  :: rays_dist, rays_tau
 integer, dimension(:), allocatable :: rays_dim
 integer, parameter :: ndim = 200

 nrays = 12*4**order ! The number of rays traced given the healpix order
 nsides = 2**order   ! The healpix nsides given the healpix order

 allocate(dirs(3, nrays))
 allocate(rays_dist(ndim, nrays))
 allocate(rays_tau(ndim, nrays))
 allocate(rays_dim(nrays))

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

 !$omp parallel default(none) &
 !$omp private(ray_dir,theta,root,sep) &
 !$omp shared(nrays,nsides,primary,kappa,xyzh,Rstar,Rcomp,rays_dist,rays_tau,rays_dim) &
 !$omp shared(uvecCompanion,normCompanion,cosphi,sinphi,theta0)
 !$omp do
 do i = 1, nrays
    !returns ray_dir, the unit vector identifying ray (index i-1 becase healpix starts counting at index 0)
    call pix2vec_nest(nsides, i-1, ray_dir)
    !rotate ray vectors by an angle = phi so the main axis points to the companion (This is because along the
    !main axis (1,0,0) rays are distributed more uniformally
    ray_dir = (/cosphi*ray_dir(1) - sinphi*ray_dir(2),sinphi*ray_dir(1) + cosphi*ray_dir(2), ray_dir(3)/)
    theta   = acos(dot_product(uvecCompanion, ray_dir))
    !the ray intersects the companion: only calculate tau up to the companion
    if (theta < theta0) then
       root  = sqrt(Rcomp**2-normCompanion**2*sin(theta)**2)
       sep   = normCompanion*cos(theta)-root
       call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,rays_tau(:,i),rays_dist(:,i),rays_dim(i), sep)
    else
       call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,rays_tau(:,i),rays_dist(:,i),rays_dim(i))
    endif
 enddo
 !$omp enddo
 !$omp end parallel

 !-----------------------------------------------
 ! DETERMINE the optical depth for each particle
 ! using the values available on the rays
 !-----------------------------------------------

 !$omp parallel default(none) &
 !$omp private(part_dir) &
 !$omp shared(npart,primary,cosphi,sinphi,nsides,xyzh,ray_dir,rays_dist,rays_tau,rays_dim,raypolation,tau)
 !$omp do
 do i = 1, npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       !vector joining the source to the particle
       part_dir = xyzh(1:3,i)-primary
       part_dir = (/cosphi*part_dir(1) + sinphi*part_dir(2),-sinphi*part_dir(1) + cosphi*part_dir(2), part_dir(3)/)
       call interpolate_tau(nsides, part_dir, rays_tau, rays_dist, rays_dim, raypolation, tau(i))
    else
       tau(i) = -99.
    endif
 enddo
 !$omp enddo
 !$omp end parallel
end subroutine get_all_tau_outwards_companion

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
 !  IN: rays_tau:        2-dimensional array containing the cumulative optical
 !                       depts along each ray
 !  IN: rays_dist:       2-dimensional array containing the distances from the
 !                       primary along each ray
 !  IN: rays_dim:        The vector containing the number of points defined along each ray
 !  IN: raypolation:     The interpolation scheme used for the ray interpolation
 !+
 !  OUT: tau:            The interpolated optical depth at the particle's location
 !+
 !--------------------------------------------------------------------------
subroutine interpolate_tau(nsides, vec, rays_tau, rays_dist, rays_dim, raypolation, tau)
 integer, intent(in) :: nsides, rays_dim(:), raypolation
 real, intent(in)    :: vec(:), rays_tau(:,:), rays_dist(:,:)
 real, intent(out)   :: tau

 integer :: rayIndex, neighbours(8), nneigh, i, k
 real    :: tautemp, ray(3), vectemp(3), weight, tempdist(8), distRay_sq, vec_norm2
 logical :: mask(8)

 ! 1 ray, no interpolation
 if (raypolation==0) then
    call vec2pix_nest(nsides, vec, rayIndex)
    rayIndex = rayIndex + 1
    call get_tau_on_ray(norm2(vec), rays_tau(:,rayIndex), rays_dist(:,rayIndex), rays_dim(rayIndex), tau)

    ! 4 rays, linear interpolation
 elseif (raypolation==1) then
    vec_norm2 = norm2(vec)
    !returns rayIndex, the index of the ray vector that points to the particle (direction vec)
    call vec2pix_nest(nsides, vec, rayIndex)
    !returns ray(3), the unit vector identifying the ray with index number rayIndex
    call pix2vec_nest(nsides, rayIndex, ray)
    vectemp       = vec - vec_norm2*ray
    distRay_sq    = norm2(vectemp)
    call get_tau_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), tautemp)
    if (distRay_sq > 0.) then
       tau    = tautemp/distRay_sq
       weight = 1./distRay_sq
    else
       ! the particle sits exactly on the ray, no need to get the neighbours
       tau    = tautemp
       return
    endif

    !returns the number nneigh and list of vectors (n) neighbouring the ray number index
    call neighbours_nest(nsides, rayIndex, neighbours, nneigh)
    !for each neighbouring ray calculate its distance to the particle
    do i=1,nneigh
       call pix2vec_nest(nsides, neighbours(i), ray)
       vectemp     = vec - vec_norm2*ray
       tempdist(i) = norm2(vectemp)
    enddo
    neighbours       = neighbours+1
    mask             = .true.
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

    ! 9 rays, linear interpolation
 elseif (raypolation==2) then
    vec_norm2 = norm2(vec)
    !returns rayIndex, the index of the ray vector that points to the particle (direction vec)
    call vec2pix_nest(nsides, vec, rayIndex)
    !returns ray(3), the unit vector identifying the ray with index number rayIndex
    call pix2vec_nest(nsides, rayIndex, ray)
    vectemp       = vec - vec_norm2*ray
    distRay_sq    = norm2(vectemp)
    call get_tau_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), tautemp)
    if (distRay_sq > 0.) then
       tau    = tautemp/distRay_sq
       weight = 1./distRay_sq
    else
       ! the particle sits exactly on the ray, no need to get the neighbours
       tau    = tautemp
       return
    endif

    !returns the number nneigh and list of vectors (n) neighbouring the ray number index
    call neighbours_nest(nsides, rayIndex, neighbours, nneigh)
    !for each neighbouring ray calculate its distance to the particle
    do i=1,nneigh
       call pix2vec_nest(nsides, neighbours(i), ray)
       vectemp     = vec - vec_norm2*ray
       tempdist(i) = norm2(vectemp)
    enddo
    neighbours       = neighbours+1
    mask             = .true.
    if (nneigh <8) mask(nneigh+1:8) = .false.
    !take tau contribution from the 3 closest rays
    do i=1,nneigh
       k       = minloc(tempdist,1,mask)
       mask(k) = .false.
       call get_tau_on_ray(vec_norm2, rays_tau(:,neighbours(k)), &
                     rays_dist(:,neighbours(k)), rays_dim(neighbours(k)), tautemp)
       tau    = tau + tautemp/tempdist(k)
       weight = weight + 1./tempdist(k)
    enddo
    tau = tau / weight

    ! 4 rays, square interpolation
 elseif (raypolation==3) then
    vec_norm2 = norm2(vec)
    !returns rayIndex, the index of the ray vector that points to the particle (direction vec)
    call vec2pix_nest(nsides, vec, rayIndex)
    !returns ray(3), the unit vector identifying the ray with index number rayIndex
    call pix2vec_nest(nsides, rayIndex, ray)
    vectemp       = vec - vec_norm2*ray
    distRay_sq    = dot_product(vectemp,vectemp)
    call get_tau_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), tautemp)
    if (distRay_sq > 0.) then
       tau    = tautemp/distRay_sq
       weight = 1./distRay_sq
    else
       ! the particle sits exactly on the ray, no need to get the neighbours
       tau    = tautemp
       return
    endif

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
       k       = minloc(tempdist,1,mask)
       mask(k) = .false.
       call get_tau_on_ray(vec_norm2, rays_tau(:,neighbours(k)), &
                     rays_dist(:,neighbours(k)), rays_dim(neighbours(k)), tautemp)
       tau    = tau + tautemp/tempdist(k)
       weight = weight + 1./tempdist(k)
    enddo
    tau = tau / weight

    ! 9 rays, square interpolation
 elseif (raypolation==4) then
    vec_norm2 = norm2(vec)
    !returns rayIndex, the index of the ray vector that points to the particle (direction vec)
    call vec2pix_nest(nsides, vec, rayIndex)
    !returns ray(3), the unit vector identifying the ray with index number rayIndex
    call pix2vec_nest(nsides, rayIndex, ray)
    vectemp       = vec - vec_norm2*ray
    distRay_sq    = dot_product(vectemp,vectemp)
    call get_tau_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), tautemp)
    if (distRay_sq > 0.) then
       tau    = tautemp/distRay_sq
       weight = 1./distRay_sq
    else
       ! the particle sits exactly on the ray, no need to get the neighbours
       tau    = tautemp
       return
    endif

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
    do i=1,nneigh
       k       = minloc(tempdist,1,mask)
       mask(k) = .false.
       call get_tau_on_ray(vec_norm2, rays_tau(:,neighbours(k)), &
                     rays_dist(:,neighbours(k)), rays_dim(neighbours(k)), tautemp)
       tau    = tau + tautemp/tempdist(k)
       weight = weight + 1./tempdist(k)
    enddo
    tau = tau / weight

    ! 4 rays, cubed interpolation
 elseif (raypolation==5) then
    vec_norm2 = norm2(vec)
    !returns rayIndex, the index of the ray vector that points to the particle (direction vec)
    call vec2pix_nest(nsides, vec, rayIndex)
    !returns ray(3), the unit vector identifying the ray with index number rayIndex
    call pix2vec_nest(nsides, rayIndex, ray)
    vectemp       = vec - vec_norm2*ray
    distRay_sq    = norm2(vectemp)**3
    call get_tau_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), tautemp)
    if (distRay_sq > 0.) then
       tau    = tautemp/distRay_sq
       weight = 1./distRay_sq
    else
       ! the particle sits exactly on the ray, no need to get the neighbours
       tau    = tautemp
       return
    endif

    !returns the number nneigh and list of vectors (n) neighbouring the ray number index
    call neighbours_nest(nsides, rayIndex, neighbours, nneigh)
    !for each neighbouring ray calculate its distance to the particle
    do i=1,nneigh
       call pix2vec_nest(nsides, neighbours(i), ray)
       vectemp     = vec - vec_norm2*ray
       tempdist(i) = norm2(vectemp)**3
    enddo
    neighbours       = neighbours+1
    mask             = .true.
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

    ! 9 rays, cubed interpolation
 elseif (raypolation==6) then
    vec_norm2 = norm2(vec)
    !returns rayIndex, the index of the ray vector that points to the particle (direction vec)
    call vec2pix_nest(nsides, vec, rayIndex)
    !returns ray(3), the unit vector identifying the ray with index number rayIndex
    call pix2vec_nest(nsides, rayIndex, ray)
    vectemp       = vec - vec_norm2*ray
    distRay_sq    = norm2(vectemp)**3
    call get_tau_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), tautemp)
    if (distRay_sq > 0.) then
       tau    = tautemp/distRay_sq
       weight = 1./distRay_sq
    else
       ! the particle sits exactly on the ray, no need to get the neighbours
       tau    = tautemp
       return
    endif

    !returns the number nneigh and list of vectors (n) neighbouring the ray number index
    call neighbours_nest(nsides, rayIndex, neighbours, nneigh)
    !for each neighbouring ray calculate its distance to the particle
    do i=1,nneigh
       call pix2vec_nest(nsides, neighbours(i), ray)
       vectemp     = vec - vec_norm2*ray
       tempdist(i) = norm2(vectemp)**3
    enddo
    neighbours       = neighbours+1
    mask             = .true.
    if (nneigh <8) mask(nneigh+1:8) = .false.
    !take tau contribution from the 3 closest rays
    do i=1,nneigh
       k       = minloc(tempdist,1,mask)
       mask(k) = .false.
       call get_tau_on_ray(vec_norm2, rays_tau(:,neighbours(k)), &
                  rays_dist(:,neighbours(k)), rays_dim(neighbours(k)), tautemp)
       tau    = tau + tautemp/tempdist(k)
       weight = weight + 1./tempdist(k)
    enddo
    tau = tau / weight
 endif
end subroutine interpolate_tau


 !--------------------------------------------------------------------------
 !+
 !  Interpolation of the optical depth for an arbitrary point on the ray,
 !  with a given distance to the starting point of the ray.
 !+
 !  IN: distance:        The distance from the staring point of the ray to a
 !                       point on the ray
 !  IN: tau_along_ray:   The vector of cumulative optical depths along the ray
 !  IN: dist_along_ray:  The vector of distances from the primary along the ray
 !  IN: len:             The length of listOfTau and listOfDist
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
    tau = 0.
 elseif (distance  >  dist_along_ray(len)) then
    tau = 99.
 else
    L = 2
    R = len-1
    !bysection search for the index of the closest ray location to the particle
    do while (L < R)
       m = (L + R)/2
       if (dist_along_ray(m) > distance) then
          R = m
       else
          L = m + 1
       endif
    enddo
    !interpolate linearly ray properties to get the particle's optical depth
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
 !                       optical depts will be calculated
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the particles opacity
 !  IN: Rstar:           The radius of the primary star
 !+
 !  OUT: taus:           The distribution of optical depths throughout the ray
 !  OUT: listOfDists:    The distribution of distances throughout the ray
 !  OUT: len:            The length of tau_along_ray and dist_along_ray
 !+
 !  OPT: maxDistance:    The maximal distance the ray needs to be traced
 !+
 !--------------------------------------------------------------------------
subroutine ray_tracer(primary, ray, xyzh, kappa, Rstar, tau_along_ray, dist_along_ray, len, maxDistance)
 use linklist, only:getneigh_pos,ifirstincell,listneigh
 use kernel,   only:radkern
 use units, only:umass,udist
 real, intent(in)     :: primary(3), ray(3), Rstar, xyzh(:,:), kappa(:)
 real, optional       :: maxDistance
 real, intent(out)    :: dist_along_ray(:), tau_along_ray(:)
 integer, intent(out) :: len

 integer, parameter :: maxcache = 0
 real, allocatable  :: xyzcache(:,:)
 real :: distance, h, dtaudr, previousdtaudr, nextdtaudr
 integer :: nneigh, inext, i

 distance = Rstar

 h = Rstar/100.
 inext=0
 do while (inext==0)
    h = h*2.
    call getneigh_pos(primary+Rstar*ray,0.,h,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
    call find_next(primary, ray, distance, xyzh, listneigh, inext, nneigh)
 enddo
 call calc_opacity(primary+Rstar*ray, xyzh, kappa, listneigh, nneigh, previousdtaudr)

 i = 1
 tau_along_ray(i)  = 0.
 distance          = Rstar
 dist_along_ray(i) = distance
 do while (hasNext(inext,tau_along_ray(i),distance,maxDistance))
    i = i + 1
    call getneigh_pos(primary + distance*ray,0.,xyzh(4,inext)*radkern, &
                              3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
    call calc_opacity(primary + distance*ray, xyzh, kappa, listneigh, nneigh, nextdtaudr)
    dtaudr            = (nextdtaudr+previousdtaudr)/2
    previousdtaudr    = nextdtaudr
    tau_along_ray(i)  = tau_along_ray(i-1)+(distance-dist_along_ray(i-1))*dtaudr
    dist_along_ray(i) = distance
    call find_next(primary, ray, distance, xyzh, listneigh, inext,nneigh)
 enddo
 len = i
 tau_along_ray = tau_along_ray*umass/(udist**2)
end subroutine ray_tracer

logical function hasNext(inext, tau, distance, maxDistance)
 integer, intent(in) :: inext
 real, intent(in)    :: distance, tau
 real, optional      :: maxDistance
 real, parameter :: tau_max = 99.
 if (present(maxDistance)) then
    hasNext = inext /= 0 .and. distance < maxDistance .and. tau < tau_max
 else
    hasNext = inext /= 0 .and. tau < tau_max
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
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: neighbors:       A list containing the indices of the neighbors of
 !                       each particle
 !  IN: kappa:           The array containing the opacity of all the SPH particles
 !  IN: Rstar:           The radius of the primary star
 !+
 !  OUT: tau:            The array of optical depths for each SPH particle
 !+
 !  OPT: companion:      The location of the companion
 !  OPT: R:              The radius of the companion
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_inwards(npart, primary, xyzh, neighbors, kappa, Rstar, tau, companion, R)
 real, intent(in)    :: primary(3), kappa(:), Rstar, xyzh(:,:)
 integer, intent(in) :: npart, neighbors(:,:)
 real, optional      :: R, companion(3)
 real, intent(out)   :: tau(:)

 if (present(companion) .and. present(R)) then
    call get_all_tau_inwards_companion(npart, primary, xyzh, neighbors, kappa, Rstar, companion, R, tau)
 else
    call get_all_tau_inwards_single(npart, primary, xyzh, neighbors, kappa, Rstar, tau)
 endif
end subroutine get_all_tau_inwards

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth of each particle, using the inwards ray-
 !  tracing scheme concerning only a single star
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: neighbors:       A list containing the indices of the neighbors of
 !                       each particle
 !  IN: kappa:           The array containing the opacity of all the SPH particles
 !  IN: Rstar:           The radius of the primary star
 !+
 !  OUT: taus:           The list of optical depths for each particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_inwards_single(npart, primary, xyzh, neighbors, kappa, Rstar, tau)
 real, intent(in)    :: primary(3), kappa(:), Rstar, xyzh(:,:)
 integer, intent(in) :: npart, neighbors(:,:)
 real, intent(out)   :: tau(:)

 integer :: i

 !$omp parallel do
 do i = 1, npart
    call get_tau_inwards(i, primary, xyzh, neighbors, kappa, Rstar, tau(i))
 enddo
 !$omp end parallel do
end subroutine get_all_tau_inwards_single

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth of each particle, using the inwards ray-
 !  tracing scheme concerning a binary system
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: neighbors:       A list containing the indices of the neighbors of
 !                       each particle
 !  IN: kappa:           The array containing the opacity of all the SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: companion:       The xyz coordinates of the companion
 !  IN: Rcomp:           The radius of the companion
 !+
 !  OUT: tau:            The array of optical depths for each SPH particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_tau_inwards_companion(npart, primary, xyzh, neighbors, kappa, Rstar, companion, Rcomp, tau)
 real, intent(in)    :: primary(3), companion(3), kappa(:), Rstar, xyzh(:,:), Rcomp
 integer, intent(in) :: npart, neighbors(:,:)
 real, intent(out)   :: tau(:)

 integer :: i
 real    :: normCompanion, theta0, uvecCompanion(3), norm, theta, root, norm0

 uvecCompanion = companion-primary
 normCompanion = norm2(uvecCompanion)
 uvecCompanion = uvecCompanion/normCompanion
 theta0        = asin(Rcomp/normCompanion)

 !$omp parallel do private(norm,theta,root,norm0)
 do i = 1, npart
    norm  = norm2(xyzh(1:3,i)-primary)
    theta = acos(dot_product(uvecCompanion, xyzh(1:3,i)-primary)/norm)
    if (theta < theta0) then
       root = sqrt(normCompanion**2*cos(theta)**2-normCompanion**2+Rcomp**2)
       norm0 = normCompanion*cos(theta)-root
       if (norm > norm0) then
          tau(i) = 99.
       else
          call get_tau_inwards(i, primary, xyzh, neighbors, kappa, Rstar, tau(i))
       endif
    else
       call get_tau_inwards(i, primary, xyzh, neighbors, kappa, Rstar, tau(i))
    endif
 enddo
 !$omp end parallel do
end subroutine get_all_tau_inwards_companion

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth for a given particle, using the inwards ray-
 !  tracing scheme
 !+
 !  IN: point:           The index of the point that needs to be calculated
 !  IN: primary:         The location of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: neighbors:       A list containing the indices of the neighbors of
 !                       each particle
 !  IN: kappa:           The array containing the opacity of all the SPH particles
 !  IN: Rstar:           The radius of the star
 !+
 !  OUT: tau:           The list of optical depth of the given particle
 !+
 !--------------------------------------------------------------------------
subroutine get_tau_inwards(point, primary, xyzh, neighbors, kappa, Rstar, tau)
 use linklist, only:getneigh_pos,ifirstincell,listneigh
 use kernel,   only:radkern
 use units, only:umass,udist
 real, intent(in)    :: primary(3), xyzh(:,:), kappa(:), Rstar
 integer, intent(in) :: point, neighbors(:,:)
 real, intent(out)   :: tau

 integer :: i, next, previous, nneigh
 integer, parameter :: nmaxcache = 0
 real  :: xyzcache(0,nmaxcache)
 real    :: ray(3), nextDist, previousDist, maxDist, dtaudr, previousdtaudr, nextdtaudr

 ray = primary - xyzh(1:3,point)
 maxDist = norm2(ray)
 ray = ray / maxDist
 maxDist=max(maxDist-Rstar,0.)
 next=point
 call getneigh_pos(xyzh(1:3,point),0.,xyzh(4,point)*radkern, &
                           3,listneigh,nneigh,xyzh,xyzcache,nmaxcache,ifirstincell)
 call calc_opacity(xyzh(1:3,point), xyzh, kappa, listneigh, nneigh, nextdtaudr)
 nextDist=0.

 tau = 0.
 i=1
 do while (nextDist < maxDist .and. next /=0)
    i = i + 1
    previous = next
    previousDist = nextDist
    call find_next(xyzh(1:3,point), ray, nextDist, xyzh, neighbors(next,:), next)
    if (nextDist  >  maxDist) then
       nextDist = maxDist
    endif
    call getneigh_pos(xyzh(1:3,point) + nextDist*ray,0.,xyzh(4,previous)*radkern, &
                              3,listneigh,nneigh,xyzh,xyzcache,nmaxcache,ifirstincell)
    previousdtaudr=nextdtaudr
    call calc_opacity(xyzh(1:3,point) + nextDist*ray, xyzh, kappa, listneigh, nneigh, nextdtaudr)
    dtaudr = (nextdtaudr+previousdtaudr)/2
    tau = tau + (nextDist-previousDist)*dtaudr
 enddo
 !fix units for tau (kappa is in cgs while rho & r are in code units)
 tau = tau*umass/(udist**2)
end subroutine get_tau_inwards

 !*********************************************************************!
 !****************************   COMMON   *****************************!
 !*********************************************************************!

 !--------------------------------------------------------------------------
 !+
 !  Find the next point on a ray
 !+
 !  IN: inpoint:         The coordinate of the initial point projected on the
 !                       ray for which the next point will be calculated
 !  IN: ray:             The unit vector of the direction in which the next
 !                       point will be calculated
 !  IN: xyzh:            The array containing the particles position+smoothing length
 !  IN: neighbors:       A list containing the indices of the neighbors of
 !                       the initial point
 !  IN: inext:           The index of the initial point
 !                       (this point will not be considered as possible next point)
 !+
 !  OPT: nneighin:       The amount of neighbors
 !+
 !  OUT: inext:          The index of the next point on the ray
 !+
 !--------------------------------------------------------------------------
subroutine find_next(inpoint, ray, dist, xyzh, neighbors, inext, nneighin)
 integer, intent(in)    :: neighbors(:)
 real, intent(in)       :: xyzh(:,:), inpoint(:), ray(:)
 integer, intent(inout) :: inext
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

 prev=inext
 inext=0
 nextDist=dist
 trace_point = inpoint + dist*ray

 i = 1
 do while (i <= nneigh .and. neighbors(i) /= 0)
    if (neighbors(i)  /=  prev) then
       vec=xyzh(1:3,neighbors(i)) - trace_point
       tempdist = dot_product(vec,ray)
       if (tempdist>0.) then
          raydist = dot_product(vec,vec) - tempdist**2
          if (raydist < dmin) then
             dmin = raydist
             inext = neighbors(i)
             nextdist = dist+tempdist
          endif
       endif
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
 !  OUT: dtaudr:         The local optical depth derivative at the given location (inpoint)
 !+
 !--------------------------------------------------------------------------
subroutine calc_opacity(r0, xyzh, opacities, listneigh, nneigh, dtaudr)
 use kernel,   only:cnormk,wkern
 use part,     only:hfact,rhoh,massoftype,igas
 use dim,      only:maxpsph
 real, intent(in)    :: r0(:), xyzh(:,:), opacities(:)
 integer, intent(in) :: listneigh(:), nneigh
 real, intent(out)   :: dtaudr

 integer :: i,j
 real    :: q

 dtaudr=0
 do i=1,nneigh
    j = listneigh(i)
    if (j > maxpsph) cycle
    q = norm2(r0 - xyzh(1:3,j))/xyzh(4,j)
    dtaudr=dtaudr+wkern(q*q,q)*opacities(j)*rhoh(xyzh(4,j), massoftype(igas))
 enddo
 dtaudr = dtaudr*cnormk/hfact**3
end subroutine calc_opacity
end module raytracer_all
