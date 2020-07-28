!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: splitpart
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: icosahedron, injectutils, io, part, random
!+
!--------------------------------------------------------------------------
module splitpart
 implicit none
 integer, parameter :: ires = 1 ! use 12 particles per sphere

contains

!--------------------------------------------------------------------------
!+
!  split particles into nchild particle
!+
!--------------------------------------------------------------------------
subroutine split_particles(nchild,npart,npartoftype,xyzh,vxyzu,massoftype,ierr)
 use injectutils, only:get_parts_per_sphere
 use icosahedron, only:pixel2vector,compute_corners,compute_matrices
 use part,        only:copy_particle
 use io,          only:error
 integer, intent(in)    :: nchild
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(inout) :: massoftype(:)
 integer, intent(out)   :: ierr
 real, allocatable, dimension(:,:) :: shifts
 integer :: npart_per_sphere,iold,inew,j,iseed
 real    :: dhfac,dx(3),sep,vel,vmax,dv(3)
 real    :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)

 ierr = 0
 if (nchild <= 0) then
    npart_per_sphere = get_parts_per_sphere(ires)
    call compute_matrices(geodesic_R)
    call compute_corners(geodesic_v)
 else
    ! initialise random number generator
    npart_per_sphere = nchild - 1
    iseed = -6542
 endif

 !--check there is enough memory
 if (size(xyzh(1,:)) < npart*(npart_per_sphere+1)) then
    call error('split_particles','not enough memory, increase MAXP and recompile')
    ierr = 1
    return
 endif

 !--divide mass of all particles by 13
 massoftype(:) = massoftype(:) / (npart_per_sphere + 1)

 !--update npartoftype
 npartoftype(:) = npartoftype*(npart_per_sphere+1)

 !--amend smoothing length of original particles
 !--also find the maximum velocity
 vmax = -1.0
 dhfac = 1./(npart_per_sphere + 1)**(1./3.)
 do iold=1,npart
    xyzh(4,iold) = xyzh(4,iold)*dhfac
    dv = vxyzu(1:3,iold)
    vel = sqrt(dot_product(dv,dv))
    if (vel > vmax) vmax = vel
 enddo

 !--find positions of the new particles
 inew = npart ! put new particles after original parts
 do iold=1,npart
    sep = xyzh(4,iold)
    ! add 12 new child particles
    do j=0,npart_per_sphere-1
       inew = inew + 1
       ! copy all properties from original particle
       call copy_particle(iold,inew)
       ! find positions of original particles
       if (nchild <= 0) then
          call pixel2vector(j,ires,geodesic_R,geodesic_v,dx)
       else
          call sample_kernel(iseed,dx)
       endif
       xyzh(1:3,inew) = xyzh(1:3,iold) + sep*dx(:)
    enddo
    ! amend smoothing length of original particle
    npart = npart + npart_per_sphere
 enddo

!-- now shift all the particles a little
allocate(shifts(3,npart))
call shift_particles(npart,xyzh,vmax,0.0001,100.0,shifts)
xyzh(1:3,1:npart) = xyzh(1:3,1:npart) + shifts(:,:)
deallocate(shifts)

end subroutine split_particles

subroutine sample_kernel(iseed,dx)
 use random, only:gauss_random,get_random_pos_on_sphere
 integer, intent(inout) :: iseed
 real :: dx(3),r

 r = 3.
 do while (abs(r) > 2.)
    r = gauss_random(iseed)
 enddo
 dx = get_random_pos_on_sphere(iseed)
 dx = r*dx

end subroutine sample_kernel

subroutine shift_particles(npart,xyzh,vmax,deltat,beta,shifts)
  real, intent(in)    :: xyzh(:,:)
  real, intent(in)    :: deltat,vmax,beta
  integer, intent(in) :: npart
  real, intent(out)   :: shifts(3,npart)
  integer             :: i,j,neighbours
  real                :: rnaught,rij2,dr3
  real                :: q2,rij(3),rsum(3)

  do i = 1,npart
    rnaught = 0.
    neighbours = 0
    rsum = 0.

    over_npart: do j = 1,npart
      if (i == j) cycle over_npart
      rij = xyzh(1:3,j) - xyzh(1:3,i)
      rij2 = dot_product(rij,rij)
      q2 = rij2/(xyzh(4,i)*xyzh(4,i))

      if (q2 < 4.0) then !assuming defaults for testing, must be updated
        neighbours = neighbours + 1
        rnaught = rnaught + sqrt(rij2)
      endif

      dr3 = 1./(rij2**1.5)
      rsum = rsum + (rij*dr3)
    enddo over_npart

    rnaught = rnaught/neighbours
    shifts(:,i) = beta*rnaught*rnaught*vmax*deltat*rsum
  enddo

end subroutine shift_particles

!-----------------------------------------------------------------------
!+
! merges nchild particles in one parent particle
!+
!-----------------------------------------------------------------------
subroutine merge_particles(nchild,children_list,mchild,npart, &
           xyzh,vxyzu,xyzh_parent,vxyzu_parent,mparent)
 use kernel,      only:get_kernel,cnormk,radkern
 integer, intent(in)    :: nchild,children_list(nchild)
 integer, intent(inout) :: npart
 real,    intent(in)    :: xyzh(:,:),vxyzu(:,:),mchild
 real,    intent(out)   :: xyzh_parent(4),vxyzu_parent(3),mparent
 integer :: i,j,ichild
 real    :: h1,h31,rij_vec(3)
 real    :: qij,rij,wchild,grkernchild,rho_parent

 !--set mass
 mparent = mchild*nchild

 !-- positions and velocities from centre of mass
 xyzh_parent = 0.
 vxyzu_parent = 0.
 do i=1,nchild
   ichild = children_list(i)
   xyzh_parent(1:3) = xyzh_parent(1:3) + mchild*xyzh(1:3,ichild)
   vxyzu_parent(:)  = vxyzu_parent(:) + mchild*vxyzu(1:3,ichild)
 enddo
 xyzh_parent = xyzh_parent/mparent
 vxyzu_parent = vxyzu_parent/mparent

!-- calculate density at the parent position from the children
rho_parent = 0.
  over_npart:  do j=1,npart
    if (xyzh(4,j) < 0.) cycle over_npart
    h1 = 1./xyzh(4,j)
    h31 = h1**3

    rij_vec = xyzh_parent(1:3) - xyzh(1:3,j)

    rij = sqrt(dot_product(rij_vec,rij_vec))
    qij = rij*h1

    wchild = 0.
    if (qij < radkern) call get_kernel(qij*qij,qij,wchild,grkernchild)

    rho_parent = rho_parent + (mchild*wchild*cnormk*h31)
  enddo over_npart

!-- smoothing length from density
xyzh_parent(4) = 1.2*(mparent/rho_parent)**(1./3.)

end subroutine merge_particles

end module splitpart
