!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module splitpart
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: icosahedron, injectutils, io, part, random
!
 implicit none

contains

!--------------------------------------------------------------------------
!+
!  split particles into nchild particle
!+
!--------------------------------------------------------------------------
subroutine split_particles(nchild,iparent,xyzh_parent,vxyzu_parent, &
           lattice_type,ires,xyzh_child,vxyzu_child)
 use icosahedron, only:pixel2vector,compute_corners,compute_matrices
 integer, intent(in)    :: nchild
 integer, intent(in)    :: iparent,lattice_type,ires
 real,    intent(in)    :: xyzh_parent(:),vxyzu_parent(:)
 real,    intent(out)   :: xyzh_child(4,nchild),vxyzu_child(3,nchild)
 integer :: j,iseed
 real    :: dhfac,dx(3),sep,geodesic_R(0:19,3,3), geodesic_v(0:11,3)

 if (lattice_type == 0) then
    call compute_matrices(geodesic_R)
    call compute_corners(geodesic_v)
 else
    ! initialise random number generator
    iseed = -6542
 endif

 ! first child sits where the parent was
 ! for now, assume all children have same velocity as parent
 xyzh_child(1:3,1) = xyzh_parent(1:3)
 vxyzu_child(:,1)  = vxyzu_parent(1:3)
 sep = xyzh_parent(4)

 do j=2,nchild
   ! find positions of original particles
   if (lattice_type == 0) then
     call pixel2vector(j-2,ires,geodesic_R,geodesic_v,dx)
   else
     call sample_kernel(iseed,dx)
   endif
      xyzh_child(1:3,j) = xyzh_parent(1:3) + sep*dx(:)
      vxyzu_child(1:3,j) = vxyzu_parent(1:3)
 enddo

 !--amend smoothing length
 dhfac = 1./(nchild)**(1./3.)
 xyzh_child(4,:) = xyzh_parent(4)*dhfac

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

subroutine shift_particles(npart,xyzh,vxyzu,deltat,beta,shifts)
  use kernel, only:radkern2
  real, intent(in)    :: xyzh(:,:),vxyzu(:,:)
  real, intent(in)    :: deltat,beta
  integer, intent(in) :: npart
  real, intent(out)   :: shifts(3,npart)
  integer             :: i,j,neighbours
  real                :: rnaught,rij2,dr3,vel2,vmax
  real                :: q2,rij(3),rsum(3)

  vmax = tiny(vmax)
  vel2 = 0.

  do i = 1,npart
    rnaught = 0.
    neighbours = 0
    rsum = 0.
    vel2 = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
    if (vel2 > vmax) vmax = vel2

    over_npart: do j = 1,npart
      if (i == j) cycle over_npart
      rij = xyzh(1:3,j) - xyzh(1:3,i)
      rij2 = dot_product(rij,rij)
      q2 = rij2/(xyzh(4,i)*xyzh(4,i))

      if (q2 < radkern2) then
        neighbours = neighbours + 1
        rnaught = rnaught + sqrt(rij2)
      endif

      dr3 = 1./(rij2**1.5)
      rsum = rsum + (rij*dr3)
    enddo over_npart

    rnaught = rnaught/neighbours
    shifts(:,i) = beta*rnaught*rnaught*deltat*rsum
  enddo

  shifts = shifts*sqrt(vel2)

end subroutine shift_particles

!-----------------------------------------------------------------------
!+
! merges nchild particles in one parent particle
!+
!-----------------------------------------------------------------------
subroutine merge_particles(nchild,children_list,mchild,npart, &
           xyzh,vxyzu,xyzh_parent,vxyzu_parent)
 use kernel, only:get_kernel,cnormk,radkern
 integer, intent(in)    :: nchild,children_list(nchild)
 integer, intent(inout) :: npart
 real,    intent(in)    :: xyzh(:,:),vxyzu(:,:),mchild
 real,    intent(out)   :: xyzh_parent(4),vxyzu_parent(3)
 integer :: i,j,ichild
 real    :: h1,h31,rij_vec(3)
 real    :: qij,rij,wchild,grkernchild,rho_parent

 !-- positions and velocities from centre of mass
 xyzh_parent = 0.
 vxyzu_parent = 0.
 do i=1,nchild
   ichild = children_list(i)
   xyzh_parent(1:3) = xyzh_parent(1:3) + xyzh(1:3,ichild)
   vxyzu_parent(:)  = vxyzu_parent(:) + vxyzu(1:3,ichild)
 enddo
 xyzh_parent = xyzh_parent/nchild
 vxyzu_parent = vxyzu_parent/nchild

!-- calculate density at the parent position from the children
!-- procedure described in Vacondio et al. 2013, around Eq 21
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
xyzh_parent(4) = 1.2*(nchild*mchild/rho_parent)**(1./3.)

end subroutine merge_particles

end module splitpart
