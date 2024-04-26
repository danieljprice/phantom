!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_sampling
!
! Contains simple routine to sample variable using specific distributions
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 public :: divide_unit_seg

contains

subroutine divide_unit_seg(lengths,mindist,nlengths)
 integer, intent(in)    :: nlengths
 real,    intent(inout) :: lengths(nlengths)
 real,    intent(in)    :: mindist
 integer :: i,j
 logical :: far
 real :: points(nlengths+1),tmp,dist
 points(nlengths+1) = 1.
 points(1)          = 0.
 tmp  = 0.

 do i=2,nlengths
    far =  .false.
    dist = huge(tmp)
    do while (far)
       tmp = rand()
       dist = min(abs(points(1)-tmp),dist)
       dist = min(abs(points(nlengths+1)-tmp),dist)
       do j=2,i-1
          dist = min(abs(points(j)-tmp),dist)
       enddo
       far = mindist<dist
    enddo
    points(i) = tmp
 enddo

 do i=2,nlengths+1
    lengths(i-1) = points(i) - points(i-1)
 enddo



end subroutine divide_unit_seg


end module utils_sampling
