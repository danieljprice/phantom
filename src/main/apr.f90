!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module apr
  !
  ! Contains everything for live adaptive particle refinement
  !
  ! :References: None
  !
  ! :Owner: Rebecca Nealon
  !
  ! :Runtime parameters: None
  !
  ! :Dependencies: None
  !
  implicit none

contains

  !--------------------------------------------------------------------
  ! Subroutine to initialise arrays
  !--------------------------------------------------------------------
  subroutine init_apr(apr_level,ierr)
    integer, intent(inout) :: ierr,apr_level(:)

    apr_level = 3
    ierr = 0

  end subroutine init_apr

  !--------------------------------------------------------------------
  ! Subroutine to check if particles need to be split or merged
  !--------------------------------------------------------------------
  subroutine update_apr(npart,xyzh,vxyzu,apr_level)
    real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
    integer, intent(inout) :: npart,apr_level(:)
    integer :: ii

    ! Do any particles need to be split?


    ! Do any particles need to be merged?

  end subroutine update_apr

end module apr
