!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Give velocity perturbation to gas particles
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: part, systemutils
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = '--perturb_factor=0.5 --perturb_sink=.false. --sink_perturb_factor=0.5'

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:vxyz_ptmass,nptmass
 use systemutils,       only:get_command_option_real,get_command_option_logical
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: perturb_factor,sink_perturb_factor
 integer                :: i
 logical                :: perturb_sink

 perturb_sink = get_command_option_logical('perturb_sink',default=.false.)
 perturb_factor = get_command_option_real('perturb_factor',default=0.5)
 sink_perturb_factor = get_command_option_real('sink_perturb_factor',default=0.5)

 print "(a,g0)", '>>> perturbing gas with factor ',perturb_factor
 do i=1,npart
    vxyzu(1:3,i) = (1. + perturb_factor) * vxyzu(1:3,i)
 enddo

 if (perturb_sink) then
    print "(a,g0)", '>>> perturbing sink particle with factor ',sink_perturb_factor
    do i=1,nptmass
       vxyz_ptmass(1:3,i) = (1. + sink_perturb_factor) * vxyz_ptmass(1:3,i)
    enddo
 endif

end subroutine modify_dump

end module moddump

