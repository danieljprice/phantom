!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! :Runtime parameters:
!   - perturb_factor      : *fractional velocity perturbation for gas*
!   - perturb_sink        : *also perturb sink particle velocities*
!   - sink_perturb_factor : *fractional velocity perturbation for sinks*
!
! :Dependencies: moddump_utils, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters
 real    :: perturb_factor      = 0.5      ! fractional velocity perturbation for gas
 real    :: sink_perturb_factor = 0.5      ! fractional velocity perturbation for sinks
 logical :: perturb_sink        = .false.  ! also perturb sink particle velocities

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:vxyz_ptmass,nptmass
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

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

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter fractional velocity perturbation for gas',perturb_factor)
 call prompt('Also perturb sink particle velocities?',perturb_sink)
 if (perturb_sink) call prompt('Enter fractional velocity perturbation for sinks',sink_perturb_factor)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(perturb_factor,'perturb_factor',db,errcount=nerr)
 call read_inopt(perturb_sink,'perturb_sink',db,errcount=nerr)
 if (perturb_sink) call read_inopt(sink_perturb_factor,'sink_perturb_factor',db,errcount=nerr)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# gas velocity perturbation parameters'
 call write_inopt(perturb_factor,'perturb_factor','fractional velocity perturbation for gas',iunit)
 call write_inopt(perturb_sink,'perturb_sink','also perturb sink particle velocities',iunit)
 if (perturb_sink) call write_inopt(sink_perturb_factor,'sink_perturb_factor','perturbation for sinks',iunit)

end subroutine write_options_moddump

end module moddump
