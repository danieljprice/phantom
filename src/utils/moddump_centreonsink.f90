!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Center positions and velocities of particles around one given sink
!
! :References: None
!
! :Owner: Antoine Alaguero
!
! :Runtime parameters:
!   - sink_ind : *index of the sink to centre on*
!
! :Dependencies: moddump_utils, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 integer :: sink_ind = 1   ! index of the sink to centre positions/velocities on

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 if (nptmass < sink_ind) then
    print*,'Selected sink index larger than number of sinks'
    return
 endif

 do i=1,npart
    xyzh(1:3,i) = xyzh(1:3,i) - xyzmh_ptmass(1:3,sink_ind)
    vxyzu(1:3,i) = vxyzu(1:3,i) - vxyz_ptmass(1:3,sink_ind)
 enddo

 do i=1,nptmass
    !skip sink_ind, because useless if put to 0 first
    if (i==sink_ind) then
       cycle
    else
       xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - xyzmh_ptmass(1:3,sink_ind)
       vxyz_ptmass(1:3,i) = vxyz_ptmass(1:3,i) - vxyz_ptmass(1:3,sink_ind)
    endif
 enddo

 xyzmh_ptmass(1:3,sink_ind) = xyzmh_ptmass(1:3,sink_ind) - xyzmh_ptmass(1:3,sink_ind)
 vxyz_ptmass(1:3,sink_ind) = vxyz_ptmass(1:3,sink_ind) - vxyz_ptmass(1:3,sink_ind)

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  set parameters interactively (when no .moddump file is found)
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter index of the sink to centre on',sink_ind,1)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# centre-on-sink parameters'
 call write_inopt(sink_ind,'sink_ind','index of the sink to centre on',iunit)

end subroutine write_options_moddump

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(sink_ind,'sink_ind',db,errcount=nerr,min=1)

end subroutine read_options_moddump

end module moddump
