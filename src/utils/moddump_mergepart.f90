!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! merges particles; input simulation with npart, get npart/nchild back
!
! :References: Vacondio et al. 2013
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters:
!   - nchild : *number of children per merged particle (>= 2)*
!
! :Dependencies: moddump_utils, part, prompting, splitpart
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 integer :: nchild = 2   ! number of children merged into each particle (>= 2)

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,     only:merge_all_particles
 use part,          only:kill_particle,delete_dead_or_accreted_particles
 use part,          only:isdead_or_accreted
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,nactive

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)
 if (nchild < 2) stop 'error nchild cannot be < 2'

 !-- how many active particles
 nactive = 0
 do i = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       nactive = nactive + 1
    else
       call kill_particle(i,npartoftype)
    endif
 enddo

 if (nactive < npart) then
    call delete_dead_or_accreted_particles(npart,npartoftype)
    print*,' discarding inactive particles'
 endif

 ! Merge 'em!
 call merge_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                          nchild,nactive)

 print*,' new npart = ',npart

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter number of children per merged particle',nchild,2)

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

 call read_inopt(nchild,'nchild',db,errcount=nerr,min=2)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# merge parameters'
 call write_inopt(nchild,'nchild','number of children per merged particle (>= 2)',iunit)

end subroutine write_options_moddump

end module moddump
