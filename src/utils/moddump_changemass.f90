!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Changes particle mass
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - disc_mass : *desired total disc mass [code units]*
!
! :Dependencies: moddump_utils, part, prompting, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter (written to / read from the prefix.moddump file)
 real :: disc_mass = 0.05   ! desired total disc mass [code units]

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:igas,isdead_or_accreted,kill_particle,shuffle_part
 use units,         only:umass
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real     :: current_disc_mass, mass_factor
 integer  :: i

 ! Remove particles that are dead or accreted
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       call kill_particle(i)
    endif
 enddo

 call shuffle_part(npart)
 npartoftype(igas) = npart

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 current_disc_mass = npartoftype(igas)*massoftype(igas)
 mass_factor = disc_mass/current_disc_mass

 massoftype(igas) = mass_factor*massoftype(igas)
 print*,'Particle mass is now ', massoftype(igas)*umass, ' g'
 print*,'Total disc mass is now ', npartoftype(igas)*massoftype(igas)*umass, ' g'
 print*,'Total disc mass is now ', npartoftype(igas)*massoftype(igas), 'Msun'

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter desired total disc mass in code units',disc_mass,0.)

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

 call read_inopt(disc_mass,'disc_mass',db,errcount=nerr,min=0.)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# changemass parameters'
 call write_inopt(disc_mass,'disc_mass','desired total disc mass [code units]',iunit)

end subroutine write_options_moddump

end module moddump
