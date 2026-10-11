!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! None
!
! :References: None
!
! :Owner: Antoine Alaguero
!
! :Runtime parameters:
!   - rho_threshold : *delete particles with density below this [code units]*
!
! :Dependencies: io, moddump_utils, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter (written to / read from the prefix.moddump file)
 real :: rho_threshold = 5e-8   ! delete particles with density below this [code units]

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:rho,kill_particle,shuffle_part
 use io,            only:fatal
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, dimension(:), intent(inout) :: npartoftype
 real, dimension(:), intent(inout) :: massoftype
 real, dimension(:,:), intent(inout) :: xyzh,vxyzu
 real   :: rhoi
 integer :: i,compt

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 compt = 0
 do i=1,npart
    rhoi = rho(i)
    ! write(*,*) rhoi      ! uncomment to have an idea of the density
    if (rhoi < rho_threshold) then
       call kill_particle(i,npartoftype)
       compt = compt+1
       !xyzh(4,i) = -abs(hi)    !call kill_particle(i,npoftype)
    endif
 enddo
 write(*,*) 'Particles deleted :', compt
 call shuffle_part(npart)
 if (npart /= sum(npartoftype)) call fatal('del_dead_part_dens','particles not conserved')

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter density threshold below which to delete particles (code units)',rho_threshold,0.)

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

 call read_inopt(rho_threshold,'rho_threshold',db,errcount=nerr,min=0.)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# remove-low-density-particles parameters'
 call write_inopt(rho_threshold,'rho_threshold','delete particles with density below this [code units]',iunit)

end subroutine write_options_moddump

end module moddump
