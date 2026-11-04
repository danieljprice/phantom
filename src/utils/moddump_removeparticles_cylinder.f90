!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Remove particles outside a cylinder
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - radius : *cylinder radius [code units]*
!   - xcen   : *cylinder centre x [code units]*
!   - ycen   : *cylinder centre y [code units]*
!   - zcen   : *cylinder centre z [code units]*
!   - zmax   : *cylinder half-height [code units]*
!
! :Dependencies: moddump_utils, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters (written to / read from the prefix.moddump file)
 real :: xcen   = 0.      ! cylinder centre x [code units]
 real :: ycen   = 0.      ! cylinder centre y [code units]
 real :: zcen   = 0.      ! cylinder centre z [code units]
 real :: radius = 1500.   ! cylinder radius [code units]
 real :: zmax   = 1500.   ! cylinder half-height [code units]

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:delete_particles_outside_cylinder
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real    :: center(3)

 print*,' Phantommoddump: Remove particles outside a cylinder'

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 center = (/xcen,ycen,zcen/)
 !
 !--removing particles
 !
 print*,'Removing particles outside the cylinder centered in ( ', center,' ), with radius ',radius,' and zmax ',zmax,' : '
 call delete_particles_outside_cylinder(center, radius, zmax, npartoftype)

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter cylinder centre x (code units)',xcen)
 call prompt('Enter cylinder centre y (code units)',ycen)
 call prompt('Enter cylinder centre z (code units)',zcen)
 call prompt('Enter cylinder radius (code units)',radius,0.)
 call prompt('Enter cylinder half-height zmax (code units)',zmax,0.)

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

 call read_inopt(xcen,'xcen',db,errcount=nerr)
 call read_inopt(ycen,'ycen',db,errcount=nerr)
 call read_inopt(zcen,'zcen',db,errcount=nerr)
 call read_inopt(radius,'radius',db,errcount=nerr,min=0.)
 call read_inopt(zmax,'zmax',db,errcount=nerr,min=0.)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# remove-outside-cylinder parameters'
 call write_inopt(xcen,'xcen','cylinder centre x [code units]',iunit)
 call write_inopt(ycen,'ycen','cylinder centre y [code units]',iunit)
 call write_inopt(zcen,'zcen','cylinder centre z [code units]',iunit)
 call write_inopt(radius,'radius','cylinder radius [code units]',iunit)
 call write_inopt(zmax,'zmax','cylinder half-height [code units]',iunit)

end subroutine write_options_moddump

end module moddump
