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
! :Owner: Daniel Mentiplay
!
! :Runtime parameters:
!   - icutinside  : *delete particles inside a given radius*
!   - icutoutside : *delete particles outside a given radius*
!   - incenterx   : *x coordinate of the centre of the inner sphere*
!   - incentery   : *y coordinate of the centre of the inner sphere*
!   - incenterz   : *z coordinate of the centre of the inner sphere*
!   - inradius    : *inward radius [au]*
!   - outcenterx  : *x coordinate of the centre of the outer sphere*
!   - outcentery  : *y coordinate of the centre of the outer sphere*
!   - outcenterz  : *z coordinate of the centre of the outer sphere*
!   - outradius   : *outward radius [au]*
!
! :Dependencies: moddump_utils, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 logical :: icutinside  = .false.
 logical :: icutoutside = .false.
 real    :: inradius    = 10.
 real    :: outradius   = 200.
 real    :: incenter(3)  = 0.
 real    :: outcenter(3) = 0.

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:delete_particles_outside_sphere
 use moddump_utils, only:get_moddump_options

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 !
 !--read parameters from the .moddump file, prompting if it does not exist
 !
 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)


 if (icutinside) then
    print*,'Phantommoddump: Remove particles inside a particular radius'
    print*,'Removing particles inside radius ',inradius
    call delete_particles_outside_sphere(incenter,inradius,npart,revert=.true.)
 endif

 if (icutoutside) then
    print*,'Phantommoddump: Remove particles outside a particular radius'
    print*,'Removing particles outside radius ',outradius
    call delete_particles_outside_sphere(outcenter,outradius,npart)
 endif

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Deleting particles inside a given radius ?',icutinside)
 call prompt('Deleting particles outside a given radius ?',icutoutside)
 if (icutinside) then
    call prompt('Enter inward radius in au',inradius,0.)
    call prompt('Enter x coordinate of the center of that sphere',incenter(1))
    call prompt('Enter y coordinate of the center of that sphere',incenter(2))
    call prompt('Enter z coordinate of the center of that sphere',incenter(3))
 endif
 if (icutoutside) then
    call prompt('Enter outward radius in au',outradius,0.)
    call prompt('Enter x coordinate of the center of that sphere',outcenter(1))
    call prompt('Enter y coordinate of the center of that sphere',outcenter(2))
    call prompt('Enter z coordinate of the center of that sphere',outcenter(3))
 endif

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(icutinside,'icutinside','delete particles inside a given radius',iunit)
 call write_inopt(icutoutside,'icutoutside','delete particles outside a given radius',iunit)
 call write_inopt(inradius,'inradius','inward radius [au]',iunit)
 call write_inopt(incenter(1),'incenterx','x coordinate of the centre of the inner sphere',iunit)
 call write_inopt(incenter(2),'incentery','y coordinate of the centre of the inner sphere',iunit)
 call write_inopt(incenter(3),'incenterz','z coordinate of the centre of the inner sphere',iunit)
 call write_inopt(outradius,'outradius','outward radius [au]',iunit)
 call write_inopt(outcenter(1),'outcenterx','x coordinate of the centre of the outer sphere',iunit)
 call write_inopt(outcenter(2),'outcentery','y coordinate of the centre of the outer sphere',iunit)
 call write_inopt(outcenter(3),'outcenterz','z coordinate of the centre of the outer sphere',iunit)

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

 call read_inopt(icutinside,'icutinside',db,errcount=nerr)
 call read_inopt(icutoutside,'icutoutside',db,errcount=nerr)
 call read_inopt(inradius,'inradius',db,min=0.,errcount=nerr)
 call read_inopt(incenter(1),'incenterx',db,errcount=nerr)
 call read_inopt(incenter(2),'incentery',db,errcount=nerr)
 call read_inopt(incenter(3),'incenterz',db,errcount=nerr)
 call read_inopt(outradius,'outradius',db,min=0.,errcount=nerr)
 call read_inopt(outcenter(1),'outcenterx',db,errcount=nerr)
 call read_inopt(outcenter(2),'outcentery',db,errcount=nerr)
 call read_inopt(outcenter(3),'outcenterz',db,errcount=nerr)

end subroutine read_options_moddump

end module moddump
