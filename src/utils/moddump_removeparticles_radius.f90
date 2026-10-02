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
! :Owner: Josh Calcino
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
! :Dependencies: infile_utils, moddump_utils, part, prompting
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 logical :: icutinside  = .false.
 logical :: icutoutside = .false.
 real    :: inradius    = 10.
 real    :: outradius   = 200.
 real    :: incenter(3)  = 0.
 real    :: outcenter(3) = 0.

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:delete_particles_outside_sphere

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: np
 !
 !--prompt only when the driver did not find a parameter file
 !
 if (prompt_for_params) call read_interactive_moddumpfile()

 np = npart

 if (icutinside) then
    print*,'Phantommoddump: Remove particles inside a particular radius'
    print*,'Removing particles inside radius ',inradius
    call delete_particles_outside_sphere(incenter,inradius,np)
 endif

 if (icutoutside) then
    print*,'Phantommoddump: Remove particles outside a particular radius'
    print*,'Removing particles outside radius ',outradius
    call delete_particles_outside_sphere(outcenter,outradius,np)
 endif

end subroutine modify_dump

!
!---Interactively set the moddump parameters--------------------------------
!
subroutine read_interactive_moddumpfile()
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

end subroutine read_interactive_moddumpfile

!
!---Write the moddump parameter file----------------------------------------
!
subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
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
 close(iunit)

end subroutine write_moddump

!
!---Read the moddump parameter file-----------------------------------------
!
subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
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
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

end module moddump
