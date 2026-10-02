!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! split every particle in a dump into nchild children
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - lattice_type : *child arrangement (0=regular lattice, 1=random)*
!   - nchild       : *number of children per particle (>= 2; forced to 13 for lattice)*
!
! :Dependencies: infile_utils, moddump_utils, part, prompting, splitpart
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters
 integer :: nchild       = 13   ! number of children per particle
 integer :: lattice_type = 0    ! 0 for regular lattice, 1 for random

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,    only:split_all_particles
 use part,         only:delete_dead_or_accreted_particles
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 if (prompt_for_params) call read_interactive_moddumpfile()

 !-- the regular lattice requires a specific number of children
 if (lattice_type == 0) nchild = 13
 if (nchild < 2) stop 'error nchild cannot be < 2'

 !-- don't split accreted particles
 call delete_dead_or_accreted_particles(npart,npartoftype)

 ! Split 'em!
 print "(/,a,i0,a)", ' >>> splitting all particles into ',nchild,' children <<<'

 if (lattice_type==0) then
    print "(a,/)", ' >>> placing children on regular lattice <<<'
 else
    print "(a,/)", ' >>> placing children using random arrangement <<<'
 endif
 call split_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                          nchild,lattice_type,1)

 print "(a,i0,/)",' new npart = ',npart

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter child arrangement (0=regular lattice, 1=random)',lattice_type,0,1)
 if (lattice_type /= 0) call prompt('Enter number of children per particle',nchild,2)

end subroutine read_interactive_moddumpfile

subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 23
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(lattice_type,'lattice_type',db,errcount=nerr,min=0,max=1)
 call read_inopt(nchild,'nchild',db,errcount=nerr,min=2)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# split parameters'
 call write_inopt(lattice_type,'lattice_type','child arrangement (0=regular lattice, 1=random)',iunit)
 call write_inopt(nchild,'nchild','number of children per particle (>= 2; forced to 13 for lattice)',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
