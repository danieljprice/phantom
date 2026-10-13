!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump_utils
!
! Shared helpers for moddumps. get_moddump_options is a thin wrapper
! around get_options that handles the prefix.moddump parameter file:
! opening/closing it and writing a header recording which dump file was
! modified. Each moddump only supplies routines that call write_inopt
! and read_inopt, in the same way as write_options_X/read_options_X in
! the main code modules. Any option can also be set on the command line
! as --option=value, overriding the value in the .moddump file
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, io, mpiutils
!
 use infile_utils, only:write_inopt,read_inopt,inopts
 implicit none

 ! set by phantommoddump, written into the .moddump file as comments
 character(len=120), public :: moddump_dumpfile_in = ''  ! input dump being modified
 real,               public :: moddump_time = 0.         ! output time of the new dump

 public :: get_moddump_options
 ! re-exported so moddumps only need to use this module for their options
 public :: write_inopt,read_inopt,inopts

 abstract interface
  subroutine read_opts_func(db,nerr)
   import :: inopts
   type(inopts), intent(inout) :: db(:)
   integer,      intent(inout) :: nerr
  end subroutine read_opts_func
 end interface

 abstract interface
  subroutine write_opts_func(iunit)
   integer, intent(in) :: iunit
  end subroutine write_opts_func
 end interface

 abstract interface
  subroutine set_func()
  end subroutine set_func
 end interface

 ! moddump-specific routines, set by get_moddump_options for use in read/write_modfile
 procedure(read_opts_func),  pointer :: read_opts  => null()
 procedure(write_opts_func), pointer :: write_opts => null()

 private

contains

!----------------------------------------------------------------------
!+
!  read/write moddump options from/to prefix.moddump via get_options:
!   - file complete   -> read it and return
!   - file incomplete -> missing options are added, then stop
!   - file absent     -> prompt with set_options (or use any options
!                        given on the command line), write the file,
!                        then stop, as for phantomsetup
!+
!----------------------------------------------------------------------
subroutine get_moddump_options(read_options,write_options,set_options)
 use infile_utils, only:get_options
 use io,           only:fileprefix,id,master
 use mpiutils,     only:barrier_mpi
 procedure(read_opts_func)     :: read_options
 procedure(write_opts_func)    :: write_options
 procedure(set_func), optional :: set_options
 character(len=len(fileprefix)+8) :: filename
 integer :: ierr

 read_opts  => read_options
 write_opts => write_options

 filename = trim(fileprefix)//'.moddump'
 if (present(set_options)) then
    call get_options(filename,id==master,ierr,read_modfile,write_modfile,set_options)
 else
    call get_options(filename,id==master,ierr,read_modfile,write_modfile)
 endif
 if (ierr /= 0) then
    if (id==master) print "(/,a,/)",' >> edit '//trim(filename)//' and rerun phantommoddump'
    stop
 endif
 ! rewrite the file so the header records the dump being modified
 ! (done before modify_dump converts any options to code units, and
 ! once all threads have finished reading it)
 call barrier_mpi()
 if (id==master) call write_modfile(filename)

end subroutine get_moddump_options

!----------------------------------------------------------------------
!+
!  write the .moddump file: header, then the moddump-specific options
!+
!----------------------------------------------------------------------
subroutine write_modfile(filename)
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing moddump options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)")    '# input file for phantommoddump'
 write(iunit,"(2a)")   '# input dump file modified: ',trim(moddump_dumpfile_in)
 write(iunit,"(a,g0)") '# output time of modified dump: ',moddump_time
 call write_opts(iunit)
 close(iunit)

end subroutine write_modfile

!----------------------------------------------------------------
!+
!  read the .moddump file, returning ierr/=0 if any option is
!  missing or invalid
!+
!----------------------------------------------------------------
subroutine read_modfile(filename,ierr)
 use infile_utils, only:open_db_from_file,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 23
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_opts(db,nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i0,a)",nerr,' error(s) during read of moddump options file: re-writing...'
    ierr = nerr
 endif

end subroutine read_modfile

end module moddump_utils
