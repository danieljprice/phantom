!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testinfile
!
! Unit tests of the reading and writing of parameter files (e.g. .setup
! and .moddump files) with get_options, including options given on the
! command line as --option=value
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, io, testutils
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_infile

 private

 ! options in the test parameter file
 integer           :: np     = 100
 real              :: omega  = 0.5
 logical           :: lflag  = .false.
 character(len=20) :: label  = 'abc'
 logical           :: prompted = .false.  ! records whether set_pars was called

 character(len=40) :: filename  ! temporary file, one per MPI task

contains

!--------------------------------------------
!+
!  tests of parameter file handling
!+
!--------------------------------------------
subroutine test_infile(ntests,npass)
 integer, intent(inout) :: ntests,npass

 ! run on all tasks (update_test_scores is collective) with one file per task
 write(filename,"(a,i0,a)") 'test_infile_tmp',id,'.setup'
 if (id==master) write(*,"(/,a)") '--> TESTING PARAMETER FILE HANDLING'

 call test_parse_option_string(ntests,npass)
 call test_override_existing_file(ntests,npass)
 call test_override_no_file(ntests,npass)
 call test_prompt_no_file(ntests,npass)
 call test_incomplete_file(ntests,npass)
 call test_unknown_override(ntests,npass)
 call test_string_default(ntests,npass)
 call delete_file(filename)

 if (id==master) write(*,"(/,a)") '<-- PARAMETER FILE TEST COMPLETE'

end subroutine test_infile

!--------------------------------------------
!+
!  parsing of --option=value strings
!+
!--------------------------------------------
subroutine test_parse_option_string(ntests,npass)
 use infile_utils, only:parse_option_string
 integer, intent(inout) :: ntests,npass
 character(len=20) :: name,val
 integer :: nfailed(10)

 if (id==master) write(*,"(/,a)") '--> testing parsing of --option=value'
 nfailed = 0
 call checkval(parse_option_string('--np=1000',name,val),.true.,nfailed(1),'--np=1000 is an option')
 call checkval(name,'np',nfailed(2),'option name')
 call checkval(val,'1000',nfailed(3),'option value')
 call checkval(parse_option_string('-np=5',name,val),.true.,nfailed(4),'-np=5 is an option')
 call checkval(name,'np',nfailed(5),'option name with single dash')
 call checkval(parse_option_string('dump_00000',name,val),.false.,nfailed(6),'positional argument')
 call checkval(parse_option_string('--flag',name,val),.false.,nfailed(7),'flag without =')
 call checkval(parse_option_string('--=x',name,val),.false.,nfailed(8),'option without name')
 call checkval(parse_option_string('np=5',name,val),.false.,nfailed(9),'option without dashes')
 call checkval(parse_option_string('',name,val),.false.,nfailed(10),'empty string')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_parse_option_string

!-------------------------------------------------------------------
!+
!  options on the command line override those in an existing file,
!  and the file is rewritten with the values used
!+
!-------------------------------------------------------------------
subroutine test_override_existing_file(ntests,npass)
 use infile_utils, only:get_options,set_overrides,clear_overrides
 integer, intent(inout) :: ntests,npass
 integer :: ierr,nfailed(5)

 if (id==master) write(*,"(/,a)") '--> testing command line options with an existing file'
 nfailed = 0
 call set_defaults()
 call write_pars(filename)
 call set_overrides([character(len=20) :: '--np=7','--omega=2.5','--maxp=1000'])
 call get_options(filename,.true.,ierr,read_pars,write_pars)
 call clear_overrides()
 call checkval(ierr,0,0,nfailed(1),'ierr')
 call checkval(np,7,0,nfailed(2),'integer option overridden')
 call checkval(omega,2.5,1.e-6,nfailed(3),'real option overridden')

 ! the file should now contain the values used
 call set_defaults()
 call get_options(filename,.true.,ierr,read_pars,write_pars)
 call checkval(np,7,0,nfailed(4),'file rewritten with integer option')
 call checkval(omega,2.5,1.e-6,nfailed(5),'file rewritten with real option')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_override_existing_file

!---------------------------------------------------------------
!+
!  with no file, options on the command line are used with the
!  defaults instead of prompting, and the file is written
!+
!---------------------------------------------------------------
subroutine test_override_no_file(ntests,npass)
 use infile_utils, only:get_options,set_overrides,clear_overrides
 integer, intent(inout) :: ntests,npass
 integer :: ierr,nfailed(4)

 if (id==master) write(*,"(/,a)") '--> testing command line options with no file'
 nfailed = 0
 call set_defaults()
 call delete_file(filename)
 call set_overrides([character(len=20) :: '--np=9'])
 call get_options(filename,.true.,ierr,read_pars,write_pars,set_pars)
 call clear_overrides()
 call checkval(ierr,1,0,nfailed(1),'ierr=1 (file written, rerun)')
 call checkval(prompted,.false.,nfailed(2),'no prompting')

 call set_defaults()
 call get_options(filename,.true.,ierr,read_pars,write_pars,set_pars)
 call checkval(ierr,0,0,nfailed(3),'file read on rerun')
 call checkval(np,9,0,nfailed(4),'file written with option value')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_override_no_file

!-------------------------------------------------------------
!+
!  with no file and no command line options, the options are
!  prompted for and the file is written
!+
!-------------------------------------------------------------
subroutine test_prompt_no_file(ntests,npass)
 use infile_utils, only:get_options
 integer, intent(inout) :: ntests,npass
 integer :: ierr,nfailed(4)

 if (id==master) write(*,"(/,a)") '--> testing prompting with no file'
 nfailed = 0
 call set_defaults()
 call delete_file(filename)
 call get_options(filename,.true.,ierr,read_pars,write_pars,set_pars)
 call checkval(ierr,1,0,nfailed(1),'ierr=1 (file written, rerun)')
 call checkval(prompted,.true.,nfailed(2),'options prompted for')

 call set_defaults()
 call get_options(filename,.true.,ierr,read_pars,write_pars,set_pars)
 call checkval(ierr,0,0,nfailed(3),'file read on rerun')
 call checkval(np,42,0,nfailed(4),'file written with prompted value')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_prompt_no_file

!---------------------------------------------------
!+
!  missing options are added to an incomplete file
!+
!---------------------------------------------------
subroutine test_incomplete_file(ntests,npass)
 use infile_utils, only:get_options
 integer, intent(inout) :: ntests,npass
 integer :: ierr,iunit,nfailed(3)

 if (id==master) write(*,"(/,a)") '--> testing incomplete file'
 nfailed = 0
 call set_defaults()
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") 'np = 5'
 close(iunit)
 call get_options(filename,.true.,ierr,read_pars,write_pars)
 call checkval(ierr /= 0,.true.,nfailed(1),'error for missing options')

 call set_defaults()
 call get_options(filename,.true.,ierr,read_pars,write_pars)
 call checkval(ierr,0,0,nfailed(2),'missing options added to file')
 call checkval(np,5,0,nfailed(3),'existing option kept')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_incomplete_file

!------------------------------------------------------------
!+
!  options on the command line that are not in the file are
!  reported and do not change the values read
!+
!------------------------------------------------------------
subroutine test_unknown_override(ntests,npass)
 use infile_utils, only:get_options,set_overrides,clear_overrides,count_unused_overrides
 integer, intent(inout) :: ntests,npass
 integer :: ierr,nfailed(3)

 if (id==master) write(*,"(/,a)") '--> testing unknown command line option'
 nfailed = 0
 call set_defaults()
 call write_pars(filename)
 call set_overrides([character(len=20) :: '--npp=3'])
 call get_options(filename,.true.,ierr,read_pars,write_pars)
 call checkval(count_unused_overrides(),1,0,nfailed(1),'unknown option detected')
 call clear_overrides()
 call checkval(ierr,0,0,nfailed(2),'ierr')
 call checkval(np,100,0,nfailed(3),'value from file unchanged')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_unknown_override

!---------------------------------------------------------------
!+
!  read_inopt for strings, with the variable itself as default
!+
!---------------------------------------------------------------
subroutine test_string_default(ntests,npass)
 use infile_utils, only:open_db_from_file,close_db,inopts,read_inopt
 integer, intent(inout) :: ntests,npass
 type(inopts), allocatable :: db(:)
 integer :: ierr,nfailed(1)

 if (id==master) write(*,"(/,a)") '--> testing string option with default'
 nfailed = 0
 call set_defaults()
 call write_pars(filename)
 label = 'xyz'
 call open_db_from_file(db,filename,23,ierr)
 call read_inopt(label,'missing_tag',db,err=ierr,default=label)
 call close_db(db)
 call checkval(label,'xyz',nfailed(1),'default kept for missing string option')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_string_default

!--------------------------------------------
!+
!  default values of the test options
!+
!--------------------------------------------
subroutine set_defaults()

 np       = 100
 omega    = 0.5
 lflag    = .false.
 label    = 'abc'
 prompted = .false.

end subroutine set_defaults

!--------------------------------------------
!+
!  "interactive" setting of the test options
!+
!--------------------------------------------
subroutine set_pars()

 prompted = .true.
 np = 42

end subroutine set_pars

!--------------------------------------------
!+
!  write the test parameter file
!+
!--------------------------------------------
subroutine write_pars(fname)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: fname
 integer :: iunit

 open(newunit=iunit,file=fname,status='replace',form='formatted')
 write(iunit,"(a)") '# test parameter file'
 call write_inopt(np,'np','number of particles',iunit)
 call write_inopt(omega,'omega','a real option',iunit)
 call write_inopt(lflag,'lflag','a logical option',iunit)
 call write_inopt(label,'label','a string option',iunit)
 close(iunit)

end subroutine write_pars

!--------------------------------------------
!+
!  read the test parameter file
!+
!--------------------------------------------
subroutine read_pars(fname,ierr)
 use infile_utils, only:open_db_from_file,close_db,inopts,read_inopt
 character(len=*), intent(in)  :: fname
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,fname,23,ierr)
 if (ierr /= 0) return
 call read_inopt(np,'np',db,errcount=nerr)
 call read_inopt(omega,'omega',db,errcount=nerr)
 call read_inopt(lflag,'lflag',db,errcount=nerr)
 call read_inopt(label,'label',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_pars

!--------------------------------------------
!+
!  delete a file if it exists
!+
!--------------------------------------------
subroutine delete_file(fname)
 character(len=*), intent(in) :: fname
 integer :: iunit,ierr

 open(newunit=iunit,file=fname,status='old',iostat=ierr)
 if (ierr == 0) close(iunit,status='delete')

end subroutine delete_file

end module testinfile
