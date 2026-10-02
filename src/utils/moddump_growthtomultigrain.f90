!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! transforms dustgrowth dump into multigrain dump for mcfost usage
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - bins_per_dex : *number of bins per dex*
!   - force_smax   : *set the maximum grain size manually*
!   - smax_user    : *maximum grain size [cm] (used if force_smax=T)*
!
! :Dependencies: dim, growth, infile_utils, io_control, moddump_utils,
!   part, prompting
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer :: bins_per_dex = 5        ! number of bins per dex
 real    :: smax_user    = 2.       ! maximum grain size in cm (used if force_smax)
 logical :: force_smax   = .false.  ! set the maximum grain size manually

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:use_dust,use_dustgrowth
 use part,           only:delete_dead_or_accreted_particles
 use io_control,     only:nmax
 use growth,         only:bin_to_multi
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 logical                :: file_exists
 character(len=20)      :: infile  = "bin_param.txt"

 if ((.not. use_dust) .or. (.not. use_dustgrowth)) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes AND DUSTGROWTH=yes'
    stop
 endif

 nmax = 0 !- deriv called once after moddump

 !- check if param file exists, created by python script growthtomcfost.py
 inquire(file=infile, exist=file_exists)

 if (file_exists .and. prompt_for_params) then
    !- file created by phantom/scripts/growthtomcfost.py module
    open(unit=420,file=infile)
    read(420,*) force_smax, smax_user, bins_per_dex
    close(unit=420)
 else
    !- otherwise read the moddump parameter file (or write a template and stop)
    if (prompt_for_params) call read_interactive_moddumpfile()
 endif

 !- delete dead or accreted particles before doing anything
 call delete_dead_or_accreted_particles(npart,npartoftype)

 !- bin dust particles into desired bins
 call bin_to_multi(bins_per_dex,force_smax,smax_user,verbose=.true.)

end subroutine modify_dump

!
!---Interactively set the moddump parameters--------------------------------
!
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Set smax manually?',force_smax)
 if (force_smax) call prompt('Enter smax in cm',smax_user,0.05)
 call prompt('Enter number of bins per dex',bins_per_dex,1)

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
 call write_inopt(force_smax,'force_smax','set the maximum grain size manually',iunit)
 call write_inopt(smax_user,'smax_user','maximum grain size [cm] (used if force_smax=T)',iunit)
 call write_inopt(bins_per_dex,'bins_per_dex','number of bins per dex',iunit)
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
 call read_inopt(force_smax,'force_smax',db,errcount=nerr)
 call read_inopt(smax_user,'smax_user',db,min=0.,errcount=nerr)
 call read_inopt(bins_per_dex,'bins_per_dex',db,min=1,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

end module moddump
