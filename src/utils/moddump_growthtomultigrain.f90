!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! transforms dustgrowth dump into multigrain dump for mcfost usage
!
! :References: None
!
! :Owner: Arnaud Vericel
!
! :Runtime parameters: None
!
! :Dependencies: dim, growth, part, prompting, timestep
!
 implicit none
 !- initialise variables
 integer, private :: bins_per_dex = 5
 real, private :: smax_user    = 2.
 logical, private :: force_smax   = .false.
 logical, private :: use_moments  = .true.

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:use_dust,use_dustgrowth
 use part,           only:delete_dead_or_accreted_particles
 use prompting,      only:prompt
 use timestep,       only:nmax
 use growth,         only:bin_to_multi,init_growth
 use io,             only:fatal
 use deriv,          only:get_derivs_global
 use dust,           only:init_drag
 use initial,        only:initialise
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: iu,ierr
 logical                :: file_exists
 character(len=20)      :: infile  = "bin_param.txt"

 if ((.not. use_dust) .or. (.not. use_dustgrowth)) then
    call fatal('moddump_growth2multi','Need to compile with DUST=yes and DUSTGROWTH=yes')
 endif

 nmax         = 0 !- deriv called once after moddump

 !- check if param file exists, created by python script growthtomcfost.py
 inquire(file=infile, exist=file_exists)

 !- if not, switch to interactive method
 if (.not.file_exists) then
    call prompt('Set smax manually?',force_smax)
    if (force_smax) call prompt('Enter smax in cm',smax_user,0.05)
    call prompt('Enter number of bins per dex',bins_per_dex,1)
    ! write the input file so don't have to ask questions again
    call write_options_moddump(infile)
 else
    call read_options_moddump(infile,ierr)
    if (ierr /= 0) then
       call write_options_moddump(infile)
       stop
    endif
 endif

 !- delete dead or accreted particles before doing anything
 call delete_dead_or_accreted_particles(npart,npartoftype)

 !- bin dust particles into desired bins
 call bin_to_multi(bins_per_dex,force_smax,smax_user,use_moments,verbose=.true.)

 !--now interpolate the dust density onto the gas particles
 call initialise()
 call init_drag(ierr)
 if (use_dustgrowth) call init_growth(ierr)
 call get_derivs_global()

end subroutine modify_dump

!------------------------------------------------------
!+
!  write the parameter file for the moddump procedure
!+
!------------------------------------------------------
subroutine write_options_moddump(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit

 open(newunit=iunit,file=filename,status='replace',form='formatted')
 call write_inopt(force_smax,'force_smax','set max grain size manually?',iunit)
 call write_inopt(smax_user,'smax_user','user-defined max grain size',iunit)
 call write_inopt(bins_per_dex,'bins_per_dex','number of grain size bins per dex',iunit)
 call write_inopt(use_moments,'use_moments','reconstruct grain size distribution from moments',iunit)
 close(iunit)

end subroutine write_options_moddump

!------------------------------------------------------
!+
!  read the parameter file for the moddump procedure
!+
!------------------------------------------------------
subroutine read_options_moddump(filename,ierr)
 use infile_utils, only:inopts,open_db_from_file,close_db,read_inopt
 use io,           only:id,master
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 !- file created by phantom/scripts/growthtomcfost.py module
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(force_smax,'force_smax',db,errcount=nerr)
 call read_inopt(smax_user,'smax_user',db,errcount=nerr)
 call read_inopt(bins_per_dex,'bins_per_dex',db,errcount=nerr)
 call read_inopt(use_moments,'use_moments',db,errcount=nerr)
 call close_db(db)

 if (nerr > 0 .and. id==master) then
    print "(1x,a,i2,a)",'read_options_moddump: ',nerr,' error(s) during read of parameter file.  Re-writing.'
    ierr = nerr
 endif

end subroutine read_options_moddump

end module moddump
