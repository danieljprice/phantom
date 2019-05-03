!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: get_mdot
!
!  DESCRIPTION: Utility to extract mdot as a function of time
!    from the .ev file. Splices the sequence of .ev files back
!    together to reconstruct a single file dataset
!    Daniel Price, March 2011, rewritten in Aug 2017
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: ev2mdot [dt_min] file01.ev file02.ev ...
!
!  DEPENDENCIES: evutils, fileutils
!+
!--------------------------------------------------------------------------
program get_mdot
 use evutils
 use fileutils, only:files_are_sequential
 implicit none
 integer, parameter :: max_files = 128
 integer            :: ncols,imacc_col,imacc1_col,imacc2_col,iend,nfiles
 character(len=20)  :: labels(max_columns)
 real, dimension(:,:), allocatable :: dat,dat_combined,temp_array
 character(len=40)  :: string,filenames(max_files),outfile,filename
 integer :: i,istart,nargs,ierr,iev,ncols1,nsteps,nsteps_prev,nsteps_keep
 logical :: combine_files
 real    :: dtmin,dt,tfile

 !
 ! get number of arguments and name of current program
 !
 nargs = command_argument_count()
 call get_command_argument(0,string)

 if (nargs <= 0) then
    print "(a)",' Usage: '//trim(string)//' [dt_min] file01.ev file02.ev ...'
    stop
 endif
 !
 ! set the time interval based on the (optional) first command line argument
 ! if not read as a sensible real number, assume dtmin = 0.
 !
 call get_command_argument(1,string)
 read(string,*,iostat=ierr) dtmin
 if (ierr /= 0) then
    dtmin = 0.
    istart = 1
 else
    print*,'Using minimum time interval = ',dtmin
    istart = 2
 endif
 !
 ! extract filenames from command line
 !
 iend = min(nargs,max_files)
 do i=istart,iend
    call get_command_argument(i,filenames(i))
 enddo
 nfiles = iend - istart + 1
 !
 ! check if filenames form a sequence file01.ev, file02.ev, file03.ev
 !
 combine_files = files_are_sequential(filenames(istart:iend))

 nsteps_prev = 0
 nsteps_keep = 0
 imacc_col = 0   ! to prevent compiler warnings
 imacc1_col = 0
 imacc2_col = 0
 ncols = 0
 over_files: do i=istart,iend
    filename = filenames(i)
    !
    ! read the .ev file
    !
    call read_evfile(filename,dat,labels,ncols1,nsteps,ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR opening '//trim(filename)//' ...skipping'
       cycle over_files
    endif
    !
    ! open the output (.mdot) file
    !
    if (.not.combine_files .or. (combine_files .and. i==istart)) then
       ncols  = ncols1
       !
       ! get labels and number of columns from the first line of the file
       !
       imacc_col = find_column(labels,'accretedmas') ! label used in standard .ev files
       if (imacc_col <= 0) imacc_col = find_column(labels,'macc') ! label used in sink particle .ev files
       if (imacc_col <= 0) imacc_col = find_column(labels,'accretedm') ! label used in old .ev files
       if (imacc_col <= 0) then
          print*,"(a)",'ERROR: could not locate accreted mass column from header information'
          if (combine_files) then
             exit over_files
          else
             cycle over_files
          endif
       endif
       !
       ! find column numbers for binary mass accretion rates
       ! (only if binary information exists in .ev file)
       !
       imacc1_col = find_column(labels,'macc1')
       imacc2_col = find_column(labels,'macc2')
       !
       ! set output file name
       !
       iev = index(filename,'.ev') - 1
       if (iev <= 0) iev = len_trim(filename)
       if (combine_files) then
          outfile = trim(filename(1:iev-2))//'.mdot'
       else
          outfile = trim(filename(1:iev))//'.mdot'
       endif
       if (combine_files) then
          allocate(dat_combined(ncols,nsteps))
          dat_combined = dat
       endif
    endif

    nsteps_keep = 0
    if (combine_files) then
       if (i > 1) then
          !
          ! exclude any times from previous files that
          ! overlap with times in the current file
          !
          tfile = dat(1,1)   ! earliest time in current file
          dt = 1.
          do while(dt > 0. .and. nsteps_keep < nsteps_prev)
             nsteps_keep = nsteps_keep + 1
             dt = tfile - dat_combined(1,nsteps_keep)
          enddo
          if (dt <= 0.) nsteps_keep = nsteps_keep - 1
          print "(a,i6,a)",trim(filename) !//' (',nsteps_prev-nsteps_keep,' steps overlap with previous file)'

          ! resize dat_combined array
          allocate(temp_array(ncols,nsteps_keep))
          temp_array(:,1:nsteps_keep) = dat_combined(:,1:nsteps_keep)
          if (allocated(dat_combined)) deallocate(dat_combined)

          ! Reassign all values to dat_combined
          allocate(dat_combined(ncols,nsteps_keep + nsteps))
          dat_combined(:,1:nsteps_keep) = temp_array(:,:)
          dat_combined(:,nsteps_keep+1:nsteps_keep+nsteps) = dat(1:ncols,1:nsteps)
          if (allocated(temp_array)) deallocate(temp_array)
       else
          print "(a)",trim(filename)
       endif
    else
       print "(a)",trim(filename)//' --> '//trim(outfile)
       call compute_and_write_mdot(outfile,dat,dtmin,nsteps,ncols,imacc_col,imacc1_col,imacc2_col)
    endif
    !
    ! close .ev file and deallocate memory
    !
    deallocate(dat)

    nsteps = nsteps + nsteps_keep
    nsteps_prev = nsteps

 enddo over_files

 if (combine_files .and. imacc_col > 0 .and. allocated(dat_combined)) then
    print "(a)",' --> '//trim(outfile)
    call compute_and_write_mdot(outfile,dat_combined,dtmin,nsteps,ncols,imacc_col,imacc1_col,imacc2_col)
 endif

 if (allocated(dat))          deallocate(dat)
 if (allocated(dat_combined)) deallocate(dat_combined)

contains

subroutine compute_and_write_mdot(outfile,dat,dtmin,nsteps,ncol,imacc_col,imacc1_col,imacc2_col)
 character(len=*), intent(in) :: outfile
 real,             intent(in) :: dat(:,:)
 real,             intent(in) :: dtmin
 integer,          intent(in) :: nsteps,ncol,imacc_col,imacc1_col,imacc2_col
 real :: maccprev,maccprev1,maccprev2,macc,macc1,macc2,mdot,mdot1,mdot2
 real :: tprev,dt
 integer :: i
 integer, parameter :: iout = 14

 !
 ! open the output file for output
 !
 open(unit=iout,file=trim(outfile),status='replace',form='formatted')
 if (imacc1_col > 0 .and. imacc2_col > 0) then
    write(iout,"('#',7(a17,2x))") 't','mdot','macc','mdot1','macc1','mdot2','macc2'
 else
    write(iout,"('#',3(a17,2x))") 't','mdot','macc'
 endif
 !
 ! Now calculate the accretion rates
 !
 tprev = 0.
 maccprev = 0.
 maccprev1 = 0.
 maccprev2 = 0.
 macc1 = 0.
 macc2 = 0.

 do i=1,nsteps
    macc = dat(imacc_col,i)
    if (imacc1_col > 0 .and. imacc1_col <= ncol) macc1 = dat(imacc1_col,i)
    if (imacc2_col > 0 .and. imacc2_col <= ncol) macc2 = dat(imacc2_col,i)
    dt = dat(1,i) - tprev
    if (dt > tiny(dt)) then
       mdot = max(0., (macc - maccprev)/dt)
       mdot1 = max(0., (macc1 - maccprev1)/dt)
       mdot2 = max(0., (macc2 - maccprev2)/dt)
    else
       mdot = 0.
       mdot1 = 0.
       mdot2 = 0.
    endif
    if (dt > dtmin) then
       if (imacc1_col > 0 .and. imacc2_col > 0) then
          write(iout,"(3(es18.10,1x))",iostat=ierr) dat(1,i),mdot,macc,mdot1,macc1,mdot2,macc2
       else
          write(iout,"(3(es18.10,1x))",iostat=ierr) dat(1,i),mdot,macc
       endif
       tprev = dat(1,i)
       maccprev = macc
       maccprev1 = macc1
       maccprev2 = macc2
    endif
 enddo

 close(unit=iout)

end subroutine compute_and_write_mdot

end program get_mdot
