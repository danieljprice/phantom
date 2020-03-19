program get_kdot
 use evutils
 use fileutils, only:files_are_sequential
 implicit none
 integer, parameter :: max_files = 128
 integer            :: ncols,ientrop_col,iend,nfiles
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
 ientrop_col = 0   ! to prevent compiler warnings
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
    ! open the output (.kdot) file
    !
    if (.not.combine_files .or. (combine_files .and. i==istart)) then
       ncols  = ncols1
       !
       ! get labels and number of columns from the first line of the file
       !
       ientrop_col = find_column(labels,'totentrop') ! label used in standard .ev files
       if (ientrop_col <= 0) then
          print*,"(a)",'ERROR: could not locate entropy column from header information'
          if (combine_files) then
             exit over_files
          else
             cycle over_files
          endif
       endif
       !
       ! set output file name
       !
       iev = index(filename,'.ev') - 1
       if (iev <= 0) iev = len_trim(filename)
       if (combine_files) then
          outfile = trim(filename(1:iev-2))//'.kdot'
       else
          outfile = trim(filename(1:iev))//'.kdot'
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
       call compute_and_write_kdot(outfile,dat,dtmin,nsteps,ncols,ientrop_col)
    endif
    !
    ! close .ev file and deallocate memory
    !
    deallocate(dat)

    nsteps = nsteps + nsteps_keep
    nsteps_prev = nsteps

 enddo over_files

 if (combine_files .and. ientrop_col > 0 .and. allocated(dat_combined)) then
    print "(a)",' --> '//trim(outfile)
    call compute_and_write_kdot(outfile,dat_combined,dtmin,nsteps,ncols,ientrop_col)
 endif

 if (allocated(dat))          deallocate(dat)
 if (allocated(dat_combined)) deallocate(dat_combined)

contains

subroutine compute_and_write_kdot(outfile,dat,dtmin,nsteps,ncol,ientrop_col)
 character(len=*), intent(in) :: outfile
 real,             intent(in) :: dat(:,:)
 real,             intent(in) :: dtmin
 integer,          intent(in) :: nsteps,ncol,ientrop_col
 real :: kprev,k,kdot
 real :: tprev,dt
 integer :: i
 integer, parameter :: iout = 14

 !
 ! open the output file for output
 !
 open(unit=iout,file=trim(outfile),status='replace',form='formatted')
 write(iout,"('#',3(a17,2x))") 't','kdot','ktot'
 !
 ! Now calculate the accretion rates
 !
 tprev = 0.
 kprev = 0.

 do i=1,nsteps
    k = dat(ientrop_col,i)
    dt = dat(1,i) - tprev
    if (dt > tiny(dt)) then
       kdot = max(0., (k - kprev)/dt)
    else
       kdot = 0.
    endif
    if (dt > dtmin) then
       write(iout,"(3(es18.10,1x))",iostat=ierr) dat(1,i),kdot,k
       tprev = dat(1,i)
       kprev = k
    endif
 enddo

 close(unit=iout)

end subroutine compute_and_write_kdot

end program get_kdot
