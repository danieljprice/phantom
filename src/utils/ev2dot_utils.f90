!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module ev2dotutils
!
! Utility routines to extract mdot as a function of time
! from the .ev file. Splices the sequence of .ev files back
! together to reconstruct a single file dataset
!
! Daniel Price, March 2011, rewritten in Aug 2017, modularised 2025
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: evutils, fileutils
!
 use evutils
 use fileutils, only:files_are_sequential
 implicit none

 private

 public :: ev2dot

contains

!-----------------------------------------------------------------------
!+
!  Calculate mdot or kdot from a sequence of .ev files and write to .mdot file
!
!  output_type = 'kdot': compute kdot
!  output_type = 'mdot': compute mdot (default if output_type not set)
!+
!-----------------------------------------------------------------------
subroutine ev2dot(nfiles,filenames,dtmin,output_type)
 integer,          intent(in) :: nfiles
 character(len=*), intent(in) :: filenames(nfiles)
 real,             intent(in) :: dtmin
 character(len=4), intent(in), optional :: output_type
 integer            :: ncols,imacc_col,imacc1_col,imacc2_col,ientrop_col
 character(len=20)  :: labels(max_columns)
 real, dimension(:,:), allocatable :: dat,dat_combined,temp_array
 character(len=120) :: outfile,filename
 character(len=5)   :: ext
 integer :: i,ierr,iev,ncols1,nsteps,nsteps_prev,nsteps_keep
 logical :: combine_files
 real    :: dt,tfile

 if (present(output_type)) then
    ext = '.'//output_type
 else
    ext = '.mdot'
 endif
 !
 ! check if filenames form a sequence file01.ev, file02.ev, file03.ev
 !
 combine_files = files_are_sequential(filenames)

 nsteps_prev = 0
 nsteps_keep = 0
 imacc_col = 0   ! to prevent compiler warnings
 imacc1_col = 0
 imacc2_col = 0
 ientrop_col = 0
 ncols = 0
 ierr = 0
 over_files: do i=1,nfiles
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
    if (.not.combine_files .or. (combine_files .and. i==1)) then
       ncols  = ncols1
       !
       ! get labels and number of columns from the first line of the file
       !
       select case(ext)
       case('.kdot')
          ientrop_col = find_column(labels,'totentrop') ! label used in standard .ev files
          if (ientrop_col <= 0) then
             print*,"(a)",'ERROR: could not locate entropy column from header information'
             if (combine_files) then
                exit over_files
             else
                cycle over_files
             endif
          endif
       case default
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
       end select
       !
       ! set output file name
       !
       iev = index(filename,'.ev') - 1
       if (iev <= 0) iev = len_trim(filename)
       if (combine_files) then
          outfile = trim(filename(1:iev-2))//ext
       else
          outfile = trim(filename(1:iev))//ext
       endif
       if (combine_files) then
          allocate(dat_combined(ncols,nsteps))
          dat_combined = dat
       endif
    else
       if (.not.allocated(dat_combined)) allocate(dat_combined(0,0)) ! to avoid compiler warning
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
       select case(ext)
       case('.kdot')
          call compute_and_write_kdot(outfile,dat,dtmin,nsteps,ncols,ientrop_col)
       case default
          call compute_and_write_mdot(outfile,dat,dtmin,nsteps,ncols,imacc_col,imacc1_col,imacc2_col)
       end select
    endif
    !
    ! close .ev file and deallocate memory
    !
    deallocate(dat)

    nsteps = nsteps + nsteps_keep
    nsteps_prev = nsteps

 enddo over_files

 if (.not.allocated(dat_combined)) allocate(dat_combined(0,0)) ! to avoid compiler warning

 if (combine_files .and. allocated(dat_combined)) then
    print "(a)",' --> '//trim(outfile)
    select case(ext)
    case('.kdot')
       call compute_and_write_kdot(outfile,dat_combined,dtmin,nsteps,ncols,ientrop_col)
    case default
       if (imacc_col > 0) call compute_and_write_mdot(outfile,dat_combined,dtmin,nsteps,ncols,&
                                                      imacc_col,imacc1_col,imacc2_col)
    end select
 endif

 if (allocated(dat))          deallocate(dat)
 if (allocated(dat_combined)) deallocate(dat_combined)

end subroutine ev2dot

!-----------------------------------------------------------------------
!+
!  Internal routine to compute and write mdot from the data
!+
!-----------------------------------------------------------------------
subroutine compute_and_write_mdot(outfile,dat,dtmin,nsteps,ncol,imacc_col,imacc1_col,imacc2_col)
 character(len=*), intent(in) :: outfile
 real,             intent(in) :: dat(:,:)
 real,             intent(in) :: dtmin
 integer,          intent(in) :: nsteps,ncol,imacc_col,imacc1_col,imacc2_col
 real :: maccprev,maccprev1,maccprev2,macc,macc1,macc2,mdot,mdot1,mdot2
 real :: tprev,dt
 integer :: i,ierr
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
    if (dt >= dtmin) then
       if (imacc1_col > 0 .and. imacc2_col > 0 .and. i > 1) then
          write(iout,"(3(es18.10,1x))",iostat=ierr) dat(1,i),mdot,macc,mdot1,macc1,mdot2,macc2
       elseif (i > 1) then
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

!-----------------------------------------------------------------------
!+
!  Internal routine to compute and write the entropy derivative from the data
!+
!-----------------------------------------------------------------------
subroutine compute_and_write_kdot(outfile,dat,dtmin,nsteps,ncol,ientrop_col)
 character(len=*), intent(in) :: outfile
 real,             intent(in) :: dat(:,:)
 real,             intent(in) :: dtmin
 integer,          intent(in) :: nsteps,ncol,ientrop_col
 real :: kprev,k,kdot
 real :: tprev,dt
 integer :: i,ierr
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
       kdot = (k - kprev)/dt
    else
       kdot = 0.
    endif
    !print*,'k is ',k,' was ',kprev,' dt = ',dt,' gradient = ',(k-kprev)/dt
    if (dt >= dtmin) then
       if (i > 1) write(iout,"(3(es18.10,1x))",iostat=ierr) dat(1,i),kdot,k
       tprev = dat(1,i)
       kprev = k
    endif
 enddo

 close(unit=iout)

end subroutine compute_and_write_kdot

end module ev2dotutils
