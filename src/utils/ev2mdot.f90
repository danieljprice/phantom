!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: get_mdot
!
!  DESCRIPTION: Utility to extract mdot as a function of time
!    from the .ev file
!    Daniel Price, March 2011
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: ev2mdot [dt_min] file1.ev file2.ev ...
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
program get_mdot
 use evutils
 implicit none
 integer, parameter :: lu1 = 11, lu2 = 13, iout = 12
 integer, parameter :: maxcol = 128
 integer :: ncols,imacc_col,imacc1_col,imacc2_col
 real    :: dat(maxcol),dat2(maxcol-1)
 character(len=20) :: labels(maxcol)
 real, dimension(:,:), allocatable :: datfull,datfull2,temp_array
 character(len=40) :: filename,outfile
 integer :: i,nargs,ierr,nlines,iev,isteps,istart,append_length,stepstart,ii,jj,stepend
 real    :: maccprev,maccprev1,maccprev2,macc,macc1,macc2,mdot,mdot1,mdot2
 real    :: tprev,dt,dtmin,thold

 !
 ! get number of arguments and name of current program
 !
 nargs = command_argument_count()
 call get_command_argument(0,filename)

 if (nargs <= 0) then
    print "(a)",' Usage: '//trim(filename)//' [dt_min] file1.ev file2.ev ...'
    stop
 endif
 !
 ! set the time interval based on the (optional) first command line argument
 ! if not read as a sensible real number, assume dtmin = 0.
 !
 call get_command_argument(1,filename)
 read(filename,*,iostat=ierr) dtmin
 if (ierr /= 0) then
    dtmin = 0.
    istart = 1
 else
    print*,'Using minimum time interval = ',dtmin
    istart = 2
 endif

 !
 ! get name of the first file
 !
 call get_command_argument(istart,filename)
 !
 ! get labels and number of columns from the first line of the file
 !
 call get_column_labels_from_ev(filename,labels,ncols,ierr)
 if (ncols < 0) then
    print "(a,i2)",' ERROR: could not determine number of columns in '//trim(filename)//' ...skipping'
    stop
 endif
 imacc_col = find_column(labels,'accretedmas') ! global .ev file
 if (imacc_col <= 0) imacc_col = find_column(labels,'macc') ! sink particle .ev file
 if (imacc_col <= 0) stop 'could not locate accreted mass column from header information'

 !
 ! find column numbers for binary mass accretion rates
 ! (only if binary information exists in .ev file)
 !
 imacc1_col = find_column(labels,'macc1')
 imacc2_col = find_column(labels,'macc2')

 !
 ! open the file for reading and set the name of the output file
 !
 open(unit=lu1,file=trim(filename),status='old',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)//' ...skipping'
 else
    ! open output file
    iev = index(filename,'.ev') - 1
    if (iev <= 0) iev = len_trim(filename)
    ! outfile = filename(1:iev)//'-mdot.ev'
    outfile = trim(filename(1:iev-2))//'.mdot'
    open(unit=iout,file=trim(outfile),status='replace',form='formatted')
    if (imacc1_col > 0 .and. imacc2_col > 0) then
       write(iout,"('#',7(a17,2x))") 't','mdot','macc','mdot1','macc1','mdot2','macc2'
    else
       write(iout,"('#',3(a17,2x))") 't','mdot','macc'
    endif
 endif
 !
 ! allocate memory
 !
 allocate(datfull2(ncols,1))
 datfull2(:,:) = 0.

 stepstart = 1
 isteps = 1
 stepend = stepstart
 imacc_col = 0

 over_args: do i=istart,nargs
    print "(a)",trim(filename)//' --> '//trim(outfile)
    call get_command_argument(i,filename)
    
    !
    ! now open the file properly
    !
    open(unit=lu1,file=trim(filename),status='old',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR opening '//trim(filename)//' ...skipping'
    else
       ! skip header lines
       ierr = 1
       nlines = 0
       do while(ierr /= 0 .and. nlines <= 200)
          nlines = nlines + 1
          read(lu1,*,iostat=ierr) dat(1:ncols)
       enddo

       if (nlines==200) then
          print*,' error: only one line or fewer in file'
          close(unit=lu1)
          cycle over_args
       endif

       !
       ! count the number of useable lines (timesteps) in the file
       !
       isteps = 0
       do while(ierr == 0)
          isteps = isteps + 1
          read(lu1,*,iostat=ierr) !dat(1:ncols)
       enddo
       !
       ! allocate enough memory to read it all
       !
       allocate(datfull(1:ncols,1:isteps))
       datfull(:,:) = 0.

       !
       ! rewind file
       !
       rewind(unit=lu1)

       !
       ! skip header lines again
       !
       do jj=1,nlines
          read(lu1,*,iostat=ierr)
       enddo
       !
       ! read the data
       !
       do jj = 1,isteps
          read(lu1,*,iostat=ierr) dat(1:ncols)
          datfull(1:ncols,jj) = dat(1:ncols)
       enddo

       ! Work out where the overlap between the files is
       call get_command_argument(i+1,filename)
       open(unit=lu2,file=trim(filename),status='old',form='formatted',iostat=ierr)
       if (ierr /= 0) then
          append_length = isteps-1
       else
          read(lu2,*,iostat=ierr)
          read(lu2,*,iostat=ierr) thold,dat2(1:ncols-1)
          append_length = minloc(abs(datfull(1,:)-thold),1)-1
       endif
       close(unit=lu2)

       ii = 1
       stepend = stepstart + append_length
       ! print*,'Loop start at',stepstart,'/',datfull(1,1),'finish at',stepend,'/',thold

       ! Resize datfull2 here
       allocate(temp_array(ncols,stepstart))
       temp_array(:,1:stepstart) = datfull2(:,:)
       deallocate(datfull2)
       ! Reassign all values to datfull2
       allocate(datfull2(ncols,1:stepend))
       datfull2(:,1:stepstart) = temp_array(:,:)
       datfull2(:,stepstart:stepend) = datfull(:,1:(append_length+1))
       deallocate(temp_array)
       close(unit=lu1)
       deallocate(datfull)
    endif
    stepstart=stepend
 enddo over_args

 !
 ! Now calculate the accretion rates
 !
 isteps = 0
 tprev = 0
 maccprev = 0.
 maccprev1 = 0.
 maccprev2 = 0.
 macc1 = 0.
 macc2 = 0.
 dat(:) = 0.
 ierr = 0

 do while(isteps  <=  (size(datfull2,2)-1))
    isteps = isteps + 1
    macc = datfull2(imacc_col,isteps) !dat(imacc_col)
    if (imacc1_col > 0) macc1 = datfull2(imacc1_col,isteps)
    if (imacc2_col > 0) macc2 = datfull2(imacc2_col,isteps)
    dt = datfull2(1,isteps) - tprev
    if (dt > tiny(dt)) then
       mdot = max(0., (macc - maccprev)/dt)
       mdot1 = max(0., (macc1 - maccprev1)/dt)
       mdot2 = max(0., (macc2 - maccprev2)/dt)
    else
       mdot = 0.
       mdot1 = 0.
       mdot2 = 0.
    endif
    if (dt > dtmin .and. ierr==0) then
       !print *,' time ',datfull2(1,isteps),' mass accreted = ',macc, ' mdot = ',mdot
       if (imacc1_col > 0 .and. imacc2_col > 0) then
          write(iout,"(3(es18.10,1x))",iostat=ierr) datfull2(1,isteps),mdot,macc,mdot1,macc1,mdot2,macc2
       else
          write(iout,"(3(es18.10,1x))",iostat=ierr) datfull2(1,isteps),mdot,macc
       endif
       tprev = datfull2(1,isteps)
       maccprev = macc
       maccprev1 = macc1
       maccprev2 = macc2
    endif
 enddo

 close(unit=iout)

end program get_mdot
