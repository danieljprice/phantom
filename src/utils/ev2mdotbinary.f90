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
!  OWNER: Chris Nixon
!
!  $Id$
!
!  USAGE: ev2mdotbinary [dt_min] file1.ev file2.ev ...
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
program get_mdot
 implicit none
 integer, parameter :: lu1 = 11, lu2 = 13, iout = 12
! integer, parameter :: nsmooth = 20
 integer, parameter :: imacc_col = 18
 integer, parameter :: imacc1_col = 20
 integer, parameter :: imacc2_col = 21
 integer, parameter :: ncols = imacc2_col
 real :: dat(ncols),dat2(ncols-1)
 real, dimension(:,:), allocatable :: datfull,datfull2,temp_array
 character(len=40) :: filename,outfile
 integer :: i,nargs,ierr,nlines,iev,isteps,istart,append_length,stepstart,ii,jj,stepend
 real :: maccprev,macc,tprev,dt,mdot,dtmin,thold
 real :: maccprev1,maccprev2,mdot1,mdot2,macc1,macc2

 allocate(datfull2(ncols,1))
 datfull2(:,:) = 7.

 nargs = command_argument_count()
 call get_command_argument(0,filename)

 if (nargs <= 0) then
    print "(a)",' Usage: '//trim(filename)//' [dt_min] file1.ev file2.ev ...'
    stop
 endif
 call get_command_argument(1,filename)
 read(filename,*,iostat=ierr) dtmin
 if (ierr /= 0) then
    dtmin = 0.
    istart = 1
 else
    print*,'Using minimum time interval = ',dtmin
    istart = 2
 endif

 call get_command_argument(istart,filename)
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
    write(iout,"('#',7(a17,2x))") 't','mdot','macc','mdot1','macc1','mdot2','macc2'
 endif

 stepstart = 1
 isteps = 1
 stepend = stepstart

 over_args: do i=istart,nargs
    print "(a)",trim(filename)//' --> '//trim(outfile)
    call get_command_argument(i,filename)
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

       isteps = 0
       do while(ierr == 0)
          isteps = isteps + 1
          read(lu1,*,iostat=ierr) !dat(1:ncols)
       enddo
       allocate(datfull(1:ncols,1:isteps))
       datfull(:,:) = 0.

       close(unit=lu1)
       open(unit=lu1,file=trim(filename),status='old',form='formatted',iostat=ierr)
       read(lu1,*,iostat=ierr)

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
!       print*,'Loop start at',stepstart,'/',datfull(1,1),'finish at',stepend,'/',thold

       ! Resize datfull2 here
       allocate(temp_array(ncols,stepend))
       temp_array(:,1:stepstart) = datfull2(:,:)
       deallocate(datfull2)
       ! Reassign all values to datfull2
       allocate(datfull2(ncols,1:stepend))
       datfull2(:,1:stepend) = temp_array(:,:)
       datfull2(:,stepstart:stepend) = datfull(:,1:(append_length+1))
       deallocate(temp_array)
       close(unit=lu1)
       deallocate(datfull)
    endif
    stepstart=stepend + 1
 enddo over_args

 ! Now calculate the accretion rates
 isteps = 0
 tprev = 0
 maccprev = 0.
 maccprev1 = 0.
 maccprev2 = 0.
 dat(:) = 0.
 ierr = 0

 do while(isteps  <=  size(datfull2,2))
    isteps = isteps + 1
    macc = datfull2(imacc_col,isteps) !dat(imacc_col)
    macc1 = datfull2(imacc1_col,isteps)
    macc2 = datfull2(imacc2_col,isteps)
    dt = datfull2(1,isteps) - tprev
    if (dt > tiny(dt)) then
       mdot = max(0., (macc - maccprev)/dt)
       mdot1 = max(0., (macc1 - maccprev1)/dt)
       mdot2 = max(0., (macc2 - maccprev2)/dt)
    else
       mdot = 0.
       mdot1 = 0.
       mdot2= 0.
    endif
    if (dt > dtmin .and. ierr==0) then
       !print *,' time ',dat(1),' mass accreted = ',macc, ' mdot = ',mdot
       write(iout,'(7(es18.10,1x))') datfull2(1,isteps),mdot,macc,mdot1,macc1,mdot2,macc2
       maccprev = macc
       maccprev1 = macc1
       maccprev2 = macc2
       tprev = dat(1)
    endif
 enddo

 close(unit=iout)

end program get_mdot
