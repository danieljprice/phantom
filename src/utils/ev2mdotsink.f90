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
!  USAGE: ev2mdotsink [dt_min] file1.ev file2.ev ...
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
program get_mdot
 implicit none
 integer, parameter :: lu1 = 11, lu2 = 13, iout = 12
! integer, parameter :: nsmooth = 20
 integer, parameter :: imacc_col = 12
 integer, parameter :: ncols = imacc_col
 real :: dat(ncols),dat2(ncols-1)
 real, dimension(:,:), allocatable :: datfull,datfull2,temp_array
 character(len=40) :: filename,outfile
 integer :: i,nargs,ierr,nlines,iev,isteps,istart,append_length,stepstart,ii,jj,stepend
 real :: maccprev,macc,tprev,dt,mdot,dtmin,thold

 allocate(datfull2(ncols,1))
 datfull2(:,:) = 0.

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
    write(iout,"('#',3(a17,2x))") 't','mdot','macc'
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

 ! Now calculate the accretion rates
 isteps = 0
 tprev = 0
 maccprev = 0.
 dat(:) = 0.
 ierr = 0

 do while(isteps  <=  (size(datfull2,2)-1))
    isteps = isteps + 1
    macc = datfull2(imacc_col,isteps) !dat(imacc_col)
    dt = datfull2(1,isteps) - tprev
    if (dt > tiny(dt)) then
       mdot = max(0., (macc - maccprev)/dt)
    else
       mdot = 0.
    endif
    if (dt > dtmin .and. ierr==0) then
       !print *,' time ',datfull2(1,isteps),' mass accreted = ',macc, ' mdot = ',mdot
       write(iout,"(3(es18.10,1x))") datfull2(1,isteps),mdot,macc
       tprev = datfull2(1,isteps)
       maccprev = macc
    endif
 enddo

 close(unit=iout)

end program get_mdot
