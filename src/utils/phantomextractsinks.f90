!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantomextractsinks
!
!  DESCRIPTION: The will create one data file per sink from the master sink
!               file generated at runtime.  In the event that there are
!               multiple files per sink, we will use the data from the most
!               recent file
!
!  REFERENCES: None
!
!  OWNER: James Wurster
!
!  $Id$
!
!  USAGE: phantomextractsinks [no arguments]
!
!  DEPENDENCIES: prompting
!+
!--------------------------------------------------------------------------
program phantomextractsinks
 integer, parameter :: maxfiles = 100
 integer            :: i,io,nargs,nfiles,isink,nsink,nsink0
 real               :: tnow,tstart(maxfiles+1),asink(17)
 logical            :: fexists
 character(len= 32) :: prefix,filename
 character(len=512) :: cdummy,header
 !
 !--If argument exists, read it in
 !
 nargs = command_argument_count()
 if (nargs > 0) then
    call get_command_argument(1,prefix)
 else
    print*, 'Enter the model prefix: '
    read(5,"(a)") prefix
 endif
 !
 !--Determine the number of files and the first time in each file
 !
 fexists = .true.
 nfiles  = 0
 do while (fexists)
    nfiles = nfiles + 1
    write(filename,'(a,I2.2,a)') trim(prefix),nfiles,'.sink'
    inquire(file=filename,exist=fexists)
 enddo
 nfiles = nfiles - 1
 if (nfiles > 0) then
    print*, 'There are ',nfiles,' sink files from which to extract data'
 else
    print*, 'There are no sink files.  Good-bye'
    stop
 endif
 !
 !--Determine the first time in each file
 !
 tstart = huge(tstart)
 do i = 1,nfiles
    write(filename,'(a,I2.2,a)') trim(prefix),i,'.sink'
    open(unit=20,file=trim(filename))
    read(20,'(a)') cdummy
    read(20,'(a)') header
    read(20,*) tstart(i),isink,asink,nsink
    close(20)
 enddo
 !
 !--Read through the master sink files and write the data to the individual files
 !
 nsink0 = 0
 do i = 1,nfiles
    write(filename,'(a,I2.2,a)') trim(prefix),i,'.sink'
    open(unit=20,file=trim(filename))
    read(20,'(a)') cdummy
    read(20,'(a)') cdummy
    io = 0
    do while (io==0)
       read(20,*,iostat=io) tnow,isink,asink,nsink
       if (io == 0) then
          if (isink > nsink0) then
             nsink0 = isink
             write(filename,'(2a,I4.4,a)') trim(prefix),'_sink',isink,'.dat'
             open(unit=20+isink,file=trim(filename))
             write(20+isink,'(a)') trim(header)
          endif
          if (tnow < tstart(i+1)) write(20+isink,"(1pe18.9,1x,I18,1x,17(1pe18.9,1x),I18)") tnow,isink,asink,nsink
       endif
    enddo
    close(20)
 enddo
 print*, 'Completed extracting sink information.  The new data files are:'
 do i = 1,nsink0
    write(filename,'(2a,I4.4,a)') trim(prefix),'_sink',i,'.dat'
    print*, '  ',trim(filename)
    close(20+i)
 enddo

!--------------------------------------------------------------------------
end program phantomextractsinks
