!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: get_ang
!
!  DESCRIPTION: Utility to extract accreted ang mom as a function of time
!    from the accretion.ev file
!    Chris Nixon, March 2014
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: acc2ang sink_ID dt_min accretion.ev
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
program get_ang
 implicit none
 integer, parameter :: lu = 11, iout = 12
 integer, parameter :: imacc_col = 17
 integer, parameter :: ncols = imacc_col
 real :: dat(ncols)
 character(len=40) :: filename,outfile
 integer :: i,nargs,ierr,nlines
 real :: dtmin
 real :: ang(3),pos(3),vel(3),pos_sink(3)
 real :: time,msink,mu,rcirc,G,mpart
 integer :: isink,sink_ID
 integer :: npart,ncount
 real :: added_ang(3)

! assuming
 G=1.0

 nargs = command_argument_count()
 call get_command_argument(0,filename)

 if (nargs <= 0) then
    print "(a)",' Usage: '//trim(filename)//' sink_ID dt_min accretion.ev'
    stop
 endif

 call get_command_argument(1,filename)
 read(filename,*,iostat=ierr) sink_ID
 if (ierr /= 0) then
    print "(a)",' Usage: '//trim(filename)//' sink_ID dt_min accretion.ev'
    stop
 endif

 call get_command_argument(2,filename)
 read(filename,*,iostat=ierr) dtmin
 if (ierr /= 0) then
    print "(a)",' Usage: '//trim(filename)//' sink_ID dt_min accretion.ev'
    stop
 endif

 call get_command_argument(3,filename)
 open(unit=lu,file=trim(filename),status='old',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)//'...'
    stop
 endif

! open output file
 outfile = 'accretion.ang'
 print "(a)",trim(filename)//' --> '//trim(outfile)
 open(unit=iout,file=trim(outfile),status='replace',form='formatted')
 write(iout,'("Angular momentum accreted with time")')
 write(iout,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
      1,'time', &
      2,'npart', &
      3,'angmx', &
      4,'angmy', &
      5,'angmz', &
      6,'totang', &
      7,'rcirc'

! skip 2 header lines and get number of accreted particles
 read(lu,*)
 read(lu,*)

 ierr = 0
 nlines = 0
 do while(ierr == 0)
    nlines = nlines + 1
    read(lu,*,iostat=ierr) dat(1:ncols)
 enddo
 nlines=nlines-1
 write(*,'(i10," accreted particles found in file")') nlines

! Close and reopen file
 close(unit=lu)
 open(unit=lu,file=trim(filename),status='old',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)//'... for a second time.'
    stop
 endif
 read(lu,*) ! skip first comment
 read(lu,*) ! skip second comment

 ncount = 1
 npart = 0
 added_ang(:) = 0.0

 do i=1,nlines
    read(lu,*) time,isink,ang(1),ang(2),ang(3),mpart,msink,mu,&
               pos(1),pos(2),pos(3),vel(1),vel(2),vel(3),&
               pos_sink(1),pos_sink(2),pos_sink(3)

    if (isink /= sink_ID) cycle

    if (time > real(ncount)*dtmin) then
       if (npart > 0) then
          added_ang(:) = added_ang(:)/real(npart)
       endif
       rcirc = dot_product(added_ang,added_ang)/(G*msink*mpart*mpart)
       write(iout,'(7(es18.10,1x))') real(ncount)*dtmin,real(npart),&
             added_ang(1),added_ang(2),added_ang(3),&
             sqrt(dot_product(added_ang,added_ang)),rcirc
       ncount = ncount+1
       npart = 0
       added_ang(:) = 0.0
    endif

    npart = npart+1
    added_ang(:) = added_ang(:) + ang(:)
 enddo

 close(unit=iout)
 close(unit=lu)

 write(*,'("Finshed conversion")')

end program get_ang
