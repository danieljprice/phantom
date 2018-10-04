!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine to print the sink masses as functions of time.
!  Each sink particle will have its own file.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, part, prompting
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'sinkmass'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 use part,     only:xyzmh_ptmass, nptmass
 use prompting,  only:prompt
 use dim, only:maxptmass
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character(len=20) :: fname
 integer, save :: j = 0
 integer :: i
 logical :: does_it_exists, cont

 ! Only does this check the first time this subroutine is called.
 if (j==0) then
    over_sinks: do i=0,maxptmass
       write(fname,"('sink_',i5.5,'.dat')") i
       inquire(file=fname,exist=does_it_exists)
       if (does_it_exists) then
          print "(a)",'It appears that you already have files for sinkmasses.'
          print "(a)",'Data will be appended to these files.'
          cont = .false.
          call prompt('Do you wish to continue?',cont)
          if (cont) then
             exit over_sinks
          else
             stop 'Bailing out!!'
          endif
       endif
    enddo over_sinks
 endif

 ! Write the mass at a specific time to a different file for each sink.
 if (nptmass > 0) then
    do i = 1,nptmass
       write(fname,"('sink_',i5.5,'.dat')") i
       inquire(file=fname,exist=does_it_exists)
       if (does_it_exists) then
          open(iunit,file=fname,position='append')
       else
          open(iunit,file=fname,status='replace')
       endif
       write(iunit,*) time, xyzmh_ptmass(4,i)
       close(iunit)
    enddo
 endif

 j=1

end subroutine do_analysis

end module
