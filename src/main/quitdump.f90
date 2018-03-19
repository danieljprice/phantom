!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: quitdump
!
!  DESCRIPTION:
! This module handles graceful exits
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, io_summary, part, readwrite_dumps, timestep
!+
!--------------------------------------------------------------------------
module quitdump
 implicit none
 public :: quit

 private

contains

subroutine quit
 use io,              only:iprint,ievfile,iscfile,ipafile,iskfile,die
 use readwrite_dumps, only:write_fulldump
 use timestep,        only:time
 use io_summary,      only:summary_printout
 use part,            only:nptmass
 integer  :: i

!$omp critical
 call summary_printout(iprint,nptmass)
 write(iprint,10)
10 format(/, &
   '  _____ _____ _  __  _ ',/, &
   ' | ____| ____| |/ / | |',/, &
   ' |  _| |  _| | '' /  | |',/, &
   ' | |___| |___| . \  |_|',/, &
   ' |_____|_____|_|\_\ (_)',/)

 write(iprint,*) 'run terminated abnormally'
 call write_fulldump(time,'phantom.debugdump')
!$omp end critical
!
! close ev, log& ptmass-related files
 close(unit=ievfile)
 close(unit=iprint)
 if (iscfile > 0) close(unit=iscfile)
 if (ipafile > 0) close(unit=ipafile)
 do i = 1,nptmass
    close(unit=iskfile+i)
 enddo
 call die

end subroutine quit

end module quitdump
