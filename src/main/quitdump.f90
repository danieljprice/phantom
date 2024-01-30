!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module quitdump
!
! This module handles graceful exits
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, io_summary, part, readwrite_dumps, timestep
!
 implicit none
 public :: quit

 private

contains

subroutine quit
 use io,              only:iprint,ievfile,iscfile,iskfile,die
 use readwrite_dumps, only:write_fulldump
 use timestep,        only:time
 use io_summary,      only:summary_printout
 use part,            only:nptmass
 integer  :: i

!$omp critical (crit_quit_summary)
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
!$omp end critical (crit_quit_summary)
!
! close ev, log& ptmass-related files
 close(unit=ievfile)
 close(unit=iprint)
 if (iscfile > 0) close(unit=iscfile)
 do i = 1,nptmass
    close(unit=iskfile+i)
 enddo
 call die

end subroutine quit

end module quitdump
