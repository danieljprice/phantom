!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine which computes neighbour lists for all particles
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: getneigbours
!+
!--------------------------------------------------------------------------
module analysis
 use getneigbours,    only:generate_neighbour_lists, read_neighbours, write_neighbours, &
                           neighcount,neighb,neighmax
 implicit none
 character(len=20), parameter, public :: analysistype = 'getneighbours'

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character(len=100)           :: neighbourfile

 !***************************************
 ! Assign neighbour lists to particles by searching shared list of host cell
 !***************************************

 call generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile,.true.)

 !**************************************
 ! Output neighbour lists to file
 !**************************************

 neighbourfile = 'neigh_'//TRIM(dumpfile)
 call write_neighbours(neighbourfile, npart)

 print*, 'Neighbour finding complete for file ', TRIM(dumpfile)

end subroutine do_analysis
!--------------------------------------------------------------------------
end module
