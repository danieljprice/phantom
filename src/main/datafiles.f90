!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: datafiles
!
!  DESCRIPTION:
!   Interface to routine to search for external data files
!   This module just provides the url and environment variable
!   settings that are specific to Phantom
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: datautils
!+
!--------------------------------------------------------------------------
module datafiles
 implicit none
 character(len=*), parameter :: data_url = &
  'http://users.monash.edu.au/~dprice/phantom/'

contains

function find_phantom_datafile(filename,loc)
 use datautils, only:find_datafile
 character(len=*), intent(in) :: filename,loc
 character(len=120) :: search_dir
 character(len=120) :: find_phantom_datafile

 search_dir = 'data/'//trim(adjustl(loc))
 find_phantom_datafile = find_datafile(filename,dir=search_dir,env_var='PHANTOM_DIR',url=data_url)

end function find_phantom_datafile

end module datafiles
