!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module datafiles
!
! Interface to routine to search for external data files
!   This module just provides the url and environment variable
!   settings that are specific to Phantom
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: datautils, io, mpiutils
!
 implicit none
 character(len=*), parameter :: data_url = &
  'https://users.monash.edu.au/~dprice/phantom/'

contains

function find_phantom_datafile(filename,loc)
 use datautils, only:find_datafile
 use io,        only:id,master
 use mpiutils,  only:barrier_mpi
 character(len=*), intent(in) :: filename,loc
 character(len=120) :: search_dir
 character(len=120) :: find_phantom_datafile

 search_dir = 'data/'//trim(adjustl(loc))
 if (id == master) then ! search for and download datafile if necessary
    find_phantom_datafile = find_datafile(filename,dir=search_dir,env_var='PHANTOM_DIR',url=data_url)
 endif
 call barrier_mpi()
 if (id /= master) then ! find datafile location, do not attempt to download it
    find_phantom_datafile = find_datafile(filename,dir=search_dir,&
                            env_var='PHANTOM_DIR',verbose=.false.)
 endif

end function find_phantom_datafile

end module datafiles
