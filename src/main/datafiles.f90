!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module datafiles
!
! Interface to routine to search for external data files
! This module just provides the url and environment variable
! settings that are specific to Phantom
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

contains

!----------------------------------------------------------------
!+
!  Find a datafile in the Phantom data directory
!+
!----------------------------------------------------------------
function find_phantom_datafile(filename,loc)
 use datautils, only:find_datafile
 use io,        only:id,master
 use mpiutils,  only:barrier_mpi
 character(len=*), intent(in) :: filename,loc
 character(len=120) :: search_dir
 character(len=120) :: find_phantom_datafile

 search_dir = 'data/'//trim(adjustl(loc))
 if (id == master) then ! search for and download datafile if necessary
    find_phantom_datafile = find_datafile(filename,dir=search_dir,env_var='PHANTOM_DIR',&
                            url=map_dir_to_web(trim(search_dir)))
 endif
 call barrier_mpi()
 if (id /= master) then ! find datafile location, do not attempt to download it
    find_phantom_datafile = find_datafile(filename,dir=search_dir,&
                            env_var='PHANTOM_DIR',verbose=.false.)
 endif

end function find_phantom_datafile

!----------------------------------------------------------------
!+
!  Find the web location for files that are not in the Phantom
!  git repo, which need to be downloaded into the data directory
!  at runtime
!+
!----------------------------------------------------------------
function map_dir_to_web(search_dir) result(url)
 character(len=*), intent(in) :: search_dir
 character(len=120) :: url

 !print*,' search_dir=',trim(search_dir)
 select case(search_dir)
 case('data/eos/mesa')
    url = 'https://zenodo.org/records/13148447/files/'
 case('data/eos/shen')
    url = 'https://zenodo.org/records/13163155/files/'
 case('data/eos/helmholtz')
    url = 'https://zenodo.org/records/13163286/files/'
 case('data/forcing')
    url = 'https://zenodo.org/records/13162225/files/'
 case('data/velfield')
    url = 'https://zenodo.org/records/13162515/files/'
 case('data/galaxy_merger')
    url = 'https://zenodo.org/records/13162815/files/'
 case('data/starcluster')
    url = 'https://zenodo.org/records/13164858/files/'
 case('data/eos/lombardi')
    url = 'https://zenodo.org/records/13842491/files/'
 case default
    url = 'https://users.monash.edu.au/~dprice/'//trim(search_dir)
 end select
 !print*,'url=',trim(new_url)

end function map_dir_to_web

end module datafiles
