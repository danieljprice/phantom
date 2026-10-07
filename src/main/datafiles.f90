!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! MESA EOS table version (hidden variable, set via input file)
 integer, public :: eosmesa_version = 2   ! 0 or 1 uses the older version of the mesa tables (Reichardt et al 2020), which are no longer used by default

 ! GitHub mirror of data files (backup if Zenodo is unreachable)
 character(len=*), parameter :: mirror_raw_base = &
    'https://raw.githubusercontent.com/phantomSPH/phantom-datafiles/main/'
 character(len=*), parameter :: mirror_release_base = &
    'https://github.com/phantomSPH/phantom-datafiles/releases/download/large-files/'

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
                            url=map_dir_to_web(trim(search_dir)),&
                            url_fallback=map_dir_to_mirror(trim(search_dir),trim(filename)))
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

 select case(search_dir)
 case('data/eos/mesa')
    ! EOS table versions:
    !   0,1 = legacy tables (Reichardt et al. 2020)
    !   2   = current tables
    select case(eosmesa_version)
    case(0,1)
       url = 'https://zenodo.org/records/13148447/files/'
    case(2)
       url = 'https://zenodo.org/records/21712204/files/'
    case default
       stop 'Unknown eosmesa_version'
    end select
 case('data/eos/mesa_opac') !!!  same link as old eos/mesa, but this is used only to download the opacity tables
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
 case('data/binarybh')
    url = 'https://zenodo.org/records/18615172/files/'
 case('data/eos/lombardi')
    url = 'https://zenodo.org/records/13842491/files/'
 case('data/star_data_files')
    url = 'https://zenodo.org/records/20738843/files/'
 case default
    url = 'https://users.monash.edu.au/~dprice/'//trim(search_dir)
 end select

end function map_dir_to_web

!----------------------------------------------------------------
!+
!  Fallback URL on the phantom-datafiles GitHub mirror.
!  Large files (>=100 MB) are stored as Release assets; others
!  are fetched from raw.githubusercontent.com under data/.
!  The returned string is a directory/prefix; retrieve_remote_file
!  appends the filename (same convention as map_dir_to_web).
!+
!----------------------------------------------------------------
function map_dir_to_mirror(search_dir,filename) result(url)
 character(len=*), intent(in) :: search_dir,filename
 character(len=200) :: url

 if (is_large_mirror_file(filename)) then
    ! release assets sit at the release root; prefix is the download base
    url = trim(mirror_release_base)
 else
    url = trim(mirror_raw_base)//trim(search_dir)//'/'
 endif

end function map_dir_to_mirror

!----------------------------------------------------------------
!+
!  files too large for normal git blobs; hosted as Release assets
!+
!----------------------------------------------------------------
logical function is_large_mirror_file(filename)
 character(len=*), intent(in) :: filename

 select case(trim(filename))
 case('galaxiesP25e5.dat','eos_binary_table.dat')
    is_large_mirror_file = .true.
 case default
    is_large_mirror_file = .false.
 end select

end function is_large_mirror_file

end module datafiles
