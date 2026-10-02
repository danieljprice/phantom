!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Reads in sphNG data without B-fields, spits out phantom data with B-fields
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: kernel, moddump_utils, part, setup_params
!
 use moddump_utils, only:init_moddump=>init_moddump_empty, &
                         read_moddump=>read_moddump_empty,write_moddump=>write_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use setup_params, only:ihavesetupB
 use part,         only:hfact,mhd
 use kernel,       only:hfact_default
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 ! Define hfact
 hfact = hfact_default

 print*,'sphNG data reformatted for phantom write'

 ! Now add B field
 if (mhd) then
    ihavesetupB = .false.
 endif

end subroutine modify_dump

end module moddump
