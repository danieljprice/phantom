!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: domain
!
!  DESCRIPTION:
!   This module performs the MPI domain decomposition
!   Since we now do the decomposition using the tree all this
!   module does is store the ibelong array for the particles
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, part
!+
!--------------------------------------------------------------------------
module domain
 use dim, only:maxp
 use io,  only:nprocs
 use part, only:ibelong
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 integer, parameter :: ndim = 3

 !--default initialisation of domainisperiodic so that
 !  in principle init_domains does not need to be called
 !  for non-MPI runs (safest to do so however)
 !
#ifdef PERIODIC
 logical, dimension(ndim), public :: isperiodic = .true.
#else
 logical, dimension(ndim), public :: isperiodic = .false.
#endif

 public :: init_domains, assign_to_domain
 public :: i_belong, i_belong_i4

! interface i_belong
!  module procedure i_belong_i4, i_belong
! end interface

 private

contains
!-----------------------------------------------------------------------
!+
!  initialise routine (does nothing for this domain decomposition)
!+
!-----------------------------------------------------------------------
subroutine init_domains(nprocs)
 integer, intent(in) :: nprocs

end subroutine init_domains

!-----------------------------------------------------------------------
!+
!  routine to initially distribute particles evenly between domains
!  i.e. during setup. The actual domain decomposition is then done
!  during the tree build
!+
!-----------------------------------------------------------------------
integer function assign_to_domain(i,id)
 integer(kind=8), intent(in) :: i
 integer,         intent(in) :: id

#ifdef MPI
 assign_to_domain = int(mod(i,int(nprocs,kind=8)),kind=kind(assign_to_domain))
#else
 assign_to_domain = id
#endif

end function assign_to_domain

logical function i_belong_i4(iparttot)
 use io, only:nprocs,id
 integer(kind=4), intent(in) :: iparttot

 i_belong_i4 = (id == int(mod(iparttot, int(nprocs, kind=kind(iparttot))), kind=kind(nprocs)) )
end function i_belong_i4

logical function i_belong(iparttot)
 use io, only:nprocs,id
 integer(kind=8), intent(in) :: iparttot

 i_belong = (id == int(mod(iparttot, int(nprocs, kind=kind(iparttot))), kind=kind(nprocs)) )

end function i_belong

end module domain
