!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: domain
!
!  DESCRIPTION:
!  This module performs the MPI domain decomposition
!  this one based on the cell ID.
!
!  REFERENCES: None
!
!  OWNER: Owen Kaluza
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module domain
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 integer, public :: cellspp = 0
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

 public :: init_domains,cellbelong

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

!----------------------------------------------------------------
!+
!  Determines if a given cell belongs on a given processor
!  to given its index and proc id
!+
!----------------------------------------------------------------
pure logical function cellbelong(icell,id,nprocs,ncells)
 integer,         intent(in) :: id,nprocs
 integer,         intent(in) :: icell
 integer(kind=8), intent(in) :: ncells

 !Strided cells: distribute over procs alternately
 if (id == MOD(icell, nprocs)) then
   cellbelong = .true.
 else
   cellbelong = .false.
 endif

end function cellbelong

end module domain
