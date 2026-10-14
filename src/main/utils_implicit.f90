!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module implicit
!
! Utility routines for implicit radiative diffusion
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: io, neighcompact
!
 implicit none
 integer              :: ncompact,ncompactlocal,icompactmax,nneigh_average
 real,    allocatable :: vari(:,:),varij(:,:),varij2(:,:),varinew(:,:)
 integer, allocatable :: ivar(:,:),ijvar(:)
 logical, allocatable :: mask(:)
 logical              :: done_allocation = .false.

contains

!---------------------------------------------------------
!+
!  allocate arrays for compacted neighbour lists
!+
!---------------------------------------------------------
subroutine allocate_memory_implicit(npart,radkern,hfact,ierr)
 use io,           only:fatal
 use neighcompact, only:compact_neigh_counts
 integer, intent(in)  :: npart
 real,    intent(in)  :: radkern,hfact
 integer, intent(out) :: ierr

 if (done_allocation) then
    return
 else
    call compact_neigh_counts(npart,radkern,hfact,nneigh_average,icompactmax)
    allocate(ivar(3,npart),stat=ierr)
    if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for ivar')
    allocate(ijvar(icompactmax),stat=ierr)
    if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for ijvar')
    allocate(vari(2,npart),varij(2,icompactmax),varij2(4,icompactmax),varinew(3,npart),mask(npart),stat=ierr)
    if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for vari, varij, varij2, varinew, mask')
    done_allocation = .true.
 endif

end subroutine allocate_memory_implicit

end module implicit
