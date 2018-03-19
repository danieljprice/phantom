!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setvfield
!
!  DESCRIPTION:
! this module contains utilities for setting up velocity fields
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: energies, io, mpiutils
!+
!--------------------------------------------------------------------------
module setvfield
 implicit none
 public :: set_vfield,normalise_vfield

 integer, parameter, public :: &
    iuniform_rotation_z = 1

 private

contains

!-------------------------------------------------------------------
!+
!  utility to set up various velocity fields
!+
!------------------------------------------------------------------
subroutine set_vfield(itype,angvel,npart,npartoftype,xyzh,massoftype,vxyzu)
 real,    intent(in)  :: angvel
 integer, intent(in)  :: itype,npart
 integer, intent(in)  :: npartoftype(:)
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: vxyzu(:,:)
 real,    intent(in)  :: massoftype(:)
 integer :: i

 select case(itype)
 case(iuniform_rotation_z)
    do i=1,npart
       vxyzu(1,i) = -angvel*xyzh(2,i)
       vxyzu(2,i) =  angvel*xyzh(1,i)
       vxyzu(3,i) = 0.
    enddo
 case default
    vxyzu(1:3,1:npart) = 0.
 end select

 return
end subroutine set_vfield

!-------------------------------------------------------------------
!+
!  utility to normalise the velocity field according to either the:
!  1) total kinetic energy
!  2) RMS velocity
!  3) RMS Mach number
!+
!------------------------------------------------------------------
subroutine normalise_vfield(npart,vxyzu,ierr,rms,ke,rmsmach)
 use io,       only:fatal,id,master
 use energies, only:compute_energies,rmsmachn=>rmsmach,vrms,ekin
 use mpiutils, only:bcast_mpi
 integer, intent(in)    :: npart
 real,    intent(inout) :: vxyzu(:,:)
 integer, intent(out)   :: ierr
 real,    intent(in), optional :: rms,ke,rmsmach
 real :: factor
 integer :: i
!
!--compute the total kinetic energy, RMS Mach number and RMS velocity
!
 call compute_energies(0.)
!
!--make sure these are same on all processors
!
 call bcast_mpi(ekin)
 call bcast_mpi(vrms)
 call bcast_mpi(rmsmachn)
!
!--work out factor needed to renormalise velocity field to correct value
!
 factor = -1.
 if (present(rms)) then
    if (vrms > 0.) factor = rms/vrms
    if (present(ke) .or. present(rmsmach)) ierr = 1
 elseif (present(ke)) then
    if (ekin > 0.) factor = sqrt(ke/ekin)
    if (present(rms) .or. present(rmsmach)) ierr = 1
 elseif (present(rmsmach)) then
    if (rmsmachn > 0.) factor = rmsmach/rmsmachn
    if (present(rms) .or. present(ke)) ierr = 1
 else
    factor = 0.
    ierr = 1
 endif
 if (factor < 0) ierr = 2
 if (ierr /= 0) return
!
!--perform renormalisation of velocity field
!
 do i=1,npart
    vxyzu(1:3,i) = vxyzu(1:3,i)*factor
 enddo
!
!--check this has worked... (not strictly necessary)
!
 call compute_energies(0.)
 if (id==master) write(*,"(1x,3(a,es10.3))") &
    'Normalising v field: Ekin = ',ekin,' rms v = ',vrms,' rms mach = ',rmsmachn

end subroutine normalise_vfield

end module setvfield
