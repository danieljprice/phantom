!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module torques
!
!  Writes diagnostics needed for analysis of torques, written
!  to compute necessary outputs for KITP binary-disc comparison
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, io, part, physcon
!
 implicit none
 public init_torques, write_torques
 integer, private :: itorque_unit
 private

contains

subroutine init_torques(evfile)
 character (len=*), intent(in) :: evfile
 character (len=len(evfile))   :: filename
 integer :: idot

 idot = index(evfile,'.ev') - 1
 filename = evfile(1:idot)//'.tq'  ! create .tq file name

 open(newunit=itorque_unit,file=filename,form='formatted',status='replace')

 !---Write header
 write(itorque_unit,"('#',9(1x,'[',i2.2,1x,a13,']'))") &
       1,'time', &
       2,'Tz', &
       3,'Tz(R>a)', &
       4,'PsiR', &
       5,'PsiI', &
       6,'<e_x>', &
       7,'<e_y>', &
       8,'Mdot1', &
       9,'Mdot2'

end subroutine init_torques

subroutine write_torques(time)
 use part, only:npart,xyzh,vxyzu,massoftype,fext,isdead_or_accreted
 use io,   only:fatal
 use vectorutils, only:cross_product3D
 use physcon,     only:twopi
 real, intent(in) :: time
 reaL :: pmassi,torque(3),torque_excised(3),rcrossa(3),r2
 real :: ri(3),vi(3),ecci(3),ecc(3),rhat(3),psix,psiy
 integer :: i

 !
 ! compute total torque
 !
 pmassi=massoftype(1)
 torque = 0.
 torque_excised = 0.
 psix = 0.
 psiy = 0.
 ecc = 0.
 !$omp parallel do default(none) &
 !$omp shared(xyzh,vxyzu,fext,npart,pmassi) &
 !$omp private(i,ri,rcrossa,r2,ecci,rhat,vi)  &
 !$omp reduction(+:torque,torque_excised,psix,psiy,ecc)
 do i = 1,npart
    ri = xyzh(1:3,i)
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       r2 = dot_product(ri(1:2),ri(1:2))
       call cross_product3D(ri,fext(1:3,i),rcrossa)
       torque(1:3) = torque(1:3) + pmassi*rcrossa(1:3)
       if (r2 > 1.) then
          torque_excised(1:3) = torque_excised(1:3) + pmassi*rcrossa(1:3)
          if (r2 < 36) then ! between r=1a and r=6a
             vi(:) = vxyzu(1:3,i)
             rhat(:) = ri(:)/sqrt(r2)
             ecci(:) = (dot_product(vi,vi)*ri - dot_product(vi,ri)*vi) - rhat
             ecc = ecc + pmassi*ecci
          endif
       endif
       psix = psix + pmassi*ri(1)
       psiy = psiy + pmassi*ri(2)
    endif
 enddo
 !$omp end parallel do

 !
 ! write the torques to file
 !
 ! This list is subject to change as we come up with new diagnostics.  At the moment the data file should have nine columns:
 ! 1   Time (in binary orbits, so be sure to divide by 2pi)
 ! 2   Total Gravitational Torque (r x F_grav).  Normalize to code units in overleaf.
 ! 3   Grav. Torque outside r>a
 ! 4   PsiR  (Mass Dipole Moment, eq 16 in overleaf)
 ! 5   PsiI  (17)
 ! 6   <e_x> (Eccentricity Vector, eq 23)
 ! 7   <e_y>
 ! 8   Mdot onto Primary
 ! 9   Mdot onto Secondary
 !
 write(itorque_unit,"(1x,9(1x,1pe18.10))") time/twopi,-torque(3),-torque_excised(3),psix,psiy,ecc(1),ecc(2),0.,0.
 call flush(itorque_unit)

end subroutine write_torques

end module torques
