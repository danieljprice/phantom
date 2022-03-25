!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Recalculate internal energy and temperature of all gas particles using a
! specified EoS.
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, part, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos,   only:equationofstate,ieos,init_eos,done_init_eos,calc_temp_and_ene,finish_eos,&
                 gmw,X_in,Z_in,irecomb,gamma,eosinfo
 use io,    only:iprint
 use part,  only:rhoh,eos_vars,itemp,igasP,igas,store_temperature
 use units, only:unit_density,unit_pressure,unit_ergg
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,ierr
 real :: densi,eni,tempi,ponrhoi,temp_residual,ene_residual
 real :: dum1,dum2,dum3,dum4

 !-SET-EOS-OF-INPUT-DUMP--------
 ieos = 12
 gamma = 5./3.
 gmw = 0.60319
 irecomb = 0
 X_in = 0.69843  ! Set X and Z. Only relevant for ieos = 10, 20
 Z_in = 0.01426
 !-------------------------------

 write(iprint,"(/,a,i2)") 'Assuming input dump has ieos = ',ieos
 if (.not. done_init_eos) call init_eos(ieos,ierr)
 call eosinfo(ieos,iprint)
 print*,'Check if this is correct. Enter to proceed'
 read*

 dum1 = 0.
 dum2 = 0.
 dum3 = 0.
 dum4 = 0.

 do i = 1,npart
    densi = rhoh(xyzh(4,i),massoftype(igas))

    ! Get pressure
    call equationofstate(ieos,ponrhoi,dum1,densi,dum2,dum3,dum4,vxyzu(4,i))
    eos_vars(igasP,i) = ponrhoi * densi
 enddo

 !-SET-EOS-OF-OUTPUT-DUMP--------
 ! Comment out to leave quantity unchanged
 ieos = 20
 gamma = 5./3.
 gmw = 0.60319
 irecomb = 3
 X_in = 0.69843  ! Set X and Z. Only relevant for ieos = 10, 20
 Z_in = 0.01426
 !-------------------------------

 write(iprint,"(/,a,i2)")'Changing to ieos = ',ieos
 call init_eos(ieos,ierr)
 call eosinfo(ieos,iprint)
 print*,'Check if this is correct. Enter to proceed'
 read*

 !!!! BEGIN CHECKING COMPOSITION
!  tempi = 0.
!  do i = 1,npart
!     densi = rhoh(xyzh(4,i),massoftype(igas))
!     call calc_temp_and_ene(20,densi*unit_density,eos_vars(igasP,i)*unit_pressure,eni,tempi,ierr)
!     temp_residual = abs( tempi / eos_vars(itemp,i) - 1. )  ! (Tnew-Told)/Told
!     ene_residual = abs( eni / unit_ergg / vxyzu(4,i) - 1. )
!     print*, 'i = ',i,' dT/T = ',temp_residual,' du/u = ',ene_residual
!  enddo
!  stop
!!!! END CHECKING COMPOSITION

 tempi = 0.
 do i = 1,npart
    densi = rhoh(xyzh(4,i),massoftype(igas))
    call calc_temp_and_ene(ieos,densi*unit_density,eos_vars(igasP,i)*unit_pressure,eni,tempi,ierr)
    vxyzu(4,i) = eni / unit_ergg
    if (store_temperature) eos_vars(itemp,i) = tempi
 enddo

 call finish_eos(ieos,ierr)

end subroutine modify_dump

end module moddump

