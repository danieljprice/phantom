!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
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
! :Dependencies: None
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos, only:equationofstate,ieos,init_eos,done_init_eos,calc_temp_and_ene,finish_eos,gmw,X_in,Z_in
 use part, only:rhoh,eos_vars,itemp,igasP,igas,store_temperature
 use units, only:unit_density,unit_pressure,unit_ergg
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,ierr
 real :: densi,eni,tempi,ponrhoi
 real :: dum,dum2

 !-SET-EOS-OF-INPUT-DUMP--------
 ieos = 12
 gmw = 0.61821
 !--------------------------
 print*,'Assuming input dump has ieos = ',ieos
 print*,'Assuming input dump has gmw = ',gmw
 if (.not. done_init_eos) call init_eos(ieos,ierr)

 dum = 0.0
 do i = 1,npart
    densi = rhoh(xyzh(4,i),massoftype(igas))

    ! Get pressure
    call equationofstate(ieos,ponrhoi,dum2,densi,dum,dum,dum,vxyzu(4,i))
    eos_vars(igasP,i) = ponrhoi * densi
 enddo

 !-SET-EOS-OF-OUTPUT-DUMP--------
 ieos = 10
 if (ieos == 10) then
    X_in = 0.69843
    Z_in = 0.01426
 endif 
 !--------------------------
 call init_eos(ieos,ierr)
 print*,'Changing to ieos = ',ieos
 do i = 1,npart
    densi = rhoh(xyzh(4,i),massoftype(igas))
    call calc_temp_and_ene(densi*unit_density,eos_vars(igasP,i)*unit_pressure,eni,tempi,ierr)
    vxyzu(4,i) = eni / unit_ergg
    if (store_temperature) eos_vars(itemp,i) = tempi
 enddo

 call finish_eos(ieos,ierr)

end subroutine modify_dump

end module moddump

