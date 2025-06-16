!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
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
! :Dependencies: eos, eos_gasradrec, io, part, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos,           only:get_pressure,ieos,init_eos,done_init_eos,calc_temp_and_ene,finish_eos,&
                         gmw,X_in,Z_in,gamma,eosinfo
 use eos_gasradrec, only:irecomb
 use io,            only:iprint
 use part,          only:rhoh,eos_vars,itemp,igasP,igas
 use units,         only:unit_density,unit_pressure,unit_ergg
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i,ierr
 real                   :: densi,eni,tempi

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

 do i = 1,npart
    densi = rhoh(xyzh(4,i),massoftype(igas))
    eos_vars(igasP,i) = get_pressure(ieos,xyzh(:,i),densi,vxyzu(:,i))
 enddo

 !-SET-EOS-OF-OUTPUT-DUMP--------
 ! Comment out to leave quantity unchanged
 ieos = 2
 gamma = 5./3.
 gmw = 0.6175
 irecomb = 3
 if (ieos == 10) then
    X_in = 0.69843  ! Set X and Z. Only relevant for ieos = 10, 20
    Z_in = 0.01426
 endif
 !-------------------------------

 write(iprint,"(/,a,i2)")'Changing to ieos = ',ieos
 call init_eos(ieos,ierr)
 call eosinfo(ieos,iprint)
 print*,'Check if this is correct. Enter to proceed'
 read*

 tempi = 0.
 do i = 1,npart
    densi = rhoh(xyzh(4,i),massoftype(igas))
    call calc_temp_and_ene(ieos,densi*unit_density,eos_vars(igasP,i)*unit_pressure,eni,tempi,ierr)
    vxyzu(4,i) = eni / unit_ergg
    eos_vars(itemp,i) = tempi
 enddo

 call finish_eos(ieos,ierr)

end subroutine modify_dump

end module moddump

