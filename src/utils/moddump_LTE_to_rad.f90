!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Convert non-radiation dump (assuming LTE, ieos=12) to radiation dump
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, eos_idealplusrad, eos_mesa, io,
!   mesa_microphysics, part, radiation_utils, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use units,            only:unit_density,unit_opacity,unit_ergg
 use dim,              only:do_radiation
 use io,               only:fatal
 use eos,              only:gmw,gamma,X_in,Z_in
 use eos_idealplusrad, only:get_idealplusrad_temp
 use eos_mesa,         only:init_eos_mesa
 use part,             only:igas,rad,iradxi,ikappa,rhoh,radprop,ithick
 use radiation_utils,  only:radiation_and_gas_temperature_equal,ugas_from_Tgas
 use mesa_microphysics,only:get_kappa_mesa
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i,ierr
 real                   :: pmass,mu,rhoi,kappa_cgs,kapt,kapr,rho_cgs,ugasi,tempi,gamma_fixed

 if (.not. do_radiation) call fatal("moddump_LTE_to_rad","Not compiled with radiation")

 X_in=0.687
 Z_in=0.0142
 mu = 0.61821
 gmw = mu
 gamma_fixed = 5/3.  ! gamma should be exactly 5/3, because that is what ieos=12 assumes
 gamma = gamma_fixed
 print*,'Assuming gmw = ',mu,' and gamma=',gamma,'X = ',X_in,'Z = ',Z_in  ! X and Z are only used for calculating opacity
 call init_eos_mesa(X_in,Z_in,ierr)

 pmass = massoftype(igas)
 do i=1,npart
    rhoi = rhoh(xyzh(4,i),pmass)
    rho_cgs = rhoi*unit_density
    call get_idealplusrad_temp(rho_cgs,vxyzu(4,i)*unit_ergg,mu,tempi,ierr)

    ! calculate u and xi
    ugasi = ugas_from_Tgas(tempi,gamma,mu)
    vxyzu(4,i) = ugasi
    rad(iradxi,i) = radiation_and_gas_temperature_equal(rhoi,ugasi,gamma,mu)

    ! calculate opacity
    call get_kappa_mesa(rho_cgs,tempi,kappa_cgs,kapt,kapr)
    radprop(ikappa,i) = kappa_cgs/unit_opacity
    radprop(ithick,i) = 1.
 enddo

end subroutine modify_dump

end module moddump

