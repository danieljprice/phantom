!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! :Runtime parameters:
!   - Xfrac : *hydrogen mass fraction*
!   - Zfrac : *metal mass fraction*
!   - mu    : *mean molecular weight*
!
! :Dependencies: dim, eos, eos_idealplusrad, eos_mesa, io,
!   mesa_microphysics, moddump_utils, part, prompting, radiation_utils,
!   units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters (composition; gamma is fixed at 5/3 by ieos=12)
 real :: Xfrac = 0.687    ! hydrogen mass fraction
 real :: Zfrac = 0.0142   ! metal mass fraction
 real :: mu    = 0.61821  ! mean molecular weight

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use units,            only:unit_density,unit_opacity,unit_ergg
 use dim,              only:do_radiation
 use io,               only:fatal
 use eos,              only:gmw,gamma,X_in,Z_in
 use eos_idealplusrad, only:get_idealplusrad_temp
 use eos_mesa,         only:init_eos_mesa
 use part,             only:rho,igas,rad,iradxi,ikappa,radprop,ithick
 use radiation_utils,  only:radiation_and_gas_temperature_equal,ugas_from_Tgas
 use mesa_microphysics,only:get_kappa_mesa
 use moddump_utils,    only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i,ierr
 real                   :: pmass,rhoi,kappa_cgs,kapt,kapr,rho_cgs,ugasi,tempi,gamma_fixed

 if (.not. do_radiation) call fatal("moddump_LTE_to_rad","Not compiled with radiation")

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 X_in = Xfrac
 Z_in = Zfrac
 gmw  = mu
 gamma_fixed = 5/3.  ! gamma should be exactly 5/3, because that is what ieos=12 assumes
 gamma = gamma_fixed
 print*,'Assuming gmw = ',mu,' and gamma=',gamma,'X = ',X_in,'Z = ',Z_in  ! X and Z are only used for calculating opacity
 call init_eos_mesa(X_in,Z_in,ierr)

 pmass = massoftype(igas)
 do i=1,npart
    rhoi = rho(i)
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

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter hydrogen mass fraction X',Xfrac,0.,1.)
 call prompt('Enter metal mass fraction Z',Zfrac,0.,1.)
 call prompt('Enter mean molecular weight mu',mu,0.)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(Xfrac,'Xfrac',db,errcount=nerr,min=0.,max=1.)
 call read_inopt(Zfrac,'Zfrac',db,errcount=nerr,min=0.,max=1.)
 call read_inopt(mu,'mu',db,errcount=nerr,min=0.)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# composition for LTE -> radiation conversion'
 call write_inopt(Xfrac,'Xfrac','hydrogen mass fraction',iunit)
 call write_inopt(Zfrac,'Zfrac','metal mass fraction',iunit)
 call write_inopt(mu,'mu','mean molecular weight',iunit)

end subroutine write_options_moddump

end module moddump
