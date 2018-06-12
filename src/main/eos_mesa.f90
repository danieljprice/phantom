!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: eos_mesa
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: mesa_microphysics
!+
!--------------------------------------------------------------------------
module eos_mesa

 use mesa_microphysics

 implicit none

contains

!----------------------------------------------------------------
!+
!  subroutine initialises the mesa eos tables
!+
!----------------------------------------------------------------
subroutine init_eos_mesa(x,z,ierr)
 real, intent(in) :: x,z
 integer, intent(out) :: ierr

 ierr=0
 if ((x + z > 1.) .or. (x < 0.) .or. (z < 0.)) then
    ierr=-1
    return
 endif

 call get_environment_variable('MESA_DATA_DIR',mesa_eos_dir)
 mesa_eos_prefix="output_DE_"
 call get_environment_variable('MESA_DATA_DIR',mesa_opacs_dir)
 mesa_opacs_suffix=""
 mesa_eos_full_output=.true.
 mesa_binary_data=.true.

!Initialize the mesa_microphysics module
 mesa_eos_z1=0.d0   ! lowest Z
 mesa_eos_h1=0.d0   ! lowest H
 mesa_eos_dz=0.02d0 ! delta Z
 mesa_eos_dh=0.2d0  ! delta H in EOS tables

 call get_eos_constants_mesa(ierr)
 if (ierr /= 0) return

 call read_eos_mesa(x,z,ierr)
 call get_opacity_constants_mesa
 call read_opacity_mesa(x,z)

end subroutine init_eos_mesa

!----------------------------------------------------------------
!+
!  subroutine closes the mesa eos tables
!+
!----------------------------------------------------------------
subroutine finish_eos_mesa

 call deallocate_arrays_mesa

end subroutine finish_eos_mesa

!----------------------------------------------------------------
!+
!  subroutine returns pressure/gamma1 as a function of
!  density/internal energy
!+
!----------------------------------------------------------------
subroutine get_eos_pressure_gamma1_mesa(den,eint,pres,gam1)
 real, intent(in) :: den, eint
 real, intent(out) :: pres, gam1

 call getvalue_mesa(den,eint,2,pres)
 call getvalue_mesa(den,eint,11,gam1)

end subroutine get_eos_pressure_gamma1_mesa

!----------------------------------------------------------------
!+
!  subroutine returns kappa as a function of
!  density, temperature and composition
!+
!----------------------------------------------------------------
subroutine get_eos_kappa_mesa(den,temp,kappa,kappat,kappar)
 real, intent(in)  :: den, temp
 real, intent(out) :: kappa, kappat, kappar

 call get_kappa_mesa(den,temp,kappa,kappat,kappar)

end subroutine get_eos_kappa_mesa

!----------------------------------------------------------------
!+
!  subroutine returns pressure/temperature as
!  a function of density/internal energy
!+
!----------------------------------------------------------------
subroutine get_eos_pressure_temp_mesa(den,eint,pres,temp)
 real, intent(in) :: den, eint
 real, intent(out) :: pres, temp

 call getvalue_mesa(den,eint,2,pres)
 call getvalue_mesa(den,eint,4,temp)

end subroutine get_eos_pressure_temp_mesa

!----------------------------------------------------------------
!+
!  subroutine returns various quantities as
!  a function of density/internal energy
!+
!----------------------------------------------------------------
subroutine get_eos_various_mesa(den,eint,pres,proint,peint,temp,troint,teint,entrop,abad,gamma1,gam)
 real, intent(in) :: den, eint
 real, intent(out) :: pres, temp, proint, peint, troint, teint, entrop, abad, gamma1, gam

 call getvalue_mesa(den,eint,2,pres)
 call getvalue_mesa(den,eint,4,temp)
 call getvalue_mesa(den,eint,5,proint)
 call getvalue_mesa(den,eint,6,peint)
 call getvalue_mesa(den,eint,7,troint)
 call getvalue_mesa(den,eint,8,teint)
 call getvalue_mesa(den,eint,9,entrop)
 call getvalue_mesa(den,eint,10,abad)
 call getvalue_mesa(den,eint,11,gamma1)
 call getvalue_mesa(den,eint,12,gam)

end subroutine get_eos_various_mesa

end module eos_mesa
