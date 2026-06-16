!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_zerotemp
!
! Implements zero temperature equation of state, e.g. for white dwarfs. Meant to 
!
! :References: Kippenhahn & Weigert, Stellar Structure and Evolution, section 15.2
!
! :Dependencies: physcon
!
 use units, only:unit_density,unit_velocity
 use physcon,  only:pi,atomic_mass_unit, mass_electron_cgs, planckh, c
 implicit none
 real, parameter :: mu_e = 2 ! mean molecular weight per free electron


 public :: f_chandra, get_zerotemp_pressure,get_zerotemp_u, get_zerotemp_spsoundi,&
           get_zerotemp_rhofrompres

 private

contains

!----------------------------------------------------------------
!+
!  Chandrasekhar's EOS function
!+
!----------------------------------------------------------------

real function f_chandra(x) result(fx)
 real, intent(in) :: x

 fx = x*(2*x**2 - 3)*sqrt(1+x**2) + 3*log(x + sqrt(1+x**2))
end function f_chandra



! ----------------------------------------------------------------
!+
!  Calculates the zero temperature pressure for a fully degenerate electron gas, i.e. a
!  white dwarf. See Kippenhahn & Weigert, Stellar Structure and Evolution, section 15.2
!  Note that this is only the electron degeneracy pressure, so does not include the ion
!+ input/output is cgs
! ----------------------------------------------------------------

subroutine get_zerotemp_pressure(rhoi,presi)
 real, intent(in)  :: rhoi
 real, intent(out) :: presi
 real :: ne, x

 ! This is the zero temperature pressure for a fully degenerate electron gas, i.e. a white dwarf. See Kippenhahn & Weigert, Stellar Structure and Evolution, section 15.2
 ! Note that this is only the electron degeneracy pressure, so does not include the ion contribution to the pressure.
 ! This is a good approximation for white dwarfs where the electrons are highly degenerate but the ions are not.
 ! Note also that this assumes a fully ionised gas, so mu is the mean molecular weight per free electron


    ne = rhoi / (mu_e * atomic_mass_unit)

    x = ( (3.0 * ne * planckh**3) / &
        (8.0 * pi * mass_electron_cgs**3 * c**3) )**(1.0/3.0)

    presi = (pi * mass_electron_cgs**4 * c**5 / (3.0 * planckh**3)) * f_chandra(x)

end subroutine get_zerotemp_pressure


!-----------------------------------------------------------------------
!+!  Inputs and outputs in cgs units
!+!  get internal energy from density for zero temperature EOS
!-----------------------------------------------------------------------
subroutine get_zerotemp_u(rhoi,u)
 real,    intent(in)    :: rhoi
 real,    intent(out)   :: u
 real :: ne, x, gx

    ne = rhoi / (mu_e * atomic_mass_unit)

    x = ( (3.0 * ne * planckh**3) / &
        (8.0 * pi * mass_electron_cgs**3 * c**3) )**(1.0/3.0)

    gx = (8*x**3)*(sqrt(x**2 + 1)-1)-f_chandra(x)
    u = (pi * mass_electron_cgs**4 * c**5 / (3.0 * planckh**3)) * f_chandra(x)

end subroutine get_zerotemp_u


!----------------------------------------------------------------
!+
!  Calculates sound speed from density (derivative of f(x), analytically found), input output in cgs
!+
!----------------------------------------------------------------

subroutine get_zerotemp_spsoundi(rhoi,spsoundi)
 real, intent(in)  :: rhoi
 real, intent(out) :: spsoundi
 real :: ne, x

 ne = rhoi / (mu_e * atomic_mass_unit)

 x = ( (3.0 * ne * planckh**3) / &
        (8.0 * pi * mass_electron_cgs**3 * c**3) )**(1.0/3.0)

 constants_cs = (mass_electron_cgs * c**2 / (24.0* atomic_mass_unit*mu_e))
 spsoundi = sqrt(constants_cs * ((6.0 * x**2 - 3.0)*sqrt(1.0 + x**2) + (2.0 * x**4 - 3.0*x**2 +3.0)/sqrt(1.0+x**2))/x**2)
end subroutine get_zerotemp_spsoundi

!----------------------------------------------------------------
!+
!  Calculates density from pressure (involves the bisection method), input output in cgs
!+
!----------------------------------------------------------------
subroutine get_zerotemp_rhofrompres(presi,densi,ierr)

   real, intent(in)  :: presi
   real, intent(out) :: densi
   integer, intent(out) :: ierr

   real :: ne, x, fx
   real :: xlo, xhi, xmid
   integer :: iter

   integer, parameter :: iter_max = 1000
   real,    parameter :: tolerance = 1.e-12

   fx = presi * (3.0*planckh**3) / &
        (pi*mass_electron_cgs**4*c**5)

   ! Initial bracket
   xlo = 0.0
   xhi = 1.0

   do while (f_chandra(xhi) < fx)
      xhi = 2.0*xhi
   end do

   ierr = 0

   do iter = 1, iter_max

      xmid = 0.5*(xlo + xhi)

      if (f_chandra(xmid) > fx) then
         xhi = xmid
      else
         xlo = xmid
      end if

      if (abs(xhi-xlo)/(xmid+1.e-300) < tolerance) exit

   end do

   if (iter == iter_max) ierr = 1

   x = 0.5*(xlo + xhi)

   ne = (8.0*pi*mass_electron_cgs**3*c**3 / &
        (3.0*planckh**3)) * x**3

   densi = ne * mu_e * atomic_mass_unit

end subroutine get_zerotemp_rhofrompres


end module eos_zerotemp