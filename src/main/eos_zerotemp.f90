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
 use physcon,  only:pi,atomic_mass_unit,mass_electron_cgs,planckh,c
 implicit none
 real :: mu_e ! mean molecular weight per free electron


 public :: eos_zerotemp_init,read_options_eos_zerotemp,write_options_eos_zerotemp,&
           f_chandra,get_zerotemp_pressure,get_zerotemp_u,get_zerotemp_spsoundi,&
           get_zerotemp_rhofrompres,eos_zerotemp_calc_mu_e,eos_zerotemp_eosinfo

 private

  ! these set the mixture of species
 ! following elements can be set by user at runtime
 real :: xh  = 0.0
 real :: xhe = 0.0
 real :: xc  = 0.5
 real :: xo  = 0.5
 real :: xne = 0.0
 real :: xmg = 0.0 
 
 integer, parameter :: speciesmax = 6
 character(len=10) :: speciesname(speciesmax)
 real :: xmass(speciesmax) ! mass fraction of species
 real :: Aion(speciesmax)  ! number of nucleons
 real :: Zion(speciesmax)  ! number of protons
contains


!----------------------------------------------------------------
!+
!  initialise zero-temperature EOS composition
!+
!----------------------------------------------------------------
subroutine eos_zerotemp_init(ierr)

 use io, only: warning, fatal
 integer, intent(out) :: ierr

 ierr = 0

 !----------------------------
 ! species definitions
 !----------------------------
 speciesname(1) = "hydrogen"
 speciesname(2) = "helium"
 speciesname(3) = "carbon"
 speciesname(4) = "oxygen"
 speciesname(5) = "neon"
 speciesname(6) = "magnesium"

 Aion(1) = 1.0   ; Zion(1) = 1.0
 Aion(2) = 4.0   ; Zion(2) = 2.0
 Aion(3) = 12.0  ; Zion(3) = 6.0
 Aion(4) = 16.0  ; Zion(4) = 8.0
 Aion(5) = 20.0  ; Zion(5) = 10.0
 Aion(6) = 24.0  ; Zion(6) = 12.0

 !----------------------------
 ! pack mass fractions
 !----------------------------
 xmass(:) = 0.0
 xmass(1) = xh
 xmass(2) = xhe
 xmass(3) = xc
 xmass(4) = xo
 xmass(5) = xne
 xmass(6) = xmg

 if (abs(sum(xmass(:)) - 1.0) > 1e-6) then
    call warning('eos_zerotemp_init','mass fractions do not sum to 1')
    ierr = 1
    return
 endif

 call eos_zerotemp_calc_mu_e()
 write(*,'(a,f10.7)') 'mean mu_e = ', mu_e
end subroutine eos_zerotemp_init

!----------------------------------------------------------------
!+
!  read options from input file
!+
!----------------------------------------------------------------
subroutine read_options_eos_zerotemp(db,nerr)

 use infile_utils, only: inopts, read_inopt
 type(inopts), intent(inout) :: db(:)
 integer, intent(inout) :: nerr

 call read_inopt(xh ,'xh',db,nerr)
 call read_inopt(xhe,'xhe',db,nerr)
 call read_inopt(xc ,'xc',db,nerr)
 call read_inopt(xo ,'xo',db,nerr)
 call read_inopt(xne,'xne',db,nerr)
 call read_inopt(xmg,'xmg',db,nerr)

end subroutine read_options_eos_zerotemp

!----------------------------------------------------------------
!+
!  write options to input file
!+
!----------------------------------------------------------------
subroutine write_options_eos_zerotemp(iunit)

 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(xh ,'xh' ,'Hydrogen mass fraction',iunit)
 call write_inopt(xhe,'xhe','Helium mass fraction',iunit)
 call write_inopt(xc ,'xc' ,'Carbon mass fraction',iunit)
 call write_inopt(xo ,'xo' ,'Oxygen mass fraction',iunit)
 call write_inopt(xne,'xne','Neon mass fraction',iunit)
 call write_inopt(xmg,'xmg','Magnesium mass fraction',iunit)

end subroutine write_options_eos_zerotemp
!----------------------------------------------------------------
!+
!  print eos information
!+
!----------------------------------------------------------------
subroutine eos_zerotemp_eosinfo(iprint)
 integer, intent(in) :: iprint
 integer :: i

 write(iprint,"(/,a)") ' Zero temperature equation of state'
 write(iprint,"(a)") '   mass fractions of each species:'
 do i=1,speciesmax
    if (xmass(i) > 0.0) then
       write(iprint,"(a,a,a,f5.3)") '     ', speciesname(i), ': ', xmass(i)
    endif
 enddo
 write(iprint,'(a,f10.7)') '   mu_e = ', mu_e
end subroutine eos_zerotemp_eosinfo

!----------------------------------------------------------------
!+
!  compute mean molecular weight per electron
!+
!----------------------------------------------------------------
subroutine eos_zerotemp_calc_mu_e()

 mu_e = 1.0 / sum(xmass(:) * Zion(:) / Aion(:))

end subroutine eos_zerotemp_calc_mu_e

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
    u = (pi * mass_electron_cgs**4 * c**5 / (3.0 * planckh**3)) * gx
    u = u/rhoi !previous equation is erg/cm^3 so here I convert it
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

 if (rhoi <= 0.0) then
    spsoundi = 0.0
    return
 endif

 ne = rhoi / (mu_e * atomic_mass_unit)

 x = ( (3.0 * ne * planckh**3) / &
        (8.0 * pi * mass_electron_cgs**3 * c**3) )**(1.0/3.0)

 spsoundi = sqrt((mass_electron_cgs * c**2 * x**2) / (3.0 * atomic_mass_unit * mu_e * sqrt(1.0 + x**2)))
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