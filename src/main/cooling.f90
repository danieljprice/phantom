!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: cooling
!
!  DESCRIPTION:
!  Interface routines and options for
!  various cooling prescriptions
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    C_cool    -- factor controlling cooling timestep
!    beta_cool -- beta factor in Gammie (2001) cooling
!    icooling  -- cooling function (0=off, 1=Gammie cooling 2=SD93)
!
!  DEPENDENCIES: h2cooling, infile_utils, io, options, part, timestep
!+
!--------------------------------------------------------------------------
module cooling
 use options,  only:icooling
 use timestep, only:C_cool
 use part,     only:h2chemistry,maxvxyzu
 implicit none
 real :: beta_cool = 3.

 public :: energ_cooling
 public :: write_options_cooling, read_options_cooling
 public :: cooling_rate_sd93

 private

contains
!----------------------------------------------------
!+
!  Implementation of various cooling prescriptions
!
!  INPUT: icool -- cooling option
!         ui    -- thermal energy per unit mass
!         xi, yi, zi -- current position
!  OUTPUT:
!         dudti -- cooling rate
!
!  The default cooling is from Gammie (2001)
!  i.e. du/dt = -u/tcool, where tcool = beta/Omega
!+
!----------------------------------------------------
!subroutine energ_cooling(icool,ui,dudti,xi,yi,zi)
subroutine energ_cooling(icool,ui,dudti,xi,yi,zi,rhoi,vxyzui)
 use units, only:utime,unit_ergg,umass,udist
 use options, only:ieos
 use eos, only:get_temperature
 use dim, only:maxvxyzu
 use physcon, only:atomic_mass_unit
 integer, intent(in)    :: icool
 !real,    intent(in)    :: ui,xi,yi,zi
 real,    intent(in)    :: ui,xi,yi,zi,rhoi
 real,         intent(in) :: vxyzui(maxvxyzu)
 real,    intent(inout) :: dudti
 real :: r2,Omegai,tcool1,temp,crate,fac

 select case(icool)
 case(2)
    !
    ! SD93 cooling
    !
    !temp = 1. !get_temperature(ieos,(/xi,yi,zi/),rhoi,ui)
    temp = get_temperature(ieos,(/xi,yi,zi/),rhoi,vxyzui)
    !
    ! convert cooling rate in cgs units to code units
    !
    !crate = cooling_rate_sd93(temp) !*unit_ergg/utime
    !crate = cooling_rate_sd93(temp)*unit_ergg/utime
    if(temp > 1.e4) then
       fac = unit_ergg/utime/umass*udist**3
       crate = cooling_rate_sd93(temp)/atomic_mass_unit**2/fac  
!WRITE(*,*) 'cooling debug: ',vxyzui(4),dudti,crate,dudti+crate
       dudti = dudti + crate*rhoi
    endif

 case default
    !
    ! Gammie (2001) cooling
    !
    r2 = xi*xi + yi*yi + zi*zi
    Omegai = r2**(-0.75)
    tcool1 = Omegai/beta_cool
    dudti  = dudti - ui*tcool1
 end select

end subroutine energ_cooling

!---------------------------------------------------------
!+
!  Cooling function from SD93 as function of temperature
!  Returns (heating rate-cooling rate)/n_H2 in cgs units
!+
!---------------------------------------------------------
pure real function cooling_rate_sd93(T)
 real, intent(in) :: T
 real, parameter :: Tbr1 = 3.e4, Tbr2 = 4.e7
 real, parameter :: fac = -6.4e-23
 real, parameter :: fac1 = 1.e-7   ! gives T/1.e7

 if (T > Tbr2) then
    cooling_rate_sd93 = fac*(Tbr2*fac1)**(-0.7)*(T/Tbr2)**0.5
 elseif (T > Tbr1) then
    cooling_rate_sd93 = fac*(T*fac1)**(-0.7)
 else
    cooling_rate_sd93 = fac*(T*fac1)**(-0.7)*(T/Tbr1)**2
 endif

end function cooling_rate_sd93

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling(iunit)
 use infile_utils, only:write_inopt
 use h2cooling,    only:write_options_h2cooling
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling cooling'
 if (h2chemistry) then
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=on)',iunit)
 else
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=Gammie cooling 2=SD93)',iunit)
 endif
 if (icooling > 0) then
    call write_inopt(C_cool,'C_cool','factor controlling cooling timestep',iunit)
 endif
 if (h2chemistry) then
    call write_options_h2cooling(iunit)
 elseif (icooling == 1) then
    call write_inopt(beta_cool,'beta_cool','beta factor in Gammie (2001) cooling',iunit)
 endif

end subroutine write_options_cooling

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling(name,valstring,imatch,igotall,ierr)
 use h2cooling, only:read_options_h2cooling
 use io,        only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 logical :: igotallh2

 imatch  = .true.
 igotall = .false.  ! cooling options are compulsory
 igotallh2 = .true.
 if (maxvxyzu < 4) igotall = .true. ! options unnecessary if isothermal

 select case(trim(name))
 case('icooling')
    read(valstring,*,iostat=ierr) icooling
    ngot = ngot + 1
 case('C_cool')
    read(valstring,*,iostat=ierr) C_cool
    ngot = ngot + 1
 case('beta_cool')
    read(valstring,*,iostat=ierr) beta_cool
    ngot = ngot + 1
    if (beta_cool < 1.) call fatal('read_options','beta_cool must be >= 1')
 case default
    imatch = .false.
    if (h2chemistry) call read_options_h2cooling(name,valstring,imatch,igotallh2,ierr)
 end select

 if (igotallh2 .and. ngot >= 1) igotall = .true.

 if (icooling > 1 .and. ngot >= 2) igotall = .true.

 if (.not.h2chemistry .and. (icooling == 1 .and. ngot < 3)) igotall= .false.

end subroutine read_options_cooling

end module cooling
