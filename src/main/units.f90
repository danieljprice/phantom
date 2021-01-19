!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module units
!
! This module contains information about physical units
!   (irrelevant if simulations are scale free except that
!   the units are printed in the dump file header
!   and log output is scaled in these units)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 implicit none
!
! Default units (note that values below are overwritten by the call to set_units
!                but this prevents errors in case set_units was never called)
!
 real(kind=8), public :: udist = 1.d0, umass = 1.d0, utime = 1.d0
 real(kind=8), public :: unit_velocity, unit_Bfield, unit_charge
 real(kind=8), public :: unit_pressure, unit_density
 real(kind=8), public :: unit_ergg, unit_energ, unit_opacity

 public :: set_units, set_units_extra, print_units
 public :: get_G_code, get_c_code, get_steboltz_code
 public :: c_is_unity, G_is_unity, in_geometric_units

contains

!------------------------------------------------------------------------------------
!+
!  Routine to set the basic units of mass,
!  length and time for the code. This should be
!  called from the setup_blah.f90 routine for your calculation.
!
!  In Phantom the units are written and read from the dump file header
!
!  Input:
!   dist : (optional) distance unit
!   mass : (optional) mass unit
!   time : (optional) time unit
!   G    : (optional) can be used to specify time unit such that G=1
!   c    : (optional) can be used to specify length unit such that c=1
!
!  Output:
!   None, sets the values of the public module variables udist,umass,utime
!   as well as derivative units for other quantities (uBfield,upres,udens,uenerg)
!+
!------------------------------------------------------------------------------------
subroutine set_units(dist,mass,time,G,c)
 use physcon, only:gg,clight=>c
 real(kind=8), intent(in), optional :: dist,mass,time,G,c

 if (present(dist)) then
    udist = dist
 else
    udist = 1.d0
 endif
 if (present(mass)) then
    umass = mass
 else
    umass = 1.d0
 endif
 if (present(time)) then
    utime = time
 else
    utime = 1.d0
 endif

 if (present(c)) then
    if (present(mass)) then
       udist = gg*umass/clight**2
       utime = udist/clight
       if (present(dist)) print "(a)",' WARNING: over-riding length unit with c=1 assumption'
    elseif (present(dist)) then
       utime = udist/clight
       umass = clight*clight*udist/gg
       if (present(time)) print "(a)",' WARNING: over-riding time unit with c=1 assumption'
    elseif (present(time)) then
       udist = utime*clight
       umass = clight*clight*udist/gg
    else
       udist = gg*umass/clight**2   ! umass is 1
       utime = udist/clight
    endif
 elseif (present(G)) then
    if (present(mass) .and. present(dist)) then
       utime = sqrt(udist**3/(gg*umass))
       if (present(time)) print "(a)",' WARNING: over-riding time unit with G=1 assumption'
    elseif (present(dist) .and. present(time)) then
       umass = udist**2/(gg*utime**2)
       if (present(mass)) print "(a)",' WARNING: over-riding mass unit with G=1 assumption'
    elseif (present(mass) .and. present(time)) then
       udist = (utime**2*(gg*umass))**(1.d0/3.d0)
       if (present(dist)) print "(a)",' WARNING: over-riding length unit with G=1 assumption'
    elseif (present(time)) then
       umass = udist**2/(gg*utime**2)     ! udist is 1
    else
       utime = sqrt(udist**3/(gg*umass))  ! udist and umass are 1
    endif
 endif

 call set_units_extra()

end subroutine set_units

!------------------------------------------------------------------------------------
!+
!  Subroutine to set the values of useful unit conversions
!  within Phantom - these derive from the basic ones set
!  in set_units.
!
!  This routine should be called after the units have been
!  read from the dump file header. Thereafter these unit
!  conversions are available throughout the code.
!+
!------------------------------------------------------------------------------------
subroutine set_units_extra()
 use physcon, only:cgsmu0
 !
 ! The way magnetic field units are set such that mu_0 = 1
 ! is described in Ref:pricemonaghan04, section 7.1.1
 !
 unit_charge   = sqrt(umass*udist/cgsmu0)
 unit_Bfield   = umass/(utime*unit_charge)

 unit_velocity = udist/utime
 unit_density  = umass/udist**3
 unit_pressure = umass/(udist*utime**2)
 unit_ergg     = unit_velocity**2
 unit_energ    = umass*unit_ergg
 unit_opacity  = udist**2/umass

end subroutine set_units_extra

!------------------------------------------------------------------------------------
!+
!  Subroutine to pretty-print unit conversion information to the log file
!+
!------------------------------------------------------------------------------------
subroutine print_units(unit)
 use physcon, only:gg,c,cgsmu0
 integer, intent(in), optional :: unit
 integer :: lu

 if (present(unit)) then
    lu = unit
 else
    lu = 6
 endif

 write(lu,"(a)") ' --- code units --- '
 write(lu,"(/,3(a,es10.3,1x),a)") '     Mass: ',umass,    'g       Length: ',udist,  'cm    Time: ',utime,'s'
 write(lu,"(3(a,es10.3,1x),a)") '  Density: ',unit_density, 'g/cm^3  Energy: ',unit_energ,'erg   En/m: ',unit_ergg,'erg/g'
 write(lu,"(2(a,es10.3,1x),a)") ' Velocity: ',unit_velocity,'cm/s    Bfield: ',unit_Bfield,'G'
 write(lu,"(3(a,es10.3,1x),/)")   '        G: ', gg*umass*utime**2/udist**3,'             c: ',c*utime/udist,&
                                 '      mu_0: ',cgsmu0*unit_charge**2/(umass*udist)

end subroutine print_units

!------------------------------------------------------------------------------------
!+
!  Subroutine to recognise mass and length units from a string
!+
!------------------------------------------------------------------------------------
subroutine select_unit(string,unit,ierr)
 use physcon
 character(len=*), intent(in)  :: string
 real(kind=8),     intent(out) :: unit
 integer,          intent(out) :: ierr
 character(len=len(string)) :: unitstr
 real(kind=8) :: fac

 ierr = 0
 call get_unit_multiplier(string,unitstr,fac,ierr)

 select case(trim(unitstr))
 case('solarr','rsun')
    unit = solarr
 case('au')
    unit = au
 case('ly','lightyear')
    unit = ly
 case('pc','parsec')
    unit = pc
 case('kpc','kiloparsec')
    unit = kpc
 case('mpc','megaparsec')
    unit = mpc
 case('km','kilometres','kilometers')
    unit = km
 case('cm','centimetres','centimeters')
    unit = 1.d0
 case('solarm','msun')
    unit = solarm
 case('earthm','mearth')
    unit = earthm
 case('jupiterm','mjup','mjupiter')
    unit = jupiterm
 case('g','grams')
    unit = 1.d0
 case default
    ierr = 1
    unit = 1.d0
 end select

 unit = unit*fac

end subroutine select_unit

!------------------------------------------------------------------------------------
!+
!  Utility routine to extract the number in front of the unit info (e.g. 100 au)
!+
!------------------------------------------------------------------------------------
subroutine get_unit_multiplier(string,unit_string,fac,ierr)
 character(len=*), intent(in)  :: string
 character(len=*), intent(out) :: unit_string
 real(kind=8),     intent(out) :: fac
 integer,          intent(out) :: ierr
 integer :: i,i1,i2

 fac = 1.d0
 i1 = 0
 i2 = 0
 ierr = 0
 !
 ! loop backwards through the string, looking
 ! for the first character that is a number
 !
 over_string: do i=len_trim(string),1,-1
    if (string(i:i)=='*') then
       i1 = i-1
       i2 = i
       exit over_string
    elseif (is_digit(string(i:i))) then
       i1 = i
       i2 = i
       exit over_string
    endif
 enddo over_string
 if (i1 > 0) then
    read(string(1:i1),*,iostat=ierr) fac
    unit_string = adjustl(string(i2+1:))
 else
    unit_string = string
 endif
 !print*,'fac = ',fac,' unitstr = ',trim(unit_string)

end subroutine get_unit_multiplier

!---------------------------------------------------------------------------
!+
!  whether or not a character is a number
!+
!---------------------------------------------------------------------------
pure logical function is_digit(ch)
 character(len=1), intent(in) :: ch

 is_digit = (iachar(ch) >= iachar('0') .and. iachar(ch) <= iachar('9'))

end function is_digit

!---------------------------------------------------------------------------
!+
!  Gravitational constant in code units
!+
!---------------------------------------------------------------------------
real(kind=8) function get_G_code() result(G_code)
 use physcon, only:gg

 G_code = gg*umass*utime**2/udist**3

end function get_G_code

!---------------------------------------------------------------------------
!+
!  speed of light in code units
!+
!---------------------------------------------------------------------------
real(kind=8) function get_c_code() result(c_code)
 use physcon, only:c

 c_code = c*utime/udist

end function get_c_code

!---------------------------------------------------------------------------
!+
!  Stefan-Boltzmann constant in code units
!+
!---------------------------------------------------------------------------
real(kind=8) function get_steboltz_code() result(steboltz_code)
 use physcon, only:steboltz

 steboltz_code = steboltz/(unit_energ/(udist**2*utime))

end function get_steboltz_code

!---------------------------------------------------------------------------
!+
!  whether or not the Gravitational constant is unity in code units
!+
!---------------------------------------------------------------------------
logical function G_is_unity()

 G_is_unity = abs(get_G_code() - 1.d0) < 1.d-12

end function G_is_unity

!---------------------------------------------------------------------------
!+
!  whether or not the speed of light is unity in code units
!+
!---------------------------------------------------------------------------
logical function c_is_unity()

 c_is_unity = abs(get_c_code() - 1.d0) < 1.d-12

end function c_is_unity
!---------------------------------------------------------------------------
!+
!  logical to check we are in geometric units (i.e. c = G = 1)
!+
!---------------------------------------------------------------------------
logical function in_geometric_units()

 in_geometric_units = c_is_unity() .and. G_is_unity()

end function in_geometric_units

end module units
