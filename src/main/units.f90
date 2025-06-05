!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module units
!
! This module contains information about physical units
!   (irrelevant if simulations are scale free except that
!   the units are printed in the dump file header
!   and log output is scaled in these units)
!
! :References:
!   Price & Monaghan (2004) MNRAS 348, 123 (section 7.1.1 on MHD units)
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
 real(kind=8), public :: unit_ergg, unit_energ, unit_opacity, unit_luminosity
 real(kind=8), public :: unit_angmom

 public :: set_units, set_units_extra, print_units, select_unit
 public :: get_G_code, get_c_code, get_radconst_code, get_kbmh_code
 public :: c_is_unity, G_is_unity, in_geometric_units, in_code_units, in_units
 public :: is_time_unit, is_length_unit, is_mdot_unit, is_density_unit
 public :: in_solarr, in_solarm, in_solarl
 public :: get_unit_multiplier

 integer, parameter :: len_utype = 10

 private

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
       umass = udist**3/(gg*utime**2)
       if (present(mass)) print "(a)",' WARNING: over-riding mass unit with G=1 assumption'
    elseif (present(mass) .and. present(time)) then
       udist = (utime**2*(gg*umass))**(1.d0/3.d0)
       if (present(dist)) print "(a)",' WARNING: over-riding length unit with G=1 assumption'
    elseif (present(time)) then
       umass = udist**3/(gg*utime**2)     ! udist is 1
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

 unit_velocity   = udist/utime
 unit_density    = umass/udist**3
 unit_pressure   = umass/(udist*utime**2)
 unit_ergg       = unit_velocity**2
 unit_energ      = umass*unit_ergg
 unit_opacity    = udist**2/umass
 unit_luminosity = unit_energ/utime
 unit_angmom     = udist*umass*unit_velocity

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

 write(lu,"(/,a)") ' --- code units --- '
 write(lu,"(/,3(a,es10.3,1x),a)") '     Mass: ',umass,    'g       Length: ',udist,  'cm    Time: ',utime,'s'
 write(lu,"(3(a,es10.3,1x),a)") '  Density: ',unit_density, 'g/cm^3  Energy: ',unit_energ,'erg   En/m: ',unit_ergg,'erg/g'
 write(lu,"(3(a,es10.3,1x),a)") ' Velocity: ',unit_velocity,'cm/s    Bfield: ',unit_Bfield,'G  opacity: ',unit_opacity,'cm^2/g'
 write(lu,"(3(a,es10.3,1x))") '        G: ', gg*umass*utime**2/udist**3,'             c: ',c*utime/udist,&
                                '      mu_0: ',cgsmu0*unit_charge**2/(umass*udist)
 write(lu,"(2(a,es10.3,1x),/)") '        a: ',get_radconst_code(),      '         kB/mH: ',get_kbmh_code()

end subroutine print_units

!------------------------------------------------------------------------------------
!+
!  Subroutine to recognise mass, length and time units from a string
!+
!------------------------------------------------------------------------------------
subroutine select_unit(string,unit,ierr,unit_type)
 use physcon
 character(len=*),  intent(in)  :: string
 real(kind=8),      intent(out) :: unit
 integer,           intent(out) :: ierr
 character(len=len_utype), intent(out), optional :: unit_type
 character(len=len(string)) :: unitstr
 character(len=len_utype)   :: utype
 real(kind=8) :: fac

 ierr = 0
 call get_unit_multiplier(string,unitstr,fac,ierr)

 select case(trim(adjustl(unitstr)))
 case('solarr','rsun')
    unit  = solarr
    utype = 'length'
 case('jupiterr','rjup','rjupiter')
    unit  = jupiterr
    utype = 'length'
 case('earthr','rearth')
    unit = earthr
    utype = 'length'
 case('au')
    unit  = au
    utype = 'length'
 case('ly','lightyear')
    unit  = ly
    utype = 'length'
 case('pc','parsec')
    unit  = pc
    utype = 'length'
 case('kpc','kiloparsec')
    unit  = kpc
    utype = 'length'
 case('mpc','megaparsec')
    unit  = mpc
    utype = 'length'
 case('km','kilometres','kilometers')
    unit  = km
    utype = 'length'
 case('cm','centimetres','centimeters')
    unit = 1.d0
    utype = 'length'
 case('solarm','msun')
    unit = solarm
    utype = 'mass'
 case('earthm','mearth')
    unit = earthm
    utype = 'mass'
 case('jupiterm','mjup','mjupiter')
    unit = jupiterm
    utype = 'mass'
 case('ceresm','mceres')
    unit = ceresm
    utype = 'mass'
 case('g','grams')
    unit = gram
    utype = 'mass'
 case('days','day')
    unit = days
    utype = 'time'
 case('Myr')
    unit = 1.d6*years
    utype = 'time'
 case('yr','year','yrs','years')
    unit = years
    utype = 'time'
 case('hr','hour','hrs','hours')
    unit = hours
    utype = 'time'
 case('min','minute','mins','minutes')
    unit = minutes
    utype = 'time'
 case('s','sec','second','seconds')
    unit = seconds
    utype = 'time'
 case('g/s','grams/second','g/second','grams/s','g/sec','grams/sec')
    unit = gram/seconds
    utype = 'mdot'
 case('Ms/yr','M_s/yr','ms/yr','m_s/yr','Msun/yr','M_sun/yr','Msolar/yr',&
      'M_solar/yr','Ms/year','M_s/year','ms/year','m_s/year','Msun/year',&
      'M_sun/year','Msolar/year','M_solar/year','msun/yr')
    unit = solarm/years
    utype = 'mdot'
 case('lsun','solarl','Lsun')
    unit =  solarl
    utype = 'luminosity'
 case('erg/s')
    unit =  1.d0
    utype = 'luminosity'
 case('cm/s')
    unit = cm/seconds
    utype = 'velocity'
 case('m/s')
    unit = 1.d2*cm/seconds
    utype = 'velocity'
 case('km/s')
    unit = km/seconds
    utype = 'velocity'
 case('km/h')
    unit = km/hours
    utype = 'velocity'
 case('au/yr')
    unit = au/years
    utype = 'velocity'
 case('c')
    unit = c
    utype = 'velocity'
 case('g/cm^3','g/cm3','g/cc','g_per_cc','g_per_cm3','g cm^-3')
    unit = gram/cm**3
    utype = 'density'
 case('kg/m^3','kg/m3','kg_per_m3','kg m^-3')
    unit = kg/metre**3
    utype = 'density'
 case default
    if (len_trim(unitstr) > 0) ierr = 1
    unit = 1.d0
    utype = 'none'
 end select

 unit = unit*fac
 if (present(unit_type)) unit_type = utype

end subroutine select_unit

!------------------------------------------------------------------------------------
!+
!  check if string is a unit of time
!+
!------------------------------------------------------------------------------------
logical function is_time_unit(string)
 character(len=*), intent(in) :: string
 character(len=len_utype) :: unit_type
 real(kind=8) :: val
 integer :: ierr

 ierr = 0
 call select_unit(string,val,ierr,unit_type)

 is_time_unit = (trim(unit_type) == 'time')

end function is_time_unit

!------------------------------------------------------------------------------------
!+
!  check if string is a unit of length
!+
!------------------------------------------------------------------------------------
logical function is_length_unit(string)
 character(len=*), intent(in) :: string
 character(len=len_utype) :: unit_type
 real(kind=8) :: val
 integer :: ierr

 ierr = 0
 call select_unit(string,val,ierr,unit_type)

 is_length_unit = (trim(unit_type) == 'length')

end function is_length_unit

!------------------------------------------------------------------------------------
!+
!  check if string is a unit of mdot
!+
!------------------------------------------------------------------------------------
logical function is_mdot_unit(string)
 character(len=*), intent(in) :: string
 character(len=len_utype) :: unit_type
 real(kind=8) :: val
 integer :: ierr

 ierr = 0
 call select_unit(string,val,ierr,unit_type)

 is_mdot_unit = (trim(unit_type) == 'mdot')

end function is_mdot_unit

!------------------------------------------------------------------------------------
!+
!  check if string is a unit of mdot
!+
!------------------------------------------------------------------------------------
logical function is_density_unit(string)
 character(len=*), intent(in) :: string
 character(len=len_utype) :: unit_type
 real(kind=8) :: val
 integer :: ierr

 ierr = 0
 call select_unit(string,val,ierr,unit_type)

 is_density_unit = (trim(unit_type) == 'density')

end function is_density_unit

!------------------------------------------------------------------------------------
!+
!  parse a string like '10.*days' or '10*au' and return the value in code units
!  if there is no recognisable units, the value is returned unscaled
!+
!------------------------------------------------------------------------------------
real function in_code_units(string,ierr,unit_type) result(rval)
 character(len=*), intent(in)  :: string
 integer,          intent(out) :: ierr
 character(len=*), intent(in), optional :: unit_type
 real(kind=8) :: val
 character(len=len_utype) :: utype

 call select_unit(string,val,ierr,unit_type=utype)

 ! return an error if incorrect dimensions (e.g. mass instead of length)
 if (present(unit_type)) then
    if ((trim(utype) /= 'none') .and. trim(utype) /= trim(unit_type)) then
       ierr = 2
       rval = real(val)
       return
    endif
 endif

 if (ierr /= 0) then
    rval = real(val)
    return
 else
    select case(trim(utype))
    case('time')
       rval = real(val/utime)
    case('length')
       rval = real(val/udist)
    case('mass')
       rval = real(val/umass)
    case('mdot')
       rval = real(val/(umass/utime))
    case('luminosity')
       rval = real(val/unit_luminosity)
    case('velocity')
       rval = real(val/unit_velocity)
    case('density')
       rval = real(val/unit_density)
    case default
       rval = real(val)  ! no unit conversion
    end select
 endif

end function in_code_units

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
 character(len=:), allocatable :: tmpstr
 integer :: i,n,stat

 ! default values
 fac = 1.d0
 ierr = 0
 tmpstr = adjustl(string)  ! remove leading spaces

 ! handle empty string
 if (len_trim(tmpstr) == 0) then
    unit_string = ''
    return
 endif

 ! find first space or *
 n = len_trim(tmpstr)
 over_string: do i=1,len_trim(tmpstr)
    if (tmpstr(i:i)==' ' .or. tmpstr(i:i)=='*') then
       n = i - 1
       exit over_string
    endif
 enddo over_string

 ! try to read the first part as a number
 read(tmpstr(1:n),*,iostat=stat) fac
 if (stat == 0) then
    ! skip exactly one separator (* or space) if present
    if (n < len_trim(tmpstr) .and. &
        (tmpstr(n+1:n+1)=='*' .or. tmpstr(n+1:n+1)==' ')) then
       n = n + 1
    endif
    ! extract the unit string, preserving internal spaces
    unit_string = trim(adjustl(tmpstr(n+1:)))
 else
    ! if no number found, use entire string as unit
    unit_string = adjustl(tmpstr)
 endif
 !print*,trim(string),': fac = ',fac,' unitstr = ',trim(unit_string),' error = ',ierr

end subroutine get_unit_multiplier

!---------------------------------------------------------------------------
!+
!  whether or not a character is a number
!+
!---------------------------------------------------------------------------
pure logical function is_digit(ch)
 character(len=1), intent(in) :: ch

 is_digit = (iachar(ch) >= iachar('0') .and. iachar(ch) <= iachar('9')) .or. (ch=='.')

end function is_digit

!---------------------------------------------------------------------------
!+
!  Gravitational constant in code units
!+
!---------------------------------------------------------------------------
real function get_G_code() result(G_code)
 use physcon, only:gg

 G_code = real(gg*umass*utime**2/udist**3)

end function get_G_code

!---------------------------------------------------------------------------
!+
!  speed of light in code units
!+
!---------------------------------------------------------------------------
real function get_c_code() result(c_code)
 use physcon, only:c

 c_code = real(c*utime/udist)

end function get_c_code

!---------------------------------------------------------------------------
!+
!  radiation constant
!+
!---------------------------------------------------------------------------
real function get_radconst_code() result(radconst_code)
 use physcon, only:radconst

 radconst_code = real(radconst/unit_energ*udist**3)

end function get_radconst_code

!---------------------------------------------------------------------------
!+
!  speed of light in code units
!+
!---------------------------------------------------------------------------
real function get_kbmh_code() result(kbmh_code)
 use physcon, only:kb_on_mh

 kbmh_code = real(kb_on_mh/unit_velocity**2)

end function get_kbmh_code

!---------------------------------------------------------------------------
!+
!  whether or not the Gravitational constant is unity in code units
!+
!---------------------------------------------------------------------------
logical function G_is_unity()

 G_is_unity = abs(get_G_code() - 1.) < 1.e-12

end function G_is_unity

!---------------------------------------------------------------------------
!+
!  whether or not the speed of light is unity in code units
!+
!---------------------------------------------------------------------------
logical function c_is_unity()

 c_is_unity = abs(get_c_code() - 1.) < 1.e-12

end function c_is_unity
!---------------------------------------------------------------------------
!+
!  logical to check we are in geometric units (i.e. c = G = 1)
!+
!---------------------------------------------------------------------------
logical function in_geometric_units()

 in_geometric_units = c_is_unity() .and. G_is_unity()

end function in_geometric_units

!---------------------------------------------------------------------------
!+
!  function to convert a mass value from code units to solar masses
!+
!---------------------------------------------------------------------------
real(kind=8) function in_solarm(val) result(rval)
 use physcon, only:solarm
 real, intent(in) :: val

 rval = val*(umass/solarm)

end function in_solarm
!---------------------------------------------------------------------------
!+
!  function to convert a distance value from code units to solar radii
!+
!---------------------------------------------------------------------------
real(kind=8) function in_solarr(val) result(rval)
 use physcon, only:solarr
 real, intent(in) :: val

 rval = val*(udist/solarr)

end function in_solarr
!---------------------------------------------------------------------------
!+
!  function to convert a luminosity value from code units to solar luminosity
!+
!---------------------------------------------------------------------------
real(kind=8) function in_solarl(val) result(rval)
 use physcon, only:solarl
 real, intent(in) :: val

 rval = val*(unit_luminosity/solarl)

end function in_solarl

!------------------------------------------------------------------------------------
!+
!  print a value in physical units, e.g. give code value of mass and
!  call this routine print*,in_units(mass,'solarm')
!+
!------------------------------------------------------------------------------------
real(kind=8) function in_units(val,unitstring) result(rval)
 real,             intent(in)  :: val
 character(len=*), intent(in)  :: unitstring
 character(len=len_utype) :: utype
 integer :: ierr
 real(kind=8) :: fac

 call select_unit(unitstring,fac,ierr,unit_type=utype)  ! handle errors silently by ignoring ierr

 select case(trim(utype))
 case('time')
    rval = val*(utime/fac)
 case('length')
    rval = val*(udist/fac)
 case('mass')
    rval = val*(umass/fac)
 case('mdot')
    rval = val*((umass/utime)/fac)
 case('luminosity')
    rval = val*(unit_luminosity/fac)
 case('velocity')
    rval = val*(unit_velocity/fac)
 case('density')
    rval = val*(unit_density/fac)
 case default
    rval = real(val)  ! no unit conversion
 end select

end function in_units

end module units
