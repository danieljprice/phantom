!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_barotropic
!
! Implements barotropic equation of state, e.g. for star formation
!
! :References: Bate, Bonnell & Bromm (2003)
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - drhocrit : *transition size between rhocrit0 & 1 (fraction of rhocrit0; barotropic eos)*
!   - gamma1   : *adiabatic index 1 (barotropic eos)*
!   - gamma2   : *adiabatic index 2 (barotropic eos)*
!   - gamma3   : *adiabatic index 3 (barotropic eos)*
!   - gamma4   : *adiabatic index 4 (barotropic eos)*
!   - rhocrit0 : *critical density 0 in g/cm^3 (barotropic eos)*
!   - rhocrit1 : *critical density 1 in g/cm^3 (barotropic eos)*
!   - rhocrit2 : *critical density 2 in g/cm^3 (barotropic eos)*
!   - rhocrit3 : *critical density 3 in g/cm^3 (barotropic eos)*
!   - rhocrit4 : *critical density 4 in g/cm^3 (barotropic eos)*
!
! :Dependencies: infile_utils, io, units
!
 use units, only:unit_density,unit_velocity
 implicit none
 !--Default initial parameters for Barotropic Eos
 real,    public :: drhocrit0   = 0.50
 real,    public :: rhocrit0cgs = 1.e-18
 real,    public :: rhocrit1cgs = 1.e-14
 real,    public :: rhocrit2cgs = 1.e-10
 real,    public :: rhocrit3cgs = 1.e-3
 real,    public :: rhocrit4cgs = 1.e-1
 real,    public :: gamma1      = 1.4
 real,    public :: gamma2      = 1.1
 real,    public :: gamma3      = 5./3.
 real,    public :: gamma4      = 1.1

 real :: rhocritT,rhocrit0,rhocrit1,rhocrit2,rhocrit3,rhocrit4
 real :: fac2,fac3,fac4,log10polyk2,log10rhocritT,rhocritT0slope

 public :: init_eos_barotropic
 public :: get_eos_barotropic,eos_info_barotropic
 public :: write_options_eos_barotropic,read_options_eos_barotropic
 public :: gamma_barotropic

 private

contains

!-----------------------------------------------------------------------
!+
!  Initialise the equation of state
!+
!-----------------------------------------------------------------------
subroutine init_eos_barotropic(polyk,polyk2,ierr)
 use io, only:warning
 real, intent(in) :: polyk,polyk2
 integer, intent(out) :: ierr
 !
 !--calculate initial variables for the barotropic equation of state
 !
 if (unit_density <= 0.) then
    ierr = 3
    return
 endif

 ! Convert to code units, and calculate constants
 rhocrit0 = real(rhocrit0cgs/unit_density)
 rhocrit1 = real(rhocrit1cgs/unit_density)
 rhocrit2 = real(rhocrit2cgs/unit_density)
 rhocrit3 = real(rhocrit3cgs/unit_density)
 rhocrit4 = real(rhocrit4cgs/unit_density)
 fac2     = polyk*(rhocrit2/rhocrit1)**(gamma1-1.)
 fac3     =  fac2*(rhocrit3/rhocrit2)**(gamma2-1.)
 fac4     =  fac3*(rhocrit4/rhocrit3)**(gamma3-1.)

 ! verify that the rhocrit's are in the correct order
 call verify_less_than(ierr,rhocrit0,rhocrit1)
 call verify_less_than(ierr,rhocrit1,rhocrit2)
 call verify_less_than(ierr,rhocrit2,rhocrit3)
 call verify_less_than(ierr,rhocrit3,rhocrit4)
 ! Calculate values for the first transition region (no transition if drhocrit0=0)
 if (polyk < tiny(polyk) .or. polyk2 < tiny(polyk2)) drhocrit0 = 0.0

 if (drhocrit0 > 0.0) then
    rhocritT       = rhocrit0*(1.0-drhocrit0)
    log10polyk2    = log10(polyk2)
    log10rhocritT  = log10(rhocritT)
    rhocritT0slope = (log10(polyk)-log10(polyk2)) /(log10(rhocritT)-log10(rhocrit0))
 else
    rhocritT       = rhocrit0  ! moving the transition boundary to rhocrit0
    rhocrit0       = 0.0       ! removing the valid threshhold to enter the transition region
    log10polyk2    = 0.0
    log10rhocritT  = 0.0
    rhocritT0slope = 0.0
 endif

 ! Reset rhocrit0 if a warm medium is not defined
 if (rhocrit0cgs > 0.0 .and. polyk2 < tiny(polyk2)) then
    call warning('init_eos','warm medium defined by critical density rho0 but not polyk2.  Resetting rho0 = 0.')
    drhocrit0   = 0.0
    rhocritT    = 0.0
    rhocrit0    = 0.0
    rhocrit0cgs = 0.0
 endif

end subroutine init_eos_barotropic

!-----------------------------------------------------------------------
!+
!  Main eos routine: calculates pressure at a given density
!+
!-----------------------------------------------------------------------
subroutine get_eos_barotropic(rhoi,polyk,polyk2,ponrhoi,spsoundi,gammai)
 real, intent(in)  :: rhoi,polyk,polyk2
 real, intent(out) :: ponrhoi,spsoundi,gammai

 ! variables calculated in the eos initialisation routine:
 !    fac2 = polyk*(rhocrit2/rhocrit1)**(gamma1-1.)
 !    fac3 =  fac2*(rhocrit3/rhocrit2)**(gamma2-1.)
 !    fac4 =  fac3*(rhocrit4/rhocrit3)**(gamma2-1.)
 !    rhocritT0slope = (log10(polyk)-log10(polyk2)) &
 !                   /(log10(rhocritT)-log10(rhocrit0)))
 !
 if (rhoi < rhocritT) then
    gammai  = 1.0
    ponrhoi = polyk2
 elseif (rhoi < rhocrit0) then
    gammai  = 1.0
    ponrhoi = 10**(log10polyk2 + rhocritT0slope*(log10rhocritT-log10(rhoi))  )
 elseif (rhoi < rhocrit1) then
    gammai  = 1.0
    ponrhoi = polyk
 elseif (rhoi < rhocrit2) then
    gammai  = gamma1
    ponrhoi = polyk*(rhoi/rhocrit1)**(gamma1-1.)
 elseif (rhoi < rhocrit3) then
    gammai  = gamma2
    ponrhoi = fac2*(rhoi/rhocrit2)**(gamma2-1.)
 elseif (rhoi < rhocrit4) then
    gammai  = gamma3
    ponrhoi = fac3*(rhoi/rhocrit3)**(gamma3-1.)
 else
    gammai  = gamma4
    ponrhoi = fac4*(rhoi/rhocrit4)**(gamma4-1.)
 endif
 spsoundi = sqrt(gammai*ponrhoi)

end subroutine get_eos_barotropic

!-----------------------------------------------------------------------
!+
!  Get gamma for thermal energy calculations when using the
!  piecewise polytrope
!+
!-----------------------------------------------------------------------
real function gamma_barotropic(rhoi) result(gammai)
 real, intent(in) :: rhoi

 if (rhoi < rhocritT) then
    gammai  = 1.0
 elseif (rhoi < rhocrit0) then
    gammai  = 1.0
 elseif (rhoi < rhocrit1) then
    gammai  = 1.0
 elseif (rhoi < rhocrit2) then
    gammai  = gamma1
 elseif (rhoi < rhocrit3) then
    gammai  = gamma2
 elseif (rhoi < rhocrit4) then
    gammai  = gamma3
 else
    gammai  = gamma4
 endif

end function gamma_barotropic

!-----------------------------------------------------------------------
!+
!  print information about the equation of state parameters
!+
!-----------------------------------------------------------------------
subroutine eos_info_barotropic(polyk,polyk2,iprint)
 real,    intent(in) :: polyk,polyk2
 integer, intent(in) :: iprint
 character(len=*), parameter :: cu   = ' code units = '
 character(len=*), parameter :: baro = ' Barotropic eq of state: '

 write(iprint,"(a)") ' '
 if (polyk2 > 0.0) then
    write(iprint,"(/,2a,2(es10.3,a))") baro, 'cs_ld            = ',sqrt(polyk2),cu,sqrt(polyk2)*unit_velocity,' cm/s'
 endif
 write(iprint,"(  2a,2(es10.3,a))")    baro, 'cs               = ',sqrt(polyk), cu,sqrt(polyk)*unit_velocity, ' cm/s'
 if (drhocrit0 > 0.0) then
    write(iprint,"(  2a,2(es10.3,a))") baro, 'rhocritT == rhoT = ',rhocritT,    cu,rhocritT*unit_density,     ' g/cm^3'
    write(iprint,"(  2a,2(es10.3,a))") baro, 'rhocrit0 == rho0 = ',rhocrit0,    cu,rhocrit0*unit_density,     ' g/cm^3'
 else
    if (rhocritT > 0.0) then
       write(iprint,"(2a,2(es10.3,a))")baro, 'rhocrit0 == rho0 = ',rhocritT,    cu,rhocritT*unit_density,     ' g/cm^3'
    endif
 endif
 write(iprint,"(  2a,2(es10.3,a))")    baro, 'rhocrit1 == rho1 = ',rhocrit1,    cu,rhocrit1*unit_density,     ' g/cm^3'
 write(iprint,"(  2a,2(es10.3,a))")    baro, 'rhocrit2 == rho2 = ',rhocrit2,    cu,rhocrit2*unit_density,     ' g/cm^3'
 write(iprint,"(  2a,2(es10.3,a))")    baro, 'rhocrit3 == rho3 = ',rhocrit3,    cu,rhocrit3*unit_density,     ' g/cm^3'
 write(iprint,"(  2a,2(es10.3,a))")    baro, 'rhocrit4 == rho4 = ',rhocrit4,    cu,rhocrit4*unit_density,     ' g/cm^3'
 write(iprint,"(a)")                   baro
 if (drhocrit0 > 0.0) then
    write(iprint,"(2a)")       baro, '        rho < rhoT: P = cs_ld^2*rho'
    write(iprint,"(2a)")       baro, 'rhoT <= rho < rho0: P = 10**(log10(cs_ld^2) + M*(log10(rhoT)-log10(rho)))'
 else
    if (polyk2 > 0.0) then
       write(iprint,"(2a)")    baro, '        rho < rho0: P = cs_ld^2*rho'
    endif
 endif
 if (polyk2 > 0.0) then
    write(iprint,"(2a)")       baro, 'rho0 <= rho < rho1: P = cs^2*rho'
 else
    write(iprint,"(2a)")       baro, '        rho < rho1: P = cs^2*rho'
 endif
 write(iprint,"(2a,f5.3)")     baro, 'rho1 <= rho < rho2: P = cs^2*rho1*(rho /rho1)^',gamma1
 write(iprint,"(a,2(a,f5.3))") baro, 'rho2 <= rho < rho3: P = cs^2*rho1*(rho2/rho1)^',gamma1,'*(rho /rho2)^',gamma2
 write(iprint,"(a,3(a,f5.3))") baro, 'rho3 <= rho < rho4: P = cs^2*rho1*(rho2/rho1)^',gamma1,'*(rho3/rho2)^',gamma2, &
                                                                      '*(rho /rho3)^',gamma3
 write(iprint,"(a,4(a,f5.3))") baro, 'rho4 <= rho:        P = cs^2*rho1*(rho2/rho1)^',gamma1,'*(rho3/rho2)^',gamma2, &
                                                                      '*(rho4/rho3)^',gamma3,'*(rho /rho4)^',gamma4

end subroutine eos_info_barotropic

!-----------------------------------------------------------------------
!+
!  internal routine to verify that val1 < val2
!+
!-----------------------------------------------------------------------
subroutine verify_less_than(ierr,val1,val2)
 use io, only:error
 integer, intent(inout) :: ierr
 real,    intent(in)    :: val1,val2

 if (val1 > val2) then
    ierr = ierr + 1
    call error('eos_barotropic','incorrect ordering of rhocrit')
 endif

end subroutine verify_less_than
!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_barotropic(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(drhocrit0,  'drhocrit','transition size between rhocrit0 & 1 (fraction of rhocrit0; barotropic eos)',iunit)
 call write_inopt(rhocrit0cgs,'rhocrit0','critical density 0 in g/cm^3 (barotropic eos)',iunit)
 call write_inopt(rhocrit1cgs,'rhocrit1','critical density 1 in g/cm^3 (barotropic eos)',iunit)
 call write_inopt(rhocrit2cgs,'rhocrit2','critical density 2 in g/cm^3 (barotropic eos)',iunit)
 call write_inopt(rhocrit3cgs,'rhocrit3','critical density 3 in g/cm^3 (barotropic eos)',iunit,exp=.true.)
 call write_inopt(rhocrit4cgs,'rhocrit4','critical density 4 in g/cm^3 (barotropic eos)',iunit,exp=.true.)
 call write_inopt(gamma1,'gamma1','adiabatic index 1 (barotropic eos)',iunit)
 call write_inopt(gamma2,'gamma2','adiabatic index 2 (barotropic eos)',iunit)
 call write_inopt(gamma3,'gamma3','adiabatic index 3 (barotropic eos)',iunit)
 call write_inopt(gamma4,'gamma4','adiabatic index 4 (barotropic eos)',iunit)

end subroutine write_options_eos_barotropic

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_eos_barotropic(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot  = 0
 character(len=30), parameter  :: label = 'eos_barotropic'

 imatch  = .true.
 select case(trim(name))
 case('drhocrit')
    read(valstring,*,iostat=ierr) drhocrit0
    if (drhocrit0 < 0.)  call fatal(label,'drhocrit0 < 0: Negative transition region is nonsense')
    if (drhocrit0 > 1.)  call fatal(label,'drhocrit0 > 1: Too large of transition region')
    ngot = ngot + 1
 case('rhocrit0')
    read(valstring,*,iostat=ierr) rhocrit0cgs
    ! if (rhocrit0cgs <= 0.) call fatal(label,'rhocrit0 <= 0')  ! This region can be 0 if the warm medium is undefined
    ngot = ngot + 1
 case('rhocrit1')
    read(valstring,*,iostat=ierr) rhocrit1cgs
    if (rhocrit1cgs <= 0.) call fatal(label,'rhocrit1 <= 0')
    ngot = ngot + 1
 case('rhocrit2')
    read(valstring,*,iostat=ierr) rhocrit2cgs
    if (rhocrit2cgs <= 0.) call fatal(label,'rhocrit2 <= 0')
    ngot = ngot + 1
 case('rhocrit3')
    read(valstring,*,iostat=ierr) rhocrit3cgs
    if (rhocrit3cgs <= 0.) call fatal(label,'rhocrit3 <= 0')
    ngot = ngot + 1
 case('rhocrit4')
    read(valstring,*,iostat=ierr) rhocrit4cgs
    if (rhocrit4cgs <= 0.) call fatal(label,'rhocrit4 <= 0')
    ngot = ngot + 1
 case('gamma1')
    read(valstring,*,iostat=ierr) gamma1
    if (gamma1 < 1.) call fatal(label,'gamma1 < 1.0')
    ngot = ngot + 1
 case('gamma2')
    read(valstring,*,iostat=ierr) gamma2
    if (gamma2 < 1.) call fatal(label,'gamma2 < 1.0')
    ngot = ngot + 1
 case('gamma3')
    read(valstring,*,iostat=ierr) gamma3
    if (gamma3 < 1.) call fatal(label,'gamma3 < 1.0')
    ngot = ngot + 1
 case('gamma4')
    read(valstring,*,iostat=ierr) gamma4
    if (gamma4 < 1.) call fatal(label,'gamma4 < 1.0')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 8)

end subroutine read_options_eos_barotropic

end module eos_barotropic
