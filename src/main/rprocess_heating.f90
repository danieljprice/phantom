!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: rprocess_heating
!
!  DESCRIPTION:
!   Implements r-process heating.
!
!  OWNER: Siva Darbha
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!
!    q_0_cgs -- heating rate coefficient, in [ergs/s/g]
!    t_b1_seconds -- time of break 1 in heating rate, in [s]'
!    exp_1 -- exponent 1 in radioactive power component
!    t_b2_seconds -- time of break 2 in heating rate, in [s]'
!    exp_2 -- exponent 2 in radioactive power component
!    t_start_seconds -- time after merger at which the simulation begins [s]
!
!  DEPENDENCIES: infile_utils, io, physcon, units
!+
!--------------------------------------------------------------------------
module rprocess_heating
 implicit none
 integer, parameter :: maxt = 1000
 real    :: temper(maxt),lambda(maxt),slope(maxt),yfunc(maxt)
 integer :: nt
 !
 ! set default values for input parameters
 !
 real :: q_0_cgs    = 0 ! heating rate coefficient [ergs/s/g]
 real :: t_b1_seconds = 0 ! time of break 1 in heating rate [s]
 real :: exp_1 = 0 ! exponent 1 in heating rate
 real :: t_b2_seconds = 0 ! time of break 2 in heating rate [s]
 real :: exp_2 = 0 ! exponent 2 in heating rate
 real :: t_start_seconds = 0 ! time after merger at which the simulation begins [s]

 public :: write_options_rprocess,read_options_rprocess
 public :: get_rprocess_heating_rate, energ_rprocess

 private

contains

!-----------------------------------------------------------------------
!+
!  get the r-process specific heating rate q at time t
!+
!-----------------------------------------------------------------------
subroutine get_rprocess_heating_rate(q,t)
 use physcon, only:days
 use units,   only:unit_ergg,utime
 real, intent(in)    :: t
 real, intent(out)   :: q
 real :: t_seconds, t_days
 real :: t_b1_days, t_b2_days
 real :: t_start_days, t_new_days
 real :: q_cgs

 t_seconds = t * utime
 t_days = t_seconds / days

 t_b1_days = t_b1_seconds / days
 t_b2_days = t_b2_seconds / days

 t_start_days = t_start_seconds / days
 t_new_days = t_days + t_start_days

 !----- Heating rate in [ergs/s/g]
 if (t_new_days < t_b1_days) then
    q_cgs = q_0_cgs
 else if ( (t_b1_days <= t_new_days) .and. (t_new_days < t_b2_days) ) then
    q_cgs = q_0_cgs * (t_days/t_b1_days)**exp_1
 else
    q_cgs = q_0_cgs * (t_b2_days/t_b1_days)**exp_1 * (t_days/t_b2_days)**exp_2
 endif

 q = q_cgs / (unit_ergg/utime)

end subroutine get_rprocess_heating_rate

!-----------------------------------------------------------------------
!+
!  get the change in the entropy variable K in the timestep t -> t+dt
!+
!-----------------------------------------------------------------------
subroutine energ_rprocess(entrop_var,rho,t,dt)
 use eos,     only:gamma
 use physcon, only:days
 use units,   only:unit_ergg,utime
 real, intent(in)    :: rho,t,dt
 real, intent(inout) :: entrop_var
 real :: dt_seconds
 real :: q, q_cgs
 real :: delta_u, delta_u_cgs
 real :: delta_entrop_var

 dt_seconds = dt * utime

 call get_rprocess_heating_rate(q, t)

 q_cgs = q * (unit_ergg/utime)

 !----- Change in the specific internal energy [ergs/g] during the time step
 delta_u_cgs = q_cgs * dt_seconds
 delta_u = delta_u_cgs / unit_ergg

 !----- Change in the entropy variable K during the time step
 !----- In Liptai & Price (2019), K is given the symbol K. -----
 delta_entrop_var =  (gamma - 1.) * rho**(1. - gamma) * delta_u
 entrop_var = entrop_var + delta_entrop_var

end subroutine energ_rprocess

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_rprocess(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(q_0_cgs,'q_0_cgs','heating rate coefficient [ergs/s/g]',iunit)
 call write_inopt(t_b1_seconds,'t_b1_seconds','time of break 1 in heating rate [s]',iunit)
 call write_inopt(exp_1,'exp_1','exponent 1 in heating rate',iunit)
 call write_inopt(t_b2_seconds,'t_b2_seconds','time of break 2 in heating rate [s]',iunit)
 call write_inopt(exp_2,'exp_2','exponent 2 in heating rate',iunit)
 call write_inopt(t_start_seconds,'t_start_seconds','time after merger at which the simulation begins [s]',iunit)

end subroutine write_options_rprocess

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_rprocess(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .false.  ! cooling options are compulsory
 select case(trim(name))
 case('q_0_cgs')
    read(valstring,*,iostat=ierr) q_0_cgs
    ngot = ngot + 1
 case('t_b1_seconds')
    read(valstring,*,iostat=ierr) t_b1_seconds
    ngot = ngot + 1
 case('exp_1')
    read(valstring,*,iostat=ierr) exp_1
    ngot = ngot + 1
 case('t_b2_seconds')
    read(valstring,*,iostat=ierr) t_b2_seconds
    ngot = ngot + 1
 case('exp_2')
    read(valstring,*,iostat=ierr) exp_2
    ngot = ngot + 1
 case('t_start_seconds')
    read(valstring,*,iostat=ierr) t_start_seconds
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if (ngot >= 6) igotall = .true.

end subroutine read_options_rprocess

end module rprocess_heating