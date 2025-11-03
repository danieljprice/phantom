!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setorbit
!
! Generic procedure for setting up two body orbits with
! different parameter sets for the orbital elements
!
! The current options are:
!   0) Campbell elements for bound or unbound orbit (aeiOwf)
!   1) Flyby parameters (periapsis, initial separation, i0we)
!   2) position and velocity for both bodies
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, orbits, physcon, setbinary, units
!
 implicit none
 public :: set_orbit
 public :: set_defaults_orbit,write_options_orbit,read_options_orbit
 public :: orbit_t
!
 ! define data types with options needed
 ! to setup an orbit
 !
 integer, parameter :: len_str = 20

 type campbell_elems
    character(len=len_str) :: a ! semi-major axis as a string because can specify units
 end type campbell_elems

 type posvel_elems
    character(len=len_str) :: x1(3)  ! position of body 1
    character(len=len_str) :: v1(3)  ! velocity of body 1
    character(len=len_str) :: x2(3)  ! position of body 2
    character(len=len_str) :: v2(3)  ! velocity of body 2
 end type posvel_elems

 type flyby_elems
    character(len=len_str) :: rp   ! pericentre distance in arbitrary units
    character(len=len_str) :: d    ! initial separation in arbitrary units
 end type flyby_elems

 type obs_elems
    character(len=len_str) :: dx(3)
    character(len=len_str) :: dv(3)
    real :: f  ! true anomaly at moment of observation
 end type obs_elems

 !
 ! generic type handling all options
 !
 type orbit_t
    integer :: input_type
    real :: a,e,i,O,w,f
    type(campbell_elems) :: elems
    type(flyby_elems)    :: flyby
    type(obs_elems)      :: obs
    type(posvel_elems)   :: posvel
    real :: period
 end type orbit_t

 private

contains

!----------------------------------------------------------------
!+
!  default parameters for orbit type
!+
!----------------------------------------------------------------
subroutine set_defaults_orbit(orbit)
 type(orbit_t), intent(out) :: orbit

 orbit%input_type = 0
 orbit%elems%a = '10.'
 orbit%a = 10.0
 orbit%e = 0.0
 orbit%i = 0.0
 orbit%O = 0.0
 orbit%w = 270.  ! argument of periapsis
 orbit%f = 180.  ! start orbit at apocentre

 orbit%flyby%rp = '200.'
 orbit%flyby%d = '2000.0'

 orbit%obs%dx(1) = '20.0'
 orbit%obs%dx(2) = '0.0'
 orbit%obs%dx(3) = '0.0'
 orbit%obs%dv(1) = '0.0'
 orbit%obs%dv(2) = '2.0'
 orbit%obs%dv(3) = '0.0'
 orbit%obs%f = 0.

 orbit%posvel%x1 = '0.0'
 orbit%posvel%v1 = '0.0'
 orbit%posvel%x2 = '0.0'
 orbit%posvel%v2 = '0.0'
 orbit%posvel%x1(1) = '10.0'
 orbit%posvel%x2(1) = '-10.0'
 orbit%posvel%v1(2) = '1.0'
 orbit%posvel%v2(2) = '-1.0'

end subroutine set_defaults_orbit

!----------------------------------------------------------------
!+
!  setup for two body orbit
!+
!----------------------------------------------------------------
subroutine set_orbit(orbit,m1,m2,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,verbose,ierr,omega_corotate)
 use units,     only:in_code_units
 use setbinary, only:set_binary
 use orbits,    only:convert_flyby_to_elements,convert_posvel_to_flyby,get_orbital_elements
 type(orbit_t), intent(inout) :: orbit
 real,          intent(in)    :: m1,m2,hacc1,hacc2
 real,          intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer,       intent(inout) :: nptmass
 logical,       intent(in)    :: verbose
 integer,       intent(out)   :: ierr
 real,          intent(out), optional :: omega_corotate
 real :: x1(3),x2(3),v1(3),v2(3)
 integer :: i

 ierr = 0
 select case(orbit%input_type)
 case(3)
    !
    ! for posvel input we do not call set_binary at all, just set x and v for both bodies.
    ! However, the result should be the same, just without the centre of mass being at zero
    !
    do i=1,3
       x1(i) = in_code_units(orbit%posvel%x1(i),ierr,unit_type='length')
       x2(i) = in_code_units(orbit%posvel%x2(i),ierr,unit_type='length')
       v1(i) = in_code_units(orbit%posvel%v1(i),ierr,unit_type='velocity')
       v2(i) = in_code_units(orbit%posvel%v2(i),ierr,unit_type='velocity')
    enddo
    ! body 1
    xyzmh_ptmass(1:3,nptmass+1) = x1
    xyzmh_ptmass(4,nptmass+1)   = m1
    xyzmh_ptmass(5,nptmass+1)   = hacc1
    vxyz_ptmass(1:3,nptmass+1)  = v1
    ! body 2
    xyzmh_ptmass(1:3,nptmass+2) = x2
    xyzmh_ptmass(4,nptmass+2)   = m2
    xyzmh_ptmass(5,nptmass+2)   = hacc2
    vxyz_ptmass(1:3,nptmass+2)  = v2

    ! set orbit elements from position and velocity, so they can be used to compute the period
    call set_orbit_elements(orbit,m1,m2,verbose=verbose)
 case default
    call set_orbit_elements(orbit,m1,m2,verbose=verbose)
    !
    !--now setup orbit using sink particles
    !
    if (present(omega_corotate)) then
       call set_binary(m1,m2,orbit%a,orbit%e,hacc1,hacc2,&
                       xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                       posang_ascnode=orbit%O,arg_peri=orbit%w,&
                       incl=orbit%i,f=orbit%f,verbose=verbose)
    else
       call set_binary(m1,m2,orbit%a,orbit%e,hacc1,hacc2,&
                       xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                       posang_ascnode=orbit%O,arg_peri=orbit%w,&
                       incl=orbit%i,f=orbit%f,verbose=verbose)
    endif
 end select

 call get_orbital_time(orbit,m1,m2,orbit%period)

end subroutine set_orbit

!----------------------------------------------------------------------
!+
!  convert input parameters to standard orbital elements (a,e,i,O,w,f)
!+
!----------------------------------------------------------------------
subroutine set_orbit_elements(orbit,m1,m2,verbose)
 use physcon, only:days
 use units,   only:in_code_units,is_time_unit,utime
 use orbits,  only:get_semimajor_axis,convert_flyby_to_elements,get_orbital_elements,&
                   get_pericentre_distance,get_true_anomaly_from_separation,rad_to_deg,&
                   get_time_to_separation,get_time_between_true_anomalies
 type(orbit_t), intent(inout) :: orbit
 real,          intent(in)    :: m1,m2
 logical,       intent(in), optional :: verbose
 integer :: ierr,i
 logical :: do_verbose
 real :: dx(3),dv(3),mu,rp,d,time_to_obs
 real :: x1(3),x2(3),v1(3),v2(3)

 ierr = 0
 do_verbose = .true.
 if (present(verbose)) do_verbose = verbose

 mu = m1+m2
 select case(orbit%input_type)
 case(3)
    ! convert input strings to code units
    do i=1,3
       x1(i) = in_code_units(orbit%posvel%x1(i),ierr,unit_type='length')
       x2(i) = in_code_units(orbit%posvel%x2(i),ierr,unit_type='length')
       v1(i) = in_code_units(orbit%posvel%v1(i),ierr,unit_type='velocity')
       v2(i) = in_code_units(orbit%posvel%v2(i),ierr,unit_type='velocity')
    enddo
    ! for posvel input, compute orbital elements from relative position and velocity
    dx = x2 - x1
    dv = v2 - v1
    call get_orbital_elements(mu,dx,dv,orbit%a,orbit%e,orbit%i,orbit%O,orbit%w,orbit%f)
 case(2)
    ! convert input strings to code units
    do i=1,3
       dx(i) = in_code_units(orbit%obs%dx(i),ierr,unit_type='length')
       dv(i) = in_code_units(orbit%obs%dv(i),ierr,unit_type='velocity')
    enddo
    d = in_code_units(orbit%flyby%d,ierr,unit_type='length')
    if (do_verbose) then
       print "(/,a,/)", ' Flyby Reconstructor^TM Inputs:'
       print*,'          separation in code units = ',dx
       print*,' velocity difference in code units = ',dv
    endif

    ! convert observed separation and velocity to rp,d,e,O,w,i,f
    call get_orbital_elements(mu,dx,dv,orbit%a,orbit%e,orbit%i,orbit%O,orbit%w,orbit%obs%f)

    if (do_verbose) then
       print "(/,a,/)", ' Flyby Reconstructor^TM Recovered Orbital Parameters (at moment of observation):'
       print "(a,g0.4)",'          rp in code units      = ',get_pericentre_distance(mu,dx,dv)
       print "(a,g0.4)",'          semi-major axis a     = ',orbit%a
       print "(a,g0.4)",'          projected separation  = ',sqrt(dot_product(dx(1:2),dx(1:2)))
       print "(a,g0.4)",'          true separation       = ',sqrt(dot_product(dx,dx))
       print "(a,g0.4)",'          eccentricity          = ',orbit%e
       print "(a,g0.4)",'          O (pos. angle, deg)   = ',orbit%O
       print "(a,g0.4)",'          w (arg. peri, deg)    = ',orbit%w
       print "(a,g0.4)",'          i (inclination, deg)  = ',orbit%i
       print "(a,g0.4)",'          f (true anomaly, deg) = ',orbit%obs%f
    endif

    ! convert desired initial separation to initial true anomaly
    orbit%f = get_true_anomaly_from_separation(orbit%a,orbit%e,d)

    ! compute time from initial true anomaly to true anomaly at moment of observation
    time_to_obs = get_time_between_true_anomalies(mu,orbit%a,orbit%e,orbit%f,orbit%obs%f)
    if (time_to_obs < 0.0) then
       orbit%f = -orbit%f
       time_to_obs = get_time_between_true_anomalies(mu,orbit%a,orbit%e,orbit%f,orbit%obs%f)
    endif
    if (do_verbose) then
       print "(a,g0.4)",'  starting true anomaly f (true anomaly, deg) = ',orbit%f
       print "(a,g0.4)",'  time to observation = ',time_to_obs
    endif
 case(1)
    ! flyby elements give pericentre distance and initial separation, convert to reals
    rp = in_code_units(orbit%flyby%rp,ierr)
    d = in_code_units(orbit%flyby%d,ierr)

    ! convert rp,e,d to semi-major axis and initial true anomaly
    call convert_flyby_to_elements(rp,d,orbit%e,orbit%a,orbit%f)
 case default
    ! convert input string for semi-major axis to code units
    orbit%a = in_code_units(orbit%elems%a,ierr)

    ! can also specify period in time units, convert to semi-major axis
    if (is_time_unit(orbit%elems%a) .and. ierr == 0 .and. orbit%e < 1.) then
       orbit%a = -abs(orbit%a)
       if (do_verbose) then
          print "(a,g0,a,g0,a)",' Using PERIOD = ',abs(orbit%a),' = ',abs(orbit%a)*utime/days,' days'
       endif
    endif
    if (orbit%a < 0. .and. orbit%e < 1.) orbit%a = get_semimajor_axis(mu,abs(orbit%a)) ! convert period to semi-major axis
 end select

end subroutine set_orbit_elements

!----------------------------------------------------------------
!+
!  get orbit time from orbital elements. This is the period
!  for a bound orbit. For unbound orbits, we compute a time
!  to reach pericentre
!+
!----------------------------------------------------------------
subroutine get_orbital_time(orbit,m1,m2,period)
 use orbits, only:get_time_between_true_anomalies,get_T_flyby_hyp,&
                  get_T_flyby_par,orbit_is_parabolic,get_orbital_period
 use units,  only:in_code_units
 type(orbit_t), intent(in) :: orbit
 real, intent(in)  :: m1,m2
 real, intent(out) :: period
 integer :: ierr
 real :: mu,flyby_d

 mu = m1+m2
 if (orbit%input_type==2) then
    ! for Flyby Reconstructor^TM input, compute time to reach observed separation
    period = get_time_between_true_anomalies(mu,orbit%a,orbit%e,orbit%f,orbit%obs%f)
 else
    if (orbit%e > 1.0) then
       period = get_T_flyby_hyp(mu,orbit%e,orbit%f,orbit%a)
    elseif (orbit_is_parabolic(orbit%e)) then
       flyby_d = in_code_units(orbit%flyby%d,ierr,unit_type='length')
       period = get_T_flyby_par(mu,orbit%a,flyby_d/orbit%a)
    else
       period = get_orbital_period(mu,orbit%a)
    endif
 endif

end subroutine get_orbital_time

!----------------------------------------------------------------
!+
!  check if separation exceeds apocentre for flyby elements
!+
!----------------------------------------------------------------
logical function sep_in_range(orbit,sep,rp,ra)
 use units,  only:in_code_units
 use orbits, only:orbit_is_parabolic
 type(orbit_t), intent(in) :: orbit
 real, intent(out) :: sep,rp,ra
 real :: a,e,d
 integer :: ierr

 ! must have already called set_orbit_elements
 a = orbit%a
 e = orbit%e
 sep = 0.
 if (orbit_is_parabolic(e)) then
    rp = a
    ra = huge(ra)
 elseif (e > 1.) then
    rp = a*(1. - e)
    ra = huge(ra)
 else
    rp = a*(1. - e)
    ra = a*(1. + e)
 endif

 sep_in_range = .true.
 select case(orbit%input_type)
 case(1,2)
    d = in_code_units(orbit%flyby%d,ierr)
    if (e < 1. .and. d > ra) sep_in_range = .false.
    if (d < rp) sep_in_range = .false.
 end select

end function sep_in_range

!----------------------------------------------------------------
!+
!  write options to .setup file
!+
!----------------------------------------------------------------
subroutine write_options_orbit(orbit,iunit,label,prefix,comment_prefix,input_type)
 use infile_utils, only:write_inopt
 type(orbit_t),    intent(in) :: orbit
 integer,          intent(in) :: iunit
 character(len=*), intent(in), optional :: label,prefix,comment_prefix
 integer,          intent(in), optional :: input_type
 character(len=10) :: c,p
 character(len=20) :: cp
 integer :: itype

 ! append optional label e.g. '1', '2' or prefix e.g. 'binary' and comment prefix e.g. 'outer binary:'
 c = ''; p = ''; cp = ''
 if (present(label)) c = trim(adjustl(label))
 if (present(prefix)) p = '_'//trim(adjustl(prefix))
 if (present(comment_prefix)) cp = trim(adjustl(comment_prefix))//':'

 write(iunit,"(/,a)") '# orbit '//trim(c)
 if (.not.present(input_type)) then
    itype = orbit%input_type
    call write_inopt(orbit%input_type,'itype'//trim(p)//trim(c),'type of orbital elements (0=aeiOwf,1=flyby,2=obs,3=posvel)',iunit)
 else
    itype = input_type
 endif
 if (present(prefix)) p = trim(adjustl(prefix))//'_'

 select case(itype)
 case(3)
    call write_inopt(orbit%posvel%x1(1),trim(p)//'x1'//trim(c),trim(cp)//'x position body 1 (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%posvel%x1(2),trim(p)//'y1'//trim(c),trim(cp)//'y position body 1 (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%posvel%x1(3),trim(p)//'z1'//trim(c),trim(cp)//'z position body 1 (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%posvel%v1(1),trim(p)//'vx1'//trim(c),trim(cp)//'x velocity body 1 (code units or e.g. 1*km/s)',iunit)
    call write_inopt(orbit%posvel%v1(2),trim(p)//'vy1'//trim(c),trim(cp)//'y velocity body 1 (code units or e.g. 1*km/s)',iunit)
    call write_inopt(orbit%posvel%v1(3),trim(p)//'vz1'//trim(c),trim(cp)//'z velocity body 1 (code units or e.g. 1*km/s)',iunit)
    call write_inopt(orbit%posvel%x2(1),trim(p)//'x2'//trim(c),trim(cp)//'x position body 2 (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%posvel%x2(2),trim(p)//'y2'//trim(c),trim(cp)//'y position body 2 (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%posvel%x2(3),trim(p)//'z2'//trim(c),trim(cp)//'z position body 2 (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%posvel%v2(1),trim(p)//'vx2'//trim(c),trim(cp)//'x velocity body 2 (code units or e.g. 1*km/s)',iunit)
    call write_inopt(orbit%posvel%v2(2),trim(p)//'vy2'//trim(c),trim(cp)//'y velocity body 2 (code units or e.g. 1*km/s)',iunit)
    call write_inopt(orbit%posvel%v2(3),trim(p)//'vz2'//trim(c),trim(cp)//'z velocity body 2 (code units or e.g. 1*km/s)',iunit)
 case(2)
    call write_inopt(orbit%obs%dx(1),trim(p)//'dx'//trim(c),trim(cp)//'observed dx at t=tmax (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%obs%dx(2),trim(p)//'dy'//trim(c),trim(cp)//'observed dy at t=tmax (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%obs%dx(3),trim(p)//'dz'//trim(c),trim(cp)//'[guessed] dz at t=tmax (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%obs%dv(1),trim(p)//'dvx'//trim(c),trim(cp)//'[guessed] dvx at t=tmax (code units or e.g. 1*km/s)',iunit)
    call write_inopt(orbit%obs%dv(2),trim(p)//'dvy'//trim(c),trim(cp)//'[guessed] dvy at t=tmax (code units or e.g. 1*km/s)',iunit)
    call write_inopt(orbit%obs%dv(3),trim(p)//'dvz'//trim(c),trim(cp)//'observed dvz at t=tmax (code units or e.g. 1*km/s)',iunit)
    call write_inopt(orbit%flyby%d,trim(p)//'d'//trim(c),trim(cp)//'separation at t=0 if unbound',iunit)
 case(1)
    call write_inopt(orbit%flyby%rp,trim(p)//'rp'//trim(c),trim(cp)//'pericentre distance (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%flyby%d,trim(p)//'d'//trim(c),trim(cp)//'initial separation (code units or e.g. 1*au)',iunit)
    call write_inopt(orbit%O,trim(p)//'O'//trim(c),trim(cp)//'position angle of the ascending node (deg)',iunit)
    call write_inopt(orbit%i,trim(p)//'i'//trim(c),trim(cp)//'inclination (deg)',iunit)
    call write_inopt(orbit%w,trim(p)//'w'//trim(c),trim(cp)//'argument of periapsis (deg)',iunit)
    call write_inopt(orbit%e,trim(p)//'e'//trim(c),trim(cp)//'eccentricity',iunit)
 case default
    call write_inopt(orbit%elems%a,trim(p)//'a'//trim(c),&
                     trim(cp)//'semi-major axis (e.g. 1 au), period (e.g. 10*days) or rp if e>=1',iunit)
    call write_inopt(orbit%e,trim(p)//'e'//trim(c),trim(cp)//'eccentricity',iunit)
    call write_inopt(orbit%i,trim(p)//'i'//trim(c),trim(cp)//'inclination (deg)',iunit)
    call write_inopt(orbit%O,trim(p)//'O'//trim(c),trim(cp)//'position angle of ascending node (deg)',iunit)
    call write_inopt(orbit%w,trim(p)//'w'//trim(c),trim(cp)//'argument of periapsis (deg)',iunit)
    call write_inopt(orbit%f,trim(p)//'f'//trim(c),trim(cp)//'initial true anomaly (180=apoastron)',iunit)
 end select

end subroutine write_options_orbit

!----------------------------------------------------------------
!+
!  read options from .setup file
!+
!----------------------------------------------------------------
subroutine read_options_orbit(orbit,m1,m2,db,nerr,label,prefix,input_type)
 use infile_utils, only:inopts,read_inopt
 type(orbit_t),             intent(out)   :: orbit
 real,                      intent(in)    :: m1,m2
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(inout) :: nerr
 character(len=*),          intent(in), optional :: label,prefix
 integer,                   intent(in), optional :: input_type
 character(len=10) :: c,p
 real :: sep,rp,ra

 ! append optional label e.g. '1', '2' or prefix e.g. 'binary' and comment prefix e.g. 'outer binary:'
 c = ''; p = ''
 if (present(label)) c = trim(adjustl(label))
 if (present(prefix)) p = '_'//trim(adjustl(prefix))

 call set_defaults_orbit(orbit)

 ! can force read of a particular input type, if provided as argument
 if (present(input_type)) then
    orbit%input_type = input_type
 else
    call read_inopt(orbit%input_type,'itype'//trim(p)//trim(c),db,errcount=nerr,min=0,max=3)
 endif
 if (present(prefix)) p = trim(adjustl(prefix))//'_'

 select case(orbit%input_type)
 case(3)
    call read_inopt(orbit%posvel%x1(1),trim(p)//'x1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x1(2),trim(p)//'y1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x1(3),trim(p)//'z1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v1(1),trim(p)//'vx1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v1(2),trim(p)//'vy1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v1(3),trim(p)//'vz1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x2(1),trim(p)//'x2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x2(2),trim(p)//'y2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x2(3),trim(p)//'z2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v2(1),trim(p)//'vx2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v2(2),trim(p)//'vy2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v2(3),trim(p)//'vz2'//trim(c),db,errcount=nerr)
 case(2)
    call read_inopt(orbit%obs%dx(1),trim(p)//'dx'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%obs%dx(2),trim(p)//'dy'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%obs%dx(3),trim(p)//'dz'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%obs%dv(1),trim(p)//'dvx'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%obs%dv(2),trim(p)//'dvy'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%obs%dv(3),trim(p)//'dvz'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%flyby%d,trim(p)//'d'//trim(c),db,errcount=nerr)
 case(1)
    call read_inopt(orbit%flyby%rp,trim(p)//'rp'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%flyby%d,trim(p)//'d'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%e,trim(p)//'e'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%i,trim(p)//'i'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%O,trim(p)//'O'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%w,trim(p)//'w'//trim(c),db,errcount=nerr)
 case default
    call read_inopt(orbit%elems%a,trim(p)//'a'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%e,trim(p)//'e'//trim(c),db,min=0.,errcount=nerr)
    call read_inopt(orbit%i,trim(p)//'i'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%O,trim(p)//'O'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%w,trim(p)//'w'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%f,trim(p)//'f'//trim(c),db,errcount=nerr)
 end select

 ! convert input parameters to standard orbital elements (a,e,i,O,w,f)
 call set_orbit_elements(orbit,m1,m2,verbose=.false.)
 if (.not.sep_in_range(orbit,sep,rp,ra)) then
    nerr = nerr + 1
    print "(4(a,1pg10.3))",' ERROR: initial distance ',sep,' out of range, need d >= ',rp,' and d <= ',ra,' for e=',orbit%e
 endif

end subroutine read_options_orbit

end module setorbit
