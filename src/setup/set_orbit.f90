!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
!   1) Flyby parameters (periapsis, initial separation, argument of periapsis, inclination)
!   2) position and velocity for both bodies
!
! While Campbell elements can be used for unbound orbits, they require
! specifying the true anomaly at the start of the simulation. This is
! not always easy to determine, so the flyby option is provided as an
! alternative. There one specifies the initial separation instead, however
! the choice of angles is more restricted
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, physcon, setbinary, setflyby, units
!
 implicit none
 public :: set_orbit
 public :: set_defaults_orbit,write_options_orbit,read_options_orbit
 public :: orbit_t
!
 ! define data types with options needed
 ! to setup an orbit
 !
 type campbell_elems
    character(len=20) :: semi_major_axis ! string because can specific units
    real :: e   ! eccentricity
    real :: i   ! inclination
    real :: O   ! position angle of the ascending node
    real :: w   ! argument of periapsis
    real :: f   ! initial true anomaly
 end type campbell_elems

 type posvel_elems
    real :: x1(3)  ! position of body 1
    real :: v1(3)  ! velocity of body 1
    real :: x2(3)  ! position of body 2
    real :: v2(3)  ! velocity of body 2
 end type posvel_elems

 type flyby_elems
    character(len=20) :: rp   ! pericentre distance in arbitrary units
    real :: d                 ! initial separation
    real :: O                 ! position angle of the ascending node
    real :: i                 ! inclination
 end type flyby_elems

 !
 ! generic type handling all options
 !
 type orbit_t
    integer :: itype
    type(campbell_elems) :: elems
    type(flyby_elems)    :: flyby
    type(posvel_elems)   :: posvel
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

 orbit%itype = 0
 orbit%elems%semi_major_axis = '10.'
 orbit%elems%e = 0.0
 orbit%elems%i = 0.0
 orbit%elems%O = 0.0
 orbit%elems%w = 270.  ! argument of periapsis
 orbit%elems%f = 180.  ! start orbit at apocentre

 orbit%flyby%rp = '10.'
 orbit%flyby%d = 100.0
 orbit%flyby%O = 0.0
 orbit%flyby%i = 0.0

 orbit%posvel%x1 = 0.0
 orbit%posvel%v1 = 0.0
 orbit%posvel%x2 = 0.0
 orbit%posvel%v2 = 0.0
 orbit%posvel%x1(1) = 10.0
 orbit%posvel%x2(1) = -10.0
 orbit%posvel%v1(2) = 1.0
 orbit%posvel%v2(2) = -1.0

end subroutine set_defaults_orbit

!----------------------------------------------------------------
!+
!  setup for two body orbit
!+
!----------------------------------------------------------------
subroutine set_orbit(orbit,m1,m2,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,verbose,ierr,omega_corotate)
 use physcon,   only:days
 use units,     only:in_code_units,is_time_unit,utime
 use setbinary, only:set_binary,get_a_from_period
 use setflyby,  only:set_flyby
 type(orbit_t), intent(in)    :: orbit
 real,          intent(in)    :: m1,m2,hacc1,hacc2
 real,          intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer,       intent(inout) :: nptmass
 logical,       intent(in)    :: verbose
 integer,       intent(out)   :: ierr
 real,          intent(out), optional :: omega_corotate
 real :: rp,a

 ierr = 0
 select case(orbit%itype)
 case(2)
    ! body 1
    xyzmh_ptmass(1:3,nptmass+1) = orbit%posvel%x1(1:3)
    xyzmh_ptmass(4,nptmass+1)   = m1
    xyzmh_ptmass(5,nptmass+1)   = hacc1
    vxyz_ptmass(1:3,nptmass+1)  = orbit%posvel%v1(1:3)
    ! body 2
    xyzmh_ptmass(1:3,nptmass+2) = orbit%posvel%x2(1:3)
    xyzmh_ptmass(4,nptmass+2)   = m2
    xyzmh_ptmass(5,nptmass+2)   = hacc2
    vxyz_ptmass(1:3,nptmass+2)  = orbit%posvel%v2(1:3)
 case(1)
    rp = in_code_units(orbit%flyby%rp,ierr)

    call set_flyby(m1,m2,rp,orbit%flyby%d,hacc1,hacc2,xyzmh_ptmass, &
                   vxyz_ptmass,nptmass,ierr,orbit%flyby%O,orbit%flyby%i,verbose=verbose)
 case default
    !
    !--if a is negative or is given time units, interpret this as a period
    !
    a = in_code_units(orbit%elems%semi_major_axis,ierr)
    if (is_time_unit(orbit%elems%semi_major_axis) .and. ierr == 0) then
       a = -abs(a)
       print "(a,g0,a,g0,a)",' Using PERIOD = ',abs(a),' = ',abs(a)*utime/days,' days'
    endif
    if (a < 0.) a = get_a_from_period(m1,m2,abs(a))
    !
    !--now setup orbit using sink particles
    !
    if (present(omega_corotate)) then
       call set_binary(m1,m2,a,orbit%elems%e,hacc1,hacc2,&
                       xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                       posang_ascnode=orbit%elems%O,arg_peri=orbit%elems%w,&
                       incl=orbit%elems%i,f=orbit%elems%f,verbose=verbose)
    else
       call set_binary(m1,m2,a,orbit%elems%e,hacc1,hacc2,&
                       xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                       posang_ascnode=orbit%elems%O,arg_peri=orbit%elems%w,&
                       incl=orbit%elems%i,f=orbit%elems%f,verbose=verbose)
    endif
 end select

end subroutine set_orbit

!----------------------------------------------------------------
!+
!  write options to .setup file
!+
!----------------------------------------------------------------
subroutine write_options_orbit(orbit,iunit,label)
 use infile_utils, only:write_inopt
 type(orbit_t), intent(in) :: orbit
 integer,       intent(in) :: iunit
 character(len=*), intent(in), optional :: label
 character(len=10) :: c

 ! append optional label e.g. '1', '2'
 c = ''
 if (present(label)) c = trim(adjustl(label))

 write(iunit,"(/,a)") '# orbit '//trim(c)
 call write_inopt(orbit%itype,'itype'//trim(c),'type of orbital elements (0=aeiOwf,1=flyby,2=posvel)',iunit)
 select case(orbit%itype)
 case(2)
    call write_inopt(orbit%posvel%x1(1),'x1'//trim(c),'x position body 1',iunit)
    call write_inopt(orbit%posvel%x1(2),'y1'//trim(c),'y position body 1',iunit)
    call write_inopt(orbit%posvel%x1(3),'z1'//trim(c),'z position body 1',iunit)
    call write_inopt(orbit%posvel%v1(1),'vx1'//trim(c),'x velocity body 1',iunit)
    call write_inopt(orbit%posvel%v1(2),'vy1'//trim(c),'y velocity body 1',iunit)
    call write_inopt(orbit%posvel%v1(3),'vz1'//trim(c),'z velocity body 1',iunit)
    call write_inopt(orbit%posvel%x2(1),'x2'//trim(c),'x position body 2',iunit)
    call write_inopt(orbit%posvel%x2(2),'y2'//trim(c),'y position body 2',iunit)
    call write_inopt(orbit%posvel%x2(3),'z2'//trim(c),'z position body 2',iunit)
    call write_inopt(orbit%posvel%v2(1),'vx2'//trim(c),'x velocity body 2',iunit)
    call write_inopt(orbit%posvel%v2(2),'vy2'//trim(c),'y velocity body 2',iunit)
    call write_inopt(orbit%posvel%v2(3),'vz2'//trim(c),'z velocity body 2',iunit)
 case(1)
    call write_inopt(orbit%flyby%rp,'rp'//trim(c),'pericentre distance',iunit)
    call write_inopt(orbit%flyby%d,'d'//trim(c),'initial separation [same units as rp]',iunit)
    call write_inopt(orbit%flyby%O,'O'//trim(c),'position angle of the ascending node',iunit)
    call write_inopt(orbit%flyby%i,'i'//trim(c),'inclination',iunit)
 case default
    call write_inopt(orbit%elems%semi_major_axis,'a'//trim(c),&
                     'semi-major axis (e.g. 1 au), period (e.g. 10*days) or rp if e=1',iunit)
    call write_inopt(orbit%elems%e,'ecc'//trim(c),'eccentricity',iunit)
    call write_inopt(orbit%elems%i,'inc'//trim(c),'inclination (deg)',iunit)
    call write_inopt(orbit%elems%O,'O'//trim(c),'position angle of ascending node (deg)',iunit)
    call write_inopt(orbit%elems%w,'w'//trim(c),'argument of periapsis (deg)',iunit)
    call write_inopt(orbit%elems%f,'f'//trim(c),'initial true anomaly (180=apoastron)',iunit)
 end select

end subroutine write_options_orbit

!----------------------------------------------------------------
!+
!  read options from .setup file
!+
!----------------------------------------------------------------
subroutine read_options_orbit(orbit,db,nerr,label)
 use infile_utils, only:inopts,read_inopt
 type(orbit_t),             intent(out)   :: orbit
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(inout) :: nerr
 character(len=*),          intent(in), optional :: label
 character(len=10) :: c

 ! append optional label e.g. '1', '2'
 c = ''
 if (present(label)) c = trim(adjustl(label))

 call read_inopt(orbit%itype,'itype'//trim(c),db,errcount=nerr,min=0,max=2)
 select case(orbit%itype)
 case(2)
    call read_inopt(orbit%posvel%x1(1),'x1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x1(2),'y1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x1(3),'z1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v1(1),'vx1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v1(2),'vy1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v1(3),'vz1'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x2(1),'x2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x2(2),'y2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%x2(3),'z2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v2(1),'vx2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v2(2),'vy2'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%posvel%v2(3),'vz2'//trim(c),db,errcount=nerr)
 case(1)
    call read_inopt(orbit%flyby%rp,'rp'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%flyby%d,'d'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%flyby%O,'O'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%flyby%i,'i'//trim(c),db,errcount=nerr)
 case default
    call read_inopt(orbit%elems%semi_major_axis,'a'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%elems%e,'ecc'//trim(c),db,min=0.,errcount=nerr)
    call read_inopt(orbit%elems%i,'inc'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%elems%O,'O'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%elems%w,'w'//trim(c),db,errcount=nerr)
    call read_inopt(orbit%elems%f,'f'//trim(c),db,errcount=nerr)
 end select

end subroutine read_options_orbit

end module setorbit
