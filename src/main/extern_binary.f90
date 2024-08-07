!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_binary
!
! This module contains routines relating to the computation
!    of a gravitational force/potential from a central binary system
!    (adapted from a routine by Graham Wynn)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - accradius1 : *accretion radius of primary*
!   - accradius2 : *accretion radius of secondary (if iexternalforce=binary)*
!   - eps_soft1  : *Plummer softening of primary*
!   - eps_soft2  : *Plummer softening of secondary*
!   - mass1      : *m1 of central binary system (if iexternalforce=binary)*
!   - mass2      : *m2 of central binary system (if iexternalforce=binary)*
!   - ramp       : *ramp up mass of secondary over first 5 orbits?*
!
! :Dependencies: dump_utils, infile_utils, io, physcon
!
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, public :: mass1 = 1.0
 real, public :: mass2 = 0.001
 real, public :: accradius1  = 0.5
 real, public :: accradius2  = 0.5
 real, public :: accretedmass1 = 0.
 real, public :: accretedmass2 = 0.
 real, public :: eps_soft1 = 0.0
 real, public :: eps_soft2 = 0.0
 logical, public :: ramp = .false.
 logical, public :: surface_force = .false.
 real,    private :: mass2_temp

 public :: binary_force, binary_posvel, binary_accreted, update_binary
 public :: write_options_externbinary, read_options_externbinary
 public :: write_headeropts_externbinary, read_headeropts_externbinary
 private

 real, private :: x1,y1,x2,y2

 ! Required for HDF5 compatibility
 real, public :: a0 = 0.
 real, public :: direction = 0.

contains

!----------------------------------------------
!+
!  update position of binary at new time
!  only do this ONCE as sin/cos expensive
!+
!----------------------------------------------
subroutine update_binary(ti)
 use physcon, only:pi
 real, intent(in) :: ti
 real :: omega,mtot,a,x,y

 if (ramp .and. ti < 10.*pi) then
    mass2_temp = mass2*(sin(ti/20.)**2)
 else
    mass2_temp = mass2
 endif

 if (surface_force) then
    omega = 0. ! fixed position
 else
    omega = 1.
 endif

 a = 1.  ! semi-major axis
 x = a*cos(omega*ti)
 y = a*sin(omega*ti)

 !--positions of primary and secondary (exact)
 mtot = mass1 + mass2_temp
 x1 = -mass2_temp/mtot*x
 y1 = -mass2_temp/mtot*y

 x2 = mass1/mtot*x
 y2 = mass1/mtot*y

end subroutine update_binary

!----------------------------------------------
!+
!  compute the force on a given particle from
!  the binary system
!+
!----------------------------------------------
subroutine binary_force(xi,yi,zi,ti,fxi,fyi,fzi,phi)
 real, intent(in)  :: xi,yi,zi,ti
 real, intent(out) :: fxi,fyi,fzi,phi
 real :: dx1,dy1,dz1,rr1,f1
 real :: dx2,dy2,dz2,rr2,f2,r2
 real :: dr1,dr2,phi1,phi2

 !--compute gravitational force on gas particle i
 !  from the binary
 dx1 = xi - x1
 dy1 = yi - y1
 dz1 = zi

 dx2 = xi - x2
 dy2 = yi - y2
 dz2 = zi

 rr1 = dx1*dx1 + dy1*dy1 + dz1*dz1
 rr2 = dx2*dx2 + dy2*dy2 + dz2*dz2

 dr1 = 1./sqrt(rr1 + eps_soft1**2)
 dr2 = 1./sqrt(rr2 + eps_soft2**2)

 if (surface_force) then
    r2 = sqrt(rr2)
    if (r2 < 2.*eps_soft2) then
       !--add surface force to keep particles outside of r_planet
       f2 = mass2_temp/(rr2*r2)*(1.-((2.*eps_soft2-r2)/eps_soft2)**4)
    else
       !--1/r potential
       f2 = mass2_temp/(rr2*r2)
    endif
 else
    !--normal softened potential
    f2 = mass2_temp*dr2*dr2*dr2
 endif
 f1 = mass1*dr1*dr1*dr1

 fxi = -dx1*f1 - dx2*f2
 fyi = -dy1*f1 - dy2*f2
 fzi = -dz1*f1 - dz2*f2

 ! Note: phi1 is the Newtonian potential and does not include then surface force above;
 !       however, this is not critical since phi is only used for timestep control
 phi1 = -mass1*dr1
 phi2 = -mass2_temp*dr2
 phi  = phi1 + phi2

end subroutine binary_force

!----------------------------------------------
!+
!  return position and velocity of two stars
!  (used to write these to the dump file)
!+
!----------------------------------------------
subroutine binary_posvel(ti,posmh,vels)
 real, intent(in)  :: ti
 real, intent(out) :: posmh(10)
 real, intent(out) :: vels(6)
 real :: mtot,omega,vx,vy

 !--positions of primary and secondary (exact)
 posmh(1) = x1
 posmh(2) = y1
 posmh(3) = 0.
 posmh(4) = mass1
 posmh(5) = accradius1

 posmh(6)  = x2
 posmh(7)  = y2
 posmh(8)  = 0.
 posmh(9)  = mass2_temp
 posmh(10) = accradius2

 if (surface_force) then
    omega = 0. ! fixed position
 else
    omega = 1.
 endif

 vx = -omega*x2
 vy = omega*y2

 mtot = mass1 + mass2_temp
 vels(1) = -mass2_temp/mtot*vx
 vels(2) = -mass2_temp/mtot*vy
 vels(3) = 0.
 vels(4) = mass1/mtot*vx
 vels(5) = mass1/mtot*vy
 vels(6) = 0.

end subroutine binary_posvel

!----------------------------------------------
!+
!  flags whether or not a particle has been
!  accreted by either star
!+
!----------------------------------------------
logical function binary_accreted(xi,yi,zi,mi,ti)
 real, intent(in) :: xi,yi,zi,mi,ti
 real :: dx1,dy1,dz1,rr1
 real :: dx2,dy2,dz2,rr2

 dx1 = xi - x1
 dy1 = yi - y1
 dz1 = zi

 dx2 = xi - x2
 dy2 = yi - y2
 dz2 = zi

 rr1 = dx1*dx1 + dy1*dy1 + dz1*dz1
 rr2 = dx2*dx2 + dy2*dy2 + dz2*dz2

 if (rr1 < accradius1**2) then
    binary_accreted = .true.
    accretedmass1 = accretedmass1 + mi
 elseif (rr2 < accradius2**2) then
    binary_accreted = .true.
    accretedmass2 = accretedmass2 + mi
 else
    binary_accreted = .false.
 endif

end function binary_accreted

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_externbinary(iunit)
 use infile_utils,   only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(mass1,'mass1','m1 of central binary system (if iexternalforce=binary)',iunit)
 call write_inopt(mass2,'mass2','m2 of central binary system (if iexternalforce=binary)',iunit)
 call write_inopt(accradius1,'accradius1','accretion radius of primary',iunit)
 call write_inopt(accradius2,'accradius2','accretion radius of secondary (if iexternalforce=binary)',iunit)
 call write_inopt(eps_soft1,'eps_soft1','Plummer softening of primary',iunit)
 call write_inopt(eps_soft2,'eps_soft2','Plummer softening of secondary',iunit)
 call write_inopt(ramp,'ramp','ramp up mass of secondary over first 5 orbits?',iunit)

end subroutine write_options_externbinary

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_externbinary(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,error
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: where = 'read_options_externbinary'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('mass1')
    read(valstring,*,iostat=ierr) mass1
    ngot = ngot + 1
    if (mass1 < 0.) then
       call fatal(where,'invalid setting for m1 (<0)')
    endif
 case('mass2')
    read(valstring,*,iostat=ierr) mass2
    ngot = ngot + 1
    if (mass2 < epsilon(mass2)) &
       call fatal(where,'invalid setting for m2 (<0)')
    if (mass2/mass1 > 1.e10)  call error(where,'binary mass ratio is huge!!!')
 case('accradius1') ! cannot be compulsory, because also handled in parent routine
    read(valstring,*,iostat=ierr) accradius1
    if (accradius1 < 0.)  call fatal(where,'negative accretion radius')
 case('accradius2')
    read(valstring,*,iostat=ierr) accradius2
    ngot = ngot + 1
    if (accradius2 < 0.)  call fatal(where,'negative accretion radius')
 case('eps_soft1')
    read(valstring,*,iostat=ierr) eps_soft1
    if (eps_soft1 < 0.)  call fatal(where,'negative eps_soft1')
 case('eps_soft2')
    read(valstring,*,iostat=ierr) eps_soft2
    if (eps_soft2 < 0.)  call fatal(where,'negative eps_soft2')
 case('ramp')
    read(valstring,*,iostat=ierr) ramp
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 2)

end subroutine read_options_externbinary

!-----------------------------------------------------------------------
!+
!  writes relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_externbinary(hdr,time,ierr)
 use dump_utils, only:lentag,dump_h,add_to_rheader
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time
 integer,      intent(out)   :: ierr
 character(len=lentag)       :: tags(18)
 real    :: rheader(18)
 integer :: i

 ierr = 0
 call binary_posvel(time,rheader(1:10),rheader(11:16))

 rheader(17) = accretedmass1
 rheader(18) = accretedmass2

 tags(1:16) = (/'x1 ','y1 ','z1 ','m1 ','h1 ','x2 ','y2 ','z2 ','m2 ', &
                'h2 ','vx1','vy1','vz1','vx2','vy2','vz2'/)
 tags(17:18) = (/'accretedmass1','accretedmass2'/)

 do i=1,18
    call add_to_rheader(rheader(i),tags(i),hdr,ierr)
 enddo

end subroutine write_headeropts_externbinary

!-----------------------------------------------------------------------
!+
!  reads relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_externbinary(hdr,ierr)
 use dump_utils, only:dump_h,extract
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 integer :: ierr1,ierr2

 ierr  = 0
 call extract('accretedmass1',accretedmass1,hdr,ierr1)
 call extract('accretedmass2',accretedmass2,hdr,ierr2)

 if (ierr1 /= 0 .or. ierr2 /= 0) then
    write(*,*) ' ERROR extracting accretedmass1 and accretedmass2 from file'
    ierr = 1
 endif

end subroutine read_headeropts_externbinary

end module extern_binary
