!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: extern_binary
!
!  DESCRIPTION:
!   This module contains routines relating to the computation
!   of a gravitational force/potential from a binary system
!   that decays due to gravitational waves
!   (2010 James Wetter & Daniel Price)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    accradius1 -- accretion radius of primary
!    accradius2 -- accretion radius of secondary (if iexternalforce=binary)
!
!  DEPENDENCIES: dump_utils, infile_utils, io
!+
!--------------------------------------------------------------------------
module extern_binary
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, public, parameter :: massp = 1.
 real, public :: massr = 0.001
 real, public :: a0 = 30.
 real, public :: accradius1  = 2.0
 real, public :: accradius2  = 2. ! accradius1 x mass ratio
 real, public :: direction = 1.
 real, public :: accretedmass1 = 0.
 real, public :: accretedmass2 = 0.

 public :: binary_force, binary_posvel, binary_accreted, update_binary
 public :: write_options_externbinary, read_options_externbinary
 public :: write_headeropts_externbinary, read_headeropts_externbinary
 private

 !real, private, save :: tset = -1.
 real, private :: xyzbin(3,2)

contains
!----------------------------------------------
!+
!  update position of binary at new time
!  only do this ONCE as sin/cos expensive
!+
!----------------------------------------------
subroutine update_binary(ti,surface_force)
 logical, optional, intent(in) :: surface_force
 real, intent(in) :: ti

 call compute_binary_pos(ti)

end subroutine update_binary

!----------------------------------------------
!+
!  Compute the position of the binary at
!   a particular time t (optimisation)
!+
!----------------------------------------------
subroutine compute_binary_pos(ti)
 real, intent(in) :: ti
 real             :: posmh(10)
 real             :: vels(6)

 call binary_posvel(ti,posmh,vels)
 xyzbin(1,1) = posmh(1)
 xyzbin(2,1) = posmh(2)
 xyzbin(3,1) = posmh(3)
 xyzbin(1,2) = posmh(6)
 xyzbin(2,2) = posmh(7)
 xyzbin(3,2) = posmh(8)
 !tset = ti

end subroutine compute_binary_pos

!----------------------------------------------
!+
!  compute the force on a given particle from
!  the binary system
!+
!----------------------------------------------
subroutine binary_force(xi,yi,zi,ti,fxi,fyi,fzi,phi,surface_force)
 logical, optional, intent(in) :: surface_force
 real, intent(in)  :: xi,yi,zi,ti
 real, intent(out) :: fxi,fyi,fzi,phi
 real :: x1,y1
 real :: dx1,dy1,dz1,dx2,dy2,dz2
 real :: rr1,rr2,f1,f2,dr1,dr2,phi1,phi2

 if (abs(massr) < tiny(massr)) then
    x1=0.
    y1=0.

    dx1 = xi - x1
    dy1 = yi - y1
    dz1 = zi

    rr1 = dx1*dx1 + dy1*dy1 + dz1*dz1
    dr1 = 1./sqrt(rr1)

    f1   = massp*dr1*dr1*dr1
    phi1 = massp*dr1

    fxi  = -dx1*f1
    fyi  = -dy1*f1
    fzi  = -dz1*f1
 else
    !if (abs(ti-tset) > epsilon(0.)) call compute_binary_pos(ti)

    !--compute gravitational force on gas particle i
    !  from the binary
    dx1 = xi - xyzbin(1,1)
    dy1 = yi - xyzbin(2,1)
    dz1 = zi - xyzbin(3,1)

    dx2 = xi - xyzbin(1,2)
    dy2 = yi - xyzbin(2,2)
    dz2 = zi - xyzbin(3,2)

    rr1 = dx1*dx1 + dy1*dy1 + dz1*dz1
    rr2 = dx2*dx2 + dy2*dy2 + dz2*dz2

    dr1 = 1./sqrt(rr1)
    dr2 = 1./sqrt(rr2)

    f1   = massp*dr1*dr1*dr1
    f2   = massp*massr*dr2*dr2*dr2

    fxi  = -dx1*f1 - dx2*f2
    fyi  = -dy1*f1 - dy2*f2
    fzi  = -dz1*f1 - dz2*f2

    phi1 = massp*dr1
    phi2 = massp*massr*dr2
    phi  = phi1 + phi2
 endif

 !write (*,*) ti,0.,0.,0.,0.,x2,y2,0.,0.,0.,0.,a,theta
 return
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
 real :: theta,a,mu,mtot,tau,x,y,omega,dadt,vx,vy

 if (abs(massr) < tiny(massr)) then
    posmh(1:3) = 0.
    posmh(4)   = massp
    posmh(5)   = accradius1
    posmh(6:8) = 0.
    posmh(9)   = massp*massr
    posmh(10)  = accradius2

    vels(1) = 0.
    vels(2) = 0.
    vels(3) = 0.
    vels(4) = 0.
    vels(5) = 0.
    vels(6) = 0.

 else

    mtot=massp+massp*massr;
    mu=(massp**2)*massr/mtot;

    !--time of coalescence
    tau=5./256.*a0**4/(mu * mtot**2)
    !--separation
    a=a0*((1.-ti/tau)**.25)

    !--relative angular position
    theta=direction*(-1./32.)*1./(mu*mtot**(2./3.))*(5./256.*1./(mu*mtot**(2./3.))*1./(tau-ti))**(-5./8.)

    !--relative position
    x=a*cos(theta)
    y=a*sin(theta)

    !--positions of primary and secondary (exact)
    posmh(1) = -massp*massr/mtot*x
    posmh(2) = -massp*massr/mtot*y
    posmh(3) = 0.
    posmh(4) = massp
    posmh(5) = accradius1

    posmh(6) = massp/mtot*x
    posmh(7) = massp/mtot*y
    posmh(8) = 0.
    posmh(9) = massp*massr
    posmh(10) = accradius2

    !-- rate of shrinking
    dadt = -64./5. * mu*mtot**2/a**3

    !-- rate of rotating
    omega = direction*(5./256. * 1./(mu*mtot**(2./3.)) *1./(tau-ti))**(3./8.)

    !--relative velocity
    vx = 1./a*dadt*x - omega*y
    vy = 1./a*dadt*y + omega*x

    !--absolute velocity
    vels(1) = -massp*massr/mtot*vx
    vels(2) = -massp*massr/mtot*vy
    vels(3) = 0.
    vels(4) = massp/mtot*vx
    vels(5) = massp/mtot*vy
    vels(6) = 0.

 endif

 return
end subroutine binary_posvel

!----------------------------------------------
!+
!  flags whether or not a particle has been
!  accreted by either star
!+
!----------------------------------------------
logical function binary_accreted(xi,yi,zi,mi,ti)
 real, intent(in) :: xi,yi,zi,mi,ti
 real :: x1,y1
 real :: dx1,dy1,dz1,dx2,dy2,dz2
 real :: rr1,rr2

 if (abs(massr) < tiny(massr)) then
    x1=0.
    y1=0.

    dx1 = xi - x1
    dy1 = yi - y1
    dz1 = zi

    rr1 = dx1*dx1 + dy1*dy1 + dz1*dz1
    rr2 = 1.e10
 else
    !if (abs(ti-tset) > epsilon(0.)) call compute_binary_pos(ti)
    dx1 = xi - xyzbin(1,1)
    dy1 = yi - xyzbin(2,1)
    dz1 = zi - xyzbin(3,1)

    dx2 = xi - xyzbin(1,2)
    dy2 = yi - xyzbin(2,2)
    dz2 = zi - xyzbin(3,2)

    rr1 = dx1*dx1 + dy1*dy1 + dz1*dz1
    rr2 = dx2*dx2 + dy2*dy2 + dz2*dz2
 endif

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
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(accradius1,'accradius1','accretion radius of primary',iunit)
 call write_inopt(accradius2,'accradius2','accretion radius of secondary (if iexternalforce=binary)',iunit)

end subroutine write_options_externbinary

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_externbinary(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,warn
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: where = 'read_options_externbinary_gw'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('accradius1') ! cannot be compulsory, because also handled in parent routine
    read(valstring,*,iostat=ierr) accradius1
    if (accradius1 < 0.)  call fatal(where,'negative accretion radius')
 case('accradius2')
    read(valstring,*,iostat=ierr) accradius2
    ngot = ngot + 1
    if (accradius2 < 0.)  call fatal(where,'negative accretion radius')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)

end subroutine read_options_externbinary

!-----------------------------------------------------------------------
!+
!  writes relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_externbinary(hdr,time,ierr)
 use dump_utils, only:dump_h,add_to_rheader,lentag
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time
 integer,      intent(out)   :: ierr
 character(len=lentag)       :: tags(20)
 real    :: rheader(20)
 integer :: i

 call binary_posvel(time,rheader(1:10),rheader(11:16))
 rheader(17) = a0
 rheader(18) = direction
 rheader(19) = accretedmass1
 rheader(20) = accretedmass2

 tags(1:16) = (/'x1 ','y1 ','z1 ','m1 ','h1 ','x2 ','y2 ','z2 ','m2 ', &
                'h2 ','vx1','vy1','vz1','vx2','vy2','vz2'/)
 tags(17:20) = (/'a0           ','direction    ',&
                 'accretedmass1','accretedmass2'/)
 do i=1,20
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
 use io,         only:error
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 real    :: time,mtot,mu,tau,m1,m2,x1,y1,z1,x2,y2,z2
 integer :: ierrs(7)

 ierr  = 0
 call extract('time',time,hdr,ierrs(1))
 call extract('m1',m1,hdr,ierrs(2))
 if (abs(m1 - massp) > epsilon(massp)) then
    call error('read_headeropts_externbinary','m1 in dump header not equal to mass of primary',var='m1',val=m1)
    ierr = 2
 endif
 call extract('m2',m2,hdr,ierrs(3))
 call extract('a0',a0,hdr,ierrs(4))
 call extract('direction',direction,hdr,ierrs(5))
 call extract('accretedmass1',accretedmass1,hdr,ierrs(6))
 call extract('accretedmass2',accretedmass2,hdr,ierrs(7))

 if (any(ierrs /= 0)) then
    write(*,*) ' ERROR: not enough parameters in dump header for external binary with gravitational waves'
    ierr = 1
    return
 endif

 ! these are for information only
 call extract('x1',x1,hdr,ierrs(1))
 call extract('y1',y1,hdr,ierrs(2))
 call extract('z1',z1,hdr,ierrs(3))
 call extract('x2',x2,hdr,ierrs(4))
 call extract('y2',y2,hdr,ierrs(5))
 call extract('z2',z2,hdr,ierrs(6))

 massr = m2/m1
 print*,'------ Inspiralling gravitational wave binary ------'
 print*,' Position of primary   : ',x1,y1,z1
 print*,' Position of secondary : ',x2,y2,z2
 print*,' Binary mass ratio     : ',massr
 print*,' Initial separation a0 : ',a0
 print*,' Accreted mass on prim : ',accretedmass1
 print*,' Accreted mass on sec. : ',accretedmass2

 mtot=massp+massp*massr;
 mu=massp**2 * massr/mtot;
 tau=5./256.*a0**4/(mu * mtot**2)
 print*,' Time of coalescence   : ',tau
 print*,' Remaining time        : ',tau - time
 print "(2x,a,2pf6.2,'%')",'Fraction complete     : ',time/tau
 if (abs(direction-1.) < tiny(direction)) then
    print*,' PROGRADE orbit'
 elseif (abs(direction+1.) < tiny(direction)) then
    print*,' RETROGRADE orbit'
 else
    write(*,*) 'ERROR: direction of orbit should be -1 or +1, but got direction = ',direction
 endif
 print "(1x,52('-'))"

end subroutine read_headeropts_externbinary

end module extern_binary
