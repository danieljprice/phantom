!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module gravwaveutils
!
! Compute gravitational wave emission using the quadrupole formula.
! Assumes a distance of 1 Mpc
!
! :References: Toscani et al. (2021) https://arxiv.org/abs/2111.05145
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - gw       : *calculate gravitational wave strain*
!   - phi_gw   : *angle betw. projection of los in xy plane and y axis (deg)*
!   - theta_gw : *rotation of xy plane (deg)*
!
! :Dependencies: infile_utils, io, mpiutils, physcon, units
!
 implicit none

 public :: calculate_strain,get_G_on_dc4,get_strain_from_circular_binary
 public :: write_rotated_strain_components
 public :: write_options_gravitationalwaves,read_options_gravitationalwaves

 logical, public :: calc_gravitwaves = .false.
 real, public  :: theta_gw = 0.  ! public as can be set by setup routine
 real, private :: phi_gw = 0.

 private

contains

!--------------------------------------------------------------------------------
!+
!  This subroutine computes the gravitational wave strain at a distance of 1Mpc
!  for an arbitrary collection of particles
!+
!--------------------------------------------------------------------------------
subroutine calculate_strain(hx,hp,pmass,ddq_xy,x0,v0,a0,npart,xyzh,vxyz,axyz,&
                            axyz1,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)
 use io,           only:master,error
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use physcon,      only:pi
 use mpiutils,     only:reduceall_mpi
 real, intent(out)             :: hx(4),hp(4),ddq_xy(3,3)
 real, intent(in)              :: xyzh(:,:), vxyz(:,:), axyz(:,:), pmass,x0(3),v0(3),a0(3)
 real, intent(inout), optional :: axyz1(:,:) !optional, only if there are external forces
 integer,intent(in),  optional :: nptmass
 real,   intent(in),  optional :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 integer, intent(in)           :: npart
 real                          :: ddq(6),x,y,z,vx,vy,vz,ax,ay,az,fac,r2
 real                          :: xp,yp,zp,vxp,vyp,vzp,axp,ayp,azp,mp
 integer                       :: i
 real                          :: eta,phi,sinphi,cosphi,sineta,coseta,cosphi2,sinphi2,&
                                  cos2phi,sin2phi,cos2eta,sin2eta,sineta2,coseta2
 real                          :: theta,lambda
 real,dimension(3,3)           :: R,Rt,quadrupole_2deriv,intermediate_result

 ! change this line if you want to start from a different value of phi and eta
 ! you need these angles for the angular distribution of the quadrupole radiation
 ! eta is the angle between the line of sight and the normal to the orbit
 ! phi is the angle between the projection of the los in the xy plane and the y axis
 ! more details and sketch in Toscani et al 2021

 eta = 0.
 phi = phi_gw*pi/180.
 ! define new functions instead of sin and cos
 sinphi = sin(phi)
 cosphi = cos(phi)
 sinphi2 = sinphi*sinphi
 cosphi2 = cosphi*cosphi
 sin2phi = sin(2*phi)
 cos2phi = cos(2*phi)
 !
 ! initialise moment of inertia and its second deriv to zero
 !
 ddq(:) = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyz,axyz,axyz1,pmass,x0,v0,a0) &
 !$omp private(i,x,y,z,vx,vy,vz,ax,ay,az,r2) &
 !$omp reduction(+:ddq)
 do i=1,npart
    if (xyzh(4,i) > tiny(xyzh)) then  !if not accreted
       x  = xyzh(1,i) - x0(1)
       y  = xyzh(2,i) - x0(2)
       z  = xyzh(3,i) - x0(3)
       vx = vxyz(1,i) - v0(1)
       vy = vxyz(2,i) - v0(2)
       vz = vxyz(3,i) - v0(3)
       ax = axyz(1,i) - a0(1)
       ay = axyz(2,i) - a0(2)
       az = axyz(3,i) - a0(3)
       if (present(axyz1)) then ! to avoid memory allocation in phantom
          ax = ax + axyz1(1,i)
          ay = ay + axyz1(2,i)
          az = az + axyz1(3,i)
       endif
       r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))

       ! calculate the components of the traceless quadrupole moment Q--not necessary but maybe useful
       !q(1) = q(1) + pmass*(x*x)!-onethird*r2) !qxx
       !q(2) = q(2) + pmass*(x*y)!-onethird*r2) !qxy
       !q(3) = q(3) + pmass*(x*z)!-onethird*r2) !qxz
       !q(4) = q(4) + pmass*(y*y)!-onethird*r2) !qyy
       !q(5) = q(5) + pmass*(y*z)!-onethird*r2) !qyz
       !q(6) = q(6) + pmass*(z*z)!-onethird*r2) !qzz

       ! second time derivative of the traceless quadrupole moment Q
       ddq(1) = ddq(1) + pmass*(2.*vx*vx+x*ax+x*ax) !ddqxx
       ddq(2) = ddq(2) + pmass*(2.*vx*vy+x*ay+y*ax) !ddqxy
       ddq(3) = ddq(3) + pmass*(2.*vx*vz+x*az+z*ax) !ddqxz
       ddq(4) = ddq(4) + pmass*(2.*vy*vy+y*ay+y*ay) !ddqyy
       ddq(5) = ddq(5) + pmass*(2.*vy*vz+y*az+z*ay) !ddqyz
       ddq(6) = ddq(6) + pmass*(2.*vz*vz+z*az+z*az) !ddqzz
    endif
 enddo
 !omp end parallel do

 ddq = reduceall_mpi('+',ddq)  ! reduce across MPI threads

 ! add contribution from sink particles
 ! note: sink particles are present on all MPI threads, so no reduceall call
 if (present(xyzmh_ptmass) .and. present(vxyz_ptmass) .and. present(nptmass) .and. present(fxyz_ptmass)) then
    do i=1,nptmass
       xp  = xyzmh_ptmass(1,i) - x0(1)
       yp  = xyzmh_ptmass(2,i) - x0(2)
       zp  = xyzmh_ptmass(3,i) - x0(3)
       mp  = xyzmh_ptmass(4,i)
       vxp = vxyz_ptmass(1,i) - v0(1)
       vyp = vxyz_ptmass(2,i) - v0(2)
       vzp = vxyz_ptmass(3,i) - v0(3)
       axp = fxyz_ptmass(1,i) - a0(1)
       ayp = fxyz_ptmass(2,i) - a0(2)
       azp = fxyz_ptmass(3,i) - a0(3)
       ddq(1) = ddq(1) + mp*(2.*vxp*vxp+xp*axp+xp*axp) !ddqxx
       ddq(2) = ddq(2) + mp*(2.*vxp*vyp+xp*ayp+yp*axp) !ddqxy
       ddq(3) = ddq(3) + mp*(2.*vxp*vzp+xp*azp+zp*axp) !ddqxz
       ddq(4) = ddq(4) + mp*(2.*vyp*vyp+yp*ayp+yp*ayp) !ddqyy
       ddq(5) = ddq(5) + mp*(2.*vyp*vzp+yp*azp+zp*ayp) !ddqyz
       ddq(6) = ddq(6) + mp*(2.*vzp*vzp+zp*azp+zp*azp) !ddqzz
    enddo
 endif

 ! define some parameters
 fac = get_G_on_dc4()

 ! rotate if the inclination angle is non-zero
 if (abs(theta_gw) > tiny(0.)) then

    theta = theta_gw*pi/180.  ! convert to radians
    lambda = theta            ! removed some old code here, hence the reassignment

    !'rotate back' the second time derivative of the mass quadrupole into the x-y plane
    ! rotation matrix
    R = transpose(reshape((/ cos(lambda), 0., sin(lambda), 0., 1., 0.,-sin(lambda), 0., cos(lambda) /), shape(R)))

    ! transpose of the rotation matrix
    Rt  = transpose(reshape((/ cos(lambda), 0., -sin(lambda), 0., 1., 0.,sin(theta), 0., cos(lambda) /), shape(Rt)))

    ! second time derivative of the moment of inertia matrix
    quadrupole_2deriv = transpose(reshape((/ ddq(1),ddq(2),ddq(3),ddq(2),ddq(4),ddq(5),ddq(3),ddq(5),ddq(6) /), &
                                  shape(quadrupole_2deriv)))
    intermediate_result = matmul(quadrupole_2deriv, R)
    ! second time derivative rotated back in the x-y plane
    ddq_xy = matmul(Rt, intermediate_result)

    ! derive the quadrupole radiation
    do i=1,4
       sineta=sin(eta)
       coseta=cos(eta)
       sineta2=sineta*sineta
       coseta2=coseta*coseta
       sin2eta=sin(2*eta)
       cos2eta=cos(2*eta)
       ! angular distribution of the quadrupole radiation
       hp(i) = fac*(ddq_xy(1,1)*(cosphi2 - sinphi2*coseta2) &
                  + ddq_xy(2,2)*(sinphi2 - cosphi2*coseta2) &
                  - ddq_xy(3,3)*sineta2-ddq_xy(1,2)*(sin2phi)*(1+coseta2) &
                  + ddq_xy(1,3)*(sinphi*sin2eta)+ddq_xy(2,3)*cosphi*sin2eta)
       hx(i) = 2.*fac*(0.5*(ddq_xy(1,1)-ddq_xy(2,2))*sin2phi*coseta &
                       + ddq_xy(1,2)*(cos2phi*coseta) &
                       - ddq_xy(1,3)*(cosphi*sineta) &
                       + ddq_xy(2,3)*sinphi*sineta)
       eta=eta+pi/6.
    enddo
 else
    ddq_xy = 0.
    !
    ! assume default values for the two angles otherwise
    !
    do i=1,4
       sineta=sin(eta)
       coseta=cos(eta)
       sineta2=sineta*sineta
       coseta2=coseta*coseta
       sin2eta=sin(2*eta)
       cos2eta=cos(2*eta)
       ! angular distribution of the quadrupole radiation
       hp(i) = fac*(ddq(1)*(cosphi2 - sinphi2*coseta2) &
                  + ddq(4)*(sinphi2 - cosphi2*coseta2) &
                  - ddq(6)*sineta2-ddq(2)*(sin2phi)*(1+coseta2) &
                  + ddq(3)*(sinphi*sin2eta)+ddq(5)*cosphi*sin2eta)
       hx(i) = 2.*fac*(0.5*(ddq(1)-ddq(4))*sin2phi*coseta &
                       + ddq(2)*(cos2phi*coseta) &
                       - ddq(3)*(cosphi*sineta) &
                       + ddq(5)*sinphi*sineta)
       eta=eta+pi/6.
    enddo
 endif

end subroutine calculate_strain

!--------------------------------------------------------------------------------
!+
!  Get prefactor for strain calculation (G/(c^4 * distance)) in code units
!  Involves some unit conversions in the code
!+
!--------------------------------------------------------------------------------
real function get_G_on_dc4() result(fac)
 use units,   only:umass,udist,utime
 use physcon, only:gg,c,Mpc
 real(kind=8) :: d,gc4

 d   = Mpc/udist                            ! 1Mpc in code units
 gc4 = (gg/c**4) / (utime**2/(umass*udist)) ! G/c^4 in code units
 fac = real(gc4/d)

end function get_G_on_dc4

!--------------------------------------------------------------------------------
!+
!   Exact solution for circular binaries (used for unit tests)
!+
!--------------------------------------------------------------------------------
subroutine get_strain_from_circular_binary(t,m1,m2,r,eta,hx,hp)
 real, intent(in)  :: t,m1,m2,r,eta
 real, intent(out) :: hx,hp
 real :: mred,term,omega2,omega,fac

 mred = m1*m2/(m1 + m2)   ! reduced mass
 omega2 = (m1 + m2)/r**3  ! Keplerian speed
 omega  = sqrt(omega2)

 ! following lines are written a bit differently to the expression
 ! in calculate_strain as a better check of the code

 fac = get_G_on_dc4()
 term = 4.*mred*r**2*omega2*fac

 ! Eqs 7-8 in Toscani et al. 2021, see e.g. Maggiore (2007)
 ! the minus sign is presumably because the binary goes
 ! anticlockwise by default ?
 hx = -term*cos(eta)*sin(2.*omega*t)
 hp = -term*0.5*(1. + cos(eta)**2)*cos(2.*omega*t)

end subroutine get_strain_from_circular_binary

!--------------------------------------------------------------------------------
!+
!   Exact solution for circular binaries (used for unit tests)
!+
!--------------------------------------------------------------------------------
subroutine write_rotated_strain_components(time,ddq_xy)
 real, intent(in) :: time,ddq_xy(3,3)
 logical, save    :: firstdump=.true.
 integer :: iuu

 ! Write a file where I append all the values of the strain wrt time
 if (firstdump) then
    firstdump = .false.
    open(newunit=iuu,file='quadrupole_plane_xy.txt',status='replace')
    write(iuu,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
          1, 'time',  &
          2, 'ddm11', &
          3, 'ddm12', &
          4, 'ddm13', &
          5, 'ddm22', &
          6, 'ddm23', &
          7, 'ddm33'
 else
    open(newunit=iuu,file='quadrupole_plane_xy.txt',position='append')
 endif
 write(iuu,'(7(es18.10,1X))') time, ddq_xy(1,1),ddq_xy(1,2),ddq_xy(1,3),&
                              ddq_xy(2,2),ddq_xy(2,3),ddq_xy(3,3)

end subroutine write_rotated_strain_components

!-----------------------------------------------------------------------
!+
!  writes options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_gravitationalwaves(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# gravitational waves'
 call write_inopt(calc_gravitwaves,'gw','calculate gravitational wave strain',iunit)
 if (calc_gravitwaves) then
    call write_inopt(theta_gw,'theta_gw','rotation of xy plane (deg)',iunit)
    call write_inopt(phi_gw,'phi_gw','angle betw. projection of los in xy plane and y axis (deg)',iunit)
 endif

end subroutine write_options_gravitationalwaves

!-----------------------------------------------------------------------
!+
!  reads options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_gravitationalwaves(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal,warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 character(len=*), parameter :: tag = 'gravwaves'
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('gw')
    read(valstring,*,iostat=ierr) calc_gravitwaves
    ngot = ngot + 1
 case('theta_gw')
    read(valstring,*,iostat=ierr) theta_gw
    ngot = ngot + 1
 case('phi_gw')
    read(valstring,*,iostat=ierr) phi_gw
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_gravitationalwaves

end module gravwaveutils
