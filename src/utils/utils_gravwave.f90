!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module gravwaveutils
!
! None
!
! :References: None
!
! :Owner: Martina Toscani
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, io, physcon, timestep, units
!
 implicit none

 public :: calculate_strain

 private

contains

! This subroutine computes the gravitational wave strain at a distance of 1Mpc
subroutine calculate_strain(hx,hp,pmass,ddq_xy,x0,v0,a0,npart,xyzh,vxyz,axyz,&
                            axyz1,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)
 use units,                    only:umass,udist,utime
 use physcon,                  only:gg,c,Mpc,pi
 use io,                       only:master,error
 use infile_utils,             only:open_db_from_file,inopts,read_inopt,close_db
 use timestep,                 only:time
 real, intent(out)             :: hx(4),hp(4),ddq_xy(3,3)
 real, intent(in)              :: xyzh(:,:), vxyz(:,:), axyz(:,:), pmass,x0(3),v0(3),a0(3)
 real, intent(inout), optional :: axyz1(:,:) !optional, only if there are external forces
 integer,intent(in),  optional :: nptmass
 real,   intent(in),  optional :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 integer, intent(in)           :: npart
 real                          :: q(6),ddq(6),x,y,z,vx,vy,vz,ax,ay,az,fac,r2,gc4,d
 real                          :: xp,yp,zp,vxp,vyp,vzp,axp,ayp,azp,mp
 real, parameter               :: onethird = 1./3.
 integer                       :: i
 real                          :: eta,phi,sinphi,cosphi,sineta,coseta,cosphi2,sinphi2,&
                                  cos2phi,sin2phi,cos2eta,sin2eta,sineta2,coseta2
 logical                       :: iexist
 character(len=120)            :: filename
 integer                       :: ierr,iuu
 integer, parameter            :: iunit = 21
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)
 real                          :: theta,ecc,lambda
 real,dimension(3,3)           :: R, Rt, quadrupole_2deriv,intermediate_result
 logical, save                 :: firstdump=.true.

!change this line if you want to start from a different value of phi and eta
!you need these angles for the angular distribution of the quadrupole radiation
!eta is the angle between the line of sight and the normal to the orbit
!phi is the angle between the projection of the los in the xy plane and the y axis
!more details and sketch in Toscani et al 2021
 phi=0.
 eta=0.

!initialise moment of inertia and its second deriv to zero
 q(:)   = 0.
 ddq(:) = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyz,axyz,axyz1,pmass,x0,v0,a0) &
 !$omp private(i,x,y,z,vx,vy,vz,ax,ay,az,r2) &
 !$omp reduction(+:q,ddq)
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
       q(1) = q(1) + pmass*(x*x)!-onethird*r2) !qxx
       q(2) = q(2) + pmass*(x*y)!-onethird*r2) !qxy
       q(3) = q(3) + pmass*(x*z)!-onethird*r2) !qxz
       q(4) = q(4) + pmass*(y*y)!-onethird*r2) !qyy
       q(5) = q(5) + pmass*(y*z)!-onethird*r2) !qyz
       q(6) = q(6) + pmass*(z*z)!-onethird*r2) !qzz

       ddq(1) = ddq(1) + pmass*(2.*vx*vx+x*ax+x*ax) !ddqxx
       ddq(2) = ddq(2) + pmass*(2.*vx*vy+x*ay+y*ax) !ddqxy
       ddq(3) = ddq(3) + pmass*(2.*vx*vz+x*az+z*ax) !ddqxz
       ddq(4) = ddq(4) + pmass*(2.*vy*vy+y*ay+y*ay) !ddqyy
       ddq(5) = ddq(5) + pmass*(2.*vy*vz+y*az+z*ay) !ddqyz
       ddq(6) = ddq(6) + pmass*(2.*vz*vz+z*az+z*az) !ddqzz
    endif
 enddo
 !omp end parallel do

 ! add contribution from sink particles
 if (present(xyzmh_ptmass) .and. present(vxyz_ptmass) .and. present(nptmass) .and. present(fxyz_ptmass)) then
      do i=1,nptmass
       xp  = xyzmh_ptmass(1,i) - x0(1)
       yp  = xyzmh_ptmass(2,i) - x0(2)
       zp  = xyzmh_ptmass(3,i) - x0(3)
       vxp = vxyz_ptmass(1,i) - v0(1)
       vyp = vxyz_ptmass(2,i) - v0(2)
       vzp = vxyz_ptmass(3,i) - v0(3)
       axp = fxyz_ptmass(1,i) - a0(1)
       ayp = fxyz_ptmass(2,i) - a0(2)
       azp = fxyz_ptmass(3,i) - a0(3)
       mp  = xyzmh_ptmass(4,i)
       ddq(1) = ddq(1) + mp*(2.*vxp*vxp+xp*axp+xp*axp) !ddqxx
       ddq(2) = ddq(2) + mp*(2.*vxp*vyp+xp*ayp+yp*axp) !ddqxy
       ddq(3) = ddq(3) + mp*(2.*vxp*vzp+xp*azp+zp*axp) !ddqxz
       ddq(4) = ddq(4) + mp*(2.*vyp*vyp+yp*ayp+yp*ayp) !ddqyy
       ddq(5) = ddq(5) + mp*(2.*vyp*vzp+yp*azp+zp*ayp) !ddqyz
       ddq(6) = ddq(6) + mp*(2.*vzp*vzp+zp*azp+zp*azp) !ddqzz
      enddo
 endif

 !define some parameters
 d   = Mpc/udist                            ! 1Mpc in code units
 gc4 = (gg/c**4) / (utime**2/(umass*udist)) ! G/c^4 in code units
 fac = (gc4/d)

 !Read inclination angle theta and eccentricity from the setup file
 filename = 'grtde'//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
  print "(a)",'reading setup options from '//trim(filename)
  nerr = 0
  ierr = 0
  call open_db_from_file(db,filename,iunit,ierr)
  call read_inopt(ecc,'ecc',db,min=0.,max=1.,errcount=nerr)
  call read_inopt(theta,'theta',db,errcount=nerr)
  call close_db(db)
  print*, "orbit inclination angle", theta
  theta = theta*pi/180       ! convert theta in radians
  print*, "cos(theta)", cos(theta) !just a check
  print*, "sin(theta)", sin(theta) !just a check
  print*, "eccentricity", ecc      !just a check

  !I do these change because for elleptic orbits david rotate the orbit wrt y axis
  !by -theta and not theta
   if (ecc==1.) then
    lambda=theta
   else
    lambda=-theta
   endif

  !'rotate back' in the x-y plane the second time derivative of the mass quadrupole
  !rotation matrix
  R = transpose(reshape((/ cos(lambda), 0., sin(lambda), 0., 1., 0.,-sin(lambda), 0., cos(lambda) /), shape(R)))
  print*, 'R elements', R(1,1),R(1,3),R(3,1),R(3,3) !just a check
  !transpose of the rotation matrix
  Rt  = transpose(reshape((/ cos(lambda), 0., -sin(lambda), 0., 1., 0.,sin(theta), 0., cos(lambda) /), shape(Rt)))
  print*, 'Rt elements', Rt(1,1),Rt(1,3),Rt(3,1),Rt(3,3) !just a check
  ! second time derivative of the mom inertia matrix
  quadrupole_2deriv =transpose(reshape((/ ddq(1),ddq(2),ddq(3),ddq(2),ddq(4),ddq(5),ddq(3),ddq(5),ddq(6) /), &
  & shape(quadrupole_2deriv)))
  intermediate_result = matmul(quadrupole_2deriv, R)
  !second time derivative rotated back in the x-y plane
  ddq_xy = matmul(Rt, intermediate_result)
  print*, ddq_xy

  ! Write a file where I append all the values of the strain wrt time
  if (firstdump) then
  firstdump = .false.
  open(newunit=iuu, file='quadrupole_plane_xy.txt',status='replace')
  write(iuu,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
  1, 'time',  &
  2, 'ddm11', &
  3, 'ddm12', &
  4, 'ddm13', &
  5, 'ddm22', &
  6, 'ddm23', &
  7, 'ddm33'
  else
  open(newunit=iuu, file='quadrupole_plane_xy.txt',position='append')
  endif
  write(iuu,'(7(es18.10,1X))') time, ddq_xy(1,1),ddq_xy(1,2),ddq_xy(1,3),&
                                ddq_xy(2,2),ddq_xy(2,3),ddq_xy(3,3)
  !maybe write an 'on fly' file with all the components of rotated second time derivative of mom inertia

  !derive the quadrupole radiation
  do i=1,4
  !define new functions instead of sin and cos
   sinphi=sin(phi)
   cosphi=cos(phi)
   sinphi2=sinphi*sinphi
   cosphi2=cosphi*cosphi
   sin2phi=sin(2*phi)
   cos2phi=cos(2*phi)
   sineta=sin(eta)
   coseta=cos(eta)
   sineta2=sineta*sineta
   coseta2=coseta*coseta
   sin2eta=sin(2*eta)
   cos2eta=cos(2*eta)

  !angular distribution of the quadrupole radiation
   hp(i) =fac*(ddq_xy(1,1)*(cosphi2 - sinphi2*coseta2)+ddq_xy(2,2)*(sinphi2 &
   -cosphi2*coseta2)-ddq_xy(3,3)*sineta2-ddq_xy(1,2)*(sin2phi)*(1+coseta2) &
   +ddq_xy(1,3)*(sinphi*sin2eta)+ddq_xy(2,3)*cosphi*sin2eta)
   hx(i)=2.*fac*(0.5*(ddq_xy(1,1)-ddq_xy(2,2))*sin2phi*coseta+ddq_xy(1,2)*(cos2phi*coseta)&
   -ddq_xy(1,3)*(cosphi*sineta)+ddq_xy(2,3)*sinphi*sineta)
   eta=eta+pi/6.
  enddo
 elseif (.not. iexist .or. ierr /= 0) then !derive the quadrupole radiation for setup different from grtde
   print*,'nosetup grtde'
    !maybe write an 'on fly'file with all the components of the second time derivative of the mom inertia
    !i am assuming that for no tde there is no need for all the rotating matrix stuff
    !if (firstdumpa) then
    !firstdumpa = .false.
    !open(unit=iu, file='quad_secondderiv_plane_xy_notde.txt',status='replace')
    !write(iu,"('#',13(1x,'[',i2.2,1x,a11,']',2x))") &
    !1, 'time',  &
    !2, 'q11', &
    !3, 'q12', &
    !4, 'q13', &
    !5, 'q22', &
    !6, 'q23', &
    !7, 'q33', &
    !8, 'ddm11', &
    !9, 'ddm12', &
    !10, 'ddm13', &
    !11, 'ddm22', &
    !12, 'ddm23', &
    !13, 'ddm33'
    !else
    !open(unit=iu, file='quad_secondderiv_plane_xy_notde.txt',position='append')
    !endif
    !write(iu,'(13(es18.10,1X))') time, q(1),q(2),q(3),q(4),q(5),q(6),&
  !                                ddq(1),ddq(2),ddq(3),ddq(4),ddq(5),ddq(6)

    do i=1,4
    !define new functions instead of sin and cos
     sinphi=sin(phi)
     cosphi=cos(phi)
     sinphi2=sinphi*sinphi
     cosphi2=cosphi*cosphi
     sin2phi=sin(2*phi)
     cos2phi=cos(2*phi)
     sineta=sin(eta)
     coseta=cos(eta)
     sineta2=sineta*sineta
     coseta2=coseta*coseta
     sin2eta=sin(2*eta)
     cos2eta=cos(2*eta)

    !angular distribution of the quadrupole radiation
      hp(i) =fac*(ddq(1)*(cosphi2 - sinphi2*coseta2)+ddq(4)*(sinphi2 &
      -cosphi2*coseta2)-ddq(6)*sineta2-ddq(2)*(sin2phi)*(1+coseta2) &
      +ddq(3)*(sinphi*sin2eta)+ddq(5)*cosphi*sin2eta)
      hx(i)=2.*fac*(0.5*(ddq(1)-ddq(4))*sin2phi*coseta+ddq(2)*(cos2phi*coseta)&
      -ddq(3)*(cosphi*sineta)+ddq(5)*sinphi*sineta)
      eta=eta+pi/6.
     enddo
 endif

end subroutine calculate_strain

end module gravwaveutils
