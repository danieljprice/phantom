!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    beta  -- penetration factor
!    mh    -- mass of black hole (code units)
!    ms    -- mass of star       (code units)
!    phi   -- stellar rotation with respect to y-axis (in degrees)
!    r0    -- starting distance
!    rs    -- radius of star     (code units)
!    theta -- stellar rotation with respect to x-axis (in degrees)
!
!  DEPENDENCIES: centreofmass, externalforces, infile_utils, io, options,
!    physcon, prompting
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

 real :: beta,   &  ! penetration factor
         mh,     &  ! BH mass
         ms,     &  ! stellar mass
         rs,     &  ! stellar radius
         theta,  &  ! stellar tilting along x
         phi,    &  ! stellar tilting along y
         r0,     &  ! starting distance
         ecc        ! eccentricity

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use centreofmass
 use externalforces, only:mass1
 use externalforces, only:accradius1
 use options,        only:iexternalforce,damp
 use prompting,      only:prompt
 use physcon,        only:pi
 integer,  intent(inout) :: npart
 integer,  intent(inout) :: npartoftype(:)
 real,     intent(inout) :: massoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=120)      :: filename
 integer                 :: i,ierr
 logical                 :: iexist
 real                    :: Lx,Ly,Lz,L,Lp,Ltot(3)
 real                    :: rp,rt
 real                    :: x,y,z,vx,vy,vz
 real                    :: x0,y0,vx0,vy0,alpha

!
!-- Default runtime parameters
!
!
 beta  = 1.     ! penetration factor
 Mh    = 1.e6   ! BH mass
 Ms    = 1.     ! stellar mass
 rs    = 1.     ! stellar radius
 theta = 0.     ! stellar tilting along x
 phi   = 0.     ! stellar tilting along y
 ecc   = 1.                     ! eccentricity

 rt = (Mh/Ms)**(1./3.) * rs         ! tidal radius
 rp = rt/beta                       ! pericenter distance
 r0 = 4.9*rt                        ! starting radius

 !
 !-- Read runtime parameters from tdeparams file
 !
 filename = 'tde'//'.tdeparams'                                ! moddump should really know about the output file prefix...
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
   call write_setupfile(filename)
   print*,' Edit '//trim(filename)//' and rerun phantommoddump'
   stop
 endif
 rt = (Mh/Ms)**(1./3.) * rs         ! tidal radius
 rp = rt/beta                       ! pericenter distance

 !--Reset center of mass
 call reset_centreofmass(npart,xyzh,vxyzu)

 phi   = 0.
 theta = 0.

 call get_angmom(ltot,npart,xyzh,vxyzu)
 Lx = ltot(1)
 Ly = ltot(2)
 Lz = ltot(3)
 Lp = sqrt(Lx**2.0+Lz**2.0)
 if (Lx > 0.) then
    phi=acos(Lz/Lp)
 elseif (Lx < 0.) then
    phi=-acos(Lz/Lp)
 endif

!
!--Rotate the star so the momentum lies in the yz plan
!
 print*,'tilting along y axis: ',(phi*180/pi),'degrees'
 do i=1,npart
    x=xyzh(1,i)
    z=xyzh(3,i)
    xyzh(1,i)=x*cos(-phi)+z*sin(-phi)
    xyzh(3,i)=-x*sin(-phi)+z*cos(-phi)
    vx=vxyzu(1,i)
    vz=vxyzu(3,i)
    vxyzu(1,i)=vx*cos(-phi)+vz*sin(-phi)
    vxyzu(3,i)=-vx*sin(-phi)+vz*cos(-phi)
 enddo

!
!--Recheck the stellar angular momentum
!
 call get_angmom(ltot,npart,xyzh,vxyzu)
 lx = ltot(1)
 ly = ltot(2)
 lz = ltot(3)
 L  = sqrt(Lx**2.0+Ly**2.0+Lz**2.0)
 if (Ly < 0.) then
    theta=acos(Lz/L)
 elseif (Ly > 0.) then
    theta=-acos(Lz/L)
 endif

!
!--Rotate the star so the momentum lies along the z axis
!
 print*, 'tilting along x axis: ',(theta*180/pi),'degrees'
 do i=1,npart
    y=xyzh(2,i)
    z=xyzh(3,i)
    xyzh(2,i)=y*cos(-theta)-z*sin(-theta)
    xyzh(3,i)=y*sin(-theta)+z*cos(-theta)
    vy=vxyzu(2,i)
    vz=vxyzu(3,i)
    vxyzu(2,i)=vy*cos(-theta)-vz*sin(-theta)
    vxyzu(3,i)=vy*sin(-theta)+vz*cos(-theta)
 enddo

!
!--Recheck the stellar angular momentum
!

 call get_angmom(ltot,npart,xyzh,vxyzu)
 print*,'Stellar spin should now be along the z axis.'

 alpha = acos((rt*(1.+ecc)/(r0*beta)-1.)/ecc)         ! starting angle anti-clockwise from positive x-axis
 x0    = r0*cos(alpha)
 y0    = r0*sin(alpha)
 vx0   = sqrt(mh*beta/((1.+ecc)*rt)) * sin(alpha)
 vy0   = -sqrt(mh*beta/((1.+ecc)*rt)) * (cos(alpha)+ecc)

 !--Set input file parameters
 mass1          = Mh
 iexternalforce = 1
 damp           = 0.
 accradius1     = (2*Mh*rs)/((6.8565e2)**2) ! R_sch = 2*G*Mh*rs/c**2

 !--Tilting the star
 theta=theta*pi/180.0
 phi=phi*pi/180.0

 if (theta  /=  0.) then
    do i=1,npart
       y=xyzh(2,i)
       z=xyzh(3,i)
       xyzh(2,i)=y*cos(theta)-z*sin(theta)
       xyzh(3,i)=y*sin(theta)+z*cos(theta)
       vy=vxyzu(2,i)
       vz=vxyzu(3,i)
       vxyzu(2,i)=vy*cos(theta)-vz*sin(theta)
       vxyzu(3,i)=vy*sin(theta)+vz*cos(theta)
    enddo
 endif
 if (phi  /=  0.) then
    do i=1,npart
       x=xyzh(1,i)
       z=xyzh(3,i)
       xyzh(1,i)=x*cos(phi)+z*sin(phi)
       xyzh(3,i)=-x*sin(phi)+z*cos(phi)
       vx=vxyzu(1,i)
       vz=vxyzu(3,i)
       vxyzu(1,i)=vx*cos(phi)+vz*sin(phi)
       vxyzu(3,i)=-vx*sin(phi)+vz*cos(phi)
    enddo
 endif

 !--Putting star into orbit
 do i = 1, npart
    xyzh(1,i)  = xyzh(1,i)  + x0
    xyzh(2,i)  = xyzh(2,i)  + y0
    vxyzu(1,i) = vxyzu(1,i) + vx0
    vxyzu(2,i) = vxyzu(2,i) + vy0
 enddo

 theta = theta*pi/180.
 phi   = phi*pi/180.

 write(*,'(a)') "======================================================================"
 write(*,'(a,Es12.5,a)') ' Pericenter distance = ',rp,' R_sun'
 write(*,'(a,Es12.5,a)') ' Tidal radius        = ',rt,' R_sun'
 write(*,'(a,Es12.5,a)') ' Radius of star      = ',rs,' R_sun'
 write(*,'(a,Es12.5,a)') ' Starting distance   = ',r0,' R_sun'
 write(*,'(a,Es12.5,a)') ' Stellar mass        = ',Ms,' M_sun'
 write(*,'(a,Es12.5,a)') ' Tilting along x     = ',theta,' degrees'
 write(*,'(a,Es12.5,a)') ' Tilting along y     = ',phi,' degrees'

 write(*,'(a)') "======================================================================"

 return
end subroutine modify_dump

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# parameters file for a TDE phantommodump'
 call write_inopt(beta,  'beta',  'penetration factor',                                  iunit)
 call write_inopt(mh,    'mh',    'mass of black hole (code units)',                     iunit)
 call write_inopt(ms,    'ms',    'mass of star       (code units)',                     iunit)
 call write_inopt(rs,    'rs',    'radius of star     (code units)',                     iunit)
 call write_inopt(theta, 'theta', 'stellar rotation with respect to x-axis (in degrees)',iunit)
 call write_inopt(phi,   'phi',   'stellar rotation with respect to y-axis (in degrees)',iunit)
 call write_inopt(r0,    'r0',    'starting distance',                                   iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(beta,  'beta',  db,min=0.,errcount=nerr)
 call read_inopt(mh,    'mh',    db,min=0.,errcount=nerr)
 call read_inopt(ms,    'ms',    db,min=0.,errcount=nerr)
 call read_inopt(rs,    'rs',    db,min=0.,errcount=nerr)
 call read_inopt(theta, 'theta', db,min=0.,errcount=nerr)
 call read_inopt(phi,   'phi',   db,min=0.,errcount=nerr)
 call read_inopt(r0,    'r0',    db,min=0.,errcount=nerr)

 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

subroutine get_angmom(ltot,npart,xyzh,vxyzu)
 real, intent(out)   :: ltot(3)
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:), vxyzu(:,:)
 integer :: i
 real    :: L

 ltot = 0.
 do i=1,npart
     ltot(1) = ltot(1)+xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i)
     ltot(2) = ltot(2)+xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i)
     ltot(3) = ltot(3)+xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i)
 enddo

 L = sqrt(dot_product(ltot,ltot))

 print*,''
 print*,'Checking angular momentum orientation and magnitude...'
 print*,'Angular momentum is L = (',ltot(1),ltot(2),ltot(3),')'
 print*,'Angular momentum modulus is |L| = ',L
 print*,''

end subroutine get_angmom

end module moddump
