!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! None
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters:
!   - Mh2                  : *mass of second black hole (code units)*
!   - a                    : *spin of SMBH*
!   - beta                 : *penetration factor*
!   - ecc                  : *eccentricity of stellar orbit (1 for parabolic)*
!   - ecc_binary           : *eccenticity of black hole binary (1 for parabolic)*
!   - incline              : *inclination (in x-z plane)*
!   - iorigin              : *0 = COM of BBH, 1 = Sink 1, 2 = Sink 2*
!   - mh                   : *mass of black hole (code units)*
!   - ms                   : *mass of star       (code units)*
!   - r0                   : *starting distance  (code units)*
!   - rs                   : *radius of star     (code units)*
!   - semimajoraxis_binary : *sepration between black hole binary(code units)*
!   - theta                : *stellar rotation with respect to y-axis (in degrees)*
!
! :Dependencies: centreofmass, dim, externalforces, infile_utils, io,
!   metric, options, orbits, part, physcon, prompting, setbinary, units,
!   vectorutils
!
 implicit none

 real :: beta,    &  ! penetration factor
         Mh1,     &  ! BH mass1
         ms,      &  ! stellar mass
         rs,      &  ! stellar radius
         theta,   &  ! stellar tilting along x
         phi,     &  ! stellar tilting along y
         r0,      &  ! starting distance
         ecc,     &  ! eccentricity
         incline, &  ! inclination (about y axis)
         spin,    &  ! spin of black hole
         Mh2,     &  ! BH mass2
         semimajoraxis_binary, & !sepration
         ecc_binary !eccentricity of the black hole

 integer, public :: iorigin  ! which black hole to use for the origin
 logical,public :: use_binary,use_sink

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use centreofmass
 use externalforces, only:mass1
 use externalforces, only:accradius1,accradius1_hard
 use options,        only:iexternalforce,damp
 use dim,            only:gr
 use prompting,      only:prompt
 use physcon,        only:pi,solarm,solarr
 use units,          only:umass,udist,get_c_code
 use metric,         only:a
 use orbits,         only:isco_kerr
 use vectorutils,    only:rotatevec
 use setbinary,      only:set_binary
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use io,             only:fatal
 integer,  intent(inout) :: npart
 integer,  intent(inout) :: npartoftype(:)
 real,     intent(inout) :: massoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=120)      :: filename
 integer                 :: i,ierr
 logical                 :: iexist
 real                    :: Ltot(3)
 real                    :: rp,rt
 real                    :: x0,y0,vx0,vy0,vz0,alpha,z0
 real                    :: x,z,vx,vz
 real                    :: c_light,m0
 real                    :: accradius2
 real                    :: xyzstar(3),vxyzstar(3)
 real                    :: semia,period,hacc1,hacc2
!
!-- Default runtime parameters
!
!
 c_light        = get_c_code()
 beta  = 1.                  ! penetration factor
 Mh1   = 1.e6*solarm/umass   ! BH mass
 ms    = 1.  *solarm/umass   ! stellar mass
 rs    = 1.  *solarr/udist   ! stellar radius
 theta = 0.                  ! stellar tilting along x
 phi   = 0.                  ! stellar tilting along y
 ecc   = 1.                  ! eccentricity
 incline = 0.                ! inclination (in x-z plane)
 semimajoraxis_binary = 1000.*solarr/udist   !separation distance
 if (.not. gr) then
    spin = 0.
 else
    spin = 1. !upper limit on Sagitarrius A*'s spin is 0.1 (Fragione and Loeb 2020)'
 endif
 Mh2        = 1.e6*solarm/umass !setting mass of a second BH
 ecc_binary = 1.                !Eccentricity of the binary system
 m0 = Mh1
 rt = (m0/ms)**(1./3.) * rs

 ! setting a default r0 value
 r0 = 10*rt

 ! default parameters for binary (overwritten from .tdeparams file)
 use_binary = .false.
 use_sink = .false.
 iorigin = 0

 filename = 'tde'//'.tdeparams'                                ! moddump should really know about the output file prefix...
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    call write_setupfile(filename)
    print*,' Edit '//trim(filename)//' and rerun phantommoddump'
    stop
 endif
 print*,"--------------------------------------------"
 print*,use_binary,"use_binary"
 print*,"--------------------------------------------"

 m0 = Mh1
 if (use_binary) then
    select case(iorigin)
    case(1)
       m0 = Mh1
    case(2)
       m0 = Mh2
    case(0)
       m0 = Mh1 + Mh2
    end select
 endif

 rt = (m0/ms)**(1./3.) * rs
 rp = rt/beta
 theta=theta*pi/180.0
 !--Reset center of mass
 call reset_centreofmass(npart,xyzh,vxyzu)
 call get_angmom(ltot,npart,xyzh,vxyzu)
 if (ecc < 1.) then
    print*, 'Elliptical orbit'

    alpha = acos((rt*(1.+ecc)/(r0*beta)-1.)/ecc)     ! starting angle anti-clockwise from positive x-axis

    print*,rt*(1.+ecc),"(rt*(1.+ecc)",(r0*beta),"(r0*beta)",(rt*(1.+ecc)/(r0*beta)-1.),"(rt*(1.+ecc)/(r0*beta)-1.)"
    print*,(rt*(1.+ecc)/(r0*beta)-1.)/ecc,"(rt*(1.+ecc)/(r0*beta)-1.)/ecc"

    semia    = rp/(1.-ecc)
    period   = 2.*pi*sqrt(semia**3/mass1)
    hacc1    = rs/1.e8    ! Something small so set_binary doesn't warn about Roche lobe
    hacc2    = hacc1
    call set_binary(m0,ms,semia,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                    posang_ascnode=0.,arg_peri=90.,incl=0.,f=-180.)
    vxyzstar = vxyz_ptmass(1:3,2)
    xyzstar  = xyzmh_ptmass(1:3,2)
    nptmass  = 0
    x0 = xyzstar(1)
    y0 = xyzstar(2)
    z0 = xyzstar(3)
    vx0 = vxyzstar(1)
    vy0 = vxyzstar(2)
    vz0 = vxyzstar(3)

 elseif (abs(ecc-1.) < tiny(1.)) then
    print*, 'Parabolic orbit',r0,"r0"
    y0    = -2.*rp + r0
    x0    = -sqrt(r0**2 - y0**2)
    z0    = 0.
    vx0   = sqrt(2*m0/r0) * 2*rp / sqrt(4*rp**2 + x0**2)
    vy0   = sqrt(2*m0/r0) * x0   / sqrt(4*rp**2 + x0**2)
    vz0   = 0.
    xyzstar = (/x0,y0,z0/)
    vxyzstar = (/vx0,vy0,vz0/)
 else
    call fatal('moddump_tidal',' Hyperbolic orbits not implemented')
    x0 = 0.; y0 = 0.; z0 = 0.; vx0 = 0.; vy0 = 0.; vz0 = 0. ! avoid compiler warning
 endif

 !--Set input file parameters
 print*,x0,y0,"x0,y0",vx0,"vx0",vy0,"vy0"
 accradius1     = (2*Mh1)/(c_light**2) ! R_sch = 2*G*Mh/c**2
 print*,"use sink", use_sink
 !--Set input file parameters
 if (gr) then
    ! single black hole in GR
    mass1          = Mh1
    a              = spin
    call isco_kerr(a,mass1,accradius1)
    accradius1_hard = accradius1
 elseif (use_binary) then
    ! binary black hole in Newtonian gravity
    accradius2     = (2*Mh2)/(c_light**2) ! R_sch = 2*G*Mh/c**2
    call set_binary(Mh1,Mh2,semimajoraxis_binary,ecc_binary, &
                    accradius1,accradius2, &
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ierr)
 elseif (use_sink) then
    ! single black hole in Newtonian gravity using a sink particle
    nptmass = nptmass + 1
    xyzmh_ptmass(:,nptmass) = 0.
    xyzmh_ptmass(4,nptmass) = m0
    xyzmh_ptmass(ihacc,nptmass) = accradius1
    xyzmh_ptmass(ihsoft,nptmass) = accradius1
    vxyz_ptmass(:,nptmass) = 0.
 else
    ! single black hole in Newtonian gravity
    mass1          = m0
    iexternalforce = 1
    damp           = 0.
 endif

 if (theta /= 0.) then
    !--Tilting the star around y axis, i.e., in xz place with angle theta
    call rotatevec(xyzstar,(/0.,1.,0./),theta)
    call rotatevec(vxyzstar,(/0.,1.,0./),theta)
    x0 = xyzstar(1)
    y0 = xyzstar(2)
    z0 = xyzstar(3)
    vx0 = vxyzstar(1)
    vy0 = vxyzstar(2)
    vz0 = vxyzstar(3)
 endif
 !--shift origin to correct black hole
 if (use_binary) then
    if (iorigin==2) then
       x0 = x0 + xyzmh_ptmass(1,2)
       y0 = y0 + xyzmh_ptmass(2,2)
       z0 = z0 + xyzmh_ptmass(3,2)
       vx0 = vx0 + vxyz_ptmass(1,2)
       vy0 = vy0 + vxyz_ptmass(2,2)
       vz0 = vz0 + vxyz_ptmass(3,2)
    elseif (iorigin==1) then
       x0 = x0 + xyzmh_ptmass(1,1)
       y0 = y0 + xyzmh_ptmass(2,1)
       z0 = z0 + xyzmh_ptmass(3,1)
       vx0 = vx0 + vxyz_ptmass(1,1)
       vy0 = vy0 + vxyz_ptmass(2,1)
       vz0 = vz0 + vxyz_ptmass(3,1)
    endif
 endif

 !--Tilting the star
 !--Putting star into orbit
 do i = 1, npart
    xyzh(1,i)  = xyzh(1,i)  + x0
    xyzh(2,i)  = xyzh(2,i)  + y0
    xyzh(3,i)  = xyzh(3,i)  + z0
    vxyzu(1,i) = vxyzu(1,i) + vx0
    vxyzu(2,i) = vxyzu(2,i) + vy0
    vxyzu(3,i) = vxyzu(3,i) + vz0
 enddo

 !check angular momentum after putting star on orbit
 call get_angmom(ltot,npart,xyzh,vxyzu)

 !--Incline the orbit
 incline = incline*pi/180
 do i = 1, npart
    x=xyzh(1,i)
    z=xyzh(3,i)
    xyzh(1,i)= x*cos(incline) - z*sin(incline)
    xyzh(3,i)= x*sin(incline) + z*cos(incline)
    vx=vxyzu(1,i)
    vz=vxyzu(3,i)
    vxyzu(1,i)= vx*cos(incline) - vz*sin(incline)
    vxyzu(3,i)= vx*sin(incline) + vz*cos(incline)
 enddo

 theta=theta*180.0/pi
 write(*,'(a)') "======================================================================"
 write(*,'(a,Es12.5,a)') ' Pericenter distance = ',rp,' code units'
 write(*,'(a,Es12.5,a)') ' Tidal radius        = ',rt,' code units'
 write(*,'(a,Es12.5,a)') ' Radius of star      = ',rs,' code units'
 write(*,'(a,Es12.5,a)') ' Starting distance   = ',r0,' code units'
 write(*,'(a,Es12.5,a)') ' Stellar mass        = ',ms,' code units'
 write(*,'(a,Es12.5,a)') ' Tilting along y     = ',theta,' degrees'
 write(*,'(a,Es12.5,a)') ' Eccentricity of stellar orbit      = ',ecc
 write(*,'(a,Es12.5,a)') ' Mass of BH =',m0,' code units'
 write(*,'(a,Es12.5,a)') ' Inclination         = ',incline,' degrees'
 if (gr) then
    write(*,'(a,Es12.5,a)') ' Spin of black hole "a"       = ',a
 endif
 if (use_binary .and. .not.gr) then
    write(*,'(a)') "Use Binary BH = True"
    write(*,'(a,Es12.5,a)')'Eccentricity of the binary black hole system = ', ecc_binary
    write(*,'(a,Es12.5,a)') 'Separation betweent the two black holes = ', semimajoraxis_binary
 elseif (use_sink .and. .not.gr) then
    write(*,'(a)') "Use Sink particle as BH = True"
 else
    write(*,'(a)') "Use Binary BH = False"
 endif

 write(*,'(a)') "======================================================================"

end subroutine modify_dump

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use dim,          only:gr
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# parameters file for a TDE phantommodump'
 call write_inopt(beta,  'beta',  'penetration factor',                                  iunit)
 call write_inopt(Mh1,    'mh',    'mass of black hole (code units)',                    iunit)
 call write_inopt(ms,    'ms',    'mass of star       (code units)',                     iunit)
 call write_inopt(rs,    'rs',    'radius of star     (code units)',                     iunit)
 call write_inopt(theta, 'theta', 'stellar rotation with respect to y-axis (in degrees)',iunit)
 call write_inopt(r0,    'r0',    'starting distance  (code units)',                     iunit)
 call write_inopt(ecc,   'ecc',   'eccentricity of stellar orbit (1 for parabolic)',                      iunit)
 call write_inopt(incline,'incline','inclination (in x-z plane)',                          iunit)
 if (gr) then
    call write_inopt(spin,   'a',   'spin of SMBH',                                       iunit)
 endif
 if (.not. gr) then
    call write_inopt(use_binary, 'use binary', 'true/false', iunit)
    call write_inopt(use_sink, 'use sink', 'true/false', iunit)
    if (use_binary) then
       call write_inopt(iorigin,'iorigin','0 = COM of BBH, 1 = Sink 1, 2 = Sink 2', iunit)
       call write_inopt(Mh2,    'Mh2',    'mass of second black hole (code units)',iunit)
       call write_inopt(ecc_binary,    'ecc_binary',    'eccenticity of black hole binary (1 for parabolic)',iunit)
       call write_inopt(semimajoraxis_binary,'semimajoraxis_binary', 'sepration between black hole binary(code units)',iunit)
    endif
 endif
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 use dim,          only:gr
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(beta,   'beta',   db,min=0.,errcount=nerr)
 call read_inopt(Mh1,    'mh',     db,min=0.,errcount=nerr)
 call read_inopt(ms,     'ms',     db,min=0.,errcount=nerr)
 call read_inopt(rs,     'rs',     db,min=0.,errcount=nerr)
 call read_inopt(theta,  'theta',  db,min=0.,errcount=nerr)
 call read_inopt(r0,     'r0',     db,min=0.,errcount=nerr)
 call read_inopt(ecc,    'ecc',    db,min=0.,max=1.,errcount=nerr)
 call read_inopt(incline,'incline',db,       errcount=nerr)

 if (gr) then
    call read_inopt(spin, 'a',    db,min=-1.,max=1.,errcount=nerr)
 endif
 if (.not. gr) then
    call read_inopt(use_binary, 'use binary', db, errcount=nerr)
    call read_inopt(use_sink, 'use sink', db, errcount=nerr)
    if (use_binary) then
       call read_inopt(iorigin,'iorigin',db,min=0,max=2,errcount=nerr)
       call read_inopt(Mh2, 'Mh2',    db,min=0.,errcount=nerr)
       call read_inopt(ecc_binary, 'ecc_binary',    db,min=0.,max=1.,errcount=nerr)
       call read_inopt(semimajoraxis_binary, 'semimajoraxis_binary',    db,min=0.05,errcount=nerr)
    endif
 endif
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
