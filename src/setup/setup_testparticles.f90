!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! setup for test particles
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - dumpsperorbit : *dumps per orbit*
!   - norbits       : *number of orbits*
!   - orbtype       : *orbit type (1=circle, 2=precession, 3=epicycle, 4=vertical-oscillation, 0=custom)*
!   - r             : *initial radius r in spherical coordinates*
!   - spin          : *black hole spin*
!   - vx0           : *initial vx velocity*
!   - vy0           : *initial vy velocity*
!   - vz0           : *initial vz velocity*
!   - x0            : *initial x position*
!   - y0            : *initial y position*
!   - z0            : *initial z position*
!
! :Dependencies: eos, externalforces, infile_utils, io, options, part,
!   physcon, prompting, timestep, units, vectorutils
!
 implicit none
 public :: setpart

 private

 ! Module variables for setup parameters
 integer :: orbtype
 integer :: dumpsperorbit
 real :: spin
 real :: r
 real :: norbits
 real :: x0, y0, z0
 real :: vx0, vy0, vz0

contains

!----------------------------------------------------------------
!+
!  Setup for a single test particles (no pressure between gas)
!
!  User has the choice to select the initial positon and velocity
!  of particle 1.
!
!   - The particle is of type gas, however it feels only the external force.
!   - As a consequence of using a gas particle, we need to also initialise a
!     few extra gas particles, in order for neighbour finding to not fail.
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use timestep,       only:dtmax,tmax
 use options,        only:iexternalforce,alpha,alphamax,alphau,beta,nfulldump
 use units,          only:set_units
 use physcon,        only:solarm
 use io,             only:master
 use externalforces, only:iext_star,a
 use eos,            only:ieos
 use physcon,        only:pi
 use prompting,      only:prompt
 use vectorutils,    only:cross_product3D
 use part,           only:gr
 use infile_utils,   only:get_options
 use kernel,         only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: i,ierr
 real    :: dr,h0,xyz0(3),rhat(3),r2,vcirc,rtan(3),period,r0,fac,omega,rkerr,z

 call set_units(mass=solarm,G=1.d0,c=1.d0)
 !
 ! default runtime parameters
 !
 orbtype = 1
 spin = 0.
 r = 10.
 norbits = 1.
 dumpsperorbit = 100
 x0 = 0.
 y0 = 0.
 z0 = 0.
 vx0 = 0.
 vy0 = 0.
 vz0 = 0.
 tmax = 500.
 dtmax = 0.5
 period = tmax
 hfact = hfact_default
 !
 ! read runtime parameters from setup file
 !
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Test particles setup'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
 !
 ! general parameters
 !
 time  = 0.
 if (gr) then
    gamma = 5./3. ! GR cannot have gamma=1
 else
    gamma = 1.
 endif
 polyk = 0.
 npart = 10
 ieos  = 11
 nfulldump = 1
 alpha     = 0.
 alphamax  = 0.
 alphau    = 0.
 beta      = 0.
 massoftype     = 1.e-10
 npartoftype(:) = 0
 npartoftype(1) = npart

 xyzh  = 0.
 vxyzu = 0.

 ! Calculate orbit parameters based on orbit type
 select case(orbtype)
 case(1,3,4) ! circular, epicycle, vertical oscillation
    omega = 1./(r**(1.5)+spin)
    x0 = sqrt(r**2 + spin**2)
    vy0 = x0*omega
    period = 2.*pi/omega
    if (orbtype == 3) then ! epicycle
       fac = 1.00001
       vy0 = fac*x0*omega
    elseif (orbtype == 4) then ! vertical oscillation
       fac = 1.00001
       z0  = fac-1.
    endif
 case(2) ! precession
    x0     = 90.
    vy0    = 0.0521157
    period = 2.*pi*sqrt((0.5*x0)**3/1.) ! approximate
 case(0) ! custom
    x0  = 10.
    vy0 = sqrt(1./x0)
 end select

 if (orbtype > 0) then
    tmax   = norbits*period
    dtmax  = period/dumpsperorbit
 endif

 print*,''
 print*,' Setup for single test particle: '
 print*,' tmax = ',tmax
 print*,' Initial  (x,y,z)   = ',x0,y0,z0
 print*,' Initial (vx,vy,vz) = ',vx0,vy0,vz0
 print*,''

 xyzh(1:3,1)   = (/x0,y0,z0/)
 vxyzu(1:3,1)  = (/vx0,vy0,vz0/)
 !
 ! Put all other particles in a radial line outwards from the origin, with their circular velocity,
 ! but only in the x-y plane. Also set smoothing lengths, and thermal energies.
 !
 xyz0          = (/xyzh(1,1),xyzh(2,1),0./)
 r0            = sqrt(dot_product(xyz0,xyz0))
 rhat          = xyz0/r0
 dr            = 0.25
 h0            = 10.*dr
 xyzh(4,:)     = h0
 vxyzu(4,:)    = 0.
 do i=2,npart
    xyzh(1:3,i) = xyz0 + (i-1)*dr*rhat
    call cross_product3D((/0.,0.,1./),xyzh(1:3,i),rtan)          ! Vector tangential to motion
    rtan  = rtan/sqrt(dot_product(rtan,rtan))                    ! Unit vector tangential to motion
    r2    = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
    z     = xyzh(3,i)
    rkerr = sqrt((r2-spin**2)/2. + sqrt((r2-spin**2)**2 + 4.*spin**2*z**2)/2.)
    omega = 1./(rkerr**(1.5)+spin)
    x0    = sqrt(rkerr**2 + spin**2)
    vcirc = x0*omega
    vxyzu(1:3,i) = rtan*vcirc
 enddo

 a = spin
 if (.not.gr) iexternalforce = iext_star

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for test particles setup'
 call write_inopt(orbtype,'orbtype','orbit type (1=circle, 2=precession, 3=epicycle, 4=vertical-oscillation, 0=custom)',iunit)
 call write_inopt(spin,'spin','black hole spin',iunit)
 select case(orbtype)
 case(1,3,4) ! circular, epicycle, vertical oscillation
    call write_inopt(r,'r','initial radius r in spherical coordinates',iunit)
 case(0) ! custom
    call write_inopt(x0,'x0','initial x position',iunit)
    call write_inopt(y0,'y0','initial y position',iunit)
    call write_inopt(z0,'z0','initial z position',iunit)
    call write_inopt(vx0,'vx0','initial vx velocity',iunit)
    call write_inopt(vy0,'vy0','initial vy velocity',iunit)
    call write_inopt(vz0,'vz0','initial vz velocity',iunit)
 end select
 if (orbtype > 0) then
    call write_inopt(norbits,'norbits','number of orbits',iunit)
    call write_inopt(dumpsperorbit,'dumpsperorbit','dumps per orbit',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(orbtype,'orbtype',db,min=0,max=4,errcount=nerr)
 call read_inopt(spin,'spin',db,min=-1.,max=1.,errcount=nerr)
 select case(orbtype)
 case(1,3,4) ! circular, epicycle, vertical oscillation
    call read_inopt(r,'r',db,min=0.,errcount=nerr)
    call read_inopt(norbits,'norbits',db,min=0.,errcount=nerr)
    call read_inopt(dumpsperorbit,'dumpsperorbit',db,min=1,errcount=nerr)
 case default ! custom
    call read_inopt(x0,'x0',db,errcount=nerr)
    call read_inopt(y0,'y0',db,errcount=nerr)
    call read_inopt(z0,'z0',db,errcount=nerr)
    call read_inopt(vx0,'vx0',db,errcount=nerr)
    call read_inopt(vy0,'vy0',db,errcount=nerr)
    call read_inopt(vz0,'vz0',db,errcount=nerr)
 end select
 if (orbtype > 0) then
    call read_inopt(norbits,'norbits',db,min=0.,errcount=nerr)
    call read_inopt(dumpsperorbit,'dumpsperorbit',db,min=1,errcount=nerr)
 endif
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Interactive setup
!+
!-----------------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 use part,      only:gr

 if (gr) call prompt('black hole spin',spin,-1.,1.)
 call prompt('select orbit type (1=circle, 2=precession, 3=epicycle, 4=vertical-oscillation, 0=custom)',orbtype,0,4)
 select case(orbtype)
 case(1,3,4) ! circular, epicycle, vertical oscillation
    call prompt('initial radius r in spherical (this is not that same as radius (x0,0,0) in Cartesian)',r)
 case(0) ! custom
    call prompt('initial x position',x0)
    call prompt('initial y position',y0)
    call prompt('initial z position',z0)
    call prompt('initial vx velocity',vx0)
    call prompt('initial vy velocity',vy0)
    call prompt('initial vz velocity',vz0)
 end select
 if (orbtype > 0) then
    call prompt('number of orbits',norbits)
    call prompt('dumps per orbit',dumpsperorbit)
 endif

end subroutine setup_interactive

end module setup
