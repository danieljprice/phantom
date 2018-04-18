module setup
 implicit none
 public :: setpart

 private

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
#ifdef GR
 use externalforces, only:iext_gr
#else
 use externalforces, only:iext_star
#endif
 use eos,            only:ieos
 use physcon,        only:pi
 use prompting,      only:prompt
 use vectorutils,    only:cross_product3D
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: i,dumpsperorbit,orbtype
 real    :: x0,y0,z0,vx0,vy0,vz0,dr,h0,xyz0(3),rhat(3),r2,vcirc,rtan(3),norbits,period,r0,fac

 call set_units(mass=solarm,G=1.d0,c=1.d0)

 !
 ! general parameters
 !
 time  = 0.
 gamma = 1.
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

#ifdef GR
 iexternalforce = iext_gr
#else
 iexternalforce = iext_star
#endif

 xyzh  = 0.
 vxyzu = 0.

 x0    = 0.
 y0    = 0.
 z0    = 0.
 vx0   = 0.
 vy0   = 0.
 vz0   = 0.

 orbtype       = 1
 norbits       = 1.
 dumpsperorbit = 100
 period        = 0.

 call prompt('select orbit type (1=cirlce, 2=precession, 3=epicycle, 4=vertical-oscillation, 0=custom)',orbtype,0,4)
 select case(orbtype)

 case(1) ! circular
    x0  = 10.
    call prompt('initial radius r0',x0)
    vy0 = sqrt(1./x0)
    period = 2.*pi*sqrt(x0**3/1.)

 case(2) ! precession
    x0     = 90.
    vy0    = 0.0521157
    period = 2.*pi*sqrt((0.5*x0)**3/1.) ! approximate

 case(3) ! epicycle
    x0  = 10.
    fac = 1.00001
    call prompt('initial radius r0',x0)
    vy0 = fac*sqrt(1./x0)
    period = 2.*pi*sqrt(x0**3/1.)

 case(4) ! vertical oscillation
    x0  = 10.
    fac = 1.00001
    z0  = fac-1.
    call prompt('initial radius r0',x0)
    vy0 = sqrt(1./x0)
    period = 2.*pi*sqrt(x0**3/1.)

 case(0) ! custom : if you just press enter a lot it gives you a circular orbit
    x0  = 10.
    vy0 = sqrt(1./x0)
    call prompt('initial x position',x0)
    call prompt('initial y position',y0)
    call prompt('initial z position',z0)
    call prompt('initial vx velocity',vx0)
    call prompt('initial vy velocity',vy0)
    call prompt('initial vz velocity',vz0)
    tmax  = 500.
    dtmax = tmax/1000.
    call prompt('enter tmax',tmax)
    call prompt('enter dtmax',dtmax)

 end select

 if (orbtype/=0) then
    call prompt('number of orbits',norbits)
    call prompt('dumps per orbit',dumpsperorbit)
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
 ! set smoothing lengths, and set thermal energies
 !
 xyz0          = xyzh(1:3,1)
 r0            = sqrt(dot_product(xyz0,xyz0))
 rhat          = xyz0/r0
 dr            = 0.25
 h0            = 10.*dr
 xyzh(4,:)     = h0
 vxyzu(4,:)    = 0.
 do i=2,npart
    xyzh(1:3,i) = xyz0 + (i-1)*dr*rhat
    call cross_product3D((/0.,0.,1./),xyzh(1:3,i),rtan)          ! Unit vector tangential to motion
    rtan  = rtan/sqrt(dot_product(rtan,rtan))
    r2    = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
    vcirc = sqrt(1./sqrt(r2))
    vxyzu(1:3,i) = rtan*vcirc
 enddo

end subroutine setpart

end module setup
