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
 integer :: i
 real    :: x0,y0,z0,vx0,vy0,vz0,dr,h0,xyz0(3),rhat(3),r2,vcirc,rtan(3)

 call set_units(mass=solarm,G=1.d0,c=1.d0)

 !
 ! general parameters
 !

 time  = 0.
 tmax  = 500.
 dtmax = tmax/1000
 gamma = 1.
 polyk = 0.
 npart = 10
 ieos  = 11
 nfulldump = 1

 alpha    = 0.    ! art. viscosity parameter
 alphamax = 0.
 alphau   = 0.    ! art. conductivity parameter
 beta     = 0.    ! beta viscosity

 massoftype     = 1.e-10
 npartoftype(:) = 0
 npartoftype(1) = npart

#ifdef GR
 iexternalforce = iext_gr
#else
 iexternalforce = iext_star
#endif


 x0 = 10.
 y0 = 0.
 z0 = 0.

 vx0 = 0.
 vy0 = sqrt(1./x0)
 vz0 = 0.

 call prompt('initial x position',x0)
 call prompt('initial y position',y0)
 call prompt('initial z position',z0)

 call prompt('initial vx velocity',vx0)
 call prompt('initial vy velocity',vy0)
 call prompt('initial vz velocity',vz0)

 print*,''
 print*,' Setup for single test particle: '
 print*,' tmax = ',tmax
 print*,' Initial  (x,y,z)   = ',x0,y0,z0
 print*,' Initial (vx,vy,vz) = ',vx0,vy0,vz0
 print*,''

 dr = 0.25
 h0 = 10.*dr

 xyzh = 0.
 vxyzu = 0.
 xyzh(4,:) = h0
 vxyzu(4,:) = 0.

 xyz0 = (/x0,y0,z0/)
 rhat = xyz0/sqrt(dot_product(xyz0,xyz0))

 xyzh(1:3,1) = xyz0
 vxyzu(1:3,1) = (/vx0,vy0,vz0/)

 !--- Put all other particles in a radial line outwards from the origin, with their circular velocity
 do i=2,npart
    xyzh(1:3,i)    = xyz0 + (i-1)*dr*rhat

    !--- Unit vector tangential to motion
    call cross_product3D((/0.,0.,1./),xyzh(1:3,i),rtan)
    rtan = rtan/sqrt(dot_product(rtan,rtan))

    r2           = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
    vcirc = sqrt(1./sqrt(r2))
    vxyzu(1:3,i) = rtan*vcirc
 enddo

end subroutine setpart

end module setup
