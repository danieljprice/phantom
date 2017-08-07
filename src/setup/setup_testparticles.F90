module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for blob problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use timestep,       only:dtmax,tmax
 use options,        only:iexternalforce,alpha,alphau,beta,nfulldump
 use units,          only:set_units
 use physcon,        only:solarm
#ifdef GR
 use externalforces, only:iext_gr
#else
 use externalforces, only:iext_star
#endif
 use eos,            only:ieos
 use physcon,        only:pi
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
 real    :: r2

 call set_units(mass=solarm,G=1.d0,c=1.d0)

 !
 ! general parameters
 !

 time  = 0.
 tmax  = 500.*5
 dtmax = 5
 gamma = 1.
 polyk = 0.
 npart = 20
 ieos  = 11
 nfulldump = 1

 alpha  = 0.000    ! art. viscosity parameter
 alphau = 0.000    ! art. conductivity parameter
 beta   = 0.000    ! beta viscosity

 massoftype     = 1.e-10
 npartoftype(:) = 0
 npartoftype(1) = npart

#ifdef GR
    iexternalforce = iext_gr
#else
    iexternalforce = iext_star
#endif

 !
 ! setup particle positions, velocities, thermal energy
 !

 xyzh = 0.
 vxyzu = 0.
 xyzh(4,:) = 10.
 vxyzu(4,:) = 0.

 do i=1,npart
    xyzh(1,i)    = 10. + (i-1)*0.001
    r2           = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
    vxyzu(1:3,i) = (/0.,sqrt(1./sqrt(r2)),0./)
 enddo

end subroutine setpart

end module setup
