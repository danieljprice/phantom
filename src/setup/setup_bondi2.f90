module setup
 implicit none

 public :: setpart

 private

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon, only:solarm
 use units,   only:set_units
 use part,    only:igas
 use timestep, only:tmax,dtmax
 use externalforces,only:accradius1,accradius1_hard
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)

 call set_units(c=1.,G=1.)

 gamma = 5./3.

 massoftype(igas) = 1.e-6 !mdot/dndt

 tmax = 1000.
 dtmax = 5.

 accradius1 = 2.5
 accradius1_hard = accradius1

end subroutine setpart

end module setup
