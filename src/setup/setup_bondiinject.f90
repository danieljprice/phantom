module setup
 implicit none
 public :: setpart

 private

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma_eos,hfact,time,fileprefix)
 use part,           only:igas,gr
 use options,        only:iexternalforce
 use units,          only:set_units
 use inject,         only:init_inject
 use timestep,       only:tmax
 use io,             only:master
 use eos,            only:gamma
 use metric,         only:imetric
 use metric_tools,   only:imet_schwarzschild
 use externalforces, only:accradius1,accradius1_hard
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma_eos,hfact
 real,              intent(inout) :: time
 character(len=*),  intent(in)    :: fileprefix
 logical :: iexist
 integer :: ierr

 if (.not.gr) call fatal('setup','This setup only works with GR on')
 if (imetric/=imet_schwarzschild) call fatal('setup','This setup is meant for use with the Schwarzschild metric')
 call set_units(G=1.,c=1.)

 time            = 0.
 polyk           = 0.
 iexternalforce  = 1

 !-- Don't overwrite these values if infile exists
 inquire(file=trim(fileprefix)//'.in',exist=iexist)
 if (.not.iexist) then
    tmax            = 360.
    accradius1      = 2.1
    accradius1_hard = accradius1
 endif

 npart            = 0
 npartoftype(:)   = 0
 massoftype(igas) = 0.012 !1.30833862e-2
 gamma            = 5./3. !Set gamma in module eos since init_inject needs to know about it.
 gamma_eos        = gamma

 call init_inject(ierr)

end subroutine setpart

end module setup
