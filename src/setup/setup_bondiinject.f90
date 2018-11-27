module setup
 implicit none
 public :: setpart

 private

 logical :: filldomain = .false.

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma_eos,hfact,time,fileprefix)
 use part,           only:igas,gr,xyzmh_ptmass,vxyz_ptmass
 use options,        only:iexternalforce
 use units,          only:set_units
 use inject,         only:init_inject,inject_particles,dtsphere,rin
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
 integer :: ierr,nspheres
 real :: dtinject,tinfall

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

 if (filldomain) then
!
!--- Inject 'real' spheres into the whole domain
!    (get the number of spheres required self-consistently from the infall time)
!
    call get_tinfall(tinfall,r1=accradius1,r2=rin,gamma=gamma)
    nspheres = int(tinfall/dtsphere) !27!100!20!
    print*,'number of "real" spheres: ',nspheres
    call inject_particles(dtsphere*nspheres,dtsphere*nspheres,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,dtinject)
 endif

end subroutine setpart

!
!-- Basic integration of the velocty to get the infall time between two points
!
subroutine get_tinfall(tinfall,r1,r2,gamma)
 use bondiexact, only:get_bondi_solution
 real, intent(in)  :: r1,r2,gamma
 real, intent(out) :: tinfall
 integer, parameter :: N=10000
 integer :: i
 real :: dr,rhor,vr,ur,r
 dr = abs(r2-r1)/N
 tinfall = 0.

 r = r1 + dr
 do i=1,N
    call get_bondi_solution(rhor,vr,ur,r,mass=1.,gamma_eos=gamma)
    tinfall = tinfall + dr/abs(vr)
    r = r + dr
 enddo

end subroutine get_tinfall

end module setup
