module inject
 use physcon, only:solarm,years
 use units,   only:utime,umass
 implicit none
 character(len=*), parameter, public :: inject_type = 'bondi'

 public :: inject_particles,write_options_inject,read_options_inject

 real, public :: dmdt = 1. ! mass injection rate
 real, public :: dndt = 5.              ! particle injection rate in particles per orbit
 real, public :: rin  = 30.0                ! injection radius

 private

 !-- Global to module, only need to be computed once since radius injection is always the same.
 real :: rhoin,vin,uin,hin

contains

! Wrapper
subroutine get_exact(rin,mass1,gamma)
 use bondiexact, only:get_bondi_solution
 use part, only:hrho
 real, intent(in)  :: rin,mass1,gamma
 call get_bondi_solution(rhoin,vin,uin,rin,mass1,gamma)
 hin = hrho(rhoin)
 ! vin = -vin
 ! ! Direction of wind
 ! if (inflow) v = -v
end subroutine get_exact

subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,      only:massoftype,igas,hrho
 use partinject,only:add_or_update_particle
 use physcon,   only:pi,twopi
 use random,    only:ran2
 use eos,       only:gamma
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    dimension(3)  :: xyz,vxyz,rhat
 real    :: theta,phi
 integer :: i,ipart,npinject,seed
 logical, save :: first = .true.

!
!-- Mass of gas particles is set by mass accretion rate and particle injection rate
!
 massoftype(igas) = dmdt/dndt

 if (first) then
    call get_exact(rin,mass1=1.,gamma=gamma)
 endif

!
!-- How many particles do we need to inject?
!   (Seems to need at least eight gas particles to not crash) <-- This statement may or may not be true...
!
 if(npartoftype(igas)<8) then
    npinject = 8-npartoftype(igas)
 else
    npinject = max(0, int(0.5 + (time*dmdt/massoftype(igas)) - npartoftype(igas) ))
 endif

 print*,'injecting ',npinject

!
!-- Randomly inject particles at the injection radius
!
 do i=1,npinject
    phi       = ran2(seed)*twopi
    theta     = ran2(seed)*pi
    xyz       = [rin*cos(phi)*sin(theta),rin*sin(phi)*sin(theta),rin*cos(theta)]
    rhat      = xyz/rin
    vxyz      = vin*rhat
    ipart     = npart + 1
    call add_or_update_particle(igas,xyz,vxyz,hin,uin,ipart,npart,npartoftype,xyzh,vxyzu)
 enddo

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 ! use infile_utils, only:write_inopt
 integer, intent(in) :: iunit
 !
 ! call write_inopt(dmdt,'dmdt','mass injection rate in grams/second',iunit)
 ! call write_inopt(dndt,'dndt','particle injection rate'            ,iunit)
 ! call write_inopt(rin ,'rin' ,'injection radius'                   ,iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 ! use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 ! integer, save :: ngot = 0
 ! character(len=30), parameter :: label = 'read_options_inject'
 !
 ! imatch  = .true.
 ! select case(trim(name))
 ! case('dmdt')
 !    read(valstring,*,iostat=ierr) dmdt
 !    ngot = ngot + 1
 !    if (dmdt  <  0.) call fatal(label,'dmdt < 0 in input options')
 ! case('dndt')
 !    read(valstring,*,iostat=ierr) dndt
 !    ngot = ngot + 1
 !    if (dndt < 0.) call fatal(label,'dndt < 0 in input options')
 ! case('rin')
 !    read(valstring,*,iostat=ierr) rin
 !    ngot = ngot + 1
 ! case default
 !    imatch = .false.
 ! end select
 !
 ! igotall = (ngot >= 1)

end subroutine read_options_inject

end module inject
