!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection module for radiation pulse setup to reset radiation energy of boundary particles
!
! :References: None
!
! :Owner:
!
! :Runtime parameters:
!
! :Dependencies:
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'radpulse'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
           set_default_options_inject,update_injected_par

 real, public :: xileft = 4.
 real, public :: xiright = 0.4
 private

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use part, only:igas,iboundary,iamtype,iphase,iradxi,rad
 use units, only:unit_ergg
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real                   :: xileftcode,xirightcode
 integer :: i

 ! fix radiation energy of boundary by resetting to original values
 xileftcode = xileft/unit_ergg
 xirightcode = xiright/unit_ergg
 do i=1,npart
    if (.not. iamtype(iphase(i))==iboundary) then
       cycle
    elseif (xyzh(1,i) < 0.) then
       rad(iradxi,i) = xileftcode
    else
       rad(iradxi,i) = xirightcode
    endif
 enddo

 dtinject = huge(0.) ! do not limit timestep

end subroutine inject_particles

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(xileft,'xileft','xileft',iunit)
 call write_inopt(xiright,'xiright','xiright',iunit)
end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only: fatal, error, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('xileft')
    read(valstring,*,iostat=ierr) xileft
    ngot = ngot + 1
 case('xiright')
    read(valstring,*,iostat=ierr) xiright
    ngot = ngot + 1
 end select

 igotall = (ngot >= 2)
end subroutine read_options_inject

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
