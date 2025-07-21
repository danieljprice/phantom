!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module ptmass_heating
!
! Heating of particles around sink particles
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: kernel, part
!

 implicit none
 public :: energ_sinkheat,heating_kernel
 real,    public :: Lnuc
 integer, public :: isink_heating

 private

contains
!-----------------------------------------------------------------------
!+
!  heating from point mass
!+
!-----------------------------------------------------------------------
subroutine energ_sinkheat(nptmass,xyzmh_ptmass,xi,yi,zi,dudtheati)
 use part,   only:ihsoft,imassenc,iLum, sink_has_heating
 use kernel, only:radkern2
 integer, intent(in) :: nptmass
 real, intent(in)    :: xi,yi,zi,xyzmh_ptmass(:,:)
 real, intent(out)   :: dudtheati
 integer             :: i
 real                :: q2,dri2

 dudtheati = 0.
 do i = 1,nptmass
    if (.not. sink_has_heating(xyzmh_ptmass(:,i))) cycle
    dri2 = (xi-xyzmh_ptmass(1,i))**2 + (yi-xyzmh_ptmass(2,i))**2 + (zi-xyzmh_ptmass(3,i))**2
    q2 = dri2/xyzmh_ptmass(ihsoft,i)**2
    if (q2 < radkern2) then
       ! this should deliberately crash if xyzmh_ptmass(imassenc,i) == 0, instead of ignoring that sink's heating
       dudtheati = dudtheati + xyzmh_ptmass(iLum,i) / xyzmh_ptmass(imassenc,i) * heating_kernel(q2,isink_heating)
    endif
 enddo

end subroutine energ_sinkheat


!-----------------------------------------------------------------------
!+
!  heating weight function (note: arbitrary normalisation)
!+
!-----------------------------------------------------------------------
real function heating_kernel(q2,kernel_type)
 use kernel, only:wkern,cnormk
 real, intent(in)    :: q2
 integer, intent(in) :: kernel_type

 select case(kernel_type)
 case(1)
    heating_kernel = cnormk*wkern(q2,sqrt(q2))
 case default
    heating_kernel = 1.
 end select

end function heating_kernel


!-----------------------------------------------------------------------
!+
!  write options to input file (not used at the moment)
!+
!-----------------------------------------------------------------------
! subroutine write_options_ptmass_heating(iunit)
!  use infile_utils, only:write_inopt
!  integer, intent(in) :: iunit

!  call write_inopt(isink_heating,'isink_heating','sink heating distirbution (0=uniform,1=kernel)',iunit)

! end subroutine write_options_ptmass_heating


!-----------------------------------------------------------------------
!+
!  read options from input file (not used at the moment)
!+
!-----------------------------------------------------------------------
! subroutine read_options_ptmass_heating(name,valstring,imatch,igotall,ierr)
!  use io,  only:fatal
!  character(len=*), intent(in)  :: name,valstring
!  logical, intent(out)          :: imatch,igotall
!  integer, intent(out)          :: ierr
!  integer, save                 :: ngot = 0
!  integer                       :: ni
!  character(len=30), parameter  :: label = 'read_options_ptmass_heating'

!  imatch  = .true.
!  igotall = .false.
!  select case(trim(name))
!  case('isink_heating')
!     read(valstring,*,iostat=ierr) isink_heating
!     ngot = ngot + 1
!     if (isink_heating < 0 .or. isink_heating > 1) call fatal(label,'invalid setting for isink_heating ([0,1])')
!  case default
!     imatch = .false.
!  end select
!  ni = 1
!  igotall = (ngot >= ni)

! end subroutine read_options_ptmass_heating


end module ptmass_heating
