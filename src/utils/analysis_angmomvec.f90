module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'angmomvec'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 integer, parameter :: iu = 1993
 logical, save      :: first = .true.
 real    :: Lhat(3),inc,rot

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 call get_angmomvec(npart,xyzh,vxyzu,Lhat,inc,rot)

 ! Write angular momentum vector information
 if (first) then
    first = .false.
    open(unit=iu, file='angmomvec.ev',status='replace')
    write(iu,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',&
          2,'Lx',  &
          3,'Ly',  &
          4,'Lz',  &
          5,'inc', &
          6,'rot'
 else
    open(unit=iu, file='angmomvec.ev',position='append')
 endif
 write(iu,'(6(es18.10,1X))') time,Lhat,inc,rot
 close(iu)

end subroutine do_analysis

!--------------------------------------------------------------------------------------------------------------------
!
!-- Actual subroutine where the analysis is done!
!
!--------------------------------------------------------------------------------------------------------------------
subroutine get_angmomvec(npart,xyzh,vxyzu,Lhat,inc,rot)
 use physcon,     only:pi
 use vectorutils, only:cross_product3D
 use part,        only:isdead_or_accreted
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:),vxyzu(:,:)
 real, intent(out)   :: Lhat(3),inc,rot
 integer :: i
 real    :: Li(3),Ltot(3)

 Ltot = 0.
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),Li)
       Ltot = Ltot + Li
    endif
 enddo

 Lhat = Ltot/sqrt(dot_product(Ltot,Ltot))
 inc  = acos(dot_product( Lhat,(/0.,0.,1./)))*180./pi ! Angle from +z axis -- should always be 0<inc<180 degrees
 rot  = atan2(Lhat(2),Lhat(1))*180./pi                ! Angle around in xy plane (from +x axis) -- should always -180<rot<180

end subroutine get_angmomvec

end module
