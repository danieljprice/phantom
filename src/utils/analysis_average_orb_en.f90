module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'average orbital energy'
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
 real    :: ekin_av,epot_av,e_av

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 call get_average_energies(npart,xyzh,vxyzu,ekin_av,epot_av,e_av)

 if (first) then
    first = .false.
    open(unit=iu, file='orbitalenergy.ev',status='replace')
    write(iu,"('#',4(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',&
          2,'ekin',&
          3,'epot',&
          4,'etot'
 else
    open(unit=iu, file='orbitalenergy.ev',position='append')
 endif
 write(iu,'(4(es18.10,1X))') time,ekin_av,epot_av,e_av
 close(iu)

end subroutine do_analysis

!--------------------------------------------------------------------------------------------------------------------
!
!-- Actual subroutine where the analysis is done!
!
!--------------------------------------------------------------------------------------------------------------------
subroutine get_average_energies(npart,xyzh,vxyzu,ekin_av,epot_av,e_av)
 use part, only:isdead_or_accreted
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:),vxyzu(:,:)
 real, intent(out)   :: ekin_av,epot_av,e_av
 integer :: i,n
 real    :: v2,r

 ekin_av = 0.
 epot_av = 0.
 n       = 0

 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       v2 = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
       r  = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
       ekin_av = ekin_av + 0.5*v2
       epot_av = epot_av - 1./r
       n = n + 1
    endif
 enddo

 ekin_av = ekin_av/n
 epot_av = epot_av/n

 e_av = ekin_av + epot_av

end subroutine get_average_energies

end module
