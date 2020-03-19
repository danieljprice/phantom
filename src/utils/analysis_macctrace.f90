module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'macc'
 public :: do_analysis

 private

 integer, allocatable :: already_accreted(:)
 integer :: iacc = 0

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 integer :: i
 logical :: fopened

 ! Allocate and initialise accretion history array
 if (.not.allocated(already_accreted)) then
    allocate(already_accreted(npart))
    already_accreted = 0
 endif

 fopened = .false.

 ! Trying to save the accretion history -- so you need to run this on all files in order, starting before anything got accreted.

 ! Check if a particle is accreted.
 ! Then check to see if the particle is already in the list of accreted things
 ! If it's not accreted, then add it to the list of accreted things

 do i=1,npart
    if (xyzh(4,i)<0.) then
       if (any(i==already_accreted)) cycle
       if (.not.fopened) call open_newfile(iunit,numfile,time,fopened)
       iacc = iacc + 1              ! Increase the counter for number of particles accreted
       already_accreted(iacc) = i   ! Save the particle to the accretion history
       write(iunit,*) i             ! Write the particle to the file of the current dump
    endif
 enddo

end subroutine do_analysis

subroutine open_newfile(iunit,numfile,time,fopened)
 integer, intent(in)  :: iunit, numfile
 real,    intent(in)  :: time
 logical, intent(out) :: fopened
 character(len=120) :: output

 write(output,"(a5,i5.5)") 'macc_',numfile
 write(*,'("Output file name is ",A)') trim(output)

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',1(1x,'[',i2.2,1x,a34,']',2x))") 1,'accreted particles since last dump'

 fopened = .true.

end subroutine open_newfile

end module
