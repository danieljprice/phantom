!!This code reads the data from the TDE dumpfile and stores the data in the format liked by kepler.
!!Flaws in this code- Does not consider bins. Instead the total npart. We need to bin data so that we
!!have grid similar to input file??

!
!--Writing the analysis module for making a kepler file.
module analysis
  implicit none
  character(len=3), parameter, public :: analysistype = 'tde'
  public :: do_analysis

private

!---- These can be changed in the params file
integer :: nbins
real :: mh     = 1.
contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

   use io,         only:warning
   use dump_utils, only:read_array_from_file
   use prompting,  only:prompt
   use readwrite_dumps, only: opened_full_dump
   character(len=*),   intent(in) :: dumpfile
   integer,            intent(in) :: numfile,npart,iunit
   real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
   real,               intent(in) :: pmass,time
   character(len=120) :: output
   character(len=20)  :: filename
   integer :: i,ierr
   logical :: iexist
   real(4) :: pressure(npart)
   real(4) :: density(npart)
   real(4) :: temperature(npart)

   if (.not.opened_full_dump) then
      write(*,'("SKIPPING FILE -- (Not a full dump)")')
      return
   endif

   !Read density of each particle caculated.
   call read_array_from_file(123,dumpfile,'rhogas',density(1:npart),ierr)
   if (ierr/=0) then
      write(*,*)
      write(*,'("WARNING: could not read density from file. It will be set to zero")')
      write(*,*)
      density = 0.
   endif
   print*, 'density found'

   !Read pressure of each particle calculated.
   call read_array_from_file(123,dumpfile,'pressure',pressure(1:npart),ierr)
   if (ierr/=0) then
      write(*,*)
      write(*,'("WARNING: could not read pressure from file. It will be set to zero")')
      write(*,*)
      pressure = 0.
   endif
   print*, 'pressure found'

   !Read temperature of each particle calculated.
   call read_array_from_file(123,dumpfile,'temperature',temperature(1:npart),ierr)
   if (ierr/=0) then
      write(*,*)
      write(*,'("WARNING: could not read temperature from file. It will be set to zero")')
      write(*,*)
      pressure = 0.
   endif
   print*, 'temperature found'

   ! Print the analysis being done
    write(*,'("Performing analysis type ",A)') analysistype
    write(*,'("Input file name is ",A)') dumpfile

    write(output,"(a4,i5.5)") 'test',numfile
    write(*,'("Output file name is ",A)') output


    ! Assuming G=1
    write(*,*)
    write(*,'("ASSUMING G == 1")')

    !Use EOS to find energy of each particle?
    open(iunit,file=output)

    print*, output,'output'
    write(iunit,'("# ",es20.12,"   # TIME")') time !Don't need time for the kepler file but just testing things.
    write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'density',      &
          2,'temperature',  &
          3,'pressure'

    do i = 1,npart
       write(iunit,'(3(es18.10,1X))') &
              density(i),   &
              temperature(i),    &
              pressure(i)
    enddo

 end subroutine do_analysis


end module analysis
