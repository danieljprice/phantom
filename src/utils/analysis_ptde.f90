!This file analysis the partial TDE. We use the mass of the star and the mass of the particles * code units?
!to find the mass lost in a PTDE.

module analysis
  implicit none

  character(len=20), parameter, public :: analysistype = 'ptde'
  logical, private :: firstcall = .true.

  public :: do_analysis

  private

  real, dimension(:), allocatable :: ebins,dmde,tbins,dmdt,rbins,dlumdr,lumcdf,lbins,dmdl,angbins,dmdang,vbins,dmdv

  !---- These can be changed in the params file
  integer :: nbins
  real :: mh     = 1.
  !--- If min=max then lims are set dynamically
  real :: emin   = 0.
  real :: emax   = 0.
  real :: trmin  = 0.
  real :: trmax  = 0.
  real :: rmin   = 0.
  real :: rmax   = 0.
  real :: lummin = 0.
  real :: lummax = 0.
  real :: angmin = 0.
  real :: angmax = 0.
  real :: vmin   = 0.
  real :: vmax   = 0.

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
  real(4) :: luminosity(npart)

  if (.not.opened_full_dump) then
     write(*,'("SKIPPING FILE -- (Not a full dump)")')
     return
  endif

  call read_array_from_file(123,dumpfile,'luminosity',luminosity(1:npart),ierr)
  if (ierr/=0) then
     write(*,*)
     write(*,'("WARNING: could not read luminosity from file. It will be set to zero")')
     write(*,*)
     luminosity = 0.
  endif

  ! Print the analysis being done
   write(*,'("Performing analysis type ",A)') analysistype
   write(*,'("Input file name is ",A)') dumpfile

   write(output,"(a4,i5.5)") 'hist',numfile
   write(*,'("Output file name is ",A)') output

   ! Assuming G=1
   write(*,*)
   write(*,'("ASSUMING G == 1")')

   ! Read black hole mass from params file
   filename = 'analysis_'//trim(analysistype)//'.params'
   inquire(file=filename,exist=iexist)
   !if (iexist) call read_tdeparams(filename,ierr)
   if (.not.iexist.or.ierr/=0) then
      nbins = int(sqrt(real(npart)))
      !call write_tdeparams(filename)
      print*,' Edit '//trim(filename)//' and rerun phantomanalysis'
      stop
   endif

   if (.not.allocated(ebins)) then
      allocate(ebins(nbins),dmde(nbins),tbins(nbins),dmdt(nbins),rbins(nbins),dlumdr(nbins),lumcdf(nbins),&
               lbins(nbins),dmdl(nbins),angbins(nbins),dmdang(nbins),vbins(nbins),dmdv(nbins))
   endif

   !call ptde_analysis(npart,pmass,xyzh,vxyzu,mass_remaining)
   open(iunit,file=output)
   write(iunit,'("# ",es20.12,"   # TIME")') time
   write(iunit, *) 'hey'


  end subroutine do_analysis


end module analysis
