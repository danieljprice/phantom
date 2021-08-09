!!This code reads the data from the TDE dumpfile and stores the data in the format liked by kepler.
!!Flaws in this code- Does not consider bins. Instead the total npart. We need to bin data so that we
!!have grid similar to input file??
!!Another problem is 0 values.
!
!--Writing the analysis module for making a kepler file.
module analysis
  implicit none
  character(len=3), parameter, public :: analysistype = 'tde'
  public :: do_analysis

private

!---- These can be changed in the params file
!integer :: nbins
!real :: mh     = 1.

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

   use io,         only:warning
   use dump_utils, only:read_array_from_file
   use units     , only:udist,umass,unit_density,unit_pressure,unit_ergg,unit_velocity !units required to convert to kepler units.
   use prompting,  only:prompt
   use readwrite_dumps, only: opened_full_dump
   use vectorutils, only:cross_product3D

   character(len=*),   intent(in) :: dumpfile
   integer,            intent(in) :: numfile,npart,iunit
   real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
   real,               intent(in) :: pmass,time
   character(len=120)             :: output
   !character(len=20)              :: filename

   integer :: i,ierr
   !logical :: iexist
   real(4) :: pressure(npart)
   real(4) :: density(npart)
   real(4) :: temperature(npart)
   real(4) :: entropy(npart)
   real(4) :: int_eng(npart)
   real(4) :: ang_vel(npart)
   real, dimension(npart) :: mass, velocity, outer_rad, tot_mass
   real    :: gamma = 5./3.
   real    :: Li(3) !defining angular momentum vector
   real , dimension(npart):: Li_mag !|Li(3)|

   !If dumpfile is not a complete dump we don't read it.
   if (.not.opened_full_dump) then
      write(*,'("SKIPPING FILE -- (Not a full dump)")')
      return
   endif

   !read the particle mass, velocity, radius and save it as an array
   do i = 1, npart
     !set the mass as the particle mass.
     mass(i) = pmass !every particle has same mass.

     !calculate the position which is the location of the particle.
     outer_rad(i) = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))

     !calculate the velocity using the magnitude of the velocity vector i.e., velocity = sqrt(v.v) where v is a vector.
     velocity(i) = sqrt(dot_product(vxyzu(1:3,i),vxyzu(1:3,i)))
   end do

   do i = 1, npart
      !Calculating angular momentum and saving it as an array.
      !Calculate the (r_x, r_y, r_z) X (v_x, v_y, v_z), ignoring mass as it gets cancelled when equating to I omega
      call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),Li)
      Li_mag(i) = sqrt(dot_product(Li, Li))
   end do

   !from angular momentum, calculate angular velocity
   do i = 1, npart
     ang_vel(i) = Li_mag(i)/outer_rad(i)**2
   end do

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

   !Read temperature of each particle calculated.
   call read_array_from_file(123,dumpfile,'entropy',entropy(1:npart),ierr)
   if (ierr/=0) then
      write(*,*)
      write(*,'("WARNING: could not read entropy from file. It will be set to zero")')
      write(*,*)
      pressure = 0.
   endif
   print*, 'entropy found'

   !Calculate internal energy
   !one problem if density is 0. fine to leave it like this??

   do i = 1, npart
     int_eng(i) = pressure(i)/((gamma - 1)*density(i))
   end do
   print*, 'energy calculated'

   ! Print the analysis being done
    write(*,'("Performing analysis type ",A)') analysistype
    write(*,'("Input file name is ",A)') dumpfile

    write(output,"(a4,i5.5)") 'test',numfile
    write(*,'("Output file name is ",A)') output


    ! Assuming G=1
    write(*,*)
    write(*,'("ASSUMING G == 1")')

    !open the output file and save the data in the format kepler likes.
    open(iunit,file=output)
    write(iunit,'("# ",es20.12,"   # TIME")') time !Don't need time for the kepler file but just testing things.
    write(iunit,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'mass[g]',               &  !mass of particle
          2,'position[cm]',          &  !position
          3,'vel[cm/s]',             &  !velocity
          4,'den[g/cm^3]',           &  !density
          5,'temp[K]',               &  !temperature
          6,'press[g/(cm s^2)]',     &  !pressure
          7,'entropy',               &  !entropy
          8,'int.eng[erg/g]' ,       &  !specific internal energy
          9,'ang vel'                  !angular velocity

    do i = 1,npart
       write(iunit,'(9(es18.10,1X))') &                             
              mass(i)*umass,                   &
              outer_rad(i)*udist,              &
              velocity(i)*unit_velocity,       &
              density(i)*unit_density,         &
              temperature(i),                  &
              pressure(i)*unit_pressure,       &
              entropy(i),                      & !no entropy unit?
              int_eng(i)*unit_ergg,            &
              ang_vel(i)
    enddo

 print*, sum(mass*umass), 'tot mass?'
 end subroutine do_analysis


end module analysis
