!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!   Analysis routine to call MCFOST code to perform live temperature
!   feedback from radiation transport
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, io
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'mcfost'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use eos, only : temperature_coef, gmw, gamma
 use io,             only:fatal
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(in)    :: particlemass,time
 real,             intent(inout) :: vxyzu(:,:)
 real(kind=4), dimension(:), allocatable :: T_SPH
 integer :: n_SPH, file_size, ierr, i, stat
 logical :: file_exists


 ! call mcfost
 call execute_command_line('rm -r data_th _voronoi.tmp')
 write(*,*) 'Calling mcfost on dumpfile: ', dumpfile
 call execute_command_line('./mcfost disc.para -phantom ' // trim(dumpfile), exitstat=stat)
 if (stat /= 0) call fatal('mcfost-phantom','error running mcfost from cmd line')
 ! note: mcfost does not return a non-zero status if it fails
 ! Christophe may change this in the future.

 ! check if file exists
 inquire(file='T_for_phantom.tmp',exist=file_exists)
 if (file_exists) then

    ! get number of SPH particles with temperatures
    inquire(file='T_for_phantom.tmp',size=file_size)
    ! file_size is size of file, 4 bytes per real
    ! with 1 at the top and 1 at the bottom
    n_SPH = file_size/4 - 2
    allocate(T_SPH(n_SPH))

    ! read values from unformatted fortran file
    open(unit=1,file='T_for_phantom.tmp',status='old',form='unformatted',&
         action='read',iostat=ierr)
    read(unit=1) T_SPH
    close(unit=1)
    write(*,*) ''
    write(*,*) 'Minimum temperature = ', minval(T_SPH, mask=(T_SPH > 0.))
    write(*,*) 'Maximum temperature = ', maxval(T_SPH)
    write(*,*) ''

    ! set thermal energy
    do i=1,n_SPH
       if (T_SPH(i) > 0.) then
          vxyzu(4,i) = T_SPH(i)/(temperature_coef*gmw*(gamma-1))
       else
          ! if mcfost doesn't return a temperature set it to 10K
          vxyzu(4,i) = 10./(temperature_coef*gmw*(gamma-1))
       endif
    enddo

 else

    write(*,*) ''
    write(*,*) 'ERROR: cannot open temperature file from mcfost'
    write(*,*) ''

 endif

end subroutine do_analysis

end module
