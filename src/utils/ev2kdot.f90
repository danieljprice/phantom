!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program get_kdot
!
! Utility to compute the time derivative of the total entropy
! from the .ev file(s), writes output to .kdot files
!
! :References: None
!
! :Owner: David Liptai
!
! :Usage: ev2kdot [dt_min] file01.ev file02.ev ...
!
! :Dependencies: ev2dotutils
!
 use ev2dotutils, only:ev2dot
 implicit none
 integer, parameter :: max_files = 128
 character(len=120) :: string,filenames(max_files)
 integer :: i,istart,iend,ierr,nargs,nfiles
 real    :: dtmin
 !
 ! get number of arguments and name of current program
 !
 nargs = command_argument_count()
 call get_command_argument(0,string)

 if (nargs <= 0) then
    print "(a)",' Usage: '//trim(string)//' [dt_min] file01.ev file02.ev ...'
    stop
 endif
 !
 ! set the time interval based on the (optional) first command line argument
 ! if not read as a sensible real number, assume dtmin = 0.
 !
 call get_command_argument(1,string)
 read(string,*,iostat=ierr) dtmin
 if (ierr /= 0) then
    dtmin = 0.
    istart = 1
 else
    print*,'Using minimum time interval = ',dtmin
    istart = 2
 endif
 !
 ! extract filenames from command line
 !
 iend = min(nargs,max_files)
 do i=istart,iend
    call get_command_argument(i,filenames(i))
 enddo
 nfiles = iend - istart + 1

 call ev2dot(nfiles,filenames(istart:iend),dtmin,output_type='kdot')

end program get_kdot
