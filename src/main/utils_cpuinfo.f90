!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: cpuinfo
!
!  DESCRIPTION:
!  This module extracts information about the system hardware
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module cpuinfo
 implicit none

 public :: print_cpuinfo, get_cpuinfo

 private

contains

!---------------------------------------------------
!+
!  Interface routine to nicely print out
!  information regarding the cpu
!+
!---------------------------------------------------
subroutine print_cpuinfo(unit)
 integer, intent(in), optional :: unit
 character(len=60) :: cpumodel
 character(len=40) :: cachesize
 integer :: ncpu,ncpureal
 real :: cpuspeed
 integer :: lu,ierr

 if (present(unit)) then
    lu = unit
 else
    lu = 6 ! stdout
 endif

 call get_cpuinfo(ncpu,ncpureal,cpuspeed,cpumodel,cachesize,ierr)
 if (ierr==0) then
    write(lu,"(/,1x,a)") 'Running on '//trim(cpumodel)
    if (ncpureal > 10) then
       write(lu,"(1x,i3,a)",advance='no') ncpureal,' cpus @'
    else
       write(lu,"(1x,i1,a)",advance='no') ncpureal,' cpus @'
    endif
    write(lu,"(f6.2,a)",advance='no') cpuspeed,' Ghz'
    write(lu,"(a)")        ', cache size '//trim(adjustl(cachesize))
 endif

end subroutine print_cpuinfo

!---------------------------------------------------
!+
!  Actual routine that finds the cpu information
!  WARNING: This is very system dependent!!
!  Works on most Linux and Mac systems
!+
!---------------------------------------------------
subroutine get_cpuinfo(ncpu,ncpureal,cpuspeed,cpumodel,cachesize,ierr)
 real,             intent(out) :: cpuspeed
 integer,          intent(out) :: ncpu,ncpureal
 character(len=*), intent(out) :: cpumodel,cachesize
 integer,          intent(out) :: ierr
 logical :: iexist
 integer, parameter :: iunit = 608
 character(len=80) :: line
 character(len=40) :: tempfile
 real :: cachesizel2,cachesizel3
!
!--on Linux, cpu info will be located in the /proc/cpuinfo file
!  So we look in this file first
!
 ierr = 0
 ncpu = 0
 inquire(file='/proc/cpuinfo',exist=iexist)
 if (iexist) then
    open(unit=iunit,file='/proc/cpuinfo',status='old',iostat=ierr)
    do while (ierr==0)
       read (iunit,"(a)",iostat=ierr) line
       if (matchname(line,'processor' )) ncpu = ncpu + 1
       if (matchname(line,'model name')) call readval(line,cpumodel)
       if (matchname(line,'cpu MHz'   )) call readvalr(line,cpuspeed)
       if (matchname(line,'cache size')) call readval(line,cachesize)
    enddo
    cpuspeed = cpuspeed*1.e-3 ! convert cpu speed to GHz
    close(unit=iunit,iostat=ierr)
    ncpureal = ncpu
 else
!
!--On a Mac, we have to use the sysctl utility
!
    tempfile='cpuinfo.tmp'
    call system('sysctl -a hw machdep > '//trim(tempfile))
    !--check to see if this file exists
    inquire(file=tempfile,exist=iexist)
    if (iexist) then
       open(unit=iunit,file=tempfile,status='old',iostat=ierr)
       do while (ierr==0)
          read (iunit,"(a)",iostat=ierr) line
          if (matchname(line,'hw.ncpu' ))        call readvali(line,ncpu)
          if (matchname(line,'hw.physicalcpu' )) call readvali(line,ncpureal)
          if (matchname(line,'machdep.cpu.brand_string' )) call readval(line,cpumodel)
          if (matchname(line,'hw.cpufrequency' )) call readvalr(line,cpuspeed)
          !if (matchname(line,'hw.cachesize' )) call readval(line,cachesize)
          !if (matchname(line,'hw.l1dcachesize' )) call readvalr(line,cachesizel1)
          if (matchname(line,'hw.l2cachesize' ))  call readvalr(line,cachesizel2)
          if (matchname(line,'hw.l3cachesize' ))  call readvalr(line,cachesizel3)
          ! print*,trim(line)
       enddo
       !--clean up by deleting the cpuinfo file
       close(unit=iunit,status='delete',iostat=ierr)
       cpuspeed = cpuspeed*1.e-9 ! convert cpu speed to GHz
       !--print cache sizes to the string
       cachesizel2 = cachesizel2/1024.  ! in Mb
       cachesizel3 = cachesizel3/1024.  ! in Mb
       write(cachesize,"(i4,' KB (L2) ',i5,' KB (L3)')") nint(cachesizel2),nint(cachesizel3)
    else
       ierr = 1
    endif
 endif


end subroutine get_cpuinfo

!------------------------------
!+
!  Utility to match substring
!+
!------------------------------
logical function matchname(line,name)
 character(len=*), intent(in) :: line,name

 if (index(line,name) > 0) then
    matchname = .true.
 else
    matchname = .false.
 endif

end function matchname

!---------------------------------------------------
!+
!  Utility to extract value from a name : val entry
!+
!---------------------------------------------------
subroutine readval(line,val)
 character(len=*), intent(in)  :: line
 character(len=*), intent(out) :: val
 integer :: j

 j = index(line,':')
 if (j <= 0) j = index(line,'=')
 if (j > 0) then
    val = trim(adjustl(line(j+1:)))
 else
    val = ''
 endif

end subroutine readval

!--------------------------------------
!+
!  Utility to extract integer value
!  from a name : val entry
!+
!--------------------------------------
subroutine readvali(line,ival)
 character(len=*), intent(in)  :: line
 integer,          intent(out) :: ival
 character(len=20) :: val
 integer :: ierr

 call readval(line,val)
 read(val,*,iostat=ierr) ival

end subroutine readvali

!--------------------------------------
!+
!  Utility to extract real value
!  from a name : val entry
!+
!--------------------------------------
subroutine readvalr(line,rval)
 character(len=*), intent(in)  :: line
 real,             intent(out) :: rval
 character(len=20) :: val
 integer :: ierr

 call readval(line,val)
 read(val,*,iostat=ierr) rval

end subroutine readvalr

end module cpuinfo
