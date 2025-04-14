!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpc
!
! Read datafiles downloaded from the IAU Minor Planet Centre
! https://minorplanetcenter.net/data
!
! For example, to read all distant solar system objects, use:
! https://www.minorplanetcenter.net/iau/MPCORB/Distant.txt
!
! The relevant code required is then::
!
!   use mpc, only:read_mpc,mpc_entry
!   type(mpc_entry), allocatable :: dat(:)
!
!   call read_mpc('Distant.txt',nbodies,dat=dat)
!
! from which orbital elements can then be extracted via::
!
!   do i=1,nbodies
!      print*,dat(i)%a,dat(i)%e
!   enddo
!
! :References:
! https://minorplanetcenter.net/Extended_Files/Extended_MPCORB_Data_Format_Manual.pdf
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

 public :: mpc_entry, mpc_comet, read_mpc

 type mpc_entry
    character(len=7) :: num
    real             :: h           ! absolute magnitude, H
    real             :: g           ! slope parameter, G
    character(len=5) :: epoch
    real             :: M           ! mean anomaly (degrees)
    real             :: w           ! argument of perihelion (degrees)
    real             :: O           ! longitude of ascending node (degrees)
    real             :: inc         ! inclination (degrees)
    real             :: ecc         ! eccentricity
    real             :: n           ! mean daily motion, degrees per day
    real             :: a           ! semi-major axis
    !character(len=1) :: uncertainty
    !character(len=9) :: reference
    !integer :: num_obs  ! number of observations
    !integer :: num_opps ! number of oppositions
    !character(len=9) :: arc_years
    !real :: rms_residual
    !character(len=3) :: perturbers   ! coarse indication of perturbers
    !character(len=3) :: perturbers_2 ! precise indication of perturbers
    !character(len=10) :: computer    ! computer name
    !character(len=4) :: flags
    !character(len=27) :: designation
    !integer(kind=8) :: last_obs
 end type mpc_entry

 type mpc_comet
    character(len=4)  :: comet_num         ! periodic coment number
    character(len=1)  :: orbit_type  ! generally 'C','P' or 'D'
    character(len=7)  :: designation ! designation
    integer           :: peri_year
    character(len=5)  :: peri_month
    real              :: peri_day
    real              :: rp          ! perihelion distance (au)
    real              :: ecc         ! eccentricity
    real              :: period_short! period in years if < 100
    integer           :: period_long ! period in years if >= 100
    real              :: omega       ! argument of perihelion
    real              :: big_omega   ! longitude of ascending node
    real              :: inc         ! inclination
    character(len=9)  :: epoch
    integer           :: num_obs     ! number of observations
    character(len=1)  :: non_grav
    character(len=5)  :: first
    character(len=5)  :: last
    character(len=8)  :: perturber_info
    integer           :: nperturbing  ! number of perturbing planets
    character(len=30) :: name

 end type mpc_comet

 private

contains

!-----------------------------------------------------------------------
!+
!  read data files from the minor planet center
!+
!-----------------------------------------------------------------------
subroutine read_mpc(filename,n,nhdr,dat,dat_c)
 character(len=*),             intent(in)  :: filename
 integer,                      intent(out) :: n
 integer,                      intent(in), optional  :: nhdr
 type(mpc_entry), allocatable, intent(out), optional :: dat(:)
 type(mpc_comet), allocatable, intent(out), optional :: dat_c(:)
 integer :: iu,ierr,i,nmax

 ! open the file
 open(newunit=iu,file=filename,status='old',iostat=ierr)
 if (ierr /= 0) then
    print*,' error opening '//trim(filename)
    return
 endif

 ! skip header lines
 if (present(nhdr)) then
    do i=1,nhdr
       read(iu,*)
    enddo
 endif

 nmax = 0
 ! count lines and allocate memory
 if (present(dat) .or. present(dat_c)) then
    n = 0
    do while(ierr==0)
       read(iu,*,iostat=ierr)
       n = n + 1
    enddo
    n = n - 1
    if (present(dat)) then
       allocate(dat(n))
       nmax = n
    endif
    if (present(dat_c)) then
       allocate(dat_c(n))
       nmax = n
    endif
    ierr = 0
    rewind(iu)
 endif

 n = 0
 ! read data from file
 do while (ierr >= 0  .and. n < nmax)
    n = n + 1
    if (present(dat_c)) then
       !
       ! comets data file format
       !
       !read(iu,"(i4,a1,a7,1x,i4,1x,a5,1x,f7.4,1x,f9.6,"// &
       !         "1x,f8.6,1x,f5.2,1x,i3,1x,3(f8.4,1x),a9,1x,i5,1x,"//&
       !         "a1,1x,a5,1x,a5,1x,a8,1x,i1,a30)",iostat=ierr) dat_c(n)
       read(iu,*,iostat=ierr) dat_c(n)
    elseif (present(dat)) then
       !
       ! standard data file format
       !
       !read(iu,"(a7,1x,2(f5.2,1x),a5,1x,4(f9.5,1x),f9.7,1x,f11.8,1x,f11.7,1x)",iostat=ierr) dat(n)
       read(iu,*,iostat=ierr) dat(n)
       if (ierr > 0) then
          print*,'ERROR reading line ',n,'err=',ierr,'object=',dat(n)%num
          n = n - 1
       endif
    else
       !
       ! just count number of lines
       !
       read(iu,*,iostat=ierr)
    endif
 enddo
 n = n - 1

 ! close
 close(iu)

end subroutine read_mpc

end module mpc
