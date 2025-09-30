!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module shock_capturing
!
! Storage of parameters and options related to shock capturing.
! Also utility routines needed for the Cullen & Dehnen
! shock detection switch.
!
! :References:
!   Cullen & Dehnen (2010), MNRAS 408, 669
!   Price et al. (2018), PASA 35, e031
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - alpha    : *shock viscosity parameter*
!   - alphaB   : *shock resistivity parameter*
!   - alphamax : *MAXIMUM shock viscosity parameter*
!   - alphau   : *shock conductivity parameter*
!   - beta     : *non-linear shock viscosity parameter*
!
! :Dependencies: dim, infile_utils, io
!
 use dim, only:maxalpha,maxp,nalpha,gr,isothermal,mhd,disc_viscosity
 implicit none

 ! shock viscosity, thermal conductivity, resistivity
 real, parameter :: avdecayconst = 0.1
 real, public :: alpha,alphau,beta
 real, public :: alphamax
 real, public :: alphaB
 integer, public :: ireconav

 public :: write_options_shock_capturing, read_options_shock_capturing

contains

!-----------------------------------------------------------------------
!+
!  set default shock capturing options
!+
!-----------------------------------------------------------------------
subroutine set_defaults_shock_capturing()

 if (maxalpha>0 .and. maxalpha==maxp) then
    alpha = 0.0 ! Cullen-Dehnen switch
 else
    alpha = 1.
 endif
 alphamax = 1.0
 beta = 2.0             ! beta viscosity term
 ireconav = -1

 ! shock thermal conductivity
 alphau = 1.
 if (gr) alphau = 0.1

 ! shock resistivity (MHD only)
 alphaB = 1.0

end subroutine set_defaults_shock_capturing

!-----------------------------------------------------------------------
!+
!  writes shock capturing options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_shock_capturing(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# shock capturing'
 if (maxalpha==maxp .and. nalpha > 0) then
    call write_inopt(alpha,'alpha','MINIMUM shock viscosity parameter',iunit)
    call write_inopt(alphamax,'alphamax','MAXIMUM shock viscosity parameter',iunit)
 else
    call write_inopt(alpha,'alpha','shock viscosity parameter',iunit)
 endif
 call write_inopt(beta,'beta','non-linear shock viscosity parameter',iunit)
 if (.not.isothermal) call write_inopt(alphau,'alphau','shock conductivity parameter',iunit)
 if (mhd) then
    call write_inopt(alphaB,'alphaB','shock resistivity parameter',iunit)
 endif
 if (gr) call write_inopt(ireconav,'ireconav',&
              'use reconstruction in shock viscosity (-1=off,0=no limiter,1=Van Leer)',iunit)

 if (disc_viscosity) call write_inopt(disc_viscosity,'disc_viscosity',&
                          'use cs, multiply by h/|rij| and apply to approaching/receding',iunit)

end subroutine write_options_shock_capturing

!-----------------------------------------------------------------------
!+
!  reads shock capturing options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_shock_capturing(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,warn
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 character(len=*), parameter :: label = 'read_options'
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .false. ! default to true for optional parameters

 select case(trim(name))
 case('alpha')
    read(valstring,*,iostat=ierr) alpha
    if (alpha < 0.) call fatal(label,'stupid choice of alpha')
    if (alpha > 10.) call warn(label,'very large alpha, need to change timestep',2)
    ngot = ngot + 1
 case('alphamax')
    read(valstring,*,iostat=ierr) alphamax
    if (alphamax < tiny(alphamax)) call warn(label,'alphamax = 0 means no shock viscosity',2)
    if (alphamax < 1.) call warn(label,'alphamax < 1 is dangerous if there are shocks: don''t publish crap',2)
    if (alphamax < 0. .or. alphamax > 100.) call fatal(label,'stupid value for alphamax (generally 0.0-1.0)')
 case('alphau')
    read(valstring,*,iostat=ierr) alphau
    if (alphau < 0.) call fatal(label,'stupid choice of alphau')
    if (alphau > 10.) call warn(label,'very large alphau, check timestep',3)
    ngot = ngot + 1
 case('alphaB')
    read(valstring,*,iostat=ierr) alphaB
    if (alphaB < 0.)  call warn(label,'stupid choice of alphaB',4)
    if (alphaB > 10.) call warn(label,'very large alphaB, check timestep',3)
    if (alphaB < 1.)  call warn(label,'alphaB < 1 is not recommended, please don''t publish rubbish',2)
 case('beta')
    read(valstring,*,iostat=ierr) beta
    if (beta < 0.) call fatal(label,'beta < 0')
    if (beta > 4.) call warn(label,'very high beta viscosity set')
    ngot = ngot + 1
 case('disc_viscosity')
    read(valstring,*,iostat=ierr) disc_viscosity
 case('ireconav')
    read(valstring,*,iostat=ierr) ireconav
 case default
    imatch = .false.
 end select

 if (ngot >= 2) igotall = .true.
 if (mhd .and. ngot >= 4) igotall = .true.

end subroutine read_options_shock_capturing

!-------------------------------------------------------------------------------
!+
!  function to return alphaloc from known values of d(divv)/dt and sound speed
!  for use in Cullen & Dehnen (2010) switch
!+
!-------------------------------------------------------------------------------
pure real function get_alphaloc(divvdti,spsoundi,hi,xi_limiter,alpha_min,alpha_max)
 real, intent(in) :: divvdti,spsoundi,hi,xi_limiter,alpha_min,alpha_max
 real :: source,temp

 source = 10.*hi**2*xi_limiter*max(-divvdti,0.)
 temp = spsoundi**2 !+ source
 if (temp > epsilon(temp)) then
    get_alphaloc = max(min(source/temp,alpha_max),alpha_min)
 else
    get_alphaloc = alpha_min
 endif

end function get_alphaloc

!-------------------------------------------------------------------------------
!+
!  return the xi_limiter function used in the Cullen & Dehnen switch
!  based on the spatial velocity gradients
!+
!-------------------------------------------------------------------------------
pure real function xi_limiter(dvdx)
 real(kind=4), intent(in) :: dvdx(9)
 real  :: dvxdx,dvxdy,dvxdz,dvydx,dvydy,dvydz,dvzdx,dvzdy,dvzdz
 real  :: fac,traceS,divv,curlvx,curlvy,curlvz

 dvxdx = dvdx(1)
 dvxdy = dvdx(2)
 dvxdz = dvdx(3)
 dvydx = dvdx(4)
 dvydy = dvdx(5)
 dvydz = dvdx(6)
 dvzdx = dvdx(7)
 dvzdy = dvdx(8)
 dvzdz = dvdx(9)
 divv = dvxdx + dvydy + dvzdz
 curlvx = dvzdy - dvydz
 curlvy = dvxdz - dvzdx
 curlvz = dvydx - dvxdy

 fac    = max(-divv,0.)**2
 traceS = curlvx**2 + curlvy**2 + curlvz**2
 if (fac + traceS > epsilon(0.)) then
    xi_limiter = fac/(fac + traceS)
 else
    xi_limiter = 1.
 endif

end function xi_limiter

end module shock_capturing
