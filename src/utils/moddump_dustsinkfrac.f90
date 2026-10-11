!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Changes the dustfrac around sink particles, prevents spurious shadows in RT
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - outer_radius : *radius within which to taper the dust fraction [code units]*
!
! :Dependencies: dim, moddump_utils, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter (written to / read from the prefix.moddump file)
 real :: outer_radius = 10.   ! radius within which to taper the dust fraction [code units]

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,           only:use_dust
 use part,          only:igas,ndusttypes,dustfrac
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,np_gas
 real    :: dust_to_gas,r_g

 if (.not. use_dust) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes'
    stop
 endif

 dust_to_gas = 0.01

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 !- grainsize and graindens already set if convert from one fluid to two fluid with growth
 np_gas = npartoftype(igas)

 do i=1,np_gas
    r_g = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
    if (r_g < outer_radius) then
       r_g = r_g/outer_radius
       dustfrac(1:ndusttypes,i) = dust_frac_func(r_g)*dustfrac(1:ndusttypes,i)
    endif
 enddo

 ! massoftype(igas) = massoftype(igas)*(1. + dust_to_gas)
 ! npart = np_gas

end subroutine modify_dump

real function dust_frac_func(x)
 real, intent(in) :: x
 real, parameter :: pi = 4.*atan(1.)

 dust_frac_func = sin(x*pi/2)

end function dust_frac_func

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter outer radius within which to taper dustfrac (code units)',outer_radius,0.)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(outer_radius,'outer_radius',db,errcount=nerr,min=0.)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# dustsinkfrac parameters'
 call write_inopt(outer_radius,'outer_radius','taper dust fraction within this radius [code units]',iunit)

end subroutine write_options_moddump

end module moddump
