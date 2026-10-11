!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! transforms dustgrowth dump into multigrain dump for mcfost usage
!
! :References: None
!
! :Owner: Arnaud Vericel
!
! :Runtime parameters:
!   - bins_per_dex : *number of bins per dex*
!   - force_smax   : *set the maximum grain size manually*
!   - smax_user    : *maximum grain size [cm] (used if force_smax=T)*
!
! :Dependencies: dim, growth, io_control,
!   moddump_utils, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer :: bins_per_dex = 5        ! number of bins per dex
 real    :: smax_user    = 2.       ! maximum grain size in cm (used if force_smax)
 logical :: force_smax   = .false.  ! set the maximum grain size manually

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:use_dust,use_dustgrowth
 use part,           only:delete_dead_or_accreted_particles
 use io_control,     only:nmax
 use growth,         only:bin_to_multi
 use moddump_utils,  only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 if ((.not. use_dust) .or. (.not. use_dustgrowth)) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes AND DUSTGROWTH=yes'
    stop
 endif

 nmax = 0 !- deriv called once after moddump

 !- read options from the .moddump file (e.g. written by scripts/growthtomcfost.py), or prompt
 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 !- delete dead or accreted particles before doing anything
 call delete_dead_or_accreted_particles(npart,npartoftype)

 !- bin dust particles into desired bins
 call bin_to_multi(bins_per_dex,force_smax,smax_user,verbose=.true.)

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Set smax manually?',force_smax)
 if (force_smax) call prompt('Enter smax in cm',smax_user,0.05)
 call prompt('Enter number of bins per dex',bins_per_dex,1)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(force_smax,'force_smax','set the maximum grain size manually',iunit)
 call write_inopt(smax_user,'smax_user','maximum grain size [cm] (used if force_smax=T)',iunit)
 call write_inopt(bins_per_dex,'bins_per_dex','number of bins per dex',iunit)

end subroutine write_options_moddump

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(force_smax,'force_smax',db,errcount=nerr)
 call read_inopt(smax_user,'smax_user',db,min=0.,errcount=nerr)
 call read_inopt(bins_per_dex,'bins_per_dex',db,min=1,errcount=nerr)

end subroutine read_options_moddump

end module moddump
