!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
! :Runtime parameters: None
!
! :Dependencies: dim, growth, part, prompting, timestep
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:use_dust,use_dustgrowth
 use part,           only:delete_dead_or_accreted_particles
 use prompting,      only:prompt
 use timestep,       only:nmax
 use growth,         only:bin_to_multi
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: bins_per_dex,iforce_smax
 real                   :: smax_user
 logical                :: force_smax,file_exists
 character(len=20)      :: infile  = "bin_param.txt"

 if ((.not. use_dust) .or. (.not. use_dustgrowth)) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes AND DUSTGROWTH=yes'
    stop
 endif

 !- initialise variables
 bins_per_dex = 5
 iforce_smax  = 2
 force_smax   = .false.
 smax_user    = 2
 nmax         = 0 !- deriv called once after moddump

 !- check if param file exists, created by python script growthtomcfost.py
 inquire(file=infile, exist=file_exists)

 !- if not, switch to interactive method
 if (.not.file_exists) then
    call prompt('Set smax manually? (1=yes, 2=no)',iforce_smax,1,2)
    if (iforce_smax == 1) then
       force_smax = .true.
    else
       force_smax = .false.
    endif
    if (force_smax) call prompt('Enter smax in cm',smax_user,0.05)
    call prompt('Enter number of bins per dex',bins_per_dex,1)
 else
    !- file created by phantom/scripts/growthtomcfost.py module
    open (unit=420, file=infile)
    read(420,*) force_smax, smax_user, bins_per_dex
    close(unit=420)
 endif

 !- delete dead or accreted particles before doing anything
 call delete_dead_or_accreted_particles(npart,npartoftype)

 !- bin dust particles into desired bins
 call bin_to_multi(bins_per_dex,force_smax,smax_user,verbose=.true.)

end subroutine modify_dump

end module moddump
