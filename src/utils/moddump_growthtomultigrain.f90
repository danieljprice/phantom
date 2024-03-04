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
 use io,             only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: bins_per_dex,iu
 real                   :: smax_user
 logical                :: force_smax,file_exists
 character(len=20)      :: infile  = "bin_param.txt"

 if ((.not. use_dust) .or. (.not. use_dustgrowth)) then
    call fatal('moddump_growth2multi','Need to compile with DUST=yes and DUSTGROWTH=yes')
 endif

 !- initialise variables
 bins_per_dex = 5
 force_smax   = .false.
 smax_user    = 2.
 nmax         = 0 !- deriv called once after moddump

 !- check if param file exists, created by python script growthtomcfost.py
 inquire(file=infile, exist=file_exists)

 !- if not, switch to interactive method
 if (.not.file_exists) then
    call prompt('Set smax manually?',force_smax)
    if (force_smax) call prompt('Enter smax in cm',smax_user,0.05)
    call prompt('Enter number of bins per dex',bins_per_dex,1)
    ! write the input file so don't have to ask questions again
    open(newunit=iu,file=infile,status='new',action='write')
    write(iu,*) force_smax,smax_user,bins_per_dex
    close(iu)
 else
    !- file created by phantom/scripts/growthtomcfost.py module
    open(newunit=iu, file=infile,status='old',action='read')
    read(iu,*) force_smax,smax_user,bins_per_dex
    close(iu)
 endif

 !- delete dead or accreted particles before doing anything
 call delete_dead_or_accreted_particles(npart,npartoftype)

 !- bin dust particles into desired bins
 call bin_to_multi(bins_per_dex,force_smax,smax_user,verbose=.true.)

end subroutine modify_dump

end module moddump
