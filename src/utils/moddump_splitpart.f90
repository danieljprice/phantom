!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! split every particle in a dump into nchild children
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: injectutils, io, part, splitpart, systemutils
!
 implicit none
 integer :: nchild = 13
 integer :: lattice_type = 0 ! 0 for lattice, 1 for random
 integer :: ires = 1         ! use 12 particles per sphere
 character(len=*), parameter :: moddump_flags = '--nchild=13 --lattice_type=0 [0=lattice,1=random]'

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,    only:split_all_particles
 use io,           only:fatal,error
 use injectutils,  only:get_parts_per_sphere
 use part,         only:delete_dead_or_accreted_particles
 use systemutils,  only:get_command_option
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ierr

 ierr = 0

 !--give nchild as command line flag
 nchild = int(get_command_option('nchild',default=nchild))
 if (nchild < 12) lattice_type = 1
 if (nchild < 2) stop 'error nchild cannot be < 2'

 !--this can be overridden by choice of lattice, if regular is requested
 lattice_type = int(get_command_option('lattice_type',default=lattice_type))

 !-- if using the regular grid, set nchild to get desired resolution
 if (lattice_type == 0) then
    nchild = get_parts_per_sphere(ires) + 1
 endif

 !-- don't split accreted particles
 call delete_dead_or_accreted_particles(npart,npartoftype)

 ! Split 'em!
 print "(/,a,i0,a)", ' >>> splitting all particles into ',nchild,' children <<<'

 if (lattice_type==0) then
    print "(a,/)", ' >>> placing children on regular lattice <<<'
 else
    print "(a,/)", ' >>> placing children using random arrangement <<<'
 endif
 call split_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                                nchild,lattice_type,ires)

 print "(a,i0,/)",' new npart = ',npart

end subroutine modify_dump

end module moddump
