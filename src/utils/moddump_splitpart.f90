!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! :Runtime parameters:
!   - lattice_type : *child arrangement (0=regular lattice, 1=random)*
!   - nchild       : *number of children per particle (13 allows a regular lattice)*
!
! :Dependencies: moddump_utils, part, prompting, splitpart
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters
 integer :: nchild       = 13   ! number of children per particle
 integer :: lattice_type = 0    ! 0 for regular lattice (nchild=13 only), 1 for random

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,     only:split_all_particles
 use part,          only:delete_dead_or_accreted_particles
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 !-- the regular lattice is only possible for 13 children
 if (nchild /= 13) lattice_type = 1

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
                          nchild,lattice_type,1)

 print "(a,i0,/)",' new npart = ',npart

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter number of children per particle',nchild,2)
 if (nchild == 13) call prompt('Enter child arrangement (0=regular lattice, 1=random)',lattice_type,0,1)

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

 call read_inopt(nchild,'nchild',db,errcount=nerr,min=2)
 if (nchild == 13) call read_inopt(lattice_type,'lattice_type',db,errcount=nerr,min=0,max=1)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# split parameters'
 call write_inopt(nchild,'nchild','number of children per particle (13 allows a regular lattice)',iunit)
 if (nchild == 13) call write_inopt(lattice_type,'lattice_type','0=regular lattice, 1=random',iunit)

end subroutine write_options_moddump

end module moddump
