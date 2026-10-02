!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! moddump routine: store temperature from cluster (HIIfeedback) setup
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters:
!   - accrad : *accrete/remove particles outside this radius [code units]*
!
! :Dependencies: HIIRegion, deriv, infile_utils, io, moddump_utils, part,
!   prompting, ptmass
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 real :: accrad = 10.   ! accrete/remove particles outside this radius [code units]

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use HIIRegion,     only:HII_feedback,initialize_H2R,update_ionrates,iH2R
 use part,          only:xyzmh_ptmass,vxyz_ptmass,nptmass,eos_vars, &
                         delete_dead_or_accreted_particles,accrete_particles_outside_sphere,rho
 use ptmass,        only:h_acc
 use deriv,         only:get_density_global
 use io,            only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,isinkdeadhead,n,nsinkdead
 integer :: ll(nptmass)

 if (prompt_for_params) call read_interactive_moddumpfile()

 ll(:) = 0
 call accrete_particles_outside_sphere(accrad)
 isinkdeadhead = -1
 nsinkdead = 0
 iH2R=1

 ! list dead sinks
 do i=1,nptmass
    if (xyzmh_ptmass(4,i) < .001 .or. xyzmh_ptmass(5,i)==h_acc) then
       xyzmh_ptmass(:,i) = 0.
       vxyz_ptmass(:,i) = 0.
       ll(i) = isinkdeadhead
       isinkdeadhead = i
       nsinkdead = nsinkdead + 1
    endif
 enddo

 ! remove dead sinks and shuffle ptmass arrays
 n = nptmass
 do while(isinkdeadhead>0)
    if (isinkdeadhead <= n) then
       if (xyzmh_ptmass(4,n) > 0.) then
          xyzmh_ptmass(:,isinkdeadhead) = xyzmh_ptmass(:,n)
          vxyz_ptmass(:,isinkdeadhead) = vxyz_ptmass(:,n)
          isinkdeadhead = ll(isinkdeadhead)
       endif
       n = n - 1
    else
       isinkdeadhead = ll(isinkdeadhead)
    endif
    if (n < 0) call fatal('shuffle','npart < 0')
 enddo

 nptmass = nptmass-nsinkdead

 print*, "number of dead sink particles :",nsinkdead

 call get_density_global(2)
 call initialize_H2R()
 call update_ionrates(nptmass,xyzmh_ptmass,h_acc)
 call HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,rho,eos_vars)

 !call delete_dead_or_accreted_particles(npart,npartoftype)

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter radius outside which to accrete particles (code units)',accrad,0.)

end subroutine read_interactive_moddumpfile

subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 23
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(accrad,'accrad',db,errcount=nerr,min=0.)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# moddump_temp parameters'
 call write_inopt(accrad,'accrad','accrete/remove particles outside this radius [code units]',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
