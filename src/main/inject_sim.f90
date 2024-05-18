!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles particle injections from another simulations (for TDE outflow only currently)
!
! :References: None
!
! :Owner: Fitz) Hu
!
! :Runtime parameters:
!   - r_inject : *radius to inject tde outflow (in cm)*
!
! :Dependencies: dump_utils, fileutils, infile_utils, io, part, partinject,
!   readwrite_dumps_common, readwrite_dumps_fortran, timestep, units
!
 use fileutils, only:getnextfilename

 implicit none
 character(len=*), parameter, public :: inject_type = 'sim'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject, &
           set_default_options_inject,update_injected_par
 private :: read_injected_par
!
!--runtime settings for this module
!

! global variables

 character(len=120) :: start_dump,final_dump,pre_dump,next_dump
 integer :: npart_sim
 real    :: r_inject,r_inject_cgs=-1,next_time!,e_inject
 real, allocatable :: xyzh_pre(:,:),xyzh_next(:,:),vxyzu_next(:,:),pxyzu_next(:,:)
 logical, allocatable :: injected(:)

 character(len=*), parameter :: label = 'inject_tdeoutflow'
 character(len=*), parameter :: injected_filename = 'injected_par'

contains

!-----------------------------------------------------------------------
!+
!  Initialize -- find the start dump to inject
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,        only:error
 use timestep,  only:time
 use units,     only:udist

 integer, intent(out) :: ierr
 integer, parameter   :: max_niter=5000, idisk=23
 integer :: niter

 !
 !--find the tde dump at the right time
 !
 next_time = -1.
 next_dump = getnextfilename(start_dump)
 call get_dump_time_npart(trim(next_dump),next_time,ierr,npart_out=npart_sim)
 ierr = 0
 niter = 0

 do while (next_time < time .and. niter < max_niter)
    niter = niter + 1
    pre_dump = next_dump
    next_dump = getnextfilename(next_dump)
    call get_dump_time_npart(trim(next_dump),next_time,ierr)
    if (ierr /= 0) then
       ierr = 0
       call error('inject','error reading time and npart from '//trim(next_dump))
       cycle
    endif
 enddo
 start_dump = next_dump

 write(*,'(a,1x,es10.2)') ' Start read sims and inject particle from '//trim(next_dump)//' at t =',next_time

 r_inject = r_inject_cgs/udist ! to code unit
 allocate(xyzh_pre(4,npart_sim),xyzh_next(4,npart_sim),vxyzu_next(4,npart_sim),pxyzu_next(4,npart_sim),injected(npart_sim))
 xyzh_pre = 0.
 injected = .false.
 call read_injected_par()
 !e_inject = -1./r_inject

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer           :: ierr
 real :: tfac

 !
 !--inject particles only if time has reached
 !
 tfac = 1.
 if (time >= next_time) then
    ! read next dump
    call read_dump(next_dump,xyzh_next,ierr,vxyzu_dump=vxyzu_next,pxyzu_dump=pxyzu_next)

    npart_old = npart
    call inject_required_part_tde(npart,npartoftype,xyzh,vxyzu,xyzh_pre,xyzh_next,vxyzu_next,pxyzu_next)

    ! copy to pre for next injection use
    pre_dump = next_dump
    xyzh_pre = xyzh_next

    call find_next_dump(next_dump,next_time,ierr)
    start_dump = next_dump

    write(*,'(i10,1x,a27,1x,a)') npart-npart_old, 'particles are injected from', trim(pre_dump)

    if (pre_dump == final_dump) then
       write(*,'(a)') ' Reach the final dumpfile. Stop injecting ...'
       next_time = huge(0.)
    endif

    tfac = 1.d-40 ! set a tiny timestep so the code has time to adjust for timestep
 endif

 ! update time to next inject
 dtinject = tfac*(next_time - time)
end subroutine inject_particles

subroutine read_dump(filename,xyzh_dump,ierr,vxyzu_dump,pxyzu_dump)
 use dump_utils, only: read_array_from_file
 character(len=*), intent(in) :: filename
 real, intent(out) :: xyzh_dump(:,:)
 integer, intent(out) :: ierr
 real, intent(out), optional :: vxyzu_dump(:,:),pxyzu_dump(:,:)
 integer, parameter :: iunit = 578
 real(kind=4) :: h(npart_sim)

 !
 !--read xyzh
 !
 call read_array_from_file(iunit,filename,'x',xyzh_dump(1,:),ierr,iprint_in=.false.)
 call read_array_from_file(iunit,filename,'y',xyzh_dump(2,:),ierr,iprint_in=.false.)
 call read_array_from_file(iunit,filename,'z',xyzh_dump(3,:),ierr,iprint_in=.false.)
 call read_array_from_file(iunit,filename,'h',h,ierr,iprint_in=.false.)
 xyzh_dump(4,:) = h

 !
 !--read vxyzu
 !
 if (present(vxyzu_dump)) then
    call read_array_from_file(iunit,filename,'vx',vxyzu_dump(1,:),ierr,iprint_in=.false.)
    call read_array_from_file(iunit,filename,'vy',vxyzu_dump(2,:),ierr,iprint_in=.false.)
    call read_array_from_file(iunit,filename,'vz',vxyzu_dump(3,:),ierr,iprint_in=.false.)
    call read_array_from_file(iunit,filename,'u',vxyzu_dump(4,:),ierr,iprint_in=.false.)
 endif

 !
 !--read vxyzu
 !
 if (present(pxyzu_dump)) then
    call read_array_from_file(iunit,filename,'px',pxyzu_dump(1,:),ierr,iprint_in=.false.)
    call read_array_from_file(iunit,filename,'py',pxyzu_dump(2,:),ierr,iprint_in=.false.)
    call read_array_from_file(iunit,filename,'pz',pxyzu_dump(3,:),ierr,iprint_in=.false.)
    call read_array_from_file(iunit,filename,'entropy',pxyzu_dump(4,:),ierr,iprint_in=.false.)
 endif

end subroutine read_dump

subroutine get_dump_time_npart(filename,time,ierr,npart_out)
 use io,                      only:iprint,id,nprocs
 use dump_utils,              only:dump_h,open_dumpfile_r,read_header,free_header
 use part,                    only:maxtypes
 use readwrite_dumps_fortran, only:unfill_header
 use readwrite_dumps_common,  only:get_options_from_fileid

 character(len=*), intent(in)   :: filename
 real, intent(out)              :: time
 integer, intent(out)           :: ierr
 integer, intent(out), optional :: npart_out
 integer, parameter :: idisk=389
 character(len=120) :: fileid
 logical :: tagged,phantomdump,smalldump,use_dustfrac
 type(dump_h) :: hdr
 integer(kind=8) :: nparttot
 integer :: nblocks,npartoftype(maxtypes),npart
 real :: hfactfile,alphafile

 call open_dumpfile_r(idisk,filename,fileid,ierr)
 call get_options_from_fileid(fileid,tagged,phantomdump,smalldump,use_dustfrac,ierr)
 call read_header(idisk,hdr,ierr,tagged=tagged)
 call unfill_header(hdr,phantomdump,tagged,nparttot, &
                    nblocks,npart,npartoftype, &
                    time,hfactfile,alphafile,iprint,id,nprocs,ierr)
 call free_header(hdr,ierr)
 close(idisk)

 if (present(npart_out)) npart_out = npart

end subroutine get_dump_time_npart

subroutine find_next_dump(next_dump,next_time,ierr)
 character(len=*), intent(inout) :: next_dump
 real, intent(out) :: next_time
 integer, intent(out) :: ierr

 next_dump = getnextfilename(next_dump)
 call get_dump_time_npart(next_dump,next_time,ierr)

end subroutine find_next_dump

subroutine inject_required_part_tde(npart,npartoftype,xyzh,vxyzu,xyzh_pre,xyzh_next,vxyzu_next,pxyzu_next)
 use part,       only:igas,pxyzu,isdead_or_accreted
 use partinject, only:add_or_update_particle
 integer, intent(inout) :: npart, npartoftype(:)
 real, intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real, intent(in) :: xyzh_pre(:,:), xyzh_next(:,:), vxyzu_next(:,:), pxyzu_next(:,:)
 integer :: i,partid
 real :: r_next,r_pre,vr_next!,e_next

 !
 !--check all the particles
 !
 do i=1,npart_sim
    if (.not. isdead_or_accreted(xyzh_next(4,i)) .and. .not. injected(i)) then
       r_next = sqrt(dot_product(xyzh_next(1:3,i),xyzh_next(1:3,i)))
       r_pre = sqrt(dot_product(xyzh_pre(1:3,i),xyzh_pre(1:3,i)))
       vr_next = (dot_product(xyzh_next(1:3,i),vxyzu_next(1:3,i)))/r_next
       !e_next = 0.5*vr_next**2 - 1./r_next

       if (r_next > r_inject .and. r_pre < r_inject .and. vr_next > 0.) then! .and. e_next > e_inject) then
          ! inject particle by copy the data into position
          partid = npart+1
          call add_or_update_particle(igas,xyzh_next(1:3,i),vxyzu_next(1:3,i),xyzh_next(4,i), &
                                    vxyzu_next(4,i),partid,npart,npartoftype,xyzh,vxyzu)
          pxyzu(:,partid) = pxyzu_next(:,i)
          injected(i) = .true.
       endif
    endif
 enddo

end subroutine inject_required_part_tde

subroutine read_injected_par()
 use io, only:fatal,warning
 integer, parameter :: iunit=242
 logical :: iexist
 integer :: nread,i

 inquire(file=trim(injected_filename),exist=iexist)

 if (iexist) then
    open(iunit,file=trim(injected_filename),status='old')
    read(iunit,*) nread

    ! check if npart in file is the same as npart_sim
    if (nread /= npart_sim) call fatal('inject_sim','npart in '//trim(injected_filename)// &
                                              ' does not match npart_sim')

    do i=1,nread
       read(iunit,*) injected(i)
    enddo
    close(iunit)
 else
    call warning('inject_sim',trim(injected_filename)//' not found, assume no particles are injected')
    injected = .false.
 endif

end subroutine

subroutine update_injected_par()
 use io, only:error
 integer, parameter :: iunit=284
 logical :: iexist
 integer :: i

 if (allocated(injected)) then
    inquire(file=trim(injected_filename),exist=iexist)
    if (iexist) then
       open(iunit,file=trim(injected_filename),status='replace')
    else
       open(iunit,file=trim(injected_filename),status='new')
    endif

    write(iunit,*) npart_sim
    do i=1,npart_sim
       write(iunit,*) injected(i)
    enddo
    close(iunit)
 endif
end subroutine

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit
 character(len=10), parameter :: start_dump_default = 'dump_00000', &
                                 final_dump_default = 'dump_02000'
 real, parameter :: r_inject_default = 5.e14

 ! write something meaningful in infile
 if (r_inject_cgs  <=  0.) then
    start_dump = start_dump_default
    r_inject_cgs = r_inject_default
    final_dump = final_dump_default
 endif

 write(iunit,"(/,a)") '# options controlling particle injection'
 call write_inopt("'"//trim(start_dump)//"'",'start_dump', &
                  'dumpfile to start for injection (with relative path if in other direc)',iunit)
 call write_inopt(r_inject_cgs,'r_inject','radius to inject tde outflow (in cm)',iunit)
 call write_inopt("'"//trim(final_dump)//"'",'final_dump', &
                  'stop injection after this dump (with relative path if in other direc)',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr
 character(len=30), parameter :: label = 'read_options_inject'
 integer, save :: ngot

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('start_dump')
    read(valstring,*,iostat=ierr) start_dump
    ngot = ngot + 1
 case('r_inject')
    read(valstring,*,iostat=ierr) r_inject_cgs
    ngot = ngot + 1
    if (r_inject_cgs < 0.) call fatal(label,'invalid setting for r_inject (<0)')
 case('final_dump')
    read(valstring,*,iostat=ierr) final_dump
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 3)

end subroutine read_options_inject

subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

end subroutine set_default_options_inject

end module inject
