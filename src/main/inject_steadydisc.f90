!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Steady disc boundary conditions
!
! When particles leave radial boundary zones they are re-injected
! into the boundary zones according to the original disc profile
! by calling the set_disc routine
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - refill_boundaries : *refill inner and outer disc?*
!
! :Dependencies: damping, eos, infile_utils, io, part, partinject, physcon,
!   setdisc
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'steadydisc'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject

 real, private :: R_ref,sig_ref
 real, private :: p_index,q_index,HoverR,M_star

 logical :: refill_boundaries = .true.
 logical, private :: read_discparams = .false.
 integer, allocatable :: listinner(:),listouter(:)
 integer, private, parameter :: maxinj = 10000
 integer, private, parameter :: inject_threshold = 100

 private

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!  Sanity check options required for this injection procedure
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use damping, only:idamp
 use io,      only:error
 integer, intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0
 allocate(listinner(maxinj),listouter(maxinj))

 if (idamp /= 3) then
    call error('init_inject','these boundary conditions require idamp=3: switching it on')
    idamp = 3
    ierr = 0
 endif

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling particle (re-)injection
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,      only:id,master,fatal
 use damping, only:r1in,r2in,r1out,r2out
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: ninner,nouter,injected

 !
 !--count the number of particles with r < R_in or r > R_out
 !
 call count_particles_outside_bounds(npart,xyzh,r1in,r2out,ninner,nouter,listinner,listouter)

 !
 !--remove those particles and place them back according to the
 !  original surface density profile
 !
 injected = 0
 if (ninner > inject_threshold) then
    if (ninner > size(listinner)) then
       call fatal('inject','too many particles outside inner boundary ',var='r',val=r1in)
    endif
    call inject_particles_in_annulus(r1in,r2in,ninner,injected,listinner)
 endif
 if (nouter > inject_threshold) then
    if (nouter > size(listouter)) then
       call fatal('inject','too many particles outside outer boundary ',var='r',val=r2out)
    endif
    call inject_particles_in_annulus(r1out,r2out,nouter,injected,listouter)
 endif

 if (id==master .and. injected > 0) print "(1x,a,es10.3,a,i0,a,g0.2,a,i0,a,g0.2)", &
      '>> t=',time,': re-injected ',ninner,' parts with R < ',r1in,' and ',nouter,' with R > ',r2out
 !
 !-- no constraint on timestep
 !
 dtinject = huge(dtinject)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine count_particles_outside_bounds(npart,xyzh,R_in,R_out,ninner,nouter,listinner,listouter)
 integer, intent(in) :: npart
 real, intent(in) :: xyzh(:,:),R_in,R_out
 integer, intent(out) :: ninner,nouter,listinner(:),listouter(:)
 integer :: i
 real :: r2

 ninner = 0
 nouter = 0
 listinner(:) = 0
 listouter(:) = 0
 do i=1,npart
    r2 = xyzh(1,i)**2 + xyzh(2,i)**2
    if (r2 < R_in**2) then
       ninner = ninner + 1
       if (ninner < size(listinner)) listinner(ninner) = i
    elseif (r2 > R_out**2) then
       nouter = nouter + 1
       if (nouter <= size(listouter)) listouter(nouter) = i
    endif
 enddo

 if (ninner == size(listinner) .or. nouter == size(listouter)) then
    print*,' WARNING: injection is limited by array sizes, increase maxinj and recompile...'
 endif

end subroutine count_particles_outside_bounds

!-----------------------------------------------------------------------
!+
!  re-inject particles in a disc annulus
!+
!-----------------------------------------------------------------------
subroutine inject_particles_in_annulus(r1,r2,ninject,injected,list)
 use eos,        only:polyk,gamma
 use io,         only:id,master
 use part,       only:igas,xyzh,vxyzu,hfact,iunknown,set_particle_type
 use setdisc,    only:set_disc
 use partinject, only:updated_particle
 use physcon,    only:pi

 real,    intent(in)    :: r1,r2
 integer, intent(inout) :: ninject,injected
 integer, intent(in) :: list(ninject)
 integer :: i,j

 real :: xyzh_inject(4,ninject)
 real :: vxyzu_inject(4,ninject)
 real :: pmass_tmp

 ! re-inject particles according to the original disc profile
 ! in the damping zones
 call set_disc(id,master  = master,       &
               npart      = ninject,      &
               rmin       = r1,           &
               rmax       = r2,           &
               rref       = R_ref,        &
               phimax     = 2.*pi,        & ! prevents CofM adjustment
               p_index    = p_index,      &
               q_index    = q_index,      &
               HoverR     = HoverR,       &
               sig_norm   = sig_ref,      &
               star_mass  = M_star,       &
               polyk      = polyk,        &
               gamma      = gamma,        &
               particle_mass = pmass_tmp, &
               hfact      = hfact,        &
               xyzh       = xyzh_inject,  &
               vxyzu      = vxyzu_inject, &
               writefile  = .false., &
               verbose = .false.)

 ! flag that particle properties have been updated
 if (ninject > 0) updated_particle = .true.

 do i=1,ninject
    j = list(i)
    xyzh(1:3,j) = xyzh_inject(1:3,i)
    xyzh(4,j) = abs(xyzh(4,j))  ! ensure not dead or accreted
    vxyzu(1:3,j) = vxyzu_inject(1:3,i)
    call set_particle_type(j,iunknown) ! this is reset in update_injected_particles
 enddo
 injected = injected + ninject

end subroutine inject_particles_in_annulus

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(refill_boundaries,'refill_boundaries','refill inner and outer disc?',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error,fatal,iunit=>imflow,lenprefix,fileprefix
 use damping,      only:r1in,r2out,r2in,r1out
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 type(inopts), allocatable     :: db(:)
 character(len=lenprefix+11)   :: filename
 integer, save :: ngot
 integer :: nerr

 imatch  = .false.
 igotall = .true.

 select case(trim(name))
 case('refill_boundaries')
    read(valstring,*,iostat=ierr) refill_boundaries
    ngot = ngot + 1
    imatch = .true.
 end select
 igotall = (ngot >= 1)

 ! read some additional options from the .discparams file
 nerr = 0
 if (.not.read_discparams) then
    filename = trim(fileprefix)//'.discparams'
    call open_db_from_file(db,filename,iunit,ierr)
    if (ierr /= 0) then
       call error('inject','could not open .discparams file')
       return
    endif
    call read_inopt(r1in,'R_in',db,errcount=nerr)
    r2in = 1.19*r1in
    call read_inopt(r2out,'R_out',db,errcount=nerr)
    r1out = 0.84*r2out
    call read_inopt(R_ref,'R_ref',db,errcount=nerr)
    call read_inopt(HoverR,'H/R_ref',db,errcount=nerr)
    call read_inopt(sig_ref,'sig_ref',db,errcount=nerr,min=tiny(0.))
    call read_inopt(M_star,'M_star',db,errcount=nerr)
    call read_inopt(p_index,'p_index',db,errcount=nerr)
    call read_inopt(q_index,'q_index',db,errcount=nerr)
    if (nerr /= 0) call fatal('inject','error with parameters in .discparams file',var='errors',ival=nerr)
    call close_db(db)
    if (ierr==0 .and. nerr==0) read_discparams = .true.
 endif

end subroutine read_options_inject

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
