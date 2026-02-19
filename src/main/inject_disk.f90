!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! Handles QPE simulations
!
! :References: None
!
! :Owner: Taj Jankovič & Aleksej Jurca
!
!
! :Dependencies: boundary, eos, infile_utils, io, part, partinject,
!   physcon, units
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'selfcrossing'
 public :: init_inject,inject_particles,write_options_inject,read_options_inject,update_injected_par

 private
 integer, dimension(:), allocatable :: list

contains





subroutine ExtendList(list, min_size)

          IMPLICIT NONE

          integer :: i, isize
          integer, intent(in) :: min_size
          integer, dimension(:), allocatable, intent(inout) :: list
          integer, dimension(:), allocatable :: clist

          if(allocated(list)) then
              isize = size(list)
              if(min_size > isize) then
                allocate(clist(max(min_size, isize*2)))
                do i=1,isize
                  clist(i) = list(i)
                end do
                deallocate(list)
                call move_alloc(clist, list)
              endif
          else
              allocate(list(max(min_size, 8)))
          end if


end subroutine ExtendList

!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use part,      only:igas,hfact,massoftype,kill_particle,shuffle_part,ihsoft
 use partinject,only:add_or_update_particle
 use io,        only:master
 use sortutils, only:indexx
 use vectorutils,   only:rotatevec
 use physcon,   only:pi,solarm,years,c
 use random,    only:ran2

 real(kind=8), allocatable :: dataold(:,:),xcyl(:),ycyl(:),zcyl(:),hcyl(:),mass(:),ycyl_shift(:)

 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 character(len=512)::namefile

 real :: mdot,vinj,dz,H0,y0,disk_h
 real :: xyzi1(3),vxyz1(3),xyzi2(3),vxyz2(3)
 real :: shift,tperiod,tcross,y0up,tprev
 real :: h,u,inc,r
 integer :: k,i_part
 integer :: unit, i, n
 integer io,ninj,nrem,nbounce
 logical :: inclined_streams,vr_coll

 real :: star_xyz(3),star_r,star_v(3)
 real :: rel_r(3),rel_v(3),dir_r(3),min_d, eps
 real :: dt,a,b,cquad,disc,sd,thit,t1,t2


 ! Ellastic collisions between sink and gas
 if (xyzmh_ptmass(ihsoft,1) > 0.) then
   nbounce=0
   ! 'bounce' gas off sink particle
   star_xyz = xyzmh_ptmass(1:3,1)
   star_r = xyzmh_ptmass(ihsoft,1)
   star_v = vxyz_ptmass(1:3,1)
   min_d = 1e+10
   eps = 0.01 ! used to shift position
   do i=1,npart
     rel_r = xyzh(1:3,i) - star_xyz
     if (min_d > norm2(rel_r)) min_d = norm2(rel_r)
     if (norm2(rel_r) < star_r) then
        nbounce = nbounce + 1

        dt  = dtlast
        rel_v = vxyzu(1:3,i) - star_v   ! relative velocity (star frame)
        rel_r = rel_r - rel_v*dt     !   (reconstruct particle position at t-dt)

        !!!! get the location on the stellar surface where particle crossed it !!!!!
        ! Solve |r0 + v t|^2 = Rstar^2 for t in [0,dt]
        a     = dot_product(rel_v,rel_v)
        b     = 2.*dot_product(rel_r,rel_v)
        cquad = dot_product(rel_r,rel_r) - star_r*star_r
        disc  = b*b - 4.*a*cquad

        thit = 0.
        if (disc > 0. .and. a > 0.) then
          sd = sqrt(disc)
          t1 = (-b - sd)/(2.*a)
          t2 = (-b + sd)/(2.*a)
          thit = dt
          if (t1 >= 0. .and. t1 <= thit) thit = t1
          if (t2 >= 0. .and. t2 <= thit) thit = t2
        endif

        ! hit point and normal
        rel_r = rel_r + rel_v*thit           ! r_hit
        dir_r = rel_r / norm2(rel_r)         ! n-hat

        ! reflect v (keep tanential component and flip normal component) and drift remaining time
        rel_v = rel_v - 2.*dot_product(rel_v,dir_r)*dir_r
        rel_r = rel_r + rel_v*(dt - thit)    ! r_new

        ! back to lab frame
        xyzh(1:3,i)  = star_xyz + rel_r
        vxyzu(1:3,i) = rel_v   + star_v
      endif

   enddo
   print *,'Did bounce:',nbounce,'min_d=',min_d, 'eps=',eps
 endif
end subroutine inject_particles




!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit


end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only: fatal, error, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 character(len=30), parameter :: label = 'read_options_inject'


end subroutine read_options_inject

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine

end module inject
