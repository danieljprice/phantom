!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! Handles star-disc collision simulations
!
! :References: Jankovič, T., et al., 2026: https://arxiv.org/abs/2602.02656
!
! :Owner: Taj Jankovič & Aleksej Jurca
!
!
! Note: Handles collision between a star (sink particle) and a local section of an accretion disk.
! The current version does not actually handle injection. This will be implemented in a new version.
!
! :Dependencies: boundary, eos, infile_utils, io, part, partinject,
!   physcon, units
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'selfcrossing'
 public :: init_inject,inject_particles,write_options_inject,read_options_inject,update_injected_par

 !private
 !integer, dimension(:), allocatable :: list ! might be used in the next version

contains


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
 use part,      only:ihsoft

 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 integer :: i,nbounce
 real :: star_xyz(3),star_r,star_v(3)
 real :: rel_r(3),rel_v(3),dir_r(3),min_d
 real :: dt,a,b,cquad,disc,sd,thit,t1,t2,vn, rmag

 ! Ellastic collisions between sink and gas
 if (xyzmh_ptmass(ihsoft,1) > 0.) then
   nbounce=0
   ! 'bounce' gas off sink particle
   star_xyz = xyzmh_ptmass(1:3,1)
   star_r = xyzmh_ptmass(ihsoft,1)
   star_v = vxyz_ptmass(1:3,1)
   min_d = 1e+10
   do i=1,npart
     rel_r = xyzh(1:3,i) - star_xyz
     if (min_d > norm2(rel_r)) min_d = norm2(rel_r)
     if (norm2(rel_r) < star_r) then
        dt  = dtlast
        rel_v = vxyzu(1:3,i) - star_v   ! relative velocity (star frame)
        rel_r = rel_r - rel_v*dt     !   (reconstruct particle position at t-dt)

        !!!! get the location on the stellar surface where particle crossed it !!!!!
        ! Solve |r0 + v t|^2 = Rstar^2 for t in [0,dt]
        a     = dot_product(rel_v,rel_v)
        b     = 2.*dot_product(rel_r,rel_v)
        cquad = dot_product(rel_r,rel_r) - star_r*star_r
        disc  = b*b - 4.*a*cquad

        thit = dt
        if (disc >= 0. .and. a > 0.) then
          sd = sqrt(disc)
          t1 = (-b - sd)/(2.*a)
          t2 = (-b + sd)/(2.*a)
          if (t1 >= 0. .and. t1 <= thit) thit = t1
          if (t2 >= 0. .and. t2 <= thit) thit = t2
        endif

        ! hit point and normal
        rel_r = rel_r + rel_v*thit           ! r_hit
        rmag  = norm2(rel_r)
        dir_r = rel_r / rmag         ! n-hat
        rel_r = dir_r * star_r               ! snap exactly to surface of the star

        ! reflect vel. only if moving inward (keep tangential component and flip normal component) and drift remaining time
        vn = dot_product(rel_v, dir_r)
        if (vn < 0.) then
           nbounce = nbounce + 1
           rel_v = rel_v - 2.*vn*dir_r
        endif

        ! drift remaining time - new method
        rel_r = rel_r + rel_v*(dt - thit)    ! r_new

        ! back to lab frame
        xyzh(1:3,i)  = star_xyz + rel_r
        vxyzu(1:3,i) = rel_v   + star_v
      endif

   enddo
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
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
end subroutine read_options_inject

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine

end module inject
