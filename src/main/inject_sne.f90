!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION:
!   Routine for injecting supernovae for test from Balsara & Kim (2004)
!
!  REFERENCES: Balsara & Kim (2004), ApJ 602, 1079
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, infile_utils, io, part, physcon
!+
!--------------------------------------------------------------------------
module inject
 implicit none
 character(len=*), parameter, public :: inject_type = 'supernovae'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 integer, parameter :: maxsn = 30
 real, parameter :: xyz_sn(3,maxsn) = &
   reshape((/ 7.825e-7, 1.315e-2, 7.556e-2, &
             -5.413e-2,-4.672e-2,-7.810e-2, &
             -3.211e-2, 6.793e-2, 9.346e-2, &
             -6.165e-2, 5.194e-2,-1.690e-2, &
              5.346e-3, 5.297e-2, 6.711e-2, &
              7.698e-4,-6.165e-2,-9.331e-2, &
              4.174e-2, 6.867e-2, 5.889e-2, &
              9.304e-2,-1.538e-2, 5.269e-2, &
              9.196e-3,-3.460e-2,-5.810e-2, &
              7.011e-2, 9.103e-2,-2.378e-2, &
             -7.375e-2, 4.746e-3,-2.639e-2, &
              3.653e-2, 2.470e-2,-1.745e-3, &
              7.268e-3,-3.683e-2, 8.847e-2, &
             -7.272e-2, 4.364e-2, 7.664e-2, &
              4.777e-2,-7.622e-2,-7.250e-2, &
             -1.023e-2,-9.079e-3, 6.056e-3, &
             -9.534e-3,-4.954e-2, 5.162e-2, &
             -9.092e-2,-5.223e-3, 7.374e-3, &
              9.138e-2, 5.297e-2,-5.355e-2, &
              9.409e-2,-9.499e-2, 7.615e-2, &
              7.702e-2, 8.278e-2,-8.746e-2, &
             -7.306e-2,-5.846e-2, 5.373e-2, &
              4.679e-2, 2.872e-2,-8.216e-2, &
              7.482e-2, 5.545e-2, 8.907e-2, &
              6.248e-2,-1.579e-2,-8.496e-2, &
             -9.090e-2, 2.745e-2,-5.857e-2, &
             -1.130e-2, 6.520e-2,-8.496e-2, &
             -3.186e-2, 3.858e-2, 3.877e-2, &
              4.997e-2,-8.524e-2, 5.871e-2, &
              8.455e-2,-4.098e-2,-4.438e-2/), shape=(/3,maxsn/))

 private

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
!  Main routine handling supernovae injection
!  Note that we actually only inject thermal energy, not kinetic energy
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast_u,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,      only:id,master
 use eos,     only:gamma
 use part,    only:rhoh,massoftype,igas
 real,    intent(in)    :: time, dtlast_u
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer            :: i,i_sn,ipart
 real    :: dx(3),uval,t_sn,r2
 logical :: inject_sn
 !
 ! parameters for supernovae injection, as in Balsara & Kim (2004)
 !
 real, parameter :: dt_sn = 0.00125
 real, parameter :: r_sn  = 0.005
 real, parameter :: pr_sn = 13649.6
 !
 ! determine if time is right for injection
 !
 t_sn = time/dt_sn
 i_sn = int(t_sn)
 inject_sn = abs(t_sn - i_sn) < 1.e-8
 print*,' time = ',time,' i_sn = ',i_sn,t_sn,' inject = ',inject_sn
 if (i_sn < 1 .or. i_sn > maxsn) return
 !
 !--inject sn by changing internal energy of particles
 !
 if (inject_sn) then
    if (id==master) print "(/,a,i2,a,3(es10.3,1x),a/)",&
       ' ===> BOOM! injecting supernova #',i_sn,' at (',xyz_sn(:,i_sn),')'
    print*,' gamma = ',gamma
    ipart = 0
    do i = 1,npart
       dx = xyzh(1:3,i) - xyz_sn(1:3,i_sn)
       r2 = dot_product(dx,dx)
       if (r2 < r_sn**2) then
          uval = pr_sn / ((gamma - 1.)*rhoh(xyzh(4,i),massoftype(igas)))
          print*,uval
          vxyzu(4,i) = uval
          ipart = ipart + 1
       endif
    enddo
    print*,' energy injected into ',ipart,' particles'
    print*,'--------'
 endif
 !
 !-timestep constraint
 !
 dtinject = dt_sn

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use physcon,      only: au, solarm, years
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 !call write_inopt(dt_sn,'dt_sn','time between supernovae injections',iunit)

end subroutine

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
    !case('dt_sn')
    !   read(valstring,*,iostat=ierr) dt_sn
    !   ngot = ngot + 1
    !   if (dt_sn < 0.)    call fatal(label,'invalid setting for time between supernovae (<0)')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 0)

end subroutine

end module inject
