!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: extern_corotate
!
!  DESCRIPTION:
!   Implementation of external forces needed to perform
!   simulations in a corotating frame
!   (i.e. coriolis & centrifugal forces)
!
!   Kinetic energy is 0.5*(v + Omega x r)^2
!                   = 0.5*v^2 + v.(Omega x r) + 0.5*(Omega x r)^2
!
!  REFERENCES: e.g. Tong (2015) classical dynamics lecture notes
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    omega_corotate -- angular speed of corotating frame
!
!  DEPENDENCIES: infile_utils, io, physcon, vectorutils
!+
!--------------------------------------------------------------------------
module extern_corotate
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, public :: omega_corotate = 1.

 public :: update_coriolis_leapfrog
 public :: get_coriolis_force,get_centrifugal_force
 public :: write_options_corotate, read_options_corotate
 private

contains

!------------------------------------------------
!+
!  Compute the spatial part of the acceleration
!  This is the centrifugal force, i.e.:
!
!   f_cent = -Omega x (Omega x r)
!          = (Omega x r) x Omega
!
!  This is associated with a "potential" energy
!   Phi_cent = 0.5*(Omega x r)^2
!+
!------------------------------------------------
subroutine get_centrifugal_force(r,fextxi,fextyi,fextzi,phi)
 use vectorutils, only:cross_product3D
 real, intent(in)    :: r(3)
 real, intent(inout) :: fextxi,fextyi,fextzi
 real, intent(inout) :: phi
 real                :: Omegavec(3),Omega_cross_r(3),f_cent(3)

 call get_omega(Omegavec)
 call cross_product3D(Omegavec,r,omega_cross_r)
 call cross_product3D(Omega_cross_r,Omegavec,f_cent)

 fextxi = fextxi + f_cent(1)
 fextyi = fextyi + f_cent(2)
 fextzi = fextzi + f_cent(3)
 phi    = phi - 0.5*dot_product(Omega_cross_r,Omega_cross_r)

end subroutine get_centrifugal_force
!---------------------------------------------------------------
!+
!  Define speed and direction of rotation
!  At present direction is hard-wired to rotation in x-y plane
!+
!---------------------------------------------------------------
subroutine get_omega(omegavec)
 real, intent(out) :: omegavec(3)

 omegavec = (/0.,0.,omega_corotate/)

end subroutine get_omega

!---------------------------------------------------------------
!+
!  subroutine to compute acceleration due to the coriolis force
!
!    f_cor = -2*(Omega x v)
!
!  (this routine computes this explicitly)
!
!  also returns potential energy term v.(Omega x r) for use
!  in energy budgeting
!+
!---------------------------------------------------------------
subroutine get_coriolis_force(r,vel,f_cor,poti)
 use vectorutils, only:cross_product3D
 real, intent(in)  :: r(3),vel(3)
 real, intent(out) :: f_cor(3)
 real, intent(inout) :: poti
 real :: Omegavec(3) !,Omega_cross_r(3)

 ! define angular speed and direction of rotation
 call get_omega(Omegavec)

 ! compute Coriolis force
 call cross_product3D(vel,Omegavec,f_cor)

 f_cor(:) = 2.*f_cor(:)

 ! potential energy associated with Coriolis term: v.(Omega x r)
 !call cross_product3D(Omegavec,r,Omega_cross_r)
 !poti = poti + dot_product(vel,Omega_cross_r)

end subroutine get_coriolis_force

!---------------------------------------------------------------
!+
! subroutine to read in the usual forces and
! add to them the coriolis force.
! Uses an implicit scheme to define v1.
! Reads in the vxyzuhalf from the leap frog half step : vhalfx y z
! The current forces on the particle                  : fxi,fyi,fzi
! and the timestep                                    : dt
! returning the forces plus the Coriolis force.
!+
!---------------------------------------------------------------
subroutine update_coriolis_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,&
                                    vcrossomega,dt)
 use vectorutils, only:cross_product3D,matrixinvert3D
 use io,          only:fatal
 real, intent(in)    :: dt
 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(out)   :: vcrossomega(3)

 integer :: ierr
 real :: dton2
 real :: A(3),v1(3),Omegap(3)
 real :: Rmat(3,3),Rinv(3,3)

 ! define angular speed and direction of rotation
 call get_omega(Omegap)
 Omegap = 2.*Omegap  ! because force is 2*(v x Omega)

!--------------------------------------------------
! Perform matrix inversion to solve the equation:
!
!  v1 = v0 + 0.5dt*(f0 + f1_sph + v1 cross Omega)
!
! for v1, where
!
!   vhalf = v0 + 0.5*dt*f0
!
! and fxi,fyi,fzi are the components of f1_sph
!--------------------------------------------------
 dton2   = 0.5*dt
 A(1)    = vhalfx + dton2*fxi
 A(2)    = vhalfy + dton2*fyi
 A(3)    = vhalfz + dton2*fzi

! This is the matrix from the equation for v1: [Rmat][v1] = [A]
 Rmat = reshape((/1.,             -dton2*Omegap(3),  dton2*Omegap(2), &
                  dton2*Omegap(3), 1.,              -dton2*Omegap(1), &
                 -dton2*Omegap(2), dton2*Omegap(1), 1.              /),(/3,3/))

! Get the inverse matrix
 call matrixinvert3D(Rmat,Rinv,ierr)
 if (ierr /= 0) then
    call fatal('extern_corotate','Error: determinant = 0 in matrix inversion')
 endif

! Comupte v1 via matrix multiplication.
 v1(:) = matmul(A,Rinv)

! Compute Coriolis force.
 call cross_product3D(v1,Omegap,vcrossomega)

! Finally add to the forces and return.
 fxi = fxi + vcrossomega(1)
 fyi = fyi + vcrossomega(2)
 fzi = fzi + vcrossomega(3)

end subroutine update_coriolis_leapfrog

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_corotate(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to corotating frame'
 call write_inopt(omega_corotate,'omega_corotate','angular speed of corotating frame',iunit)

end subroutine write_options_corotate

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_corotate(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use physcon, only:pi
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_corotate'

 igotall = .false.
 imatch = .true.
 select case(trim(name))
 case('omega_corotate')
    read(valstring,*,iostat=ierr) omega_corotate
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_corotate

end module extern_corotate

