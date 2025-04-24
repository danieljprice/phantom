!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_corotate
!
! Implementation of external forces needed to perform
!   simulations in a corotating frame
!   (i.e. coriolis & centrifugal forces)
!
!   Kinetic energy is 0.5*(v + Omega x r)^2
!                   = 0.5*v^2 + v.(Omega x r) + 0.5*(Omega x r)^2
!
! :References: e.g. Tong (2015) classical dynamics lecture notes
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - companion_mass    : *mass of companion*
!   - companion_xpos    : *x-position of companion*
!   - hsoft             : *softening radius of companion gravity*
!   - icompanion_grav   : *1=add companion potential, 2=add companion and primary core potential*
!   - omega_corotate    : *angular speed of corotating frame*
!   - primarycore_hsoft : *softening radius of primary core*
!   - primarycore_mass  : *mass of primary*
!   - primarycore_xpos  : *x-position of primary*
!
! :Dependencies: infile_utils, io, kernel, physcon, vectorutils
!
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, public    :: omega_corotate = 1.,hsoft = 1.,primarycore_hsoft = 1.
 real, public    :: companion_xpos = 1.,companion_mass = 1.
 real, public    :: primarycore_xpos = 1., primarycore_mass = 1.
 integer, public :: icompanion_grav = 0

 public          :: update_coriolis
 public          :: get_coriolis_force,get_centrifugal_force,get_companion_force
 public          :: write_options_corotate, read_options_corotate
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
subroutine update_coriolis(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,&
                                    vcrossomega,dkdt)
 use vectorutils, only:cross_product3D,matrixinvert3D
 use io,          only:fatal
 real, intent(in)    :: dkdt
 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(out)   :: vcrossomega(3)

 integer :: ierr
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
 A(1)    = vhalfx + dkdt*fxi
 A(2)    = vhalfy + dkdt*fyi
 A(3)    = vhalfz + dkdt*fzi

! This is the matrix from the equation for v1: [Rmat][v1] = [A]
 Rmat = reshape((/1.,             -dkdt*Omegap(3),  dkdt*Omegap(2), &
                  dkdt*Omegap(3), 1.,              -dkdt*Omegap(1), &
                 -dkdt*Omegap(2), dkdt*Omegap(1), 1.              /),(/3,3/))

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

end subroutine update_coriolis

!-----------------------------------------------------------------------
!+
!  Calculate softened gravitational force due to a companion and/or
!  primary core
!+
!-----------------------------------------------------------------------
subroutine get_companion_force(r,fextxi,fextyi,fextzi,phi)
 use kernel, only:kernel_softening
 real, intent(in)    :: r(3)
 real, intent(inout) :: fextxi,fextyi,fextzi,phi
 real                :: disp_from_companion(3),sep_from_companion,&
                        disp_from_primary(3),sep_from_primary,fmag,&
                        fmag_on_sep,phigrav,q

 disp_from_companion = (/companion_xpos,0.,0./) - r
 sep_from_companion = sqrt(dot_product(disp_from_companion,disp_from_companion))
 q = sep_from_companion / hsoft
 call kernel_softening(q**2,q,phigrav,fmag)
 fmag = fmag * companion_mass / hsoft**2
 phigrav = phigrav * companion_mass / hsoft
 fmag_on_sep = fmag / sep_from_companion
 fextxi = fextxi + disp_from_companion(1) * fmag_on_sep
 fextyi = fextyi + disp_from_companion(2) * fmag_on_sep
 fextzi = fextzi + disp_from_companion(3) * fmag_on_sep
 phi = phi + phigrav

 if (icompanion_grav == 2) then ! Get gravity from primary core
    disp_from_primary = (/primarycore_xpos,0.,0./) - r
    sep_from_primary = sqrt(dot_product(disp_from_primary,disp_from_primary))
    q = sep_from_primary / primarycore_hsoft
    call kernel_softening(q**2,q,phigrav,fmag)
    fmag = fmag * primarycore_mass / primarycore_hsoft**2
    phigrav = phigrav * primarycore_mass / primarycore_hsoft
    fmag_on_sep = fmag / sep_from_primary
    fextxi = fextxi + disp_from_primary(1) * fmag_on_sep
    fextyi = fextyi + disp_from_primary(2) * fmag_on_sep
    fextzi = fextzi + disp_from_primary(3) * fmag_on_sep
    phi = phi + phigrav
 endif

end subroutine get_companion_force

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
 call write_inopt(icompanion_grav,'icompanion_grav','1=add companion potential, 2=add companion and primary core potential',iunit)

 if ( (icompanion_grav == 1) .or. (icompanion_grav == 2) ) then
    call write_inopt(companion_mass,'companion_mass','mass of companion',iunit)
    call write_inopt(companion_xpos,'companion_xpos','x-position of companion',iunit)
    call write_inopt(hsoft,'hsoft','softening radius of companion gravity',iunit)
    if (icompanion_grav == 2) then
       call write_inopt(primarycore_mass,'primarycore_mass','mass of primary',iunit)
       call write_inopt(primarycore_xpos,'primarycore_xpos','x-position of primary',iunit)
       call write_inopt(primarycore_hsoft,'primarycore_hsoft','softening radius of primary core',iunit)
    endif
 endif
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
 case('companion_mass')
    read(valstring,*,iostat=ierr) companion_mass
    ngot = ngot + 1
 case('companion_xpos')
    read(valstring,*,iostat=ierr) companion_xpos
    ngot = ngot + 1
 case('primarycore_mass')
    read(valstring,*,iostat=ierr) primarycore_mass
    ngot = ngot + 1
 case('primarycore_xpos')
    read(valstring,*,iostat=ierr) primarycore_xpos
    ngot = ngot + 1
 case('icompanion_grav')
    read(valstring,*,iostat=ierr) icompanion_grav
    ngot = ngot + 1
 case('hsoft')
    read(valstring,*,iostat=ierr) hsoft
    ngot = ngot + 1
 case('primarycore_hsoft')
    read(valstring,*,iostat=ierr) primarycore_hsoft
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 2)

end subroutine read_options_corotate

end module extern_corotate

