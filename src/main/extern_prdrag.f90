!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_prdrag
!
! This module contains routines relating to the computation
! of radial and transverse radiation forces from a central sink particle
! centered at the origin.
!
! This module is intended to represent radiation drag as
! generically as possible. Drag is parametrized using "beta"
! which is defined as the ratio of radiation force to gravitational
! force. You need to specify your own routine for calculating
! beta, and include it in this module. This is the exhaustive list
! of lines you need to change:
!
! subroutine get_prdrag_spatial_force-- use beta_module, only:beta
! subroutine get_prdrag_vdependent_force-- use beta_module, only:beta
! subroutine update_prdrag_leapfrog-- use beta_module, only:beta
! subroutine write_options_prdrag-- use beta_module, only:write_options_beta
! subroutine read_options_prdrag-- use beta_module, only:read_options_beta
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - beta : *beta parameter*
!
! :Dependencies: eos, infile_utils, io, units, vectorutils
!
 use eos, only:qfacdisc

 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, private    :: k2 = 1.        ! transverse drag
 real, private    :: k0 = 1.        ! radiation pressure
 real, private    :: k1 = 1.        ! redshift
 real, private    :: beta = 0.01

 public  :: get_prdrag_spatial_force, get_prdrag_vdependent_force
 public  :: update_prdrag_leapfrog
 public  :: read_options_prdrag, write_options_prdrag

 private

contains

!------------------------------------------------
!+
!  compute the spatial part of the acceleration
!+
!------------------------------------------------
subroutine get_prdrag_spatial_force(xi,yi,zi,MStar,fextxi,fextyi,fextzi,phi)
 use units,        only:get_G_code
 real, intent(in)    :: xi,yi,zi,Mstar
 real, intent(inout) :: fextxi,fextyi,fextzi
 real, intent(out)   :: phi
 real                   :: r2,dr,dr3,betai,Mbdr3,rbetai,gcode

 gcode = get_G_code()

 r2 = xi*xi + yi*yi + zi*zi
 betai = beta
 rbetai = k0*betai
 if (r2 > epsilon(r2)) then
    dr = 1./sqrt(r2)
    dr3 = dr**3
    Mbdr3 = Mstar*gcode*(1.-rbetai)*dr3
    fextxi = fextxi - xi*Mbdr3
    fextyi = fextyi - yi*Mbdr3
    fextzi = fextzi - zi*Mbdr3
    phi    = -Mstar*gcode*dr*(1.-rbetai)
 endif

end subroutine get_prdrag_spatial_force

!-----------------------------------------------------------------------
!+
!  Routine to return velocity-dependent part of PR drag force
!+
!-----------------------------------------------------------------------
subroutine get_prdrag_vdependent_force(xyzi,vel,Mstar,fexti)
 use units, only:get_c_code,get_G_code
 real, intent(in)  :: xyzi(3), vel(3)
 real, intent(in)  :: Mstar
 real, intent(out) :: fexti(3)
 real :: rhat(3)
 real :: betai, r, r2, r3, vr, gcode, ccode

 ccode = get_c_code()
 gcode = get_G_code()

 r2     = dot_product( xyzi, xyzi )
 r      = sqrt(r2)
 r3     = r*r2

 rhat = xyzi/r
 vr = dot_product(vel, rhat)

 betai  = beta

 fexti = (-betai*Mstar*gcode/ccode)* &
            ( (vr/r3)*xyzi*k1 + vel/r2*k2 )

end subroutine get_prdrag_vdependent_force

!-----------------------------------------------------------------------
!+
!  solve for the velocity update in the leapfrog corrector step
!  i.e. v^n+1 = vhalf + 0.5*dt*f_sph + 0.5*dt*f_pr(x,v^n+1)
!+
!-----------------------------------------------------------------------
subroutine update_prdrag_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dt,xi,yi,zi,Mstar)
 use units,         only:get_c_code,get_G_code
 use io,            only:warn,fatal
 use vectorutils,   only:matrixinvert3D
 real, intent(in)    :: dt,xi,yi,zi,Mstar
 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(inout) :: fexti(3)
 integer :: ierr
 real :: dton2,r2,dr,rx,ry,rz
 real :: gcode,ccode,betai,bterm,b,vr
 real :: rhat(3),vel(3),A(3),Rmat(3,3),Rinv(3,3)
 character(len=30), parameter :: label = 'update_prdrag_leapfrog'

 ccode = get_c_code()
 gcode = get_G_code()

 ! we are solving half a timestep forwards in time
 dton2 = 0.5*dt

 r2     = xi*xi + yi*yi + zi*zi
 dr     = 1./sqrt(r2)
 rx     = xi*dr
 ry     = yi*dr
 rz     = zi*dr
 rhat   = (/rx,ry,rz/)

 ! solve for v^1 using matrix inversion of [Rmat][v1] = [A]
 A(1) = vhalfx + dton2*fxi
 A(2) = vhalfy + dton2*fyi
 A(3) = vhalfz + dton2*fzi

 betai = beta
 bterm = betai*gcode*Mstar/(ccode*r2)
 b = dton2*bterm

 ! This is the matrix from the equation for v1: [Rmat][v1] = [A]
 Rmat = reshape((/1. + b*(k2 + k1*rx*rx), b*k1*ry*rx,             b*k1*rz*rx, &
                  b*k1*rx*ry,             1. + b*(k2 + k1*ry*ry), b*k1*rz*ry, &
                  b*k1*rx*rz,             b*k1*ry*rz,             1. + b*(k2 + k1*rz*rz)/),(/3,3/))

! Get the inverse matrix
 call matrixinvert3D(Rmat,Rinv,ierr)
 if (ierr /= 0) then
    call fatal('extern_prdrag','Error: determinant = 0 in matrix inversion')
 endif

! Compute v1 via matrix multiplication.
 vel(:) = matmul(A,Rinv)

 vr = dot_product(vel,rhat)

 ! velocity dependent part of the P-R drag force (e.g. equation 142 of Klacka 1992)
 fexti(:) = -bterm*(vr*rhat*k1 + vel*k2)

 !v1check(:) = A(:) + dton2*fexti(:)  ! this should match expression for v1

 fxi = fxi + fexti(1)
 fyi = fyi + fexti(2)
 fzi = fzi + fexti(3)

end subroutine update_prdrag_leapfrog

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_prdrag(iunit)
 use infile_utils,         only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to Poynting-Robertson drag'

 call write_inopt(beta,'beta','beta parameter',iunit)
 call write_inopt(k0, 'RadiationPressure', &
                  'Radiation pressure multiplier', iunit)
 call write_inopt(k2, 'TransverseDrag', &
                  'Transverse multiplier', iunit)
 call write_inopt(k1, 'Redshift', &
                  'Redshift multiplier', iunit)

end subroutine write_options_prdrag

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_prdrag(name,valstring,imatch,igotall,ierr)
 use io, only:fatal, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_prdrag'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('beta')
    read(valstring,*,iostat=ierr) beta
    ngot = ngot + 1
 case('RadiationPressure')
    read(valstring,*,iostat=ierr) k0
    ngot = ngot + 1
 case('TransverseDrag')
    read(valstring,*,iostat=ierr) k2
    ngot = ngot + 1
 case('Redshift')
    read(valstring,*,iostat=ierr) k1
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_prdrag

end module extern_prdrag
