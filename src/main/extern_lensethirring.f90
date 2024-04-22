!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_lensethirring
!
! This module contains routines relating to the computation
! of a Lense--Thirring precession torque on the gas particles
! (See Lodato & Pringle 2006 and references therein)
! CJN 22/06/11
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - blackhole_spin : *spin of central black hole (-1 to 1)*
!
! :Dependencies: infile_utils, io, physcon, units, vectorutils
!
 implicit none
 real, parameter, private :: clight = 1.0
 real, parameter, private :: bigG   = 1.0
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, public :: blackhole_spin = 1.
 real, public :: blackhole_spin_angle = 0.
 real, public :: cos_spinangle = 1., sin_spinangle = 0.

 public :: update_ltforce
 public :: get_lense_thirring_force,check_lense_thirring_settings
 public :: write_options_ltforce, read_options_ltforce
 private

contains

!----------------------------------------------
!+
!  compute the precession frequency `Omegap'
!  note this is both azimuthal AND radial
!  precession.
!  (`h' in Nelson & Papaloizou 2000)
!+
!----------------------------------------------
subroutine calc_omegap(xi,yi,zi,bh_mass,Omegap)
 real, intent(in)  :: xi,yi,zi      ! Particles position
 real, intent(in)  :: bh_mass       ! Black hole mass
 real, intent(out) :: Omegap(3)        ! Precession frequency

 real :: xsink,ysink,zsink,modr,wpdotr,m2
 real :: Jsink(3),wp(3),rsinkgas(3)

 m2 = bigG*bh_mass*bh_mass/clight

! Assuming one sink for now and no backreaction - fix later...
 Jsink(1) = m2*blackhole_spin*sin_spinangle  ! eq. (6) in Nelson/Pap (2000)
 Jsink(2) = 0.
 Jsink(3) = m2*blackhole_spin*cos_spinangle

 xsink = 0.0
 ysink = 0.0
 zsink = 0.0

! Feed in the spin and mass of the ptmass
 wp(:) = (bigG/clight**2) * Jsink(:)

 rsinkgas(1) = xi - xsink
 rsinkgas(2) = yi - ysink
 rsinkgas(3) = zi - zsink

 modr = sqrt(dot_product(rsinkgas,rsinkgas))

 wpdotr = dot_product(wp,rsinkgas)

 Omegap(:) = (2.0/modr**3)*wp(:) - (6.0*wpdotr/modr**5)*rsinkgas(:)

end subroutine calc_omegap

!---------------------------------------------------------------
!+
!  subroutine to compute the Lense-Thirring force
!  (explicitly)
!+
!---------------------------------------------------------------
subroutine get_lense_thirring_force(r,vel,bh_mass,vcrossomega)
 use vectorutils, only:cross_product3D
 real, intent(in)  :: r(3),vel(3)
 real, intent(in)  :: bh_mass
 real, intent(out) :: vcrossomega(3)
 real :: Omegap(3)

! First get Omega:
 call calc_omegap(r(1),r(2),r(3),bh_mass,Omegap)
! Compute Lense--Thirring force.
 call cross_product3D(vel,Omegap,vcrossomega)

end subroutine get_lense_thirring_force

!---------------------------------------------------------------
!+
! subroutine to read in the usual forces and
! add to them the Lense--Thirring force.
! Uses an implicit scheme to define v1.
! Reads in the vxyzuhalf from the leap frog half step : vhalfx y z
! The current forces on the particle                  : fxi,fyi,fzi
! and the timestep                                    : dt
! returning the forces plust the Lense-Thirring force.
!+
!---------------------------------------------------------------
subroutine update_ltforce(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,&
                                   vcrossomega,dkdt,xi,yi,zi,bh_mass)

 use vectorutils, only : cross_product3D,matrixinvert3D
 use io,          only : fatal,warning

 real, intent(in)    :: dkdt,xi,yi,zi,bh_mass
 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(out)   :: vcrossomega(3)

 integer :: ierr
 real :: A(3),v1(3),Omegap(3) !,v1check
 real :: Rmat(3,3),Rinv(3,3)

! Half the timestep and compute its square

! Equation we are solving is: v1 = v0 + 0.5dt*(f0 + f1_sph + v1 cross Omega)
! vhalf = v0 + 0.5*dt*f0
! fxi,fyi,fzi are the components of f1_sph
! Then rearranging the above equation to get v1 = ?

! First get Omega:
 call calc_omegap(xi,yi,zi,bh_mass,Omegap)

!----------------------------------------------------------------------
! Third attempt with matrix inversion.
!----------------------------------------------------------------------
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
    call fatal('extern_lensethirring','Error: determinant = 0 in matrix inversion')
 endif

! Comupte v1 via matrix multiplication.
 v1(:) = matmul(A,Rinv)

! Compute Lense--Thirring force.
 call cross_product3D(v1,Omegap,vcrossomega)

! Finally add to the forces and return.
 fxi = fxi + vcrossomega(1)
 fyi = fyi + vcrossomega(2)
 fzi = fzi + vcrossomega(3)

! Sanity check for v1
! v1check(:) = A(:) + dton2*vcrossomega(:)
! if (abs(v1(1)-v1check(1))>1.0d-12 .or. abs(v1(2)-v1check(2))>1.0d-12 .or. abs(v1(3)-v1check(3)) > 1.0d-12) then
!    write(*,'("v1x = ",ES12.5," vlxcheck = ",ES12.5," abs(difference) = ",ES12.5)') &
!         v1(1),v1check(1),abs(v1(1) - v1check(1))
!    write(*,'("v1y = ",ES12.5," vlycheck = ",ES12.5," abs(difference) = ",ES12.5)') &
!         v1(2),v1check(2),abs(v1(2) - v1check(2))
!    write(*,'("v1z = ",ES12.5," vlzcheck = ",ES12.5," abs(difference) = ",ES12.5)') &
!         v1(3),v1check(3),abs(v1(3) - v1check(3))
!    call fatal('extern_lensethirring','incorrect velocity calculation')
! endif

! Check that the Lense-Thirring force is only a small fraction of the total force
! f2   = fxi*fxi + fyi*fyi + fzi*fzi
! flt2 = dot_product(vcrossomega,vcrossomega)

!--warn if the L-T force is greater than 10% of the total force
! if ((flt2 > 0.01*f2) .and. (f2 > 1.e-9)) then
!    call warning('extern_lensethirring',' lense-thirring force > 10% of total force')
! endif

end subroutine update_ltforce

!---------------------------------------------------------------
!+
!  checks settings for the Lense-Thirring force
!+
!---------------------------------------------------------------
subroutine check_lense_thirring_settings(ierr,accradius1)
 use units,   only:c_is_unity,G_is_unity,get_c_code,get_G_code
 use io,      only:error,warning
 integer, intent(out) :: ierr
 real,    intent(in)  :: accradius1
 !
 !--check that c=1 in code units
 !
 ierr = 0
 if (.not.c_is_unity()) then
    ierr = ierr + 1
    call error('Lense-Thirring','c /= 1 in code units, but assumed for Lense-Thirring force',&
         var='c',val=real(get_c_code()))
 endif
 !
 !--check that G=1 in code units
 !
 if (.not.G_is_unity()) then
    ierr = ierr + 1
    call error('Lense-Thirring','G /= 1 in code units, but assumed for Lense-Thirring force',&
         var='G',val=real(get_G_code()))
 endif
 !
 !--check that the inner disc edge is >= the Schwarzschild radius
 !
 if (accradius1 < 1.0) then
    print*, 'accradius1 = ',accradius1
    print*, 'R_g        = ',1.0
    call error('Lense-Thirring',&
        'Accretion radius is < R_g',&
         var='accradius1',val=accradius1)
    ierr = ierr + 1
 elseif (accradius1 < 30.) then
    call warning('Lense-Thirring',&
        'Accretion radius is < 15 R_Schwarszchild, but Lense-Thirring not good this close to the BH',&
         var='accradius1',val=accradius1)
 endif

end subroutine check_lense_thirring_settings

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_ltforce(iunit)
 use infile_utils, only:write_inopt
 use physcon, only:pi
 integer, intent(in) :: iunit

 blackhole_spin_angle = blackhole_spin_angle*(180.0/pi)
 write(iunit,"(/,a)") '# options relating to Lense-Thirring precession'
 call write_inopt(blackhole_spin,'blackhole_spin','spin of central black hole (-1 to 1)',iunit)
 call write_inopt(blackhole_spin_angle, &
                 'blackhole_spin_angle','black hole spin angle w.r.t. x-y plane (0 to 180)',iunit)
 blackhole_spin_angle = blackhole_spin_angle*(pi/180.0)

end subroutine write_options_ltforce

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_ltforce(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use physcon, only:pi
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_ltforce'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('blackhole_spin')
    read(valstring,*,iostat=ierr) blackhole_spin
    if (blackhole_spin > 1 .or. blackhole_spin < -1.) then
       call fatal(label,'invalid spin parameter for black hole')
    endif
    ngot = ngot + 1
 case('blackhole_spin_angle')
    read(valstring,*,iostat=ierr) blackhole_spin_angle
    if (blackhole_spin_angle > 180. .or. blackhole_spin_angle < 0.) then
       call fatal(label,'invalid spin angle for black hole (should be between 0 and 180 degrees)')
    else
       blackhole_spin_angle = blackhole_spin_angle*(pi/180.0)
       sin_spinangle = sin(blackhole_spin_angle)
       cos_spinangle = cos(blackhole_spin_angle)
    endif
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_ltforce

end module extern_lensethirring
