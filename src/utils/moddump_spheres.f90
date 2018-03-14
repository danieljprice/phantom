!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  This will modify a sphere:
!  1) This option will add a radial perturbation of fac*r to a sphere
!  2) This option will add a rotational velocity to the sphere
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, options, part, physcon, prompting, timestep
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass
 use physcon,      only: pi
 use centreofmass, only: reset_centreofmass
 use timestep,     only: tmax
 use options,      only: alphamax,beta,damp
 use prompting,    only: prompt
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, parameter     :: ivr       = 1 ! option 1: radial pulsation
 integer, parameter     :: ivphi     = 2 ! option 2: rotational velocity
 integer, parameter     :: nmod_opts = 2 ! max number of options
 integer                :: i,choice
 real                   :: rad,theta,phi,vr,fac
 real                   :: omega_inner,omega_outer,domega,rad2,rad2i,rstar
 character(len=30)      :: mod_opt(nmod_opts)
 !
 !
 !--Determine modification to make
 mod_opt(:)     = 'none'
 mod_opt(ivr)   = 'Add radial pulsation' ! just make this
 mod_opt(ivphi) = 'Add rotational velocity'
 !
 write(*,"(a)") 'Velocity options: '
 do i = 1, nmod_opts
    if (trim(mod_opt(i)) /= 'none') write(*,"(a5,i2,1x,a30)") 'Case ', i, mod_opt(i)
 enddo
 !
 choice = 1
 call prompt('Enter choice',choice,1,nmod_opts)

 !
 !--Reset centre of mass
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 !
 !--Zero velocity
 vxyzu(1:3,:) = 0.0
 !
 !--Add new velocity profile
 select case(choice)
    !
 case (ivr)
    ! Author: James Wurster
    !
    !--Determine v_r as a function of radius
    fac = 0.2
    call prompt('Enter fac, where v_r = fac*r:',fac,0.)
    !
    !--Add radial velocity
!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,fac) &
!$omp private(i,rad,theta,phi,vr)
!$omp do
    do i = 1,npart
       rad   = sqrt( dot_product(xyzh(1:3,i),xyzh(1:3,i)) )
       theta = atan(xyzh(2,i)/xyzh(1,i))
       if (xyzh(1,i) < 0.0) theta = theta + pi
       phi   = acos(xyzh(3,i)/rad)
       vr    = fac*rad
       vxyzu(1,i) = vr*cos(theta)*sin(phi)
       vxyzu(2,i) = vr*sin(theta)*sin(phi)
       vxyzu(3,i) = vr*cos(phi)
    enddo
!$omp enddo
!$omp end parallel
    !--Reset .in parameters
    tmax     = 2.0*tmax
    alphamax = 0.0
    beta     = 0.0
    write(*,'(a)') "moddump: Radial velocity added"
    !
 case (ivphi)
    ! Author: Bernard Field under supervision of James Wurster
    !
    !--determine rotational profile
    omega_inner = 0.025
    omega_outer = 0.014
    call prompt('Enter angular velocity at center: ',omega_inner)
    call prompt('Enter angular velocity at surface:',omega_outer)
    domega = omega_outer - omega_inner
    !--calculate stellar radius
    rad2 = 0.
!$omp parallel default(none) &
!$omp shared(npart,xyzh) &
!$omp private(i,rad2i) &
!$omp reduction(max:rad2)
!$omp do
    do i = 1,npart
       rad2i = xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i)
       rad2  = max(rad2,rad2i)
    enddo
!$omp enddo
!$omp end parallel
    rstar = sqrt(rad2)
    !
    !--Add rotational profile
!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,rstar,domega,omega_inner) &
!$omp private(i,rad) &
!$omp reduction(max:rad2)
!$omp do
    do i=1,npart
       rad        = sqrt( xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) )
       vxyzu(1,i) = -xyzh(2,i) * (domega*rad/rstar + omega_inner)
       vxyzu(2,i) =  xyzh(1,i) * (domega*rad/rstar + omega_inner)
    enddo
!$omp enddo
!$omp end parallel
    !
    !--Reset .in parameters
    damp = 0.0
    !
    write(*,'(a)') "moddump: Rotational profile added"
    write(*,'(a)') "moddump: Supervisor's note: The output may not be in hydrostatic equilibrium, thus may be unstable"
    write(*,'(a)') "moddump:                    Further testing is required"
 end select
 !
 return
end subroutine modify_dump

end module moddump

