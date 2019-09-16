!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: ptmass_radiation
!
!  DESCRIPTION:
!   Implementation of radiation from sink particles
!   Contains routines to compute dust temperature assuming radiative equilibrium
!   Also routine to compute radiative acceleration based on sink particle luminosity
!
!  REFERENCES: None
!
!  OWNER: Lionel Siess
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    alpha_rad       -- fraction of the gravitational acceleration imparted to the gas
!    iget_tdust      -- method for computing dust temperature (0=none 1=T(R) 2=Lucy 3=MCFOST)
!    isink_radiation -- sink radiation pressure method (0=off,1=alpha,2=dust)
!
!  DEPENDENCIES: dim, dust_formation, infile_utils, io, kernel, part
!+
!--------------------------------------------------------------------------
module ptmass_radiation
 use part, only:iTeff,iLum
 implicit none
 integer, public  :: isink_radiation = 0
 integer, private :: iget_tdust = 0
 real,    private :: alpha_rad = 0.

 public :: get_rad_accel_from_ptmass,read_options_ptmass_radiation,write_options_ptmass_radiation
 public :: get_dust_temperature_from_ptmass

contains
!-----------------------------------------------------------------------
!+
!  Initialisation (if needed)
!+
!-----------------------------------------------------------------------
subroutine init_radiation_ptmass(ierr)
  integer, intent(out) :: ierr

  ierr = 0

end subroutine init_radiation_ptmass

!-----------------------------------------------------------------------
!+
!  compute radiative acceleration from ALL sink particles
!+
!-----------------------------------------------------------------------
subroutine get_rad_accel_from_ptmass(nptmass,npart,xyzh,xyzmh_ptmass,fext)
  use part,    only:isdead_or_accreted
#ifdef NUCLEATION
 use part,  only:nucleation
#endif
 use part,  only:dust_temp
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:)
 real,     intent(in)    :: xyzmh_ptmass(:,:)
 real,     intent(inout) :: fext(:,:)
 real                    :: dx,dy,dz,xa,ya,za,r,pmassj,plumj,pmlossj,ax,ay,az
 integer                 :: i,j

 print *,'@@@@@@@@@@@',isink_radiation,alpha_rad
 do j=1,nptmass
    pmassj  = xyzmh_ptmass(4,j)
    plumj   = xyzmh_ptmass(12,j)
    pmlossj = xyzmh_ptmass(14,j)
    xa = xyzmh_ptmass(1,j)
    ya = xyzmh_ptmass(2,j)
    za = xyzmh_ptmass(3,j)
    !$omp parallel  do default(none) &
#ifdef NUCLEATION
    !$omp shared(nucleation)&
#endif
    !$omp shared(dust_temp) &
    !$omp shared(npart,xa,ya,za,pmassj,plumj,pmlossj,xyzh,fext) &
    !$omp private(i,dx,dy,dz,ax,ay,az,r)
    do i=1,npart
      if (.not.isdead_or_accreted(xyzh(4,i))) then
        dx = xyzh(1,i) - xa
        dy = xyzh(2,i) - ya
        dz = xyzh(3,i) - za
        r = sqrt(dx**2 + dy**2 + dz**2)
#ifdef NUCLEATION
        call get_radiative_acceleration_from_star(r,dx,dy,dz,pmassj,plumj,pmlossj,ax,ay,az,K3=nucleation(5:i))
#elif BOWEN
        call get_radiative_acceleration_from_star(r,dx,dy,dz,pmassj,plumj,pmlossj,ax,ay,az,Tdust=dust_temp(i))
#else
        call get_radiative_acceleration_from_star(r,dx,dy,dz,pmassj,plumj,pmlossj,ax,ay,az)
#endif
        fext(1,i) = fext(1,i) + ax
        fext(2,i) = fext(2,i) + ay
        fext(3,i) = fext(3,i) + az
      endif
    end do
    !$omp end parallel do
 enddo

end subroutine get_rad_accel_from_ptmass

!-----------------------------------------------------------------------
!+
!  compute radiative acceleration from a SINGLE sink particle
!  based on sink particle luminosity and computed opacities / column depth
!+
!-----------------------------------------------------------------------
 subroutine get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar,Lstar,Mdot,ax,ay,az,K3,Tdust)
#ifdef NUCLEATION
  use dust_formation, only:calc_alpha_dust
#elif BOWEN
  use dust_formation, only:calc_alpha_bowen
#endif
  real, intent(in)    :: r,dx,dy,dz,Mstar,Lstar,Mdot
  real, optional, intent(in)  :: K3,Tdust
  real, intent(out)   :: ax,ay,az
  real :: fac,alpha_dust

  select case (isink_radiation)
  ! alpha wind
  case (1)
     fac = alpha_rad*Mstar/r**3
  case(2)
#ifdef NUCLEATION
     call calc_alpha_dust(Mstar,Lstar,Mdot,K3,alpha_dust)
#elif BOWEN
     call calc_alpha_bowen(Mstar,Lstar,Tdust,alpha_dust)
#endif
     fac = alpha_dust*Mstar/r**3
  case default
     fac = 0.
  end select
  ax = fac*dx
  ay = fac*dy
  az = fac*dz
end subroutine get_radiative_acceleration_from_star

!-----------------------------------------------------------------------
!+
!  compute dust temperature by ray tracing from sink particles
!  through the gas, or by simpler approximations
!+
!-----------------------------------------------------------------------
subroutine get_dust_temperature_from_ptmass(npart,xyzh,nptmass,xyzmh_ptmass,dust_temp)
 use part,    only:isdead_or_accreted,iLum,iTeff,iReff
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:)
 real,     intent(in)    :: xyzmh_ptmass(:,:)
! real,     intent(in)    :: opacity(:)
 real,     intent(out)   :: dust_temp(:)
 real                    :: r,pLumj,pTeffj,pReffj,xa,ya,za
 integer                 :: i,j

 !
 ! sanity check, return zero if no sink particles
 !
 if (nptmass < 1) then
    dust_temp = 0.
    return
 endif

 select case (iget_tdust)
 case (1)
  ! simple T(r) relation
    j = 1
    pTeffj = xyzmh_ptmass(iTeff,j)
    plumj  = xyzmh_ptmass(iLum,j)
    pReffj = xyzmh_ptmass(iReff,j) !sqrt(pLumj/(4.*pi*steboltz*utime**3/umass*pTeffj**4))
    xa = xyzmh_ptmass(1,j)
    ya = xyzmh_ptmass(2,j)
    za = xyzmh_ptmass(3,j)
    !$omp parallel  do default(none) &
    !$omp shared(npart,xa,ya,za,pReffj,pTeffj,xyzh,dust_temp) &
    !$omp private(i,r)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          r = sqrt((xyzh(1,i)-xa)**2 + (xyzh(2,i)-ya)**2 + (xyzh(3,i)-za)**2)
          dust_temp(i) = pTeffj*sqrt(pReffj/r)
       endif
    enddo
    !$omp end parallel do
 case(2)
    !call Lucy_approximation(npart,xyzh,nptmass,xyzmh_ptmass,opacity,dust_temp)
 end select

 end subroutine get_dust_temperature_from_ptmass

!-----------------------------------------------------------------------
!+
!  Interpolates a quantity computed on the discretized line of sight for all SPH particles
!  (spherical symmetry assumed)
!+
!-----------------------------------------------------------------------
subroutine interpolate_on_particles(npart, d, h, N, dmax, quantity, output)
 use part,    only:isdead_or_accreted
 integer, intent(in)  :: npart
 real,    intent(in)  :: d(npart)
 integer, intent(in)  :: N
 real,    intent(in)  :: h(N), dmax, quantity(N)
 real,    intent(out) :: output(npart)

 real :: r, dr
 integer :: i, j

 dr = dmax / N
 !should start at nwall
 do i=1,npart
    if (.not.isdead_or_accreted(h(i))) then
       r = d(i)
       j = min(int(r/dr),N-1)
       output(i) = (r-dr*j)*(quantity(j+1)-quantity(j))/dr + quantity(j)
    endif
 enddo
end subroutine interpolate_on_particles

!-----------------------------------------------------------------------
!+
!  Calculates the radiative equilibrium temperature along the half-line
!+
!-----------------------------------------------------------------------
subroutine calculate_Teq(N, R_star, Teff, tau_prime, OR2, Teq)
 integer, intent(in)  :: N
 real,    intent(in)  :: R_star, Teff, tau_prime(N), OR2(N)
 real,    intent(out) :: Teq(N)

 where(OR2  <  R_star**2)
    Teq = Teff
 elsewhere
    Teq = Teff*(0.5*(1.-sqrt(1.-(R_star**2/OR2))) + 0.75*tau_prime)**(1./4.)
 end where
end subroutine calculate_Teq

!-----------------------------------------------------------------------
!+
!  Calculates tau_prime along the discretized line of sight
!+
!-----------------------------------------------------------------------
subroutine calculate_tau_prime(N, dmax, R_star, kappa, kap, rho_over_r2, OR, tau_prime)
 integer, intent(in)  :: N
 real,    intent(in)  :: dmax, R_star, kappa, kap(N), rho_over_r2(2*N+1)
 real,    intent(out) :: OR(N), tau_prime(N)

 real :: dr, fact(N)
 integer :: i

 tau_prime(N) = 0.
 dr = dmax/N
 fact = dr/4. * R_star**2 * (kappa+kap)
 do i=1,N-1
    OR(i) = i*dr
    tau_prime(N-i) = tau_prime(N-i+1) +&
                       fact(i)*(rho_over_r2(i)+rho_over_r2(i+1)+rho_over_r2(2*N-i+1)+rho_over_r2(2*N-i+2))
 enddo
 OR(N) = dmax
end subroutine calculate_tau_prime

!-----------------------------------------------------------------------
!+
!  compute radiative equilibrium structure (Lucy approximation)
!+
!-----------------------------------------------------------------------
subroutine density(npart, position, distance2, part_h, part_mass, nfound, found,&
                     rmin, rmax, r_star, N, rho, rho_over_r2)
 use kernel, only:cnormk,wkern
 integer, intent(in)  :: npart
 real,    intent(in)  :: position(npart), distance2(npart), part_h(npart), part_mass
 integer, intent(in)  :: nfound, found(npart)
 real,    intent(in)  :: rmin, rmax, r_star
 integer :: N
 real,    intent(out) :: rho(N), rho_over_r2(N)

 real :: OH, PH2, OR(N), OR2, HR, q2, q, fact0, fact, h, h2
 real :: delta_r, rmin_o, rmin_p, rmax_p, dr, r_star2
 integer :: i, np, j, j_min, j_max

 ! Discretization of the line :
 dr = (rmax-rmin)/(N-1.)
 !     OR(i) = dr*(i-1)+rmin
 !     i = (OR(i)-rmin)/dr + 1

 ! (Very) little optimization :
 rmin_o = rmin - dr
 !     OR(i) = dr*i+rmin_o
 !     i = (OR(i)-rmin_o)/dr

 rho(:) = 0.
 fact0 = part_mass*cnormk
 do i=1,N
    OR(i) = dr*i+rmin_o
 enddo
 do i=1,nfound
    np = found(i)
    OH = position(np)
    PH2 = distance2(np)
    h = part_h(np)
    h2 = h*h
    delta_r = sqrt(4.*h2 - PH2)
    ! rmin_p and rmax_p are the positions on the line of the two intersections between the line and the interaction sphere
    rmin_p = OH-delta_r
    rmax_p = OH+delta_r
    j_min = ceiling((rmin_p-rmin_o)/dr)
    j_max = floor((rmax_p-rmin_o)/dr)
    j_min = max(1, j_min)
    j_max = min(N, j_max)
    ! Adds the contribution of the particle i on the density at all the discretized locations in the interaction sphere
    fact = fact0/h**3
    do j=j_min, j_max
       HR = OR(j) - OH
       q2 = (PH2+HR**2)/h2
       q = sqrt(q2)
       rho(j) = rho(j) + fact*wkern(q2,q)
    enddo
 enddo

 ! rho_over_r2 = 0 inside the star so that we do not divide by zero!
 r_star2 = r_star**2
 do i=1,N
    OR2 = OR(i)**2
    if (OR2  <  r_star2) then
       rho_over_r2(i) = 0.
    else
       rho_over_r2(i) = rho(i)/OR2
    endif
 enddo

end subroutine density

!-----------------------------------------------------------------------------------
!+
!  Select interesting particles (whose interaction zone intersects the line of sight)
!+
!-----------------------------------------------------------------------------------
subroutine select_particles(npart, d2_axis, part_h, nfound, found)
 integer, intent(in)  :: npart
 real,    intent(in)  :: d2_axis(npart), part_h(npart)
 integer, intent(out) :: nfound, found(npart)

 integer :: i
 real :: h

 nfound = 0
 !should start at nwall
 do i=1,npart
    h = part_h(i)
    if (d2_axis(i)  <  4.*h**2) then
       nfound = nfound + 1
       found(nfound) = i
    endif
 enddo
end subroutine select_particles

!-----------------------------------------------------------------------
!+
!  Project SPH particles on the z-axis and calculate distance to axis
!+
!-----------------------------------------------------------------------
subroutine project_on_z(npart, r, z_coord, d2_axis)
 integer, intent(in)  :: npart
 real,    intent(in)  :: r(3,npart)
 real,    intent(out) :: z_coord(npart), d2_axis(npart)

 z_coord = r(3,:)
 d2_axis = r(1,:)**2+r(2,:)**2
end subroutine project_on_z

!-----------------------------------------------------------------------
!+
!  Project SPH particles on a given axis specified by the angles theta and phi
!  and calculate distance to axis
!+
!-----------------------------------------------------------------------
!  subroutine project_on_line(npart, r, theta, phi, z_coord, d2_axis)
!    integer, intent(in)  :: npart
!    real,    intent(in)  :: r(3,npart), theta, phi
!    real,    intent(out) :: z_coord(npart), d2_axis(npart)

!    integer :: i
!    real :: u(3), OP(3), HP(3), p

 ! Direction vector
!    u(1) = sin(theta)*cos(phi)
!    u(2) = sin(theta)*sin(phi)
!    u(3) = cos(theta)

!    do i=1,npart
!      OP = r(:,i)
!      p = dot_product(OP,u)
!      z_coord(i) = p
!      HP = OP - p*u
!      d2_axis(i) = dot_product(HP,HP)
!    enddo
!  end subroutine


!-----------------------------------------------------------------------
!+
! Computes the coordinates and radial distance of SPH particles in the
! reference frame attached to the wind emitting star
!+
!-----------------------------------------------------------------------
subroutine center_star(npart, xyzh, xzy_origin, r, d, dmin, dmax)
 use dim,   only:maxp
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(4,maxp), xzy_origin(3)
 real,    intent(out) :: r(3,npart), d(npart), dmin, dmax

 integer :: i,nwall_particles

 do i=1,npart
    r(:,i) = xyzh(1:3,i)-xzy_origin
 enddo
 nwall_particles = 1
 d = sqrt(r(1,:)**2 + r(2,:)**2 + r(3,:)**2)
 dmin = minval(d(nwall_particles+1:npart))
 dmax = maxval(d(nwall_particles+1:npart))
end subroutine center_star

!-----------------------------------------------------------------------
!+
!  write options to input file
!+
!-----------------------------------------------------------------------
subroutine write_options_ptmass_radiation(iunit)
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(isink_radiation,'isink_radiation','sink radiation pressure method (0=off,1=alpha,2=dust)',iunit)
 if (isink_radiation == 1) then
    call write_inopt(alpha_rad,'alpha_rad','fraction of the gravitational acceleration imparted to the gas',iunit)
 elseif (isink_radiation == 2) then
    call write_inopt(iget_tdust,'iget_tdust','method for computing dust temperature (0=none 1=T(R) 2=Lucy 3=MCFOST)',iunit)
 endif

end subroutine write_options_ptmass_radiation

!-----------------------------------------------------------------------
!+
!  read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_ptmass_radiation(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_ptmass_radiation'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('alpha_rad')
    read(valstring,*,iostat=ierr) alpha_rad
    ngot = ngot + 1
    if (alpha_rad < 0.) call fatal(label,'invalid setting for alpha_rad (must be >= 0)')
 case('isink_radiation')
    read(valstring,*,iostat=ierr) isink_radiation
    ngot = ngot + 1
    if (isink_radiation < 0 .or. isink_radiation > 2) call fatal(label,'invalid setting for isink_radiation ([0,2])')
 case('iget_tdust')
    read(valstring,*,iostat=ierr) iget_tdust
    ngot = ngot + 1
    if (iget_tdust < 0 .or. iget_tdust > 2) call fatal(label,'invalid setting for iget_tdust ([0,2])')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)
 if (isink_radiation > 0) igotall = (ngot >= 2)

end subroutine read_options_ptmass_radiation

end module ptmass_radiation
