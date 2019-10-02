!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!          1D-SPH Example Programme: Supplementary Subroutines         !
!                                                                      !
!                 Copyright (c) 2015-2019 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!
! This is an example 1D SPH programme to outline the general procedure
! of how to implement NICIL into a pre-existing SPH code.
!
! These are supplementary codes required for use in nicil_ex_sph_main.
! These should be simplistic versions of subroutines already existing
! in the user's code.
!
! WARNING! This example code is intented to be an example and not used
! for actual computations.
! Comments indicated with '!**' indicate code that is required for NICIL
!
!----------------------------------------------------------------------!
module sphsup
 implicit none
 !
 ! unit number for writing to the file/screen
 integer, public, parameter :: iprint     =  6             ! Unit to write to the screen
 integer, public, parameter :: iprintdat  = 17             ! Unit to write to the .dat file
 integer, public, parameter :: iprintwarn = 18             ! Unit to write warnings to file

 private                    :: kernel,dwdh
 public                     :: dkernel
 public                     :: initialise_particles,calculate_rhoh,calculate_jcurrent
 public                     :: fatal
 private                    :: fill_xyzh
 private
 !
contains
!+
!----------------------------------------------------------------------!
!+
! Initialisation of particle positions, densities, and magnetic fields
!+
!----------------------------------------------------------------------!
subroutine initialise_particles(nmax,npart,ixmin,ixmax,umass,udist,unit_density   &
                               ,unit_Bfield,unit_energy,unit_velocity &
                               ,xyzh,vxyz,rho,omega,B,u,cs,mass)
 integer, intent(in)  :: nmax
 real,    intent(in)  :: udist,umass,unit_density,unit_Bfield,unit_energy,unit_velocity
 real,    intent(in)  :: cs
 integer, intent(out) :: npart,ixmin,ixmax
 real,    intent(out) :: xyzh(4,nmax),vxyz(3,nmax),rho(nmax),omega(nmax),B(3,nmax),u(nmax),mass
 integer, parameter   :: nboundary=2   ! number of boundary particles
 integer              :: i,j,k
 real                 :: dx1,dx0,dx,dxrat,B0,rho0,gam
 !
 gam   = 5.0/3.0
 dx    = 0.001*udist
 dx0   = dx
 dxrat = 0.998
 !--Set positions (CGS)
 npart = 0
 do k = 1,2*nboundary + 1
    do j = 1,2*nboundary + 1
       dx = dx0
       do i=1,8*nboundary
          call fill_xyzh(npart,nmax,i,j,k,dx,dxrat,dx0,xyzh)
       enddo
       if (k==nboundary + 1 .and. j==k ) ixmin = npart + 1
       dx1 = dx*dxrat**2
       do while (xyzh(1,npart)-xyzh(1,npart-1) > 0.5*dx1 )
          call fill_xyzh(npart,nmax,i,j,k,dx,dxrat,dx0,xyzh)
       enddo
       if (k==nboundary + 1 .and. j==k ) ixmax = npart
       do i=1,10*nboundary
          call fill_xyzh(npart,nmax,i,j,k,dx,dxrat,dx0,xyzh)
       enddo
    enddo
 enddo
 if (npart==nmax) then
    write(*,*) "NICIL: SPH_TEST: nmax < npart. Aborting."
    call fatal(1)
 endif
 mass  = umass/(ixmax-ixmin+1) * 8.17d-10
 !
 !--Set density, smoothing length & magnetic field (CGS)
 B         = 0.0
 B0        = 1.08d-4  ! [G]
 rho0      = 7.43d-18 ! [g cm^-3]
!$omp parallel default(none) &
!$omp shared(npart,rho,omega,B,B0,rho0,xyzh,mass,nmax,dx0) &
!$omp private(i)
!$omp do schedule(runtime)
 do i = 1,npart
    call calculate_rhoh(nmax,npart,xyzh,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),rho(i),omega(i),mass)
    B(3,i) = B0 * rho(i)/rho0
 enddo
!$omp enddo
!$omp end parallel
 !
 !--Convert to code units
 rho  = rho /unit_density
 B    = B   /unit_Bfield
 xyzh = xyzh/udist
 mass = mass/umass
 vxyz = 0.0
 !
 !--Set internal energy (in code units)
 do i = 1,npart
    u(i) = cs**2/(gam-1.0)
 enddo
 !
end subroutine initialise_particles
!----------------------------------------------------------------------!
!+
! Fills in the npart'th element of the xyzh array
!+
!----------------------------------------------------------------------!
subroutine fill_xyzh(npart,nmax,i,j,k,dx,dxrat,dx0,xyzh)
 integer, intent(inout) :: npart
 integer, intent(in)    :: nmax,i,j,k
 real,    intent(inout) :: dx
 real,    intent(in)    :: dxrat,dx0
 real,    intent(out)   :: xyzh(:,:)
 !
 npart           = min(npart + 1, nmax)
 dx              = dx*dxrat
 if (npart==1) then
    xyzh(1,npart) = 0.0
 else
    xyzh(1,npart) = xyzh(1,npart-1) + dx
 endif
 xyzh(2,npart)   = dx0*(j-1)
 xyzh(3,npart)   = dx0*(k-1)
 xyzh(4,npart)   = 3.6*dx
 !
end subroutine fill_xyzh
!----------------------------------------------------------------------!
!+
! Iteratively calculates h & rho, where h = hfac*(mass/rho)^(1/3)
!+
!----------------------------------------------------------------------!
subroutine calculate_rhoh(nmax,npart,xyzh,xi,yi,zi,hi,rhoi,omegai,mass)
 integer, parameter     :: ctrmax = 200
 real,    parameter     :: tol    = 1.0d-8
 real,    parameter     :: hfac   = 1.2
 integer, intent(in)    :: nmax,npart
 real,    intent(in)    :: xyzh(4,nmax),xi,yi,zi,mass
 real,    intent(out)   :: rhoi,omegai
 real,    intent(inout) :: hi
 integer                :: k,ctr
 real                   :: rhosum,omegatmp,dx,dy,dz,rik2,rik,q,f1,f2,hrat,hnew
 !
 !--Initialise variables
 ctr     = 0
 hrat    = 2.0*tol
 !
 do while (hrat > tol .and. ctr < ctrmax)
    ctr      = ctr + 1
    rhoi     = mass*(hfac/hi)**3
    rhosum   = 0.0
    omegatmp = 0.0
    do k = 1,npart
       dx   = xi - xyzh(1,k)
       dy   = yi - xyzh(2,k)
       dz   = zi - xyzh(3,k)
       rik2 = dx*dx + dy*dy + dz*dz
       if (rik2 < 4.0*hi**2) then
          rik      = sqrt(rik2)
          q        = rik/hi
          rhosum   = rhosum   + mass*kernel(q) /hi**3
          omegatmp = omegatmp + mass*dwdh(q,hi)/hi**4
       endif
    enddo
    omegai = 1.0 + hi/(3.0*rhoi)*omegatmp
    f1     = rhoi - rhosum
    f2     = -3.0*rhoi/hi*omegai
    hnew   = hi - f1/f2
    if (hnew.le.0.0) then
       write(*,*) "NICIL: SPH_TEST: h < 0. Aborting."
       call fatal(1)
    endif
    hrat = abs(hnew-hi)/hi
    hi   = hnew
 enddo
 if (ctr >= ctrmax) then
    write(*,*) "NICIL: SPH_TEST: rho-h failed to converge. Aborting."
    call fatal(1)
 endif
 !
end subroutine calculate_rhoh
!----------------------------------------------------------------------!
!+
! Calculates the magnetic current, J, using the symmetric operator
!+
!----------------------------------------------------------------------!
subroutine calculate_jcurrent(nmax,npart,mass,xyzh,rho,B,jcurrent,omega)
 integer, intent(in)  :: nmax,npart
 real,    intent(in)  :: mass,xyzh(4,nmax),rho(nmax),B(3,nmax),omega(nmax)
 real,    intent(out) :: jcurrent(3,nmax)
 integer              :: i,k,io
 real                 :: dx,dy,dz,rik2,rik,r1,dxr1,dyr1,dzr1,q
 real                 :: dB(3),BcrossDr(3)
 !
 jcurrent = 0.0
!$omp parallel default(none) &
!$omp shared(npart,xyzh,rho,omega,jcurrent,B,mass) &
!$omp private(i,k,dx,dy,dz,rik2,rik,r1,dxr1,dyr1,dzr1,q,dB,BcrossDr)
!$omp do schedule(runtime)
 do i = 1,npart
    do k = 1,npart
       if (i/=k) then
          dx     = xyzh(1,i) - xyzh(1,k)
          dy     = xyzh(2,i) - xyzh(2,k)
          dz     = xyzh(3,i) - xyzh(3,k)
          rik2   = dx*dx + dy*dy + dz*dz
          if (rik2 < 4.0*xyzh(4,i)**2) then
             rik  = sqrt(rik2)
             r1   = 1.0/rik
             dxr1 = dx*r1
             dyr1 = dy*r1
             dzr1 = dz*r1
             q    = rik/xyzh(4,i)
             dB   = B(:,i) - B(:,k)
             BcrossDr(1)   = dB(2)*dzr1 - dB(3)*dyr1
             BcrossDr(2)   = dB(3)*dxr1 - dB(1)*dzr1
             BcrossDr(3)   = dB(1)*dyr1 - dB(2)*dxr1
             jcurrent(:,i) = jcurrent(:,i) + BcrossDr(:)*dkernel(q)/xyzh(4,i)**4
          endif
       endif
    enddo
    jcurrent(:,i) = jcurrent(:,i)*mass/(rho(i)*omega(i))
 enddo
!$omp enddo
!$omp end parallel
 !
end subroutine calculate_jcurrent
!----------------------------------------------------------------------!
!+
! Smoothing kernel and its derivatives
!+
!----------------------------------------------------------------------!
!The M4 cubic spline softening kernel
real function kernel(q)
 real, parameter   :: fourpi = 12.5663706144d0
 real, intent(in)  :: q
 !
 if (q.lt.0.0) then
    write(*,*) 'invalid kernel'
    call fatal(1)
 else if (q.lt.1.0) then
    kernel = 4.0 - 6.0*q**2 + 3*q**3
 else if (q.lt.2.0) then
    kernel = (2.0 - q)**3
 else
    kernel = 0.0
 endif
 kernel = kernel / fourpi
 !
end function kernel
!
! The derivative of the M4 Cubic spline softening kernel
real function dkernel(q)
 real, parameter   :: fourpi = 12.5663706144d0
 real, intent(in)  :: q
 !
 if (q.lt.0.0) then
    print*, 'invalid kernel gradient'
    call fatal(1)
 else if (q.lt.1.0) then
    dkernel = q*(9.0*q - 12.0)
 else if (q.lt.2.0) then
    dkernel = -3.0*(q - 2.)**2
 else
    dkernel = 0.0
 endif
 dkernel = dkernel/fourpi
end function dkernel
!
! Derivative of the smoothing kernel with respect to smoothing length
real function dwdh(q,h)
 real, intent(in)  :: q,h
 !
 dwdh = - ( 3.0*kernel(q) + q*dkernel(q) )
end function dwdh
!----------------------------------------------------------------------!
!+
! Terminates the sph test programme, or prints a warning
!+
!----------------------------------------------------------------------!
subroutine fatal(ierr,rho,temperature)
 integer,           intent(in) :: ierr
 real,    optional, intent(in) :: rho,temperature
 integer                       :: i
 !
 if (ierr > 0) then
    ! Fatal error encountered.  Aborting.
    write(iprint,'(a)') "NICIL: SPH TEST: error encountered in NICIL.  Aborting"
    close(iprintdat)
    close(iprintwarn)
    stop
 else
    ! Warning error encountered.  Print and continue.
    if (present(rho) .and. present(temperature)) then
       write(iprintwarn,'(2(a,Es10.3),a)') "For the above error, rho = ",rho," g/cm^3 and T = ",temperature," K"
    endif
 endif
 !
end subroutine
!----------------------------------------------------------------------!
end module sphsup
