!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for torii by RN
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, io, part, physcon
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'torus'
 public :: do_analysis, torus_analysis

 integer, parameter :: ntheta = 45
 integer, parameter :: nr = 300

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunitone)
 use io,      only:fatal
 use physcon, only:pi
 use eos,     only:get_spsound
 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunitone,numfile
 integer          :: ninbinangle(ntheta),ninbinr(nr)
 character(len=9) :: outputone,outputtwo
 integer :: i,iunittwo
 real :: rmin,rmax,radius(nr),H(nr),tvisc(nr),sigma(nr),h_smooth(nr),alpha_ss(nr)
 real :: theta(ntheta),rho_dev(ntheta),Limag(nr),ecc(nr)
 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24

 iunittwo = iunitone + 1

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(outputone,"(a4,i5.5)") 'angm',numfile
 write(*,'("Output file name is ",A)') outputone

 write(outputtwo,"(a4,i5.5)") 'azim',numfile
 write(*,'("Output file name is ",A)') outputtwo

! Setup rmin and rmax for the analysis
! Tori spreads quickly, so would need to set this manually anyway
 rmin = 0.001
 rmax = 5.0

 call torus_analysis(xyzh,vxyz,npart,pmass,time,ntheta,nr,rmin,rmax,radius,H,&
                     theta,rho_dev,ninbinangle,ninbinr,tvisc,sigma,h_smooth,alpha_ss,Limag,ecc)

 open(iunitone,file=outputone)
 write(iunitone,'("# Analysis data at t = ",es20.12)') time
 write(iunitone,"('#',7(1x,'[',i2.2,1x,a11,']',3x))") &
       1,'radius', &
       2,'sigma', &
       3,'H/R', &
       4,'<h>/H', &
       5,'alpha_ss', &
       6,'Spec l', &
       7,'ecc'

 do i = 1,nr
    if (ninbinr(i) > 0) then
       write(iunitone,'(7(es18.10,1X))') radius(i),sigma(i),(H(i)/radius(i)),h_smooth(i),alpha_ss(i),Limag(i),ecc(i)
    else
       write(iunitone,'(7(es18.10,1X))') radius(i),0.0,0.0,0.0,0.0,0.0,0.0
    endif
 enddo

 close(iunitone)

 open(iunittwo,file=outputtwo)
 write(iunittwo,'("# Analysis data at t = ",es20.12)') time
 write(iunittwo,"('#',2(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'theta', &
       2,'rho_dev'

 do i = 1,ntheta
    if (ninbinangle(i) > 0) then
       write(iunittwo,'(2(es18.10,1X))') theta(i),rho_dev(i)
    else
       write(iunittwo,'(2(es18.10,1X))') theta(i), 0.0
    endif
 enddo

 close(iunittwo)

end subroutine do_analysis

!----------------------------------------------------------------
!+
!  Torus analysis routine - so that this can be called externally
!  Adapted from a disc analysis routine written by Chris Nixon
!+
!----------------------------------------------------------------

subroutine torus_analysis(xyzh,vxyz,npart,pmass,time,ntheta,nr,rmin,rmax,radius,H,&
                     theta,rho_dev,ninbinangle,ninbinr,tvisc,sigma,h_smooth,alpha_ss,Limag,ecc)
 use physcon, only:pi
 use eos,     only:get_spsound
 use part,    only:alphaind,igas,rhoh
 real, intent(inout) :: xyzh(:,:)
 real, intent(inout) :: vxyz(:,:)
 real, intent(inout) :: pmass,time
 integer, intent(in) :: ntheta,nr,npart
 real, intent(in)    :: rmin,rmax
 real, intent(out)   :: theta(ntheta),rho_dev(ntheta)
 real                :: rho_ave,rad,dtheta,theta_i,dr,cs
 real                :: H(nr),z(npart,nr),meanz(nr),radius(nr),tvisc(nr),sigma(nr)
 real                :: alphaav,area,h_smooth(nr),alpha_ss(nr),Limag(nr),Li(3)
 real                :: ecc(nr),R_torus,Limag_mag,E,term,mu,vel(3)
 integer             :: i,ii,jj,ninbinr(nr),ninbinangle(ntheta)

 R_torus = 1.

! Set up the arrays
 dtheta = 2*pi/real(ntheta)  ! Gives segments in degrees
 do i=1,ntheta
    theta(i)= real(i-1)*dtheta
 enddo

 dr = (rmax - rmin)/real(nr-1)
 do i = 1,nr
    radius(i) = rmin + real(i-1)*dr
 enddo

! Initialise arrays to zero
 ninbinangle(:)=0
 ninbinr(:) = 0
 rho_dev(:) = 0.
 rho_ave = 0.0
 z = 0.
 meanz(:) = 0.
 H(:) = 0.
 tvisc(:) = 0.
 sigma(:) = 0.
 h_smooth(:) = 0.
 alpha_ss(:) = 0.
 Limag(:) = 0.
 Li(:) = 0.
 ecc(:) = 0.

! Loop over particles putting properties into the correct bin
 do i = 1,npart
    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       rad = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))

       if (rad > 0.95*R_torus .and. rad < 1.05*R_torus .and. abs(xyzh(3,i)) < 0.05) then
          theta_i = atan2(xyzh(2,i),xyzh(1,i))
          if (theta_i < 0.) theta_i = theta_i + 2.*pi
          ii = int(theta_i/dtheta) + 1

          if (ii > ntheta) cycle
          if (ii < 1)  cycle

          rho_dev(ii) = rho_dev(ii) + rhoh(xyzh(4,i),pmass)
          rho_ave = rho_ave + rhoh(xyzh(4,i),pmass)
          ninbinangle(ii) = ninbinangle(ii) + 1
       endif

       jj = int((rad - radius(1))/dr + 1)
       if (jj > nr) cycle
       if (jj < 1) cycle

       ninbinr(jj) = ninbinr(jj) + 1
       area = (pi*((radius(jj)+dr/2.)**2-(radius(jj)- dr/2.)**2))
       h_smooth(jj) = h_smooth(jj) + xyzh(4,i)
       sigma(jj) = sigma(jj) + pmass/area
       cs = get_spsound(2,xyzh(1:3,i),rhoh(xyzh(4,i),pmass),vxyz(:,i))
       alphaav = alphaind(1,i)
       alpha_ss(jj) = alpha_ss(jj) + alphaav
       z(ninbinr(jj),jj) = xyzh(3,i)
       if (alphaav > tiny(alphaav)) tvisc(jj) = tvisc(jj) + 10.*rad**2/(alphaav*cs*xyzh(4,i))

       ! Specific angular momentum contributions
       Li(1) = pmass*(xyzh(2,i)*vxyz(3,i)-xyzh(3,i)*vxyz(2,i))
       Li(2) = pmass*(xyzh(3,i)*vxyz(1,i)-xyzh(1,i)*vxyz(3,i))
       Li(3) = pmass*(xyzh(1,i)*vxyz(2,i)-xyzh(2,i)*vxyz(1,i))
       Limag_mag = sqrt(dot_product(Li,Li))/pmass
       mu = 10.E6
       vel = vxyz(1:3,i)
       E = 0.5*dot_product(vel,vel) - mu/rad
       term = 2.*E*Limag_mag**2/(mu**2)
       ecc(jj) = ecc(jj) + sqrt(1. + term)

       Limag(jj) = Limag(jj) + Limag_mag

    endif
 enddo

 rho_ave = rho_ave/sum(ninbinangle)

! Smooth across bins
 do i = 1,ntheta
    if (ninbinangle(i) > 0) then
       rho_dev(i) = (rho_dev(i)/ninbinangle(i) - rho_ave)
    endif
 enddo

 do i = 1,nr
    if (ninbinr(i) > 1) then
       meanz(i) = sum(z(1:ninbinr(i),i))/real(ninbinr(i))
       H(i) = sqrt(sum(((z(1:ninbinr(i),i) - meanz(i))**2)/real(ninbinr(i)-1)))
       tvisc(i) = tvisc(i)/ninbinr(i)
       h_smooth(i) = h_smooth(i)/(ninbinr(i)*H(i))
       alpha_ss(i) = 0.1*alpha_ss(i)*H(i)/(h_smooth(i)*ninbinr(i))
       Limag(i) = Limag(i)/real(ninbinr(i))
       ecc(i) = ecc(i)/real(ninbinr(i))
    else
       h_smooth(i) = 0.
       H(i) = 0.
       alpha_ss(i) = 0.
    endif
    !print*,radius(i),(H(i)/radius(i)),tvisc(i)
 enddo

end subroutine torus_analysis

end module analysis

