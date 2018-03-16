!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, part, physcon, setup_params, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,         only:Bxyz,rhoh,hrho,mhd
 use part,         only:iphase,iamtype,igas,maxphase,get_partinfo
 use dim,          only:maxp,maxvxyzu
 use setup_params, only:ihavesetupB
 use physcon,      only:pi,au,solarm,solarr
 use units,        only:set_units
 integer,           intent(in)  :: id
 integer,           intent(out) :: npart
 integer,           intent(out) :: npartoftype(:)
 real,              intent(out) :: xyzh(:,:)
 real,              intent(out) :: polyk,gamma,hfact
 real,              intent(out) :: vxyzu(:,:)
 real,              intent(out) :: massoftype(:)
 real,              intent(in)  :: time
 character(len=20), intent(in)  :: fileprefix
 integer :: i,nrings,nlayers,iring,iz,ipart,npartphi,ii,iamtypei
 real :: Rtorus,dfac,Mstar,Mtorus,zmax,deltaz,bigG
 real :: massp,r_in,r_out,deltar,polyn,sumA,np
 real :: ri,zi,rhofac,deltaphi,densi,pri
 real :: deltartemp,denstemp,rtemp,deltar0,dens0
 real :: omegai,v2onr,rcyl2,rcyl,rsph,rhosum,pmassi,pmassii
 real :: beta,Bzi,dbeta,densmax,densmin
 real, parameter :: dndim = 1./3.
 logical :: iactivei,iamdusti

! call set_units(dist=au,mass=solarm,G=1.d0)
 bigG = 1.d0
 call set_units(dist=100.*solarr,mass=10.E6*solarm,G=bigG)

 Rtorus = 1.0
 dfac = 1.01
 Mstar = 1.0
 Mtorus = 4.e-4
 densmax = 0.25e-8 !1.e-3 ! set maximum density in torus
 hfact = 1.2
 gamma = 5./3.

 print*,'Setup for beer-related torus thing (Owen, Laure)'
 print*,'maxvxyzu value: ',maxvxyzu
 print*,'particles: ',npart

!
!--integrate to get total mass (in order to set A in P = A*rho^gamma)
!  (Owen was particularly fussy about using A not K)
!
 r_in = 0.01
 r_out = 2.5
 nrings = 200
 nlayers = 40
 deltar = (r_out - r_in)/real(nrings-1)
 zmax = 1.5
 deltaz = 2.*zmax/real(nlayers-1)
 polyn = 1./(gamma-1.)
! sumA = 0.
! do iring=1,nrings
!    ri = r_in + (iring-1)*deltar
!    zi = 0.
!    sumA = sumA + 2.*pi*ri*deltar*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
! enddo
!--rearrange to get A (which I call polyk)
! polyk = ((sumA/Mtorus)**(1./polyn))/(polyn + 1.)
!
!--alternatively determine polyk from maximum density
!
 polyk = bigG*Mstar/(Rtorus*(polyn + 1.)*densmax**(gamma-1.))*(dfac - 1.)/(2.*dfac)
 print*,' polyk = ',polyk
!
!--rhofac is the factor which multiplies rhofunc to give rho
!  this is 1/B in Owen lingo
!
 rhofac = 1./(polyk*(polyn+1.))**polyn
!
!--work out total mass in torus and set particle mass
!
 Mtorus = 0.
 do iring=1,nrings
    ri = r_in + (iring-1)*deltar
    do iz = 1,nlayers
       zi = (iz - 1)*deltaz - zmax
       Mtorus = Mtorus + 2.*pi*ri*deltar*deltaz*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
!       print*,'ri = ',ri,'zi = ',zi, ' rho = ',rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    enddo
 enddo
 Mtorus = Mtorus*rhofac
 print*,'torus mass (by integration) = ',Mtorus
 if (Mtorus > Mstar) stop 'error Mtorus > Mstar'

 print*,'tentative npart = ',maxp
 np = 0.5*maxp      ! Decrease this factor to allow for higher d values
 massp = Mtorus/real(np/1.0)
 massoftype(1) = massp
 densmin = 1.e-20 !densmax/100.

!
!--setup first ring at r=Rtorus
!
 ipart = 0
 ri = Rtorus
 zi = 0.
 densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
!--calculate delta r and delta phi from rho
 deltar = (massp/densi)**dndim
 deltaphi = deltar/ri
 npartphi = int(2.*pi/deltaphi)
!--correct deltaphi to integer number of particles
 deltaphi = 2.*pi/npartphi
!--correct deltar to match new deltaphi
 deltar = ri*deltaphi
!--setup particles in z and phi at this r
 call setring(npartphi,ipart,ri,deltar,deltaphi,densi)
 deltar0 = deltar
 dens0 = densi

!
!--setup rings from Rtorus outwards until dens < 0
!
 ri = Rtorus
 do while (densi > densmin)
    zi = 0.
    !--take half step in r using current deltar
    rtemp = ri + 0.5*deltar
    !--get density
    denstemp = rhofac*rhofunc(rtemp,zi,polyn,dfac,Mstar,Rtorus)
    !--calculate delta r and delta phi from rho
    if (denstemp > tiny(denstemp)) then
       deltartemp = (massp/denstemp)**dndim
       !--corrector step on r using midpoint density
       ri = ri + deltartemp
       !--get density
       densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    else
       densi = 0.
    endif
    if (densi > densmin) then
       zi = 0.
       !--get new deltar at new position
       deltar = (massp/densi)**dndim
       deltaphi = deltar/ri
       npartphi = int(2.*pi/deltaphi)
       !--correct deltaphi to integer number of particles
       deltaphi = 2.*pi/npartphi

       if (ri*deltaphi < 2.0*deltartemp) then
          !--correct deltar to match new deltaphi
          deltar = ri*deltaphi
          call setring(npartphi,ipart,ri,deltar,deltaphi,densi)
          !print*,' ri = ',ri,' npart = ',ipart
       else
          deltar = ri*deltaphi
          print*,'skipping ring, too few particles : ',npartphi
       endif
    endif

 enddo
!
!--setup rings from Rtorus inwards until dens < 0
!
 ri = Rtorus
 deltar = deltar0
 densi = dens0
 do while (densi > densmin)
    zi = 0.
    !--take half step in r using current deltar
    rtemp = ri - 0.5*deltar
    !--get density
    denstemp = rhofac*rhofunc(rtemp,zi,polyn,dfac,Mstar,Rtorus)
    if (denstemp > densmin) then
       !--calculate delta r and delta phi from rho
       deltartemp = (massp/denstemp)**dndim
       !--corrector step on r using midpoint density
       ri = ri - deltartemp
       !--get density
       densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    else
       densi = 0.
    endif
    if (densi > densmin) then
       !--get new deltar at new position
       deltar = (massp/densi)**dndim
       deltaphi = deltar/ri
       npartphi = int(2.*pi/deltaphi)
       !--correct deltaphi to integer number of particles
       deltaphi = 2.*pi/npartphi
       if (ri*deltaphi < 2.0*deltartemp) then
          !--correct deltar to match new deltaphi
          deltar = ri*deltaphi
          call setring(npartphi,ipart,ri,deltar,deltaphi,densi)
          !print*,' ri = ',ri,' npart = ',npart
       else
          deltar = ri*deltaphi
          print*,'skipping ring, too few particles : ',npartphi
       endif
    endif
 enddo

 npart = ipart
 npartoftype(1) = ipart
 print*,'Mtorus = ',Mtorus,massp*npart
 massoftype(1) = Mtorus/real(npart)
 print*,'Number of particles used',npart

 !
 !--balance pressure forces with centripedal acceleration
 !
 !--set magnetic field using plasma beta
 beta = 100. !1.e4
 rhosum = 0.0
 pmassii = massoftype(igas)
 do ii=1,npart
    if (maxphase==maxp) pmassii = massoftype(iamtype(iphase(ii)))
    rhosum = rhosum + rhoh(xyzh(4,ii),pmassii)
 enddo
 rhosum = rhosum/npart

 Bzi = sqrt(2.*polyk*rhosum**gamma/beta)
 print*,' using beta = ',beta
 print*,' using rho = ',rhosum
 if (mhd) then
    print*,'initial Bz field = ',Bzi
 else
    Bzi = 0.
    print*,'No magnetic field for now.'
 endif
 dbeta = 1./beta
 ihavesetupB=.true.

 print*,'setting v to balance pressure gradients'
 !
 !--analytic velocities
 !
 pmassi   = massoftype(igas)
 iamtypei = igas
 do i=1,npart
    if (maxphase==maxp) then
       call get_partinfo(iphase(i),iactivei,iamdusti,iamtypei)
       pmassi = massoftype(iamtypei)
    endif
    densi = rhoh(xyzh(4,i),pmassi)
    Bzi   = sqrt(2.*polyk*densi**gamma/beta)
    rcyl2 = dot_product(xyzh(1:2,i),xyzh(1:2,i))
    rcyl  = sqrt(rcyl2)
    rsph  = sqrt(rcyl2 + xyzh(3,i)*xyzh(3,i))
    ! Velcoity derived from Daniel's torus paper, Equation 124
    v2onr = bigG*Mstar*Rtorus/(rcyl*rcyl2)
    !--compare to frad from SPH forces
    !frad = abs(dot_product(force(1:2,i),x(1:2,i)/rcyl))
    !omegai = sqrt(frad/rpart)

    omegai = sqrt(v2onr/rcyl)
    vxyzu(1,i) = 0.0 !-omegai*xyzh(2,i)
    vxyzu(2,i) = 0.0 !omegai*xyzh(1,i)
    vxyzu(3,i) = 0.

! Poloidal Field
    if (mhd .and. iamtypei==igas) then
       Bxyz(3,i) = dbeta*(densi**2/rcyl &
                      + 2.*Mstar/(polyk*gamma)*densi**(3.-gamma) &
                      *(Rtorus/rcyl**3 - rcyl/rsph**3))
       Bxyz(1,i) = dbeta*((2.*Mstar/(polyk*gamma))*densi**(3.-gamma) &
                       *(xyzh(3,i)/rsph**3))*xyzh(1,i)/rcyl
       Bxyz(2,i) = dbeta*((2.*Mstar/(polyk*gamma))*densi**(3.-gamma) &
                       *(xyzh(3,i)/rsph**3))*xyzh(2,i)/rcyl
       !     print*,'B ',i,' = ',Bxyz(:,i)
    endif
 enddo

! Toroidal field



 return

contains

real function rhofunc(rr,zz,polyn,dd,Mstar,R0)
 real, intent(in) :: rr,zz,polyn,dd,Mstar,R0
 real :: term
!
!--functional form of rho/(A(n+1))^n
!
 term = Mstar/R0*(R0/sqrt(rr**2 + zz**2) - 0.5*(R0/rr)**2 - 1./(2.*dd))
 if (term > tiny(term)) then
    rhofunc = term**polyn
 else
    rhofunc = 0.
 endif

end function rhofunc

!
!--sets several ring of particles above and below the midplane
!  and around in phi
!
subroutine setring(npartphi,ipart,ri,deltar,deltaphi,densi)
 integer, intent(in)    :: npartphi
 integer, intent(inout) :: ipart
 real,    intent(in)    :: ri,deltar,deltaphi,densi
 real :: deltaz,phii,zi,ztemp,denszi,denstemp,deltaztemp
 integer :: i

 deltaz = deltar
 denszi = densi
 zi = 0.
!--upwards from and including the midplane
 do while (denszi > densmin)
    !--now setup ring of particles
    do i = 1,npartphi
       ipart = ipart + 1
       if (ipart > maxp) stop 'error: ipart > maxp'
       phii = (i-1)*deltaphi
       xyzh(1,ipart) = ri*cos(phii)
       xyzh(2,ipart) = ri*sin(phii)
       xyzh(3,ipart) = zi
       xyzh(4,ipart) = hrho(denszi)
       pri = polyk*denszi**gamma
       if (maxvxyzu >= 4) vxyzu(4,ipart) = pri/(denszi*(gamma-1.))
    enddo
    !--take half step in r using current deltaz
    ztemp = zi + 0.5*deltaz
    denstemp = rhofac*rhofunc(ri,ztemp,polyn,dfac,Mstar,Rtorus)
    if (denstemp > densmin) then
       deltaztemp = massp/(denstemp*deltar**2)
       zi = zi + deltaztemp
       !--get density
       denszi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
       deltaz = massp/(denstemp*deltar**2)
    else
       denszi = 0.
    endif
    !print*,'zi = ',ri,zi,denszi,ipart
 enddo
 deltaz = deltar
 denszi = densi
 zi = 0.
!--downwards from midplane
 do while (denszi > densmin)
    !--take half step in r using current deltaz
    ztemp = zi - 0.5*deltaz
    denstemp = rhofac*rhofunc(ri,ztemp,polyn,dfac,Mstar,Rtorus)
    if (denstemp > tiny(denstemp)) then
       deltaztemp = massp/(denstemp*deltar**2)
       zi = zi - deltaztemp
       !--get density
       denszi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
       deltaz = massp/(denstemp*deltar**2)
    else
       denszi = 0.
    endif
    if (denszi > densmin) then
       !--now setup ring of particles
       do i = 1,npartphi
          ipart = ipart + 1
          if (ipart > maxp) stop 'error: ipart > maxp'
          phii = (i-1)*deltaphi
          xyzh(1,ipart) = ri*cos(phii)
          xyzh(2,ipart) = ri*sin(phii)
          xyzh(3,ipart) = zi
          xyzh(4,ipart) = hrho(denszi)
          pri = polyk*denszi**gamma
          if (maxvxyzu >= 4) vxyzu(4,ipart) = pri/(denszi*(gamma-1.))
       enddo
    endif
    !print*,'zi = ',ri,zi,denszi,ipart
 enddo

! print*,'set ring: r = ',ri,' npart = ',ipart

end subroutine setring

end subroutine setpart

end module setup
