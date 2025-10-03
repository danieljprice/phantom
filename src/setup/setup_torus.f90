!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for a rotating torus, for simulations of
! the magnetorotational instability and/or the
! Papaloizoi-Pringle instbaility
!
! Originally written over a few beers in collaboration
! with Owen Matthews and Laure Fouchet. This shows.
!
! Most recently used in Nealon et al. (2018)
!
! :References:
!  - Papaloizoi & Pringle (1984), MNRAS 208, 721
!  - Zurek & Benz (1986), ApJ 308, 123
!  - Nealon et al. (2018), MNRAS 474, 1737
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Mstar   : *central star mass*
!   - Mtorus  : *torus mass*
!   - Rtorus  : *torus radius*
!   - beta    : *plasma beta parameter*
!   - densmax : *maximum density in torus*
!   - densmin : *minimum density cutoff*
!   - dfac    : *density factor*
!   - nlayers : *number of layers for integration*
!   - nrings  : *number of rings for integration*
!   - r_in    : *inner radius for integration*
!   - r_out   : *outer radius for integration*
!   - zmax    : *maximum z for integration*
!
! :Dependencies: dim, infile_utils, io, part, physcon, setup_params, units
!
 use physcon, only:pi
 use dim,     only:maxvxyzu,mhd
 use part,    only:hrho,rhoh,igas
 implicit none
 public :: setpart

 private

 ! Setup parameters
 real :: Rtorus = 1.0
 real :: dfac = 1.01
 real :: Mstar = 1.0
 real :: Mtorus = 4.e-4
 real :: densmax = 0.25e-8
 real :: densmin = 1.e-20
 real :: beta = 100.0
 real :: r_in = 0.01
 real :: r_out = 2.5
 integer :: nrings = 200
 integer :: nlayers = 40
 real :: zmax = 1.5

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,&
           gamma,hfact,time,fileprefix)
 use part,         only:Bxyz
 use physcon,      only:au,solarm,solarr
 use units,        only:set_units,get_G_code
 use infile_utils, only:get_options
 use io,           only:master
 integer,           intent(in)  :: id
 integer,           intent(out) :: npart
 integer,           intent(out) :: npartoftype(:)
 real,              intent(out) :: xyzh(:,:)
 real,              intent(out) :: polyk,gamma,hfact
 real,              intent(out) :: vxyzu(:,:)
 real,              intent(out) :: massoftype(:)
 real,              intent(in)  :: time
 character(len=20), intent(in)  :: fileprefix
 integer :: ipart,npartphi,ierr
 real :: massp,deltar,polyn,np
 real :: ri,zi,rhofac,deltaphi,densi,bigG
 real :: deltar0,dens0,maxp
 real, parameter :: dndim = 1./3.

 ! Read setup parameters from file
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

! call set_units(dist=au,mass=solarm,G=1.d0)
 call set_units(dist=100.*solarr,mass=1.e7*solarm,G=1.d0)

 bigG = get_G_code()
 hfact = 1.2
 gamma = 5./3.

 print*,'Setup for beer-related torus thing (Owen, Laure)'
 print*,'particles: ',npart

!
!--calculate polyk, rhofac and total mass
!
 call calculate_polyk_and_rhofac(polyk,rhofac,Mtorus,bigG,gamma)

 maxp = 100000
 print*,'tentative npart = ',maxp
 np = 0.5*maxp      ! Decrease this factor to allow for higher d values
 massp = Mtorus/real(np)
 massoftype(igas) = massp

!
!--setup first ring at r=Rtorus
!
 ipart = 0
 ri = Rtorus
 zi = 0.
 polyn = 1./(gamma-1.)
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
 call setring(npartphi,ipart,ri,deltar,deltaphi,densi,xyzh,vxyzu,polyk,gamma,massp,rhofac,polyn)
 deltar0 = deltar
 dens0 = densi

!
!--setup rings outward and inward from Rtorus
!
 call setup_rings_outward(ipart,ri,deltar0,dens0,massp,rhofac,polyn,bigG,gamma,xyzh,vxyzu,polyk)
 call setup_rings_inward(ipart,ri,deltar0,dens0,massp,rhofac,polyn,bigG,gamma,xyzh,vxyzu,polyk)

 npart = ipart
 npartoftype(igas) = ipart
 print*,'Mtorus = ',Mtorus,massp*npart
 massoftype(igas) = Mtorus/real(npart)
 print*,'Number of particles used',npart

 !
 !--setup velocities and magnetic fields
 !
 call setup_velocities_and_Bfields(npart,polyk,gamma,bigG,xyzh,vxyzu,Bxyz,massoftype)

end subroutine setpart

!---------------------------------------------
!+
!  density function needed for Torus
!+
!---------------------------------------------
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

!---------------------------------------------
!+
!  setup several rings of particles above
!  and below the midplane and around in phi
!+
!---------------------------------------------
subroutine setring(npartphi,ipart,ri,deltar,deltaphi,densi,xyzh,vxyzu,polyk,gamma,massp,rhofac,polyn)
 integer, intent(in)    :: npartphi
 integer, intent(inout) :: ipart
 real,    intent(in)    :: ri,deltar,deltaphi,densi,polyk,gamma,massp,rhofac,polyn
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: deltaz,phii,zi,ztemp,denszi,denstemp,deltaztemp,pri
 integer :: i

 deltaz = deltar
 denszi = densi
 zi = 0.
!--upwards from and including the midplane
 do while (denszi > densmin)
    !--now setup ring of particles
    do i = 1,npartphi
       ipart = ipart + 1
       if (ipart > size(xyzh,2)) stop 'error: ipart > maxp: use ./phantomsetup --maxp=10000000'
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
          if (ipart > size(xyzh,2)) stop 'error: ipart > maxp: use ./phantomsetup --maxp=10000000'
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

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for torus setup'
 write(iunit,"(/,a)") '# torus parameters'
 call write_inopt(Rtorus,'Rtorus','torus radius',iunit)
 call write_inopt(dfac,'dfac','density factor',iunit)
 call write_inopt(Mstar,'Mstar','central star mass',iunit)
 call write_inopt(Mtorus,'Mtorus','torus mass',iunit)
 call write_inopt(densmax,'densmax','maximum density in torus',iunit)
 call write_inopt(densmin,'densmin','minimum density cutoff',iunit)
 if (mhd) call write_inopt(beta,'beta','plasma beta parameter',iunit)
 write(iunit,"(/,a)") '# integration parameters'
 call write_inopt(r_in,'r_in','inner radius for integration',iunit)
 call write_inopt(r_out,'r_out','outer radius for integration',iunit)
 call write_inopt(nrings,'nrings','number of rings for integration',iunit)
 call write_inopt(nlayers,'nlayers','number of layers for integration',iunit)
 call write_inopt(zmax,'zmax','maximum z for integration',iunit)
 close(iunit)

end subroutine write_setupfile

!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(Rtorus,'Rtorus',db,min=0.,errcount=nerr)
 call read_inopt(dfac,'dfac',db,min=1.,errcount=nerr)
 call read_inopt(Mstar,'Mstar',db,min=0.,errcount=nerr)
 call read_inopt(Mtorus,'Mtorus',db,min=0.,errcount=nerr)
 call read_inopt(densmax,'densmax',db,min=0.,errcount=nerr)
 call read_inopt(densmin,'densmin',db,min=0.,errcount=nerr)
 if (mhd) call read_inopt(beta,'beta',db,min=0.,errcount=nerr)
 call read_inopt(r_in,'r_in',db,min=0.,errcount=nerr)
 call read_inopt(r_out,'r_out',db,min=0.,errcount=nerr)
 call read_inopt(nrings,'nrings',db,min=1,errcount=nerr)
 call read_inopt(nlayers,'nlayers',db,min=1,errcount=nerr)
 call read_inopt(zmax,'zmax',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Calculate polyk and rhofac from torus parameters
!+
!-----------------------------------------------------------------------
subroutine calculate_polyk_and_rhofac(polyk,rhofac,Mtorus_calc,bigG,gamma)
 real, intent(out) :: polyk,rhofac,Mtorus_calc
 real, intent(in)  :: bigG,gamma
 real :: deltar,deltaz,polyn
 integer :: iring,iz
 real :: ri,zi

 polyn = 1./(gamma-1.)
 deltar = (r_out - r_in)/real(nrings-1)
 deltaz = 2.*zmax/real(nlayers-1)

 !--determine polyk from maximum density
 polyk = bigG*Mstar/(Rtorus*(polyn + 1.)*densmax**(gamma-1.))*(dfac - 1.)/(2.*dfac)
 print*,' polyk = ',polyk

 !--rhofac is the factor which multiplies rhofunc to give rho
 rhofac = 1./(polyk*(polyn+1.))**polyn

 !--work out total mass in torus
 Mtorus_calc = 0.
 do iring=1,nrings
    ri = r_in + (iring-1)*deltar
    do iz = 1,nlayers
       zi = (iz - 1)*deltaz - zmax
       Mtorus_calc = Mtorus_calc + 2.*pi*ri*deltar*deltaz*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    enddo
 enddo
 Mtorus_calc = Mtorus_calc*rhofac
 print*,'torus mass (by integration) = ',Mtorus_calc
 if (Mtorus_calc > Mstar) stop 'error Mtorus > Mstar'

end subroutine calculate_polyk_and_rhofac

!-----------------------------------------------------------------------
!+
!  Setup rings outward from Rtorus
!+
!-----------------------------------------------------------------------
subroutine setup_rings_outward(ipart,ri,deltar0,dens0,massp,rhofac,polyn,bigG,gamma,xyzh,vxyzu,polyk)
 integer, intent(inout) :: ipart
 real,    intent(inout) :: ri
 real,    intent(in)    :: deltar0,dens0,massp,rhofac,polyn,bigG,gamma,polyk
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: deltar,deltaphi,densi,denstemp,rtemp,deltartemp
 integer :: npartphi
 real, parameter :: dndim = 1./3.

 ri = Rtorus
 deltar = deltar0
 densi = dens0
 deltartemp = deltar0  ! Initialize to avoid uninitialized variable warning
 do while (densi > densmin)
    !--take half step in r using current deltar
    rtemp = ri + 0.5*deltar
    !--get density
    denstemp = rhofac*rhofunc(rtemp,0.,polyn,dfac,Mstar,Rtorus)
    !--calculate delta r and delta phi from rho
    if (denstemp > tiny(denstemp)) then
       deltartemp = (massp/denstemp)**dndim
       !--corrector step on r using midpoint density
       ri = ri + deltartemp
       !--get density
       densi = rhofac*rhofunc(ri,0.,polyn,dfac,Mstar,Rtorus)
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
          call setring(npartphi,ipart,ri,deltar,deltaphi,densi,xyzh,vxyzu,polyk,gamma,massp,rhofac,polyn)
       else
          deltar = ri*deltaphi
          print*,'skipping ring, too few particles : ',npartphi
       endif
    endif
 enddo

end subroutine setup_rings_outward

!-----------------------------------------------------------------------
!+
!  Setup rings inward from Rtorus
!+
!-----------------------------------------------------------------------
subroutine setup_rings_inward(ipart,ri,deltar0,dens0,massp,rhofac,polyn,bigG,gamma,xyzh,vxyzu,polyk)
 integer, intent(inout) :: ipart
 real,    intent(inout) :: ri
 real,    intent(in)    :: deltar0,dens0,massp,rhofac,polyn,bigG,gamma,polyk
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: deltar,deltaphi,densi,denstemp,rtemp,deltartemp
 integer :: npartphi
 real, parameter :: dndim = 1./3.

 ri = Rtorus
 deltar = deltar0
 densi = dens0
 deltartemp = deltar0  ! Initialize to avoid uninitialized variable warning
 do while (densi > densmin)
    !--take half step in r using current deltar
    rtemp = ri - 0.5*deltar
    !--get density
    denstemp = rhofac*rhofunc(rtemp,0.,polyn,dfac,Mstar,Rtorus)
    if (denstemp > densmin) then
       !--calculate delta r and delta phi from rho
       deltartemp = (massp/denstemp)**dndim
       !--corrector step on r using midpoint density
       ri = ri - deltartemp
       !--get density
       densi = rhofac*rhofunc(ri,0.,polyn,dfac,Mstar,Rtorus)
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
          call setring(npartphi,ipart,ri,deltar,deltaphi,densi,xyzh,vxyzu,polyk,gamma,massp,rhofac,polyn)
       else
          deltar = ri*deltaphi
          print*,'skipping ring, too few particles : ',npartphi
       endif
    endif
 enddo

end subroutine setup_rings_inward

!-----------------------------------------------------------------------
!+
!  Setup velocities and magnetic fields
!+
!-----------------------------------------------------------------------
subroutine setup_velocities_and_Bfields(npart,polyk,gamma,bigG,xyzh,vxyzu,Bxyz,massoftype)
 use setup_params, only:ihavesetupB
 integer, intent(in) :: npart
 real,    intent(in) :: polyk,gamma,bigG
 real,    intent(in) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:),Bxyz(:,:),massoftype(:)
 integer :: i
 real :: densi,rcyl2,rcyl,rsph,omegai,v2onr,Bzi,dbeta,rhosum,pmassi

 !--set magnetic field using plasma beta
 rhosum = 0.0
 pmassi = massoftype(igas)
 do i=1,npart
    rhosum = rhosum + rhoh(xyzh(4,i),pmassi)
 enddo
 rhosum = rhosum/npart

 print*,' using rho = ',rhosum
 if (mhd) then
    print*,'initial Bz field = ',Bzi
    Bzi = sqrt(2.*polyk*rhosum**gamma/beta)
    print*,' using beta = ',beta
 else
    Bzi = 0.
    print*,'No magnetic field for now.'
 endif
 ihavesetupB=.true.

 print*,'setting v to balance pressure gradients'
 !
 !--analytic velocities
 !
 do i=1,npart
    densi = rhoh(xyzh(4,i),pmassi)
    rcyl2 = dot_product(xyzh(1:2,i),xyzh(1:2,i))
    rcyl  = sqrt(rcyl2)
    rsph  = sqrt(rcyl2 + xyzh(3,i)*xyzh(3,i))
    ! Velcoity derived from Daniel's torus paper, Equation 124
    v2onr = bigG*Mstar*Rtorus/(rcyl*rcyl2)
    omegai = sqrt(v2onr/rcyl)
    vxyzu(1,i) = -omegai*xyzh(2,i)
    vxyzu(2,i) = omegai*xyzh(1,i)
    vxyzu(3,i) = 0.

! Poloidal Field
    if (mhd) then
       dbeta = 1./beta
       Bxyz(3,i) = dbeta*(densi**2/rcyl &
                      + 2.*Mstar/(polyk*gamma)*densi**(3.-gamma) &
                      *(Rtorus/rcyl**3 - rcyl/rsph**3))
       Bxyz(1,i) = dbeta*((2.*Mstar/(polyk*gamma))*densi**(3.-gamma) &
                       *(xyzh(3,i)/rsph**3))*xyzh(1,i)/rcyl
       Bxyz(2,i) = dbeta*((2.*Mstar/(polyk*gamma))*densi**(3.-gamma) &
                       *(xyzh(3,i)/rsph**3))*xyzh(2,i)/rcyl
    endif
 enddo

end subroutine setup_velocities_and_Bfields

end module setup
