!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module rho_profile
!
! This computes several radial density profiles useful for stars
! and gravitational collapse calculations, including:
!
!  1. uniform
!  2. polytrope
!  3. piecewise polytrope
!  4. Evrard
!  5. Bonnor-Ebert sphere
!
! :References: Evrard (1988), MNRAS 235, 911-934
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon, prompting, units
!
 use physcon, only:pi,fourpi
 implicit none

 public  :: rho_uniform,rho_polytrope,rho_piecewise_polytrope, &
            rho_evrard,rho_bonnorebert,prompt_BEparameters
 public  :: calc_mass_enc
 private :: integrate_rho_profile

 abstract interface
  real function func(x)
   real, intent(in) :: x
  end function func
 end interface

contains

!-----------------------------------------------------------------------
!+
!  Uniform density sphere
!+
!-----------------------------------------------------------------------
subroutine rho_uniform(ng,mass,radius,rtab,rhotab)
 integer, intent(in)  :: ng
 real,    intent(in)  :: mass,radius
 real,    intent(out) :: rtab(:),rhotab(:)
 integer              :: i
 real                 :: dr,density

 density = 3.0*mass/(fourpi*radius**3)
 dr      = radius/real(ng)
 do i=1,ng
    rtab(i)   = i*dr
    rhotab(i) = density
 enddo

end subroutine rho_uniform

!-----------------------------------------------------------------------
!+
!  Density profile for a polytrope (assumes G==1)
!+
!-----------------------------------------------------------------------
subroutine rho_polytrope(gamma,polyk,Mstar,rtab,rhotab,npts,rhocentre,set_polyk,Rstar)
 integer, intent(out)             :: npts
 real,    intent(in)              :: gamma
 real,    intent(in)              :: Mstar
 real,    intent(inout)           :: rtab(:),polyk
 real,    intent(out)             :: rhotab(size(rtab))
 real,    intent(inout), optional :: Rstar
 real,    intent(out),   optional :: rhocentre
 logical, intent(in),    optional :: set_polyk
 integer                          :: i,j
 real                             :: r(size(rtab)),v(size(rtab)),den(size(rtab))
 real                             :: dr,an,rhs,Mstar_f,rhocentre0
 real                             :: fac,rfac

 dr   = 0.001
 an   = 1./(gamma-1.)
 v(1) = 0.0
 v(2) = dr*(1.0 - dr*dr/6. )
 r(1) = 0.

 i = 2
 do while (v(i) >= 0.)
    r(i)    = (i-1)*dr
    rhs    = - r(i)*(v(i)/r(i))**an
    v(i+1) = 2*v(i) - v(i-1) + dr*dr*rhs
    i      = i + 1
    if (i+1 > size(rtab)) then ! array is not large enough; restart with larger dr
       dr   = dr*2.
       r(2) = dr
       v(2) = dr*(1.0 - dr*dr/6. )
       i = 2
    endif
 enddo
 npts = i-1
 !
 !--Calculate the mass, Mstar_f, out to radius r using the density without
 !  the central density multiplier.
 !
 den(1) = 1.0
 Mstar_f = 0.
 do j = 2,npts
    den(j)   = (v(j)/r(j))**an
    Mstar_f  = Mstar_f + fourpi*r(j)*r(j)*den(j)*dr
 enddo
 !
 !--Rescale the central density to give desired mass, Mstar
 !  This is using the incorrect polyk
 !
 fac        = (gamma*polyk)/(fourpi*(gamma - 1.))
 rhocentre0 = ((Mstar/Mstar_f)/fac**1.5)**(2./(3.*gamma - 4.))
 rfac       = sqrt(fac*rhocentre0**(gamma - 2.))

 if (present(set_polyk) .and. present(Rstar) ) then
    if ( set_polyk ) then
       !--Rescale radius to get polyk
       rfac      = Rstar/(r(npts)*rfac)
       polyk     = polyk*rfac
       !
       !--Re-rescale central density to give desired mass (using the correct polyk)
       fac        = (gamma*polyk)/(fourpi*(gamma - 1.))
       rhocentre0 = ((Mstar/Mstar_f)/fac**1.5)**(2./(3.*gamma - 4.))
       rfac       = sqrt(fac*rhocentre0**(gamma - 2.))
    endif
 endif

 rtab   = r * rfac
 rhotab = rhocentre0 * den
 if (present(Rstar))     Rstar     = r(npts)*rfac
 if (present(rhocentre)) rhocentre = rhocentre0

end subroutine rho_polytrope

!-----------------------------------------------------------------------
!+
!  Calculate the density profile for a piecewise polytrope
!  Original Authors: Madeline Marshall & Bernard Field
!  Supervisors: James Wurster & Paul Lasky
!+
!-----------------------------------------------------------------------
subroutine rho_piecewise_polytrope(rtab,rhotab,rhocentre,mstar_in,get_dPdrho,npts,ierr)
 integer, intent(out)   :: npts,ierr
 real,    intent(in)    :: mstar_in
 real,    intent(out)   :: rhocentre,rtab(:),rhotab(:)
 integer, parameter     :: itermax = 1000
 integer                :: iter,lastsign
 real                   :: dr,drho,mstar
 logical                :: iterate,bisect
 procedure(func), pointer :: get_dPdrho
 !
 !--initialise variables
 iter      = 0
 ierr      = 0
 drho      = 0.0
 dr        = 30.0/size(rtab)
 rhocentre = 1.0
 lastsign  = 1
 iterate   = .true.
 bisect    = .false.
 !
 !--Iterate to get the correct density profile
 do while ( iterate )
    call integrate_rho_profile(rtab,rhotab,rhocentre,get_dPdrho,dr,npts,ierr)
    if (ierr > 0) then
       !--did not complete the profile; reset dr
       dr   = 2.0*dr
       ierr = 0
    else
       call calc_mass_enc(npts,rtab,rhotab,mstar=mstar)
       !--iterate to get the correct mass
       if (iter==0) then
          rhocentre = rhocentre * (mstar_in/mstar)**(1./3.)
          lastsign  = int( sign(1.0,mstar_in-mstar) )
       elseif (iter==1) then
          drho      = 0.1*rhocentre*lastsign
          lastsign  = int( sign(1.0,mstar_in-mstar) )
       else
          if (bisect) then
             drho = 0.5*drho*lastsign*sign(1.0,mstar_in-mstar)
          else
             if (lastsign /= int( sign(1.0,mstar_in-mstar) ) ) then
                bisect = .true.
                drho   = -0.5*drho
             endif
          endif
       endif
       rhocentre = rhocentre + drho
       lastsign  = int( sign(1.0,mstar_in-mstar) )
       iter      = iter + 1
       !--Converged: exit
       if (abs(mstar_in-mstar) < epsilon(mstar_in)*1.0d4) iterate = .false.
       !--Did not converge: abort
       if (iter > itermax) then
          ierr    = 2
          iterate = .false.
       endif
    endif
 enddo

end subroutine rho_piecewise_polytrope
!-----------------------------------------------------------------------
!  Calculate the density profile using an arbitrary EOS and
!  given a central density
!-----------------------------------------------------------------------
subroutine integrate_rho_profile(rtab,rhotab,rhocentre,get_dPdrho,dr,npts,ierr)
 integer, intent(out) :: npts,ierr
 real,    intent(out) :: rtab(:),rhotab(:)
 real,    intent(in)  :: rhocentre,dr
 integer              :: i
 real                 :: drhodr,dPdrho,dPdrho_prev
 logical              :: iterate
 procedure(func), pointer :: get_dPdrho

 !
 !--Initialise variables
 !
 i         = 1
 ierr      = 0
 rtab      = 0.0
 rhotab    = 0.0
 drhodr    = 0.0
 rtab(1)   = 0.0
 rhotab(1) = rhocentre
 iterate   = .true.
 dPdrho_prev = 0.0
 !
 do while ( iterate )
    i = i + 1
    rhotab(i) = rhotab(i-1) + dr*drhodr
    rtab(i)   = rtab(i-1)   + dr
    dPdrho    = get_dPdrho(rhotab(i))
    if (i==2) then
       drhodr = drhodr - fourpi*rhotab(i-1)**2*dr/dPdrho
    else
       drhodr = drhodr + dr*(drhodr**2/rhotab(i-1) &
              - fourpi*rhotab(i)**2/dPdrho &
              - (dPdrho-dPdrho_prev)/(dr*dPdrho)*drhodr - 2.0*drhodr/rtab(i) )
    endif
    dPdrho_prev = dPdrho
    if (rhotab(i) < 0.0) iterate = .false.
    if (i >=size(rtab)) then
       ierr    = 1
       iterate = .false.
    endif
 enddo

 npts         = i
 rhotab(npts) = 0.0

end subroutine integrate_rho_profile

!-----------------------------------------------------------------------
!  Calculate the enclosed mass of a star
!-----------------------------------------------------------------------
subroutine calc_mass_enc(npts,rtab,rhotab,mtab,mstar)
 integer, intent(in)            :: npts
 real,    intent(in)            :: rtab(:),rhotab(:)
 real,    intent(out), optional :: mtab(:),mstar
 integer                        :: i
 real                           :: ri,ro,menc(npts)

 ro      = 0.5*( rtab(1) + rtab(2) )
 menc(1) = ro**3*rhotab(1)/3.0
 do i = 2,npts-1
    ri      = 0.5*(rtab(i) + rtab(i-1))
    ro      = 0.5*(rtab(i) + rtab(i+1))
    menc(i) = menc(i-1) + rhotab(i)*rtab(i)**2*(ro - ri)
 enddo
 ri         = 0.5*(rtab(npts) + rtab(npts-1))
 menc(npts) = menc(npts-1) + rhotab(npts)*rtab(npts)**2*(rtab(npts) - ri)
 menc       = menc*fourpi

 if (present(mtab))  mtab  = menc
 if (present(mstar)) mstar = menc(npts)

end subroutine calc_mass_enc

!-----------------------------------------------------------------------
!+
!  Calculate the density profile for the Evrard Collapse
!+
!-----------------------------------------------------------------------
subroutine rho_evrard(ng,mass,radius,rtab,rhotab)
 integer, intent(in)  :: ng
 real,    intent(in)  :: mass,radius
 real,    intent(out) :: rtab(:),rhotab(:)
 integer              :: i
 real                 :: dr

 dr = radius/real(ng)
 do i=1,ng
    rtab(i)   = i*dr
    rhotab(i) = mass/(2.0*pi*radius*radius*rtab(i))
 enddo

end subroutine rho_evrard

!-----------------------------------------------------------------------
!+
!  Calculates a Bonnor-Ebert sphere
!
!  Examples:
!  To reproduce the sphere in Wurster & Bate (2019):
!     iBEparam = 5, normalised radius = 7.45; physical mass = 1Msun; fac = 1.0
!  To reproduce the sphere in Saiki & Machida (2020):
!     iBEparam = 4, normalised radius = 12.9; physical radius = 5300au; fac = 6.98
!     cs_sphere = 18900cm/s (this is 10K, assuming gamma = 1)
!     density_contrast = 4.48
!  To define both physical radius & mass, the overdensity factor is automatically changed
!+
!-----------------------------------------------------------------------
subroutine rho_bonnorebert(iBEparam,central_density,edge_density,rBE,xBE,mBE,facBE,csBE,gmw,npts,iBElast,rtab,rhotab,ierr)
 use physcon, only:au,pc,mass_proton_cgs,solarm
 use units,   only:umass,udist
 integer, intent(in)    :: iBEparam,npts
 integer, intent(out)   :: iBElast,ierr
 real,    intent(in)    :: csBE,gmw
 real,    intent(inout) :: rBE,mBE,xBE,facBE,central_density
 real,    intent(out)   :: edge_density,rtab(:),rhotab(:)
 integer                :: j,iu
 real                   :: xi,phi,func,containedmass,dxi,dfunc,rho,dphi
 real                   :: rBE0,fac_close,facBEm,facBEr
 real                   :: mtab(npts)
 logical                :: write_BE_profile = .true.
 logical                :: override_critical = .false.  ! if true, will not error out if the density ratio is too small

 !--Initialise variables
 xi             = 0.0
 phi            = 0.0
 func           = 0.0
 containedmass  = 0.0
 dxi            = 5.01*6.45/float(npts)
 dfunc          = (-exp(phi))*dxi
 rtab           = 0.  ! array of radii
 mtab           = 0.  ! array of enclosed masses
 rhotab         = 0.  ! array of densities
 rhotab(1)      = 1.  ! initial normalised density
 rho            = 1.
 ierr           = 0

 ! initialise variables not required for chosen iBEparam (to avoid errors)
 if (iBEparam/=1 .and. iBEparam/=2 .and. iBEparam/=3) central_density = 3.8d-18
 if (iBEparam/=1 .and. iBEparam/=4 .and. iBEparam/=6) rBE   = 7000.*au/udist
 if (iBEparam/=2 .and. iBEparam/=4 .and. iBEparam/=5) xBE   = 7.45
 if (iBEparam/=3 .and. iBEparam/=5 .and. iBEparam/=6) mBE   = 1.0*solarm/umass
 if (iBEparam/=4 .and. iBEparam/=5)                   facBE = 1.0

 !--Calculate a normalised BE profile out to 5 critical radii
 do j = 2,npts
    xi    = (j-1)*dxi
    func  = func + dfunc
    dphi  = func*dxi
    phi   = phi + dphi
    dfunc = (-exp(phi) - 2.0*func/xi)*dxi
    rho   = exp(phi)
    containedmass = containedmass + fourpi*xi*xi*rho*dxi
    rtab(j)       = xi
    mtab(j)       = containedmass
    rhotab(j)     = rho
 enddo
 iBElast = npts

 !--Determine scaling factors for the BE
 fac_close = 1000.
 if (iBEparam==4) then
    central_density = (csBE*xBE/rBE)**2/fourpi
 elseif (iBEparam==5 .or. iBEparam==6) then
    do j = 1, npts
       if (rtab(j) < xBE) iBElast = j
    enddo
    central_density = (csBE**3*mtab(iBElast)*facBE/mBE)**2/fourpi**3
 endif
 rBE0 = csBE/sqrt(fourpi*central_density)

 !--Scale the entire profile to match the input parameters
 do j = 1, npts
    if (iBEparam == 2 .and. rtab(j) < xBE) iBElast = j

    rtab(j)   = rBE0 * rtab(j)
    mtab(j)   = mtab(j) * central_density*rBE0**3
    rhotab(j) = central_density * rhotab(j)

    if ((iBEparam == 1 .or. iBEparam == 4) .and. rtab(j) < rBE) then
       iBElast = j
    elseif (iBEparam == 3 .and. mtab(j) < mBE) then
       iBElast = j
    endif
 enddo
 !--Set the remaining properties
 if (iBEparam==4) then
    central_density = central_density*facBE
    mtab            = mtab*facBE
    rhotab          = rhotab*facBE
 endif
 if (iBEparam==5) then
    central_density = central_density/sqrt(facBE)
    mtab            = mtab*facBE
    rhotab          = rhotab/sqrt(facBE)
 endif
 if (iBEparam==6) then
    facBEr          = rBE/rtab(iBElast)
    facBEm          = mBE/mtab(iBElast)
    rtab            = rtab*facBEr
    mtab            = mtab*facBEm
    facBE           = facBEm/facBEr**3
    rhotab          = rhotab*facBE
    central_density = central_density*facBE
 endif
 mBE = mtab(iBElast)
 rBE = rtab(iBElast)
 xBE = rBE/rBE0
 edge_density = rhotab(iBElast)

 print*, '------ BE sphere properties --------'
 print*, ' Value of central density (code units) = ',central_density
 print*, ' Value of central density (g/cm^3)     = ',central_density*umass/udist**3
 print*, ' Value of central density (1/cm^3)     = ',central_density*umass/(gmw*mass_proton_cgs*udist**3)
 print*, ' Radius (dimensionless) = ',xBE
 print*, ' Radius (code)          = ',rBE
 print*, ' Radius (cm)            = ',rBE*udist
 print*, ' Radius (au)            = ',rBE*udist/au
 print*, ' Radius (pc)            = ',rBE*udist/pc
 print*, ' Total mass (Msun)      = ',mBE*umass/solarm
 print*, ' Overdensity factor     = ',facBE
 print*, ' rho_c/rho_outer             = ',central_density/edge_density
 print*, ' Equilibrium temperature (K) = ',mBE*umass*pc/(rBE*udist*solarm*2.02)
 print*, '------------------------------------'

 !--Error out if required
 if (central_density/rhotab(iBElast) < 14.1) then
    print*, 'The density ratio between the central and edge densities is too low and the sphere will not collapse.'
    if (.not. override_critical) then
       print*, 'Aborting.'
       ierr = 1
       return
    endif
 endif
 !--Sanity check on enclosed mass
 containedmass = 0.
 do j = 1,iBElast
    if (j == 1) then
       containedmass = containedmass + 4.0*pi/3.0*rhotab(j)*rtab(j)**3
    else
       containedmass = containedmass + 4.0*pi*rhotab(j)*rtab(j)**2*(rtab(j)-rtab(j-1))
    endif
 enddo
 print*, 'By density, the contained mass is ',containedmass
 if (abs(containedmass-mBE)/mBE > 0.05) then
    print*, 'WARNING! The defined mass and input mass are not the same! Aborting.'
    ierr = 1
    return
 endif

 !--Write the scaled BE profile that is to be used
 if (write_BE_profile) then
    open(newunit=iu,file='BonnorEbert.txt')
    write(iu,'(a)') "# [01     r(code)]   [02 M_enc(code)]   [03   rho(code)]"
    do j = 1,iBElast
       write(iu,'(3(1pe18.10,1x))') rtab(j),mtab(j),rhotab(j)
    enddo
    close(iu)
 endif

end subroutine rho_bonnorebert
!-----------------------------------------------------------------------
!  Prompts for the BE sphere
!  (see setup_sphereinbox for read/write_setup commands)
!-----------------------------------------------------------------------
subroutine prompt_BEparameters(iBEparam,rho_cen,rad_phys,rad_norm,mass_phys,fac,umass,udist,au,solarm)
 use prompting, only:prompt
 integer,      intent(out) :: iBEparam
 real,         intent(out) :: rho_cen,rad_phys,rad_norm,mass_phys,fac
 real,         intent(in)  :: au,solarm
 real(kind=8), intent(in)  :: umass,udist

 print*, 'Please select parameters used to fit the BE sphere:'
 print*, 'The pairs are: '
 print*, '  1: central density & physical radius'
 print*, '  2: central density & normalised radius'
 print*, '  3: central density & physical mass'
 print*, '  4: normalised radius & physical radius & overdensity factor'
 print*, '  5: normalised radius & physical mass   & overdensity factor'
 print*, '  6: physical radius & physical mass'
 iBEparam = 5
 call prompt('Please enter your choice now: ',iBEparam,1,6)

 !--Default values
 rho_cen   = 3.8d-18
 rad_phys  = 7000.*au/udist
 rad_norm  = 7.45
 mass_phys = 1.0*solarm/umass
 fac       = 1.0 ! This might need to be removed

 !--Ask for the values depending on iBEparam
 if (iBEparam==1 .or. iBEparam==2 .or. iBEparam==3) call prompt('Enter the central density [cgs]: ',rho_cen,0.)
 if (iBEparam==1 .or. iBEparam==4 .or. iBEparam==6) call prompt('Enter the physical radius [code]: ',rad_phys,0.)
 if (iBEparam==2 .or. iBEparam==4 .or. iBEparam==5) call prompt('Enter the normalised radius (critical==6.45): ',rad_norm,0.)
 if (iBEparam==3 .or. iBEparam==5 .or. iBEparam==6) call prompt('Enter the physical mass [code]: ',mass_phys,0.)
 if (iBEparam==4 .or. iBEparam==5) call prompt('Enter density enhancement factor (for mass = fac*mBE): ',fac,1.)
 rho_cen = rho_cen * udist**3/umass ! convert to code units

end subroutine prompt_BEparameters

end module rho_profile
