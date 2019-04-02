!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: rho_profile
!
!  DESCRIPTION: This contains several density profiles, including
!               1) uniform
!               2) polytrope
!               3) piecewise polytrope
!               4) Evrard
!               5) Read data from MESA file
!               6) Read data from KEPLER file
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: datafiles, eos, physcon, units
!+
!--------------------------------------------------------------------------
module rho_profile
 use physcon, only: pi,fourpi
 implicit none

 public  :: rho_uniform,rho_polytrope,rho_piecewise_polytrope, &
            rho_evrard,read_mesa_file,read_kepler_file
 public  :: calc_mass_enc
 private :: integrate_rho_profile,get_dPdrho

contains

!-----------------------------------------------------------------------
!+
!  Calculate a uniform density profile
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
!  Calculate the density profile for a polytrope (recall G==1)
!+
!-----------------------------------------------------------------------
subroutine rho_polytrope(gamma,polyk,Mstar,rtab,rhotab,npts,rhocentre,set_polyk,Rstar)
 implicit none
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
subroutine rho_piecewise_polytrope(rtab,rhotab,rhocentre,mstar_in,npts,ierr)
 integer, intent(out)   :: npts,ierr
 real,    intent(in)    :: mstar_in
 real,    intent(out)   :: rhocentre,rtab(:),rhotab(:)
 integer, parameter     :: itermax = 1000
 integer                :: iter,lastsign
 real                   :: dr,drho,mstar
 logical                :: iterate,bisect
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
    call integrate_rho_profile(rtab,rhotab,rhocentre,dr,npts,ierr)
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
       else if (iter==1) then
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
!+
!  Calculate the density profile using an arbitrary EOS and
!  given a central density
!+
!-----------------------------------------------------------------------
subroutine integrate_rho_profile(rtab,rhotab,rhocentre,dr,npts,ierr)
 integer, intent(out) :: npts,ierr
 real,    intent(out) :: rtab(:),rhotab(:)
 real,    intent(in)  :: rhocentre,dr
 integer              :: i
 real                 :: drhodr,dPdrho,dPdrho_prev
 logical              :: iterate
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
!+
!  Calculates pressure at a given density
!+
!-----------------------------------------------------------------------
real function get_dPdrho(rho)
 use units, only: unit_density,unit_pressure
 use eos,   only: rhocrit0pwpcgs,rhocrit1pwpcgs,rhocrit2pwpcgs,p1pwpcgs, &
                  gamma0pwp,gamma1pwp,gamma2pwp,gamma3pwp
 real, intent(in)  :: rho
 real              :: rhocrit0pwp,rhocrit1pwp,rhocrit2pwp,presscrit
 real              :: polyk0,polyk1,polyk2,polyk3
 real              :: gamma,polyk

 rhocrit0pwp = rhocrit0pwpcgs/unit_density
 rhocrit1pwp = rhocrit1pwpcgs/unit_density
 rhocrit2pwp = rhocrit2pwpcgs/unit_density
 presscrit   = p1pwpcgs/unit_pressure
 polyk1      = presscrit/rhocrit1pwp**gamma1pwp
 polyk2      = presscrit/rhocrit1pwp**gamma2pwp
 polyk3      = polyk2*rhocrit2pwp**(gamma2pwp-gamma3pwp)
 polyk0      = polyk1*rhocrit0pwp**(gamma1pwp-gamma0pwp)

 if (rho < rhocrit0pwp) then
    gamma = 5./3.
    polyk = polyk0
 elseif (rho < rhocrit1pwp) then
    gamma = gamma1pwp
    polyk = polyk1
 else if (rho < rhocrit2pwp) then
    gamma = gamma2pwp
    polyk = polyk2
 else
    gamma = gamma3pwp
    polyk = polyk3
 endif
 get_dPdrho = gamma * polyk * rho**(gamma-1.0)

end function get_dPdrho

!-----------------------------------------------------------------------
!+
!  Calculate the enclosed mass of a star
!+
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
!  Read in data output by the MESA stellar evolution code
!+
!-----------------------------------------------------------------------
subroutine read_mesa_file(filepath,ng_max,n,rtab,rhotab,ptab,temperature,&
                               enitab,totmass,ierr,mcut,rcut)
 use units,     only:udist,umass,unit_density,unit_pressure,unit_ergg
 use datafiles, only:find_phantom_datafile
 integer,          intent(in)  :: ng_max
 integer,          intent(out) :: ierr,n
 real,             intent(out) :: rtab(:),rhotab(:),ptab(:),temperature(:),enitab(:),totmass
 real,             intent(out), optional :: rcut
 real,             intent(in), optional :: mcut
 character(len=*), intent(in)  :: filepath
 character(len=120)            :: fullfilepath
 integer                       :: i,iread,aloc,iunit
 integer, parameter            :: maxstardatacols = 6
 real                          :: stardata(ng_max,maxstardatacols)
 logical                       :: iexist,n_too_big
 !
 !--Get path name
 !
 ierr = 0
 fullfilepath = find_phantom_datafile(filepath,'star_data_files')
 inquire(file=trim(fullfilepath),exist=iexist)
 if (.not.iexist) then
    ierr = 1
    return
 endif
 !
 !--Read data from file
 !
 n = 0
 stardata(:,:) = 0.
 n_too_big = .false.
 do iread=1,2
    !--open
    open(newunit=iunit, file=trim(fullfilepath), status='old',iostat=ierr)
    if (.not. n_too_big) then
       !--skip two header lines
       read(iunit,*)
       read(iunit,*)
       if (iread==1) then
          !--first reading
          n = 0
          do while (ierr==0 .and. n < size(stardata(:,1)))
             n = n + 1
             read(iunit,*,iostat=ierr)
             if (ierr /= 0) n = n - 1
          enddo
          if (n >= size(stardata(:,1))) n_too_big = .true.
       else
          !--Second reading
          do i=1,n
             read(iunit,*,iostat=ierr) stardata(n-i+1,:)
          enddo
       endif
    endif
    close(iunit)
 enddo
 if (n < 1) then
    ierr = 2
    return
 endif
 if (n_too_big) then
    ierr = 3
    return
 endif
 !
 !--convert relevant data from CGS to code units
 !
 !radius
 stardata(1:n,4)  = stardata(1:n,4)/udist
 rtab(1:n)        = stardata(1:n,4)
 !density
 stardata(1:n,5)  = stardata(1:n,5)/unit_density
 rhotab(1:n)      = stardata(1:n,5)
 !mass
 stardata(1:n,1)  = stardata(1:n,1)/umass
 totmass          = stardata(n,1)
 !pressure
 stardata(1:n,2)  = stardata(1:n,2)/unit_pressure
 ptab(1:n)        = stardata(1:n,2)
 !temp
 temperature(1:n) = stardata(1:n,3)
 !specific internal energy
 stardata(1:n,6)  = stardata(1:n,6)/unit_ergg
 enitab(1:n)      = stardata(1:n,6)

 if (present(rcut) .and. present(mcut)) then
    aloc = minloc(abs(stardata(1:n,1) - mcut),1)
    rcut = rtab(aloc)
    print*, 'rcut = ', rcut
 endif
end subroutine read_mesa_file

!-----------------------------------------------------------------------
!+
!  Read in datafile from the KEPLER stellar evolution code
!+
!-----------------------------------------------------------------------
subroutine read_kepler_file(filepath,ng_max,n,rtab,rhotab,ptab,temperature,&
                               enitab,totmass,ierr,mcut,rcut)
 use units,     only:udist,umass,unit_density,unit_pressure,unit_ergg
 use datafiles, only:find_phantom_datafile
 integer,          intent(in)  :: ng_max
 integer,          intent(out) :: ierr,n
 real,             intent(out) :: rtab(:),rhotab(:),ptab(:),temperature(:),enitab(:)
 real,             intent(out) :: totmass
 real,             intent(out), optional :: rcut
 real,             intent(in), optional :: mcut
 character(len=*), intent(in)  :: filepath
 character(len=120)            :: fullfilepath
 integer                       :: i,iread,aloc,iunit
 integer, parameter            :: maxstardatacols = 9
 real                          :: stardata(ng_max,maxstardatacols)
 logical                       :: iexist,n_too_big
 !
 !--Get path name
 !
 ierr = 0
 fullfilepath = find_phantom_datafile(filepath,'star_data_files')
 inquire(file=trim(fullfilepath),exist=iexist)
 if (.not.iexist) then
    ierr = 1
    return
 endif
 !
 !--Read data from file
 !
 n = 0
 stardata(:,:) = 0.
 n_too_big = .false.
 do iread=1,2
    !--open
    open(newunit=iunit, file=trim(fullfilepath), status='old',iostat=ierr)
    if (.not. n_too_big) then
       !--skip 23 header lines
       do i=1,23
          read(iunit,*)
       enddo
       if (iread==1) then
          !--first reading
          n = 0
          do while (ierr==0 .and. n < size(stardata(:,1)))
             n = n + 1
             read(iunit,*,iostat=ierr)
             if (ierr /= 0) n = n - 1
          enddo
          if (n >= size(stardata(:,1))) n_too_big = .true.
       else
          !--Second reading
          do i=1,n
             read(iunit,*,iostat=ierr) stardata(i,:)
          enddo
          ! fills hole in center of star
          ! copy first row from second row
          !stardata(1,:) = stardata(2,:)
          ! setting mass, radius, velocity to zero
          !stardata(1,3:5) = 0
       endif
    endif
    close(iunit)
 enddo
 if (n < 1) then
    ierr = 2
    return
 endif
 if (n_too_big) then
    ierr = 3
    return
 endif
 !
 !--convert relevant data from CGS to code units
 !
 !radius
 stardata(1:n,4)  = stardata(1:n,4)/udist
 rtab(1:n)        = stardata(1:n,4)

 !density
 stardata(1:n,6)  = stardata(1:n,6)/unit_density
 rhotab(1:n)      = stardata(1:n,6)

 !mass
 stardata(1:n,3)  = stardata(1:n,3)/umass
 totmass          = stardata(n,3)

 !pressure
 stardata(1:n,8)  = stardata(1:n,8)/unit_pressure
 ptab(1:n)        = stardata(1:n,8)

 !temperature
 temperature(1:n) = stardata(1:n,7)

 !specific internal energy
 stardata(1:n,9)  = stardata(1:n,9)/unit_ergg
 enitab(1:n)      = stardata(1:n,9)

 if (present(rcut) .and. present(mcut)) then
    aloc = minloc(abs(stardata(1:n,1) - mcut),1)
    rcut = rtab(aloc)
    print*, 'rcut = ', rcut
 endif

end subroutine read_kepler_file

end module rho_profile
