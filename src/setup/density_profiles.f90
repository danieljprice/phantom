!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module rho_profile
!
! This contains several density profiles, including
!               1) uniform
!               2) polytrope
!               3) piecewise polytrope
!               4) Evrard
!               5) Read data from MESA file
!               6) Read data from KEPLER file
!               7) Bonnor-Ebert sphere
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: datafiles, eos, physcon, prompting, units
!
 use physcon, only: pi,fourpi
 implicit none

 public  :: rho_uniform,rho_polytrope,rho_piecewise_polytrope, &
            rho_evrard,read_mesa_file,read_mesa,read_kepler_file, &
            rho_bonnorebert,prompt_BEparameters
 public  :: write_softened_profile,calc_mass_enc
 private :: integrate_rho_profile,get_dPdrho

contains

!-----------------------------------------------------------------------
!+
!  Option 1:
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
!  Option 2:
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
!  Option 3:
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
!  Calculates pressure at a given density
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
 elseif (rho < rhocrit2pwp) then
    gamma = gamma2pwp
    polyk = polyk2
 else
    gamma = gamma3pwp
    polyk = polyk3
 endif
 get_dPdrho = gamma * polyk * rho**(gamma-1.0)

end function get_dPdrho

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
!  Option 4:
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
!  Option 5:
!  Read in data output by the MESA stellar evolution code
!+
!-----------------------------------------------------------------------
subroutine read_mesa_file(filepath,ng_max,n,rtab,rhotab,ptab,temperature,&
                               enitab,mtab,totmass,ierr,mcut,rcut)
 use units,     only:udist,umass,unit_density,unit_pressure,unit_ergg
 use datafiles, only:find_phantom_datafile
 integer,          intent(in)  :: ng_max
 integer,          intent(out) :: ierr,n
 real,             intent(out) :: rtab(:),rhotab(:),ptab(:),temperature(:),enitab(:),mtab(:),totmass
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
 mtab(1:n)        = stardata(1:n,1)
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
!  Read quantities from MESA profile or from profile in the format of 
!  the P12 star (phantom/data/star_data_files/P12_Phantom_Profile.data) 
!+
!-----------------------------------------------------------------------
subroutine read_mesa(filepath,rho,r,pres,m,ene,temp,Xfrac,Yfrac,Mstar,ierr,cgsunits)
 use physcon,   only:solarm
 use eos,       only:X_in,Z_in
 use fileutils, only:get_nlines,get_ncolumns,string_delete,lcase
 use datafiles, only:find_phantom_datafile
 use units,     only:udist,umass,unit_density,unit_pressure,unit_ergg
 integer                                    :: lines,rows,i,ncols,nheaderlines
 character(len=*), intent(in)               :: filepath
 logical, intent(in), optional              :: cgsunits
 integer, intent(out)                       :: ierr
 character(len=10000)                       :: dumc
 character(len=120)                         :: fullfilepath
 character(len=24),allocatable              :: header(:),dum(:)
 logical                                    :: iexist,usecgs
 real,allocatable,dimension(:,:)            :: dat
 real,allocatable,dimension(:),intent(out)  :: rho,r,pres,m,ene,temp,Xfrac,Yfrac
 real, intent(out)                          :: Mstar

 rows = 0
 usecgs = .false.
 if (present(cgsunits)) usecgs = cgsunits
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
 lines = get_nlines(fullfilepath) ! total number of lines in file

 open(unit=40,file=fullfilepath,status='old')
 call get_ncolumns(40,ncols,nheaderlines)
 if (nheaderlines == 6) then ! Assume file is a MESA profile, and so it has 6 header lines, and (row=3, col=2) = number of zones 
    read(40,'()')
    read(40,'()')
    read(40,*) lines,lines
    read(40,'()')
    read(40,'()')
 else
    lines = lines - nheaderlines
    do i=1,nheaderlines-1
       read(40,'()')
    enddo
 endif

 read(40,'(a)') dumc! counting rows
 call string_delete(dumc,'[')
 call string_delete(dumc,']')
 allocate(dum(500)) ; dum = 'aaa'
 read(dumc,*,end=101) dum
 101 do i = 1,500
    if (dum(i)=='aaa') then
       rows = i-1
       exit
    endif
 enddo

 allocate(header(rows),dat(lines,rows))
 header(1:rows) = dum(1:rows)
 deallocate(dum)

 do i = 1,lines
    read(40,*) dat(lines-i+1,1:rows)
 enddo

 allocate(m(lines),r(lines),pres(lines),rho(lines),ene(lines), &
             temp(lines),Xfrac(lines),Yfrac(lines))

 ! Set mass fractions to default in eos module if not in file
 Xfrac = X_in
 Yfrac = 1. - X_in - Z_in
 do i = 1, rows
    select case(trim(lcase(header(i))))
    case('mass_grams','mass')
       m = dat(1:lines,i)
    case('rho','density')
       rho = dat(1:lines,i)
    case('energy','e_int')
       ene = dat(1:lines,i)
    case('radius','radius_cm')
       r = dat(1:lines,i)
    case('pressure')
       pres = dat(1:lines,i)
    case('temperature')
       temp = dat(1:lines,i)
    case('x_mass_fraction_h')
       Xfrac = dat(1:lines,i)
    case('y_mass_fraction_he')
       Yfrac = dat(1:lines,i)
    end select
 enddo

 if (nheaderlines == 6) m = m * solarm

 if (.not. usecgs) then
    m = m / umass
    r = r / udist
    pres = pres / unit_pressure
    rho = rho / unit_density
    ene = ene / unit_ergg
 endif

 Mstar = m(lines)
end subroutine read_mesa

!----------------------------------------------------------------
!  Write stellar profile in format readable by read_mesa_file;
!  used in star setup to write softened stellar profile.
!----------------------------------------------------------------
subroutine write_softened_profile(outputpath, m, pres, temp, r, rho, ene, Xfrac, Yfrac, csound)
 real, intent(in)                :: m(:),rho(:),pres(:),r(:),ene(:),temp(:)
 real, intent(in), optional      :: Xfrac(:),Yfrac(:),csound(:)
 character(len=120), intent(in)  :: outputpath
 integer                         :: i

 open(1, file = outputpath, status = 'new')

 if (present(Xfrac) .and. present(Yfrac)) then
    if (present(csound)) then
       write(1,'(a)') '[    Mass   ]  [  Pressure ]  [Temperature]  [   Radius  ]  &
       &[  Density  ]  [   E_int   ]  [   Xfrac   ]  [   Yfrac   ]  [Sound speed]'
       write(1,101) (m(i),pres(i),temp(i),r(i),rho(i),ene(i),Xfrac(i),Yfrac(i),csound(i),i=1,size(r))
101    format (es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,&
       2x,es13.6,2x,es13.6)
    else
       write(1,'(a)') '[    Mass   ]  [  Pressure ]  [Temperature]  [   Radius  ]  &
       &[  Density  ]  [   E_int   ]  [   Xfrac   ]  [   Yfrac   ]'
       write(1,102) (m(i),pres(i),temp(i),r(i),rho(i),ene(i),Xfrac(i),Yfrac(i),i=1,size(r))
102    format (es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,&
       2x,es13.6)
    endif
 else
    write(1,'(a)') '[    Mass   ]  [  Pressure ]  [Temperature]  [   Radius  ]  &
    &[  Density  ]  [   E_int   ]'
    write(1,103) (m(i),pres(i),temp(i),r(i),rho(i),ene(i),i=1,size(r))
103 format (es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6)
 endif

 close(1, status = 'keep')

end subroutine write_softened_profile

!-----------------------------------------------------------------------
!+
!  Option 6:
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
!-----------------------------------------------------------------------
!+
!  Option 7:
!  Calculates a Bonnor-Ebert sphere
!  An error will be returned if the user is request a normalised
!  radius > 5x the critical radius or if the density ratio between
!  centre and edge values is not large enough
!
!  Examples:
!  To reproduce the sphere in Wurster & Bate (2019):
!     iBEparam = 5, normalised radius = 7.45; physical mass = 1Msun; fac = 1.0
!  To reproduce the sphere in Saiki & Machida (2020):
!     iBEparam = 4, normalised radius = 12.9; physical radius = 5300au; fac = 6.98
!     cs_sphere = 18900cm/s (this is 10K, assuming gamma = 1)
!     density_contrast = 4.48
!+
!-----------------------------------------------------------------------
subroutine rho_bonnorebert(iBEparam,central_density,edge_density,rBE,xBE,mBE,facBE,csBE,npts,iBElast,rtab,rhotab,ierr)
 use physcon, only:au,pc,mass_proton_cgs,solarm
 use units,   only:umass,udist
 use eos,     only:gmw
 integer, intent(in)    :: iBEparam,npts
 integer, intent(out)   :: iBElast,ierr
 real,    intent(in)    :: csBE
 real,    intent(inout) :: rBE,mBE,xBE,facBE,central_density
 real,    intent(out)   :: edge_density,rtab(:),rhotab(:)
 integer                :: j
 real                   :: xi,phi,func,containedmass,dxi,dfunc,rho,dphi
 real                   :: rBE0,rho1,rho2,fac_close
 real                   :: mtab(npts)
 logical                :: write_BE_profile = .true.
 logical                :: debug = .true.
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
 if (iBEparam==4) central_density = (csBE*xBE/rBE)**2/fourpi
 if (iBEparam==5) then
    do j = 1, npts
       if (rtab(j) < xBE) iBElast = j
    enddo
    central_density = (csBE**3*mtab(iBElast)*facBE/mBE)**2/fourpi**3
 endif
 if (iBEparam==6) then
    do j = 1,npts
       rho1 = (csBE*rtab(j)/rBE)**2/fourpi
       rho2 = (csBE**3*mtab(j)/mBE)**2/fourpi**3
       if (debug) print*, j,rtab(j),rho1,rho2,rho1/rho2
       if (abs(rho1/rho2 - 1.) < fac_close) then
          fac_close = abs(rho1/rho2 - 1.)
          iBElast = j
       endif
    enddo
    central_density = (csBE**3*mtab(iBElast)/mBE)**2/fourpi**3
    !--Error out if required
    if (fac_close > 0.1) then
       print*, 'A BE sphere with the requested mass and radius cannot be constructed.  Aborting.'
       ierr = 1
       return
    endif
 endif
 rBE0 = csBE/sqrt(fourpi*central_density)

 !--Scale the entire profile to match the input parameters
 do j = 1, npts
    if (iBEparam == 2 .and. rtab(j) < xBE) iBElast = j

    rtab(j)   = rBE0 * rtab(j)
    mtab(j)   = mtab(j) * central_density*rBE0**3
    rhotab(j) = central_density * rhotab(j)

    if (iBEparam == 1 .and. rtab(j) < rBE) iBElast = j
    if (iBEparam == 4 .and. rtab(j) < rBE) iBElast = j
    if (iBEparam == 3 .and. mtab(j) < mBE) iBElast = j
 enddo
 !--Set the remaining properties
 if (iBEparam==4) then
    central_density = central_density*facBE
    mtab(iBElast)   = mtab(iBElast)*facBE
    rhotab          = rhotab*facBE
 endif
 if (iBEparam==5) then
    central_density = central_density/sqrt(facBE)
    mtab(iBElast)   = mtab(iBElast)*facBE
    rhotab          = rhotab/sqrt(facBE)
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
 print*, ' rho_c/rho_outer             = ',central_density/edge_density
 print*, ' Equilibrium temperature (K) = ',mBE*umass*pc/(rBE*udist*solarm*2.02)
 print*, '------------------------------------'

 !--Error out if required
 if (iBEparam==6 .and. fac_close > 0.1) then
    print*, 'A BE sphere with the requested mass and radius cannot be constructed.  Aborting.'
    ierr = 1
    return
 endif
 if (central_density/rhotab(iBElast) < 14.1) then
    print*, 'The density ratio between the central and edge densities is too low and the sphere will not collapse.'
    if (.not. override_critical) then
       print*, 'Aborting.'
       ierr = 1
       return
    endif
 endif

 !--Write the scaled BE profile that is to be used
 if (write_BE_profile) then
    open(unit = 26393,file='BonnorEbert.txt')
    write(26393,'(a)') "# [01     r(code)]   [02 M_enc(code)]   [03   rho(code)]"
    do j = 1,iBElast
       write(26393,'(3(1pe18.10,1x))') rtab(j),mtab(j),rhotab(j)
    enddo
 endif

end subroutine rho_bonnorebert
!-----------------------------------------------------------------------
!  Prompts for the BE sphere
!  (see setup_sphereinbox for read/write_setup commands)
!-----------------------------------------------------------------------
subroutine prompt_BEparameters(iBEparam,rho_cen,rad_phys,rad_norm,mass_phys,fac,umass,udist,au,solarm)
 use prompting,    only:prompt
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
 print*, '  6: physical mass & physical radius'
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
!-----------------------------------------------------------------------
end module rho_profile
