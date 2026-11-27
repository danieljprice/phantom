!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setstar_utils
!
! Utility routines for mapping 1D stars into 3D spheres
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, eos_piecewise, extern_densprofile, io, kernel,
!   mpiutils, part, physcon, radiation_utils, readwrite_aton,
!   readwrite_kepler, readwrite_mesa, rho_profile, setsoftenedcore,
!   sortutils, spherical, table_utils, unifdis, units
!
 use extern_densprofile, only:nrhotab
 use readwrite_kepler,   only:write_kepler_comp
 use readwrite_aton,     only:read_aton,write_aton
 implicit none
 !
 ! Index of setup options
 !
 integer, parameter, public :: nprofile_opts =  8 ! maximum number of initial configurations
 integer, parameter, public :: ipointmass = 0
 integer, parameter, public :: iuniform   = 1
 integer, parameter, public :: ipoly      = 2
 integer, parameter, public :: ifromfile  = 3
 integer, parameter, public :: ikepler    = 4
 integer, parameter, public :: imesa      = 5
 integer, parameter, public :: ibpwpoly   = 6
 integer, parameter, public :: ievrard    = 7
 integer, parameter, public :: iaton      = 8

 character(len=*), parameter, public :: profile_opt(0:nprofile_opts) = &
    (/'Sink particle/point mass    ', &
      'Uniform density sphere      ', &
      'Polytrope                   ', &
      'Density vs r from ascii file', &
      'KEPLER star from file       ', &
      'MESA star from file         ', &
      'Piecewise polytrope         ', &
      'Evrard collapse             ', &
      'ATON star from file         '/)

 public :: read_star_profile
 public :: set_star_density
 public :: set_star_composition
 public :: set_star_thermalenergy
 public :: set_stellar_core
 public :: write_kepler_comp
 public :: need_inputprofile,need_polyk,need_rstar,need_mu
 public :: get_mass_coord
 public :: check_convective_stability
 public :: get_softened_profile_filename
 public :: get_core_softening_hints,print_core_softening_hints

 private

 integer, parameter :: ng_max = nrhotab

contains

!-------------------------------------------------------------------------------
!+
!  read stellar profile (density, pressure, temperature, composition) from file
!+
!-------------------------------------------------------------------------------
subroutine read_star_profile(iprofile,ieos,input_profile,gamma,polyk,ui_coef,&
                             r,den,pres,temp,en,mtab,X_in,Z_in,Xfrac,Yfrac,mu,&
                             npts,rmin,Rstar,Mstar,rhocentre,&
                             isoftcore,isofteningopt,rcore,mcore,hsoft,outputfilename,&
                             composition,comp_label,columns_compo,np,hfact)
 use extern_densprofile, only:read_rhotab_wrapper
 use eos_piecewise,      only:get_dPdrho_piecewise
 use kernel,             only:radkern
 use part,               only:do_radiation
 use rho_profile,        only:rho_uniform,rho_polytrope,rho_piecewise_polytrope,rho_evrard,func,calc_mass_enc
 use eos,                only:gmw,get_mean_molecular_weight
 use readwrite_mesa,     only:read_mesa,write_mesa
 use readwrite_kepler,   only:read_kepler_file
 use setsoftenedcore,    only:set_softened_core
 use io,                 only:fatal
 use units,              only:udist,umass
 use physcon,            only:solarr,solarm
 integer,                        intent(in)    :: iprofile,ieos
 character(len=*),               intent(in)    :: input_profile,outputfilename
 real,                           intent(in)    :: ui_coef
 real,                           intent(inout) :: gamma,polyk,hsoft
 real,                           intent(in)    :: X_in,Z_in
 real, allocatable,              intent(out)   :: r(:),den(:),pres(:),temp(:),en(:),mtab(:)
 real, allocatable,              intent(out)   :: Xfrac(:),Yfrac(:),mu(:),composition(:,:)
 integer,                        intent(out)   :: npts
 real,                           intent(inout) :: rmin,Rstar,Mstar,rhocentre
 integer,                        intent(in)    :: isoftcore,isofteningopt
 real,                           intent(inout) :: rcore,mcore
 integer,                        intent(out)   :: columns_compo
 character(len=20), allocatable, intent(out)   :: comp_label(:)
 integer,                        intent(in)    :: np
 real,                           intent(in)    :: hfact
 integer :: ierr,eos_type,i
 logical :: calc_polyk,iexist,regrid_core
 procedure(func), pointer :: get_dPdrho
 !
 ! set up tabulated density profile
 !
 character(len=120) :: profile_filename
 calc_polyk = .true.
 allocate(r(ng_max),den(ng_max),pres(ng_max),temp(ng_max),en(ng_max),mtab(ng_max))
 temp = 0.  ! Initialize temperature array for non-file-based profiles
 mtab = 0.  ! mass coordinate is computed from the density profile if not read from file
 print "(/,a,/)",' Using '//trim(profile_opt(iprofile))
 select case(iprofile)
 case(ipoly)
    polyk = 1.  ! overwrite initial value of polyk
    call rho_polytrope(gamma,polyk,Mstar,r,den,npts,rhocentre,calc_polyk,Rstar)
    rmin = r(1)
    pres(1:npts) = polyk*den(1:npts)**gamma
 case(ifromfile)
    call read_rhotab_wrapper(trim(input_profile),ng_max,r,den,npts,&
                             polyk,gamma,rhocentre,Mstar,iexist,ierr)
    if (.not.iexist) call fatal('setup','density file does not exist')
    if (ierr > 0)    call fatal('setup','error in reading density file')
    rmin = r(1)
    pres(1:npts) = polyk*den(1:npts)**gamma
 case(ibpwpoly)
    get_dPdrho => get_dPdrho_piecewise
    call rho_piecewise_polytrope(r,den,rhocentre,Mstar,get_dPdrho,npts,ierr)
    if (ierr == 1) call fatal('setup','ng_max is too small')
    if (ierr == 2) call fatal('setup','failed to converge to a self-consistent density profile')
    rmin  = r(1)
    Rstar = r(npts)
    pres(1:npts) = polyk*den(1:npts)**gamma
 case(imesa,ikepler,iaton)
    deallocate(r,den,pres,temp,en,mtab)
    if (isoftcore > 0) then
       if (iprofile == imesa) then
          call read_mesa(input_profile,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr,cgsunits=.true.)
       elseif (iprofile == iaton) then
          call read_aton(input_profile,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr,cgsunits=.true.)
       else
          call read_kepler_file(trim(input_profile),ng_max,npts,r,den,pres,mtab,temp,en,&
                           Mstar,composition,comp_label,Xfrac,Yfrac,columns_compo,ierr,cgsunits=.true.)
       endif
       if (ierr /= 0) call fatal('setup','error in reading stellar profile from '//trim(input_profile))
       ! print advice on choosing the core radius and mass, and check the values supplied
       call print_core_softening_hints(r,mtab,den,np,hfact,rcore*udist,mcore*umass)
       if (do_radiation) then
          eos_type = 12
       else
          eos_type = ieos
       endif
       regrid_core = (isoftcore==2)  ! uniform grid for the fixed-entropy softened region
       mcore = mcore * umass / solarm ! mcore and rcore needed in solar units for set softened star
       rcore = rcore * udist / solarr
       call set_softened_core(eos_type,isoftcore,isofteningopt,regrid_core,rcore,mcore,r,den,pres,mtab,Xfrac,Yfrac,ierr)
       mcore = mcore * solarm / umass
       rcore = rcore * solarr / udist
       hsoft = rcore/radkern

       call solve_uT_profiles(eos_type,r,den,pres,Xfrac,Yfrac,regrid_core,temp,en,mu)
       if (iprofile == iaton) then
          call write_aton(outputfilename,mtab,pres,temp,r,den,en,Xfrac,Yfrac,mu=mu)
          ! now read the softened profile instead
          call read_aton(outputfilename,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)
       else
          call write_mesa(outputfilename,mtab,pres,temp,r,den,en,Xfrac,Yfrac,mu=mu)
          ! now read the softened profile instead
          call read_mesa(outputfilename,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)
       endif
       profile_filename = outputfilename
    else
       if (iprofile == imesa) then
          call read_mesa(input_profile,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)
       elseif (iprofile == iaton) then
          call read_aton(input_profile,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)
       else
          call read_kepler_file(trim(input_profile),ng_max,npts,r,den,pres,mtab,temp,en,&
                Mstar,composition,comp_label,Xfrac,Yfrac,columns_compo,ierr)
       endif
       profile_filename = input_profile
    endif
    if (ierr==1) call fatal('set_star',trim(profile_filename)//' does not exist')
    if (ierr==2) call fatal('set_star','insufficient data points read from file')
    if (ierr==3) call fatal('set_star','too many data points; increase ng')
    if (ierr /= 0) call fatal('set_star','error in reading stellar profile from '//trim(profile_filename))
    npts = size(den)
    rmin  = r(1)
    Rstar = r(npts)
 case(ievrard)
    call rho_evrard(ng_max,Mstar,Rstar,r,den)
    npts = ng_max
    polyk = ui_coef*Mstar/Rstar
    rmin = r(1)
    pres = polyk*den**gamma
    print*,' Assuming polyk = ',polyk
 case default  ! set up uniform sphere by default
    call rho_uniform(ng_max,Mstar,Rstar,r,den) ! use this array for continuity of call to set_sphere
    npts = ng_max
    rmin = r(1)
    pres = polyk*den**gamma
    print*,' Assuming polyk = ',polyk
 end select
 !
 ! always return the mass coordinate and mean molecular weight tables,
 ! computed from the density profile or the composition if not read from file
 !
 if (all(mtab(1:npts) < tiny(0.))) call calc_mass_enc(npts,r,den,mtab=mtab(1:npts))
 if (.not.allocated(mu)) then
    allocate(mu(size(r)))
    mu = gmw
 endif
 if (.not.allocated(Xfrac)) then
    allocate(Xfrac(size(r)),Yfrac(size(r)))
    Xfrac = X_in
    Yfrac = 1. - X_in - Z_in
 endif
 do i=1,size(mu)
    if (mu(i) <= 0.) mu(i) = get_mean_molecular_weight(Xfrac(i),1.-Xfrac(i)-Yfrac(i))
 enddo

end subroutine read_star_profile

!-------------------------------------------------------------------------------
!+
!  query function for whether a particular profile needs to read an input file
!+
!-------------------------------------------------------------------------------
logical function need_inputprofile(iprofile)
 integer, intent(in) :: iprofile

 select case(iprofile)
 case(imesa,ikepler,iaton,ifromfile)
    need_inputprofile = .true.
 case default
    need_inputprofile = .false.
 end select

end function need_inputprofile

!-------------------------------------------------------------------------------
!+
!  query function for whether a particular profile needs to read the
!  polytropic constant
!+
!-------------------------------------------------------------------------------
logical elemental function need_polyk(iprofile)
 integer, intent(in) :: iprofile

 select case(iprofile)
 case(iuniform,ibpwpoly)
    need_polyk = .true.
 case default
    need_polyk = .false.
 end select

end function need_polyk

!-------------------------------------------------------------------------------
!+
!  query function for whether mean molecular weight is needed
!+
!-------------------------------------------------------------------------------
logical elemental function need_mu(isoftcore)
 integer, intent(in) :: isoftcore

 need_mu = (isoftcore <= 0)

end function need_mu

!-------------------------------------------------------------------------------
!+
!  query function for whether a particular profile needs to read the
!  stellar radius
!+
!-------------------------------------------------------------------------------
logical function need_rstar(iprofile)
 integer, intent(in) :: iprofile

 select case(iprofile)
 case(ibpwpoly,:0) ! no Rstar for iprofile <= 0
    need_rstar = .false.
 case default
    need_rstar = .true.
 end select

end function need_rstar

!-------------------------------------------------------------------------------
!+
!  set up particles according to the desired density profile
!+
!-------------------------------------------------------------------------------
subroutine set_star_density(lattice,id,master,rmin,Rstar,Mstar,hfact,&
                            npts,den,r,npart,npartoftype,massoftype,xyzh,&
                            use_exactN,np,rhozero,npart_total,mask)
 use part,         only:set_particle_type,igas
 use spherical,    only:set_sphere
 use unifdis,      only:mask_prototype
 use physcon,      only:pi
 use mpiutils,     only:reduceall_mpi
 character(len=*), intent(in)    :: lattice
 integer,          intent(in)    :: id,master,npts,np
 integer(kind=8),  intent(out)   :: npart_total
 real,             intent(in)    :: rmin,Rstar,Mstar,hfact
 real,             intent(in)    :: den(npts),r(npts)
 real,             intent(out)   :: rhozero
 integer,          intent(inout) :: npart,npartoftype(:)
 real,             intent(inout) :: xyzh(:,:),massoftype(:)
 logical,          intent(in)    :: use_exactN
 procedure(mask_prototype) :: mask
 integer :: nx,i,ntot,npart_old,n
 integer(kind=8) :: npart_old_total
 real :: vol_sphere,psep
 logical :: mass_is_set

 !
 ! if this is the first call to set_star_density (npart_old=0),
 ! set the number of particles to the requested np and set
 ! the particle mass accordingly, otherwise assume particle mass
 ! is already set and use Mstar to determine np
 !
 npart_old = npart
 n = np
 mass_is_set = .false.
 if (massoftype(igas) > tiny(0.)) then
    n = nint(Mstar/massoftype(igas))
    mass_is_set = .true.
    print "(a,i0)",' WARNING: particle mass is already set, using np = ',n
 endif
 !
 ! place particles in sphere
 !
 vol_sphere  = 4./3.*pi*Rstar**3
 nx          = int(n**(1./3.))
 psep        = vol_sphere**(1./3.)/real(nx)
 ! set_sphere adds the new particles to the total number across all MPI threads
 npart_old_total = reduceall_mpi('+',int(npart_old,kind=8))
 npart_total     = npart_old_total
 call set_sphere(lattice,id,master,rmin,Rstar,psep,hfact,npart,xyzh, &
                 rhotab=den(1:npts),rtab=r(1:npts),nptot=npart_total, &
                 exactN=use_exactN,np_requested=n,mask=mask)
 !
 ! get total particle number across all MPI threads
 ! and use this to set the particle mass
 !
 ntot = int(npart_total - npart_old_total,kind=(kind(ntot)))  ! number of particles in this star
 if (.not.mass_is_set) then
    massoftype(igas) = Mstar/ntot
 endif
 !
 ! set particle type as gas particles
 !
 npartoftype(igas) = npartoftype(igas) + npart - npart_old   ! npart is number on this thread only
 do i=npart_old+1,npart
    call set_particle_type(i,igas)
 enddo
 !
 ! mean density
 !
 rhozero = Mstar/vol_sphere

end subroutine set_star_density

!-----------------------------------------------------------------------
!+
!  Add a sink particle as a stellar core
!+
!-----------------------------------------------------------------------
subroutine set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,ihsoft,mcore,&
                            hsoft,ilum,lcore,iptmass_core,ierr)
 integer, intent(out) :: nptmass,ierr,iptmass_core
 real,    intent(out) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(in)  :: mcore,hsoft,lcore
 integer, intent(in)  :: ihsoft,ilum
 integer              :: n

 ierr = 0
 ! Check for sensible values
 if (mcore < tiny(mcore)) then
    ierr = 1
    return
 endif
 if (hsoft < tiny(hsoft)) then
    ierr = 2
    return
 endif
 if (lcore < 0.) then
    ierr = 3
    return
 endif

 nptmass                = 1
 n                      = nptmass
 xyzmh_ptmass(:,n)      = 0. ! zero all quantities by default
 xyzmh_ptmass(4,n)      = mcore
 xyzmh_ptmass(ihsoft,n) = hsoft
 xyzmh_ptmass(ilum,n)   = lcore
 vxyz_ptmass(:,n)       = 0.
 iptmass_core = n

end subroutine set_stellar_core

!----------------------------------------------------------------
!+
!  get the mass coordinate of a particle m(r)
!  this gives the mass enclosed EXCLUSIVE of self, i.e. m(<r)
!+
!----------------------------------------------------------------
subroutine get_mass_coord(i1,npart,xyzh,mass_enclosed_r,x0)
 use dim,       only:use_apr
 use part,      only:igas,apr_level,massoftype,aprmassoftype,sync_aprmassoftype
 use sortutils, only:sort_by_radius
 integer,           intent(in)  :: i1,npart
 real,              intent(in)  :: xyzh(:,:),x0(3)
 real, allocatable, intent(out) :: mass_enclosed_r(:)
 integer, allocatable :: iorder(:)
 real :: massri,mass_at_r,pmassi,r2,r2prev
 integer :: i,j,iprev

 ! allocate memory
 allocate(mass_enclosed_r(npart-i1),iorder(npart-i1))

 ! massoftype may have just been set by the calling routine
 call sync_aprmassoftype()

 ! sort particles by radius
 call sort_by_radius(npart-i1,xyzh(:,i1+1:npart),iorder,x0)

 ! calculate cumulative mass
 massri = 0.
 mass_at_r = 0.
 r2prev = huge(0.)
 iprev = i1+1
 do i=i1+1,npart
    j = i1 + iorder(i-i1)
    if (use_apr) then
       pmassi = aprmassoftype(igas,apr_level(j))
    else
       pmassi = massoftype(igas)
    endif
    r2 = dot_product(xyzh(1:3,j)-x0,xyzh(1:3,j)-x0)
    !
    ! key point here is to handle the situation where particles are at the same
    ! radius, in which case they should get the same mass coordinate so that
    ! they interpolate identical stellar properties from the table
    !
    if (abs(r2 - r2prev) > tiny(0.)) then
       massri = massri + mass_at_r
       mass_at_r = pmassi
       r2prev = r2
    else
       mass_at_r = mass_at_r + pmassi
    endif
    mass_enclosed_r(j-i1) = massri
 enddo

 ! clean up
 deallocate(iorder)

end subroutine get_mass_coord

!-----------------------------------------------------------------------
!+
!  Set the composition, if variable composition is used
!+
!-----------------------------------------------------------------------
subroutine set_star_composition(eos_outputs_mu,npart,xyzh,Xfrac,Yfrac,&
           mu,mtab,eos_vars,npin,x0)
 use part,        only:iX,iZ,imu
 use table_utils, only:yinterp
 logical, intent(in)  :: eos_outputs_mu
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(in)  :: Xfrac(:),Yfrac(:),mu(:),mtab(:)
 real,    intent(out) :: eos_vars(:,:)
 real,    intent(in), optional :: x0(3)
 integer, intent(in), optional :: npin
 real, allocatable :: mass_enclosed_r(:)
 real :: massri,xorigin(3)
 integer :: i,i1

 i1 = 0
 if (present(npin)) i1 = npin  ! starting position in particle array

 xorigin = 0.
 if (present(x0)) xorigin = x0

 ! this does NOT work with MPI
 call get_mass_coord(i1,npart,xyzh,mass_enclosed_r,xorigin)

 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(i1,npart,mass_enclosed_r,eos_outputs_mu) &
 !$omp shared(Xfrac,Yfrac,mtab,mu,eos_vars) &
 !$omp private(i,massri)
 do i = i1+1,npart
    massri = mass_enclosed_r(i-i1)
    eos_vars(iX,i) = yinterp(Xfrac,mtab,massri)
    eos_vars(iZ,i) = 1. - eos_vars(iX,i) - yinterp(Yfrac,mtab,massri)
    if (.not. eos_outputs_mu) eos_vars(imu,i) = yinterp(mu,mtab,massri)
 enddo
 !$omp end parallel do

end subroutine set_star_composition

!-----------------------------------------------------------------------
!+
!  Set the thermal energy profile
!+
!-----------------------------------------------------------------------
subroutine set_star_thermalenergy(ieos,den,pres,temp,r,npts,npart,xyzh,vxyzu,rad,eos_vars,&
                                  relaxed,use_var_comp,initialtemp,polyk_in,npin,x0)
 use part,            only:rho,itemp,igasP,iX,iZ,imu,iradxi,icv,radprop
 use eos,             only:equationofstate,calc_temp_and_ene,eos_outputs_mu,get_cv,gmw
 use radiation_utils, only:radxi_from_Trad
 use table_utils,     only:yinterp
 use units,           only:unit_density,unit_ergg,unit_pressure
 use dim,             only:do_radiation
 use physcon,         only:Rg,radconst
 use io,              only:fatal
 integer, intent(in)    :: ieos,npart,npts
 real,    intent(in)    :: den(:),pres(:),temp(:),r(:)  ! density, pressure, and temperature tables
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:),eos_vars(:,:),rad(:,:)
 logical, intent(in)    :: relaxed,use_var_comp
 real,    intent(in)    :: initialtemp,polyk_in
 integer, intent(in), optional :: npin
 real,    intent(in), optional :: x0(3)
 integer :: eos_type,cv_type,i,ierr
 real    :: presi,densi,tempi,eni,ri,egasrad,eint,mu
 real    :: rho_cgs,p_cgs,u_gasrec,xorigin(3),dum
 logical :: do_radiation_local
 integer :: i1

 i1  = 0
 eni = 0. ! to prevent compiler warning
 tempi = 0.
 if (present(npin)) i1 = npin  ! starting position in particle array

 xorigin = 0.
 if (present(x0)) xorigin = x0

 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(i1,npart,xyzh,vxyzu,rad,eos_vars,den,pres,temp,r,npts) &
 !$omp shared(relaxed,use_var_comp,ieos,initialtemp,polyk_in,rho) &
 !$omp shared(xorigin,unit_density,unit_ergg,unit_pressure) &
 !$omp shared(radprop,gmw) &
 !$omp private(i,densi,presi,ri,tempi,eni,rho_cgs,p_cgs) &
 !$omp private(egasrad,eint,mu,u_gasrec) &
 !$omp private(dum,eos_type,cv_type,ierr,do_radiation_local)
 do i = i1+1,npart
    ri = sqrt(dot_product(xyzh(1:3,i)-xorigin,xyzh(1:3,i)-xorigin))
    if (relaxed .and. rho(i) > tiny(rho)) then
       densi = rho(i)
       presi = eos_vars(igasP,i)  ! retrieve pressure from relax_star calculated with the fake (ieos=2) internal energy
    else
       ! interpolate density and pressure from table
       densi = yinterp(den(1:npts),r(1:npts),ri)
       presi = yinterp(pres(1:npts),r(1:npts),ri)
    endif

    rho_cgs = densi*unit_density
    p_cgs = presi*unit_pressure
    if (ieos==15 .and. temp(1) > 0.) then       ! should really be a check if we actually have the temperature table
       tempi = yinterp(temp(1:npts),r(1:npts),ri)  ! use MESA temperature as initial guess for Helmholtz
    elseif (rho_cgs > tiny(rho_cgs)) then
       tempi = min((3.*p_cgs/radconst)**0.25, p_cgs/(rho_cgs*Rg))  ! temperature guess
    else
       tempi = (3.*p_cgs/radconst)**0.25
    endif

    if (do_radiation) then
       select case(ieos)
       case(2)
          eos_type = 12  ! Calculate temperature from both gas and radiation pressure
          cv_type = 0
          if (use_var_comp) then
             mu = eos_vars(imu,i)
          else
             mu = gmw
          endif
          call calc_temp_and_ene(eos_type,rho_cgs,p_cgs,egasrad,tempi,ierr,mu_local=mu)
          dum = 0.
          radprop(icv,i) = get_cv(cv_type,densi,dum,mu)

       case(20)
          do_radiation_local = .false.  ! so that p_cgs is interpreted as total (gas + radiation)
          ! pressure and eint also contains gas, rad., and ionisation components
          if (use_var_comp) then
             call calc_temp_and_ene(ieos,rho_cgs,p_cgs,eint,tempi,ierr,mu_local=eos_vars(imu,i),&
                                    X_local=eos_vars(iX,i),Z_local=eos_vars(iZ,i),radhydro=do_radiation_local)
          else
             call calc_temp_and_ene(ieos,rho_cgs,p_cgs,eint,tempi,ierr,mu_local=eos_vars(imu,i),radhydro=do_radiation_local)
          endif
          u_gasrec = (eint-radconst*tempi**4/rho_cgs)/unit_ergg
          radprop(icv,i) = u_gasrec/tempi
       case default
          call fatal('setstar_utils','only ieos=2,20 are supported with radiation')
       end select

       rad(iradxi,i) = radxi_from_Trad(densi,tempi)
       eos_vars(itemp,i) = tempi
       vxyzu(4,i) = radprop(icv,i)*tempi
    else
       select case(ieos)
       case(23) ! Tillotson
          vxyzu(4,i) = 1.5*polyk_in
       case(16) ! Shen EoS
          vxyzu(4,i) = initialtemp
       case default ! Recalculate eint and temp for each particle according to EoS
          if (use_var_comp) then
             call calc_temp_and_ene(ieos,rho_cgs,p_cgs,eni,tempi,ierr,&
                                    mu_local=eos_vars(imu,i),X_local=eos_vars(iX,i),Z_local=eos_vars(iZ,i))
          elseif (eos_outputs_mu(ieos)) then
             call calc_temp_and_ene(ieos,rho_cgs,p_cgs,eni,tempi,ierr,mu_local=eos_vars(imu,i))
          else
             call calc_temp_and_ene(ieos,rho_cgs,p_cgs,eni,tempi,ierr)
          endif
          eos_vars(itemp,i) = tempi
          vxyzu(4,i) = eni / unit_ergg
       end select
    endif
 enddo
 !$omp end parallel do

end subroutine set_star_thermalenergy

!-----------------------------------------------------------------------
!+
!  Solve for u, T profiles given rho, P
!+
!-----------------------------------------------------------------------
subroutine solve_uT_profiles(eos_type,r,den,pres,Xfrac,Yfrac,regrid_core,temp,en,mu)
 use eos,     only:get_mean_molecular_weight,calc_temp_and_ene
 use physcon, only:radconst,Rg
 integer,           intent(in)    :: eos_type
 real,              intent(in)    :: r(:),den(:),pres(:),Xfrac(:),Yfrac(:)
 logical,           intent(in)    :: regrid_core
 real, allocatable, intent(inout) :: temp(:),en(:),mu(:)
 integer             :: i,ierr
 real                :: guessene,tempi,eni

 if (regrid_core) then  ! lengths of rho, P arrays have changed, so need to resize temp, en, and mu
    deallocate(temp,en,mu)
    allocate(temp(size(r)),en(size(r)),mu(size(r)))
 endif

 do i=1,size(r)
    mu(i) = get_mean_molecular_weight(Xfrac(i),1.-Xfrac(i)-Yfrac(i))  ! only used in u, T calculation if ieos==2,12
    if (i==1) then
       guessene = 1.5*pres(i)/den(i)  ! initial guess
       tempi = min((3.*pres(i)/radconst)**0.25, pres(i)*mu(i)/(den(i)*Rg)) ! guess for temperature
    else
       guessene = en(i-1)
       tempi = temp(i-1)
    endif
    call calc_temp_and_ene(eos_type,den(i),pres(i),eni,tempi,ierr,guesseint=guessene,mu_local=mu(i))  ! for ieos==20, mu is outputted here
    en(i) = eni
    temp(i) = tempi
 enddo

end subroutine solve_uT_profiles

!-----------------------------------------------------------------------
!+
!  default filename for the softened stellar profile, from the name of
!  the input profile (without any directory): profile.data gives
!  profile-softened.data, other names have .softened appended
!+
!-----------------------------------------------------------------------
function get_softened_profile_filename(input_profile) result(filename)
 use fileutils, only:basename
 character(len=*), intent(in) :: input_profile
 character(len=120) :: filename
 integer :: n

 filename = basename(input_profile)
 n = len_trim(filename)
 if (n > 5 .and. filename(max(n-4,1):n) == '.data') then
    filename = filename(1:n-5)//'-softened.data'
 else
    filename = trim(filename)//'.softened'
 endif

end function get_softened_profile_filename

!-----------------------------------------------------------------------
!+
!  Find the edge of the dense (hydrogen-depleted) core of a giant star
!  from the steepest density gradient, dln(rho)/dln(r), inside 0.1 R,
!  taking the outer edge of the steep region where the gradient rises
!  back above dln(rho)/dln(r) = -5. The gradient is measured on a uniform
!  grid in ln(r) so that it does not depend on the zoning of the profile.
!  A core is only identified if the steepest gradient is below -5
!  (cf. ~ -0.3 inside 0.1 R for a main sequence star), otherwise
!  r_core = m_core = 0. This agrees with the X < 0.1 boundary for
!  profiles with composition.
!+
!-----------------------------------------------------------------------
subroutine find_dense_core(r,m,den,r_core,m_core)
 use table_utils, only:yinterp
 real, intent(in)  :: r(:),m(:),den(:)
 real, intent(out) :: r_core,m_core
 integer, parameter :: ngrid = 200
 real,    parameter :: grad_crit = -5.
 real, allocatable  :: lnr(:),lnrho(:)
 real    :: lnr1,dlnr,grad(ngrid-1),lnrho_grid(ngrid)
 integer :: i,imin,n

 r_core = 0.
 m_core = 0.
 n = size(r)
 lnr   = log(max(r,tiny(0.)))
 lnrho = log(max(den,tiny(0.)))
 lnr1  = log(max(r(1),1.e-6*r(n)))
 dlnr  = (log(0.1*r(n)) - lnr1)/(ngrid - 1)
 do i=1,ngrid
    lnrho_grid(i) = yinterp(lnrho,lnr,lnr1 + (i-1)*dlnr)
 enddo
 grad = (lnrho_grid(2:ngrid) - lnrho_grid(1:ngrid-1))/dlnr
 imin = minloc(grad,1)
 if (grad(imin) >= grad_crit) return

 ! move to the outer edge of the steep region
 i = imin
 do while (i < ngrid-1 .and. grad(i) < grad_crit)
    i = i + 1
 enddo
 r_core = exp(lnr1 + (i-1)*dlnr)
 m_core = yinterp(m,r,r_core)

end subroutine find_dense_core

!-----------------------------------------------------------------------
!+
!  Suggest values for the radius and mass of a softened stellar core,
!  given a stellar profile (any consistent units, ordered centre to
!  surface) and the number of particles used to resolve the envelope:
!
!   r_core,m_core : radius and mass of the dense core (see find_dense_core),
!                   which should be contained within the sink particle
!   rcore_min     : smallest rcore for which the sink softening length
!                   hsoft = rcore/radkern exceeds the gas smoothing length
!                   at rcore, i.e. rcore > radkern*h(rcore)
!   rcore_rec     : as above but with hsoft > 2h(rcore)
!
!  The smoothing length is estimated from the density of the original
!  profile, h = hfact*(mpart/rho)^(1/3). The local density is used because
!  almost all of the mass inside rcore is given to the sink when the core
!  is softened, so the softened region holds only of order
!  4/3 pi rcore^3 rho(rcore) of gas. rcore/h first decreases outwards from
!  the dense core then increases through the envelope, so the radii are
!  found by moving inwards from where rcore/h peaks.
!+
!-----------------------------------------------------------------------
subroutine get_core_softening_hints(r,m,den,np,hfact,r_core,m_core,rcore_min,rcore_rec)
 use kernel,  only:radkern
 real,    intent(in)  :: r(:),m(:),den(:),hfact
 integer, intent(in)  :: np
 real,    intent(out) :: r_core,m_core,rcore_min,rcore_rec
 real, allocatable :: ratio(:)
 real    :: mpart
 integer :: i,n,ipeak

 n = size(r)
 call find_dense_core(r,m,den,r_core,m_core)

 ! particle mass if np particles are used to resolve the gas outside the core
 mpart = (m(n) - m_core)/max(np,1)

 ! rcore/(radkern*h) at each radius outside the dense core
 allocate(ratio(n))
 ratio = 0.
 where (r > r_core .and. den > 0.) ratio = r/(radkern*hfact*(mpart/den)**(1./3.))
 ipeak = maxloc(ratio,1)

 rcore_min = r(ipeak)
 rcore_rec = r(ipeak)
 do i=ipeak,1,-1
    if (ratio(i) < 1.) exit
    rcore_min = r(i)
 enddo
 do i=ipeak,1,-1
    if (ratio(i) < 2.) exit
    rcore_rec = r(i)
 enddo
 rcore_min = max(rcore_min,r_core)
 rcore_rec = max(rcore_rec,rcore_min)

end subroutine get_core_softening_hints

!-----------------------------------------------------------------------
!+
!  print advice on choosing the core softening radius and core mass,
!  checking the values supplied by the user (input in cgs units)
!+
!-----------------------------------------------------------------------
subroutine print_core_softening_hints(r,m,den,np,hfact,rcore,mcore)
 use physcon, only:solarr,solarm
 use io,      only:num
 real,    intent(in) :: r(:),m(:),den(:),hfact,rcore,mcore
 integer, intent(in) :: np
 real    :: r_core,m_core,rcore_min,rcore_rec

 call get_core_softening_hints(r,m,den,np,hfact,r_core,m_core,rcore_min,rcore_rec)

 write(*,"(/,a,i0,a)",advance='no') ' Core softening (np = ',np,'):'
 if (m_core > 0.) then
    write(*,"(a)") ' dense core inside r = '//num(r_core/solarr)//' Rsun, try mcore ~ '//num(m_core/solarm)
    write(*,"(a)") '  need rcore > '//num(rcore_min/solarr)//' Rsun for sink softening > gas smoothing length h'
 else
    write(*,"(a)") ' no dense core found'
 endif

 if (rcore > 0. .and. rcore < r_core) write(*,"(/,a)") &
    ' WARNING: rcore = '//num(rcore/solarr)//' Rsun is inside dense core (r < '//num(r_core/solarr)//' Rsun)'
 if (rcore > 0. .and. rcore < rcore_min) write(*,"(/,a,i0,a)") &
    ' WARNING: rcore = '//num(rcore/solarr)//' Rsun unresolved with np = ',np,' particles, need rcore > '//num(rcore_min/solarr)
 if (mcore > m(size(m))) write(*,"(/,a)") ' WARNING: mcore > stellar mass'

end subroutine print_core_softening_hints

!-----------------------------------------------------------------------
!+
!  Check the convective stability of the stellar profile, given the
!  density and pressure (in code units) that the star is set up with,
!  where the thermal energy is obtained by inverting the equation of state.
!
!  A fluid parcel displaced adiabatically outwards from zone i to i+1 keeps
!  its composition and changes density according to dln(rho) = dln(P)/Gamma_1.
!  It is buoyant (unstable) if it ends up less dense than its surroundings,
!  i.e. if
!
!     dln(rho)/dln(P) < 1/Gamma_1
!
!  This is the Ledoux criterion, since any composition gradient is included
!  in the density of the surroundings, reducing to the Schwarzschild
!  criterion for uniform composition. Gamma_1 = cs^2 rho/P is obtained from
!  the sound speed returned by the equation of state.
!
!  Unstable regions will convect. Since nothing removes heat at the stellar
!  surface, the convected energy can only be released by expanding the star.
!+
!-----------------------------------------------------------------------
subroutine check_convective_stability(ieos,r,den,pres,mu,Xfrac,Yfrac,mtab,mcore,nunstable,rinner)
 use units,   only:udist,umass
 use physcon, only:solarr,solarm
 use eos,     only:get_gamma1
 use io,      only:num
 integer, intent(in)  :: ieos
 real,    intent(in)  :: r(:),den(:),pres(:),mu(:),Xfrac(:),Yfrac(:)
 real,    intent(in)  :: mtab(:)  ! mass coordinate, used for reporting only
 real,    intent(in)  :: mcore    ! mass of sink particle core (not included in mtab)
 integer, intent(out) :: nunstable
 real,    intent(out) :: rinner
 real, parameter   :: tol = 1.e-3  ! ignore regions where the buoyancy of a rising parcel is < tol
 real, allocatable :: gamma1(:),buoy(:),superad(:)
 real    :: buoyreg,underdense,superad_max
 integer :: i,n,istart,ierr,nskip,nempty
 character(len=13) :: criterion

 nunstable = 0
 rinner = 0.
 n = size(r)
 if (n < 3) return

 ! get Gamma_1 in each zone, skipping empty zones and zones where this cannot be computed
 allocate(gamma1(n),buoy(n),superad(n))
 nskip  = 0
 nempty = 0
 do i=1,n
    if (den(i) <= 0. .or. pres(i) <= 0.) then
       gamma1(i) = -1.
       nempty = nempty + 1
       cycle
    endif
    gamma1(i) = get_gamma1(ieos,den(i),pres(i),ierr,mu_local=mu(i),&
                           X_local=Xfrac(i),Z_local=1.-Xfrac(i)-Yfrac(i))
    if (ierr /= 0) nskip = nskip + 1
 enddo
 if (nskip + nempty==n) then
    write(*,"(/,a,i0,/)") ' Convective stability check not available for ieos=',ieos
    return
 endif

 ! log density contrast acquired by a parcel displaced adiabatically from zone i to i+1:
 ! positive if the parcel is less dense than its surroundings (buoyant). Dividing by the
 ! change in ln P gives the superadiabaticity, 1/Gamma_1 - dln rho/dln P, which equals
 ! nabla - nabla_ad for an ideal gas of uniform composition
 buoy = 0.
 superad = 0.
 do i=1,n-1
    if (min(gamma1(i),gamma1(i+1),den(i),den(i+1),pres(i),pres(i+1)) <= 0.) cycle
    if (pres(i+1)==pres(i)) cycle
    buoy(i) = log(den(i+1)/den(i)) - log(pres(i+1)/pres(i))/(0.5*(gamma1(i) + gamma1(i+1)))
    superad(i) = buoy(i)/abs(log(pres(i+1)/pres(i)))
    if (r(i+1) < r(i)) buoy(i) = -buoy(i)  ! profile ordered surface to centre
    if (r(i+1) < r(i)) superad(i) = -superad(i)
 enddo

 criterion = 'Schwarzschild'
 if (maxval(mu(1:n)) - minval(mu(1:n)) > epsilon(0.)*maxval(mu(1:n))) criterion = 'Ledoux'
 write(*,"(/,a,i0,a)") ' Checking convective stability of stellar profile with ieos=',ieos,&
                       ' ('//trim(criterion)//' criterion)'
 if (nempty > 0) write(*,"(a,i0,a)") ' WARNING: ',nempty,' zones have zero density or pressure'
 if (nskip > 0) write(*,"(a,i0,a)") ' WARNING: Ignoring invalid Gamma_1 in ',nskip,' zones'

 ! find contiguous unstable regions and report the significant ones
 istart = 0
 do i=1,n
    if (i < n .and. buoy(min(i,n-1)) > 0.) then
       if (istart==0) istart = i
    elseif (istart > 0) then
       buoyreg = sum(buoy(istart:i-1))
       if (buoyreg > tol) then
          nunstable = nunstable + 1
          if (nunstable==1) then
             rinner = min(r(istart),r(i))
             write(*,"(/,a)") ' WARNING: star is convectively unstable:'
          endif
          ! fraction by which a parcel rising through the whole region is less dense than its surroundings
          underdense  = 1. - exp(-buoyreg)
          superad_max = maxval(superad(istart:i-1))
          write(*,"(a)") '   '//num(min(r(istart),r(i))*udist/solarr,5)//' < r < '//&
             num(max(r(istart),r(i))*udist/solarr,5)//' Rsun ('//&
             num((min(mtab(istart),mtab(i))+mcore)*umass/solarm,5)//' < m < '//&
             num((max(mtab(istart),mtab(i))+mcore)*umass/solarm,5)//' Msun): '//&
             num(100.*underdense,4)//'% underdense -> '//trim(adjustl(instability_label(underdense)))
       endif
       istart = 0
    endif
 enddo

 if (nunstable > 0) then
    write(*,"(/,a)") ' WARNING: star will expand after relaxation for r > '//num(rinner*udist/solarr)//' Rsun'
    write(*,"(a)") '          (unstable regions convect, and since nothing removes heat at the stellar surface', &
                   '           the convected energy can only be released by expanding the star)'
    write(*,"(a)")
 else
    write(*,"(a,/)") ' Stellar profile is convectively stable'
 endif

end subroutine check_convective_stability

!-----------------------------------------------------------------------
!+
!  describe the strength of a convective instability by the fraction
!  by which a parcel rising through the unstable region becomes less
!  dense than its surroundings
!+
!-----------------------------------------------------------------------
function instability_label(underdense) result(label)
 real, intent(in) :: underdense
 character(len=19) :: label  ! right-aligned so that the columns line up

 if (underdense < 0.05) then
    label = 'mildly unstable'
 elseif (underdense < 0.3) then
    label = 'moderately unstable'
 else
    label = 'strongly unstable'
 endif
 label = adjustr(label)

end function instability_label

end module setstar_utils
