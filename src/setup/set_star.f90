!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setstar
!
! Utility routines for mapping 1D stars into 3D spheres
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: eos, eos_piecewise, extern_densprofile, io, part, physcon,
!   radiation_utils, rho_profile, setsoftenedcore, setup_params, sortutils,
!   spherical, table_utils, unifdis, units
!
 use extern_densprofile, only:nrhotab
 implicit none
 !
 ! Index of setup options
 !
 integer, parameter, public :: nprofile_opts =  7 ! maximum number of initial configurations
 integer, parameter, public :: iuniform   = 1
 integer, parameter, public :: ipoly      = 2
 integer, parameter, public :: ifromfile  = 3
 integer, parameter, public :: ikepler    = 4
 integer, parameter, public :: imesa      = 5
 integer, parameter, public :: ibpwpoly   = 6
 integer, parameter, public :: ievrard    = 7

 character(len=*), parameter, public :: profile_opt(nprofile_opts) = &
    (/'Uniform density profile     ', &
      'Polytrope                   ', &
      'Density vs r from ascii file', &
      'KEPLER star from file       ', &
      'MESA star from file         ', &
      'Piecewise polytrope         ', &
      'Evrard collapse             '/)

 public :: read_star_profile
 public :: set_star_density
 public :: set_star_composition
 public :: set_star_thermalenergy
 public :: set_stellar_core

 private

 integer, parameter :: ng_max = nrhotab

contains

!-------------------------------------------------------------------------------
!+
!  read stellar profile (density, pressure, temperature, composition) from file
!+
!-------------------------------------------------------------------------------
subroutine read_star_profile(iprofile,ieos,input_profile,gamma,polyk,ui_coef,r,den,pres,temp,en,mtab,&
                             Xfrac,Yfrac,mu,npts,rmin,Rstar,Mstar,rhocentre,&
                             isoftcore,isofteningopt,rcore,hsoft,outputfilename)
 use extern_densprofile, only:read_rhotab_wrapper
 use eos_piecewise,      only:get_dPdrho_piecewise
 use eos,                only:get_mean_molecular_weight,calc_temp_and_ene,init_eos
 use rho_profile,        only:rho_uniform,rho_polytrope,rho_piecewise_polytrope,rho_evrard,&
                              read_mesa,read_kepler_file,write_profile,func
 use setsoftenedcore,    only:set_softened_core
 use io,                 only:fatal
 use physcon,            only:kb_on_mh,radconst
 integer,           intent(in)    :: iprofile,ieos
 character(len=*),  intent(in)    :: input_profile,outputfilename
 real,              intent(in)    :: ui_coef
 real,              intent(inout) :: gamma,polyk
 real, allocatable, intent(out)   :: r(:),den(:),pres(:),temp(:),en(:),mtab(:)
 real, allocatable, intent(out)   :: Xfrac(:),Yfrac(:),mu(:)
 integer,           intent(out)   :: npts
 real,              intent(inout) :: rmin,Rstar,Mstar,rhocentre,hsoft
 integer,           intent(in)    :: isoftcore,isofteningopt
 real,              intent(in)    :: rcore
 integer :: ierr,i
 logical :: calc_polyk,iexist
 real    :: eni,tempi,guessene
 procedure(func), pointer :: get_dPdrho
 !
 ! set up tabulated density profile
 !
 calc_polyk = .true.
 allocate(r(ng_max),den(ng_max),pres(ng_max),temp(ng_max),en(ng_max),mtab(ng_max))

 print "(/,a,/)",' Using '//trim(profile_opt(iprofile))
 select case(iprofile)
 case(ipoly)
    call rho_polytrope(gamma,polyk,Mstar,r,den,npts,rhocentre,calc_polyk,Rstar)
    rmin = r(1)
    pres = polyk*den**gamma
 case(ifromfile)
    call read_rhotab_wrapper(trim(input_profile),ng_max,r,den,npts,&
                             polyk,gamma,rhocentre,Mstar,iexist,ierr)
    if (.not.iexist) call fatal('setup','density file does not exist')
    if (ierr > 0)    call fatal('setup','error in reading density file')
    rmin = r(1)
    pres = polyk*den**gamma
 case(ibpwpoly)
    get_dPdrho => get_dPdrho_piecewise
    call rho_piecewise_polytrope(r,den,rhocentre,Mstar,get_dPdrho,npts,ierr)
    if (ierr == 1) call fatal('setup','ng_max is too small')
    if (ierr == 2) call fatal('setup','failed to converge to a self-consistent density profile')
    rmin  = r(1)
    Rstar = r(npts)
    pres = polyk*den**gamma
 case(imesa)
    deallocate(r,den,pres,temp,en,mtab)
    if (isoftcore > 0) then
       call read_mesa(input_profile,den,r,pres,mtab,en,temp,Xfrac,Yfrac,Mstar,ierr,cgsunits=.true.)
       allocate(mu(size(den)))
       mu = 0.
       if (ierr /= 0) call fatal('setup','error in reading stellar profile from'//trim(input_profile))
       call set_softened_core(isoftcore,isofteningopt,r,den,pres,mtab,Xfrac,Yfrac,ierr) ! sets mcore, rcore
       hsoft = 0.5 * rcore
       ! solve for temperature and energy profile
       do i=1,size(r)
          mu(i) = get_mean_molecular_weight(Xfrac(i),1.-Xfrac(i)-Yfrac(i))  ! only used in u, T calculation if ieos==2,12
          if (i==1) then
             guessene = 1.5*pres(i)/den(i)  ! initial guess
             tempi = min((3.*pres(i)/radconst)**0.25, pres(i)*mu(i)/(den(i)*kb_on_mh)) ! guess for temperature
          else
             guessene = en(i-1)
             tempi = temp(i-1)
          endif
          call calc_temp_and_ene(ieos,den(i),pres(i),eni,tempi,ierr,guesseint=guessene,mu_local=mu(i))  ! for ieos==20, mu is outputted here
          en(i) = eni
          temp(i) = tempi
       enddo
       call write_profile(outputfilename,mtab,pres,temp,r,den,en,Xfrac,Yfrac,mu=mu)
       ! now read the softened profile instead
       call read_mesa(outputfilename,den,r,pres,mtab,en,temp,Xfrac,Yfrac,Mstar,ierr)
    else
       call read_mesa(input_profile,den,r,pres,mtab,en,temp,Xfrac,Yfrac,Mstar,ierr)
    endif
    if (ierr==1) call fatal('set_star',trim(input_profile)//' does not exist')
    if (ierr==2) call fatal('set_star','insufficient data points read from file')
    if (ierr==3) call fatal('set_star','too many data points; increase ng')
    if (ierr /= 0) call fatal('set_star','error in reading stellar profile from'//trim(input_profile))
    npts = size(den)
    rmin  = r(1)
    Rstar = r(npts)
 case(ikepler)
    call read_kepler_file(trim(input_profile),ng_max,npts,r,den,pres,temp,en,Mstar,ierr)
    if (ierr==1) call fatal('set_star',trim(input_profile)//' does not exist')
    if (ierr==2) call fatal('set_star','insufficient data points read from file')
    if (ierr==3) call fatal('set_star','too many data points; increase ng')
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

end subroutine read_star_profile

!-------------------------------------------------------------------------------
!+
!  set up particles according to the desired density profile
!+
!-------------------------------------------------------------------------------
subroutine set_star_density(lattice,id,master,rmin,Rstar,Mstar,hfact,&
                            npts,den,r,npart,npartoftype,massoftype,xyzh,use_exactN,np,mask)
 use part,         only:set_particle_type,igas
 use spherical,    only:set_sphere
 use setup_params, only:rhozero,npart_total
 use unifdis,      only:mask_prototype
 use physcon,      only:pi
 character(len=*), intent(in)    :: lattice
 integer,          intent(in)    :: id,master,npts,np
 real,             intent(in)    :: rmin,Rstar,Mstar,hfact
 real,             intent(in)    :: den(npts),r(npts)
 integer,          intent(inout) :: npart,npartoftype(:)
 real,             intent(inout) :: xyzh(:,:),massoftype(:)
 logical,          intent(in)    :: use_exactN
 procedure(mask_prototype) :: mask
 integer :: nx,i,ntot
 real :: vol_sphere,psep
 !
 ! place particles in sphere
 !
 vol_sphere  = 4./3.*pi*Rstar**3
 nx          = int(np**(1./3.))
 psep        = vol_sphere**(1./3.)/real(nx)
 call set_sphere(lattice,id,master,rmin,Rstar,psep,hfact,npart,xyzh, &
                 rhotab=den(1:npts),rtab=r(1:npts),nptot=npart_total, &
                 exactN=use_exactN,np_requested=np,mask=mask)
 !
 ! get total particle number across all MPI threads
 ! and use this to set the particle mass
 !
 ntot = int(npart_total,kind=(kind(ntot)))
 massoftype(igas) = Mstar/ntot
 !
 ! set particle type as gas particles
 !
 npartoftype(igas) = npart   ! npart is number on this thread only
 do i=1,npart
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
subroutine set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,ihsoft,mcore,hsoft)
 integer, intent(out) :: nptmass
 real, intent(out)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(in)     :: mcore,hsoft
 integer              :: n,ihsoft

 nptmass                = 1
 n                      = nptmass
 xyzmh_ptmass(:,n)      = 0. ! zero all quantities by default
 xyzmh_ptmass(4,n)      = mcore
 xyzmh_ptmass(ihsoft,n) = hsoft
 vxyz_ptmass(:,n)       = 0.

end subroutine set_stellar_core

!-----------------------------------------------------------------------
!+
!  Set the composition, if variable composition is used
!+
!-----------------------------------------------------------------------
subroutine set_star_composition(use_var_comp,use_mu,npart,xyzh,Xfrac,Yfrac,mu,mtab,Mstar,eos_vars)
 use part,        only:iorder=>ll,iX,iZ,imu  ! borrow the unused linklist array for the sort
 use sortutils,   only:find_rank,r2func
 use table_utils, only:yinterp
 logical, intent(in)  :: use_var_comp,use_mu
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(in)  :: Xfrac(:),Yfrac(:),mu(:),mtab(:),Mstar
 real,    intent(out) :: eos_vars(:,:)
 real :: ri,massri
 integer :: i

 ! this does NOT work with MPI
 call find_rank(npart,r2func,xyzh(1:3,:),iorder)
 do i = 1,npart
    ri  = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    massri = Mstar * real(iorder(i)-1) / real(npart) ! mass coordinate of particle i
    if (use_var_comp) then
       eos_vars(iX,i) = yinterp(Xfrac,mtab,massri)
       eos_vars(iZ,i) = 1. - eos_vars(iX,i) - yinterp(Yfrac,mtab,massri)
    endif
    if (use_mu) eos_vars(imu,i) = yinterp(mu,mtab,massri)
 enddo

end subroutine set_star_composition

!-----------------------------------------------------------------------
!+
!  Set the thermal energy profile
!+
!-----------------------------------------------------------------------
subroutine set_star_thermalenergy(ieos,den,pres,r,npart,xyzh,vxyzu,rad,eos_vars,relaxed,use_var_comp,initialtemp)
 use part,            only:do_radiation,rhoh,massoftype,igas,itemp,igasP,iX,iZ,imu,iradxi
 use eos,             only:equationofstate,calc_temp_and_ene,gamma,gmw
 use radiation_utils, only:ugas_from_Tgas,radE_from_Trad
 use table_utils,     only:yinterp
 use units,           only:unit_density,unit_ergg,unit_pressure
 integer, intent(in)    :: ieos,npart
 real,    intent(in)    :: den(:), pres(:), r(:)  ! density and pressure tables
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:),eos_vars(:,:),rad(:,:)
 logical, intent(in)    :: relaxed,use_var_comp
 real,    intent(in)    :: initialtemp
 integer :: eos_type,i,ierr
 real    :: xi,yi,zi,hi,presi,densi,tempi,eni,ri,p_on_rhogas,spsoundi
 real    :: rho_cgs,p_cgs

 if (do_radiation) then
    eos_type=12  ! Calculate temperature from both gas and radiation pressure
 else
    eos_type=ieos
 endif
 do i = 1,npart
    if (relaxed) then
       hi = xyzh(4,i)
       densi = rhoh(hi,massoftype(igas))
       presi = eos_vars(igasP,i)  ! retrieve pressure from relax_star calculated with the fake (ieos=2) internal energy
    else
       !  Interpolate density and pressure from table
       ri    = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
       densi = yinterp(den,r,ri)
       presi = yinterp(pres,r,ri)
    endif

    select case(ieos)
    case(16) ! Shen EoS
       vxyzu(4,i) = initialtemp
    case(15) ! Helmholtz EoS
       xi    = xyzh(1,i)
       yi    = xyzh(2,i)
       zi    = xyzh(3,i)
       tempi = initialtemp
       call equationofstate(ieos,p_on_rhogas,spsoundi,densi,xi,yi,zi,tempi,eni)
       vxyzu(4,i) = eni
       eos_vars(itemp,i) = initialtemp
    case default ! Recalculate eint and temp for each particle according to EoS
       rho_cgs = densi*unit_density
       p_cgs = presi*unit_pressure
       if (use_var_comp) then
          call calc_temp_and_ene(eos_type,rho_cgs,p_cgs,eni,tempi,ierr,&
                                 mu_local=eos_vars(imu,i),X_local=eos_vars(iX,i),Z_local=eos_vars(iZ,i))
       else
          call calc_temp_and_ene(eos_type,rho_cgs,p_cgs,eni,tempi,ierr)
       endif
       if (do_radiation) then
          vxyzu(4,i) = ugas_from_Tgas(tempi,gamma,gmw)
          rad(iradxi,i) = radE_from_Trad(tempi)/densi
       else
          vxyzu(4,i) = eni / unit_ergg
       endif
       eos_vars(itemp,i) = tempi
    end select
 enddo

end subroutine set_star_thermalenergy

end module setstar
