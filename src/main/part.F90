!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module part
!
! This module contains main particle data for the code and
!  functions used to query particle properties not stored.
!
!  Basically this module defines any quantity that is
!  stored on the particles and defines how that storage
!  is implemented. Thus any routine which requires knowledge
!  of the specifics of this storage should be placed here.
!  (for example, any routine that copies all of the variables
!   stored on a given particle).
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, dim, dtypekdtree, io, krome_user, mpiutils
!
 use dim, only:ndim,maxp,maxsts,ndivcurlv,ndivcurlB,maxvxyzu,maxalpha,&
               maxptmass,maxdvdx,nsinkproperties,mhd,maxmhd,maxBevol,&
               maxp_h2,maxindan,nabundances,periodic,ind_timesteps,&
               maxgrav,ngradh,maxtypes,gravity,maxp_dustfrac,&
               use_dust,use_dustgrowth,lightcurve,maxlum,nalpha,maxmhdni, &
               maxp_growth,maxdusttypes,maxdustsmall,maxdustlarge, &
               maxphase,maxgradh,maxan,maxdustan,maxmhdan,maxneigh,maxprad,maxp_nucleation,&
               maxTdust,store_dust_temperature,use_krome,maxp_krome, &
               do_radiation,gr,maxgr,maxgran,n_nden_phantom,do_nucleation,&
               inucleation,itau_alloc,itauL_alloc,use_apr,apr_maxlevel,maxp_apr,maxptmassgr
 use dtypekdtree, only:kdnode
#ifdef KROME
 use krome_user, only: krome_nmols
#endif
 implicit none
!
!--basic storage needed for read/write of particle data
!

 real,         allocatable :: xyzh(:,:)
 real,         allocatable :: xyzh_soa(:,:)
 real,         allocatable :: vxyzu(:,:)
 real(kind=4), allocatable :: alphaind(:,:)
 real(kind=4), allocatable :: divcurlv(:,:)
 real(kind=4), allocatable :: divcurlB(:,:)
 real,         allocatable :: Bevol(:,:)
 real,         allocatable :: Bxyz(:,:)
 character(len=*), parameter :: xyzh_label(4) = (/'x','y','z','h'/)
 character(len=*), parameter :: vxyzu_label(4) = (/'vx','vy','vz','u '/)
 character(len=*), parameter :: Bxyz_label(3) = (/'Bx','By','Bz'/)
 character(len=*), parameter :: Bevol_label(4) = (/'Bx/rho','By/rho','Bz/rho','psi   '/)
 character(len=*), parameter :: alphaind_label(3) = (/'alpha   ','alphaloc','div_a   '/)

!
!--tracking particle IDs
!
 integer(kind=8)              :: norig
 integer(kind=8), allocatable :: iorig(:)
!
!--storage of dust properties
!
 real :: grainsize(maxdusttypes)
 real :: graindens(maxdusttypes)
!
!--storage of dust growth properties
!
 real, allocatable :: dustprop(:,:)    !- mass and intrinsic density
 real, allocatable :: dustgasprop(:,:) !- gas related quantites interpolated on dust particles (see Force.F90)
 real, allocatable :: VrelVf(:)
 character(len=*), parameter :: dustprop_label(2) = (/'grainmass','graindens'/)
 character(len=*), parameter :: dustgasprop_label(4) = (/'csound','rhogas','St    ','dv    '/)
 character(len=*), parameter :: VrelVf_label = 'Vrel/Vfrag'

 !- porosity
 integer, allocatable :: dragreg(:)    !- drag regime
 real, allocatable    :: mprev(:)      !- previous mass
 real, allocatable    :: filfac(:)     !- filling factor
 real, allocatable    :: filfacprev(:) !- previous filling factor needed for minimum St condition
 real, allocatable    :: probastick(:) !-probabily of sticking, when bounce is on
 character(len=*), parameter :: filfac_label = 'filfac'
 !- options
 logical, public             :: this_is_a_test  = .false.
 logical, public             :: this_is_a_flyby = .false.
!
!--storage in divcurlv
!
 integer, parameter :: idivv = 1
 integer, parameter :: icurlvx = 2
 integer, parameter :: icurlvy = 3
 integer, parameter :: icurlvz = 4
 character(len=*), parameter :: divcurlv_label(4) = &
   (/'divv  ','curlvx','curlvy','curlvz'/)
!
!--storage in divcurlB
!
 integer, parameter :: idivB = 1
 integer, parameter :: icurlBx = 2
 integer, parameter :: icurlBy = 3
 integer, parameter :: icurlBz = 4
 character(len=*), parameter :: divcurlB_label(4) = &
   (/'divB  ','curlBx','curlBy','curlBz'/)
!
!--velocity gradients
!
 real(kind=4), allocatable :: dvdx(:,:)
 character(len=*), parameter :: dvdx_label(9) = &
   (/'dvxdx','dvxdy','dvxdz', &
     'dvydx','dvydy','dvydz', &
     'dvzdx','dvzdy','dvzdz'/)
!
!--H2 chemistry
!
 integer, parameter :: ih2ratio  = 1 ! ratio of H2 to H
 integer, parameter :: iHI       = 2 ! HI abundance
 integer, parameter :: iproton   = 3 ! proton abundance
 integer, parameter :: ielectron = 4 ! electron abundance
 integer, parameter :: iCO       = 5 ! CO abundance
 real, allocatable :: abundance(:,:)
!
!--KROME chemistry
!
#ifdef KROME
 character(len=16)  :: abundance_label(krome_nmols)
#else
 character(len=*), parameter :: abundance_label(nabundances) = &
   (/'h2ratio','  abHIq','  abhpq','   abeq','   abco'/)
#endif
 character(len=*), parameter :: abundance_meaning(nabundances) = &
      (/'ratio of molecular to atomic Hydrogen       ',&
        'nHI/nH:  fraction of neutral atomic Hydrogen',&
        'nHII/nH: fraction of ionised Hydrogen (HII) ',&
        'ne/nH:   fraction of electrons              ',&
        'nCO/nH:  fraction of Carbon Monoxide        '/)

!
!--make a public krome_nmols variable to avoid ifdefs elsewhere
!
#ifndef KROME
 integer, parameter :: krome_nmols = 0
#endif

!
!--eos_variables
!
 real, allocatable  :: eos_vars(:,:)
 integer, parameter :: igasP = 1, &
                       ics   = 2, &
                       itemp = 3, &
                       imu   = 4, &
                       iX    = 5, &
                       iZ    = 6, &
                       igamma = 7, &
                       maxeosvars = 7
 character(len=*), parameter :: eos_vars_label(maxeosvars) = &
    (/'pressure   ','sound speed','temperature','mu         ','H fraction ','metallicity','gamma      '/)

!
!--energy_variables
!
 integer, public :: ien_type
 integer, public, parameter :: ien_entropy = 1, &
                               ien_etotal  = 2, &
                               ien_entropy_s = 3
!
!--one-fluid dust (small grains)
!
 real, allocatable :: dustfrac(:,:)
 real, allocatable :: dustevol(:,:)
 real, allocatable :: deltav(:,:,:)
 character(len=*), parameter :: dustfrac_label(maxdusttypes) = 'dustfrac'
 character(len=*), parameter :: tstop_label(maxdusttypes) = 'tstop'
 character(len=*), parameter :: deltav_label(3) = &
   (/'deltavx','deltavy','deltavz'/)
!
!--General relativity
!
 real, allocatable :: pxyzu(:,:) !pxyzu(maxvxyzu,maxgr)
 character(len=*), parameter :: pxyzu_label(4) = (/'px     ','py     ','pz     ','entropy'/)
 real, allocatable :: dens(:) !dens(maxgr)
 real, allocatable :: metrics(:,:,:,:) !metrics(0:3,0:3,2,maxgr)
 real, allocatable :: metricderivs(:,:,:,:) !metricderivs(0:3,0:3,3,maxgr)
 real, allocatable :: tmunus(:,:,:) !tmunus(0:3,0:3,maxgr)
 real, allocatable :: sqrtgs(:) ! sqrtg(maxgr)
!
!--sink particles in General relativity
!
 real, allocatable :: pxyzu_ptmass(:,:) !pxyz_ptmass(maxvxyzu,maxgr)
 real, allocatable :: metrics_ptmass(:,:,:,:) !metrics(0:3,0:3,2,maxgr)
 real, allocatable :: metricderivs_ptmass(:,:,:,:) !metricderivs(0:3,0:3,3,maxgr)
!
!--sink particles
!
 integer, parameter :: ihacc  = 5  ! accretion radius
 integer, parameter :: ihsoft = 6  ! softening radius
 integer, parameter :: imacc  = 7  ! accreted mass
 integer, parameter :: ispinx = 8  ! spin angular momentum x
 integer, parameter :: ispiny = 9  ! spin angular momentum y
 integer, parameter :: ispinz = 10 ! spin angular momentum z
 integer, parameter :: i_tlast = 11 ! time of last injection
 integer, parameter :: ilum   = 12 ! luminosity
 integer, parameter :: iTeff  = 13 ! effective temperature
 integer, parameter :: iReff  = 14 ! effective radius
 integer, parameter :: imloss = 15 ! mass loss rate
 integer, parameter :: imdotav = 16 ! accretion rate average
 integer, parameter :: i_mlast = 17 ! accreted mass of last time
 integer, parameter :: imassenc = 18 ! mass enclosed in sink softening radius
 integer, parameter :: iJ2 = 19      ! 2nd gravity moment due to oblateness
 integer, parameter :: irstrom = 20  ! Stromgren radius of the stars (icreate_sinks == 2)
 integer, parameter :: irateion = 21 ! Ionisation rate of the stars (log)(icreate_sinks == 2)
 integer, parameter :: itbirth = 22  ! birth time of the new sink
 integer, parameter :: ndptmass = 13 ! number of properties to conserve after accretion phase or merge
 integer, allocatable :: sf_ptmass(:,:) ! star form prop 1 : type (1 sink ,2 star, 3 dead sink ), 2 : number of seeds
 real,    allocatable :: xyzmh_ptmass(:,:)
 real,    allocatable :: vxyz_ptmass(:,:)
 real,    allocatable :: fxyz_ptmass(:,:),fxyz_ptmass_sinksink(:,:),fsink_old(:,:)
 real,    allocatable :: dsdt_ptmass(:,:),dsdt_ptmass_sinksink(:,:)
 real,    allocatable :: dptmass(:,:)
 integer :: nptmass = 0   ! zero by default
 real    :: epot_sinksink
 character(len=*), parameter :: xyzmh_ptmass_label(nsinkproperties) = &
  (/'x        ','y        ','z        ','m        ','h        ',&
    'hsoft    ','maccreted','spinx    ','spiny    ','spinz    ',&
    'tlast    ','lum      ','Teff     ','Reff     ','mdotloss ',&
    'mdotav   ','mprev    ','massenc  ','J2       ','Rstrom   ',&
    'rate_ion ','tbirth   '/)
 character(len=*), parameter :: vxyz_ptmass_label(3) = (/'vx','vy','vz'/)
 character(len=*), parameter :: sf_ptmass_label(2) = (/'type  ','nseed '/)
!
!--self-gravity
!
 real(kind=4), allocatable :: poten(:)
!
!--Non-ideal MHD
!
 real, allocatable :: nden_nimhd(:,:)
 real, allocatable :: eta_nimhd(:,:)
 integer, parameter :: iohm  = 1 ! eta_ohm
 integer, parameter :: ihall = 2 ! eta_hall
 integer, parameter :: iambi = 3 ! eta_ambi
 integer, parameter :: iion  = 4 ! ionisation fraction
 character(len=*), parameter :: eta_nimhd_label(4) = (/'eta_{OR}','eta_{HE}','eta_{AD}','ne/n    '/)
!
!-- Ray tracing : optical depth
!
 real, allocatable :: tau(:)
 real, allocatable :: tau_lucy(:)
!
!--Dust formation - theory of moments
!
 real, allocatable :: dust_temp(:)
 integer, parameter :: n_nucleation = 10
 real, allocatable :: nucleation(:,:)
 ! please note that in nucleation, we save the *normalized* moments, i.e. \hat{K}_i = K_i/n<H>
 character(len=*), parameter :: nucleation_label(n_nucleation) = &
       (/'Jstar','K0   ','K1   ','K2   ','K3   ',&
         'mu   ','gamma','S    ','kappa','alphw'/)
 integer, parameter :: idJstar = 1, &
                       idK0    = 2, &
                       idK1    = 3, &
                       idK2    = 4, &
                       idK3    = 5, &
                       idmu    = 6, &
                       idgamma = 7, &
                       idsat   = 8, & !for logging
                       idkappa = 9, & !for logging
                       idalpha = 10   !for logging
!
!--KROME variables
!
 real, allocatable :: T_gas_cool(:)
!
!--radiation hydro, evolved quantities (which have time derivatives)
!
 real, allocatable  :: rad(:,:)
 integer, parameter :: iradxi = 1, &
                       maxirad= 1
 character(len=*), parameter :: rad_label(maxirad) = (/'xi'/)
!
!--additional quantities required for radiation hydro (not evolved)
!
 real, allocatable  :: radprop(:,:)
 integer, parameter :: ifluxx = 1, &
                       ifluxy = 2, &
                       ifluxz = 3, &
                       ikappa = 4, &
                       ithick = 5, &
                       inumph = 6, &
                       ivorcl = 7, &
                       iradP  = 8, &
                       ilambda = 9, &
                       iedd   = 10, &
                       icv    = 11, &
                       maxradprop = 11
 character(len=*), parameter :: radprop_label(maxradprop) = &
    (/'radFx ','radFy ','radFz ','kappa ','thick ','numph ','vorcl ','radP  ','lambda','edd   ','cv    '/)
!
!--lightcurves
!
 real(kind=4), allocatable :: luminosity(:)
!

!--APR - we need these arrays whether we use apr or not
!
 integer(kind=1), allocatable :: apr_level(:)
 integer(kind=1), allocatable :: apr_level_soa(:)
!

!-- Regularisation algorithm allocation
!
 integer(kind=1), allocatable :: nmatrix(:,:)    ! adjacency matrix used to construct each groups

 integer, allocatable :: group_info(:,:) ! array storing group id/idx of each group comp/boundary idx/bin comp id
 integer, parameter   :: igarg  = 1 ! idx of the particle member of a group
 integer, parameter   :: igcum  = 2 ! cumulative sum of the indices to find the starting and ending point of a group
 integer, parameter   :: igid   = 3 ! id of the group, correspond to the root of the group in the dfs/union find construction
 integer, parameter   :: icomp  = 4 ! id of the binary companion if it exists, otherwise equal to the id

 real, allocatable    :: bin_info(:,:) ! array storing important orbital parameters and quantities of each binary
 integer, parameter   :: isemi = 1 ! semi major axis
 integer, parameter   :: iecc  = 2 ! eccentricity
 integer, parameter   :: iapo  = 3 ! apocenter
 integer, parameter   :: iorb  = 4 ! orbital period
 integer, parameter   :: ipert = 5 ! perturbation
 integer, parameter   :: ikap  = 6 ! kappa slow down


 ! needed for group identification and sorting
 integer  :: n_group = 0
 integer  :: n_ingroup = 0
 integer  :: n_sing = 0
 ! Gradient of the time transformation function
 real, allocatable :: gtgrad(:,:)
 !
!-- Regularisation algorithm allocation
!
 logical, allocatable :: isionised(:)
!
!--derivatives (only needed if derivs is called)
!
 real, allocatable         :: fxyzu(:,:)
 real, allocatable         :: fxyz_drag(:,:)
 real, allocatable         :: fxyz_dragold(:,:)
 real, allocatable         :: dBevol(:,:)
 real(kind=4), allocatable :: divBsymm(:)
 real, allocatable         :: fext(:,:)
 real, allocatable         :: ddustevol(:,:)
 real, allocatable         :: ddustprop(:,:) !--grainmass is the only prop that evolves for now
 real, allocatable         :: drad(:,:)
 character(len=*), parameter :: ddustprop_label(2) = (/' dm/dt ','drho/dt'/)
!
!--storage associated with/dependent on timestepping
!
 real, allocatable   :: vpred(:,:)
 real, allocatable   :: ppred(:,:)
 real, allocatable   :: dustpred(:,:)
 real, allocatable   :: Bpred(:,:)
 real, allocatable   :: dustproppred(:,:)
 real, allocatable   :: filfacpred(:)
 real, allocatable   :: radpred(:,:)
 integer(kind=1), allocatable :: ibin(:)
 integer(kind=1), allocatable :: ibin_old(:)
 integer(kind=1), allocatable :: ibin_wake(:)
 real(kind=4),    allocatable :: dt_in(:)
 real,            allocatable :: twas(:)

 integer(kind=1), allocatable    :: iphase(:)
 integer(kind=1), allocatable    :: iphase_soa(:)
 logical, public    :: all_active = .true.

 real(kind=4), allocatable :: gradh(:,:)
 real, allocatable         :: tstop(:,:)
!
!--storage associated with link list
!  (used for dead particle list also)
!
 integer, allocatable :: ll(:)
 real    :: dxi(ndim) ! to track the extent of the particles
!
!--particle belong
!
 integer, allocatable :: ibelong(:)
!
!--super time stepping
!
 integer(kind=1), allocatable :: istsactive(:)
 integer(kind=1), allocatable :: ibin_sts(:)
!
!--size of the buffer required for transferring particle
!  information between MPI threads
!
 integer, parameter, private :: usedivcurlv = min(ndivcurlv,1)
 integer, parameter :: ipartbufsize = 129

 real            :: hfact,Bextx,Bexty,Bextz
 integer         :: npart
 integer(kind=8) :: ntot
 integer         :: ideadhead = 0

 integer         :: npartoftype(maxtypes)
 integer(kind=8) :: npartoftypetot(maxtypes)
 real            :: massoftype(maxtypes),aprmassoftype(maxtypes,apr_maxlevel)

 integer :: ndustsmall,ndustlarge,ndusttypes
!
!--labels for each type
!  NOTE: If new particle is added, and it is allowed to be accreted onto
!        a sink particle, add it to the list in 'is_accretable' below
!  NOTE: set_boundaries_to_active = .true., but will be set to .false. at the
!        end of densityiterate.  This will allow boundary particles to always be
!         initialised (even on restarts where not all arrays, e.g. gradh,
!         are not saved)
!
 integer, parameter :: igas        = 1
 integer, parameter :: iboundary   = 3
 integer, parameter :: istar       = 4
 integer, parameter :: idarkmatter = 5
 integer, parameter :: ibulge      = 6
 integer, parameter :: idust       = 7
 integer, parameter :: idustlast   = idust + maxdustlarge - 1
 integer, parameter :: idustbound  = idustlast + 1
 integer, parameter :: idustboundl = idustbound + maxdustlarge - 1
 integer, parameter :: iunknown    = 0
 logical            :: set_boundaries_to_active = .true.
 integer :: i
 character(len=7), dimension(maxtypes), parameter :: &
   labeltype = (/'gas    ','empty  ','bound  ','star   ','darkm  ','bulge  ', &
                 ('dust   ', i=idust,idustlast),('dustbnd',i=idustbound,idustboundl)/)
!
!--generic interfaces for routines
!
 interface hrho
  module procedure hrho4,hrho8,hrho4_pmass,hrho8_pmass,hrhomixed_pmass
 end interface hrho
 interface get_ntypes
  module procedure get_ntypes_i4, get_ntypes_i8
 end interface get_ntypes

 private :: hrho4,hrho8,hrho4_pmass,hrho8_pmass,hrhomixed_pmass
 private :: get_ntypes_i4,get_ntypes_i8

contains

subroutine allocate_part
 use allocutils, only:allocate_array

 call allocate_array('xyzh', xyzh, 4, maxp)
 call allocate_array('xyzh_soa', xyzh_soa, maxp, 4)
 call allocate_array('vxyzu', vxyzu, maxvxyzu, maxp)
 call allocate_array('alphaind', alphaind, nalpha, maxalpha)
 call allocate_array('divcurlv', divcurlv, ndivcurlv, maxp)
 call allocate_array('dvdx', dvdx, 9, maxp)
 call allocate_array('divcurlB', divcurlB, ndivcurlB, maxp)
 call allocate_array('Bevol', Bevol, maxBevol, maxmhd)
 call allocate_array('apr_level',apr_level,maxp_apr)
 call allocate_array('apr_level_soa',apr_level_soa,maxp_apr)
 call allocate_array('Bxyz', Bxyz, 3, maxmhd)
 call allocate_array('iorig', iorig, maxp)
 call allocate_array('dustprop', dustprop, 2, maxp_growth)
 call allocate_array('dustgasprop', dustgasprop, 4, maxp_growth)
 call allocate_array('VrelVf', VrelVf, maxp_growth)
 call allocate_array('eosvars', eos_vars, maxeosvars, maxan)
 call allocate_array('dustfrac', dustfrac, maxdusttypes, maxp_dustfrac)
 call allocate_array('dustevol', dustevol, maxdustsmall, maxp_dustfrac)
 call allocate_array('ddustevol', ddustevol, maxdustsmall, maxdustan)
 call allocate_array('ddustprop', ddustprop, 2, maxp_growth)
 call allocate_array('dragreg', dragreg, maxp_growth)
 call allocate_array('filfac', filfac, maxp_growth)
 call allocate_array('mprev', mprev, maxp_growth)
 call allocate_array('probastick', probastick, maxp_growth)
 call allocate_array('filfacprev', filfacprev, maxp_growth)
 call allocate_array('deltav', deltav, 3, maxdustsmall, maxp_dustfrac)
 call allocate_array('pxyzu', pxyzu, maxvxyzu, maxgr)
 call allocate_array('dens', dens, maxgr)
 call allocate_array('metrics', metrics, 4, 4, 2, maxgr)
 call allocate_array('metricderivs', metricderivs, 4, 4, 3, maxgr)
 call allocate_array('tmunus', tmunus, 4, 4, maxgr)
 call allocate_array('sqrtgs', sqrtgs, maxgr)
 call allocate_array('pxyzu_ptmass', pxyzu_ptmass, maxvxyzu, maxptmassgr)
 call allocate_array('metrics_ptmass', metrics_ptmass, 4, 4, 2, maxptmassgr)
 call allocate_array('metricderivs_ptmass', metricderivs_ptmass, 4, 4, 3, maxptmassgr)
 call allocate_array('xyzmh_ptmass', xyzmh_ptmass, nsinkproperties, maxptmass)
 call allocate_array('vxyz_ptmass', vxyz_ptmass, 3, maxptmass)
 call allocate_array('fxyz_ptmass', fxyz_ptmass, 4, maxptmass)
 call allocate_array('fxyz_ptmass_sinksink', fxyz_ptmass_sinksink, 4, maxptmass)
 call allocate_array('fsink_old', fsink_old, 4, maxptmass)
 call allocate_array('dptmass', dptmass, ndptmass,maxptmass)
 call allocate_array('sf_ptmass', sf_ptmass, 2, maxptmass)
 call allocate_array('dsdt_ptmass', dsdt_ptmass, 3, maxptmass)
 call allocate_array('dsdt_ptmass_sinksink', dsdt_ptmass_sinksink, 3, maxptmass)
 call allocate_array('poten', poten, maxgrav)
 call allocate_array('nden_nimhd', nden_nimhd, n_nden_phantom, maxmhdni)
 call allocate_array('eta_nimhd', eta_nimhd, 4, maxmhdni)
 call allocate_array('luminosity', luminosity, maxlum)
 call allocate_array('fxyzu', fxyzu, maxvxyzu, maxan)
 call allocate_array('fxyz_drag', fxyz_drag, 3, maxdustan)
 call allocate_array('fxyz_dragold', fxyz_dragold, 3, maxdustan)
 call allocate_array('dBevol', dBevol, maxBevol, maxmhdan)
 call allocate_array('divBsymm', divBsymm, maxmhdan)
 call allocate_array('fext', fext, 3, maxan)
 call allocate_array('vpred', vpred, maxvxyzu, maxan)
 call allocate_array('ppred', ppred, maxvxyzu, maxgran)
 call allocate_array('dustpred', dustpred, maxdustsmall, maxdustan)
 call allocate_array('Bpred', Bpred, maxBevol, maxmhdan)
 call allocate_array('dustproppred', dustproppred, 2, maxp_growth)
 call allocate_array('filfacpred', filfacpred, maxp_growth)
 call allocate_array('dust_temp',dust_temp,maxTdust)
 call allocate_array('rad', rad, maxirad, maxprad)
 call allocate_array('radpred', radpred, maxirad, maxprad)
 call allocate_array('drad', drad, maxirad, maxprad)
 call allocate_array('radprop', radprop, maxradprop, maxprad)
 call allocate_array('ibin', ibin, maxindan)
 call allocate_array('ibin_old', ibin_old, maxindan)
 call allocate_array('ibin_wake', ibin_wake, maxindan)
 call allocate_array('dt_in', dt_in, maxindan)
 call allocate_array('twas', twas, maxindan)
 call allocate_array('iphase', iphase, maxphase)
 call allocate_array('iphase_soa', iphase_soa, maxphase)
 call allocate_array('gradh', gradh, ngradh, maxgradh)
 call allocate_array('tstop', tstop, maxdusttypes, maxan)
 call allocate_array('ll', ll, maxan)
 call allocate_array('ibelong', ibelong, maxp)
 call allocate_array('istsactive', istsactive, maxsts)
 call allocate_array('ibin_sts', ibin_sts, maxsts)
 call allocate_array('nucleation', nucleation, n_nucleation, maxp_nucleation*inucleation)
 call allocate_array('tau', tau, maxp*itau_alloc)
 call allocate_array('tau_lucy', tau_lucy, maxp*itauL_alloc)
 if (use_krome) then
    call allocate_array('abundance', abundance, krome_nmols, maxp_krome)
 else
    call allocate_array('abundance', abundance, nabundances, maxp_h2)
 endif
 call allocate_array('T_gas_cool', T_gas_cool, maxp_krome)
 call allocate_array('group_info', group_info, 4, maxptmass)
 call allocate_array('bin_info', bin_info, 6, maxptmass)
 call allocate_array("nmatrix", nmatrix, maxptmass, maxptmass)
 call allocate_array("gtgrad", gtgrad, 3, maxptmass)
 call allocate_array('isionised', isionised, maxp)


end subroutine allocate_part

subroutine deallocate_part

 if (allocated(xyzh))     deallocate(xyzh)
 if (allocated(xyzh_soa)) deallocate(xyzh_soa)
 if (allocated(vxyzu))    deallocate(vxyzu)
 if (allocated(alphaind)) deallocate(alphaind)
 if (allocated(divcurlv)) deallocate(divcurlv)
 if (allocated(dvdx))     deallocate(dvdx)
 if (allocated(divcurlB)) deallocate(divcurlB)
 if (allocated(Bevol))    deallocate(Bevol)
 if (allocated(Bxyz))     deallocate(Bxyz)
 if (allocated(iorig))    deallocate(iorig)
 if (allocated(dustprop)) deallocate(dustprop)
 if (allocated(dustgasprop))  deallocate(dustgasprop)
 if (allocated(VrelVf))       deallocate(VrelVf)
 if (allocated(abundance))    deallocate(abundance)
 if (allocated(eos_vars))     deallocate(eos_vars)
 if (allocated(dustfrac))     deallocate(dustfrac)
 if (allocated(dustevol))     deallocate(dustevol)
 if (allocated(ddustevol))    deallocate(ddustevol)
 if (allocated(ddustprop))    deallocate(ddustprop)
 if (allocated(dragreg))      deallocate(dragreg)
 if (allocated(filfac))       deallocate(filfac)
 if (allocated(mprev))        deallocate(mprev)
 if (allocated(filfacprev))   deallocate(filfacprev)
 if (allocated(probastick))   deallocate(probastick)
 if (allocated(deltav))       deallocate(deltav)
 if (allocated(pxyzu))        deallocate(pxyzu)
 if (allocated(dens))         deallocate(dens)
 if (allocated(metrics))      deallocate(metrics)
 if (allocated(metricderivs)) deallocate(metricderivs)
 if (allocated(tmunus))       deallocate(tmunus)
 if (allocated(sqrtgs))       deallocate(sqrtgs)
 if (allocated(pxyzu_ptmass)) deallocate(pxyzu_ptmass)
 if (allocated(metrics_ptmass))  deallocate(metrics_ptmass)
 if (allocated(metricderivs_ptmass))  deallocate(metricderivs_ptmass)
 if (allocated(xyzmh_ptmass)) deallocate(xyzmh_ptmass)
 if (allocated(vxyz_ptmass))  deallocate(vxyz_ptmass)
 if (allocated(fxyz_ptmass))  deallocate(fxyz_ptmass)
 if (allocated(fxyz_ptmass_sinksink)) deallocate(fxyz_ptmass_sinksink)
 if (allocated(fsink_old))    deallocate(fsink_old)
 if (allocated(dptmass))      deallocate(dptmass)
 if (allocated(sf_ptmass)) deallocate(sf_ptmass)
 if (allocated(dsdt_ptmass))  deallocate(dsdt_ptmass)
 if (allocated(dsdt_ptmass_sinksink)) deallocate(dsdt_ptmass_sinksink)
 if (allocated(poten))        deallocate(poten)
 if (allocated(nden_nimhd))   deallocate(nden_nimhd)
 if (allocated(eta_nimhd))    deallocate(eta_nimhd)
 if (allocated(luminosity))   deallocate(luminosity)
 if (allocated(fxyzu))        deallocate(fxyzu)
 if (allocated(fxyz_drag))    deallocate(fxyz_drag)
 if (allocated(fxyz_dragold)) deallocate(fxyz_dragold)
 if (allocated(dBevol))       deallocate(dBevol)
 if (allocated(divBsymm))     deallocate(divBsymm)
 if (allocated(fext))         deallocate(fext)
 if (allocated(vpred))        deallocate(vpred)
 if (allocated(ppred))        deallocate(ppred)
 if (allocated(dustpred))     deallocate(dustpred)
 if (allocated(Bpred))        deallocate(Bpred)
 if (allocated(dustproppred)) deallocate(dustproppred)
 if (allocated(filfacpred))   deallocate(filfacpred)
 if (allocated(ibin))         deallocate(ibin)
 if (allocated(ibin_old))     deallocate(ibin_old)
 if (allocated(ibin_wake))    deallocate(ibin_wake)
 if (allocated(dt_in))        deallocate(dt_in)
 if (allocated(twas))         deallocate(twas)
 if (allocated(nucleation))   deallocate(nucleation)
 if (allocated(tau))          deallocate(tau)
 if (allocated(tau_lucy))     deallocate(tau_lucy)
 if (allocated(T_gas_cool))   deallocate(T_gas_cool)
 if (allocated(dust_temp))    deallocate(dust_temp)
 if (allocated(rad))          deallocate(rad,radpred,drad,radprop)
 if (allocated(iphase))       deallocate(iphase)
 if (allocated(iphase_soa))   deallocate(iphase_soa)
 if (allocated(gradh))        deallocate(gradh)
 if (allocated(tstop))        deallocate(tstop)
 if (allocated(ll))           deallocate(ll)
 if (allocated(ibelong))      deallocate(ibelong)
 if (allocated(istsactive))   deallocate(istsactive)
 if (allocated(ibin_sts))     deallocate(ibin_sts)
 if (allocated(apr_level))    deallocate(apr_level)
 if (allocated(apr_level_soa)) deallocate(apr_level_soa)
 if (allocated(group_info))   deallocate(group_info)
 if (allocated(bin_info))     deallocate(bin_info)
 if (allocated(nmatrix))      deallocate(nmatrix)
 if (allocated(gtgrad))       deallocate(gtgrad)
 if (allocated(isionised))    deallocate(isionised)

end subroutine deallocate_part

!----------------------------------------------------------------
!+
!  initialise all particle arrays to default values
!  i.e. reset the simulation
!+
!----------------------------------------------------------------
subroutine init_part
 integer :: i

 npart = 0
 nptmass = 0
 npartoftype(:) = 0
 npartoftypetot(:) = 0
 massoftype(:)  = 0.
 isionised(:) = .false.
!--initialise point mass arrays to zero
 xyzmh_ptmass = 0.
 vxyz_ptmass  = 0.
 dsdt_ptmass  = 0.
 sf_ptmass = 0

 ! initialise arrays not passed to setup routine to zero
 if (mhd) then
    Bevol = 0.
    Bxyz  = 0.
    Bextx = 0.
    Bexty = 0.
    Bextz = 0.
 endif
 if (maxphase > 0) iphase = 0 ! phases not set
 if (maxalpha==maxp)  alphaind = 0.
 if (ndivcurlv > 0) divcurlv = 0.
 if (maxdvdx==maxp) dvdx = 0.
 if (ndivcurlB > 0) divcurlB = 0.
 if (maxgrav > 0) poten = 0.
 if (use_dust) then
    dustfrac = 0.
    dustevol = 0.
 endif
 ndustsmall = 0
 ndustlarge = 0
 if (lightcurve) luminosity = 0.
 if (use_apr) apr_level = 1 ! this is reset if the simulation is to derefine
 if (do_radiation) then
    rad(:,:) = 0.
    radprop(:,:) = 0.
    radprop(ikappa,:) = huge(0.) ! set opacity to infinity
    radprop(ithick,:) = 1.       ! optically thick, i.e. use diffusion approximation
 endif
!
!--initialise chemistry arrays if this has been compiled
!  (these may be altered by the specific setup routine)
!
 if (maxp_h2==maxp) then
    abundance(:,:)   = 0.
    abundance(iHI,:) = 1.  ! assume all atomic hydrogen initially
 endif
 if (store_dust_temperature) dust_temp = -1.
 !-- Initialise dust properties to zero
 if (use_dustgrowth) then
    dustprop(:,:)    = 0.
    dustgasprop(:,:) = 0.
    VrelVf(:)        = 0.
 endif
 if (ind_timesteps) then
    ibin(:)       = 0
    ibin_old(:)   = 0
    ibin_wake(:)  = 0
    dt_in(:)      = 0.
    twas(:)       = 0.
 endif

 ideadhead = 0
!
!--Initialise particle id's
!
!$omp parallel do default(none) &
!$omp shared(iorig,maxp) &
!$omp private(i)
 do i = 1,maxp
    iorig(i) = i
 enddo
!$omp end parallel do
 norig = maxp

end subroutine init_part

!----------------------------------------------------------------
!+
!  this function determines the mass of the particle
!  where use_gas == .not. (maxphase==maxp)
!  currently used only in readwrite_dump
!+
!----------------------------------------------------------------
real function get_pmass(i,use_gas)
 integer, intent(in) :: i
 logical, intent(in) :: use_gas

 if (use_gas) then
    get_pmass = massoftype(igas)
 else
    if (iphase(i) /= 0) then
       get_pmass = massoftype(iamtype(iphase(i)))
    else
       get_pmass = massoftype(igas)
    endif
 endif

end function get_pmass

!
!----------------------------------------------------------------
!+
!  this function gives rho as a function of h
!  (used by output routines to get rho given that we only
!   store h)
!+
!----------------------------------------------------------------
pure real function rhoh(hi,pmassi)
 real, intent(in) :: hi,pmassi

 rhoh = pmassi*(hfact/abs(hi))**3

end function rhoh

!----------------------------------------------------------------
!+
!  this function gives dh/drho as a function of h
!+
!----------------------------------------------------------------
pure real function dhdrho(hi,pmassi)
 real, intent(in) :: hi,pmassi
 real :: rhoi

 rhoi   = rhoh(hi,pmassi)
 dhdrho = -hi/(3.*rhoi)

end function dhdrho
!----------------------------------------------------------------
!+
!  this subroutine does both of the above
!  (routine is an optimisation - also returns divisions)
!+
!----------------------------------------------------------------
pure subroutine rhoanddhdrho(hi,hi1,rhoi,rho1i,dhdrhoi,pmassi)
 real,         intent(in)  :: hi, pmassi
 real(kind=8), intent(out) :: hi1
 real,         intent(out) :: rhoi
 real,         intent(out) :: rho1i,dhdrhoi
 real, parameter :: third = 1./3.

 hi1 = 1./abs(hi)
 rhoi = real(pmassi*(hfact*hi1)**3)
 rho1i = 1./rhoi
 dhdrhoi = -third*hi*rho1i

end subroutine rhoanddhdrho

!----------------------------------------------------------------
!+
!  this function gives h as a function of rho
!+
!----------------------------------------------------------------
real(kind=4) function hrho4(rhoi)
 real(kind=4), intent(in) :: rhoi

 hrho4 = real(hfact*(massoftype(1)/abs(rhoi))**(1./3.),kind=4)

end function hrho4

real(kind=8) function hrho8(rhoi)
 real(kind=8), intent(in) :: rhoi

 hrho8 = hfact*(massoftype(1)/abs(rhoi))**(1.d0/3.d0)

end function hrho8

real(kind=4) function hrho4_pmass(rhoi,pmassi)
 real(kind=4), intent(in) :: rhoi,pmassi

 hrho4_pmass = real(hfact*(pmassi/abs(rhoi))**(1./3.),kind=4)

end function hrho4_pmass

real(kind=8) function hrho8_pmass(rhoi,pmassi)
 real(kind=8), intent(in) :: rhoi,pmassi

 hrho8_pmass = hfact*(pmassi/abs(rhoi))**(1.d0/3.d0)

end function hrho8_pmass

real(kind=8) function hrhomixed_pmass(rhoi,pmassi)
 real(kind=8), intent(in) :: rhoi
 real(kind=4), intent(in) :: pmassi

 hrhomixed_pmass = hfact*(pmassi/abs(rhoi))**(1.d0/3.d0)

end function hrhomixed_pmass

!------------------------------------------------------------------------
!+
!  Query function to see if any sink particles have a non-zero luminosity
!+
!------------------------------------------------------------------------
logical function sinks_have_luminosity(nptmass,xyzmh_ptmass)
 integer, intent(in) :: nptmass
 real, intent(in) :: xyzmh_ptmass(:,:)

 sinks_have_luminosity = any(xyzmh_ptmass(iTeff,1:nptmass) > 0. .and. &
                              xyzmh_ptmass(ilum,1:nptmass) > 0.)

end function sinks_have_luminosity

!------------------------------------------------------------------------
!+
!  Query function to see if any sink particles have heating
!+
!------------------------------------------------------------------------
logical function sinks_have_heating(nptmass,xyzmh_ptmass)
 integer, intent(in) :: nptmass
 real, intent(in) :: xyzmh_ptmass(:,:)

 sinks_have_heating = any(xyzmh_ptmass(iTeff,1:nptmass) <= 0. .and. &
                              xyzmh_ptmass(ilum,1:nptmass) > 0.)

end function sinks_have_heating

!------------------------------------------------------------------------
!+
!  Query function to see if any sink particles have heating
!+
!------------------------------------------------------------------------
logical function sink_has_heating(xyzmh_ptmassi)
 real, intent(in) :: xyzmh_ptmassi(:)

 sink_has_heating = xyzmh_ptmassi(iTeff) <= 0. .and. &
                              xyzmh_ptmassi(ilum) > 0.

end function sink_has_heating

!----------------------------------------------------------------
!+
!  query function returning whether or not a particle is dead
!  (currently indicated by having a negative smoothing length)
!+
!----------------------------------------------------------------
logical function isdead(i)
 integer, intent(in) :: i

 ! h = 0 indicates a dead particle
 if (abs(xyzh(4,i)) < tiny(xyzh)) then
    isdead = .true.
 else
    isdead = .false.
 endif

end function isdead

pure logical function isdeadh(hi)
 real, intent(in) :: hi

 ! h = 0 indicates a dead particle
 if (abs(hi) < tiny(hi)) then
    isdeadh = .true.
 else
    isdeadh = .false.
 endif

end function isdeadh

pure logical function isdead_or_accreted(hi)
 real, intent(in) :: hi

 ! h <= 0 indicates either dead or accreted
 if (hi < tiny(hi)) then
    isdead_or_accreted = .true.
 else
    isdead_or_accreted = .false.
 endif

end function isdead_or_accreted

!----------------------------------------------------------------
!+
! routine which kills a particle and adds it to the dead list
!+
!----------------------------------------------------------------
subroutine kill_particle(i,npoftype)
 integer, intent(in) :: i
 integer, intent(inout), optional :: npoftype(:)

 if (i < 1 .or. i > npart) return ! do nothing
 !
 ! WARNING : this routine is *NOT THREAD SAFE *
 !
 ! do not kill particles that are already dead
 ! because this causes endless loop in shuffle_part
 if (.not.isdeadh(xyzh(4,i))) then
    xyzh(4,i) = 0.
    if (present(npoftype)) call remove_particle_from_npartoftype(i,npoftype)
    ll(i) = ideadhead
    ideadhead = i
 endif

end subroutine kill_particle

!----------------------------------------------------------------
!+
!  decrement npartoftype when a particle is killed, according
!  to the type of the particle that was destroyed
!+
!----------------------------------------------------------------
subroutine remove_particle_from_npartoftype(i,npoftype)
 integer, intent(in)    :: i
 integer, intent(inout) :: npoftype(:)
 integer :: itype

 ! get the type so we know how to decrement npartoftype
 if (maxphase==maxp) then
    itype = iamtype(iphase(i))
    if (itype <= 0) itype = igas ! safety check
 else
    itype = igas
 endif
 npoftype(itype) = npoftype(itype) - 1

end subroutine remove_particle_from_npartoftype

!----------------------------------------------------------------
!+
!  recount particle types, useful after particles have moved
!  between MPI tasks
!+
!----------------------------------------------------------------
subroutine recount_npartoftype
 integer :: itype

 npartoftype(:) = 0
 do i=1,npart
    itype = iamtype(iphase(i))
    npartoftype(itype) = npartoftype(itype) + 1
 enddo

end subroutine recount_npartoftype

!----------------------------------------------
!+
!  functions to deconstruct iphase
!  abs(iphase) is the particle type
!  sign(iphase) gives whether it is active/inactive
!+
!----------------------------------------------
pure integer(kind=1) function isetphase(itype,iactive)
 integer, intent(in) :: itype
 logical, intent(in) :: iactive

 if ((set_boundaries_to_active .and.      iamboundary(itype))   .or. &
     (iactive                  .and. .not.iamboundary(itype)) ) then
    isetphase = int(itype,kind=1)
 else
    isetphase = -int(abs(itype),kind=1)
 endif

end function isetphase

pure subroutine get_partinfo(iphasei,isactive,isgas,isdust,itype)
 integer(kind=1), intent(in)  :: iphasei
 logical,         intent(out) :: isactive,isgas,isdust
 integer,         intent(out) :: itype

! isactive = iactive(iphasei)
! itype = iamtype(iphasei)
! isdust = itype==idust

!--inline versions of above (for speed)
 if (iphasei > 0) then
    isactive = .true.
    itype    = iphasei
 else
    isactive = .false.
    itype    = -iphasei
 endif
 isgas = (itype==igas .or. itype==iboundary)
#ifdef DUST
 isdust = ((itype>=idust) .and. (itype<=idustlast))  .or. &
          ((itype>=idustbound) .and. (itype<=idustboundl))
#else
 isdust = .false.
#endif
 !
 ! boundary particles (always inactive unless set to active)
 !
 if (itype==iboundary) then
    if (set_boundaries_to_active) then
       isactive = .true.
       itype = igas
    else
       isactive = .false.
    endif
 elseif (itype>= idustbound .and. itype <= idustboundl) then
    if (set_boundaries_to_active) then
       isactive = .true.
       itype = idust + itype - idustbound
    else
       isactive = .false.
    endif
 endif

end subroutine get_partinfo

pure logical function iactive(iphasei)
 integer(kind=1), intent(in) :: iphasei

 iactive = (iphasei > 0)

end function iactive

pure elemental integer function iamtype(iphasei)
 integer(kind=1), intent(in) :: iphasei

 iamtype = abs(iphasei)

end function iamtype

pure elemental function iamtype_int1(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer(kind=1) :: iamtype_int1

 iamtype_int1 = abs(iphasei)

end function iamtype_int1

pure function iamtype_int11(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer(kind=1) :: iamtype_int11

 iamtype_int11 = abs(iphasei)

end function iamtype_int11

pure elemental logical function iamgas(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer :: itype

 itype = iamtype(iphasei)
 iamgas = int(itype)==igas

end function iamgas

pure elemental logical function iamboundary(itype)
 integer, intent(in) :: itype

 !itype = abs(itype) unnecessary as always called with type, not iphase
 iamboundary = itype==iboundary .or. (itype>=idustbound .and. itype<=idustboundl)

end function iamboundary

pure elemental integer function ibasetype(itype)
 integer, intent(in) :: itype
 !integer :: itype

 ! return underlying (base) type for particle
 !itype = abs(itype)
 if (itype==iboundary) then
    ibasetype = igas                       ! boundary particles are gas
 elseif (itype>=idustbound .and. itype<=idustboundl) then
    ibasetype = idust + (itype-idustbound) ! dust boundaries are dust
 else
    ibasetype = itype                      ! otherwise same as current type
 endif

end function ibasetype

pure elemental logical function iamdust(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer :: itype

 itype = iamtype(iphasei)
 iamdust = ((itype>=idust) .and. (itype<=idustlast))

end function iamdust

pure elemental integer function idusttype(iphasei)
 integer(kind=1), intent(in) :: iphasei

 if (iamdust(iphasei)) then
    idusttype = iamtype(iphasei) - idust + 1
 else
    idusttype = 1
 endif

end function idusttype

pure function get_ntypes_i4(noftype) result(get_ntypes)
 integer, intent(in) :: noftype(:)
 integer :: get_ntypes
 integer :: i

 get_ntypes = 0
 do i=1,size(noftype)
    if (noftype(i) > 0) get_ntypes = i
 enddo

end function get_ntypes_i4

pure function get_ntypes_i8(noftype) result(get_ntypes)
 integer(kind=8), intent(in) :: noftype(:)
 integer :: get_ntypes
 integer :: i

 get_ntypes = 0
 do i=1,size(noftype)
    if (noftype(i) > 0) get_ntypes = i
 enddo

end function get_ntypes_i8

!-----------------------------------------------------------------------
!+
!  Determine if particle is of a type that is accretable
!    Modify the if-statement to include all the types of particles
!    that can be accreted onto the sink
!+
!-----------------------------------------------------------------------
pure logical function is_accretable(itype)
 integer, intent(in)  :: itype

 if (itype==igas .or. (itype>=idust .and. itype<=idustlast)) then
    is_accretable = .true.
 else
    is_accretable = .false.
 endif

end function is_accretable

!----------------------------------------------------------------
!+
!  utility function for setup routines to set initial value
!  of iphase (assumes particle is active)
!+
!----------------------------------------------------------------
subroutine set_particle_type(i,itype)
 use io, only:fatal
 integer, intent(in) :: i,itype

 if (maxphase==maxp) then
    iphase(i) = isetphase(itype,iactive=.true.)
 elseif (itype /= igas) then
    call fatal('set_particle_type','attempt to setup a particle of type > 1, but iphase not allocated')
 endif

end subroutine set_particle_type

!----------------------------------------------------------------
!+
!  utility function to retrieve particle type
!  This routine accesses iphase from the global arrays:
!  hence it is NOT safe to use in parallel loops
!+
!----------------------------------------------------------------
subroutine get_particle_type(i,itype)
 integer, intent(in)  :: i
 integer, intent(out) :: itype

 if (maxphase==maxp) then
    itype = iamtype(iphase(i))
 else
    itype = igas
 endif

end subroutine get_particle_type

!----------------------------------------------------------------
!+
!  utility function to get strain tensor from dvdx array
!+
!----------------------------------------------------------------
pure function strain_from_dvdx(dvdxi) result(strain)
 real, intent(in) :: dvdxi(9)
 real :: strain(6)

 strain(1) = 2.*dvdxi(1)
 strain(2) = dvdxi(2) + dvdxi(4)
 strain(3) = dvdxi(3) + dvdxi(7)
 strain(4) = 2.*dvdxi(5)
 strain(5) = dvdxi(6) + dvdxi(8)
 strain(6) = 2.*dvdxi(9)

end function strain_from_dvdx

!----------------------------------------------------------------
!+
! routine which copies a particle from one location to another
! (prior to a derivs evaluation - so no derivs required)
!+
!----------------------------------------------------------------
subroutine copy_particle(src,dst,new_part)
 integer, intent(in) :: src, dst
 logical, intent(in) :: new_part

 xyzh(:,dst)  = xyzh(:,src)
 vxyzu(:,dst) = vxyzu(:,src)
 fext(:,dst)  = fext(:,src)
 if (mhd) then
    Bevol(:,dst) = Bevol(:,src)
    Bxyz(:,dst)  = Bxyz(:,src)
 endif
 if (do_radiation) then
    rad(:,dst) = rad(:,src)
    radprop(:,dst) = radprop(:,src)
 endif
 if (gr) pxyzu(:,dst) = pxyzu(:,src)
 if (ndivcurlv  > 0) divcurlv(:,dst)  = divcurlv(:,src)
 if (maxalpha ==maxp) alphaind(:,dst) = alphaind(:,src)
 if (maxgradh ==maxp) gradh(:,dst)    = gradh(:,src)
 if (maxphase ==maxp) iphase(dst)   = iphase(src)
 if (maxgrav  ==maxp) poten(dst) = poten(src)
 if (ind_timesteps) then
    ibin(dst)       = ibin(src)
    ibin_old(dst)   = ibin_old(src)
    ibin_wake(dst)  = ibin_wake(src)
    dt_in(dst)      = dt_in(src)
    twas(dst)       = twas(src)
 endif
 if (use_dust) then
    dustfrac(:,dst) = dustfrac(:,src)
    dustevol(:,dst) = dustevol(:,src)
 endif
 if (use_apr) apr_level(dst) = apr_level(src)
 if (maxp_h2==maxp .or. maxp_krome==maxp) abundance(:,dst) = abundance(:,src)
 eos_vars(:,dst) = eos_vars(:,src)
 if (store_dust_temperature) dust_temp(dst) = dust_temp(src)
 if (do_nucleation) nucleation(:,dst) = nucleation(:,src)
 if (itau_alloc == 1) tau(dst) = tau(src)
 if (itauL_alloc == 1) tau_lucy(dst) = tau_lucy(src)

 if (new_part) then
    norig      = norig + 1
    iorig(dst) = norig      ! we are creating a new particle; give it the new ID
 else
    iorig(dst) = iorig(src) ! we are moving the particle within the list; maintain ID
 endif

 return
end subroutine copy_particle

!----------------------------------------------------------------
!+
! routine which copies a particle from one location to another
! (copies everything which is stored on a particle)
!
! Note that link list information CANNOT be copied so link list
! must be rebuilt after a copy operation.
!+
!----------------------------------------------------------------
subroutine copy_particle_all(src,dst,new_part)
 integer, intent(in) :: src,dst
 logical, intent(in) :: new_part

 xyzh(:,dst)  = xyzh(:,src)
 xyzh_soa(dst,:)  = xyzh_soa(src,:)
 vxyzu(:,dst) = vxyzu(:,src)
 if (maxan==maxp) then
    vpred(:,dst) = vpred(:,src)
    fxyzu(:,dst) = fxyzu(:,src)
    fext(:,dst)  = fext(:,src)
 endif
 if (mhd) then
    Bevol(:,dst)  = Bevol(:,src)
    if (maxmhdan==maxp) then
       Bpred(:,dst)  = Bpred(:,src)
       dBevol(:,dst) = dBevol(:,src)
       divBsymm(dst) = divBsymm(src)
    endif
    Bxyz(:,dst)   = Bxyz(:,src)
    if (maxmhdni==maxp) then
       nden_nimhd(:,dst) = nden_nimhd(:,src)
       eta_nimhd(:,dst)  = eta_nimhd(:,src)
    endif
 endif
 if (do_radiation) then
    rad(:,dst) = rad(:,src)
    radpred(:,dst) = radpred(:,src)
    radprop(:,dst) = radprop(:,src)
    drad(:,dst) = drad(:,src)
 endif
 if (gr) then
    pxyzu(:,dst) = pxyzu(:,src)
    if (maxgran==maxp) then
       ppred(:,dst) = ppred(:,src)
    endif
    dens(dst) = dens(src)
 endif

 if (ndivcurlv > 0) divcurlv(:,dst) = divcurlv(:,src)
 if (ndivcurlB > 0) divcurlB(:,dst) = divcurlB(:,src)
 if (maxdvdx ==maxp)  dvdx(:,dst) = dvdx(:,src)
 if (maxalpha ==maxp) alphaind(:,dst) = alphaind(:,src)
 if (maxgradh ==maxp) gradh(:,dst) = gradh(:,src)
 if (maxphase ==maxp) iphase(dst) = iphase(src)
 if (maxphase ==maxp) iphase_soa(dst) = iphase_soa(src)
 if (maxgrav  ==maxp) poten(dst) = poten(src)
 if (maxlum   ==maxp) luminosity(dst) = luminosity(src)
 if (maxindan==maxp) then
    ibin(dst)       = ibin(src)
    ibin_old(dst)   = ibin_old(src)
    ibin_wake(dst)  = ibin_wake(src)
    dt_in(dst)      = dt_in(src)
    twas(dst)       = twas(src)
 endif
 if (use_dust) then
    if (maxp_dustfrac==maxp) dustfrac(:,dst)  = dustfrac(:,src)
    dustevol(:,dst)  = dustevol(:,src)
    if (maxdustan==maxp) then
       dustpred(:,dst)  = dustpred(:,src)
       ddustevol(:,dst) = ddustevol(:,src)
       if (maxdusttypes > 0) tstop(:,dst) = tstop(:,src)
    endif
    deltav(:,:,dst)  = deltav(:,:,src)
    if (maxp_growth==maxp) then
       dustprop(:,dst) = dustprop(:,src)
       ddustprop(:,dst) = ddustprop(:,src)
       dustgasprop(:,dst) = dustgasprop(:,src)
       VrelVf(dst) = VrelVf(src)
       dustproppred(:,dst) = dustproppred(:,src)
       filfacpred(dst) = filfacpred(src)
    endif
    fxyz_drag(:,dst) = fxyz_drag(:,src)
    fxyz_dragold(:,dst) = fxyz_dragold(:,src)

 endif
 if (maxp_h2==maxp .or. maxp_krome==maxp) abundance(:,dst) = abundance(:,src)
 eos_vars(:,dst) = eos_vars(:,src)
 if (store_dust_temperature) dust_temp(dst) = dust_temp(src)
 if (do_nucleation) nucleation(:,dst) = nucleation(:,src)
 if (itau_alloc == 1) tau(dst) = tau(src)
 if (itauL_alloc == 1) tau_lucy(dst) = tau_lucy(src)

 if (use_krome) then
    T_gas_cool(dst)       = T_gas_cool(src)
 endif
 ibelong(dst) = ibelong(src)
 if (maxsts==maxp) then
    istsactive(dst) = istsactive(src)
    ibin_sts(dst) = ibin_sts(src)
 endif
 if (use_apr) then
    apr_level(dst)      = apr_level(src)
    apr_level_soa(dst)  = apr_level_soa(src)
 endif

 if (new_part) then
    norig      = norig + 1
    iorig(dst) = norig      ! we are creating a new particle; give it the new ID
 else
    iorig(dst) = iorig(src) ! we are moving the particle within the list; maintain ID
 endif

 return
end subroutine copy_particle_all

!------------------------------------------------------------------
!+
! routine which reorders the particles according to an input list
!+
!------------------------------------------------------------------
subroutine reorder_particles(iorder,np)
 integer, intent(in) :: iorder(:)
 integer, intent(in) :: np

 integer :: isrc,nbuf
 real    :: xtemp(ipartbufsize)

 do i=1,np
    isrc = iorder(i)

    ! If particle has already been moved
    do while (isrc < i)
       isrc = iorder(isrc)
    enddo

    ! Swap particles around
    call fill_sendbuf(i,xtemp,nbuf)
    call copy_particle_all(isrc,i,.false.)
    call unfill_buffer(isrc,xtemp)

 enddo

end subroutine reorder_particles

!-----------------------------------------------------------------------
!+
!  routine to compactify the list of particles by removing dead
!  particles from the list
!  (could be openMP parallel if we sent in ndead)
!+
!-----------------------------------------------------------------------
subroutine shuffle_part(np)
 use io,  only:fatal
 use dim, only: mpi
 integer, intent(inout) :: np
 integer :: newpart

 do while (ideadhead /= 0)
    newpart = ideadhead
    if (newpart <= np) then
       if (.not.isdead(np)) then
          ! move particle to new position
          call copy_particle_all(np,newpart,.false.)
          ! move ibelong to new position
          if (mpi) ibelong(newpart) = ibelong(np)
          ! update deadhead
          ideadhead = ll(newpart)
       endif
       np = np - 1
    else
       ideadhead = ll(newpart)
    endif
    if (np < 0) call fatal('shuffle','npart < 0')
 enddo

end subroutine shuffle_part

integer function count_dead_particles()
 integer :: i

 i = ideadhead
 count_dead_particles = 0
 do while (i > 0)
    count_dead_particles = count_dead_particles + 1
    i = ll(i)
 enddo

end function count_dead_particles

!-----------------------------------------------------------------------
!+
!  routine to completely remove dead or accreted particles
!  uses the routines above for efficiency
!+
!-----------------------------------------------------------------------
subroutine delete_dead_or_accreted_particles(npart,npoftype)
 integer, intent(inout) :: npart,npoftype(:)
 integer :: i

 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) call kill_particle(i,npoftype)
 enddo
 call shuffle_part(npart)

end subroutine delete_dead_or_accreted_particles

!----------------------------------------------------------------
!+
!   change the position and status of a dead particle
!
!+
!----------------------------------------------------------------
subroutine change_status_pos(npart,x,y,z,h,vx,vy,vz)
 integer, intent(in) :: npart
 real, intent (in) :: x,y,z,h
 real, intent (in) :: vx,vy,vz
 integer  :: i,ix

 ix=0
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       ix=i
       exit
    endif
 enddo

 xyzh(1,ix)=x
 xyzh(2,ix)=y
 xyzh(3,ix)=z
 xyzh(4,ix)=h
 vxyzu(1,ix)=vx
 vxyzu(2,ix)=vy
 vxyzu(3,ix)=vz

end subroutine change_status_pos

!----------------------------------------------------------------
!+
!  pack particle information into a contiguous buffer
!  to send to another processor
!+
!----------------------------------------------------------------
subroutine fill_sendbuf(i,xtemp,nbuf)
 use io,       only:fatal
 use mpiutils, only:fill_buffer
 integer, intent(in)  :: i
 real,    intent(out) :: xtemp(ipartbufsize)
 integer, intent(out) :: nbuf
!
!--package particle information into one simple wrapper
!
 nbuf = 0
!--NB: could use MPI_PACK here...
 if (i > 0) then
    call fill_buffer(xtemp,xyzh(:,i),nbuf)
    call fill_buffer(xtemp,vxyzu(:,i),nbuf)
    call fill_buffer(xtemp,vpred(:,i),nbuf)
    if (gr) then
       call fill_buffer(xtemp,pxyzu(:,i),nbuf)
       call fill_buffer(xtemp,ppred(:,i),nbuf)
       call fill_buffer(xtemp,dens(i),nbuf)
    endif
    call fill_buffer(xtemp,fxyzu(:,i),nbuf)
    call fill_buffer(xtemp,fext(:,i),nbuf)
    if (ndivcurlv > 0) then
       call fill_buffer(xtemp,divcurlv(1,i),nbuf)
    endif
    if (maxalpha==maxp) then
       call fill_buffer(xtemp,alphaind(:,i),nbuf)
    endif
    if (maxgradh==maxp) then
       call fill_buffer(xtemp,gradh(:,i),nbuf)
    endif
    if (mhd) then
       call fill_buffer(xtemp,Bevol(:,i),nbuf)
       call fill_buffer(xtemp,Bpred(:,i),nbuf)
    endif
    if (do_radiation) then
       call fill_buffer(xtemp,rad(:,i),nbuf)
       call fill_buffer(xtemp,radpred(:,i),nbuf)
       call fill_buffer(xtemp,drad(:,i),nbuf)
       call fill_buffer(xtemp,radprop(:,i),nbuf)
    endif
    if (maxphase==maxp) then
       call fill_buffer(xtemp,iphase(i),nbuf)
    endif
    if (use_dust) then
       call fill_buffer(xtemp, dustfrac(:,i),nbuf)
       call fill_buffer(xtemp, dustevol(:,i),nbuf)
       call fill_buffer(xtemp, dustpred(:,i),nbuf)
       if (use_dustgrowth) then
          call fill_buffer(xtemp, dustprop(:,i),nbuf)
          call fill_buffer(xtemp, dustproppred(:,i),nbuf)
          call fill_buffer(xtemp, dustgasprop(:,i),nbuf)
       endif
       call fill_buffer(xtemp,fxyz_drag(:,i),nbuf)
       call fill_buffer(xtemp,fxyz_dragold(:,i),nbuf)
    endif
    if (maxp_h2==maxp .or. maxp_krome==maxp) then
       call fill_buffer(xtemp, abundance(:,i),nbuf)
    endif
    call fill_buffer(xtemp, eos_vars(:,i),nbuf)
    if (store_dust_temperature) then
       call fill_buffer(xtemp, dust_temp(i),nbuf)
    endif
    if (do_nucleation) then
       call fill_buffer(xtemp, nucleation(:,i),nbuf)
    endif
    if (itau_alloc == 1)  call fill_buffer(xtemp, tau(i),nbuf)
    if (itauL_alloc == 1) call fill_buffer(xtemp, tau_lucy(i),nbuf)

    if (maxgrav==maxp) then
       call fill_buffer(xtemp, poten(i),nbuf)
    endif
    if (ind_timesteps) then
       call fill_buffer(xtemp,ibin(i),nbuf)
       call fill_buffer(xtemp,ibin_old(i),nbuf)
       call fill_buffer(xtemp,ibin_wake(i),nbuf)
       call fill_buffer(xtemp,dt_in(i),nbuf)
       call fill_buffer(xtemp,twas(i),nbuf)
    endif
    call fill_buffer(xtemp,iorig(i),nbuf)
    if (use_apr) call fill_buffer(xtemp,apr_level(i),nbuf)
 endif
 if (nbuf > ipartbufsize) call fatal('fill_sendbuf','error: send buffer size overflow',var='nbuf',ival=nbuf)

end subroutine fill_sendbuf

!----------------------------------------------------------------
!+
!  unpack particle information from the send buffer
!  after receiving from another processor
!+
!----------------------------------------------------------------
subroutine unfill_buffer(ipart,xbuf)
 use mpiutils, only:unfill_buf
 integer, intent(in) :: ipart
 real,    intent(in) :: xbuf(ipartbufsize)
 integer :: j

 j = 0
 xyzh(:,ipart)          = unfill_buf(xbuf,j,4)
 vxyzu(:,ipart)         = unfill_buf(xbuf,j,maxvxyzu)
 vpred(:,ipart)         = unfill_buf(xbuf,j,maxvxyzu)
 if (gr) then
    pxyzu(:,ipart)       = unfill_buf(xbuf,j,maxvxyzu)
    ppred(:,ipart)       = unfill_buf(xbuf,j,maxvxyzu)
    dens(ipart)          = unfill_buf(xbuf,j)
 endif
 fxyzu(:,ipart)         = unfill_buf(xbuf,j,maxvxyzu)
 fext(:,ipart)          = unfill_buf(xbuf,j,3)
 if (ndivcurlv > 0) then
    divcurlv(1,ipart)  = real(unfill_buf(xbuf,j),kind=kind(divcurlv))
 endif
 if (maxalpha==maxp) then
    alphaind(:,ipart)   = real(unfill_buf(xbuf,j,nalpha),kind(alphaind))
 endif
 if (maxgradh==maxp) then
    gradh(:,ipart)      = real(unfill_buf(xbuf,j,ngradh),kind(gradh))
 endif
 if (mhd) then
    Bevol(:,ipart)      = real(unfill_buf(xbuf,j,maxBevol),kind=kind(Bevol))
    Bpred(:,ipart)      = real(unfill_buf(xbuf,j,maxBevol),kind=kind(Bevol))
 endif
 if (do_radiation) then
    rad(:,ipart)     = real(unfill_buf(xbuf,j,maxirad))
    radpred(:,ipart) = real(unfill_buf(xbuf,j,maxirad))
    drad(:,ipart)    = real(unfill_buf(xbuf,j,maxirad))
    radprop(:,ipart) = real(unfill_buf(xbuf,j,maxradprop))
 endif
 if (maxphase==maxp) then
    iphase(ipart)       = nint(unfill_buf(xbuf,j),kind=1)
 endif
 if (use_dust) then
    dustfrac(:,ipart)   = unfill_buf(xbuf,j,maxdusttypes)
    dustevol(:,ipart)   = unfill_buf(xbuf,j,maxdustsmall)
    dustpred(:,ipart)   = unfill_buf(xbuf,j,maxdustsmall)
    if (use_dustgrowth) then
       dustprop(:,ipart)       = unfill_buf(xbuf,j,2)
       dustproppred(:,ipart)   = unfill_buf(xbuf,j,2)
       dustgasprop(:,ipart)    = unfill_buf(xbuf,j,4)
    endif
    fxyz_drag(:,ipart)   = unfill_buf(xbuf,j,3)
    fxyz_dragold(:,ipart)   = unfill_buf(xbuf,j,3)
 endif
 if (maxp_h2==maxp .or. maxp_krome==maxp) then
    abundance(:,ipart)  = unfill_buf(xbuf,j,nabundances)
 endif
 eos_vars(:,ipart) = unfill_buf(xbuf,j,maxeosvars)
 if (store_dust_temperature) then
    dust_temp(ipart)    = unfill_buf(xbuf,j)
 endif
 if (do_nucleation) then
    nucleation(:,ipart) = unfill_buf(xbuf,j,n_nucleation)
 endif
 if (itau_alloc == 1)  tau(ipart) = unfill_buf(xbuf,j)
 if (itauL_alloc == 1) tau_lucy(ipart) = unfill_buf(xbuf,j)
 if (maxgrav==maxp) then
    poten(ipart)        = real(unfill_buf(xbuf,j),kind=kind(poten))
 endif
 if (ind_timesteps) then
    ibin(ipart)         = nint(unfill_buf(xbuf,j),kind=1)
    ibin_old(ipart)     = nint(unfill_buf(xbuf,j),kind=1)
    ibin_wake(ipart)    = nint(unfill_buf(xbuf,j),kind=1)
    dt_in(ipart)        = real(unfill_buf(xbuf,j),kind=kind(dt_in))
    twas(ipart)         = unfill_buf(xbuf,j)
 endif
 iorig(ipart)           = nint(unfill_buf(xbuf,j),kind=8)
 if (use_apr) apr_level(ipart) = nint(unfill_buf(xbuf,j),kind=kind(apr_level))

!--just to be on the safe side, set other things to zero
 if (mhd) then
    divBsymm(ipart) = 0.
 endif

end subroutine unfill_buffer

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (rank 2 arrays)
!+
!----------------------------------------------------------------

subroutine copy_array(array,ilist)
 real,    intent(inout) :: array(:,:)
 integer, intent(in)    :: ilist(:)
 real :: arraytemp(size(array(1,:)))
 integer :: i

 do i=1,size(array(:,1))
    arraytemp(:) = array(i,ilist(:))
    array(i,:) = arraytemp
 enddo

 return
end subroutine copy_array

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (real4, rank 2 arrays)
!+
!----------------------------------------------------------------

subroutine copy_arrayr4(array,ilist)
 real(kind=4), intent(inout) :: array(:,:)
 integer,      intent(in)    :: ilist(:)
 real(kind=4) :: arraytemp(size(array(1,:)))
 integer :: i

 do i=1,size(array(:,1))
    arraytemp(:) = array(i,ilist(:))
    array(i,:) = arraytemp
 enddo

 return
end subroutine copy_arrayr4

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (real4, rank 1 arrays)
!+
!----------------------------------------------------------------

subroutine copy_array1(array,ilist)
 real(kind=4), intent(inout) :: array(:)
 integer,      intent(in)    :: ilist(:)
 real(kind=4) :: arraytemp(size(array(:)))

 arraytemp(:) = array(ilist(:))
 array = arraytemp

 return
end subroutine copy_array1

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (int1, rank 1 arrays)
!+
!----------------------------------------------------------------

subroutine copy_arrayint1(iarray,ilist)
 integer(kind=1), intent(inout) :: iarray(:)
 integer,         intent(in)    :: ilist(:)
 integer(kind=1) :: iarraytemp(size(iarray(:)))

 iarraytemp(:) = iarray(ilist(:))
 iarray = iarraytemp

 return
end subroutine copy_arrayint1

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (int8, rank 1 arrays)
!+
!----------------------------------------------------------------

subroutine copy_arrayint8(iarray,ilist)
 integer(kind=8), intent(inout) :: iarray(:)
 integer,         intent(in)    :: ilist(:)
 integer(kind=8) :: iarraytemp(size(iarray(:)))

 iarraytemp(:) = iarray(ilist(:))
 iarray = iarraytemp

 return
end subroutine copy_arrayint8

!----------------------------------------------------------------
!+
!  Delete particles outside of a defined box
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_box(xmin, xmax, ymin, ymax, zmin, zmax)
 real, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax

 integer :: i
 real :: x, y, z, h

 do i=1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    h = xyzh(4,i)
    if (x  <  xmin .or. x  >  xmax .or. y  <  ymin .or. y  >  ymax .or. z  <  zmin .or. z  >  zmax) then
       xyzh(4,i) = -abs(h)
    endif
 enddo
end subroutine delete_particles_outside_box

!----------------------------------------------------------------
!+
!  Delete particles outside (or inside) of a defined sphere
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_sphere(center,radius,np,revert,mytype)
 use io, only:fatal
 real,    intent(in)    :: center(3), radius
 integer, intent(inout) :: np
 logical, intent(in), optional :: revert
 integer, intent(in), optional :: mytype

 integer :: i
 real    :: r(3), radius_squared
 logical :: use_revert

 if (present(revert)) then
    use_revert = revert
 else
    use_revert = .false.
 endif

 radius_squared = radius**2

 if (present(mytype)) then
    do i=1,np
       r = xyzh(1:3,i) - center
       if (use_revert) then
          if (dot_product(r,r) < radius_squared .and. iamtype(iphase(i)) == mytype) call kill_particle(i,npartoftype)
       else
          if (dot_product(r,r) > radius_squared .and. iamtype(iphase(i)) == mytype) call kill_particle(i,npartoftype)
       endif
    enddo
 else
    do i=1,np
       r = xyzh(1:3,i) - center
       if (use_revert) then
          if (dot_product(r,r) < radius_squared) call kill_particle(i,npartoftype)
       else
          if (dot_product(r,r) > radius_squared) call kill_particle(i,npartoftype)
       endif
    enddo
 endif
 call shuffle_part(np)
 if (np /= sum(npartoftype)) call fatal('del_part_outside_sphere','particles not conserved')

end subroutine delete_particles_outside_sphere

!----------------------------------------------------------------
!+
!  Delete particles outside of a defined cylinder
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_cylinder(center, radius, zmax)
 real, intent(in) :: center(3), radius, zmax
 integer :: i
 real :: x, y, z, rcyl

 do i=1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    rcyl=sqrt((x-center(1))**2 + (y-center(2))**2)
    if (rcyl > radius .or. abs(z) > zmax) call kill_particle(i,npartoftype)
 enddo

end subroutine delete_particles_outside_cylinder

!----------------------------------------------------------------
!+
!  Delete particles within radius
!+
!----------------------------------------------------------------
subroutine delete_dead_particles_inside_radius(center,radius,np)
 use io, only:fatal
 real, intent(in) :: center(3), radius
 integer, intent(inout) :: np
 integer :: i
 real :: r(3), radius_squared

 radius_squared = radius**2
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       r = xyzh(1:3,i) - center
       if (dot_product(r,r) > radius_squared) call kill_particle(i,npartoftype)
    endif
 enddo
 call shuffle_part(np)
 if (np /= sum(npartoftype)) call fatal('del_dead_part_outside_sphere','particles not conserved')

end subroutine delete_dead_particles_inside_radius

!----------------------------------------------------------------
!+
!  Delete particles within radius
!+
!----------------------------------------------------------------
subroutine delete_particles_inside_radius(center,radius,npart,npoftype)
 real, intent(in) :: center(3), radius
 integer, intent(inout) :: npart,npoftype(:)
 integer :: i
 real :: x,y,z,r

 do i=1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    r=sqrt((x-center(1))**2+(y-center(2))**2+(z-center(3))**2)
    if (r < radius) call kill_particle(i,npoftype)
 enddo
 call shuffle_part(npart)

end subroutine delete_particles_inside_radius

!----------------------------------------------------------------
 !+
 ! Accrete particles outside a given radius
 !+
 !----------------------------------------------------------------
subroutine accrete_particles_outside_sphere(radius)
 real, intent(in) :: radius
 integer :: i
 real :: r2
 !
 ! accrete particles outside some outer radius
 !
 !$omp parallel default(none) &
 !$omp shared(npart,nptmass,xyzh,xyzmh_ptmass,radius) &
 !$omp private(i,r2)
 !$omp do
 do i=1,npart
    r2 = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
    if (r2 > radius**2) xyzh(4,i) = -abs(xyzh(4,i))
 enddo
 !$omp enddo

 !$omp do
 do i=1,nptmass
    r2 = xyzmh_ptmass(1,i)**2 + xyzmh_ptmass(2,i)**2 + xyzmh_ptmass(3,i)**2
    if (r2 > radius**2) xyzmh_ptmass(4,i) = -abs(xyzmh_ptmass(4,i))
 enddo
!$omp enddo

 !$omp end parallel

end subroutine accrete_particles_outside_sphere

!----------------------------------------------------------------
!+
!  Returns Keplerian frequency of particle i
!  USE WITH EXTREME CAUTION
!+
!----------------------------------------------------------------
real function Omega_k(i)
 integer, intent(in)  :: i
 real                 :: m_star,r
 integer              :: j

 m_star = 0.

 !-- if flyby: r should be centered around the primary, else around center of mass
 if (this_is_a_flyby) then
    r = sqrt((xyzh(1,i)-xyzmh_ptmass(1,1))**2 + (xyzh(2,i)-xyzmh_ptmass(2,1))**2 + (xyzh(3,i)-xyzmh_ptmass(3,1))**2)
 else
    r = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
 endif

 if (this_is_a_flyby) then
    m_star = xyzmh_ptmass(4,1)
 else
    do j=1,nptmass
       if (xyzmh_ptmass(4,j) > 0.) m_star = m_star + xyzmh_ptmass(4,j)
    enddo
 endif

 if (r > 0. .and. m_star > 0.) then
    Omega_k = sqrt(m_star/r) / r
 elseif (this_is_a_test) then
    Omega_k = 1/(r**1.5)
 else
    Omega_k = 0.
 endif

end function Omega_k

subroutine update_npartoftypetot
 use mpiutils, only:reduceall_mpi

 npartoftypetot = reduceall_mpi('+',npartoftype)

end subroutine update_npartoftypetot

end module part
