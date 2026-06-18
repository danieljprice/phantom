!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Convert a dust-nucleation file into a dust-growth file. Initial conditions of the dust-growth simulation are provided by the dust-nucleation simulation.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
use dim,          only:use_dust,maxdusttypes,maxdustlarge,maxdustsmall,use_dustgrowth,&
                       update_max_sizes,isothermal,do_nucleation
use part,         only:igas,idust,ndusttypes,ndustsmall,&
                       grainsize,graindens,dustfrac,&
                       nucleation
use dust_formation, only:mass_per_H
use set_dust,     only:set_dustfrac,set_dustbinfrac
use options,      only:use_dustfrac
use growth,       only:set_dustprop,convert_to_twofluid,iporosity,init_growth,ifrag,gsizemincgs,&
                       ivrelkin
use prompting,    only:prompt
use dust,         only:grainsizecgs,graindenscgs
use table_utils,  only:logspace
use units,   only:umass,udist
use physcon,  only:rho_Cdust,a0,mc
use cooling_solver, only:excitation_HI
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,dust_method,np_ratio,np_gas,maxdust
 real    :: dust_to_gas,smincgs,smaxcgs,sindex,udens
 real    :: pwl_sizedistrib,R_ref,H_R_ref,q_index
 logical :: sizedistrib


! Add flag for dust nucleation
 if (.not. (do_nucleation .and. use_dustgrowth .and. (.not. isothermal))) then
    print*,' DOING NOTHING: COMPILE WITH DUSTGROWTH=yes DUST_NUCLEATION=yes ISOTHERMAL=no'
    stop
 endif

! units
 udens = umass/(udist**3)

! parameters
 ndusttypes = 1              ! growth can only work with 1 dust population

! cooling
 excitation_HI = 1

! dust
 dust_method = 1             ! 1:dust as a mxiture ; 2:dust as particles (not working yet)
 np_ratio = 5                ! np_dust/np_gas if dust_method=2
 dust_to_gas = 0.01          ! only to inialise the variable, changed later

 sizedistrib = .false.       ! initial grain size distribution
 smincgs = 1.e-5
 smaxcgs = 1.
 sindex = 3.5

 pwl_sizedistrib = -2        ! spatial distribution of dust, if sizedistrib = .true.
 R_ref = 100
 H_R_ref = 0.0895
 q_index = 0.25

 gsizemincgs = 1.e-8         ! minimum allowed grainsize in cm
 ivrelkin = 0                ! keep 0 until relative dust motions from 1-fluid are available

 grainsizecgs = 1.           !dummies but need definition to be used by dust.f90
 graindenscgs = 1.


! Interactive prompts
if (dust_method==1) then
    maxdust = maxdustsmall
elseif (dust_method==2) then
    print*,' DOING NOTHING: NOT YET READY WITH DUST AS PARTICLES'
    stop
      !--We do not care if modulo(npart,np_ratio) is stricly zero, since npart can
      !  be a weird value depdending on the simulation it comes from.
   endif

ndusttypes = 1
if (use_dustgrowth) then
    ifrag = 1
    iporosity = 0 !turn on use_porosity if you change that
    grainsize(1) = grainsizecgs/udist
    graindens(1) = graindenscgs/udens
endif

np_gas = npartoftype(igas)

if (dust_method == 1) then

  use_dustfrac = .true.
  ndustsmall = ndusttypes

  do i=1,np_gas
     if (ndusttypes > 1) then
        print*,' DOING NOTHING: GROWTH NOT WORKING FOR MORE THAN 1 DUST POPULATION'
     else
        if (nucleation(5,i)<tiny(nucleation(5,i))) then
            dust_to_gas = tiny(dust_to_gas)
        else
            dust_to_gas = nucleation(5,i) * mc / mass_per_H
        endif
        call set_dustfrac(dust_to_gas,dustfrac(:,i))
     endif
  enddo

  massoftype(igas) = massoftype(igas)*(1. + dust_to_gas)
  npart = np_gas

   if (use_dustgrowth) then
      call set_dustprop(npart,xyzh,sizedistrib,pwl_sizedistrib,R_ref,H_R_ref,q_index)
   endif
endif





end subroutine modify_dump

end module moddump

