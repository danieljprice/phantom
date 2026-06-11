!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! default moddump routine: does not make any modifications
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
                       update_max_sizes,isothermal
use part,         only:igas,idust,ndusttypes,ndustsmall,ndustlarge,&
                       grainsize,graindens,dustfrac,&
                       nucleation
use dust_formation, only:mass_per_H
use set_dust,     only:set_dustfrac,set_dustbinfrac
use options,      only:use_dustfrac,use_porosity
use growth,       only:set_dustprop,convert_to_twofluid,iporosity,init_growth,ifrag,gsizemincgs
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
 integer :: i,ierr,j,itype,ipart,iloc,dust_method,np_ratio,np_gas,np_dust,maxdust
 real    :: dust_to_gas,smincgs,smaxcgs,sindex,dustbinfrac(maxdusttypes),udens
 real    :: inradius,outradius,pwl_sizedistrib,R_ref,H_R_ref,q_index
 logical :: sizedistrib


! cooling values
 excitation_HI = 1

! dust default values
 udens = umass/(udist**3)
 dust_method = 1
 np_ratio = 5
 dust_to_gas = 0.01
 ndusttypes = 1
 smincgs = 1.e-5
 smaxcgs = 1.
 sindex = 3.5
 dustbinfrac = 0.
 pwl_sizedistrib = -2
 R_ref = 100
 H_R_ref = 0.0895
 q_index = 0.25

 sizedistrib = .false.

 grainsizecgs = 1. !dummy
 graindenscgs = 1.

 gsizemincgs = 1.e-8


! Interactive prompts
!call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
dust_method = 1
if (dust_method==1) then
    maxdust = maxdustsmall
elseif (dust_method==2) then
    maxdust = maxdustlarge
    call prompt('Enter ratio between number of gas particles and dust particles',np_ratio,1)
      !--We do not care if modulo(npart,np_ratio) is stricly zero, since npart can
      !  be a weird value depdending on the simulation it comes from.
   endif

! call prompt('Enter total dust to gas ratio',dust_to_gas,0.) computed from nuc.
   
!call prompt('How many grain sizes do you want?',ndusttypes,1,maxdust)
ndusttypes = 1
!   if (ndusttypes > 1) then  !uncomment if several dust populations
!      !--grainsizes
!      call prompt('Enter minimum grain size in cm',smincgs,0.)
!      call prompt('Enter maximum grain size in cm',smaxcgs,0.)
!      !--mass distribution
!      call prompt('Enter power-law index, e.g. MRN',sindex)
!      call set_dustbinfrac(smincgs/udist,smaxcgs/udist,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
!      !--grain density
!      call prompt('Enter grain density in g/cm^3',graindens(1),0.)
!      graindens = graindens(1)/udens
!   else
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
        dustfrac(1:ndusttypes,i) = dust_to_gas*dustbinfrac(1:ndusttypes)
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

