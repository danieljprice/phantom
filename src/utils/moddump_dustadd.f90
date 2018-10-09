!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  adds dust particles to pre-existing gas particle dump
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, dust, growth, options, part, prompting, set_dust,
!    table_utils
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,          only:use_dust,maxdusttypes,maxdustlarge,maxdustsmall,use_dustgrowth
 use part,         only:igas,idust,set_particle_type,ndusttypes,ndustsmall,ndustlarge,&
                        grainsize,graindens,dustfrac
 use set_dust,     only:set_dustfrac,set_dustbinfrac
 use options,      only:use_dustfrac
 use growth,       only:set_dustprop
 use prompting,    only:prompt
 use dust,         only:grainsizecgs,graindenscgs
 use table_utils,  only:logspace
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,j,itype,ipart,iloc,dust_method,np_ratio,np_gas,np_dust,maxdust
 real    :: dust_to_gas,smincgs,smaxcgs,sindex,dustbinfrac(maxdusttypes)

 if (.not. use_dust) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes'
    stop
 endif

 dust_method = 1
 np_ratio = 5
 dust_to_gas = 0.01
 ndusttypes = 1
 smincgs = 1.e-5
 smaxcgs = 1.
 sindex = 3.5
 dustbinfrac = 0.
 grainsize = 1.
 graindens = 3.
 grainsizecgs = 1.
 graindenscgs = 3.

 call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)

 if (dust_method==1) then
    maxdust = maxdustsmall
 elseif (dust_method==2) then
    maxdust = maxdustlarge
    call prompt('Enter ratio between number of gas particles and dust particles',np_ratio,1,10)
    !--We do not care if modulo(npart,np_ratio) is stricly zero, since npart can
    !  be a weird value depdending on the simulation it comes from.
 endif

 call prompt('Enter total dust to gas ratio',dust_to_gas,0.)
 call prompt('How many grain sizes do you want?',ndusttypes,1,maxdust)

 if (ndusttypes > 1) then
    !--grainsizes
    call prompt('Enter minimum grain size in cm',smincgs,0.)
    call prompt('Enter maximum grain size in cm',smaxcgs,0.)
    call logspace(grainsize(1:ndusttypes),smincgs,smaxcgs)
    !--mass distribution
    call prompt('Enter power-law index, e.g. MRN',sindex)
    call set_dustbinfrac(smincgs,smaxcgs,sindex,dustbinfrac(1:ndusttypes))
    !--grain density
    call prompt('Enter grain density in g/cm^3',graindens(1),0.)
    graindens = graindens(1)
 else
    if (use_dustgrowth) then
       call prompt('Enter initial grain size in cm',grainsizecgs,0.)
    else
       call prompt('Enter grain size in cm',grainsizecgs,0.)
    endif
    call prompt('Enter grain density in g/cm^3',graindenscgs,0.)
    grainsize(1) = grainsizecgs
    graindens(1) = graindenscgs
 endif

 np_gas = npartoftype(igas)

 if (dust_method == 1) then

    use_dustfrac = .true.
    ndustsmall = ndusttypes

    do i=1,np_gas
       if (ndusttypes > 1) then
          dustfrac(1:ndusttypes,i) = dust_to_gas*dustbinfrac(1:ndusttypes)
       else
          call set_dustfrac(dust_to_gas,dustfrac(:,i))
       endif
    enddo

    massoftype(igas) = massoftype(igas)*(1. + dust_to_gas)
    npart = np_gas

 elseif (dust_method == 2) then

    use_dustfrac = .false.
    ndustlarge = ndusttypes
    np_dust = np_gas/np_ratio

    do i=1,ndustlarge

       itype = idust + i - 1
       npartoftype(itype) = np_dust
       massoftype(itype)  = massoftype(igas)*dust_to_gas*np_ratio

       do j=1,np_dust
          ipart = np_gas + (i-1)*np_dust + j
          iloc  = np_ratio*j

          xyzh(1,ipart) = xyzh(1,iloc)
          xyzh(2,ipart) = xyzh(2,iloc)
          xyzh(3,ipart) = xyzh(3,iloc)
          xyzh(4,ipart) = xyzh(4,iloc)

          vxyzu(1,ipart) = vxyzu(1,iloc)
          vxyzu(2,ipart) = vxyzu(2,iloc)
          vxyzu(3,ipart) = vxyzu(3,iloc)

          call set_particle_type(ipart,itype)
       enddo

    enddo

    npart = np_gas + np_dust*ndustlarge

    if (use_dustgrowth) then
       call set_dustprop(npart)
    endif

 endif

end subroutine modify_dump

end module moddump

