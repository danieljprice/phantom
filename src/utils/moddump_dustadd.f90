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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, dust, options, part
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,     only:use_dust,ndusttypes
 use part,    only:dustfrac,igas,idust,set_particle_type
 use dust,    only:set_dustfrac,smincgs,smaxcgs,sindex
 use options, only:use_dustfrac
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 real    :: dust_to_gas,dustfrac_temp(ndusttypes)

 if (use_dust) then
    dust_to_gas = 0.01
    print*,' SETTING DUST-TO-GAS RATIO = ',dust_to_gas
    if (use_dustfrac) then
       if (ndusttypes>1) then
          call set_dustfrac(dust_to_gas,dustfrac_temp,smincgs,smaxcgs,sindex)
       else
          call set_dustfrac(dust_to_gas,dustfrac_temp)
       endif
       do i=1,ndusttypes
          dustfrac(i,:) = dustfrac_temp(i)
       enddo
       massoftype(igas) = massoftype(igas)*(1. + dust_to_gas)
    else
       npart = npartoftype(igas)
       npartoftype(idust) = npart
       massoftype(idust)  = massoftype(igas)*dust_to_gas

       do i=npart+1,2*npart
          xyzh(1,i) = xyzh(1,i-npart)
          xyzh(2,i) = xyzh(2,i-npart)
          xyzh(3,i) = xyzh(3,i-npart)
          xyzh(4,i) = xyzh(4,i-npart)

          vxyzu(1,i) = vxyzu(1,i-npart)
          vxyzu(2,i) = vxyzu(2,i-npart)
          vxyzu(3,i) = vxyzu(3,i-npart)
          call set_particle_type(i,idust)
       enddo
       npart=2*npart
    endif
 else
    print*,' DOING NOTHING: COMPILE WITH DUST=yes'
 endif

end subroutine modify_dump

end module moddump

