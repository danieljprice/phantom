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
!  DEPENDENCIES: dim, part, readwrite_dust
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:use_dust,ndusttypes
 use part,           only:igas,idust,set_particle_type
 use readwrite_dust, only:write_temp_grains_file,set_dustfrac_from_inopts
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,dust_method
 real    :: dust_to_gas
 real    :: dustfrac_percent(ndusttypes) = 0.

 if (use_dust) then
    call write_temp_grains_file(dust_to_gas,dustfrac_percent,imethod=dust_method)
    if (dust_method == 1) then
       call set_dustfrac_from_inopts(dust_to_gas,percent=dustfrac_percent)

       massoftype(igas) = massoftype(igas)*(1. + dust_to_gas)
    elseif (dust_method == 2) then
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

