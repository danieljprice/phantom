!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Change the units of a dumpfile.
!
! :References: None
!
! :Owner: josh
!
! :Runtime parameters: None
!
! :Dependencies: prompting, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use prompting,    only:prompt
 use units, only:umass,udist,utime
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: utime_bool,udist_bool,umass_bool,utime_fixed,udist_fixed,umass_fixed,fixed_tot,total_units_to_change
 real    :: umass_factor,utime_factor,udist_factor,umass_tmp,utime_tmp,udist_tmp,grav_const

 grav_const = udist**3/(utime**2*umass)

 utime_bool = 0
 udist_bool = 0
 umass_bool = 0

 utime_fixed = 0
 udist_fixed = 0
 umass_fixed = 0

 utime_factor = 1.0
 udist_factor = 1.0
 umass_factor = 1.0

 print*,'Current time unit is ',utime,'.'
 call prompt('Would you like time unit to be adjusted',utime_bool,0,1)
 if (utime_bool==1) then
    call prompt('Would you like time unit to be adjusted by a fixed value',utime_fixed,0,1)
    if (utime_fixed==1) then
       call prompt('Enter in value you want to scale time unit by',utime_factor,0.)
    endif
 endif

 print*,'Current length unit is ',udist,'.'
 call prompt('Would you like length unit to be adjusted',udist_bool,0,1)
 if (udist_bool==1) then
    call prompt('Would you like length unit to be adjusted by a fixed value',udist_fixed,0,1)
    if (udist_fixed==1) then
       call prompt('Enter in value you want to scale lenth unit by',udist_factor,0.)
    endif
 endif

 fixed_tot = udist_fixed + utime_fixed

 print*,'Current length unit is ',umass,'.'
 call prompt('Would you like mass unit to be adjusted', umass_bool,0,1)
 if (umass_bool==1) then
    if (fixed_tot<2) then
       call prompt('Would you like mass unit to be adjusted by a fixed value',udist_fixed,0,1)
    else
       print*,'Cannot change the mass unit by a fixed value since two other units are being adjusted by fixed value.'
    endif

    if (udist_fixed==1) then
       call prompt('Enter in value you want to scale lenth unit by',udist_factor,0.)
    endif
 endif

 ! Check to make sure that more than one unit is being changed
 total_units_to_change = utime_bool+udist_bool+umass_bool

 if (total_units_to_change<=1) then
    print*,'You must allow more than one unit to be changed'
    stop
 endif

 ! Begin writing tmp code units
 utime_tmp = utime_factor*utime
 udist_tmp = udist_factor*udist
 umass_tmp = umass_factor*umass

 print*,'Temporary time unit is ',utime_tmp,'.'
 print*,'Temporary length unit is ',udist_tmp,'.'
 print*,'Temporary mass unit is ',umass_tmp,'.'

 ! Begin changing code units that are allowed to vary, but are not fixed
 if ((utime_fixed==0) .AND. (utime_bool==1)) then
    utime_tmp = (udist_tmp**3/(grav_const*umass_tmp))**(1.0/2.0)
 endif

 if ((udist_fixed==0) .AND. (udist_bool==1)) then
    udist_tmp = (grav_const*utime_tmp**2*umass_tmp)**(1.0/3.0)
 endif

 if ((umass_fixed==0) .AND. (umass_bool==1)) then
    umass_tmp = udist_tmp**3/(grav_const*utime_tmp**2)
 endif

 print*,'Temporary time unit is ',utime_tmp,'.'
 print*,'Temporary length unit is ',udist_tmp,'.'
 print*,'Temporary mass unit is ',umass_tmp,'.'

 ! Check that everything has been changed properly
 print*,'Gravitational constant is ',grav_const,'.'
 print*,'With new units, it is now ', udist_tmp**3/(utime_tmp**2*umass_tmp),'.'

 ! Write new code units to header
 umass = umass_tmp
 udist = udist_tmp
 utime = utime_tmp

end subroutine modify_dump

end module moddump

