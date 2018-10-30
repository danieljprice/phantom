!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: checkoptions
!
!  DESCRIPTION:
!  this module performs checks of the compile time options
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, part
!+
!--------------------------------------------------------------------------
module checkoptions
 implicit none
 public :: check_compile_time_settings

 private

contains

!-------------------------------------------------------------------
!
! This subroutine checks that the compile time options are sensible
! and mutually compatible with each other
!
!-------------------------------------------------------------------
subroutine check_compile_time_settings(ierr)
 use part,  only:mhd,maxBevol,gravity,ngradh,h2chemistry,maxvxyzu,use_dust,gr
 use dim,   only:use_dustgrowth,maxtypes,maxsts
 use io,    only:error,id,master,fatal,warning
#ifdef GR
 use metric_tools, only:icoordinate,icoord_cartesian
#endif
 integer, intent(out) :: ierr
 character(len=16), parameter :: string = 'compile settings'

 ierr = 0
!
!--check MHD dimension settings are OK
!
#ifdef MHD
 if (.not.mhd) then
    if (id==master) call error(string,'-DMHD but mhd=.false.')
    ierr = 1
 endif
#endif
 if (mhd) then
    if (maxBevol < 3 .or. maxBevol > 4) then
       if (id==master) call error(string,'must have maxBevol=3 (no cleaning) or maxBevol=4 (cleaning)')
       ierr = 1
    endif
 endif
#ifdef NONIDEALMHD
 if (.not.mhd) then
    if (id==master) call error(string,'-DNONIDEALMHD requires -DMHD')
    ierr = 1
 endif
#ifdef USE_CMAC_IONISE
 if (id==master) call error(string,'can not use both -DNONIDEALMHD and -DUSE_CMAC_IONISE')
 ierr = 1
#endif
#endif

#ifdef GRAVITY
 if (.not.gravity) call error(string,'-DGRAVITY but gravity=.false.')
#endif
 if (gravity) then
    if (ngradh < 2) then
       if (id==master) call error(string,'gravity requires ngradh=2 for gradsoft storage')
       ierr = 2
    endif
 endif
 if (h2chemistry) then
    if (maxvxyzu /= 4) then
       if (id==master) call error(string,'must store thermal energy (maxvxyzu=4 in dim file) if using H2 chemistry')
       ierr = 3
    endif
 endif
 if (maxtypes > 64) then
    if (id==master) call error(string,'cannot use more than 64 particle types' // &
       ' unless iphase is changed to int*2')
    ierr = 4
 endif


#ifdef DISC_VISCOSITY
#ifdef CONST_AV
 if (id==master) call error(string,'should not use both -DCONST_AV and -DDISC_VISCOSITY')
 ierr = 4
#endif
#endif

#ifdef DUST
#ifdef MHD
 if (id==master) call error(string,'-DDUST currently not compatible with magnetic fields (-DMHD)')
#endif
#endif

#ifdef GR
 if (mhd) then
    call error(string,'General relativity not compatible with MHD.')
    ierr = 6
 endif
 if (use_dust) then
    call error(string,'General relativity not compatible with dust.')
    ierr = 7
 endif
 if (gravity) then
    call warning(string,'You are using SELF GRAVITY in GENERAL RELATIVITY. Proceed with caution...!')
 endif
 if (h2chemistry) then
    call error(string,'General relativity not compatible with chemistry.')
    ierr = 8
 endif
 if (maxsts > 1) then
    call error(string,'General relativity not compatible with super-timestepping.')
    ierr = 10
 endif
#ifdef DRIVING
 call error(string,'General relativity not compatible with turbulent driving.')
 ierr = 11
#endif
 if (icoordinate /= icoord_cartesian) then
    call fatal('checkoptions (GR)',&
   "You must use Cartesian-like coordinates in PHANTOM! Please change to Cartesian in metric_tools!'")
    ierr = 12
 endif
 if (.not.gr) then
    call error(string,'-DGR but gr=.false.')
    ierr = 13
 endif
#ifdef DISC_VISCOSITY
 call error(string,'General relativity not compatible with disc viscosity.')
 ierr = 14
#endif
#ifndef CONST_AV
 call error(string,'General relativity should have CONST_AV=yes.')
 ierr = 15
#endif
#endif

#ifdef DUSTGROWTH
 if (.not. use_dustgrowth) call error(string,'-DDUSTGROWTH but use_dustgrowth = .false.')
#endif

 return
end subroutine check_compile_time_settings

end module checkoptions
