!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module checkoptions
!
! this module performs checks of the compile time options
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, metric_tools, mpiutils, part
!
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
 use part,     only:mhd,gravity,ngradh,maxvxyzu,use_dust,gr
 use dim,      only:use_dustgrowth,maxtypes,mpi,inject_parts,h2chemistry,driving,disc_viscosity
 use io,       only:error,id,master,warning
 use mpiutils, only:barrier_mpi
 use metric_tools, only:icoordinate,icoord_cartesian
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
#ifdef NONIDEALMHD
 if (.not.mhd) then
    if (id==master) call error(string,'-DNONIDEALMHD requires -DMHD')
    ierr = 1
 endif
#endif
!
!--check additional dimension settings are OK
!
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
!
!--check gravity flags are OK
!
#ifdef GRAVITY
 if (.not.gravity) call error(string,'-DGRAVITY but gravity=.false.')
#endif
 if (gravity) then
    if (ngradh < 2) then
       if (id==master) call error(string,'gravity requires ngradh=2 for gradsoft storage')
       ierr = 2
    endif
 endif

 if (use_dust .and. mhd .and. id==master) call error(string,'-DDUST currently not compatible with magnetic fields (-DMHD)')

 if (gr .and. mhd) then
    if (id==master) call error(string,'General relativity not compatible with MHD.')
    ierr = 6
 endif
 if (gr .and. use_dust) then
    if (id==master) call error(string,'General relativity not compatible with dust.')
    ierr = 7
 endif
 if (gr .and. gravity) then
    if (id==master) call warning(string,'You are using SELF GRAVITY in GENERAL RELATIVITY. Proceed with caution...!')
 endif
 if (gr .and. h2chemistry) then
    if (id==master) call error(string,'General relativity not compatible with chemistry.')
    ierr = 8
 endif
 if (gr .and. driving) then
    if (id==master) call error(string,'General relativity not compatible with turbulent driving.')
    ierr = 11
 endif
 if (gr .and. icoordinate /= icoord_cartesian) then
    if (id==master) call error('checkoptions (GR)','Must use Cartesian-like coordinates in PHANTOM! Please change in metric_tools!')
    ierr = 12
 endif
 if (gr .and. disc_viscosity) then
    if (id==master) call warning(string,'General relativity not properly tested with disc viscosity.')
 endif
#ifndef CONST_AV
 if (gr) then
    if (id==master) call error(string,'General relativity should have CONST_AV=yes.')
    ierr = 14
 endif
#endif

#ifdef DUSTGROWTH
 if (.not. use_dustgrowth) then
    if (id==master) call error(string,'-DDUSTGROWTH but use_dustgrowth = .false.')
    ierr = 15
 endif
#endif

 if (mpi .and. inject_parts) then
    if (id==master) call error(string,'MPI currently not compatible with particle injection')
    ierr = 16
 endif

 call barrier_mpi

end subroutine check_compile_time_settings

end module checkoptions
