!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: params
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------

!*******************************************************************************
module params
 integer, parameter :: mfile=60,mlabel=60                                        ! string lengths
 integer, parameter :: txt_unit=10,in_unit=11,out_unit=12,sink_unit=13, &        ! i/o unit numbers
                      power_unit=14
 logical, save :: do_seq_read=.false.,do_seq_write=.false.,do_exp=.true., &      ! flags
                 do_massflux=.false.
 real offset                                                                   ! vector component offsets
 integer x_index(2)                                                            ! vector x-component (U & B)
 character(len=70), save :: hl="----------------------------------------------------------------------"
END module
