!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mcfost_utils
!
! This module contains options for the interface to the MCFOST
!  Monte Carlo Radiation Transfer code when called live from Phantom
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - use_mcfost : *use the mcfost library*
!
! :Dependencies: dim, infile_utils
!
 implicit none
 logical, public :: use_mcfost,use_Voronoi_limits_file,use_mcfost_stellar_parameters,mcfost_computes_Lacc
 logical, public :: mcfost_uses_PdV,mcfost_dust_subl
 integer, public :: ISM
 integer, parameter :: sp = 4
 real(kind=sp), public :: mcfost_keep_part
 character(len=80), public :: Voronoi_limits_file

 public :: set_defaults_mcfost,write_options_mcfost,read_options_mcfost

 private

contains

!-----------------------------------------------------------------
!+
!  set defaults for mcfost options
!+
!-----------------------------------------------------------------
subroutine set_defaults_mcfost()

 use_mcfost = .false.
 use_mcfost_stellar_parameters = .false.
 mcfost_computes_Lacc = .false.
 mcfost_dust_subl = .false.
 mcfost_uses_PdV = .true.
 mcfost_keep_part = 0.999_sp
 ISM = 0
 use_Voronoi_limits_file = .false.

end subroutine set_defaults_mcfost

!-----------------------------------------------------------------
!+
!  write options for mcfost
!+
!-----------------------------------------------------------------
subroutine write_options_mcfost(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(use_mcfost,'use_mcfost','use the mcfost library',iunit)
 if (use_Voronoi_limits_file) call write_inopt(Voronoi_limits_file,'Voronoi_limits_file',&
      'Limit file for the Voronoi tesselation',iunit)
 call write_inopt(use_mcfost_stellar_parameters,'use_mcfost_stars',&
      'Fix the stellar parameters to mcfost values or update using sink mass',iunit)
 call write_inopt(mcfost_computes_Lacc,'mcfost_computes_Lacc',&
      'Should mcfost compute the accretion luminosity',iunit)
 call write_inopt(mcfost_uses_PdV,'mcfost_uses_PdV',&
      'Should mcfost use the PdV work and shock heating?',iunit)
 call write_inopt(mcfost_keep_part,'mcfost_keep_part',&
      'Fraction of particles to keep for MCFOST',iunit)
 call write_inopt(ISM,'ISM',&
      'ISM heating : 0 -> no ISM radiation field, 1 -> ProDiMo, 2 -> Bate & Keto',iunit)
 call write_inopt(mcfost_dust_subl,'mcfost_dust_subl',&
      'Should mcfost do dust sublimation (experimental!)',iunit)

end subroutine write_options_mcfost

!-----------------------------------------------------------------
!+
!  read options for mcfost
!+
!-----------------------------------------------------------------
subroutine read_options_mcfost(db,nerr)
 use dim,          only:track_lum
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 real :: tmp

 call read_inopt(use_mcfost,'use_mcfost',db,errcount=nerr) ! compulsory
 if (use_Voronoi_limits_file) call read_inopt(Voronoi_limits_file,'Voronoi_limits_file',db,errcount=nerr)
 call read_inopt(use_mcfost_stellar_parameters,'use_mcfost_stars',db,errcount=nerr,default=use_mcfost_stellar_parameters)
 call read_inopt(mcfost_computes_Lacc,'mcfost_computes_Lacc',db,errcount=nerr,default=mcfost_computes_Lacc)
 call read_inopt(mcfost_uses_PdV,'mcfost_uses_PdV',db,errcount=nerr) ! compulsory

 tmp = real(mcfost_keep_part) ! single precision but read as real
 call read_inopt(tmp,'mcfost_keep_part',db,errcount=nerr,default=tmp)
 mcfost_keep_part = real(tmp,kind=sp)

 call read_inopt(ISM,'ISM',db,errcount=nerr,default=ISM)
 call read_inopt(mcfost_dust_subl,'mcfost_dust_subl',db,errcount=nerr,default=mcfost_dust_subl)

 if (mcfost_uses_PdV) track_lum = .true.

end subroutine read_options_mcfost

end module mcfost_utils
