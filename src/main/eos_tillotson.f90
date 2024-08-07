!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_tillotson
!
! EoS from Tillotson et al. (1962)
!
! :References: https://apps.dtic.mil/sti/pdfs/AD0486711.pdf
!
! Implementation from Brundage, A. 2013
! :DOI: 10.1016/j.proeng.2013.05.053
!
! :Owner:
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, io
!
 implicit none

contains
!-----------------------------------------------------------------------
!+
!  EoS from Tillotson et al. (1962)
!+
!-----------------------------------------------------------------------
subroutine equationofstate_tillotson(rho,energy,temperature,mu,X,Z,pressure,spsound,gamma)
 real, intent(in) :: rho,energy,temperature,mu,X,Z
 real, intent(out) :: pressure,spsound,gamma

end subroutine equationofstate_tillotson

!----------------------------------------------------------------
!+
!  print eos information
!+
!----------------------------------------------------------------
subroutine eos_info_tillotson(iprint)
 integer, intent(in) :: iprint

 write(iprint,"(/,a)") ' Tillotson EoS'

end subroutine eos_info_tillotson

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
 subroutine read_options_eos_tillotson(name,valstring,imatch,igotall,ierr)
  use io, only:fatal
  character(len=*),  intent(in)  :: name,valstring
  logical,           intent(out) :: imatch,igotall
  integer,           intent(out) :: ierr
  integer,           save        :: ngot  = 0
  character(len=30), parameter   :: label = 'eos_tillotson'
 
  imatch  = .true.
  select case(trim(name))
!   case('irecomb')
!      read(valstring,*,iostat=ierr) irecomb
!      if ((irecomb < 0) .or. (irecomb > 3)) call fatal(label,'irecomb = 0,1,2,3')
!      ngot = ngot + 1
  case default
     imatch = .false.
  end select
 
  igotall = (ngot >= 1)
 
 end subroutine read_options_eos_tillotson

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_tillotson(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

!  call write_inopt(irecomb,'irecomb','recombination energy to include. 0=H2+H+He, 1=H+He, 2=He, 3=none',iunit)

end subroutine write_options_eos_tillotson

end module eos_tillotson