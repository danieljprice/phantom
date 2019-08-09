!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: extern_staticsine
!
!  DESCRIPTION:
! This module contains routines relating to the computation
! of a static sinusoid potential (in 1D), i.e.
! phi = A cos(k(x+B))
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    amplitude   -- Amplitude of sine perturbation
!    inclination -- Orientation angle of perturbation (rad, 0.0=aligned with x axis)
!    phase       -- Phase of perturbation
!    wavek       -- Wavenumber of perturbation
!
!  DEPENDENCIES: infile_utils, io, physcon
!+
!--------------------------------------------------------------------------
module extern_staticsine
 use physcon,   only: pi,twopi
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !

 real(kind=8), public :: amplitude = 100.0
 real(kind=8), public :: wavek = pi/4.0
 real(kind=8), public :: phase = 2.0
 real(kind=8), public :: inclination = 0.0

 public :: staticsine_force
 public :: write_options_staticsine, read_options_staticsine
 private

contains

!----------------------------------------------
!+
!  compute the force on a given particle from
!  the static sinusoid
!+
!----------------------------------------------
subroutine staticsine_force(xi,yi,fxi,fyi,fzi,phi)
 real, intent(in)  :: xi,yi
 real, intent(out) :: fxi,fyi,fzi,phi
 real :: xforce, fi

 !--Rotate in x-y plane to align to wavefront
 xforce = xi*cos(inclination) - yi*sin(inclination)

 !--Potential
 phi = amplitude*cos(wavek*(xforce+phase))

 !--Magnitude of Force (along axis aligned with sinusoid)
 fi = wavek*amplitude*sin(wavek*(xforce+phase))

 !--Separate into x and y components according to the inclination of the sinusoidal perturbation
 fxi = fi*cos(inclination)
 fyi = -fi*sin(inclination)
 fzi = 0.0

 return
end subroutine staticsine_force

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_staticsine(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(amplitude,'amplitude','Amplitude of sine perturbation',iunit)
 call write_inopt(wavek,'wavek','Wavenumber of perturbation',iunit)
 call write_inopt(phase,'phase','Phase of perturbation',iunit)
 call write_inopt(inclination,'inclination','Orientation angle of perturbation (rad, 0.0=aligned with x axis)',iunit)

end subroutine write_options_staticsine

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_staticsine(name,valstring,imatch,igotall,ierr)
 use physcon, only: pi
 use io,      only:fatal,error
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: where = 'read_options_externstaticsine'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('amplitude')
    read(valstring,*,iostat=ierr) amplitude
    ngot = ngot+1
 case('wavek')
    read(valstring,*,iostat=ierr) wavek
    ngot = ngot + 1
 case('phase')
    read(valstring,*,iostat=ierr) phase
    ngot = ngot+1
 case('inclination')
    read(valstring,*,iostat=ierr) inclination
    ngot = ngot + 1

    ! Prevent overly large values of inclination
    if (inclination > twopi) then
       inclination = mod(inclination,twopi)
       print"(a,es10.3)", 'Inclination> 2pi, set to mod(inclination,2pi)= ',inclination
    endif
 end select

 igotall = (ngot >= 4)

end subroutine read_options_staticsine

end module extern_staticsine
