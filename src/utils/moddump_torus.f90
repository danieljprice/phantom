!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! add velocity and density perturbation to a torus after it has been relaxed
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Rtorus : *torus radius [code units]*
!   - ampl   : *amplitude of the m-mode density perturbation*
!   - mmode  : *azimuthal mode number m*
!
! :Dependencies: moddump_utils, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters
 real :: ampl   = 0.05   ! amplitude of the m-mode density perturbation
 real :: mmode  = 3.0    ! azimuthal mode number m
 real :: Rtorus = 1.0    ! torus radius [code units]

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ii
 real    :: phi,x,y,r
 real    :: rcyl2,rcyl,rsph,v2onr,omegai

 ! Implementing the m-mode density perturbation from Price & Bate 2007
 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 ! Adding velocities
 do ii = 1,npart
    rcyl2 = dot_product(xyzh(1:2,ii),xyzh(1:2,ii))
    rcyl  = sqrt(rcyl2)
    rsph  = sqrt(rcyl2 + xyzh(3,ii)*xyzh(3,ii))
    v2onr = 1./(Rtorus)*(-Rtorus*rcyl/rsph**3 + Rtorus**2/(rcyl2*rcyl)) + rcyl/rsph**3

    omegai = sqrt(v2onr/rcyl)
    vxyzu(1,ii) = -omegai*xyzh(2,ii)
    vxyzu(2,ii) = omegai*xyzh(1,ii)
    vxyzu(3,ii) = 0.
 enddo
 print*,'Velocities added.'

 do ii=1,npart
    x=xyzh(1,ii)
    y=xyzh(2,ii)
    r = sqrt(x**2 + y**2)
    phi = atan2(y,x)
    phi = phi - 0.5*ampl*sin(mmode*phi)
    xyzh(1,ii) = r*cos(phi)
    xyzh(2,ii) = r*sin(phi)
 enddo

 print*,'Density perturbation added.'

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter amplitude of the density perturbation',ampl,0.)
 call prompt('Enter azimuthal mode number m',mmode,0.)
 call prompt('Enter torus radius (code units)',Rtorus,0.)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(ampl,'ampl',db,errcount=nerr,min=0.)
 call read_inopt(mmode,'mmode',db,errcount=nerr,min=0.)
 call read_inopt(Rtorus,'Rtorus',db,errcount=nerr,min=0.)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# torus perturbation parameters'
 call write_inopt(ampl,'ampl','amplitude of the m-mode density perturbation',iunit)
 call write_inopt(mmode,'mmode','azimuthal mode number m',iunit)
 call write_inopt(Rtorus,'Rtorus','torus radius [code units]',iunit)

end subroutine write_options_moddump

end module moddump
