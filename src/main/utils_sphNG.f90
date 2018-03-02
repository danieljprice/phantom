!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: sphNGutils
!
!  DESCRIPTION:
!  This module contains routines allowing read/write compatibility
!  with sphNG dump files
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part
!+
!--------------------------------------------------------------------------
module sphNGutils
 implicit none
 ! labels for sphNG types, used when converting dumps (these cannot duplicate current itypes)
 integer, parameter :: isphNG_accreted  = 18
 integer, parameter :: isphNG_sink_temp = 19


 public

contains
!--------------------------------------------------------------
!+
!  utility functions to convert phantom iphase to/from sphNG iphase
!+
!--------------------------------------------------------------
elemental function iphase_sphNG(iphasein)
 use part, only:iamtype,idust,iunknown
 integer(kind=1), intent(in) :: iphasein
 integer(kind=1) :: iphase_sphNG
 integer :: itype

 itype = iamtype(iphasein)
 select case(itype)
 case(idust)
    iphase_sphNG = 20
 case(iunknown)
    iphase_sphNG = -1
 case default
    iphase_sphNG = 0
 end select

end function iphase_sphNG

function itype_from_sphNG_iphase(iphasein)
 use part, only:idust,igas,iunknown
 integer(kind=1), intent(in) :: iphasein
 integer(kind=1) :: itype_from_sphNG_iphase

 select case(iphasein)
 case(20)
    itype_from_sphNG_iphase = idust
 case(1:10)
    itype_from_sphNG_iphase = isphNG_sink_temp
 case(0)
    itype_from_sphNG_iphase = igas
 case(:-1)
    itype_from_sphNG_iphase = isphNG_accreted
 case default
    itype_from_sphNG_iphase = iunknown
 end select

end function itype_from_sphNG_iphase

logical function is_sphNG_sink(iphasein)
 integer(kind=1), intent(in) :: iphasein

 if (iphasein > 1 .and. iphasein <= 10) then
    is_sphNG_sink = .true.
 else
    is_sphNG_sink = .false.
 endif

end function is_sphNG_sink

!--------------------------------------------------------------
!+
!  convert sink particle arrays/info to Phantom format
!+
!--------------------------------------------------------------
subroutine convert_sinks_sphNG(npart,nptmass,iphase,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,ierr)
 use part, only:iamtype,ihacc,ihsoft,set_particle_type,igas,kill_particle
 integer :: i,nsink
 integer, intent(in)    :: npart
 integer, intent(inout) :: nptmass
 integer(kind=1), intent(inout) :: iphase(:)
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer,         intent(out)   :: ierr

 nsink = 0
 ierr = 0
 do i=1,npart
    if (iamtype(iphase(i))==isphNG_sink_temp) then
       nsink = nsink + 1
       if (nsink > nptmass) then
          write(*,*) 'ERROR converting sinks from sphNG->Phantom, nsinks from iphase > nptmass'
          ierr = 66
       else
          xyzmh_ptmass(1,nsink) = xyzh(1,i)
          xyzmh_ptmass(2,nsink) = xyzh(2,i)
          xyzmh_ptmass(3,nsink) = xyzh(3,i)
          xyzmh_ptmass(ihacc,nsink) = xyzh(4,i)
          xyzmh_ptmass(ihsoft,nsink) = xyzh(4,i)
          vxyz_ptmass(1,nsink) = vxyzu(1,i)
          vxyz_ptmass(2,nsink) = vxyzu(2,i)
          vxyz_ptmass(3,nsink) = vxyzu(3,i)
          if (nptmass < 100) then
             print "(1x,a,i2,a,es13.6,a,es10.3,a)",'[CONVERTING SINK #',nsink,' sphNG->Phantom, M=',&
                   xyzmh_ptmass(4,nsink),' h= ',xyzmh_ptmass(ihacc,nsink),']'
          endif
       endif
       !
       !--remove the particle from the SPH data arrays
       !
       call set_particle_type(i,igas) ! Phantom will not run if any particles have iphase=iunknown
       call kill_particle(i)          ! particle is killed, so iphase is irrelevant
    endif
 enddo
 if (nsink < nptmass) then
    write(*,*) 'ERROR converting sinks from sphNG->Phantom, nsinks from iphase < nptmass'
    write(*,*) 'using nptmass = ',nptmass
    nptmass = nsink
    ierr = 67
 endif

end subroutine convert_sinks_sphNG

end module sphNGutils
