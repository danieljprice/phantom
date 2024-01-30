!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module sphNGutils
!
! This module contains routines allowing read/write compatibility
!  with sphNG dump files
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, part
!
 implicit none
 ! labels for sphNG types, used when converting dumps (these cannot duplicate current itypes)
 integer, parameter :: isphNG_accreted  = 18
 integer, parameter :: isphNG_sink_temp = 19
 real,allocatable,public :: mass_sphng(:),spin_sphng(:,:)
 logical,public :: got_mass=.false.,got_spin(3)=.false.

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
 use part, only:iamtype,ihacc,ihsoft,set_particle_type,igas,kill_particle,ispinx,ispiny,ispinz,&
       npartoftype,iunknown,isdead,shuffle_part
 integer :: i,nsink
 integer, intent(inout)    :: npart
 integer, intent(inout) :: nptmass
 integer(kind=1), intent(inout) :: iphase(:)
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer,         intent(out)   :: ierr
 integer :: gascount,unkncount,othercount

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
          xyzmh_ptmass(4,nsink) = mass_sphng(i)
          xyzmh_ptmass(ihacc,nsink) = xyzh(4,i)
          xyzmh_ptmass(ihsoft,nsink) = xyzh(4,i)
          vxyz_ptmass(1,nsink) = vxyzu(1,i)
          vxyz_ptmass(2,nsink) = vxyzu(2,i)
          vxyz_ptmass(3,nsink) = vxyzu(3,i)
          if (nptmass < 100) then
             print "(1x,a,i2,a,es13.6,a,es10.3,a,es13.6,a)",'[CONVERTING SINK #',nsink,' sphNG->Phantom, M=',&
                   xyzmh_ptmass(4,nsink),' h= ',xyzmh_ptmass(ihacc,nsink),' spinx= ',xyzmh_ptmass(ispinx,nsink),']'
          endif
       endif
       !
       !--remove the particle from the SPH data arrays
       !
       call set_particle_type(i,igas) ! Phantom will not run if any particles have iphase=iunknown
       call kill_particle(i)          ! particle is killed, so iphase is irrelevant
       print *, "npartoftype 1",npartoftype(1)
       print *, "killed",i,iamtype(iphase(i)),xyzh(4,i)
       print *, "npartoftype 1",npartoftype(1)
    endif
 enddo
 if (nsink < nptmass) then
    write(*,*) 'ERROR converting sinks from sphNG->Phantom, nsinks from iphase < nptmass'
    write(*,*) 'using nptmass = ',nptmass
    nptmass = nsink
    ierr = 67
 endif

 gascount = 0
 unkncount = 0
 othercount = 0
 do i=1,npart
    if (iamtype(iphase(i)) == igas) then
       gascount = gascount+1
    elseif (iamtype(iphase(i)) == iunknown) then
       unkncount = unkncount + 1
    else
       othercount = othercount + 1
    endif
    if (isdead(i)) then
       print *, i," is dead"
    endif
 enddo

 call shuffle_part(npart) !to remove killed particles
 gascount = 0
 unkncount = 0
 othercount = 0
 do i=1,npart
    if (iamtype(iphase(i)) == igas) then
       gascount = gascount+1
    elseif (iamtype(iphase(i)) == iunknown) then
       unkncount = unkncount + 1
    else
       othercount = othercount + 1
    endif
    if (isdead(i)) then
       print *, i," is dead"
    endif
 enddo
 print *, "gas:",gascount,"unkncount:",unkncount,"othercount:",&
      othercount, "total=", gascount+unkncount+othercount
end subroutine convert_sinks_sphNG

subroutine set_gas_particle_mass(mass_sphng)
 use part, only: massoftype,igas,iphase,iamtype,hfact
 use dim, only: maxp
 real,intent(in) :: mass_sphng(maxp)
 integer :: i

 if (hfact < 1) then
    hfact = 1.2
    print *, "Set hfact=1.2"
 endif
 print *, "Setting gas particle mass"
 do i=1,maxp
    if (iamtype(iphase(i)) == igas) then
       massoftype(igas) = mass_sphng(i)
       print *, "Mass of gas particles =", massoftype(igas)
       exit
    endif
 enddo
 if (i == maxp) then
    print *, "gas part mass not set!"
 else
    print *, "in set_gas...: i= ",i
 endif
end subroutine set_gas_particle_mass

end module sphNGutils
