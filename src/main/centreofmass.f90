!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: centreofmass
!
!  DESCRIPTION:
!   Utilities for computing the centre of mass on the particles
!   and correcting the bulk motion (used for turbulent driving)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, mpiutils, part, vectorutils
!+
!--------------------------------------------------------------------------
module centreofmass
 implicit none
 public :: reset_centreofmass,get_centreofmass,correct_bulk_motion,get_total_angular_momentum

 private

contains

!----------------------------------------------------------------
!+
!  routine to reset the centre of mass in the initial conditions
!  (assuming equal mass particles)
!+
!----------------------------------------------------------------
subroutine reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 use io, only:iprint
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in),    optional :: nptmass
 real,    intent(inout), optional :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 real :: xcom(3),vcom(3),xcomold(3),vcomold(3)
 integer :: i

 if (present(xyzmh_ptmass) .and. present(vxyz_ptmass) .and. present(nptmass)) then
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 else
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu)
 endif

 xcomold = xcom
 vcomold = vcom
 do i=1,npart
    xyzh(1:3,i) = xyzh(1:3,i) - xcom(1:3)
    vxyzu(1:3,i) = vxyzu(1:3,i) - vcom(1:3)
 enddo
 if (present(xyzmh_ptmass) .and. present(vxyz_ptmass) .and. present(nptmass)) then
    do i=1,nptmass
       xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - xcom(1:3)
       vxyz_ptmass(1:3,i)  = vxyz_ptmass(1:3,i) - vcom(1:3)
    enddo
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 else
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu)
 endif
 write(iprint,"(' reset CofM: (',3(es9.2,1x),') -> (',3(es9.2,1x),')')") xcomold,xcom

 return
end subroutine reset_centreofmass

!----------------------------------------------------------------
!+
! Routine returns the centre of mass and centre of mass velocity
! Accounting for sink particles is optional
!+
!----------------------------------------------------------------
subroutine get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,mass)
 use io,       only:id,master
 use part,     only:massoftype,iamtype,iphase,igas,maxphase,maxp,isdead_or_accreted
 use mpiutils, only:reduceall_mpi
 real,         intent(out) :: xcom(3),vcom(3)
 integer,      intent(in)  :: npart
 real,         intent(in)  :: xyzh(:,:),vxyzu(:,:)
 integer,      intent(in),  optional :: nptmass
 real,         intent(in),  optional :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 real,         intent(out), optional :: mass
 integer :: i,itype
 real :: xi,yi,zi,hi
 real(kind=8) :: xpos,ypos,zpos,vxpos,vypos,vzpos
 real(kind=8) :: dm,pmassi,totmass

 xpos  = 0.d0
 ypos  = 0.d0
 zpos  = 0.d0
 vxpos = 0.d0
 vypos = 0.d0
 vzpos = 0.d0
 totmass = 0.d0
!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,iphase,massoftype) &
!$omp private(i,itype,xi,yi,zi,hi,pmassi) &
!$omp reduction(+:xpos,ypos,zpos,vxpos,vypos,vzpos,totmass)
!$omp do
 do i=1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
          if (itype > 0) then ! avoid problems if called from ICs
             pmassi = massoftype(itype)
          else
             pmassi = massoftype(igas)
          endif
       else
          pmassi = massoftype(igas)
       endif
       totmass = totmass + pmassi
       xpos    = xpos  + pmassi*xi
       ypos    = ypos  + pmassi*yi
       zpos    = zpos  + pmassi*zi
       vxpos   = vxpos + pmassi*vxyzu(1,i)
       vypos   = vypos + pmassi*vxyzu(2,i)
       vzpos   = vzpos + pmassi*vxyzu(3,i)
    endif
 enddo
!$omp enddo
!$omp end parallel
 if (id==master .and. present(xyzmh_ptmass) .and. present(vxyz_ptmass) .and. present(nptmass)) then
    do i=1,nptmass
       pmassi = xyzmh_ptmass(4,i)
       totmass = totmass + pmassi
       xpos  = xpos  + pmassi*xyzmh_ptmass(1,i)
       ypos  = ypos  + pmassi*xyzmh_ptmass(2,i)
       zpos  = zpos  + pmassi*xyzmh_ptmass(3,i)
       vxpos = vxpos + pmassi*vxyz_ptmass(1,i)
       vypos = vypos + pmassi*vxyz_ptmass(2,i)
       vzpos = vzpos + pmassi*vxyz_ptmass(3,i)
    enddo
 endif
 xcom = (/xpos,ypos,zpos/)
 vcom = (/vxpos,vypos,vzpos/)
 xcom = reduceall_mpi('+',xcom)
 vcom = reduceall_mpi('+',vcom)
 totmass = reduceall_mpi('+',totmass)

 if (totmass > tiny(totmass)) then
    dm = 1.d0/totmass
 else
    dm = 0.d0
 endif
 xcom = xcom*dm
 vcom = vcom*dm

 if (present(mass)) mass = totmass

 return
end subroutine get_centreofmass

!----------------------------------------------------------------
!+
!  Subroutine to compute and correct the bulk motion
!
!  This is really for use with turbulent driving routines
!  which give a net motion to the flow.
!+
!----------------------------------------------------------------
subroutine correct_bulk_motion()
 use dim,      only:maxp
 use part,     only:npart,xyzh,vxyzu,fxyzu,iamtype,maxphase,igas,iphase,&
                    nptmass,xyzmh_ptmass,vxyz_ptmass,isdead_or_accreted,&
                    massoftype
 use mpiutils, only:reduceall_mpi
 use io,       only:iprint,iverbose,id,master
 real    :: totmass,pmassi,hi,xmom,ymom,zmom
 real    :: fmeanx,fmeany,fmeanz,totmass1
 integer :: i

 xmom   = 0.
 ymom   = 0.
 zmom   = 0.
 fmeanx = 0.
 fmeany = 0.
 fmeanz = 0.
 totmass = 0.

!$omp parallel default(none) &
!$omp shared(xyzh,vxyzu,fxyzu,npart) &
!$omp shared(massoftype,iphase) &
!$omp shared(xyzmh_ptmass,vxyz_ptmass,nptmass) &
!$omp private(i,pmassi,hi) &
!$omp reduction(+:fmeanx,fmeany,fmeanz) &
!$omp reduction(+:xmom,ymom,zmom,totmass)
!$omp do
 do i=1,npart
    hi = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       if (maxphase==maxp) then
          pmassi = massoftype(iamtype(iphase(i)))
       else
          pmassi = massoftype(igas)
       endif
       totmass = totmass + pmassi

       xmom = xmom + pmassi*vxyzu(1,i)
       ymom = ymom + pmassi*vxyzu(2,i)
       zmom = zmom + pmassi*vxyzu(3,i)

       fmeanx = fmeanx + pmassi*fxyzu(1,i)
       fmeany = fmeany + pmassi*fxyzu(2,i)
       fmeanz = fmeanz + pmassi*fxyzu(3,i)
    endif
 enddo
!$omp enddo
!
!--add contribution from sink particles
!
!$omp do
 do i=1,nptmass
    pmassi  = xyzmh_ptmass(4,i)
    totmass = totmass + pmassi
    xmom    = xmom + pmassi*vxyz_ptmass(1,i)
    ymom    = ymom + pmassi*vxyz_ptmass(2,i)
    zmom    = zmom + pmassi*vxyz_ptmass(3,i)
 enddo
!$omp enddo
!$omp end parallel

 xmom = reduceall_mpi('+',xmom)
 ymom = reduceall_mpi('+',ymom)
 zmom = reduceall_mpi('+',zmom)
 fmeanx = reduceall_mpi('+',fmeanx)
 fmeany = reduceall_mpi('+',fmeany)
 fmeanz = reduceall_mpi('+',fmeanz)
 totmass = reduceall_mpi('+',totmass)

 totmass1 = 1./totmass
 xmom     = xmom*totmass1
 ymom     = ymom*totmass1
 zmom     = zmom*totmass1
 fmeanx   = fmeanx*totmass1
 fmeany   = fmeany*totmass1
 fmeanz   = fmeanz*totmass1

 if (id==master .and. iverbose >= 2) then
    write(iprint,*) ' correcting bulk motion ',xmom,ymom,zmom
 endif
 !$omp parallel default(none) &
 !$omp shared(npart,xyzh,vxyzu,fxyzu,xmom,ymom,zmom,fmeanx,fmeany,fmeanz) &
 !$omp shared(nptmass,vxyz_ptmass) &
 !$omp private(i)
 !$omp do schedule(static)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       vxyzu(1,i) = vxyzu(1,i) - xmom
       vxyzu(2,i) = vxyzu(2,i) - ymom
       vxyzu(3,i) = vxyzu(3,i) - zmom
       fxyzu(1,i) = fxyzu(1,i) - fmeanx
       fxyzu(2,i) = fxyzu(2,i) - fmeany
       fxyzu(3,i) = fxyzu(3,i) - fmeanz
    endif
 enddo
 !$omp enddo
 !$omp do
 do i=1,nptmass
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) - xmom
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) - ymom
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) - zmom
 enddo
 !$omp enddo
 !$omp end parallel

end subroutine correct_bulk_motion

!------------------------------------------------------------------------
!
! Small routine to calculate the total angular momentum vector of
! the whole system (particles + sinks)
!
!------------------------------------------------------------------------
subroutine get_total_angular_momentum(xyzh,vxyz,npart,L_tot,xyzmh_ptmass,vxyz_ptmass,npart_ptmass)
 use vectorutils, only:cross_product3D
 use part,        only:iphase,iamtype,massoftype,isdead_or_accreted
 use mpiutils,    only:reduceall_mpi
 real, intent(in)  :: xyzh(:,:),vxyz(:,:)
 real, optional, intent(in):: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: npart
 integer, optional, intent(in) :: npart_ptmass
 real, intent(out) :: L_tot(3)
 integer           :: ii,itype
 real              :: temp(3),pmassi

 L_tot(:) = 0.

 ! Calculate the angular momentum from all the particles
 ! Check if particles are dead or have been accreted first
!$omp parallel default(none) &
!$omp shared(xyzh,vxyz,npart) &
!$omp shared(massoftype,iphase) &
!$omp shared(xyzmh_ptmass,vxyz_ptmass,npart_ptmass) &
!$omp private(ii,itype,pmassi,temp) &
!$omp reduction(+:L_tot)
!$omp do
 do ii = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,ii))) then
       itype = iamtype(iphase(ii))
       pmassi = massoftype(itype)
       call cross_product3D(xyzh(1:3,ii),vxyz(1:3,ii),temp)
       L_tot = L_tot + temp*pmassi
    endif
 enddo
!$omp enddo

 ! Calculate from the sinks
 if (present(npart_ptmass)) then
    !$omp do
    do ii = 1,npart_ptmass
       call cross_product3D(xyzmh_ptmass(1:3,ii),vxyz_ptmass(1:3,ii),temp)
       L_tot = L_tot + temp*xyzmh_ptmass(4,ii)
    enddo
    !$omp enddo
 endif
!$omp end parallel

 L_tot = reduceall_mpi('+',L_tot)

end subroutine get_total_angular_momentum

end module centreofmass
