!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION: Traces the centre of (stellar) mass each galaxy that
!  that is part of a major merger.  Stars initially assigned to each
!  galaxy are assumed to remain in that galaxy for all time.
!  Also output max and average densities for gas, stars and dark matter
!
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, physcon, units
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'GalMerger'
 logical, private :: firstcall = .true.

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,         only: iphase,istar,idarkmatter,igas, &
                         massoftype,maxp,maxphase,isdead_or_accreted,iamtype,rhoh
 use units,        only: utime,udist,unit_density
 use physcon,      only: years,mpc
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i,itype,nstar,ngas,ndm
 real                         :: xi,yi,zi,hi,pmassi,time0,rhoi
 real                         :: xpos1,ypos1,zpos1,xpos2,ypos2,zpos2,totmass1,totmass2,r
 real                         :: rhogasX,rhostarX,rhodmX,rhogasA,rhostarA,rhodmA
 logical                      :: iexist
 character(len=200)           :: fileout
 !
 !--Open file (appendif exists)
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_stellarCoM.dat'
 inquire(file=fileout,exist=iexist)
 if ( .not.iexist .or. firstcall ) then
    firstcall = .false.
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',14(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',       &
          2,'x1',         &
          3,'y1',         &
          4,'z1',         &
          5,'x2',         &
          6,'y2',         &
          7,'z2',         &
          8,'r',          &
          9,'rho_gas max',&
         10,'rho_str max',&
         11,'rho_dm  max',&
         12,'rho_gas ave',&
         13,'rho_str ave',&
         14,'rho_dm  ave'
 else
    open(iunit,file=fileout,position='append')
 endif
 !
 xpos1    = 0.d0
 ypos1    = 0.d0
 zpos1    = 0.d0
 totmass1 = 0.d0
 xpos2    = 0.d0
 ypos2    = 0.d0
 zpos2    = 0.d0
 totmass2 = 0.d0
 rhogasX   = 0.d0
 rhostarX  = 0.d0
 rhodmX    = 0.d0
 rhogasA   = 0.d0
 rhostarA  = 0.d0
 rhodmA    = 0.d0
 nstar     = 0
 ngas      = 0
 ndm       = 0
!$omp parallel default(none) &
!$omp shared(npart,xyzh,iphase,massoftype) &
!$omp private(i,itype,xi,yi,zi,hi,pmassi,rhoi) &
!$omp reduction(+:xpos1,ypos1,zpos1,xpos2,ypos2,zpos2,totmass1,totmass2) &
!$omp reduction(+:rhogasA,rhostarA,rhodmA,nstar,ngas,ndm) &
!$omp reduction(max:rhogasX,rhostarX,rhodmX)
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
       rhoi = rhoh(hi,pmassi)
       if (itype==istar) then
          rhostarX = max(rhostarX,rhoi)
          rhostarA =     rhostarA+rhoi
          nstar    = nstar + 1
          if (i <= npart/2 ) then
             totmass1 = totmass1 + pmassi
             xpos1    = xpos1    + pmassi*xi
             ypos1    = ypos1    + pmassi*yi
             zpos1    = zpos1    + pmassi*zi
          else
             totmass2 = totmass2 + pmassi
             xpos2    = xpos2    + pmassi*xi
             ypos2    = ypos2    + pmassi*yi
             zpos2    = zpos2    + pmassi*zi
          endif
       else if (itype==idarkmatter) then
          rhodmX  = max(rhodmX,rhoi)
          rhodmA  =     rhodmA+rhoi
          ndm     = ndm + 1
       else if (itype==igas) then
          rhogasX = max(rhogasX,rhoi)
          rhogasA =     rhogasA+rhoi
          ngas    = ngas + 1
       endif
    endif
 enddo
!$omp enddo
!$omp end parallel
 !
 xpos1  = xpos1/totmass1 * udist/mpc + 0.5
 ypos1  = ypos1/totmass1 * udist/mpc + 0.5
 zpos1  = zpos1/totmass1 * udist/mpc + 0.5
 xpos2  = xpos2/totmass2 * udist/mpc + 0.5
 ypos2  = ypos2/totmass2 * udist/mpc + 0.5
 zpos2  = zpos2/totmass2 * udist/mpc + 0.5
 r      = sqrt( (xpos1-xpos2)**2 + (ypos1-ypos2)**2 + (zpos1-zpos2)**2 )
 time0  = time*utime/years*1.0d-10
 rhogasX = rhogasX *unit_density
 rhostarX= rhostarX*unit_density
 rhodmX  = rhodmX  *unit_density
 rhogasA = rhogasA *unit_density/real(ngas)
 rhostarA= rhostarA*unit_density/real(nstar)
 rhodmA  = rhodmA  *unit_density/real(ndm)
 !
 write(iunit,'(14(es18.10,1x))') time0,xpos1,ypos1,zpos1,xpos2,ypos2,zpos2,r, &
                                 rhogasX,rhostarX,rhodmX,rhogasA,rhostarA,rhodmA
 close(iunit)
 !
end subroutine do_analysis

end module

