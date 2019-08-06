!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine calculating to determine the radial profile of a sphere
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, dim, part, physcon, units
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'sphere'
 logical :: firstcall = .true.
 real, parameter :: rhothresh = 2.126d-03*0.0 ! density threshhold in code units
 public  :: do_analysis

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass
 use physcon,      only: pi,gg,years
 use part,         only: igas,iamtype,iphase,maxphase,rhoh
 use units,        only: umass,udist,utime
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time
 integer, parameter :: nbins = 256
 real,    parameter :: rmin  =   0.01
 real,    parameter :: rmax  =   2.00
 logical            :: logr  = .true.
 integer            :: i,j,itype
 integer            :: ibins(nbins)
 real               :: dr,total_mass,density,vr,rad,velocity
 real               :: xi,yi,zi,hi,vxi,vyi,vzi,ui
 real               :: angx,angy,angz,angtot,angtota,angxa,angya,angza
 real               :: angp,angpa,angm,angma,rvphi
 real               :: rbins(nbins),ubins(nbins),vbins(nbins),hbins(nbins)
 logical            :: iexist
 character(len=200) :: fileout
 !
 !-- Initialise parameters
 !
 ubins = 0.0
 vbins = 0.0
 hbins = 0.0
 rbins = 0.0
 angx  = 0.0
 angy  = 0.0
 angz  = 0.0
 angxa = 0.0
 angya = 0.0
 angza = 0.0
 angp  = 0.0
 angpa = 0.0
 angm  = 0.0
 angma = 0.0
 ibins = 0
 !
 !--Calculate the angular momentum about the origin
 !
 aparts: do i = 1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (hi < tiny (hi)) cycle aparts      ! dead particle
    if (maxphase==maxp) then
       itype = iamtype(iphase(i))
    else
       itype = igas
    endif
    if (itype/=igas) cycle aparts         ! not gas
    vxi = vxyzu(1,i)
    vyi = vxyzu(2,i)
    vzi = vxyzu(3,i)
    rvphi = xi*vyi - yi*vxi
    ! angular momentum above a given density threshhold
    if (rhoh(hi,particlemass) > rhothresh) then
       angxa = angxa + particlemass*(yi*vzi - zi*vyi)
       angya = angya + particlemass*(zi*vxi - xi*vzi)
       angza = angza + particlemass*(xi*vyi - yi*vxi)
       if (rvphi > 0.0) then
          angpa = angpa + particlemass*rvphi
       else
          angma = angma + particlemass*rvphi
       endif
    endif
 enddo aparts
 angtota = sqrt(angxa*angxa + angya*angya + angza*angza)
 !
 !--Set bins (log or linear)
 if ( logr ) then
    dr   = rmax/float(nbins)
    do i = 1,nbins
       rbins(i) = float(i)*dr
    enddo
 else
    dr = (log10(rmax)-log10(rmin))/float(nbins)
    do i = 1,nbins
       rbins(i) = 10**(log10(rmin) +(i-1)*dr )
    enddo
 endif
 !
 !--Reset centre of mass
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 !--Bin the data
 parts: do i = 1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (hi < tiny (hi)) cycle parts      ! dead particle
    if (maxphase==maxp) then
       itype = iamtype(iphase(i))
    else
       itype = igas
    endif
    if (itype/=igas) cycle parts         ! not gas
    vxi = vxyzu(1,i)
    vyi = vxyzu(2,i)
    vzi = vxyzu(3,i)
    ui  = vxyzu(4,i)
    rvphi = xi*vyi - yi*vxi
    ! angular momentum above a given density threshhold
    if (rhoh(hi,particlemass) > rhothresh) then
       angx = angx + particlemass*(yi*vzi - zi*vyi)
       angy = angy + particlemass*(zi*vxi - xi*vzi)
       angz = angz + particlemass*(xi*vyi - yi*vxi)
       if (rvphi > 0.0) then
          angp = angp + particlemass*rvphi
       else
          angm = angm + particlemass*rvphi
       endif
    endif
    !--Calculate properties of the particle
    rad      = sqrt(xi*xi + yi*yi + zi*zi)
    velocity = sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
    vr       = (xi*vxi + yi*vyi + zi*vzi)/rad
    !--Find the bin
    j = 1
    do while (rad > rbins(j) .and. j<nbins)
       j = j + 1
    enddo
    if (j < nbins) then
       ibins(j) = ibins(j) + 1
       ubins(j) = ubins(j) + ui
       vbins(j) = vbins(j) + vr
       hbins(j) = hbins(j) + hi
    endif
 enddo parts
 angtot = sqrt(angx*angx + angy*angy + angz*angz)
 !
 !--Write results to file
 write(fileout,'(3a)') 'analysisout_',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'r',   &
       2,'M(r)',&
       3,'rho', &
       4,'u',   &
       5,'v_r', &
       6,'h'
 total_mass = 0.0
 do i = 1,nbins
    total_mass = total_mass + ibins(i)*particlemass
    density    = rhoh(hbins(i)/ibins(i),particlemass)
    write(iunit,'(6(1pe18.10,1x))') rbins(i),total_mass,density, &
    ubins(i)/ibins(i),vbins(i)/ibins(i),hbins(i)/ibins(i)
 enddo
 close(iunit)
 !
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_AM.dat'
 inquire(file=fileout,exist=iexist)
 if ( .not.iexist .or. firstcall ) then
    firstcall = .false.
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',14(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time (code)',&
          2,'time (kyr)', &
          3,'L   (com)',  &
          4,'L_x (com)',  &
          5,'L_y (com)',  &
          6,'L_z (com)',  &
          7,'L+ (com)',   &
          8,'L- (com)',   &
          9,'L   (act)',  &
         10,'L_x (act)',  &
         11,'L_y (act)',  &
         12,'L_z (act)',  &
         13,'L+ (act)',   &
         14,'L- (act)'
 else
    open(iunit,file=fileout,position='append')
 endif
 angtot  = angtot  * umass*udist**2 / utime
 angx    = angx    * umass*udist**2 / utime
 angy    = angy    * umass*udist**2 / utime
 angz    = angz    * umass*udist**2 / utime
 angtota = angtota * umass*udist**2 / utime
 angxa   = angxa   * umass*udist**2 / utime
 angya   = angya   * umass*udist**2 / utime
 angza   = angza   * umass*udist**2 / utime

 angpa   = angpa   * umass*udist**2 / utime
 angma   = angma   * umass*udist**2 / utime
 angp    = angp    * umass*udist**2 / utime
 angm    = angm    * umass*udist**2 / utime

 write(iunit,'(14(es18.10,1x))') time,time*utime/years/1000.0, &
                                 angtot,angx,angy,angz,angp,angm, &
                                 angtota,angxa,angya,angza,angpa,angma
 close(iunit)
 !
end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
