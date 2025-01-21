!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine calculating to determine the radial profile of a sphere
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, dim, part, physcon, units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'sphere'
 logical            :: firstcall     = .true.
 real,    parameter :: rhothresh_cgs = 2.126d-03*0.0 ! cgs value above which we will calculate the values
 integer, parameter :: nbins         = 512           ! number of bins for radial profile
 logical            :: logr          = .true.        ! log radial profile; else linear
 logical            :: centre_rhomax = .true.        ! centre on the densest particle, otherwise reset to the centre of mass
 real,    parameter :: rmin_pc       =  1.d-6        ! minimum radial bin
 real,    parameter :: rmax_pc       =  1.100        ! maximum radial bin

 public  :: do_analysis

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass
 use physcon,      only: pi,gg,years,pc
 use part,         only: igas,iamtype,iphase,maxphase,rhoh
 use units,        only: umass,udist,utime,unit_density
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time
 integer            :: i,j,itype,icom
 integer            :: ibins(3,nbins)
 real               :: dr,dr2,vr,rad,rad2d,total_mass,total_ang,total_sang
 real               :: xi,yi,zi,hi,vxi,vyi,vzi,vphi,ui,rhoi
 real               :: angx,angy,angz,angtot,angi,angxi,angyi,angzi,angp,angm
 real               :: rhothresh,rmin,rmax,rhomax,hmin,twoh2
 real               :: rbins(nbins),vbins(8,nbins),xcom(3),xdense(3)
 logical            :: iexist
 character(len=200) :: fileout

 !--Initialise parameters
 rbins = 0.0
 vbins = 0.0
 angx  = 0.0
 angy  = 0.0
 angz  = 0.0
 angp  = 0.0
 angm  = 0.0
 ibins = 0
 rhothresh = rhothresh_cgs/unit_density
 rmin      = rmin_pc*pc/udist
 rmax      = rmax_pc*pc/udist

 !--Set bins (linear or log)
 if (logr ) then
    dr = (log10(rmax)-log10(rmin))/float(nbins)
    do i = 1,nbins
       rbins(i) = 10**(log10(rmin) +(i-1)*dr )
    enddo
 else
    dr   = rmax/float(nbins)
    do i = 1,nbins
       rbins(i) = float(i)*dr
    enddo
 endif

 !--Reset centre
 if (centre_rhomax) then
    hmin = huge(hmin)
    xdense = 0.
    ! find densest particle
    do i = 1,npart
       hi = xyzh(4,i)
       if (hi > tiny(hi) .and. hi < hmin) then
          hmin = hi
          xdense = xyzh(1:3,i)
       endif
    enddo
    ! find centre of mass of nearby dense particles
    xcom   = 0.
    icom   = 0
    twoh2  = 4.0*hmin*hmin
    rhomax = rhoh(hmin,particlemass)
    print*, 'rhomax = ',rhomax*unit_density
    do i = 1,npart
       hi = xyzh(4,i)
       if (hi > tiny(hi)) then
          rhoi = rhoh(hi,particlemass)
          if (rhoi > 0.1*rhomax) then
             dr2 = dot_product((xdense-xyzh(1:3,i)),(xdense-xyzh(1:3,i)))
             if (dr2 < twoh2 .and. .true.) then
                xcom(1:3) = xcom(1:3) + xyzh(1:3,i)
                icom      = icom + 1
             endif
          endif
       endif
    enddo
    xcom = xcom/icom
    print*, 'Using ',icom,' particles to move the centre to ',xcom
    ! adjust centres
    do i = 1,npart
       xyzh(1:3,i) = xyzh(1:3,i) - xcom
    enddo
 else
    ! move to centre of mass
    call reset_centreofmass(npart,xyzh,vxyzu)
 endif

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
    rhoi = rhoh(hi,particlemass)
    if (rhoi > rhothresh) then
       rad   = sqrt(xi*xi + yi*yi + zi*zi)
       rad2d = sqrt(xi*xi + yi*yi)
       vxi   = vxyzu(1,i)
       vyi   = vxyzu(2,i)
       vzi   = vxyzu(3,i)
       ui    = vxyzu(4,i)

       angxi = yi*vzi - zi*vyi
       angyi = zi*vxi - xi*vzi
       angzi = xi*vyi - yi*vxi
       angi  = particlemass*sqrt(angxi**2 + angyi**2 + angzi**2)
       if (rad2d > 0.) then
          vphi = angzi/sqrt(xi*xi + yi*yi)
       else
          vphi = 0.0
       endif
       if (rad > 0.) then
          vr = (xi*vxi + yi*vyi + zi*vzi)/rad
       else
          vr = 0.0
       endif
       ! sum for total values
       angx = angx + angxi
       angy = angy + angyi
       angz = angz + angzi
       if (vphi > 0.0) then
          angp = angp + vphi*rad2d
       else
          angm = angm + vphi*rad2d
       endif

       ! find the bin
       j = 1
       do while (rad > rbins(j) .and. j<nbins)
          j = j + 1
       enddo
       if (j < nbins) then
          ibins(1,j) = ibins(1,j) + 1
          vbins(1,j) = vbins(1,j) + angi
          vbins(2,j) = vbins(2,j) + rhoi
          vbins(3,j) = vbins(3,j) + hi
          vbins(4,j) = vbins(4,j) + ui
          vbins(5,j) = vbins(5,j) + vr
          vbins(6,j) = vbins(6,j) + vphi
          if (vphi > 0.) then
             ibins(2,j) = ibins(2,j) + 1
             vbins(7,j) = vbins(7,j) + vphi
          elseif (vphi < 0.) then
             ibins(3,j) = ibins(3,j) + 1
             vbins(8,j) = vbins(8,j) + vphi
          endif
       endif
    endif
 enddo parts
 angtot = sqrt(angx*angx + angy*angy + angz*angz)

 !--Write results to file
 write(fileout,'(3a)') 'analysisout_',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'r',       &
       2,'M(<r)',   &
       3,'L(<r)',   &
       4,'l(<r)',   &
       5,'rho',     &
       6,'h',       &
       7,'u',       &
       8,'v_r',     &
       9,'v_phi',   &
      10,'v_phi>0', &
      11,'v_phi<0'

 total_mass = 0.0
 total_ang  = 0.0
 total_sang = 0.0
 do i = 1,nbins-1
    total_mass = total_mass + ibins(1,i)*particlemass
    total_ang  = total_ang  + vbins(1,i)
    if (ibins(1,i) > 0) total_sang = total_sang + vbins(1,i)/(ibins(1,i)*particlemass)
    if (ibins(1,i) > 0) vbins(2:6,i) = vbins(2:6,i)/ibins(1,i)
    if (ibins(2,i) > 0) vbins(  7,i) = vbins(  7,i)/ibins(2,i)
    if (ibins(3,i) > 0) vbins(  8,i) = vbins(  8,i)/ibins(3,i)
    write(iunit,'(11(1pe18.10,1x))') 0.5*(rbins(i)+rbins(i+1)),total_mass,total_ang,total_sang,vbins(2:8,i)

 enddo
 close(iunit)

 fileout = trim(dumpfile(1:index(dumpfile,'_')-1))//'_AM.dat'
 inquire(file=fileout,exist=iexist)
 if ( .not.iexist .or. firstcall ) then
    firstcall = .false.
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time (code)',&
          2,'time (kyr)', &
          3,'L    (com)', &
          4,'L_x  (com)', &
          5,'L_y  (com)', &
          6,'L_z  (com)', &
          7,'L_z+ (com)', &
          8,'L_z- (com)'
 else
    open(iunit,file=fileout,position='append')
 endif
 angtot = particlemass*angtot * umass*udist**2 / utime
 angx   = particlemass*angx   * umass*udist**2 / utime
 angy   = particlemass*angy   * umass*udist**2 / utime
 angz   = particlemass*angz   * umass*udist**2 / utime
 angp   = particlemass*angp   * umass*udist**2 / utime
 angm   = particlemass*angm   * umass*udist**2 / utime

 write(iunit,'(14(es18.10,1x))') time,time*utime/years/1000.0,angtot,angx,angy,angz,angp,angm

 close(iunit)
end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module analysis

