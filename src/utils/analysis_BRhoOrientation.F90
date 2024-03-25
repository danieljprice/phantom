!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to determine orientation between the magnetic field
! vector & density gradients
! Written by Boden Simpson, under the guidance of James Wurster
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: boundary, centreofmass, kernel, part, physcon, sortutils,
!   units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'Orientation'
 public  :: do_analysis

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use centreofmass, only: reset_centreofmass
 use physcon,      only: pi,gg,years
 use part,         only: rhoh,Bxyz
 use units,        only: unit_density,unit_Bfield,unit_velocity
 use kernel,       only: grkern,cnormk
 use sortutils,    only: indexx
#ifdef PERIODIC
 use boundary,     only:dxbound,dybound,dzbound
#endif
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time
 integer, parameter :: nbins      = 100
 real,    parameter :: costmin    =   0.00
 real,    parameter :: costmax    =   1.00
 real,    parameter :: vtmin      =   -1.00
 real,    parameter :: vtmax      =   1.00
 real,    parameter :: rhomin_cgs =   1.0d-25
 real,    parameter :: rhomax_cgs =   1.0d-20
 real,    parameter :: vmax_cgs   =   1.0d7
 real,    parameter :: vmin_cgs   =   1.0d3
 real,    parameter :: Bmin_cgs   =   1.0d-7
 real,    parameter :: Bmax_cgs   =   1.0d-4
 real,    parameter :: omega      = 1.00
 integer            :: lst(npart),ipos(npart),ii,jj,ikount
 integer            :: i,j,k,p,t,b,r,o,l,vt
 integer            :: thetB(nbins,nbins),thetrho(nbins,nbins)
 integer            :: perpi(nbins,nbins),mixedi(nbins,nbins),parali(nbins,nbins)
 integer            :: vcost(nbins),vvt(nbins,nbins), bvt(nbins,nbins),cost(nbins)
 real               :: dpos(npart)
 real               :: xi,yi,zi,hi,Bxi,Byi,Bzi,vxi,vyi,vzi
 real               :: rhxi,rhyi,rhzi,xj,yj,zj,hj,dxi,dyi,dzi,dvt
 real               :: rhomin,rhomax,Bmin,Bmax,vmax,vmin
 real               :: q,dri,rhoi,rhoi1,rhoj,twohi
 real               :: dB,dcost,drho,dv
 real               :: grki,grkxi,grkzi,grkyi
 real               :: absB,absrho,costheta,vtheta,absV
 real               :: rhobins(nbins),Bbins(nbins),costbins(nbins),vbins(nbins),vtbins(nbins)
 real               :: mixedavg(nbins,nbins),paralavg(nbins,nbins),perpavg(nbins,nbins)
 logical            :: keep_searching
 character(len=200) :: fileout
 !
 !-- Initialise parameters
 !-- Converting cgs units to code units
 !
 rhomin = rhomin_cgs/unit_density
 rhomax = rhomax_cgs/unit_density
 Bmin   = Bmin_cgs/unit_Bfield
 Bmax   = Bmax_cgs/unit_Bfield
 vmax   = vmax_cgs/unit_velocity
 vmin   = vmin_cgs/unit_velocity
 rhobins   = 0.
 Bbins     = 0.
 costbins  = 0.
 vtbins    = 0.
 vbins     = 0.
 thetB     = 0.
 thetrho   = 0.
 vcost     = 0.
 mixedavg  = 0.
 paralavg  = 0.
 perpavg   = 0.
 parali    = 0.
 perpi     = 0.
 mixedi    = 0.
 vvt       = 0.
 bvt       = 0.
 vcost     = 0.
 cost      = 0.

 !--Set bin ranges, log for B/rho and linear for costheta
 drho  = (log10(rhomax) - log10(rhomin))/float(nbins)
 dB    = (log10(Bmax)   - log10(Bmin))/float(nbins)
 dv    = (log10(vmax)   - log10(vmin))/float(nbins)
 dcost = (costmax - costmin)/float(nbins)
 dvt   = (vtmax   - vtmin  )/float(nbins)

 do i = 1,nbins
    rhobins(i)  = 10**(log10(rhomin) + (i-1)*drho)
    Bbins(i)    = 10**(log10(Bmin) + (i-1)*dB)
    vbins(i)    = 10**(log10(vmin) + (i-1)*dv)
    costbins(i) = (costmin + (i-1)*dcost)
    vtbins(i)   = (vtmin   + (i-1)*dvt)
 enddo

 !--Sorting all particles into list ordered by ascending Z position
 !  Used to find the neighbouring particles without the full neighbour-finding process
 ikount = 0
 do i = 1,npart
    if (xyzh(4,i) > 0) then
       ikount = ikount + 1
       ipos(ikount) = i
       dpos(ikount) = xyzh(3,i)
    endif
 enddo
 call indexx(ikount,dpos,lst)

 !$omp parallel default(none) &
 !$omp shared(npart,xyzh,particlemass,Bxyz,costbins,Bbins,rhobins,Bmin,rhomin,unit_density,ikount,vxyzu) &
 !$omp shared(ipos,lst,vbins,vtbins) &
#ifdef PERIODIC
!$omp shared(dxbound,dybound,dzbound) &
#endif
 !$omp private(i,xi,yi,zi,hi,rhoi,rhoi1,Bxi,Byi,Bzi,rhxi,rhyi,rhzi,xj,yj,zj,hj,j,dxi,dyi,dzi,dri,q,rhoj,l,vt) &
 !$omp private(grki,grkxi,grkyi,grkzi,absB,absrho,costheta,k,p,t,b,r,o,ii,jj,twohi,vxi,vyi,vzi,vtheta,absV) &
 !$omp private(keep_searching) &
 !$omp reduction(+:thetB,thetrho,vvt,bvt,vcost,cost) &
 !$omp reduction(+:paralavg,perpavg,mixedavg,perpi,parali,mixedi)
 !$omp do schedule(runtime)
 aparts: do ii = 1,ikount
    ! properties of particle i
    i  = ipos(lst(ii))
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (hi < tiny (hi)) cycle aparts      ! skip dead particles
    twohi = 2.0*hi
    rhoi  = rhoh(hi, particlemass)
    rhoi1 = 1.0/rhoi
    Bxi   = Bxyz(1,i)
    Byi   = Bxyz(2,i)
    Bzi   = Bxyz(3,i)
    vxi   = vxyzu(1,i)
    vyi   = vxyzu(2,i)
    vzi   = vxyzu(3,i)
    rhxi  = 0.0
    rhyi  = 0.0
    rhzi  = 0.0

    ! Quick algorithm to find neighouring particles within 2h.
    ! Particles have been sorted by z, so start at the location of i and search in both directions
    ! until abs(dz) > 2h.  This is much faster than the N^2 loop, but easier to code since it does
    ! not require the full neighbour finding algorithm
    bparts: do p = 1,2
       jj = ii

       keep_searching = .true.
       do while (keep_searching)
          if (p==1) then
             jj = jj-1
             if (jj == 0) jj = ikount    ! if we go off the bottom of the list, loop back to the top (required for periodicity)
          else
             jj = jj+1
             if (jj == ikount+1) jj = 1  ! if we go off the top of the list, loop back to the bottom (required for periodicity)
          endif
          j   = ipos(lst(jj))
          dzi = zi - xyzh(3,j)
#ifdef PERIODIC
          if (abs(dzi) > 0.5*dzbound) dzi = dzi - dzbound*SIGN(1.0,dzi)
#endif
          if (abs(dzi) > twohi) then
             keep_searching = .false.
          else
             dxi = xi-xyzh(1,j)
             dyi = yi-xyzh(2,j)
#ifdef PERIODIC
             if (abs(dxi) > 0.5*dzbound) dxi = dxi - dzbound*SIGN(1.0,dxi)
             if (abs(dyi) > 0.5*dzbound) dyi = dyi - dzbound*SIGN(1.0,dyi)
#endif
             dri = sqrt(dxi*dxi + dyi*dyi + dzi*dzi)

             if (twohi > dri) then
                hj   = xyzh(4,j)
                q    = hi/dri
                rhoj = rhoh(hj, particlemass)

                !Grad of the kernel
                grki  = cnormk * grkern(q*q,q) / (dri*hi**4)
                grkxi = grki * dxi
                grkyi = grki * dyi
                grkzi = grki * dzi

                !x,y,z components of grad density
                rhxi = rhxi + particlemass * rhoi1 * (rhoj - rhoi) * grkxi
                rhyi = rhyi + particlemass * rhoi1 * (rhoj - rhoi) * grkyi
                rhzi = rhzi + particlemass * rhoi1 * (rhoj - rhoi) * grkzi
             endif
          endif
       enddo
    enddo bparts
    ! Absolute value of grad density and B field
    absB   = sqrt(Bxi*Bxi   + Byi*Byi   + Bzi*Bzi)
    absV   = sqrt(vxi*vxi   + vyi*vyi   + vzi*vzi)
    absrho = sqrt(rhxi*rhxi + rhyi*rhyi + rhzi*rhzi)

    ! Angle between grad density and B field
    costheta = abs((Bxi*rhxi + Byi*rhyi + Bzi*rhzi) / (absB*absrho))
    vtheta   = (Bxi*vxi  + Byi*vyi  + Bzi*vzi) / (absB*absV)

    ! Finding bins
    t  = 1 ! cosTheta/angle
    b  = 1 ! mag/B field
    r  = 1 ! rho/density
    l  = 1 ! velocity
    vt = 1 ! psi angle
    do while (costheta > costbins(t) .and. t < nbins)
       t = t + 1
    enddo
    do while (vtheta > vtbins(vt) .and. vt < nbins)
       vt = vt + 1
    enddo
    do while (absB > Bbins(b) .and. b < nbins)
       b = b + 1
    enddo
    do while (rhoi > rhobins(r) .and. r < nbins)
       r = r + 1
    enddo
    do while (absV > vbins(l) .and. l < nbins)
       l = l + 1
    enddo

    ! Binning particles (B-costheta, rho-costheta, rho-B plane and by orientation)
    if (t < nbins) then
       cost(t) = cost(t) + 1
       if (b < nbins) then
          thetB(t,b) = thetB(t,b) + 1
          if (r < nbins) then
             !-- Binning by orientation, perpendicular, parallel and mixed - not being used
             if (costheta > 0.0 .and. costheta < 0.4) then
                perpavg(b,r) = perpavg(b,r) + costheta
                perpi(b,r)   = perpi(b,r)   + 1
             endif
             if (costheta > 0.4 .and. costheta < 0.6) then
                mixedavg(b,r) = mixedavg(b,r) + costheta
                mixedi(b,r)   = mixedi(b,r)   + 1
             endif
             if (costheta > 0.6 .and. costheta < 1.0) then
                paralavg(b,r) = paralavg(b,r) + costheta
                parali(b,r)   = parali(b,r)   + 1
             endif
          endif
       endif
       if (r < nbins) then
          thetrho(t,r) = thetrho(t,r) + 1
       endif
    endif

    if (vt < nbins) then
       if (b < nbins) then
          bvt(vt,b) = bvt(vt,b) + 1
       endif
       if (l < nbins) then
          vvt(vt,l) = vvt(vt,l) + 1
       endif
       vcost(vt) = vcost(vt) + 1
    endif
 enddo aparts
 !$omp enddo
 !$omp end parallel

 !--calculate averages
 do i = 1,nbins
    do j = 1,nbins
       if (perpi(i,j) > 0) then
          perpavg(i,j) = perpavg(i,j)/real(perpi(i,j))
       endif
       if (mixedi(i,j) > 0) then
          mixedavg(i,j) = mixedavg(i,j)/real(mixedi(i,j))
       endif
       if (parali(i,j) > 0) then
          paralavg(i,j) = paralavg(i,j)/real(parali(i,j))
       endif
    enddo
 enddo

 !--Write results to file
 write(fileout,'(3a)') 'analysisout_costB',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
      1,'cost', &
      2,'B',    &
      3,'freq'
 do i = 1,nbins
    do j = 1,nbins
       write(iunit,'(2(1pe18.10,1x),(I18,1x))') costbins(i),Bbins(j)*unit_Bfield,thetB(i,j)
    enddo
 enddo
 close(iunit)

 write(fileout,'(3a)') 'analysisout_costrho',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
      1,'cost',  &
      2,'rho',   &
      3,'freq'
 do i = 1,nbins
    do j = 1,nbins
       write(iunit,'(2(1pe18.10,1x),(I18,1x))') costbins(i),rhobins(j)*unit_density,thetrho(i,j)
    enddo
 enddo
 close(iunit)

 write(fileout,'(3a)') 'analysisout_mixedBrho',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
      1,'B',     &
      2,'rho',   &
      3,'perp',  &
      4,'mixed', &
      5,'paral'
 do i = 1,nbins
    do j = 1,nbins
       write(iunit,'(2(1pe18.10,1x),3(I18,1x))') Bbins(i)*unit_Bfield,rhobins(j)*unit_density,perpi(i,j),mixedi(i,j),parali(i,j)
    enddo
 enddo
 close(iunit)

 write(fileout,'(3a)') 'analysisout_vcostB',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
      1,'vcost', &
      2,'B',     &
      3,'freq'
 do i = 1,nbins
    do j = 1,nbins
       write(iunit,'(2(1pe18.10,1x),(I18,1x))') vtbins(i),Bbins(j)*unit_Bfield,bvt(i,j)
    enddo
 enddo
 close(iunit)

 write(fileout,'(3a)') 'analysisout_vcostv',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
      1,'vcost', &
      2,'v',     &
      3,'freq'
 do i = 1,nbins
    do j = 1,nbins
       write(iunit,'(2(1pe18.10,1x),(I18,1x))') vtbins(i),vbins(j)*unit_velocity,vvt(i,j)
    enddo
 enddo
 close(iunit)

 write(fileout,'(3a)') 'analysisout_hist',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',4(1x,'[',i2.2,1x,a11,']',2x))") &
      1,'vcost',   &
      2,'freq' ,   &
      3,'cost' ,   &
      4,'freq'
 do i = 1,nbins
    write(iunit,'((1pe18.10,1x),(I18,1x),(1pe18.10,1x),(I18,1x))') vtbins(i), vcost(i), costbins(i), cost(i)
 enddo
 close(iunit)
end subroutine do_analysis
!-----------------------------------------------------------------------

end module analysis
