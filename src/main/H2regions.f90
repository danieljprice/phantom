!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module HIIRegion
!
! HIIRegion
! Contains routines to model HII region expansion due to ionization (and radiation pressure)..
! Routine made by Hopkins et al. (2012),reused by Fujii et al. (2021)
! Routine made by Dale et al. (2014)
! Ionizing rate computed using observation fit made by Fujii et al. (2021)
! Adapted by Yann Bernard in Phantom
!
! :References: Fujii et al. (2021), Hopkins et al. (2012)
!
! :Owner: Yann Bernard
!
! :Runtime parameters: iH2R (1:Hopkins like,2:Dale like)
!
! :Dependencies: dim, eos, infile_utils, io, linklist, part, physcon,
!   sortutils, timing, units
!
 use eos_HIIR,   only:csion,uIon,Tion,Tcold,muion
 use eos,        only:gmw
 implicit none

 public :: update_ionrates,update_ionrate, HII_feedback,HII_feedback_ray,&
           initialize_H2R,read_options_H2R,write_options_H2R

 integer, parameter, public    :: HIIuprate   = 8 ! update rate when IND_TIMESTEPS=yes (Hopkins like method)
 integer, public               :: iH2R = 0
 real   , public               :: Rmax = 15 ! Maximum HII region radius (pc) to avoid artificial expansion...
 real   , public               :: Mmin = 8  ! Minimum mass (Msun) to produce HII region
 integer, public               :: nHIIsources = 0
 real   , public               :: ar
 real   , public               :: mH
 real   , public               :: mH1
 real   , public               :: DNcoeff

 integer, parameter :: maxcache  = 1024

 real, parameter   :: a = -39.3178 !
 real, parameter   :: b =  221.997 !  fitted parameters to compute
 real, parameter   :: c = -227.456 !  ionisation rate for massive
 real, parameter   :: d =  117.410 !  extracted from Fujii et al. (2021).
 real, parameter   :: e = -30.1511 !  (expressed in function of log(solar masses) and s)
 real, parameter   :: f =  3.06810 !
 real, parameter   :: ar_cgs = 2.7d-13
 real, parameter   :: sigd_cgs = 1.d-21
 real              :: sigd
 real              :: hv_on_c
 real              :: Rst_max
 real              :: Minmass

 private

contains

!-----------------------------------------------------------------------
!+
!  Initialise stellar feedbacks
!+
!-----------------------------------------------------------------------
subroutine initialize_H2R
 use io,      only:iprint,iverbose,id,master
 use part,    only:massoftype
 use units,   only:udist,umass,utime
 use physcon, only:mass_proton_cgs,kboltz,pc,eV,solarm,c
 use eos,     only:gmw

 !calculate the useful constant in code units
 ar = ar_cgs*(utime/udist**3)
 sigd = sigd_cgs*udist**2
 hv_on_c = ((18.6*eV)/c)*(utime/(udist*umass))
 Rst_max = sqrt(((Rmax*pc)/udist)**2)
 Minmass = (Mmin*solarm)/umass
 mH = gmw*mass_proton_cgs/umass
 mH1 = 1./mH
 DNcoeff = log10(massoftype(1)*ar*mH1**2/utime)

 if (id == master .and. iverbose > 1) then
    write(iprint,"(a,es18.10,es18.10)") " feedback constants mH,uIon      : ", mH,uIon
    write(iprint,"(a,es18.10,es18.10)") " Max strögrem radius (code/pc)   : ", Rst_max, Rmax
    write(iprint,"(a,es18.10,es18.10)") " Min feedback mass   (code/Msun) : ", Minmass, Mmin
 endif

end subroutine initialize_H2R

!-----------------------------------------------------------------------
!+
!  Calculation of the the ionizing photon rate of all stars (Only for restart)
!+
!-----------------------------------------------------------------------
subroutine update_ionrates(nptmass,xyzmh_ptmass,h_acc)
 use io,     only:iprint,iverbose
 use units,  only:umass
 use part,   only:irateion,ihacc,irstrom
 use physcon,only:solarm
 integer, intent(in)    :: nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 real,    intent(in)    :: h_acc
 real    :: logmi,log_Q,mi,hi
 integer :: i

 nHIIsources = 0
 !$omp parallel do default(none) &
 !$omp shared(xyzmh_ptmass,iprint,iverbose,umass)&
 !$omp shared(Minmass,h_acc,nptmass)&
 !$omp private(logmi,log_Q,mi,hi)&
 !$omp reduction(+:nHIIsources)
 do i=1,nptmass
    mi = xyzmh_ptmass(4,i)
    hi = xyzmh_ptmass(ihacc,i)
    if (mi > Minmass .and. hi < h_acc) then
       logmi = log10(mi*(umass/solarm))
       ! caluclation of the ionizing photon rate  of each sources
       ! this calculation uses Fujii's formula derived from OSTAR2002 databases
       log_Q = (a+b*logmi+c*logmi**2+d*logmi**3+e*logmi**4+f*logmi**5)
       xyzmh_ptmass(irateion,i) = log_Q
       xyzmh_ptmass(irstrom,i)  = -1.
       nHIIsources = nHIIsources + 1
       if (iverbose >= 0) then
          write(iprint,"(/a,es18.10,es18.10/)") "Massive stars detected : Log Q, Mass : ",log_Q,mi
       endif
    else
       xyzmh_ptmass(irateion,i) = -1.
       xyzmh_ptmass(irstrom,i)  = -1.
    endif
 enddo
 !$omp end parallel do
 if (iverbose > 1) then
    write(iprint,"(/a,i8/)") "nb_feedback sources : ",nHIIsources
 endif

end subroutine update_ionrates

!-----------------------------------------------------------------------
!+
!  update the ionizing photon rate
!+
!-----------------------------------------------------------------------
subroutine update_ionrate(i,xyzmh_ptmass,h_acc)
 use io,     only:iprint,iverbose
 use units,  only:umass
 use part,   only:irateion,ihacc,irstrom
 use physcon,only:solarm
 integer, intent(in)    :: i
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 real,    intent(in)    :: h_acc
 real    :: logmi,log_Q,mi,hi

 mi = xyzmh_ptmass(4,i)
 hi = xyzmh_ptmass(ihacc,i)
 if (mi > Minmass .and. hi < h_acc) then
    logmi = log10(mi*(umass/solarm))
    ! caluclation of the ionizing photon rate  of each sources
    ! this calculation uses Fujii's formula derived from OSTAR2002 databases
    log_Q = (a+b*logmi+c*logmi**2+d*logmi**3+e*logmi**4+f*logmi**5)
    xyzmh_ptmass(irateion,i) =  log_Q
    xyzmh_ptmass(irstrom,i)  = -1.
    nHIIsources = nHIIsources + 1
    if (iverbose >= 0) then
       write(iprint,"(/a,es18.10,es18.10/)") "Massive stars detected : Log Q, Mass : ",log_Q,mi
    endif
 else
    xyzmh_ptmass(irateion,i) = -1.
    xyzmh_ptmass(irstrom,i)  = -1.
 endif

 if (iverbose > 1) then
    write(iprint,"(/a,i8/)") "nb_feedback sources : ",nHIIsources
 endif

end subroutine update_ionrate

!-----------------------------------------------------------------------
!+
! HII feedback routine using K nearest neighbors prescription
!+
!-----------------------------------------------------------------------
subroutine HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,eos_vars,dt)
 use part,       only:rhoh,massoftype,ihsoft,igas,irateion,isdead_or_accreted,&
                      irstrom,get_partinfo,iphase,itemp,imu
 use linklist,   only:listneigh,getneigh_pos,ifirstincell
 use sortutils,  only:Knnfunc,set_r2func_origin,r2func_origin
 use physcon,    only:pc,pi
 use timing,     only:get_timings,increment_timer,itimer_HII
 use dim,        only:maxvxyzu
 use io,         only:iverbose,iprint,warning
 integer,          intent(in)    :: nptmass,npart
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 real,             intent(inout) :: eos_vars(:,:)
 real,             intent(in)    :: dt
 integer, parameter :: maxc  = 128
 real, save :: xyzcache(maxc,3)
 integer            :: i,k,j,npartin,nneigh,itypej
 real(kind=4)       :: t1,t2,tcpu1,tcpu2
 real               :: pmass,Ndot,DNdot,DNratio,logNdiff,taud,mHII,r,r_in,hcheck
 real               :: xi,yi,zi,log_Qi,stromi,xj,yj,zj,dx,dy,dz,vkx,vky,vkz
 logical            :: momflag,isactive,isgas,isdust,converged

 momflag = .false.
 r = 0.
 r_in = 0.
 pmass = massoftype(igas)

 !if (present(dt)) momflag = .true.
 call get_timings(t1,tcpu1)
 eos_vars(imu,:) = gmw
 if (maxvxyzu < 4) eos_vars(itemp,:) = Tcold  ! cooling in isothermal case

 !
 !-- Rst derivation and thermal feedback
 !
 if (nHIIsources > 0) then
    do i=1,nptmass
       converged = .false.
       npartin=0
       log_Qi = xyzmh_ptmass(irateion,i)
       if (log_Qi <=0.) cycle
       Ndot = log_Qi
       xi = xyzmh_ptmass(1,i)
       yi = xyzmh_ptmass(2,i)
       zi = xyzmh_ptmass(3,i)
       stromi = xyzmh_ptmass(irstrom,i)
       if (stromi >= 0. ) then
          hcheck = stromi + 2.*csion*dt
       else
          hcheck = Rmax
       endif

       call getneigh_pos((/xi,yi,zi/),0.,hcheck,3,listneigh,nneigh,xyzcache,maxc,ifirstincell)
       call Knnfunc(nneigh,(/xi,yi,zi/),xyzh,listneigh)

       if (nneigh > 0) then
          do k=1,nneigh
             j = listneigh(k)
             call get_partinfo(iphase(j),isactive,isgas,isdust,itypej)
             if (isgas) then
                if (.not. isdead_or_accreted(xyzh(4,j))) then
                   ! ionising photons needed to fully ionise the current particle
                   DNdot = DNcoeff + log10(rhoh(xyzh(4,j),pmass))
                   if (DNdot < Ndot) then
                      if (eos_vars(imu,j) > muion ) then ! is ionised ?
                         logNdiff = DNdot - Ndot
                         Ndot = Ndot + log10(1.-10.**(logNdiff))
                         eos_vars(itemp,j)  = Tion
                         eos_vars(imu,j) = muion
                         if (maxvxyzu >= 4) vxyzu(4,j) = uIon
                      endif
                   else
                      if (k > 1) then ! end of the HII region
                         if (maxvxyzu >= 4) then
                            DNratio = Ndot/DNdot
                            eos_vars(itemp,j)  = DNratio * Tion
                            eos_vars(imu,j)    = DNratio * muion
                            vxyzu(4,j)         = DNratio * uIon
                            Ndot = 0.
                         endif
                         r = sqrt((xi-xyzh(1,j))**2 + (yi-xyzh(2,j))**2 + (zi-xyzh(3,j))**2)
                         j = listneigh(1)
                      else ! unresolved case
                         r = 0.
                      endif
                      npartin = k
                      xyzmh_ptmass(irstrom,i) = r
                      converged = .true.
                      exit
                   endif
                endif
             endif
          enddo
       endif
       if (iverbose == 2) write(iprint,*)'Rstrom from sink ',i,' = ',r," with N = ",k,&
                                         ' ionised particle, remaining Nphot : ',Ndot,DNdot
       if (.not.converged) call warning('HII_feedback','Photon march did not converge...',var='Ndot',val=Ndot)
       !
       !-- Momentum feedback
       !
       if (momflag .and. npartin > 3) then
          j = listneigh(1)
          r_in = sqrt((xi-xyzh(1,j))**2 + (yi-xyzh(2,j))**2 + (zi-xyzh(3,j))**2)
          mHII = ((4.*pi*(r**3-r_in**3)*rhoh(xyzh(4,j),pmass))/3)
          if (mHII>3*pmass) then
             !$omp parallel do default(none) &
             !$omp shared(mHII,xyzh,sigd,dt) &
             !$omp shared(mH,vxyzu,log_Qi,hv_on_c,npartin,pmass,xi,yi,zi) &
             !$omp private(j,dx,dy,dz,vkx,vky,vkz,xj,yj,zj,r,taud)
             do k=1,npartin
                j = listneigh(1)
                xj = xyzh(1,j)
                yj = xyzh(2,j)
                zj = xyzh(3,j)
                dx = xj - xi
                dy = yj - yi
                dz = zj - zi
                r = dx**2 + dy**2 + dz**2
                taud = (rhoh(xyzh(4,j),pmass)/mH)*sigd*r
                if (taud > 1.97) taud=1.97
                vkz = (1.+1.5*exp(-taud))*((10**log_Qi)/mHII)*hv_on_c*(dz/r)
                vkx = (1.+1.5*exp(-taud))*((10**log_Qi)/mHII)*hv_on_c*(dx/r)
                vky = (1.+1.5*exp(-taud))*((10**log_Qi)/mHII)*hv_on_c*(dy/r)
                vxyzu(1,j) = vxyzu(1,j) +  vkx*dt
                vxyzu(2,j) = vxyzu(2,j) +  vky*dt
                vxyzu(3,j) = vxyzu(3,j) +  vkz*dt
             enddo
             !$omp end parallel do
          endif
       endif
    enddo
 endif
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_HII,t2-t1,tcpu2-tcpu1)

end subroutine HII_feedback

!-----------------------------------------------------------------------
!+
! HII feedback routine using inversed ray tracing method
!+
!-----------------------------------------------------------------------
subroutine HII_feedback_ray(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,eos_vars)
 use part,     only:massoftype,igas,irateion,isdead_or_accreted,&
                    iphase,get_partinfo,noverlap,rhoh,itemp,imu
 use timing,   only:get_timings,increment_timer,itimer_HII
 use linklist, only:getneigh_pos,ifirstincell,listneigh
 use dim,      only:maxvxyzu
 use io,         only:iverbose,iprint
 integer,          intent(in)    :: nptmass,npart
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 real,             intent(inout) :: eos_vars(:,:)
 real, save         :: xyzcache(maxcache,4)
 real(kind=4)       :: t1,t2,tcpu1,tcpu2
 real,allocatable   :: rhosrc(:)
 real               :: pmass,log_Qi,fluxi,xj,yj,zj
 logical            :: isactive,isgas,isdust
 integer            :: i,j,itypei,noverlapi,nneighi,k
 !$omp threadprivate(xyzcache)

 k=0

 if (nHIIsources > 0) then
    if (iH2R == 3) then
       allocate(rhosrc(nptmass))
       rhosrc = 0.
    endif
    pmass = massoftype(igas)


    call get_timings(t1,tcpu1)
    !
    !-- Rst derivation and thermal feedback
    !
!$omp parallel default(none)&
!$omp shared(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,eos_vars,noverlap)&
!$omp shared(pmass,iphase,ifirstincell,rhosrc,iH2R,Tcold,uIon,gmw)&
!$omp private(log_Qi,i,j,xj,yj,zj,fluxi,itypei,isgas,isdust)&
!$omp private(isactive,noverlapi,nneighi)&
!$omp reduction(+:k)
    if (iH2R == 3) then
!$omp do
       do i=1,nptmass ! if rough approx used, then we must find the density close to the sources
          if (xyzmh_ptmass(irateion,i) <= 0.) cycle
          call getneigh_pos((/xyzmh_ptmass(1,i),xyzmh_ptmass(2,i),xyzmh_ptmass(3,i)/),0.,&
                           0.19,3,listneigh,nneighi,xyzcache,maxcache,ifirstincell)
          if (nneighi > 0.) then
             j = listneigh(1)
             rhosrc(i) = rhoh(xyzh(4,j),pmass)
          endif
       enddo
!$omp enddo
    endif
!$omp do
    do i=1,npart
       if (.not. isdead_or_accreted(xyzh(4,i))) then
          call get_partinfo(iphase(i),isactive,isgas,isdust,itypei)
          if (.not. isactive) cycle
          if (isgas) then
             if (maxvxyzu < 4) eos_vars(itemp,i) = Tcold
             eos_vars(imu,i) = gmw
             noverlapi = 0
             do j=1,nptmass
                log_Qi = xyzmh_ptmass(irateion,j)
                if (log_Qi <= 0.) cycle
                xj = xyzmh_ptmass(1,j)
                yj = xyzmh_ptmass(2,j)
                zj = xyzmh_ptmass(3,j)
                if (iH2R == 3) then
                   call inversed_raytracing(i,(/xj,yj,zj/),xyzh,xyzcache,noverlap,&
                                            pmass,log_Qi,fluxi,rhosrc(j))
                else
                   call inversed_raytracing(i,(/xj,yj,zj/),xyzh,xyzcache,noverlap,&
                                            pmass,log_Qi,fluxi)
                endif
                if (fluxi > 0)then
                   eos_vars(itemp,i) = Tion
                   eos_vars(imu,i) = muion
                   noverlapi = noverlapi + 1
                   if (maxvxyzu >= 4) vxyzu(4,i) = uIon
                   k = k + 1
                endif
             enddo
             noverlap(i) = noverlapi
          endif
       endif
    enddo
!$omp enddo
!$omp end parallel
    call get_timings(t2,tcpu2)
    call increment_timer(itimer_HII,t2-t1,tcpu2-tcpu1)
    if (iverbose == 2) write(iprint,*) 'N ionised = ',k
    if (iH2R == 3) deallocate(rhosrc)
 endif


end subroutine HII_feedback_ray

!-----------------------------------------------------------------------
!+
! find nearest particles along the ray from target to source
!+
!-----------------------------------------------------------------------
subroutine inversed_raytracing(itarg,srcpos,xyzh,xyzcache,noverlap,pmass,log_Q,flux,rhosrcj)
 use part,     only:rhoh,iphase,get_partinfo
 use linklist, only:getneigh_pos,ifirstincell,listneigh
 use kernel,   only:radkern
 use physcon,  only:fourpi
 use units,    only:utime
 integer,        intent(in)    :: itarg
 integer,        intent(in)    :: noverlap(:)
 real,           intent(in)    :: xyzh(:,:)
 real,           intent(inout) :: xyzcache(:,:)

 real,           intent(in)    :: srcpos(3),pmass,log_Q
 real,           intent(out)   :: flux
 real, optional, intent(in)    :: rhosrcj
 real,    parameter :: thetalim = 0.5
 real            :: rayx,rayy,rayz,drisrc1,drisrc,dr,hpmass1,lumS,recvol
 real            :: xij,yij,zij,xi,yi,zi,hi,hj,drproj,dr2toray,icoeff
 real            :: dr2ij,dr2toraymin,nmean,intensity,hi2,drprojmin,hcheck
 integer         :: nneigh,k,i,j,knext,inext,itypei,noverlapj
 logical         :: isactive,isgas,isdust,unreached,rough

 rough = .false.
 if (present(rhosrcj)) rough = .true.

 hpmass1 = 0.5*mH1
 if (rough .and. rhosrcj == 0.) hpmass1 = 2*hpmass1

 lumS    = ((10**log_Q)*utime)/fourpi

 xi      = xyzh(1,itarg)
 yi      = xyzh(2,itarg)
 zi      = xyzh(3,itarg)
 hi      = xyzh(4,itarg)
 rayx    = srcpos(1)-xi
 rayy    = srcpos(2)-yi
 rayz    = srcpos(3)-zi
 drisrc  = sqrt(rayx*rayx + rayy*rayy + rayz*rayz)
 drisrc1 = 1./drisrc
 rayx    = rayx*drisrc1
 rayy    = rayy*drisrc1
 rayz    = rayz*drisrc1
 dr      = drisrc
 hcheck  = hi*radkern
 intensity = 0.
 icoeff = 0.5
 i = itarg
 inext = i
 knext = i
 drprojmin  = huge(dr2ij)
 unreached   = .true.
 if (rough) then
    unreached = .false.
    icoeff = 0.57735026919
 endif

 reachsrc:do while (unreached)

    call getneigh_pos((/xi,yi,zi/),0.,hcheck,3,listneigh,&
                      nneigh,xyzcache,maxcache,ifirstincell)
    hi2 = (hi*radkern)**2
    dr2toraymin = huge(dr2ij)
    drprojmin  = huge(dr2ij)

    !
    !-- identify neighbor closest to the line of sight
    !
    neighsearch:do k=1,nneigh
       j = listneigh(k)
       if (j==i) cycle
       call get_partinfo(iphase(j),isactive,isgas,isdust,itypei)
       if (isgas) then
          if(k <= maxcache) then
             xij = xyzcache(k,1) - xi
             yij = xyzcache(k,2) - yi
             zij = xyzcache(k,3) - zi
          else
             xij = xyzh(1,j) - xi
             yij = xyzh(2,j) - yi
             zij = xyzh(3,j) - zi
          endif
          dr2ij   = xij*xij + yij*yij + zij*zij
          if( dr2ij < hi2) then
             drproj = rayx*xij + rayy*yij + rayz*zij
             if (drproj > 0.) then
                dr2toray = dr2ij - drproj*drproj
                if (dr2toray < dr2toraymin) then
                   dr2toraymin = dr2toray
                   drprojmin = drproj
                   knext = k
                   inext = j
                   if((drprojmin*drprojmin)<acos(thetalim)**2*dr2ij) exit neighsearch
                endif
             endif
          endif
       endif
    enddo neighsearch

    unreached = (dr > drprojmin)
    if (unreached) then
       !
       !-- add projected distance to the distance between sample point and target part
       !
       dr = dr - (0.5*drprojmin)
       if(k<=maxcache) then
          hj  = 1./xyzcache(knext,4)
       else
          hj  = xyzh(4,inext)
       endif
       !
       !-- compute partial recombination volume if overlapped region
       !
       nmean = (rhoh(hi,pmass)+rhoh(hj,pmass))*hpmass1
       noverlapj = noverlap(inext)
       if (noverlapj>0) then
          recvol = (ar*nmean**2)/noverlap(inext)
       else
          recvol = ar*nmean**2
       endif
       !
       !-- integrate intensity
       !
       intensity = intensity + dr**2*drprojmin*recvol

       if (lumS < intensity) exit reachsrc! exit if already neutral
       !
       !-- prepare the next search
       !
       if(knext <= maxcache) then
          xi = xyzcache(knext,1)
          yi = xyzcache(knext,2)
          zi = xyzcache(knext,3)
       else
          xi = xyzh(1,inext)
          yi = xyzh(2,inext)
          zi = xyzh(3,inext)
       endif
       hi = hj
       hcheck = hi*radkern
       i  = inext
       dr = dr - (0.5*drprojmin)
    endif
 enddo reachsrc
 if (lumS > intensity) then
    !
    !-- Last iteration before reaching the source
    !
    drprojmin  = dr
    if (rough) then
       nmean    = (rhoh(hi,pmass)+rhosrcj)*(hpmass1)
    else
       nmean    = rhoh(hi,pmass)*2*hpmass1
    endif
    noverlapj = noverlap(inext)
    if (noverlapj>0) then
       recvol = (ar*(nmean)**2)/noverlap(inext)
    else
       recvol = ar*(nmean)**2
    endif
    intensity = intensity + (icoeff*dr)**2*drprojmin*recvol
    flux = (lumS-intensity)
 else
    flux = -1.
 endif


end subroutine inversed_raytracing


!-----------------------------------------------------------------------
!+
!  write options to input file
!+
!-----------------------------------------------------------------------
subroutine write_options_H2R(iunit)
 use infile_utils, only:write_inopt
 use physcon,      only:solarm
 integer, intent(in) :: iunit

 if (iH2R > 0) then
    write(iunit,"(/,a)") '# options controlling HII region expansion feedback'
    call write_inopt(iH2R, 'iH2R',&
    "HII region expansion feedback in star forming region(1 Knn method,2 inv RT method)", iunit)
    call write_inopt(Mmin, 'Mmin', "Minimum star mass to trigger HII region (MSun)", iunit)
    call write_inopt(Rmax, 'Rmax', "Maximum radius for HII region (pc)(only Knn)", iunit)
 endif

end subroutine write_options_H2R

!-----------------------------------------------------------------------
!+
!  read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_H2R(name,valstring,imatch,igotall,ierr)
 use io,         only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter  :: label = 'read_options_H2R'

 imatch = .true.
 select case(trim(name))
 case('iH2R')
    read(valstring,*,iostat=ierr) iH2R
    if (iH2R < 0) call fatal(label,'HII region option out of range')
    ngot = ngot + 1
 case('Mmin')
    read(valstring,*,iostat=ierr) Mmin
    if (Mmin < 8.) call fatal(label,'Minimimum mass can not be inferior to 8 solar masses')
    ngot = ngot + 1
 case('Rmax')
    read(valstring,*,iostat=ierr) Rmax
    if (Rmax < 10.) call fatal(label,'Maximum radius can not be inferior to 10 pc')
    ngot = ngot + 1
 case default
    imatch = .true.
 end select
 igotall = (ngot >= 3)

end subroutine read_options_H2R

end module HIIRegion
