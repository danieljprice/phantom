!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module HIIRegion
!
! HIIRegion
! contains routines to model HII region expansion due to ionization and radiation pressure..
! routine originally made by Hopkins et al. (2012),reused by Fujii et al. (2021)
! and adapted in Phantom by Yann Bernard
!
! :References: Fujii et al. (2021), Hopkins et al. (2012)
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, infile_utils, io, linklist, part, physcon,
!   sortutils, timing, units
!
 implicit none

 public :: update_ionrates,update_ionrate, HII_feedback,initialize_H2R,read_options_H2R,write_options_H2R

 integer, parameter, public    :: HIIuprate   = 8 ! update rate when IND_TIMESTEPS=yes
 integer, public               :: iH2R = 0
 real   , public               :: Rmax = 15 ! Maximum HII region radius (pc) to avoid artificial expansion...
 real   , public               :: Mmin = 8  ! Minimum mass (Msun) to produce HII region
 integer, public               :: nHIIsources = 0
 real   , public               :: ar
 real   , public               :: mH

 real, parameter   :: a = -39.3178 !
 real, parameter   :: b =  221.997 !  fitted parameters to compute
 real, parameter   :: c = -227.456 !  ionisation rate for massive
 real, parameter   :: d =  117.410 !  extracted from Fujii et al. (2021).
 real, parameter   :: e = -30.1511 ! (Expressed in function of log(solar masses) and s)
 real, parameter   :: f =  3.06810 !
 real, parameter   :: ar_cgs = 2.7d-13
 real, parameter   :: sigd_cgs = 1.d-21
 real              :: sigd
 real              :: hv_on_c
 real              :: Tion
 real              :: Rst_max
 real              :: Minmass
 real              :: uIon

 private

contains

!-----------------------------------------------------------------------
!+
!  Initialise stellar feedbacks
!+
!-----------------------------------------------------------------------
subroutine initialize_H2R
 use io,      only:iprint,iverbose,id,master
 use part,    only:isionised
 use units,   only:udist,umass,utime
 use physcon, only:mass_proton_cgs,kboltz,pc,eV,solarm
 use eos,     only:gmw,gamma

 isionised(:)=.false.
 !calculate the useful constant in code units
 mH = gmw*mass_proton_cgs
 Tion = 1.e4
 ar = ar_cgs*(utime/udist**3)
 sigd = sigd_cgs*udist**2
 hv_on_c = ((18.6*eV)/2.997924d10)*(utime/(udist*umass))
 Rst_max = sqrt(((Rmax*pc)/udist)**2)
 Minmass = (Mmin*solarm)/umass
 if (gamma>1.) then
    uIon = kboltz*Tion/(mH*(gamma-1.))*(utime/udist)**2
 else
    uIon = 1.5*(kboltz*Tion/(mH))*(utime/udist)**2
 endif

 mH = mH/umass

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
!  Main subroutine : Application of the HII feedback using Hopkins's like prescription
!+
!-----------------------------------------------------------------------
subroutine HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,isionised,dt)
 use part,       only:rhoh,massoftype,ihsoft,igas,irateion,isdead_or_accreted,&
                      irstrom
 use linklist,   only:listneigh=>listneigh_global,getneigh_pos,ifirstincell
 use sortutils,  only:Knnfunc,set_r2func_origin,r2func_origin
 use physcon,    only:pc,pi
 use timing,     only:get_timings,increment_timer,itimer_HII
 use dim,        only:maxvxyzu,maxpsph
 use units,      only:utime
 integer,          intent(in)    :: nptmass,npart
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 logical,          intent(inout) :: isionised(:)
 real,   optional, intent(in)    :: dt
 integer, parameter :: maxcache      = 12000
 real, save :: xyzcache(maxcache,3)
 integer            :: i,k,j,npartin,nneigh
 real(kind=4)       :: t1,t2,tcpu1,tcpu2
 real               :: pmass,Ndot,DNdot,logNdiff,taud,mHII,r,r_in,hcheck
 real               :: xi,yi,zi,log_Qi,stromi,xj,yj,zj,dx,dy,dz,vkx,vky,vkz
 logical            :: momflag

 momflag = .false.
 r = 0.
 r_in = 0.

 if (present(dt)) momflag = .true.

 ! at each new kick we reset all the particles status
 isionised(:) = .false.
 pmass = massoftype(igas)

 call get_timings(t1,tcpu1)
 !
 !-- Rst derivation and thermal feedback
 !
 if (nHIIsources > 0) then
    do i=1,nptmass
       npartin=0
       log_Qi = xyzmh_ptmass(irateion,i)
       if (log_Qi <=0.) cycle
       Ndot = log_Qi ! instead of working with very large number, we'll work in logspace now
       xi = xyzmh_ptmass(1,i)
       yi = xyzmh_ptmass(2,i)
       zi = xyzmh_ptmass(3,i)
       stromi = xyzmh_ptmass(irstrom,i)
       if (stromi >= 0. ) then
          hcheck = 1.4*stromi + 0.01*Rmax
       else
          hcheck = Rmax
       endif
       do while(hcheck <= Rmax)
          call getneigh_pos((/xi,yi,zi/),0.,hcheck,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
          call set_r2func_origin(xi,yi,zi)
          call Knnfunc(nneigh,r2func_origin,xyzh,listneigh) !! Here still serial version of the quicksort. Parallel version in prep..
          if (nneigh > 0) exit
          hcheck = hcheck + 0.01*Rmax  ! additive term to allow unresolved case to open
       enddo
       do k=1,nneigh
          j = listneigh(k)
          if (j > maxpsph) cycle
          if (.not. isdead_or_accreted(xyzh(4,j))) then
             ! ionising photons needed to fully ionise the current particle
             DNdot = log10((((pmass*ar*rhoh(xyzh(4,j),pmass))/(mH**2))/utime))
             if (Ndot>DNdot) then
                if (.not.(isionised(j))) then
                   logNdiff = DNdot -Ndot
                   Ndot = Ndot + log10(1-10**(logNdiff))
                   isionised(j)=.true.
                   if (maxvxyzu >= 4) vxyzu(4,j) = uIon
                endif
             else
                if (k > 1) then
                   ! end of the HII region
                   r = sqrt((xi-xyzh(1,j))**2 + (yi-xyzh(2,j))**2 + (zi-xyzh(3,j))**2)
                   j = listneigh(1)
                else
                   ! unresolved case
                   r = 0.
                endif
                exit
             endif
          endif
       enddo
       npartin = k
       xyzmh_ptmass(irstrom,i) = r
       !
       !-- Momentum feedback
       !
       if (momflag .and. npartin > 3) then
          j = listneigh(1)
          r_in = sqrt((xi-xyzh(1,j))**2 + (yi-xyzh(2,j))**2 + (zi-xyzh(3,j))**2)
          mHII = ((4.*pi*(r**3-r_in**3)*rhoh(xyzh(4,j),pmass))/3)
          if (mHII>3*pmass) then
             !$omp parallel do default(none) &
             !$omp shared(mHII,listneigh,xyzh,sigd,dt) &
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
!  write options to input file
!+
!-----------------------------------------------------------------------
subroutine write_options_H2R(iunit)
 use infile_utils, only:write_inopt
 use physcon,      only:solarm
 integer, intent(in) :: iunit

 if (iH2R > 0) then
    write(iunit,"(/,a)") '# options controlling HII region expansion feedback'
    call write_inopt(iH2R, 'iH2R', "enable the HII region expansion feedback in star forming reigon", iunit)
    call write_inopt(Mmin, 'Mmin', "Minimum star mass to trigger HII region (MSun)", iunit)
    call write_inopt(Rmax, 'Rmax', "Maximum radius for HII region (pc)", iunit)
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
