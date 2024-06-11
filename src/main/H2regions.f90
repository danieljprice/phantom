!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
! Module proposed By Yann BERNARD to implement stellar feedbacks in cluster!
! simulations                                                              !
!--------------------------------------------------------------------------!

module HIIRegion
 !
 !
 ! contains routine for Stromgren radius calculation and Radiative pressure velocity kick
 ! routine originally made by Hopkins et al. (2012) Fujii et al. (2021)
 ! adapted in Phantom by Yann BERNARD
 ! reference : Fujii et al. 2021 SIRIUS Project Paper III
 !
 !

 implicit none

 public :: update_ionrates,update_ionrate, HII_feedback,initialize_H2R,read_options_H2R,write_options_H2R

 integer, public               :: iH2R = 0
 real   , public               :: Rmax = 15 ! Maximum HII region radius (pc) to avoid artificial expansion...
 real   , public               :: Mmin = 8  ! Minimum mass (Msun) to produce HII region
 real   , public               :: nHIIsources = 0

 real,    private, parameter   :: a = -39.3178 !
 real,    private, parameter   :: b =  221.997 !  fitted parameters to compute
 real,    private, parameter   :: c = -227.456 !  ionisation rate for massive
 real,    private, parameter   :: d =  117.410 !  extracted from Fujii et al. (2021).
 real,    private, parameter   :: e = -30.1511 ! (Expressed in function of log(solar masses) and s)
 real,    private, parameter   :: f =  3.06810 !
 real,    private, parameter   :: ar_cgs = 2.7d-13
 real,    private, parameter   :: sigd_cgs = 1.d-21
 real,    private              :: ar
 real,    private              :: sigd
 real,    private              :: hv_on_c
 real,    private              :: mH
 real,    private              :: T_ion
 real,    private              :: u_to_t
 real,    private              :: Rst2_max
 real,    private              :: Rst_max

 private

contains

 !-----------------------------------------------------------------------
 !+
 !  Initialise stellar feedbacks
 !+
 !-----------------------------------------------------------------------
subroutine initialize_H2R
 use io,      only:iprint,iverbose
 use part,    only:isionised
 use units,   only:udist,umass,utime
 use physcon, only:mass_proton_cgs,kboltz,pc,eV,solarm
 use eos    , only:gmw
 isionised(:)=.false.
 !calculate the useful constant in code units
 mH = gmw*mass_proton_cgs
 u_to_t = (3./2)*(kboltz/mH)*(utime/udist)**2
 mH = mH/umass
 T_ion = 1.d4
 ar = ar_cgs*(utime/udist**3)
 sigd = sigd_cgs*udist**2
 hv_on_c = ((18.6*eV)/2.997924d10)*(utime/(udist*umass))
 Rst2_max = ((Rmax*pc)/udist)**2
 Rst_max = sqrt(Rst_max)
 Mmin = (Mmin*solarm)/umass
 if (iverbose > 1) then
    write(iprint,"(/a,es18.10,es18.10/)") "feedback constants mH, u_to_t : ", mH, u_to_t
 endif
 return
end subroutine initialize_H2R

!-----------------------------------------------------------------------
!+
!  Calculation of the the ionizing photon rate of all stars (Only for restart)
!+
!-----------------------------------------------------------------------

subroutine update_ionrates(nptmass,xyzmh_ptmass,h_acc)
 use io,     only:iprint,iverbose
 use units,  only:utime
 use part,   only:irateion,ihacc
 integer, intent(in)    :: nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 real,    intent(in)    :: h_acc
 real    :: logmi,log_Q,mi,hi,Q
 integer :: i
 nHIIsources = 0
 !$omp parallel do default(none) &
 !$omp shared(xyzmh_ptmass,nptmass,iprint,iverbose,utime,Mmin,h_acc)&
 !$omp private(logmi,log_Q,Q,mi,hi)&
 !$omp reduction(+:nHIIsources)
 do i=1,nptmass
    mi = xyzmh_ptmass(4,i)
    hi = xyzmh_ptmass(ihacc,i)
    if(mi > Mmin .and. hi < h_acc)then
       logmi = log10(mi)
       ! caluclation of the ionizing photon rate  of each sources
       ! this calculation uses Fujii's formula derived from OSTAR2002 databases
       log_Q = (a+b*logmi+c*logmi**2+d*logmi**3+e*logmi**4+f*logmi**5)
       Q = (10.**log_Q)*utime
       xyzmh_ptmass(irateion,i) = Q
       nHIIsources = nHIIsources + 1
       if (iverbose > 0) then
          write(iprint,"(/a,es18.10/)")"Massive stars detected : Log Q : ",log_Q
       endif
    else
       xyzmh_ptmass(irateion,i) = -1.
    endif
 enddo
 !$omp end parallel do
 if (iverbose > 1) then
    write(iprint,"(/a,i8/)") "nb_feedback sources : ",nHIIsources
 endif
 return
end subroutine update_ionrates

subroutine update_ionrate(i,xyzmh_ptmass,h_acc)
 use io,     only:iprint,iverbose
 use units,  only:utime
 use part,   only:irateion,ihacc
 integer, intent(in)    :: i
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 real,    intent(in)    :: h_acc
 real    :: logmi,log_Q,mi,hi,Q
 mi = xyzmh_ptmass(4,i)
 hi = xyzmh_ptmass(ihacc,i)
 if(mi > Mmin .and. hi < h_acc)then
    logmi = log10(mi)
    ! caluclation of the ionizing photon rate  of each sources
    ! this calculation uses Fujii's formula derived from OSTAR2002 databases
    log_Q = (a+b*logmi+c*logmi**2+d*logmi**3+e*logmi**4+f*logmi**5)
    Q = (10.**log_Q)*utime
    xyzmh_ptmass(irateion,i) = Q
    nHIIsources = nHIIsources + 1
    if (iverbose > 0) then
       write(iprint,"(/a,es18.10/)")"(HII region) Massive stars detected : Log Q : ",log_Q
    endif
 else
    xyzmh_ptmass(irateion,i) = -1.
 endif

 if (iverbose > 1) then
    write(iprint,"(/a,i8/)") "nb_feedback sources : ",nHIIsources
 endif
 return
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
 use timing,     only: get_timings
 integer,          intent(in)    :: nptmass,npart
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 logical,          intent(inout) :: isionised(:)
 real,   optional, intent(in)    :: dt
 integer, parameter :: maxcache      = 12000
 real, save :: xyzcache(maxcache,3)
 integer            :: i,k,j,npartin,nneigh
 real(kind=4)       :: t1,t2,tcpu1,tcpu2
 real               :: pmass,Ndot,DNdot,taud,mHII,r,r_in,hcheck
 real               :: xi,yi,zi,Qi,stromi,xj,yj,zj,dx,dy,dz,vkx,vky,vkz
 logical            :: momflag

 momflag = .false.
 r = 0.
 r_in = 0.

 if (present(dt)) momflag = .true.

 ! at each new kick we reset all the particles status
 isionised(:) = .false.
 pmass = massoftype(igas)
 !
 !-- Rst derivation and thermal feedback
 !
 if (nHIIsources > 0) then
    do i=1,nptmass
       npartin=0
       Qi = xyzmh_ptmass(irateion,i)
       if (Qi <=0.) cycle
       Ndot = Qi
       xi = xyzmh_ptmass(1,i)
       yi = xyzmh_ptmass(2,i)
       zi = xyzmh_ptmass(3,i)
       stromi = xyzmh_ptmass(irstrom,i)
       ! for each source we compute the distances of each particles and sort to have a Knn list
       ! Patch : We need to be aware of dead particles that will pollute the scheme if not taking into account.
       ! The simpliest way is to put enormous distance for dead particle to be at the very end of the knn list.
       if(stromi > 0 ) then
          hcheck = 2.*stromi
          if (hcheck > Rmax) hcheck = Rmax
       else
          hcheck = Rmax
       endif
       call get_timings(t1,tcpu1)
       call getneigh_pos((/xi,yi,zi/),0.,hcheck,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
       call set_r2func_origin(xi,yi,zi)
       call Knnfunc(nneigh,r2func_origin,xyzh,listneigh)
       do k=1,npart
          j = listneigh(k)
          if (.not. isdead_or_accreted(xyzh(4,j))) then
             ! calculation of the ionised mass
             DNdot = (pmass*ar*rhoh(xyzh(4,j),pmass))/(mH**2)
             if (Ndot>DNdot) then
                ! iteration on the Knn until we used all the source photons
                if (.not.(isionised(j))) then
                   Ndot = Ndot - DNdot
                   isionised(j)=.true.
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
       if(momflag .and. npartin > 3) then
          j = listneigh(1)
          r_in = sqrt((xi-xyzh(1,j))**2 + (yi-xyzh(2,j))**2 + (zi-xyzh(3,j))**2)
          mHII = ((4.*pi*(r**3-r_in**3)*rhoh(xyzh(4,j),pmass))/3)
          if (mHII>3*pmass) then
!$omp parallel do default(none) &
!$omp shared(mHII,listneigh,xyzh,sigd,dt) &
!$omp shared(mH,vxyzu,Qi,hv_on_c,npartin,pmass,xi,yi,zi) &
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
                vkx = (1.+1.5*exp(-taud))*(QI/mHII)*hv_on_c*(dx/r)
                vky = (1.+1.5*exp(-taud))*(QI/mHII)*hv_on_c*(dy/r)
                vkz = (1.+1.5*exp(-taud))*(QI/mHII)*hv_on_c*(dz/r)
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
 !print*, "HII feedback CPU time : ",t2-t1
 return
end subroutine HII_feedback

subroutine write_options_H2R(iunit)
 use infile_utils, only:write_inopt
 use physcon,      only:solarm
 use units,        only:umass
 integer, intent(in) :: iunit
 write(iunit,"(/,a)") '# options controlling HII region expansion feedback'
 if(iH2R>0) then
    call write_inopt(iH2R, 'iH2R', "enable the HII region expansion feedback in star forming reigon", iunit)
    call write_inopt((Mmin*umass)/solarm, 'Mmin', "Minimum star mass to trigger HII region (MSun)", iunit)
    call write_inopt(Rmax, 'Rmax', "Maximum radius for HII region (pc)", iunit)
 endif
end subroutine write_options_H2R

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
