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
 ! routine originally made by Fujii et al 2021
 ! adapted in Phantom by Yann BERNARD
 ! reference : Fujii et al. 2021 SIRIUS Project Paper III
 !
 !
 use part, only:nbpart,npart

 implicit none

 public :: update_fbsource, update_Q_list, HII_feedback,initialize_fb,search_connected_HII,check_ionized_sinks
 public :: allocate_fb,deallocate_fb

 integer, public               :: nbfbmax = 300
 integer, public               :: nbfbs = 0
 integer, public               :: iFb = 0
 real,    private, parameter   :: a = -39.3178
 real,    private, parameter   :: b =  221.997
 real,    private, parameter   :: c = -227.456
 real,    private, parameter   :: d =  117.410
 real,    private, parameter   :: e = -30.1511
 real,    private, parameter   :: f =  3.06810
 real,    private, parameter   :: ar_cgs = 2.6d-13
 real,    private, parameter   :: sigd_cgs = 1.d-21
 real,    private              :: ar
 real,    private              :: sigd
 real,    private              :: hv_on_c
 real,    private              :: mu = 2.38
 real,    private              :: mH
 real,    private              :: T_ion
 real,    private              :: u_to_t
 real,    private              :: Rst2_max
 logical, private              :: overlapping =.false.
 integer,          allocatable :: source_id(:)
 real,             allocatable :: Qsource(:)
 real,             allocatable :: dxyz (:,:,:)
 real,             allocatable :: r2(:)
 real,             allocatable :: overlap_e(:)
 real,             allocatable :: Rst_source(:)
 integer,          allocatable :: arg_r2(:,:)
 logical, public,  allocatable :: isionised(:)
 private

contains

 !-----------------------------------------------------------------------
 !+
 !  Initialise stellar feedbacks
 !+
 !-----------------------------------------------------------------------
subroutine initialize_fb
 use units,   only:udist,umass,utime
 use physcon, only:mass_proton_cgs,kboltz,atomic_mass_unit,pc,eV
 use eos    , only:gamma,gmw
 call allocate_fb(nbfbmax)
 isionised(:)=.false.
 source_id(:)= 0
 !calculate the useful constant in code units
 mH = gmw*mass_proton_cgs
 u_to_t = (3./2)*(kboltz/mH)*(utime/udist)**2
 mH = mH/umass
 T_ion = 1.d4
 ar = ar_cgs*utime/udist**3
 sigd = sigd_cgs*udist**2
 hv_on_c = ((18.6*eV)/2.997924d10)*(utime/(udist*umass))
 Rst2_max = ((15*pc)/udist)**2
 print*,"feedback constants mH,u_to_t,T_ion,u_to_t*T_ion,gmw  : ",mH,u_to_t,T_ion,u_to_t*T_ion,gmw
 !open(20,file="Rst.dat")
 return
end subroutine initialize_fb

 !-----------------------------------------------------------------------
 !+
 !  subroutine that gives the number of sources
 !+
 !-----------------------------------------------------------------------

subroutine update_fbsource(pmass,i)
 use part, only: nbpart
 real,    intent(in) :: pmass
 integer, intent(in) :: i
 ! select feedback source with a minimum mass of 8 Msun
 if(pmass>8) then
    nbfbs = nbfbs + 1
    source_id(nbfbs) = i
    call update_Q_list(pmass)
 endif
 return
end subroutine update_fbsource

 !-----------------------------------------------------------------------
 !+
 !  Calculation of the the ionizing photon rate
 !+
 !-----------------------------------------------------------------------

subroutine update_Q_list(pmass)
 use units, only:utime
 use part,  only:xyzmh_bpart
 real, intent(in) :: pmass
 real    :: log_pmassj,log_Q
 ! caluclation of the ionizing photon rate  of each sources
 ! this calculation uses Fujii's formula derived from OSTAR2002 databases
 log_pmassj = log10(pmass)
 log_Q = (a+b*log_pmassj+c*log_pmassj**2+d*log_pmassj**3+e*log_pmassj**4+f*log_pmassj**5)
 Qsource(nbfbs) = (10.**log_Q)*utime
 print*,"New source detected : Log Q : ",log_Q
 print*,"nb_feedback sources : ",nbfbs
 return
end subroutine update_Q_list

 !-----------------------------------------------------------------------
 !+
 !  Main subroutine : Application of the HII feedback using Hopkins's like prescription
 !+
 !-----------------------------------------------------------------------

subroutine HII_feedback(dt)
 use part,   only:xyzh,xyzmh_bpart,vxyzu,rhoh,massoftype
 use units,  only:unit_density,udist,umass
 use physcon,only:pc,pi
 use utils_stellarfb,only:merge_argsort,print_fblog_time
 use timing, only: get_timings
 real(kind=4)    :: t1,t2,tcpu1,tcpu2
 real, intent(in) :: dt
 integer :: i,j,l,k,n
 real    :: pmass,Ndot,DNdot,R_stop,eps,taud_on_r,taud,mHII,v_kick,r
 pmass = massoftype(1)
 ! at each new kick we reset all the particles status
 isionised(:) = .false.
 Rst_source(:) = 0.
 overlap_e(:) = 0.
 eps = xyzmh_bpart(5,1)
 !
 !!!!!!! Rst derivation and thermal feedback
 !
 call get_timings(t1,tcpu1)
 do i=1,nbfbs
    n=size(r2)
    j=source_id(i)
    ! for each source we compute the distances of each particles and sort to have a Knn list
    ! Patch : We need to be aware of dead particles that will pollute the scheme if not taking into account.
    ! The simpliest way is to put enormous distance for dead particle to be at the very end of the knn list.
    do l=1,npart
       if (xyzh(4,l)<0.) then
          dxyz(:,l,i) = huge(pmass)
       else
          dxyz(1,l,i) = xyzh(1,l)-xyzmh_bpart(1,j)
          dxyz(2,l,i) = xyzh(2,l)-xyzmh_bpart(2,j)
          dxyz(3,l,i) = xyzh(3,l)-xyzmh_bpart(3,j)
       endif
    enddo
    r2(:) = dxyz(1,:,i)**2+dxyz(2,:,i)**2+dxyz(3,:,i)**2
    call merge_argsort(r2,arg_r2(:,i))
    k = arg_r2(n,i)
    ! calculation of the ionised mass
    Ndot = Qsource(i)
    DNdot = (pmass*ar*rhoh(xyzh(4,k),pmass))/(mH**2)
    !print*,"Ndot : DNdot : local rho : ",Ndot,DNdot,rhoh(xyzh(4,k),pmass)*unit_density
    if (Ndot>DNdot) then
       ! iteration on the Knn until we used all the source photons
       if (r2(k)<Rst2_max) then
          do while (Ndot>DNdot .and. n/=0 .and. r2(k)<Rst2_max)
             if (.not.(isionised(k))) then
                Ndot = Ndot - DNdot
                vxyzu(4,k) = u_to_t*T_ion
                isionised(k)=.true.
             else !If already ionised : the ionization energy is stored for regularization
                Ndot = Ndot - DNdot
                overlap_e(i) = overlap_e(i) + DNdot
                overlapping = .true.
             endif
             n = n-1
             k = arg_r2(n,i)
             DNdot = (pmass*ar*rhoh(xyzh(4,k),pmass))/(mH**2)
          enddo
          Rst_source(i) = sqrt(r2(k))
       else
          Rst_source(i) = sqrt(Rst2_max)
       endif
       !write(20,*) Rst_source(i),Rst
       !print*,"remaining ionization Energy: ",(Ndot)/Qsource(i),(Ndot/DNdot)*T_ion
       !print*,"ionised particles, Rst  : ",npart-n,Rst_source(i)
    else
       ! unresolved cased...
       !print*,"source unresolved"
       vxyzu(4,k) = (Ndot/DNdot)*u_to_t*T_ion
       Rst_source(i) = 0.
    endif
 enddo
 !
 !!!!!!!! overlap regulartization
 !
 if(overlapping) then
    !print*,"overlapping detection : regularization of the ionization front"
    call search_connected_HII(nbfbs)
    do i=1,nbfbs
       n=size(r2)
       k = arg_r2(n,i)
       Ndot = overlap_e(i)
       DNdot = (pmass*ar*rhoh(xyzh(4,k),pmass))/(mH**2)
       ! iteration on the Knn until we used all the source photons
       do while (Ndot>DNdot .and. n/=0 .and. r2(k)<Rst2_max)
          if (.not.(isionised(k))) then
             Ndot = Ndot - DNdot
             vxyzu(4,k) = u_to_t*T_ion
             isionised(k)=.true.
          endif
          n = n-1
          k = arg_r2(n,i)
          DNdot = (pmass*ar*rhoh(xyzh(4,k),pmass))/(mH**2)
       enddo
    enddo
 endif
 !
 !!!!!!!! momentum feedback
 !
 !print*, "Adding momentum feedback in HII region"
 do i=1, nbfbs
    if (Rst_source(i)< 2*eps) then
       R_stop = 2*eps
    else
       R_stop = Rst_source(i)
    endif
    n=size(r2)
    k = arg_r2(n,i)
    r = sqrt(dxyz(1,k,i)**2 + dxyz(2,k,i)**2 + dxyz(3,k,i)**2)
    taud_on_r = (rhoh(xyzh(4,k),pmass)/mH)*sigd
    mHII = ((4.*pi*(R_stop**3-r**3)*rhoh(xyzh(4,k),pmass))/3)
    !print*,"MHII and Momentum prefactor :",mHII,(Qsource(i)/mHII)*hv_on_c,R_stop,Rst_source(i)
    if (mHII>3*pmass) then
       do while (r<R_stop)
          do j=1,3
             taud = taud_on_r*abs(dxyz(j,k,i)/r)
             if (taud > 1.97) taud=1.97
             v_kick = (1.+1.5*exp(-taud))*(Qsource(i)/mHII)*hv_on_c*(dxyz(j,k,i)/r)
             vxyzu(j,k) = vxyzu(j,k) +  v_kick*dt
          enddo
          n=n-1
          k = arg_r2(n,i)
          r = sqrt(dxyz(1,k,i)**2 + dxyz(2,k,i)**2 + dxyz(3,k,i)**2)
       enddo
       !print*, "real MHII", (size(r2)-n)*pmass
    endif
 enddo



 ! resetting overlap flag for the next step.
 overlapping = .false.
 call get_timings(t2,tcpu2)
 call print_fblog_time(nbfbs,tcpu2-tcpu1)
 return
end subroutine HII_feedback


subroutine search_connected_HII(nb)
 use part,   only: xyzmh_bpart
 use utils_stellarfb, only:jacobi_eigenvalue
 integer, intent(in) :: nb
 real                :: LapMatrix(nb,nb)
 real                :: EigenVec(nb,nb)
 real                :: EigenV(nb)
 integer            :: i,j,k,l,nb_region,nb_node
 real                :: dist,dx,dy,dz
 real                :: region_ov_e
 ! construct laplacian matrix to identify unconnected components of the graph
 LapMatrix = 0.
 do i=1,nbfbs
    do j=1,nbfbs
       k = source_id(i)
       l = source_id(j)
       dx =(xyzmh_bpart(1,k)-xyzmh_bpart(1,l))
       dy =(xyzmh_bpart(2,k)-xyzmh_bpart(2,l))
       dz =(xyzmh_bpart(3,k)-xyzmh_bpart(3,l))
       dist = sqrt(dx**2+dy**2+dz**2)
       !print*,dist,Rst_source(i)+Rst_source(j)
       ! ici il faut résoudre le soucis de la connexion sur des sources non résolues ou coupées car trop loin
       if (Rst_source(i)/=0.0 .and. Rst_source(j)/=0.0) then
          if (dist< Rst_source(i)+Rst_source(j) .and. i/=j) then
             LapMatrix(i,j)= - 1
          endif
       endif
    enddo
    LapMatrix(i,i) = abs(sum(LapMatrix(i,:)))
    !print*,LapMatrix(i,:)
 enddo
 ! compute egeinvalues and vectors of the Laplacian matrix
 call jacobi_eigenvalue(nb,LapMatrix,1000,EigenVec,EigenV)
 nb_region = count(EigenV<0.000001)
 !print*,EigenV
 do i=1, nb_region
    !print*,"region : ",i
    region_ov_e=0
    nb_node=0
    do j=1, nbfbs
       if (EigenVec(j,i)/=0.)then
          region_ov_e = region_ov_e + overlap_e(j)
          nb_node = nb_node + 1
          !print*,"member : ",j
       endif
    enddo
    do j=1, nbfbs
       if (EigenVec(j,i)/=0.)then
          overlap_e(j) = region_ov_e/nb_node
       endif
    enddo
 enddo

end subroutine search_connected_HII

 !-----------------------------------------------------------------------
 !+
 !  The aim of this subroutine is to warned if a sink has been ionized by
 !  a massive star. One star can only ionized one sink, beacause it can be only in one sink.
 !  Thus these check need to collect m_acc t_mean to compute the ionizing criterion.
 !+
 !-----------------------------------------------------------------------

subroutine check_ionized_sinks(msflag,merged_ptmass)
 use part,    only: nptmass,xyzmh_ptmass,xyzmh_bpart,ihacc,icmpast,itcreate,t_acc
 use physcon, only: pi
 logical, intent(out) :: msflag(nptmass)
 integer, intent(in) :: merged_ptmass(:)
 integer :: i,j,k,kmax,nmerged
 real :: rsquare,h2,tmean,macc,rhos,rst
 msflag(:)=.false.
 ! Check if any MS are in the vinicity of a sink
 do i = 1, nbfbs
    if (xyzmh_bpart(4,i)<15.)cycle
    kmax  = 0
    tmean = 0.
    macc  = 0.
    do j = 1, nptmass
       if (xyzmh_ptmass(ihacc,j)<0.) cycle
       if (xyzmh_ptmass(4,j)<0.)cycle
       if (msflag(j)) cycle
       h2 = xyzmh_ptmass(ihacc,i)**2
       rsquare = (xyzmh_ptmass(1,j)-xyzmh_bpart(1,i))**2+(xyzmh_ptmass(2,j)-xyzmh_bpart(2,i))**2&
                +(xyzmh_ptmass(3,j)-xyzmh_bpart(3,i))**2
       if (rsquare<h2) then !check passed only for one sink posssibe
          kmax = j
          exit
       endif
    enddo
    if (kmax/=0) then
       do k=1,nptmass
          if (merged_ptmass(k)==kmax) then
             macc = macc + xyzmh_ptmass(icmpast,k)
             tmean = tmean + xyzmh_ptmass(itcreate,k)
             nmerged = nmerged + 1
          endif
       enddo
       tmean = tmean/(nmerged*t_acc)
       macc = xyzmh_ptmass(4,kmax)-macc
       rhos = (tmean*macc)/(4*pi*(xyzmh_ptmass(ihacc,kmax)**3)/3)
       rst = ((3*Qsource(i))/(4*pi*(rhos/mH)**2))**(1./3)
       if (rst>0.5*(xyzmh_ptmass(ihacc,kmax))) then
          msflag(kmax) = .true.
       endif
    endif
 enddo
end subroutine check_ionized_sinks

subroutine allocate_fb(n)
 use allocutils, only:allocate_array
 integer, intent(in) :: n
 call allocate_array("Qsource"  , Qsource  , n)
 call allocate_array("Rst_source"  , Rst_source  , n)
 call allocate_array("overlap_e", overlap_e, n)
 call allocate_array("dxyz", dxyz, 3,npart, n)
 call allocate_array("r2", r2, npart)
 call allocate_array("arg_r2", arg_r2, npart,n)
 call allocate_array("isionised", isionised, npart)
 call allocate_array('source_id', source_id, n)
end subroutine allocate_fb

subroutine deallocate_fb
 deallocate(source_id)
 deallocate(Qsource)
 deallocate(overlap_e)
 deallocate(Rst_source)
 deallocate(dxyz)
 deallocate(r2)
 deallocate(arg_r2)
 deallocate(isionised)
end subroutine deallocate_fb


end module HIIRegion
