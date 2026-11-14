!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! :Runtime parameters: None
!
! :Dependencies: dim, eos, eos_HIIR, infile_utils, io, kernel, neighkdtree,
!   part, physcon, sortutils, units
!
 use eos_HIIR,   only:csion,uIon,Tion,Tcold,muion
 use eos,        only:gmw
 implicit none

 public :: update_ionrates,update_ionrate, HII_feedback,&
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
 logical, public               :: HIIupdateflag
 real   , public               :: HIIdt

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
 ! photons needed to ionise a particle are this times its mass and density
 DNcoeff = log10(ar*mH1**2/utime)
 HIIupdateflag = .true.
 HIIdt = 0.

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
 if (nptmass > 0) then
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
! HII feedback interface for switching algorithms
!+
!-----------------------------------------------------------------------
subroutine HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyz_ptmass,vxyzu,rho,eos_vars)
 integer, intent(in)    :: nptmass,npart
 real,    intent(in)    :: xyzh(:,:),rho(:),vxyz_ptmass(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 real,    intent(inout) :: eos_vars(:,:)
 select case(iH2R)
 case(1)
    call HII_Knn(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,rho,eos_vars)
 case(2,3)
    call HII_ray(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,rho,eos_vars)
 case(4)
    call HII_bound(nptmass,npart,xyzh,xyzmh_ptmass,vxyz_ptmass,vxyzu,rho,eos_vars)
 end select

end subroutine HII_feedback

!-----------------------------------------------------------------------
!+
! HII feedback routine using boundaries to spare the full knn
!+
!-----------------------------------------------------------------------
subroutine HII_bound(nptmass,npart,xyzh,xyzmh_ptmass,vxyz_ptmass,vxyzu,rho,eos_vars)
 use part,        only:ihacc,irateion,irstrom,itemp,imu
 use neighkdtree, only:listneigh,getneigh_pos,leaf_is_active
 use sortutils,   only:indexx
 use dim,         only:maxvxyzu,mpi
 use io,          only:iverbose,iprint,warning,id,nprocs
 use mpiutils,    only:reduceall_mpi,allgatherv_mpi,allgather_mpi,get_displs_mpi
 integer, intent(in)    :: nptmass,npart
 real,    intent(in)    :: xyzh(:,:),rho(:),vxyz_ptmass(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 real,    intent(inout) :: eos_vars(:,:)
 integer, parameter   :: maxc = 128, maxits = 50
 real, save           :: xyzcache(3,maxc)
 real,    allocatable :: dnpart(:,:),dnall(:),dr2all(:),rprev2(:)
 integer, allocatable :: id_shell(:)
 integer              :: nshell(nprocs),displs(nprocs)
 integer              :: i,k,j,nneigh,nneigh_shell,nneigh_inner,nshell_tot,its,ifirst,ilast,kstart
 real                 :: xi,yi,zi,log_Qi,coef,rs,r_lo,r_hi,deltar,vnorm
 real                 :: dn_lo,dn_hi,dntot,dnj,phi
 logical              :: converged,escaped,query

 !
 !-- reset every local particle to neutral before the sources ionise them again
 !
 eos_vars(imu,:) = gmw
 if (maxvxyzu < 4) eos_vars(itemp,:) = Tcold  ! cooling in isothermal case

 !
 !-- squared radii of the previous update (-1 if none): source s reached particle j
 !   last time if d_sj < r_S(s). Used for the flux sharing of overlapping regions;
 !   ptmass arrays are synced, so every rank gets the same shares
 !
 allocate(rprev2(nptmass))
 do i=1,nptmass
    if (xyzmh_ptmass(irateion,i) > epsilon(0.) .and. xyzmh_ptmass(irstrom,i) > 0.) then
       rprev2(i) = xyzmh_ptmass(irstrom,i)**2
    else
       rprev2(i) = -1.
    endif
 enddo

 sources: do i=1,nptmass
    log_Qi = xyzmh_ptmass(irateion,i)
    if (log_Qi < epsilon(log_Qi)) cycle sources
    xi   = xyzmh_ptmass(1,i)
    yi   = xyzmh_ptmass(2,i)
    zi   = xyzmh_ptmass(3,i)
    rs   = xyzmh_ptmass(irstrom,i)
    !-- photons used by particle j as a fraction of Q: y_j = coef*m_j*rho_j (eq. 3)
    coef = 10.**(DNcoeff - log_Qi)

    !
    !-- starting range [r_lo,r_hi] (eq. 7 and 8), identical on every rank
    !
    if (rs > 0.) then  ! D-type front: moves less than (2 c_ion + v_sink) dt
       vnorm  = sqrt(vxyz_ptmass(1,i)**2 + vxyz_ptmass(2,i)**2 + vxyz_ptmass(3,i)**2)
       !-- floor: the front is only known to about one smoothing length (stand-in for 2 h_edge)
       deltar = max((2.*csion + vnorm)*HIIdt, xyzmh_ptmass(ihacc,i), 0.01*rs)
       r_lo   = max(rs - deltar, 0.)
       r_hi   = min(rs + deltar, Rst_max)
    else               ! new source: continuous Stromgren radius from the nearest gas
       rs     = stromgren_guess(xi,yi,zi,xyzmh_ptmass(ihacc,i),log_Qi,xyzh,rho)
       deltar = rs
       r_lo   = 0.
       r_hi   = min(2.*rs, Rst_max)
    endif

    !
    !-- grow or shrink the range until W(r_lo) < 1 <= W(r_hi) (eq. 9)
    !
    converged = .false.
    escaped   = .false.
    query     = .true.
    its       = 0
    deltafind: do while (.not.converged)
       its = its + 1
       if (its > maxits) then
          call warning('HII_bound','range search did not converge',ival=i)
          cycle sources
       endif
       if (query) then  ! only needed when r_hi grows
          call getneigh_pos((/xi,yi,zi/),0.,r_hi,listneigh,nneigh,xyzcache,maxc,leaf_is_active)
          if (allocated(dnpart)) deallocate(dnpart)
          allocate(dnpart(2,nneigh))
       endif
       call compute_dnspheres(i,nptmass,xyzmh_ptmass,rprev2,coef,r_lo,r_hi,xyzh,rho,&
                              listneigh,nneigh,nneigh_shell,nneigh_inner,dn_lo,dn_hi,dnpart)
       if (mpi) then
          dn_lo = reduceall_mpi('+',dn_lo)
          dn_hi = reduceall_mpi('+',dn_hi)
       endif

       if (dn_hi < 1.) then        ! photons left over at r_hi: the region is larger
          if (r_hi >= Rst_max) then
             escaped   = .true.    ! a fraction 1 - dn_hi of the photons escapes
             converged = .true.
          else
             r_lo   = r_hi
             deltar = 2.*deltar
             r_hi   = min(r_hi + deltar, Rst_max)
             query  = .true.
          endif
       elseif (dn_lo >= 1.) then   ! the region shrank inside r_lo: no new query needed
          r_hi   = r_lo
          deltar = 2.*deltar
          r_lo   = max(r_lo - deltar, 0.)
          query  = .false.
       else
          converged = .true.
       endif
    enddo deltafind

    !
    !-- ionise the local particles inside r_lo (inside r_hi if photons escape);
    !   listneigh = [ shell | inner | rejected ], so no distance is recomputed
    !
    if (escaped) then
       kstart = 1
    else
       kstart = nneigh_shell + 1
    endif
    do k=kstart,nneigh_shell+nneigh_inner
       call ionise_particle(eos_vars,vxyzu,listneigh(k),1.)
    enddo
    if (escaped) then
       xyzmh_ptmass(irstrom,i) = Rst_max
       if (iverbose >= 2) write(iprint,*) 'HII_bound: photons escape from sink ',i,', fraction ',1.-dn_hi
       cycle sources
    endif

    !
    !-- gather the shell [r_lo,r_hi) on every rank: entry j of rank p sits at displs(p+1)+j
    !
    if (mpi) then
       call allgather_mpi(nneigh_shell,nshell)
       displs     = get_displs_mpi(nshell)
       nshell_tot = sum(nshell)
    else
       nshell(1)  = nneigh_shell
       displs(1)  = 0
       nshell_tot = nneigh_shell
    endif
    allocate(dnall(2*nshell_tot),dr2all(nshell_tot),id_shell(nshell_tot))
    if (mpi) then
       call allgatherv_mpi(reshape(dnpart(:,1:nneigh_shell),(/2*nneigh_shell/)),dnall,&
                           counts=2*nshell,displs=2*displs)
    else
       dnall = reshape(dnpart(:,1:nneigh_shell),(/2*nneigh_shell/))
    endif
    dr2all = dnall(1:2*nshell_tot-1:2)

    !
    !-- every rank sorts an index array (the gathered data keep their order, so displs
    !   still gives the owner) and marches the shell; the boundary particle gets eq. 12
    !
    call indexx(nshell_tot,dr2all,id_shell)
    ifirst = displs(id+1) + 1
    ilast  = displs(id+1) + nshell(id+1)
    dntot  = dn_lo
    do k=1,nshell_tot
       j   = id_shell(k)
       dnj = dnall(2*j)
       if (dntot + dnj <= 1.) then
          phi = 1.
       else
          phi = (1. - dntot)/dnj
       endif
       if (j >= ifirst .and. j <= ilast) &
          call ionise_particle(eos_vars,vxyzu,listneigh(j-displs(id+1)),phi)
       rs    = sqrt(dr2all(j))
       dntot = dntot + dnj
       if (phi < 1.) exit
    enddo
    xyzmh_ptmass(irstrom,i) = rs
    if (iverbose >= 2) write(iprint,*) 'HII_bound: Rstrom from sink ',i,' = ',rs,&
                                       ' shell size ',nshell_tot,' in ',its,' iterations'
    deallocate(dnall,dr2all,id_shell)
 enddo sources
 if (allocated(dnpart)) deallocate(dnpart)
 deallocate(rprev2)

end subroutine HII_bound


!-----------------------------------------------------------------------
!+
!  photon demand of source i, as a fraction of its Q, inside r_lo and
!  inside r_hi (eq. 6). A particle that other sources also reached at the
!  previous update costs source i only its flux share (eq. 14):
!    f_ij = F_ij / sum_s F_sj,  F_sj = Q_s/(4 pi d_sj^2),
!  the sum running over i and every source s with d_sj < r_S(s) last update.
!  listneigh is partitioned in place into three blocks:
!    1 : nneigh_shell                          shell [r_lo,r_hi), same order as dnpart
!    nneigh_shell+1 : nneigh_shell+nneigh_inner inside r_lo
!    the rest                                   beyond r_hi, not gas, or dead
!+
!-----------------------------------------------------------------------
subroutine compute_dnspheres(i,nptmass,xyzmh_ptmass,rprev2,coef,r_lo,r_hi,xyzh,rho,&
                             listneigh,nneigh,nneigh_shell,nneigh_inner,dn_lo,dn_hi,dnpart)
 use part, only:get_pmass,iphase,get_partinfo,isdead_or_accreted,irateion
 integer, intent(in)    :: i,nptmass,nneigh
 integer, intent(inout) :: listneigh(:)
 integer, intent(out)   :: nneigh_shell,nneigh_inner
 real,    intent(in)    :: xyzmh_ptmass(:,:),rprev2(:),coef,r_lo,r_hi
 real,    intent(in)    :: xyzh(:,:),rho(:)
 real,    intent(out)   :: dn_lo,dn_hi,dnpart(:,:)
 integer :: j,n,s,ia,ib,itypej
 real    :: dr2,r_hi2,r_lo2,dnj,xi,yi,zi,d2s,fsum
 logical :: isactive,isgas,isdust

 xi = xyzmh_ptmass(1,i)
 yi = xyzmh_ptmass(2,i)
 zi = xyzmh_ptmass(3,i)
 dn_lo = 0.
 dn_hi = 0.
 nneigh_shell = 0
 nneigh_inner = 0
 r_hi2 = r_hi**2
 r_lo2 = r_lo**2

 do n=1,nneigh
    j = listneigh(n)
    call get_partinfo(iphase(j),isactive,isgas,isdust,itypej)
    if (.not.isgas .or. isdead_or_accreted(xyzh(4,j))) cycle
    dr2 = (xyzh(1,j)-xi)**2 + (xyzh(2,j)-yi)**2 + (xyzh(3,j)-zi)**2
    if (dr2 >= r_hi2) cycle  ! the tree returns whole leaves, beyond r_hi too
    !-- sum_s F_sj / F_ij, with Q_s/Q_i = 10**(logQ_s - logQ_i) to keep numbers O(1)
    fsum = 1.
    do s=1,nptmass
       if (s == i .or. rprev2(s) < 0.) cycle  ! rprev2 < 0: not a source, or no radius yet
       d2s = (xyzh(1,j)-xyzmh_ptmass(1,s))**2 + (xyzh(2,j)-xyzmh_ptmass(2,s))**2 &
           + (xyzh(3,j)-xyzmh_ptmass(3,s))**2
       if (d2s < rprev2(s)) fsum = fsum + 10.**(xyzmh_ptmass(irateion,s)-xyzmh_ptmass(irateion,i)) &
                                          *dr2/max(d2s,tiny(d2s))
    enddo
    dnj   = coef*get_pmass(j,.true.)*rho(j)/fsum
    dn_hi = dn_hi + dnj
    ib = nneigh_shell + nneigh_inner + 1  ! first rejected slot (or n itself)
    if (dr2 < r_lo2) then
       dn_lo = dn_lo + dnj
       nneigh_inner = nneigh_inner + 1
       listneigh(n)  = listneigh(ib)        ! swap j into the end of the inner block
       listneigh(ib) = j
    else
       ia = nneigh_shell + 1                ! first inner slot (= ib if no inner yet)
       listneigh(n)  = listneigh(ib)        ! rotate: first inner -> end of inner block,
       listneigh(ib) = listneigh(ia)        !         rejected -> n, j -> end of shell
       listneigh(ia) = j
       nneigh_shell  = nneigh_shell + 1
       dnpart(1,nneigh_shell) = dr2
       dnpart(2,nneigh_shell) = dnj
    endif
 enddo

end subroutine compute_dnspheres

!-----------------------------------------------------------------------
!+
!  first guess of the Stromgren radius of a new source (eq. 8), from the
!  density of the nearest gas particle over all ranks
!+
!-----------------------------------------------------------------------
real function stromgren_guess(xi,yi,zi,hacc,log_Qi,xyzh,rho) result(rs)
 use neighkdtree, only:listneigh,getneigh_pos,leaf_is_active
 use mpiutils,    only:reduceall_mpi
 use dim,         only:mpi
 use units,       only:utime
 use physcon,     only:pi
 use part,        only:iphase,get_partinfo,isdead_or_accreted
 real, intent(in) :: xi,yi,zi,hacc,log_Qi
 real, intent(in) :: xyzh(:,:),rho(:)
 integer, parameter :: maxc = 128
 real    :: xyzcache(3,maxc)
 integer :: k,j,nneigh,itypej
 real    :: rq,dr2,d2loc,d2min,rholoc,rhonear,nH,Q
 logical :: isactive,isgas,isdust

 rs = Rst_max
 rq = max(hacc,1.e-3*Rst_max)
 do while (rq < Rst_max)
    call getneigh_pos((/xi,yi,zi/),0.,rq,listneigh,nneigh,xyzcache,maxc,leaf_is_active)
    d2loc  = huge(d2loc)
    rholoc = 0.
    do k=1,nneigh
       j = listneigh(k)
       call get_partinfo(iphase(j),isactive,isgas,isdust,itypej)
       if (.not.isgas .or. isdead_or_accreted(xyzh(4,j))) cycle
       dr2 = (xyzh(1,j)-xi)**2 + (xyzh(2,j)-yi)**2 + (xyzh(3,j)-zi)**2
       if (dr2 < d2loc) then
          d2loc  = dr2
          rholoc = rho(j)
       endif
    enddo
    if (mpi) then
       d2min = reduceall_mpi('min',d2loc)
    else
       d2min = d2loc
    endif
    if (d2min < rq**2) then
       if (mpi) then  ! density of the nearest particle, held by whichever rank owns it
          if (d2loc > d2min) rholoc = 0.
          rhonear = reduceall_mpi('max',rholoc)
       else
          rhonear = rholoc
       endif
       nH = rhonear*mH1                 ! number density, code units
       Q  = 10.**log_Qi*utime           ! photons per code time
       rs = (3.*Q/(4.*pi*ar*nH**2))**(1./3.)
       return
    endif
    rq = 2.*rq
 enddo

end function stromgren_guess

!-----------------------------------------------------------------------
!+
!  ionise particle j with ionised fraction phi (1 = fully ionised),
!  never lowering an ionisation already set by another source
!+
!-----------------------------------------------------------------------
subroutine ionise_particle(eos_vars,vxyzu,j,phi)
 use part, only:itemp,imu
 use dim,  only:maxvxyzu
 integer, intent(in)    :: j
 real,    intent(in)    :: phi
 real,    intent(inout) :: eos_vars(:,:),vxyzu(:,:)
 real :: phi_old

 phi_old = (gmw - eos_vars(imu,j))/(gmw - muion)
 if (phi <= phi_old) return

 eos_vars(itemp,j) = phi*Tion  + (1. - phi)*Tcold
 eos_vars(imu,j)   = phi*muion + (1. - phi)*gmw
 if (maxvxyzu >= 4) vxyzu(4,j) = phi*uIon + (1. - phi)*vxyzu(4,j)

end subroutine ionise_particle

!-----------------------------------------------------------------------
!+
! HII feedback routine using K nearest neighbors prescription
!+
!-----------------------------------------------------------------------
subroutine HII_Knn(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,rho,eos_vars)
 use part,       only:massoftype,ihsoft,igas,irateion,isdead_or_accreted,&
                      irstrom,get_partinfo,iphase,itemp,imu,get_pmass
 use neighkdtree,   only:listneigh=>listneigh_global,getneigh_pos,leaf_is_active
 use sortutils,  only:Knnfunc,set_r2func_origin,r2func_origin
 use physcon,    only:pc,pi
 use dim,        only:maxvxyzu
 use io,         only:iverbose,iprint,warning,id,master
 integer, intent(in)    :: nptmass,npart
 real,    intent(in)    :: xyzh(:,:),rho(:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 real,    intent(inout) :: eos_vars(:,:)
 integer, parameter :: maxc  = 128
 real, save :: xyzcache(3,maxc)
 integer            :: i,k,j,npartin,nneigh,itypej,its
 real               :: pmass,Ndot,DNdot,DNratio,logNdiff,taud,mHII,r,r_in,hcheck
 real               :: xi,yi,zi,log_Qi,xj,yj,zj,dx,dy,dz,vkx,vky,vkz
 logical            :: momflag,isactive,isgas,isdust,converged,max_search

 momflag = .false.
 r = 0.
 r_in = 0.
 pmass = massoftype(igas)

 !
 !-- Rst derivation and thermal feedback
 !
 if (nHIIsources > 0 .and. id==master) then
    eos_vars(imu,:) = gmw
    if (maxvxyzu < 4) eos_vars(itemp,:) = Tcold  ! cooling in isothermal case
    do i=1,nptmass
       converged  = .false.
       max_search = .false.
       its        = 0
       npartin    = 0
       log_Qi     = xyzmh_ptmass(irateion,i)
       if (log_Qi <epsilon(log_Qi)) cycle
       Ndot   = log_Qi
       xi     = xyzmh_ptmass(1,i)
       yi     = xyzmh_ptmass(2,i)
       zi     = xyzmh_ptmass(3,i)
       hcheck = xyzmh_ptmass(irstrom,i)
       iterations:do while (.not.max_search)
          if (hcheck > 0. ) then
             hcheck = hcheck + 6.*csion*HIIdt
             if (hcheck > Rmax) then
                hcheck = Rmax
                max_search = .true.
             endif
          else
             hcheck = Rmax
             max_search = .true.
          endif

          its = its + 1

          call getneigh_pos((/xi,yi,zi/),0.,hcheck,listneigh,nneigh,xyzcache,maxc,leaf_is_active)
          call Knnfunc(nneigh,(/xi,yi,zi/),xyzh,listneigh)

          if (nneigh > 0) then
             neigh_cycle:do k=1,nneigh
                j = listneigh(k)
                call get_partinfo(iphase(j),isactive,isgas,isdust,itypej)
                if (isgas) then
                   if (.not. isdead_or_accreted(xyzh(4,j))) then
                      ! ionising photons needed to fully ionise the current particle
                      DNdot = DNcoeff + log10(get_pmass(j,.true.)*rho(j))
                      if (DNdot < Ndot) then
                         if (eos_vars(imu,j) - muion > epsilon(muion) ) then ! is ionised ?
                            logNdiff = DNdot - Ndot
                            Ndot = Ndot + log10(1.-10.**(logNdiff))
                            eos_vars(itemp,j)  = Tion
                            eos_vars(imu,j) = muion
                            npartin = npartin + 1
                            if (maxvxyzu >= 4) vxyzu(4,j) = uIon
                         endif
                      else
                         if (k > 1) then ! end of the HII region
                            if (maxvxyzu >= 4) then
                               DNratio = Ndot/DNdot
                               eos_vars(itemp,j)  = DNratio * Tion + (1.0 - DNratio) * Tcold
                               eos_vars(imu,j)    = DNratio * muion + (1.0 - DNratio) * gmw
                               vxyzu(4,j)         = DNratio * uIon + (1.0 - DNratio) * vxyzu(4,j)
                               Ndot = 0.
                            endif
                            r = sqrt((xi-xyzh(1,j))**2 + (yi-xyzh(2,j))**2 + (zi-xyzh(3,j))**2)
                            j = listneigh(1)
                         else ! unresolved case
                            r = epsilon(r)
                         endif
                         xyzmh_ptmass(irstrom,i) = r
                         converged = .true.
                         max_search = .true.
                         exit neigh_cycle
                      endif
                   endif
                endif
             enddo neigh_cycle
          endif
       enddo iterations
       if (iverbose == 2) write(iprint,*)'Rstrom from sink ',i,' = ',r," with N = ",npartin,&
                                         ' ionised particle, remaining Nphot : ',Ndot,DNdot,"in ",its," iterations."
       if (.not.converged) call warning('HII_feedback','Photon march did not converge...',var='Ndot',&
                                         val=((log_Qi-Ndot)/log_Qi),ival=its)
       !
       !-- Momentum feedback
       !
       if (momflag .and. npartin > 3) then
          j = listneigh(1)
          r_in = sqrt((xi-xyzh(1,j))**2 + (yi-xyzh(2,j))**2 + (zi-xyzh(3,j))**2)
          mHII = ((4.*pi*(r**3-r_in**3)*rho(j))/3)
          if (mHII>3*pmass) then
             !$omp parallel do default(none) &
             !$omp shared(mHII,xyzh,listneigh,sigd,HIIdt,rho) &
             !$omp shared(mH,vxyzu,log_Qi,hv_on_c,npartin,pmass,xi,yi,zi) &
             !$omp private(j,dx,dy,dz,vkx,vky,vkz,xj,yj,zj,r,taud)
             do k=1,npartin
                j = listneigh(k)
                xj = xyzh(1,j)
                yj = xyzh(2,j)
                zj = xyzh(3,j)
                dx = xj - xi
                dy = yj - yi
                dz = zj - zi
                r = dx**2 + dy**2 + dz**2
                taud = (rho(j)/mH)*sigd*r
                if (taud > 1.97) taud=1.97
                vkz = (1.+1.5*exp(-taud))*((10**log_Qi)/mHII)*hv_on_c*(dz/r)
                vkx = (1.+1.5*exp(-taud))*((10**log_Qi)/mHII)*hv_on_c*(dx/r)
                vky = (1.+1.5*exp(-taud))*((10**log_Qi)/mHII)*hv_on_c*(dy/r)
                vxyzu(1,j) = vxyzu(1,j) +  vkx*HIIdt
                vxyzu(2,j) = vxyzu(2,j) +  vky*HIIdt
                vxyzu(3,j) = vxyzu(3,j) +  vkz*HIIdt
             enddo
             !$omp end parallel do
          endif
       endif
    enddo
 endif

end subroutine HII_Knn

!-----------------------------------------------------------------------
!+
! HII feedback routine using inversed ray tracing method
!+
!-----------------------------------------------------------------------
subroutine HII_ray(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,rho,eos_vars)
 use part,     only:massoftype,igas,irateion,irstrom,isdead_or_accreted,&
                    iphase,get_partinfo,noverlap,itemp,imu
 use dim,      only:maxvxyzu
 use io,         only:iverbose,iprint,id,master
 integer, intent(in)    :: nptmass,npart
 real,    intent(in)    :: xyzh(:,:),rho(:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 real,    intent(inout) :: eos_vars(:,:)
 real, save         :: xyzcache(4,maxcache)
 real               :: pmass,log_Qi,fluxi,xj,yj,zj,rsrci
 logical            :: isactive,isgas,isdust
 integer            :: i,j,itypei,noverlapi,k
 !$omp threadprivate(xyzcache)

 if (nHIIsources > 0 .and. id==master) then
    pmass = massoftype(igas)
    xyzmh_ptmass(irstrom,:) = -1.
    k=0
    !
    !-- Rst derivation and thermal feedback
    !
!$omp parallel default(none)&
!$omp shared(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,eos_vars,noverlap,rho)&
!$omp shared(pmass,iphase,iH2R,Tcold,uIon,gmw)&
!$omp private(log_Qi,i,j,xj,yj,zj,fluxi,itypei,rsrci,isgas,isdust)&
!$omp private(isactive,noverlapi)&
!$omp reduction(+:k)
!$omp do
    do i=1,npart
       if (.not. isdead_or_accreted(xyzh(4,i))) then
          call get_partinfo(iphase(i),isactive,isgas,isdust,itypei)
          ! if (.not. isactive) cycle
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
                call inversed_raytracing(i,(/xj,yj,zj/),xyzh,rho,xyzcache,noverlap,&
                                            pmass,log_Qi,fluxi,rsrci)
                if (fluxi > 0.) then
                   eos_vars(itemp,i) = Tion
                   eos_vars(imu,i) = muion
                   noverlapi = noverlapi + 1
                   if (maxvxyzu >= 4) vxyzu(4,i) = uIon
                   k = k + 1
                   !$omp atomic update
                   xyzmh_ptmass(irstrom,j) = max(rsrci,xyzmh_ptmass(irstrom,j))
                   !$omp end atomic
                endif
             enddo
             noverlap(i) = noverlapi
          endif
       endif
    enddo
!$omp enddo
!$omp end parallel
    if (iverbose == 2) write(iprint,*) 'N ionised = ',k
 endif

end subroutine HII_ray

!-----------------------------------------------------------------------
!+
! find nearest particles along the ray from target to source
!+
!-----------------------------------------------------------------------
subroutine inversed_raytracing(itarg,srcpos,xyzh,rho,xyzcache,noverlap,pmass,log_Q,flux,rsrc)
 use part,     only:iphase,get_partinfo
 use neighkdtree, only:getneigh_pos,leaf_is_active,listneigh
 use kernel,   only:radkern
 use physcon,  only:fourpi
 use units,    only:utime
 integer, intent(in)    :: itarg
 integer, intent(in)    :: noverlap(:)
 real,    intent(in)    :: xyzh(:,:),rho(:)
 real,    intent(inout) :: xyzcache(:,:)

 real, intent(in)  :: srcpos(3),pmass,log_Q
 real, intent(out) :: flux,rsrc
 real, parameter :: thetalim2 = 0.0025
 real            :: rayx,rayy,rayz,drisrc1,drisrc,dr,hpmass1,lumS,recvol
 real            :: xij,yij,zij,xi,yi,zi,hi,hj,drproj,dr2toray,icoeff
 real            :: dr2ij,dr2toraymin,intensity,hi2,drprojmin,hcheck
 integer         :: nneigh,k,i,j,knext,inext,itypei,noverlapj
 logical         :: isactive,isgas,isdust,unreached

 hpmass1 = 0.5*mH1

 lumS    = ((10**log_Q)*utime)/fourpi
 flux = -1.
 xi      = xyzh(1,itarg)
 yi      = xyzh(2,itarg)
 zi      = xyzh(3,itarg)
 hi      = xyzh(4,itarg)
 rayx    = srcpos(1)-xi
 rayy    = srcpos(2)-yi
 rayz    = srcpos(3)-zi
 drisrc  = sqrt(rayx*rayx + rayy*rayy + rayz*rayz)
 if (drisrc < epsilon(drisrc)) return
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
 rsrc = drisrc

 reachsrc:do while (unreached)
    call getneigh_pos((/xi,yi,zi/),0.,hcheck,listneigh,&
                      nneigh,xyzcache,maxcache,leaf_is_active)
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
          if (k <= maxcache) then
             xij = xyzcache(1,k) - xi
             yij = xyzcache(2,k) - yi
             zij = xyzcache(3,k) - zi
          else
             xij = xyzh(1,j) - xi
             yij = xyzh(2,j) - yi
             zij = xyzh(3,j) - zi
          endif
          dr2ij   = xij*xij + yij*yij + zij*zij
          if ( dr2ij < hi2) then
             drproj = rayx*xij + rayy*yij + rayz*zij
             if (drproj > 0.) then
                dr2toray = (xij - drproj*rayx)**2 + (yij - drproj*rayy)**2 + (zij - drproj*rayz)**2
                if (dr2toray < dr2toraymin) then
                   dr2toraymin = dr2toray
                   drprojmin = drproj
                   knext = k
                   inext = j
                   if ( dr2toraymin < (drprojmin*drprojmin)*thetalim2) exit neighsearch
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
       if (k<=maxcache) then
          hj  = 1./xyzcache(4,knext)
       else
          hj  = xyzh(4,inext)
       endif
       !
       !-- compute partial recombination volume if overlapped region
       !
       recvol = ar*((rho(i)+rho(inext))*hpmass1)**2
       noverlapj = noverlap(inext)
       if (noverlapj>0) recvol = recvol*(1./noverlap(inext))
       !
       !-- integrate intensity
       !
       intensity = intensity + drprojmin*recvol*dr**2

       if (lumS < intensity) exit reachsrc! exit if already neutral
       !
       !-- prepare the next search
       !
       if (knext <= maxcache) then
          xi = xyzcache(1,knext)
          yi = xyzcache(2,knext)
          zi = xyzcache(3,knext)
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
    drprojmin = dr
    recvol = ar*(2.*rho(i)*hpmass1)**2

    noverlapj = noverlap(inext)
    if (noverlapj>0) recvol = recvol/noverlap(inext)

    intensity = intensity + drprojmin*recvol*(icoeff*dr)**2
    flux = (lumS-intensity)
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
    "HII region expansion feedback in star forming region(1 Knn method,2 inv RT method,4 bounded Knn)", iunit)
    call write_inopt(Mmin, 'Mmin', "Minimum star mass to trigger HII region (MSun)", iunit)
    call write_inopt(Rmax, 'Rmax', "Maximum radius for HII region (pc)(only Knn)", iunit)
 endif

end subroutine write_options_H2R

!-----------------------------------------------------------------------
!+
!  read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_H2R(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(iH2R,'iH2R',db,errcount=nerr,min=0,default=0)
 if (iH2R > 0) then
    call read_inopt(Mmin,'Mmin',db,errcount=nerr,min=8.)
    call read_inopt(Rmax,'Rmax',db,errcount=nerr,min=10.)
 endif

end subroutine read_options_H2R

end module HIIRegion
