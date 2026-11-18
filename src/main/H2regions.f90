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
! Two phase thermal model (ieos=21): the fully ionised gas is held at uIon
! and the neutral gas away from the front at Tcold (its hydro du/dt is
! ignored in force). The gas of the front (partially ionised, or neutral
! with an ionised neighbour) is heated at once to its target
! T = phi*Tion + (1-phi)*Tcold and, above it, cools over Crelax*h/c_s(u)
!
! :References: Fujii et al. (2021), Hopkins et al. (2012)
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, eos_HIIR, healpix, infile_utils, io, kdtree,
!   kernel, mpiutils, neighkdtree, part, physcon, sortutils, units
!
 use eos_HIIR,   only:csion,uIon,muion
 use eos,        only:gmw
 implicit none

 public :: update_ionrates,update_ionrate, HII_feedback,&
           initialize_H2R,read_options_H2R,write_options_H2R
 public :: relax_u_HIIR,heated_u_HIIR,is_ionised,is_fully_ionised

 integer, parameter, public    :: HIIuprate   = 8 ! update rate when IND_TIMESTEPS=yes (Hopkins like method)
 integer, public               :: iH2R = 0
 real   , public               :: Crelax = 1. ! relaxation time of the gas of the front in units of h/c_s
 real   , public               :: Rmax = 15 ! Maximum HII region radius (pc) to avoid artificial expansion...
 real   , public               :: Mmin = 8  ! Minimum mass (Msun) to produce HII region
 integer, public               :: nHIIsources = 0
 real   , public               :: ar
 real   , public               :: mH
 real   , public               :: mH1
 real   , public               :: DNcoeff
 logical, public               :: HIIupdateflag
 real   , public               :: HIIdt
 integer, public               :: nside_H2R = 1   ! HEALPix Nside of the sectors (iH2R=4): 12*Nside**2 sectors
 integer, public               :: nsec_H2R  = 12
 real   , allocatable, public  :: rstrom_sec(:,:) ! Stromgren radius of each sector of each sink (nsec,maxptmass)
 logical, allocatable, public  :: sectors_on(:)   ! sink solved in sectors (else isotropic) (maxptmass)
 real   , allocatable          :: secdir(:,:)     ! unit vector of the centre of each sector (3,nsec)
 real   , allocatable          :: secbeta(:)      ! angular radius of each sector: centre to farthest corner

 integer, parameter :: maxcache  = 1024

 !-- (H2R == 4)
 integer, parameter :: nsec_on  = 4 ! sector on if nsec_on*nsec = n_hot
 integer, parameter :: nsec_off = 2 ! sector off if nsec_off*nsec = n_hot

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
 real, parameter   :: tolphase = 1.e-3 ! relative excess of u below which the gas is back in its phase

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
 use dim,     only:maxptmass
 use healpix, only:pix2vec_nest,vec2pix_nest
 use part,    only:eos_vars,imu
 integer :: k,iv,ipix
 real    :: vertex(3,4)

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
 eos_vars(imu,:) = gmw  ! neutral until the first update

 !
 !-- HEALPix sectors: centres and angular radii
 !
 if (iH2R == 4) then
    nsec_H2R = 12*nside_H2R**2

    !-- sector allocation
    if (allocated(secdir)) deallocate(secdir)
    if (allocated(secbeta)) deallocate(secbeta)
    if (allocated(rstrom_sec)) deallocate(rstrom_sec)
    if (allocated(sectors_on)) deallocate(sectors_on)

    allocate(secdir(3,nsec_H2R))
    allocate(secbeta(nsec_H2R))
    allocate(rstrom_sec(nsec_H2R,maxptmass))
    allocate(sectors_on(maxptmass))

    !-- construct HEALPix grid
    do k=1,nsec_H2R
       call pix2vec_nest(nside_H2R,k-1,secdir(:,k),vertex)
       secbeta(k) = 0.
       do iv=1,4
          secbeta(k) = max(secbeta(k),acos(min(1.,dot_product(secdir(:,k),vertex(:,iv)))))
       enddo
    enddo
    rstrom_sec = -1.
    sectors_on = .true.
    !-- vec2pix_nest fills its lookup table on the first call: do it here, outside OpenMP
    call vec2pix_nest(nside_H2R,secdir(:,1),ipix)
 endif

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
 use part,        only:ihacc,irateion,irstrom,imu
 use neighkdtree, only:listneigh=>listneigh_global,getneigh_sectors,getneigh_pos,&
                       leaf_is_active,node
 use kdtree,      only:irootnode
 use sortutils,   only:indexx
 use dim,         only:mpi
 use io,          only:iverbose,iprint,warning,id,nprocs
 use mpiutils,    only:reduceall_mpi,allgatherv_mpi,allgather_mpi,get_displs_mpi
 integer, intent(in)    :: nptmass,npart
 real,    intent(in)    :: xyzh(:,:),rho(:),vxyz_ptmass(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyzu(:,:)
 real,    intent(inout) :: eos_vars(:,:)
 integer, parameter   :: maxits = 50, maxc = 128
 real,    allocatable :: dnpart(:,:),dnall(:),rprev2(:,:),rprevmax2(:),dr2sec(:)
 real,    allocatable :: dn_lo(:),dn_hi(:),r_lo(:),r_hi(:),deltar(:),rs(:)
 integer, allocatable :: id_shell(:),isort(:),nsecall(:),offsec(:)
 logical, allocatable :: converged(:),escaped(:)
 integer              :: nshell(nprocs),displs(nprocs)
 integer              :: i,k,j,m,n,nsec,nneigh,nneigh_shell,nneigh_inner,nshell_tot,its,ifirst,ilast,nk
 integer              :: nhot,nhot_loc
 real                 :: xsrc(3),log_Qi,coef,rs0,vnorm,ddyn,rbound,dntot,dnj,phi
 real                 :: xyzcache(3,maxc)
 logical              :: query,warm

 !
 !-- reset every local particle to neutral before the sources ionise them again
 !
 eos_vars(imu,:) = gmw

 ! squared sector radii of the previous update (-1 if none)
 ! used for the flux sharing of overlapping regions
 allocate(rprev2(nsec_H2R,nptmass),rprevmax2(nptmass))
 do i=1,nptmass
    if (xyzmh_ptmass(irateion,i) > epsilon(0.) .and. xyzmh_ptmass(irstrom,i) > 0. &
        .and. all(rstrom_sec(:,i) > 0.)) then
       rprev2(:,i)  = rstrom_sec(:,i)**2
       rprevmax2(i) = maxval(rprev2(:,i))
    else
       rprev2(:,i)  = -1.
       rprevmax2(i) = -1.
    endif
 enddo

 sources: do i=1,nptmass

    log_Qi = xyzmh_ptmass(irateion,i)
    if (log_Qi < epsilon(log_Qi)) cycle sources !-- not bright enough...
    xsrc = xyzmh_ptmass(1:3,i)

    warm = (xyzmh_ptmass(irstrom,i) > 0. .and. all(rstrom_sec(:,i) > 0.))
    if (.not.warm) sectors_on(i) = .true.

    if (sectors_on(i)) then
       nsec = nsec_H2R
    else
       nsec = 1
    endif

    if (allocated(dn_lo))    deallocate(dn_lo)
    if (allocated(dn_hi))    deallocate(dn_hi)
    if (allocated(r_lo))     deallocate(r_lo)
    if (allocated(r_hi))     deallocate(r_hi)
    if (allocated(deltar))   deallocate(deltar)
    if (allocated(rs))       deallocate(rs)
    if (allocated(converged)) deallocate(converged)
    if (allocated(escaped))  deallocate(escaped)
    if (allocated(nsecall))  deallocate(nsecall)
    if (allocated(offsec))   deallocate(offsec)

    allocate(dn_lo(nsec))
    allocate(dn_hi(nsec))
    allocate(r_lo(nsec))
    allocate(r_hi(nsec))
    allocate(deltar(nsec))
    allocate(rs(nsec))
    allocate(converged(nsec))
    allocate(escaped(nsec))
    allocate(nsecall(nsec))
    allocate(offsec(nsec+1))

    !-- photons used by particle j as a fraction of its sector's budget Q/nsec (eq. 3, 13)
    coef = real(nsec)*10.**(DNcoeff - log_Qi)

    !-- distance beyond which there is no gas: a sector reaching it lets its photons escape
    if (npart > 0) then
       rbound = sqrt(sum((node(irootnode)%xcen - xsrc)**2)) + node(irootnode)%size
    else
       rbound = 0.
    endif
    if (mpi) rbound = reduceall_mpi('max',rbound)  ! same on every rank, also those without gas
    !
    !-- starting range [r_lo,r_hi] of each sector (eq. 7 and 8), identical on every rank
    !
    if (warm) then
       !-- D-type front: moves less than (2 c_ion + v_sink) dt; floor: the front is only
       !   known to about one smoothing length (stand-in for 2 h_edge)
       vnorm = sqrt(vxyz_ptmass(1,i)**2 + vxyz_ptmass(2,i)**2 + vxyz_ptmass(3,i)**2)
       ddyn  = (2.*csion + vnorm)*HIIdt
       if (nsec == 1) then  ! isotropic: the range encloses every sector of the last update
          deltar(1) = max(ddyn, xyzmh_ptmass(ihacc,i), 0.01*maxval(rstrom_sec(:,i)))
          r_lo(1)   = max(minval(rstrom_sec(:,i)) - deltar(1), 0.)
          r_hi(1)   = min(maxval(rstrom_sec(:,i)) + deltar(1), rbound)
       else
          do k=1,nsec
             deltar(k) = max(ddyn, xyzmh_ptmass(ihacc,i), 0.01*rstrom_sec(k,i))
             r_lo(k)   = max(rstrom_sec(k,i) - deltar(k), 0.)
             r_hi(k)   = min(rstrom_sec(k,i) + deltar(k), rbound)
          enddo
       endif
    else  ! new source: continuous Stromgren radius from the nearest gas, in every sector
       rs0    = stromgren_guess(xsrc(1),xsrc(2),xsrc(3),xyzmh_ptmass(ihacc,i),log_Qi,xyzh,rho,rbound)
       deltar = rs0
       r_lo   = 0.
       r_hi   = min(2.*rs0, rbound)
    endif

    !
    !-- grow or shrink the range of each sector until W(r_lo) < 1 <= W(r_hi) (eq. 9);
    !   the tree is walked again whenever one sector's r_hi grows
    !
    converged = .false.
    escaped   = .false.
    query     = .true.
    its       = 0
    deltafind: do while (.not.all(converged))
       its = its + 1
       if (its > maxits) then
          call warning('HII_bound','range search did not converge',ival=i)
          cycle sources
       endif
       if (query) then

          if (nsec > 1) then ! solid angle sector around the star
             call getneigh_sectors(xsrc,nsec,secdir,secbeta,r_hi,listneigh,nneigh)
          else  ! isotropic: a sphere of radius r_hi
             call getneigh_pos(xsrc,0.,r_hi(1),listneigh,nneigh,xyzcache,maxc,leaf_is_active)
          endif

          if (allocated(dnpart)) deallocate(dnpart)
          allocate(dnpart(3,nneigh))
       endif

       !-- Loop over neighbours to find r_lo dn_lo and particles in the shell to cycle to
       call compute_dnspheres(i,nptmass,xyzmh_ptmass,rprev2,rprevmax2,coef,r_lo,r_hi,xyzh,rho,&
                              listneigh,nneigh,nneigh_shell,nneigh_inner,dn_lo,dn_hi,dnpart)

       if (mpi) then
          dn_lo = reduceall_mpi('+',dn_lo)
          dn_hi = reduceall_mpi('+',dn_hi)
       endif

       query = .false.
       do k=1,nsec
          if (converged(k)) cycle
          if (dn_hi(k) < 1.) then        ! photons left over at r_hi: the sector reaches further
             if (r_hi(k) >= rbound) then  ! no gas left: a fraction 1 - dn_hi escapes
                if (r_lo(k) < r_hi(k)) then
                   r_lo(k) = r_hi(k)     ! one more pass: the whole sector moves to the inner block
                else
                   rs(k) = r_hi(k)
                   escaped(k)   = .true.
                   converged(k) = .true.
                endif
             else
                r_lo(k)   = r_hi(k)
                deltar(k) = 2.*deltar(k)
                r_hi(k)   = min(r_hi(k) + deltar(k), rbound)
                query     = .true.
             endif
          elseif (dn_lo(k) >= 1.) then   ! the sector shrank inside r_lo: no new walk needed
             r_hi(k)   = r_lo(k)
             deltar(k) = 2.*deltar(k)
             r_lo(k)   = max(r_lo(k) - deltar(k), 0.)
          else
             converged(k) = .true.
          endif
       enddo
    enddo deltafind

    !
    !-- ionise the local particles inside the r_lo of their sector;
    !
    nhot_loc = nneigh_inner
    !$omp parallel do default(none) &
    !$omp shared(nneigh_shell,nneigh_inner,listneigh,eos_vars,vxyzu) &
    !$omp private(n)
    do n=nneigh_shell+1,nneigh_shell+nneigh_inner
       call ionise_particle(eos_vars,vxyzu,listneigh(n),1.)
    enddo
    !$omp end parallel do

    !
    !-- gather the shell on every rank: (d^2, y_j, sector) of entry j of rank p sits at displs(p+1)+j
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
    allocate(dnall(3*nshell_tot),dr2sec(nshell_tot),id_shell(nshell_tot),isort(nshell_tot))
    if (mpi) then
       call allgatherv_mpi(reshape(dnpart(:,1:nneigh_shell),(/3*nneigh_shell/)),dnall,&
                           counts=3*nshell,displs=3*displs)
    else
       dnall = reshape(dnpart(:,1:nneigh_shell),(/3*nneigh_shell/))
    endif

    !
    !-- bucket the gathered entries by sector (the gathered data keep their order,
    !   so displs still gives the owner of entry j)
    !
    nsecall = 0
    do j=1,nshell_tot
       k = nint(dnall(3*j))
       nsecall(k) = nsecall(k) + 1
    enddo
    offsec(1) = 0
    do k=1,nsec
       offsec(k+1) = offsec(k) + nsecall(k)
    enddo
    nsecall = 0
    do j=1,nshell_tot
       k = nint(dnall(3*j))
       nsecall(k) = nsecall(k) + 1
       isort(offsec(k) + nsecall(k)) = j
    enddo

    !
    !-- every rank marches each sector independently from its own W(r_lo)
    !
    ifirst = displs(id+1) + 1
    ilast  = displs(id+1) + nshell(id+1)
    nhot   = 0
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(nsec,escaped,nsecall,offsec,isort,dnall,dr2sec,id_shell,dn_lo,r_lo,rs) &
    !$omp shared(ifirst,ilast,displs,id,listneigh,eos_vars,vxyzu) &
    !$omp private(k,nk,m,j,dnj,dntot,phi) &
    !$omp reduction(+:nhot)
    do k=1,nsec

       if (escaped(k)) cycle
       nk = nsecall(k)
       do m=1,nk
          dr2sec(offsec(k)+m) = dnall(3*isort(offsec(k)+m)-2)
       enddo

       call indexx(nk,dr2sec(offsec(k)+1:offsec(k)+nk),id_shell(offsec(k)+1:offsec(k)+nk))
       dntot = dn_lo(k)
       rs(k) = r_lo(k)

       do m=1,nk
          j   = isort(offsec(k) + id_shell(offsec(k)+m))
          dnj = dnall(3*j-1)
          if (dntot + dnj <= 1.) then
             phi  = 1.
             nhot = nhot + 1
          else
             phi = (1. - dntot)/dnj
          endif
          if (j >= ifirst .and. j <= ilast) &
             call ionise_particle(eos_vars,vxyzu,listneigh(j-displs(id+1)),phi)
          rs(k) = sqrt(dnall(3*j-2))
          dntot = dntot + dnj
          if (phi < 1.) exit
       enddo

    enddo
    !$omp end parallel do

    if (nsec == 1) then
       rstrom_sec(:,i) = rs(1)
    else
       rstrom_sec(:,i) = rs
    endif
    !-- stored value is the mean over sectors
    xyzmh_ptmass(irstrom,i) = sum(rs)/real(nsec)

    !
    !-- resolution of the bubble: sectors only once each holds a few hot particles
    !
    if (mpi) nhot_loc = int(reduceall_mpi('+',nhot_loc))
    nhot = nhot + nhot_loc
    if (sectors_on(i)) then
       if (nhot < nsec_off*nsec_H2R) sectors_on(i) = .false.
    else
       if (nhot >= nsec_on*nsec_H2R) sectors_on(i) = .true.
    endif

    if (iverbose >= 2) write(iprint,*) 'HII_bound: sink ',i,' Rstrom mean/min/max = ',&
                       sum(rs)/real(nsec),minval(rs),maxval(rs),' escaped sectors ',count(escaped),&
                       ' shell size ',nshell_tot,' in ',its,' iterations, ',nhot,' hot particles, ',&
                       nsec,' sectors'
    deallocate(dnall,dr2sec,id_shell,isort)
 enddo sources

 if (allocated(dnpart)) deallocate(dnpart)
 if (allocated(dn_lo)) deallocate(dn_lo,dn_hi,r_lo,r_hi,deltar,rs,converged,escaped,nsecall,offsec)
 deallocate(rprev2,rprevmax2)

end subroutine HII_bound

!-----------------------------------------------------------------------
!+
!  photon demand of source i inside r_lo and inside r_hi of each sector,
!  as a fraction of the sector's budget Q/nsec (eq. 6, 13). A particle that
!  other sources also reached at the previous update costs source i only
!  its flux share (eq. 14):
!    f_ij = F_ij / sum_s F_sj,  F_sj = Q_s/(4 pi d_sj^2),
!  the sum running over i and every source s with d_sj < r_S(s) last update,
!  r_S(s) taken in the sector of s that contains j.
!  listneigh is reordered into three blocks:
!    1 : nneigh_shell                          shell [r_lo,r_hi) of the particle's sector,
!                                              same order as dnpart(1:3,:) = (d^2, y_j, sector)
!    nneigh_shell+1 : nneigh_shell+nneigh_inner inside r_lo of the particle's sector
!    the rest                                   beyond r_hi, not gas, or dead
!+
!-----------------------------------------------------------------------
subroutine compute_dnspheres(i,nptmass,xyzmh_ptmass,rprev2,rprevmax2,coef,r_lo,r_hi,xyzh,rho,&
                             listneigh,nneigh,nneigh_shell,nneigh_inner,dn_lo,dn_hi,dnpart)
 use part,    only:get_pmass,iphase,get_partinfo,isdead_or_accreted,irateion
 use healpix, only:vec2pix_nest
 integer, intent(in)    :: i,nptmass,nneigh
 integer, intent(inout) :: listneigh(:)
 integer, intent(out)   :: nneigh_shell,nneigh_inner
 real,    intent(in)    :: xyzmh_ptmass(:,:),rprev2(:,:),rprevmax2(:),coef,r_lo(:),r_hi(:)
 real,    intent(in)    :: xyzh(:,:),rho(:)
 real,    intent(out)   :: dn_lo(:),dn_hi(:),dnpart(:,:)
 integer, allocatable :: iblock(:),listold(:)
 integer :: j,n,s,k,ks,ipix,itypej,ninner,nrej
 real    :: dx(3),dxs(3),dr2,dnj,d2s,fsum
 logical :: isactive,isgas,isdust

 dn_lo = 0.
 dn_hi = 0.
 !
 !-- per particle: sector, distance and demand into dnpart(:,n), and its block
 !   (0 rejected, 1 shell, 2 inner). This is the costly part: run it in parallel
 !
 allocate(iblock(nneigh))
 !$omp parallel do default(none) &
 !$omp shared(nneigh,listneigh,iphase,xyzh,xyzmh_ptmass,i,nptmass,rprev2,rprevmax2) &
 !$omp shared(coef,rho,r_lo,r_hi,nside_H2R,dnpart,iblock) &
 !$omp private(n,j,isactive,isgas,isdust,itypej,dx,dr2,ipix,k,fsum,s,dxs,d2s,ks,dnj)
 do n=1,nneigh

    iblock(n) = 0
    j = listneigh(n)
    call get_partinfo(iphase(j),isactive,isgas,isdust,itypej)
    if (.not.isgas .or. isdead_or_accreted(xyzh(4,j))) cycle

    dx  = xyzh(1:3,j) - xyzmh_ptmass(1:3,i)
    dr2 = dot_product(dx,dx)

    if (size(r_hi) > 1 .and. dr2 > 0.) then  ! size 1: isotropic, a single sector
       call vec2pix_nest(nside_H2R,dx,ipix)
       k = ipix + 1
    else
       k = 1
    endif

    if (dr2 >= r_hi(k)**2) cycle  ! the tree returns whole leaves, beyond r_hi too

    !-- sum_s F_sj / F_ij, with Q_s/Q_i = 10**(logQ_s - logQ_i) to keep numbers O(1)
    fsum = 1.
    do s=1,nptmass
       if (s == i .or. rprevmax2(s) < 0.) cycle  ! < 0: not a source, or no radius yet
       dxs = xyzh(1:3,j) - xyzmh_ptmass(1:3,s)
       d2s = dot_product(dxs,dxs)
       if (d2s >= rprevmax2(s)) cycle             ! beyond every sector of s
       if (d2s > 0.) then
          call vec2pix_nest(nside_H2R,dxs,ipix)
          ks = ipix + 1
       else
          ks = 1
       endif
       if (d2s < rprev2(ks,s)) fsum = fsum + 10.**(xyzmh_ptmass(irateion,s)-xyzmh_ptmass(irateion,i)) &
                                             *dr2/max(d2s,tiny(d2s))
    enddo

    dnj = coef*get_pmass(j,.true.)*rho(j)/fsum
    dnpart(1,n) = dr2
    dnpart(2,n) = dnj
    dnpart(3,n) = real(k)

    if (dr2 < r_lo(k)**2) then
       iblock(n) = 2
    else
       iblock(n) = 1
    endif

 enddo
 !$omp end parallel do
 !
 !-- sums and partition in serial, in list order: the sums do not depend on the
 !   number of threads. The shell keeps dnpart compacted in its own order (m <= n)
 !
 nneigh_shell = count(iblock == 1)
 nneigh_inner = count(iblock == 2)
 allocate(listold(nneigh))
 listold = listneigh(1:nneigh)
 ninner  = nneigh_shell
 nrej    = nneigh_shell + nneigh_inner
 nneigh_shell = 0

 do n=1,nneigh
    select case(iblock(n))
    case(1) !-- shell
       k = nint(dnpart(3,n))
       dn_hi(k) = dn_hi(k) + dnpart(2,n)
       nneigh_shell = nneigh_shell + 1
       listneigh(nneigh_shell) = listold(n)
       dnpart(:,nneigh_shell)  = dnpart(:,n)
    case(2) !-- inner
       k = nint(dnpart(3,n))
       dn_hi(k) = dn_hi(k) + dnpart(2,n)
       dn_lo(k) = dn_lo(k) + dnpart(2,n)
       ninner = ninner + 1
       listneigh(ninner) = listold(n)
    case default !--rejected
       nrej = nrej + 1
       listneigh(nrej) = listold(n)
    end select
 enddo
 deallocate(iblock,listold)

end subroutine compute_dnspheres

!-----------------------------------------------------------------------
!+
!  first guess of the Stromgren radius of a new source (eq. 8), from the
!  density of the nearest gas particle over all ranks; rbound if no gas
!+
!-----------------------------------------------------------------------
real function stromgren_guess(xi,yi,zi,hacc,log_Qi,xyzh,rho,rbound) result(rs)
 use neighkdtree, only:listneigh,getneigh_pos,leaf_is_active
 use mpiutils,    only:reduceall_mpi
 use dim,         only:mpi
 use units,       only:utime
 use physcon,     only:pi
 use part,        only:iphase,get_partinfo,isdead_or_accreted
 real, intent(in) :: xi,yi,zi,hacc,log_Qi,rbound
 real, intent(in) :: xyzh(:,:),rho(:)
 integer, parameter :: maxc = 128
 real    :: xyzcache(3,maxc)
 integer :: k,j,nneigh,itypej
 real    :: rq,dr2,d2loc,d2min,rholoc,rhonear,nH,Q
 logical :: isactive,isgas,isdust

 rs = rbound
 rq = max(hacc,1.e-3*rbound)
 do while (rq < rbound)
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
 use part,         only:imu
 integer, intent(in)    :: j
 real,    intent(in)    :: phi
 real,    intent(inout) :: eos_vars(:,:),vxyzu(:,:)
 real :: phi_old

 phi_old = (gmw - eos_vars(imu,j))/(gmw - muion)
 if (phi <= phi_old) return

 eos_vars(imu,j) = phi*muion + (1. - phi)*gmw
 !-- heated at once for the forces (vpred in derivs); the step relaxes u to the same target
 vxyzu(4,j) = heated_u_HIIR(vxyzu(4,j),eos_vars(imu,j))

end subroutine ionise_particle

!-----------------------------------------------------------------------
!+
! HII feedback routine using K nearest neighbors prescription
!+
!-----------------------------------------------------------------------
subroutine HII_Knn(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,rho,eos_vars)
 use part,       only:massoftype,ihsoft,igas,irateion,isdead_or_accreted,&
                      irstrom,get_partinfo,iphase,imu,get_pmass
 use neighkdtree,   only:listneigh=>listneigh_global,getneigh_pos,leaf_is_active
 use sortutils,  only:Knnfunc,set_r2func_origin,r2func_origin
 use physcon,    only:pc,pi
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
                            call ionise_particle(eos_vars,vxyzu,j,1.)
                            npartin = npartin + 1
                         endif
                      else
                         if (k > 1) then ! end of the HII region
                            DNratio = 10**(Ndot)/10**(DNdot)
                            call ionise_particle(eos_vars,vxyzu,j,DNratio)
                            Ndot = 0.
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
                    iphase,get_partinfo,noverlap,imu
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
!$omp shared(pmass,iphase,iH2R,uIon,gmw)&
!$omp private(log_Qi,i,j,xj,yj,zj,fluxi,itypei,rsrci,isgas,isdust)&
!$omp private(isactive,noverlapi)&
!$omp reduction(+:k)
!$omp do
    do i=1,npart
       if (.not. isdead_or_accreted(xyzh(4,i))) then
          call get_partinfo(iphase(i),isactive,isgas,isdust,itypei)
          ! if (.not. isactive) cycle
          if (isgas) then
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
                   call ionise_particle(eos_vars,vxyzu,i,1.)
                   noverlapi = noverlapi + 1
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
!  relaxation of the thermal energy to the target of its phase over dt,
!  exact for any dt:
!    u_new = u_target + (u - u_target)*exp(-dt/t_relax),  t_relax = Crelax*h/c_s(u)
!  uIon for the fully ionised gas, u_target at once below the target.
!  It does not commute with the kicks, so the step applies it once per
!  interval of time of the particle, not in each kick
!+
!-----------------------------------------------------------------------
real function relax_u_HIIR(ui,hi,mui,dt) result(unew)
 use eos, only:gamma
 real, intent(in) :: ui,hi,mui,dt
 real :: utarget,trelax

 unew = ui
 if (dt <= 0.) return

 if (is_fully_ionised(mui)) then
    unew = uIon
 else
    utarget = get_utarget(mui)
    if (unew <= utarget) then
       unew = utarget  ! heated by the stars
    else
       trelax = Crelax*hi/sqrt(gamma*(gamma - 1.)*unew)
       if (trelax > 0.) then
          unew = utarget + (unew - utarget)*exp(-dt/trelax)
          if (unew - utarget < tolphase*utarget) unew = utarget  ! back in its phase
       else
          unew = utarget  ! Crelax = 0: two isothermal phases
       endif
    endif
 endif

end function relax_u_HIIR

!-----------------------------------------------------------------------
!+
!  thermal energy given to the forces when a particle is ionised:
!  uIon for the fully ionised gas, at least its target otherwise
!+
!-----------------------------------------------------------------------
real function heated_u_HIIR(ui,mui) result(uheated)
 real, intent(in) :: ui,mui

 if (is_fully_ionised(mui)) then
    uheated = uIon
 else
    uheated = max(ui,get_utarget(mui))
 endif

end function heated_u_HIIR

!-----------------------------------------------------------------------
!+
!  gas (partially) ionised by the HII routines
!+
!-----------------------------------------------------------------------
logical function is_ionised(mui)
 real, intent(in) :: mui

 is_ionised = (mui > 0. .and. mui < gmw)

end function is_ionised

logical function is_fully_ionised(mui)
 real, intent(in) :: mui

 is_fully_ionised = (mui > 0. .and. mui <= muion)

end function is_fully_ionised

!-----------------------------------------------------------------------
!+
!  thermal energy of the phase, T_target = phi*Tion + (1-phi)*Tcold
!+
!-----------------------------------------------------------------------
real function get_utarget(mui) result(utarget)
 use eos,      only:gamma,temperature_coef
 use eos_HIIR, only:Tion,Tcold
 real, intent(in) :: mui
 real :: mu,phi

 mu = mui
 if (mu <= 0.) mu = gmw
 phi = min(max((gmw - mu)/(gmw - muion),0.),1.)
 utarget = (phi*Tion + (1. - phi)*Tcold)/(temperature_coef*mu*(gamma - 1.))

end function get_utarget

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
    call write_inopt(Crelax, 'Crelax', "relaxation time of the gas of the HII fronts in units of h/c_s", iunit)
    if (iH2R == 4) call write_inopt(nside_H2R, 'nsideH2R', &
       "HEALPix Nside of the sectors, 12*Nside**2 sectors (power of 2, only bounded Knn)", iunit)
 endif

end subroutine write_options_H2R

!-----------------------------------------------------------------------
!+
!  read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_H2R(db,nerr)
 use infile_utils, only:inopts,read_inopt
 use io,           only:error
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(iH2R,'iH2R',db,errcount=nerr,min=0,default=0)
 if (iH2R > 0) then
    call read_inopt(Mmin,'Mmin',db,errcount=nerr,min=8.)
    call read_inopt(Rmax,'Rmax',db,errcount=nerr,min=10.)
    call read_inopt(Crelax,'Crelax',db,errcount=nerr,min=0.,default=Crelax)
    if (iH2R == 4) then
       call read_inopt(nside_H2R,'nsideH2R',db,errcount=nerr,min=1,default=1)
       if (iand(nside_H2R,nside_H2R-1) /= 0) then
          call error('read_options_H2R','nsideH2R must be a power of 2')
          nerr = nerr + 1
       endif
    endif
 endif

end subroutine read_options_H2R

end module HIIRegion
