!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module energies
!
! Calculates global quantities on the particles
!   To Developer: See instructions in evwrite.F90 about adding values
!                 to the .ev file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary_dyn, centreofmass, dim, dust, eos, eos_piecewise,
!   externalforces, gravwaveutils, io, kernel, metric_tools, mpiutils,
!   nicil, options, part, ptmass, subgroup, timestep, units, utils_gr,
!   vectorutils, viscosity
!
 use dim,   only:maxdusttypes,maxdustsmall
 use units, only:utime
 implicit none

 logical,         public    :: gas_only,track_mass,track_lum
 real,            public    :: ekin,etherm,emag,epot,etot,eacc,totmom,angtot,mtot,xyzcom(3)
 real,            public    :: ekinacc,ethermacc,emagacc,epotacc,eradacc,etotall
 real,            public    :: hdivBonB_ave,hdivBonB_max
 real,            public    :: vrms,rmsmach,accretedmass,mdust(maxdusttypes),mgas
 real,            public    :: xcom,ycom,zcom,xmom,ymom,zmom,angx,angy,angz
 real,            public    :: totlum,angxall,angyall,angzall,angall
 real,            public    :: hx(4),hp(4),ddq_xy(3,3)
 integer,         public    :: iquantities
 integer(kind=8), public    :: ndead,npartall,np_cs_eq_0,np_e_eq_0
 integer,         public    :: iev_time,iev_ekin,iev_etherm,iev_emag,iev_epot,iev_etot,iev_totmom,iev_com(3),&
                               iev_angmom,iev_rho,iev_dt,iev_dtx,iev_entrop,iev_rmsmach,iev_vrms,iev_rhop(6),&
                               iev_alpha,iev_B,iev_divB,iev_hdivB,iev_beta,iev_temp,iev_etao,iev_etah(2),&
                               iev_etaa,iev_vel,iev_vhall,iev_vion,iev_n(7),&
                               iev_dtg,iev_ts,iev_dm(maxdusttypes),iev_momall,iev_angall,iev_maccsink(2),&
                               iev_macc,iev_eacc,iev_totlum,iev_erot(4),iev_viscrat,iev_gws(8),&
                               iev_mass,iev_bdy(3,2)
 integer,         public    :: iev_erad
 real,            public    :: erad
 integer,         parameter :: inumev  = 150  ! maximum number of quantities to be printed in .ev
 integer,         parameter :: iev_sum = 1    ! array index of the sum of the quantity
 integer,         parameter :: iev_max = 2    ! array index of the maximum of the quantity
 integer,         parameter :: iev_min = 3    ! array index of the minimum of the quantity
 integer,         parameter :: iev_ave = 4    ! array index of the average of the quantity
 ! Subroutines
 public  :: compute_energies,ev_data_update
 private :: get_erot,initialise_ev_data,collate_ev_data,finalise_ev_data
 ! Arrays
 real,             public :: ev_data(4,0:inumev),erot_com(6)

contains

!----------------------------------------------------------------
!+
!  Subroutine to compute global quantities on the particles
!+
!----------------------------------------------------------------
subroutine compute_energies(t)
 use dim,            only:maxp,maxvxyzu,maxalpha,maxtypes,mhd_nonideal,maxp_hard,&
                          lightcurve,use_dust,maxdusttypes,do_radiation,gr,use_krome,&
                          use_apr,use_sinktree,maxpsph
 use part,           only:rhoh,xyzh,vxyzu,massoftype,npart,maxphase,iphase,&
                          alphaind,Bevol,divcurlB,iamtype,igamma,&
                          igas,idust,iboundary,istar,idarkmatter,ibulge,&
                          nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,isdeadh,&
                          isdead_or_accreted,epot_sinksink,imacc,ispinx,ispiny,&
                          ispinz,mhd,gravity,poten,dustfrac,eos_vars,itemp,igasP,ics,&
                          nden_nimhd,eta_nimhd,iion,ndustsmall,graindens,grainsize,&
                          iamdust,ndusttypes,rad,iradxi,gtgrad,group_info,bin_info,n_group
 use part,           only:pxyzu,fxyzu,fext,apr_level,aprmassoftype,pxyzu_ptmass
 use gravwaveutils,  only:calculate_strain,calc_gravitwaves
 use centreofmass,   only:get_centreofmass_accel
 use eos,            only:polyk,gamma,eos_is_non_ideal,eos_outputs_gasP
 use eos_piecewise,  only:gamma_pwp
 use io,             only:id,fatal,master
 use externalforces, only:externalforce,externalforce_vdependent,was_accreted,accradius1
 use options,        only:iexternalforce,calc_erot,alpha,ieos,use_dustfrac
 use mpiutils,       only:reduceall_mpi
 use ptmass,         only:get_accel_sink_gas,use_regnbody
 use subgroup,       only:get_pot_subsys
 use viscosity,      only:irealvisc,shearfunc
 use nicil,          only:nicil_update_nimhd,nicil_get_halldrift,nicil_get_ambidrift, &
                     use_ohm,use_hall,use_ambi,n_data_out,n_warn,eta_constant
 use boundary_dyn,   only:dynamic_bdy,find_dynamic_boundaries
 use kernel,         only:radkern
 use timestep,       only:dtmax
 use part,           only:metrics,metrics_ptmass
 use metric_tools,   only:unpack_metric
 use utils_gr,       only:dot_product_gr,get_geodesic_accel
 use vectorutils,    only:cross_product3D
 use part,           only:luminosity
 use dust,           only:get_ts,idrag
 real, intent(in) :: t
 integer :: iregime,idusttype,ierr
 real    :: ev_data_thread(4,0:inumev)
 real    :: xi,yi,zi,hi,vxi,vyi,vzi,v2i,Bxi,Byi,Bzi,Bi,B2i,rhoi
 real    :: xmomacc,ymomacc,zmomacc,angaccx,angaccy,angaccz,dm
 real    :: epoti,pmassi,dnptot,dnpgas,tsi
 real    :: xmomall,ymomall,zmomall,rho1i,vsigi
 real    :: ponrhoi,spsoundi,dumx,dumy,dumz,gammai
 real    :: divBi,hdivBonBi,alphai,valfven2i,betai
 real    :: n_total,n_total1,n_ion,shearparam_art,shearparam_phys,ratio_phys_to_av
 real    :: gasfrac,rhogasi,dustfracisum,dustfraci(maxdusttypes),dust_to_gas(maxdusttypes)
 real    :: etaohm,etahall,etaambi,vhall,vion
 real    :: curlBi(3),vhalli(3),vioni(3),data_out(n_data_out)
 real    :: erotxi,erotyi,erotzi,fdum(3),x0(3),v0(3),a0(3),xyz_x_all(3),xyz_n_all(3)
 real    :: ekini,ethermi,epottmpi,eradi,emagi
 real    :: pdotv,bigvi(1:3),alpha_gr,beta_gr_UP(1:3),lorentzi,pxi,pyi,pzi
 real    :: gammaijdown(1:3,1:3),angi(1:3),fourvel_space(3)
 integer :: i,j,itype,iu
 integer :: ierrlist(n_warn)
 integer(kind=8) :: np,npgas,nptot,np_rho(maxtypes),np_rho_thread(maxtypes)
 logical :: was_not_accreted

 ! initialise values
 itype  = igas
 pmassi = massoftype(igas)
 ekin   = 0.
 etherm = 0.
 if (maxvxyzu < 4 .and. gamma < 1.0001 .and. ieos/=9) etherm = 1.5*polyk
 epot = 0.
 emag = 0.
 etot = 0.
 erad = 0.
 xcom = 0.
 ycom = 0.
 zcom = 0.
 mtot = 0.
 dm   = 0.
 xmom = 0.
 ymom = 0.
 zmom = 0.
 angx = 0.
 angy = 0.
 angz = 0.
 iu   = 4
 np   = 0
 vrms = 0.
 rmsmach = 0.
 npgas   = 0
 xmomacc = 0.
 ymomacc = 0.
 zmomacc = 0.
 angaccx = 0.
 angaccy = 0.
 angaccz = 0.
 ekinacc = 0.
 ethermacc = 0.
 emagacc = 0.
 epotacc = 0.
 eradacc = 0.
 mgas    = 0.
 mdust   = 0.
 mgas    = 0.
 np_cs_eq_0 = 0
 np_e_eq_0  = 0
 ierrlist = 0
 if (maxalpha==maxp) then
    alphai = 0.
 else
    alphai = alpha
 endif
 np_rho      = 0
 call initialise_ev_data(ev_data)

!$omp parallel default(none) &
!$omp shared(maxp,maxphase,maxalpha,maxpsph) &
!$omp shared(xyzh,vxyzu,pxyzu,rad,iexternalforce,npart,t,id) &
!$omp shared(alphaind,massoftype,irealvisc,iu,aprmassoftype) &
!$omp shared(ieos,gamma,nptmass,xyzmh_ptmass,vxyz_ptmass,xyzcom) &
!$omp shared(Bevol,divcurlB,iphase,poten,dustfrac,use_dustfrac) &
!$omp shared(use_ohm,use_hall,use_ambi,nden_nimhd,eta_nimhd,eta_constant) &
!$omp shared(ev_data,np_rho,erot_com,calc_erot,gas_only,track_mass) &
!$omp shared(calc_gravitwaves) &
!$omp shared(iev_erad,iev_rho,iev_dt,iev_entrop,iev_rhop,iev_alpha) &
!$omp shared(iev_B,iev_divB,iev_hdivB,iev_beta,iev_temp,iev_etao,iev_etah) &
!$omp shared(iev_etaa,iev_vel,iev_vhall,iev_vion,iev_n) &
!$omp shared(iev_dtg,iev_ts,iev_macc,iev_totlum,iev_erot,iev_viscrat) &
!$omp shared(eos_vars,grainsize,graindens,ndustsmall,metrics,metrics_ptmass,pxyzu_ptmass) &
!$omp private(i,j,xi,yi,zi,hi,rhoi,vxi,vyi,vzi,Bxi,Byi,Bzi,Bi,B2i,epoti,vsigi,v2i) &
!$omp private(ponrhoi,spsoundi,gammai,dumx,dumy,dumz,valfven2i,divBi,hdivBonBi,curlBi) &
!$omp private(rho1i,shearparam_art,shearparam_phys,ratio_phys_to_av,betai) &
!$omp private(gasfrac,rhogasi,dustfracisum,dustfraci,dust_to_gas,n_total,n_total1,n_ion) &
!$omp private(etaohm,etahall,etaambi,vhalli,vhall,vioni,vion,data_out) &
!$omp private(ekini,ethermi,emagi,eradi,epottmpi) &
!$omp private(erotxi,erotyi,erotzi,fdum) &
!$omp private(ev_data_thread,np_rho_thread) &
!$omp firstprivate(alphai,itype,pmassi) &
!$omp private(pxi,pyi,pzi,gammaijdown,alpha_gr,beta_gr_UP,bigvi,lorentzi,pdotv,angi,fourvel_space) &
!$omp shared(idrag) &
!$omp private(tsi,iregime,idusttype,was_not_accreted) &
!$omp shared(luminosity,track_lum,apr_level) &
!$omp reduction(+:np,npgas,np_cs_eq_0,np_e_eq_0) &
!$omp reduction(+:xcom,ycom,zcom,mtot,xmom,ymom,zmom,angx,angy,angz,mdust,mgas) &
!$omp reduction(+:xmomacc,ymomacc,zmomacc,angaccx,angaccy,angaccz) &
!$omp reduction(+:ekinacc,ethermacc,emagacc,epotacc,eradacc) &
!$omp reduction(+:ekin,etherm,emag,epot,erad,vrms,rmsmach,ierrlist)
 call initialise_ev_data(ev_data_thread)
 np_rho_thread  = 0
!$omp do
 do i=1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    was_not_accreted = .not.was_accreted(iexternalforce,hi)
    if (.not.isdead_or_accreted(hi) .or. .not. was_not_accreted) then
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
          if (itype <= 0) call fatal('energies','particle type <= 0')
          if (use_apr) then
             pmassi = aprmassoftype(itype,apr_level(i))
          else
             pmassi = massoftype(itype)
          endif
       else
          if (use_apr) then
             pmassi = aprmassoftype(igas,apr_level(i))
          else
             pmassi = massoftype(igas)
          endif
       endif

       rhoi = rhoh(hi,pmassi)
       if (was_not_accreted) then
          call ev_data_update(ev_data_thread,iev_rho,rhoi)
          if (.not.gas_only) then
             select case(itype)
             case(igas)
                call ev_data_update(ev_data_thread,iev_rhop(1), rhoi)
                np_rho_thread(igas) =  np_rho_thread(igas) + 1
             case(idust)
                call ev_data_update(ev_data_thread,iev_rhop(2),rhoi)
                np_rho_thread(idust) =  np_rho_thread(idust) + 1
             case(iboundary)
                call ev_data_update(ev_data_thread,iev_rhop(3), rhoi)
                np_rho_thread(iboundary) =  np_rho_thread(iboundary) + 1
             case(istar)
                call ev_data_update(ev_data_thread,iev_rhop(4),rhoi)
                np_rho_thread(istar) =  np_rho_thread(istar) + 1
             case(idarkmatter)
                call ev_data_update(ev_data_thread,iev_rhop(5),  rhoi)
                np_rho_thread(idarkmatter) =  np_rho_thread(idarkmatter) + 1
             case(ibulge)
                call ev_data_update(ev_data_thread,iev_rhop(6), rhoi)
                np_rho_thread(ibulge) =  np_rho_thread(ibulge) + 1
             end select
          endif
          np = np + 1
       endif

       vxi  = vxyzu(1,i)
       vyi  = vxyzu(2,i)
       vzi  = vxyzu(3,i)

       if (gr) then
          pxi  = pxyzu(1,i)
          pyi  = pxyzu(2,i)
          pzi  = pxyzu(3,i)

          call unpack_metric(metrics(:,:,:,i),betaUP=beta_gr_UP,alpha=alpha_gr,gammaijdown=gammaijdown)
          bigvi    = (vxyzu(1:3,i)+beta_gr_UP)/alpha_gr
          v2i      = dot_product_gr(bigvi,bigvi,gammaijdown)
          lorentzi = 1./sqrt(1.-v2i)
          pdotv    = pxi*vxi + pyi*vyi + pzi*vzi

          ! angular momentum
          fourvel_space = (lorentzi/alpha_gr)*vxyzu(1:3,i)
          call cross_product3D(xyzh(1:3,i),fourvel_space,angi) ! position cross with four-velocity

          ! kinetic energy
          ekini = pmassi*(pdotv + alpha_gr/lorentzi - 1.) ! The 'kinetic term' in total specific energy, minus rest mass
       else
          pxi = vxi
          pyi = vyi
          pzi = vzi

          ! centre of mass
          xcom = xcom + pmassi*xi
          ycom = ycom + pmassi*yi
          zcom = zcom + pmassi*zi

          ! angular momentum
          angi(1) = (yi*vzi - zi*vyi)
          angi(2) = (zi*vxi - xi*vzi)
          angi(3) = (xi*vyi - yi*vxi)

          ! kinetic energy and rms velocity
          v2i   = vxi*vxi + vyi*vyi + vzi*vzi
          ekini = pmassi*v2i
       endif

       if (was_not_accreted) then
          ! total mass
          mtot = mtot + pmassi

          ! linear momentum
          xmom = xmom + pmassi*pxi
          ymom = ymom + pmassi*pyi
          zmom = zmom + pmassi*pzi

          ! angular momentum
          angx = angx + pmassi*angi(1)
          angy = angy + pmassi*angi(2)
          angz = angz + pmassi*angi(3)

          ! kinetic energy & rms velocity
          ekin = ekin + ekini
          vrms = vrms + v2i
       else
          call ev_data_update(ev_data_thread,iev_macc,pmassi)

          ! linear momentum (accreted particles)
          xmomacc = xmomacc + pmassi*pxi
          ymomacc = ymomacc + pmassi*pyi
          zmomacc = zmomacc + pmassi*pzi

          ! angular momentum (accreted particles)
          angaccx = angaccx + pmassi*angi(1)
          angaccy = angaccy + pmassi*angi(2)
          angaccz = angaccz + pmassi*angi(3)

          ! kinetic energy (accreted particles
          ekinacc = ekinacc + ekini
       endif

       ! rotational energy around each axis through the Centre of mass
       ! note: for efficiency, centre of mass is from the previous time energies was called
       if (calc_erot .and. was_not_accreted) then
          call get_erot(xi,yi,zi,vxi,vyi,vzi,xyzcom,pmassi,erotxi,erotyi,erotzi)
          call ev_data_update(ev_data_thread,iev_erot(1),erotxi)
          call ev_data_update(ev_data_thread,iev_erot(2),erotyi)
          call ev_data_update(ev_data_thread,iev_erot(3),erotzi)
       endif

       ! potential energy
       epoti = 0.
       if (iexternalforce > 0 .and. .not.gr) then
          dumx = 0.
          dumy = 0.
          dumz = 0.
          epottmpi = 0.
          call externalforce(iexternalforce,xi,yi,zi,hi,t,dumx,dumy,dumz,epottmpi,ii=i)
          call externalforce_vdependent(iexternalforce,xyzh(1:3,i),vxyzu(1:3,i),fdum,epottmpi)
          epoti = pmassi*epottmpi
       endif
       if (nptmass > 0 .and. .not.use_sinktree) then ! No need to compute if sink in tree
          dumx = 0.
          dumy = 0.
          dumz = 0.
          epottmpi = 0.
          call get_accel_sink_gas(nptmass,xi,yi,zi,hi,xyzmh_ptmass,dumx,dumy,dumz,epottmpi)
          epoti = epoti + pmassi*epottmpi
       endif
       if (gravity) epoti = epoti + poten(i)
       if (was_not_accreted) then
          epot = epot + epoti
       else
          epotacc = epotacc + epoti
       endif

       !
       ! total dust mass for each species
       !
       if (use_dust .and. was_not_accreted) then
          if (iamdust(iphase(i))) then
             idusttype = ndustsmall + itype - idust + 1
             mdust(idusttype) = mdust(idusttype) + pmassi
          endif
       endif

       if (do_radiation) then
          eradi = pmassi*rad(iradxi,i)
          if (was_not_accreted) then
             erad = erad + eradi
          else
             eradacc = eradacc + eradi
          endif
       endif

       !
       ! the following apply ONLY to gas particles
       !
       isgas: if (itype==igas) then

          if (use_dustfrac) then
             dustfraci    = dustfrac(:,i)
             dustfracisum = sum(dustfraci)
             gasfrac      = 1. - dustfracisum
             dust_to_gas  = dustfraci(:)/gasfrac
             if (was_not_accreted) then
                do j=1,ndustsmall
                   call ev_data_update(ev_data_thread,iev_dtg,dust_to_gas(j))
                enddo
                mdust(1:ndustsmall) = mdust(1:ndustsmall) + pmassi*dustfraci(1:ndustsmall)
             endif
          else
             dustfraci    = 0.
             dustfracisum = 0.
             gasfrac      = 1.
          endif
          if (was_not_accreted) then
             npgas = npgas + 1
             mgas = mgas + pmassi*gasfrac
          endif

          ! thermal energy
          ponrhoi  = eos_vars(igasP,i)/rhoi
          spsoundi = eos_vars(ics,i)
          gammai   = eos_vars(igamma,i)
          if (maxvxyzu >= 4) then
             ethermi = pmassi*vxyzu(4,i)*gasfrac
             if (gr) ethermi = (alpha_gr/lorentzi)*ethermi

             if (was_not_accreted) then
                if (vxyzu(iu,i) < tiny(vxyzu(iu,i))) np_e_eq_0 = np_e_eq_0 + 1
                if (spsoundi < tiny(spsoundi) .and. vxyzu(iu,i) > 0. ) np_cs_eq_0 = np_cs_eq_0 + 1
             endif
          else
             if ((ieos==2 .or. ieos == 5 .or. ieos == 17) .and. gammai > 1.001) then
                !--thermal energy using polytropic equation of state
                ethermi = pmassi*ponrhoi/(gammai-1.)*gasfrac
             elseif (ieos==9) then
                !--thermal energy using piecewise polytropic equation of state
                ethermi = pmassi*ponrhoi/(gamma_pwp(rhoi)-1.)*gasfrac
             else
                ethermi = 0.
             endif
             if (spsoundi < tiny(spsoundi) .and. was_not_accreted) np_cs_eq_0 = np_cs_eq_0 + 1
          endif
          vsigi = spsoundi

          if (was_not_accreted) then
             etherm = etherm + ethermi
          else
             ethermacc = ethermacc + ethermi
          endif

          if (was_not_accreted) then
             ! entropy
             call ev_data_update(ev_data_thread,iev_entrop,pmassi*ponrhoi*rhoi**(1.-gammai))

             ! gas temperature
             if (eos_is_non_ideal(ieos) .or. eos_outputs_gasP(ieos)) then
                call ev_data_update(ev_data_thread,iev_temp,eos_vars(itemp,i))
             endif

             ! min and mean stopping time
             if (use_dustfrac) then
                rhogasi = rhoi*gasfrac
                do j=1,ndustsmall
                   call get_ts(idrag,j,grainsize(j),graindens(j),rhogasi,rhoi*dustfracisum,spsoundi,0.,tsi,iregime)
                   call ev_data_update(ev_data_thread,iev_ts,tsi)
                enddo
             endif

             if (track_lum .and. lightcurve) call ev_data_update(ev_data_thread,iev_totlum,real(luminosity(i)))

             ! rms mach number
             if (spsoundi > 0.) rmsmach = rmsmach + v2i/spsoundi**2

             ! max of dissipation parameters
             if (maxalpha==maxp) then
                alphai = alphaind(1,i)
                call ev_data_update(ev_data_thread,iev_alpha,alphai)
             endif

             ! physical viscosity
             if (irealvisc /= 0) then
                shearparam_art  = 0.1*alphai*hi*vsigi
                shearparam_phys = shearfunc(xi,yi,zi,spsoundi)
                if (shearparam_art > 0.) then
                   ratio_phys_to_av = shearparam_phys/shearparam_art
                else
                   ratio_phys_to_av = 0.
                endif
                call ev_data_update(ev_data_thread,iev_viscrat,ratio_phys_to_av)
             endif
          endif

          ! mhd parameters
          if (mhd) then
             Bxi       = Bevol(1,i)*rhoi
             Byi       = Bevol(2,i)*rhoi
             Bzi       = Bevol(3,i)*rhoi
             B2i       = Bxi*Bxi + Byi*Byi + Bzi*Bzi
             Bi        = sqrt(B2i)
             rho1i     = 1./rhoi
             valfven2i = B2i*rho1i
             vsigi     = sqrt(valfven2i + spsoundi*spsoundi)
             emagi     = pmassi*B2i*rho1i

             if (was_not_accreted) then
                emag      = emag + emagi
                divBi     = abs(divcurlB(1,i))
                if (B2i > 0.) then
                   hdivBonBi = hi*divBi/Bi
                   betai     = 2.0*ponrhoi*rhoi/B2i ! plasma beta
                else
                   hdivBonBi = 0.
                   betai     = 0.
                endif
                call ev_data_update(ev_data_thread,iev_B,    Bi       )
                call ev_data_update(ev_data_thread,iev_divB, divBi    )
                call ev_data_update(ev_data_thread,iev_hdivB,hdivBonBi)
                call ev_data_update(ev_data_thread,iev_beta, betai    )

                if ( mhd_nonideal ) then
                   call nicil_update_nimhd(0,etaohm,etahall,etaambi,Bi,rhoi, &
                                        eos_vars(itemp,i),nden_nimhd(:,i),ierrlist,data_out)
                   curlBi = divcurlB(2:4,i)
                   if (use_ohm) then
                      call ev_data_update(ev_data_thread,iev_etao,   etaohm      )
                   endif
                   if (use_hall) then
                      call nicil_get_halldrift(etahall,Bxi,Byi,Bzi,curlBi,vhalli)
                      vhall = sqrt( dot_product(vhalli,vhalli) )
                      call ev_data_update(ev_data_thread,iev_etah(1),etahall     )
                      call ev_data_update(ev_data_thread,iev_etah(2),abs(etahall))
                      call ev_data_update(ev_data_thread,iev_vhall  ,vhall       )
                   endif
                   if (use_ambi) then
                      call nicil_get_ambidrift(etaambi,Bxi,Byi,Bzi,curlBi,vioni)
                      vion   = sqrt( dot_product(vioni,  vioni  ) )
                      call ev_data_update(ev_data_thread,iev_etaa,   etaambi     )
                      call ev_data_update(ev_data_thread,iev_vel,    sqrt(v2i)   )
                      call ev_data_update(ev_data_thread,iev_vion,   vion        )
                   endif
                   if (.not.eta_constant) then
                      n_ion = 0
                      do j = 9,21
                         n_ion = n_ion + data_out(j)
                      enddo
                      n_total = data_out(5)
                      if (n_total > 0.) then
                         n_total1 = 1.0/n_total
                      else
                         n_total1 = 0.0         ! only possible if eta_constant = .true.
                      endif
                      eta_nimhd(iion,i) = n_ion*n_total1    ! Save ionisation fraction for the dump file
                      call ev_data_update(ev_data_thread,iev_n(1),n_ion*n_total1)
                      call ev_data_update(ev_data_thread,iev_n(2),data_out( 8)*n_total1)
                      call ev_data_update(ev_data_thread,iev_n(3),data_out( 8))
                      call ev_data_update(ev_data_thread,iev_n(4),n_total-n_ion)
                      call ev_data_update(ev_data_thread,iev_n(5),data_out(24))
                      call ev_data_update(ev_data_thread,iev_n(6),data_out(23))
                      call ev_data_update(ev_data_thread,iev_n(7),data_out(22))
                   endif
                endif
             else
                emagacc = emagacc + emagi
             endif
          endif
       endif isgas
    endif
 enddo
!$omp enddo
!
!--add contribution from sink particles
!
 if (id==master) then

    if (.not. gr) then
       !$omp do
       do i=1,nptmass
          xi     = xyzmh_ptmass(1,i)
          yi     = xyzmh_ptmass(2,i)
          zi     = xyzmh_ptmass(3,i)
          pmassi = xyzmh_ptmass(4,i)
          if (pmassi < 0.) cycle

          vxi    = vxyz_ptmass(1,i)
          vyi    = vxyz_ptmass(2,i)
          vzi    = vxyz_ptmass(3,i)

          !phii   = fxyz_ptmass(4,i)

          xcom = xcom + pmassi*xi
          ycom = ycom + pmassi*yi
          zcom = zcom + pmassi*zi
          mtot = mtot + pmassi

          xmom   = xmom + pmassi*vxi
          ymom   = ymom + pmassi*vyi
          zmom   = zmom + pmassi*vzi

          angx   = angx + pmassi*(yi*vzi - zi*vyi)
          angy   = angy + pmassi*(zi*vxi - xi*vzi)
          angz   = angz + pmassi*(xi*vyi - yi*vxi)

          if (use_sinktree) epot = epot + poten(i+maxpsph)
          angx   = angx + xyzmh_ptmass(ispinx,i)
          angy   = angy + xyzmh_ptmass(ispiny,i)
          angz   = angz + xyzmh_ptmass(ispinz,i)

          v2i    = vxi*vxi + vyi*vyi + vzi*vzi
          ekin   = ekin + pmassi*v2i


          ! rotational energy around each axis through the origin
          if (calc_erot) then
             call get_erot(xi,yi,zi,vxi,vyi,vzi,xyzcom,pmassi,erotxi,erotyi,erotzi)
             call ev_data_update(ev_data_thread,iev_erot(1),erotxi)
             call ev_data_update(ev_data_thread,iev_erot(2),erotyi)
             call ev_data_update(ev_data_thread,iev_erot(3),erotzi)
          endif
       enddo
       !$omp enddo
    else
       !$omp do
       do i=1,nptmass
          ! calculate Kinetic and thermal energy for the GR-sink case.
          xi     = xyzmh_ptmass(1,i)
          yi     = xyzmh_ptmass(2,i)
          zi     = xyzmh_ptmass(3,i)
          pmassi = xyzmh_ptmass(4,i)
          if (pmassi < 0.) cycle

          vxi    = vxyz_ptmass(1,i)
          vyi    = vxyz_ptmass(2,i)
          vzi    = vxyz_ptmass(3,i)

          pxi    = pxyzu_ptmass(1,i)
          pyi    = pxyzu_ptmass(2,i)
          pzi    = pxyzu_ptmass(3,i)

          mtot = mtot + pmassi

          call unpack_metric(metrics_ptmass(:,:,:,i),betaUP=beta_gr_UP,alpha=alpha_gr,gammaijdown=gammaijdown)
          bigvi    = (vxyz_ptmass(1:3,i)+beta_gr_UP)/alpha_gr
          v2i      = dot_product_gr(bigvi,bigvi,gammaijdown)
          lorentzi = 1./sqrt(1.-v2i)
          pdotv    = pxi*vxi + pyi*vyi + pzi*vzi

          ! angular momentum
          fourvel_space = (lorentzi/alpha_gr)*vxyz_ptmass(1:3,i)
          call cross_product3D(xyzmh_ptmass(1:3,i),fourvel_space,angi) ! position cross with four-velocity

          ! kinetic energy
          ekini = pmassi*(pdotv + alpha_gr/lorentzi - 1.) ! The 'kinetic term' in total specific energy, minus rest mass

          ! kinetic energy & rms velocity
          ekin = ekin + ekini
          vrms = vrms + v2i

          ! linear momentum
          xmom = xmom + pmassi*pxi
          ymom = ymom + pmassi*pyi
          zmom = zmom + pmassi*pzi

          ! angular momentum
          angx = angx + pmassi*angi(1)
          angy = angy + pmassi*angi(2)
          angz = angz + pmassi*angi(3)

          ! rotational energy around each axis through the origin
          if (calc_erot) then
             call get_erot(xi,yi,zi,vxi,vyi,vzi,xyzcom,pmassi,erotxi,erotyi,erotzi)
             call ev_data_update(ev_data_thread,iev_erot(1),erotxi)
             call ev_data_update(ev_data_thread,iev_erot(2),erotyi)
             call ev_data_update(ev_data_thread,iev_erot(3),erotzi)
          endif
       enddo
       !$omp enddo
    endif
 endif

!$omp critical(collatedata)
 call collate_ev_data(ev_data_thread,ev_data)
 if (.not.gas_only) then
    do i = 1,maxtypes
       np_rho(i) = np_rho(i) + np_rho_thread(i)
    enddo
 endif
!$omp end critical(collatedata)
!$omp end parallel

 !--Determing the number of active gas particles
 nptot    = reduceall_mpi('+',np)
 npgas    = reduceall_mpi('+',npgas)
 npartall = reduceall_mpi('+',npart)
 ndead    = npartall - nptot
 if (nptot > 0) then
    dnptot = 1./real(nptot)
 else
    dnptot = 0.
 endif
 if (npgas > 0) then
    dnpgas = 1./real(npgas)
 else
    dnpgas = 0.
 endif
 !--Number of gas particles without a sound speed or energy
 np_cs_eq_0 = reduceall_mpi('+',np_cs_eq_0)
 np_e_eq_0  = reduceall_mpi('+',np_e_eq_0)
 !
 !--Finalise the arrays & correct as necessary;
 !  Almost all of the average quantities are over gas particles only
 !
 call finalise_ev_data(ev_data,dnpgas)

 if (.not.gr) ekin = 0.5*ekin
 emag = 0.5*emag
 ekin = reduceall_mpi('+',ekin)
 !LS I don't know what to do here ? gamma should be replaced by gammai ?
 if (maxvxyzu >= 4 .or. gamma >= 1.0001) etherm = reduceall_mpi('+',etherm)
 emag = reduceall_mpi('+',emag)
 epot = reduceall_mpi('+',epot)
 erad = reduceall_mpi('+',erad)
 if (nptmass > 1) then
    if (use_regnbody) then
       call get_pot_subsys(n_group,group_info,bin_info,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad,epot_sinksink)
    endif
    epot = epot + epot_sinksink
 endif

 etot = ekin + etherm + emag + epot + erad
 etotall = etot

 xcom = reduceall_mpi('+',xcom)
 ycom = reduceall_mpi('+',ycom)
 zcom = reduceall_mpi('+',zcom)
 mtot = reduceall_mpi('+',mtot)
 if (mtot > 0.0) dm = 1.0 / mtot
 xcom = xcom * dm
 ycom = ycom * dm
 zcom = zcom * dm

 xmom = reduceall_mpi('+',xmom)
 ymom = reduceall_mpi('+',ymom)
 zmom = reduceall_mpi('+',zmom)
 totmom = sqrt(xmom*xmom + ymom*ymom + zmom*zmom)

 angx = reduceall_mpi('+',angx)
 angy = reduceall_mpi('+',angy)
 angz = reduceall_mpi('+',angz)
 angtot = sqrt(angx*angx + angy*angy + angz*angz)

 vrms    = reduceall_mpi('+',vrms)
 rmsmach = reduceall_mpi('+',rmsmach)
 vrms    = sqrt(vrms*dnptot)
 rmsmach = sqrt(rmsmach*dnpgas)

 !--fill in the relevant array elements for energy & momentum
 ev_data(iev_sum,iev_time  ) = t
 ev_data(iev_sum,iev_ekin  ) = ekin
 ev_data(iev_sum,iev_etherm) = etherm
 ev_data(iev_sum,iev_emag  ) = emag
 ev_data(iev_sum,iev_epot  ) = epot
 ev_data(iev_sum,iev_etot  ) = etot
 ev_data(iev_sum,iev_erad  ) = erad
 ev_data(iev_sum,iev_totmom) = totmom
 ev_data(iev_sum,iev_angmom) = angtot
 ev_data(iev_sum,iev_com(1)) = xcom
 ev_data(iev_sum,iev_com(2)) = ycom
 ev_data(iev_sum,iev_com(3)) = zcom
 do i=1,ndusttypes
    ev_data(iev_sum,iev_dm(i)) = mdust(i)
 enddo
 ev_data(iev_sum,iev_vrms   ) = vrms
 ev_data(iev_sum,iev_rmsmach) = rmsmach
 xyzcom(1) = xcom
 xyzcom(2) = ycom
 xyzcom(3) = zcom

 if (calc_erot) then
    ev_data(iev_sum,iev_erot(1)) = 0.5*ev_data(iev_sum,iev_erot(1))
    ev_data(iev_sum,iev_erot(2)) = 0.5*ev_data(iev_sum,iev_erot(2))
    ev_data(iev_sum,iev_erot(3)) = 0.5*ev_data(iev_sum,iev_erot(3))
    ev_data(iev_sum,iev_erot(4)) = sqrt(ev_data(iev_sum,iev_erot(1))**2 &
                                 +      ev_data(iev_sum,iev_erot(2))**2 &
                                 +      ev_data(iev_sum,iev_erot(3))**2)
 endif

 if (use_dust) then
    mgas  = reduceall_mpi('+',mgas)
    mdust = reduceall_mpi('+',mdust)
 endif

 if (.not. gas_only) then
    do i = 1,maxtypes
       np_rho(i) = reduceall_mpi('+',np_rho(i))
    enddo
    ! correct the average densities so that division is by n_p and not n_gas
    ev_data(iev_ave,iev_rho) = ev_data(iev_ave,iev_rho)*real(npgas)*dnptot
    if (np_rho(idust)       > 0) ev_data(iev_ave,iev_rhop(2)) = ev_data(iev_ave,iev_rhop(2))*real(npgas)/real(np_rho(idust))
    if (np_rho(iboundary)   > 0) ev_data(iev_ave,iev_rhop(3)) = ev_data(iev_ave,iev_rhop(3))*real(npgas)/real(np_rho(iboundary))
    if (np_rho(istar)       > 0) ev_data(iev_ave,iev_rhop(4)) = ev_data(iev_ave,iev_rhop(4))*real(npgas)/real(np_rho(istar))
    if (np_rho(idarkmatter) > 0) ev_data(iev_ave,iev_rhop(5)) = ev_data(iev_ave,iev_rhop(5))*real(npgas)/real(np_rho(idarkmatter))
    if (np_rho(ibulge)      > 0) ev_data(iev_ave,iev_rhop(6)) = ev_data(iev_ave,iev_rhop(6))*real(npgas)/real(np_rho(ibulge))
 endif

 if (iexternalforce > 0) then
    xmomacc   = reduceall_mpi('+',xmomacc)
    ymomacc   = reduceall_mpi('+',ymomacc)
    zmomacc   = reduceall_mpi('+',zmomacc)

    xmomall   = xmom + xmomacc
    ymomall   = ymom + ymomacc
    zmomall   = zmom + zmomacc
    ev_data(iev_sum,iev_momall) = sqrt(xmomall*xmomall + ymomall*ymomall + zmomall*zmomall)

    angaccx = reduceall_mpi('+',angaccx)
    angaccy = reduceall_mpi('+',angaccy)
    angaccz = reduceall_mpi('+',angaccz)
    angxall = angx + angaccx
    angyall = angy + angaccy
    angzall = angz + angaccz
    angall  = sqrt(angxall*angxall + angyall*angyall + angzall*angzall)
    ev_data(iev_sum,iev_angall) = angall

    ekinacc   = reduceall_mpi('+',ekinacc)
    epotacc   = reduceall_mpi('+',epotacc)
    ethermacc = reduceall_mpi('+',ethermacc)
    emagacc   = reduceall_mpi('+',emagacc)
    eradacc   = reduceall_mpi('+',eradacc)
    eacc      = ekinacc + ethermacc + emagacc + epotacc + eradacc
    etotall   = etotall + eacc
 endif

 if (track_mass) then
    accretedmass = ev_data(iev_sum,iev_macc)
    if (accradius1 > 0.) then
       !eacc = accretedmass/accradius1
       ev_data(iev_sum,iev_eacc) = eacc ! total accretion energy
    endif
 endif
 if (track_lum) totlum = ev_data(iev_sum,iev_totlum)

 if (calc_gravitwaves) then
    if (use_apr) then
       pmassi = aprmassoftype(igas,apr_level(i))
    else
       pmassi = massoftype(igas)
    endif
    x0 = 0.; v0 = 0.; a0 = 0.  ! use the origin by default
    if (gr) then
       !call get_geodesic_accel(axyz,npart,vxyzu(1:3,:),metrics,metricderivs)
       !call calculate_strain(hx,hp,pmassi,x0,v0,a0,npart,xyzh,vxyzu,axyz)
       call calculate_strain(hx,hp,pmassi,ddq_xy,x0,v0,a0,npart,xyzh,vxyzu(1:3,:),fxyzu,&
              fext,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)
    else
       if (iexternalforce==0) then  ! if no external forces, use centre of mass of particles
          x0 = (/xcom,ycom,zcom/)
          v0 = (/xmom,ymom,zmom/)
          call get_centreofmass_accel(a0,npart,xyzh,fxyzu,fext,nptmass,xyzmh_ptmass,fxyz_ptmass)
       endif
       call calculate_strain(hx,hp,pmassi,ddq_xy,x0,v0,a0,npart,xyzh,vxyzu(1:3,:),fxyzu,&
              fext,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)
    endif
    ev_data(iev_sum,iev_gws(1)) = hx(1)
    ev_data(iev_sum,iev_gws(2)) = hp(1)
    ev_data(iev_sum,iev_gws(3)) = hx(2)
    ev_data(iev_sum,iev_gws(4)) = hp(2)
    ev_data(iev_sum,iev_gws(5)) = hx(3)
    ev_data(iev_sum,iev_gws(6)) = hp(3)
    ev_data(iev_sum,iev_gws(7)) = hx(4)
    ev_data(iev_sum,iev_gws(8)) = hp(4)
 endif

 if (mhd) then
    hdivBonB_max = ev_data(iev_max,iev_hdivB)
    hdivBonB_ave = ev_data(iev_ave,iev_hdivB)
 endif

 if (maxp==maxp_hard) then
    ev_data(iev_sum,iev_mass) = mtot
 endif

 if (dynamic_bdy) then
    call find_dynamic_boundaries(npart,nptmass,dtmax,xyz_n_all,xyz_x_all,ierr)
    ev_data(iev_sum,iev_bdy(1,1)) = xyz_n_all(1)
    ev_data(iev_sum,iev_bdy(1,2)) = xyz_x_all(1)
    ev_data(iev_sum,iev_bdy(2,1)) = xyz_n_all(2)
    ev_data(iev_sum,iev_bdy(2,2)) = xyz_x_all(2)
    ev_data(iev_sum,iev_bdy(3,1)) = xyz_n_all(3)
    ev_data(iev_sum,iev_bdy(3,2)) = xyz_x_all(3)
    if (ierr==1) call fatal('energies','there is no high density gas for the dynamic boundaries')
 endif

end subroutine compute_energies
!----------------------------------------------------------------
!+
!  calculates rotational energy
!+
!----------------------------------------------------------------
subroutine get_erot(xi,yi,zi,vxi,vyi,vzi,xyzcom,pmassi,erotxi,erotyi,erotzi)
 real, intent(in)  :: xi,yi,zi,vxi,vyi,vzi,pmassi,xyzcom(3)
 real, intent(out) :: erotxi,erotyi,erotzi
 real              :: dx,dy,dz,dvx,dvy,dvz
 real              :: rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2

 erotxi = 0.0
 erotyi = 0.0
 erotzi = 0.0

 dx  = xi  - xyzcom(1)
 dy  = yi  - xyzcom(2)
 dz  = zi  - xyzcom(3)
 dvx = vxi              ! results are less reliable if subtracting vcom
 dvy = vyi
 dvz = vzi

 rcrossvx = (dy*dvz - dz*dvy)
 rcrossvy = (dz*dvx - dx*dvz)
 rcrossvz = (dx*dvy - dy*dvx)

 radxy2 = dx*dx + dy*dy
 radyz2 = dy*dy + dz*dz
 radxz2 = dx*dx + dz*dz

 if (radyz2 > 0.) erotxi = pmassi*rcrossvx*rcrossvx/radyz2
 if (radxz2 > 0.) erotyi = pmassi*rcrossvy*rcrossvy/radxz2
 if (radxy2 > 0.) erotzi = pmassi*rcrossvz*rcrossvz/radxy2

end subroutine get_erot
!----------------------------------------------------------------
!+
!  initialise the ev_data array
!+
!----------------------------------------------------------------
subroutine initialise_ev_data(evdata)
 real,    intent(inout) :: evdata(4,0:inumev)

 evdata            = 0.0
 evdata(iev_max,:) = -huge(evdata(iev_max,:))
 evdata(iev_min,:) =  huge(evdata(iev_min,:))

end subroutine initialise_ev_data
!----------------------------------------------------------------
!+
!  update the ev_data_array
!+
!----------------------------------------------------------------
subroutine ev_data_update(evdata,itag,val)
 integer, intent(in)    :: itag
 real,    intent(in)    :: val
 real,    intent(inout) :: evdata(4,0:inumev)

 evdata(iev_sum,itag) =     evdata(iev_sum,itag)+val
 evdata(iev_max,itag) = max(evdata(iev_max,itag),val)
 evdata(iev_min,itag) = min(evdata(iev_min,itag),val)

end subroutine ev_data_update
!----------------------------------------------------------------
!+
!  combines the ev_data from the various threads
!+
!----------------------------------------------------------------
subroutine collate_ev_data(evdata_thread,evdata)
 real,            intent(in)    :: evdata_thread(4,0:inumev)
 real,            intent(inout) :: evdata(4,0:inumev)
 integer                        :: i

 do i = 1,iquantities
    evdata(iev_sum,i) =     evdata(iev_sum,i)+evdata_thread(iev_sum,i)
    evdata(iev_max,i) = max(evdata(iev_max,i),evdata_thread(iev_max,i))
    evdata(iev_min,i) = min(evdata(iev_min,i),evdata_thread(iev_min,i))
 enddo

end subroutine collate_ev_data
!----------------------------------------------------------------
!+
!  Performs final generic housekeeping on the ev_data array
!+
!----------------------------------------------------------------
subroutine finalise_ev_data(evdata,dnptot)
 use mpiutils, only:reduceall_mpi
 real,            intent(inout) :: evdata(4,0:inumev)
 real,            intent(in)    :: dnptot
 integer                        :: i

 do i = 1,iquantities
    evdata(iev_sum,i) = reduceall_mpi('+',  evdata(iev_sum,i))
    evdata(iev_max,i) = reduceall_mpi('max',evdata(iev_max,i))
    evdata(iev_min,i) = reduceall_mpi('min',evdata(iev_min,i))
    evdata(iev_ave,i) = evdata(iev_sum,i)*dnptot
 enddo

end subroutine finalise_ev_data

end module energies
