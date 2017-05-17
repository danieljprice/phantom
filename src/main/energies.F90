!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: energies
!
!  DESCRIPTION:
!   Calculates global quantities on the particles
!   To Developer: See instructions in evwrite.F90 about adding values
!                 to the .ev file
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, dim, dust, eos, externalforces, io,
!    mpiutils, nicil, options, part, ptmass, viscosity
!+
!--------------------------------------------------------------------------

#define reduce_fn(a,b) reduceall_mpi(a,b)
module energies
 use dim, only: calc_erot, calc_erot_com
 implicit none

 logical,         public    :: gas_only,track_mass,track_lum
 real,            public    :: ekin,etherm,emag,epot,etot,totmom,angtot
 real,            public    :: vrms,rmsmach,accretedmass,mdust,mgas
 real,            public    :: xmom,ymom,zmom
 real,            public    :: totlum
 integer,         public    :: ielements,iquantities
 integer,         parameter :: inumev = 150
 ! Subroutines
 public  :: compute_energies,get_erot_com,ev_data_update
 private :: get_erot,initialise_ev_data,collate_ev_data,finalise_ev_data,ev_data_correction
 ! Functions
 public  :: ev_get_value
 ! Arrays
 real,             public :: ev_data(0:inumev),erot_com(6)
 character(len=11),public :: ev_tag(inumev)    ! the tag of the various quantities to print
 character(len= 3),public :: ev_action(inumev) ! the actions to be performed on the given quantity
 integer,          public :: ev_istart(inumev) ! the first array element in ev_data

contains

!----------------------------------------------------------------
!+
!  Subroutine to compute global quantities on the particles
!+
!----------------------------------------------------------------
subroutine compute_energies(t)
 use dim,  only:maxp,maxvxyzu,maxalpha,maxtypes,use_dustfrac,mhd_nonideal,lightcurve
 use part, only:rhoh,xyzh,vxyzu,massoftype,npart,maxphase,iphase,npartoftype, &
                alphaind,Bxyz,Bevol,divcurlB,iamtype,igas,idust,iboundary,istar,idarkmatter,ibulge, &
                nptmass,xyzmh_ptmass,vxyz_ptmass,isdeadh,isdead_or_accreted,epot_sinksink,&
                imacc,ispinx,ispiny,ispinz,mhd,maxvecp,divBsymm,gravity,poten,dustfrac,&
                n_R,n_electronT,ionfrac_eta
 use eos,            only:polyk,utherm,gamma,equationofstate,get_temperature_from_ponrho,gamma_pwp
 use io,             only:id,fatal,master
 use externalforces, only:externalforce,externalforce_vdependent,was_accreted,accradius1
 use options,        only:iexternalforce,alpha,alphaB,ieos
 use mpiutils,       only:reduceall_mpi
 use ptmass,         only:get_accel_sink_gas
 use viscosity,      only:irealvisc,shearfunc
 use nicil,          only:nicil_get_eta,nicil_get_vion,use_ohm,use_hall,use_ambi,ion_rays,ion_thermal, &
                     nelements_max,nelements,nlevels
#ifdef DUST
 use dust,           only:get_ts,graindens,grainsize,idrag
 integer :: iregime
 real    :: tsi
#endif
#ifdef LIGHTCURVE
 use part,         only:luminosity
#endif
 real, intent(in) :: t
 real    :: ev_data_thread(0:inumev)
 real    :: xi,yi,zi,hi,vxi,vyi,vzi,v2i,Bxi,Byi,Bzi,rhoi,angx,angy,angz
 real    :: xmomacc,ymomacc,zmomacc,angaccx,angaccy,angaccz
 real    :: epoti,pmassi,acci,dnptot
 real    :: xmomall,ymomall,zmomall,angxall,angyall,angzall,rho1i,vsigi
 real    :: ponrhoi,spsoundi,B2i,dumx,dumy,dumz,divBi,hdivBonBi,alphai,valfven2i,betai
 real    :: n_total,n_ion,shearparam_art,shearparam_phys,ratio_phys_to_av
 real    :: gasfrac,dustfraci,dust_to_gas
 real    :: temperature,etaart,etaart1,etaohm,etahall,etaambi,vion,vdrift
 real    :: vioni(3),data_out(17+nelements_max*nlevels-3)
 real    :: erotxi,erotyi,erotzi,fdum(3)
 integer :: i,j,itype,ierr
 integer(kind=8) :: np,nptot,np_rho(maxtypes),np_rho_thread(maxtypes)

 ! initialise values
 itype  = igas
 pmassi = massoftype(igas)
 ekin   = 0.
 etherm = 0.
 if (maxvxyzu < 4 .and. gamma < 1.0001 .and. ieos/=9) etherm = 1.5*polyk
 epot = 0.
 emag = 0.
 etot = 0.
 xmom = 0.
 ymom = 0.
 zmom = 0.
 angx = 0.
 angy = 0.
 angz = 0.
 np = 0
 xmomacc = 0.
 ymomacc = 0.
 zmomacc = 0.
 angaccx = 0.
 angaccy = 0.
 angaccz = 0.
 mgas    = 0.
 mdust   = 0.
 if (maxalpha==maxp) then
    alphai = 0.
 else
    alphai = alpha
 endif
 ionfrac_eta = 0.
 np_rho      = 0
 call initialise_ev_data(ev_data)
!
!$omp parallel default(none) &
!$omp shared(xyzh,vxyzu,iexternalforce,npart,t,id,npartoftype) &
!$omp shared(alphaind,massoftype,irealvisc) &
!$omp shared(ieos,gamma,nptmass,xyzmh_ptmass,vxyz_ptmass) &
!$omp shared(Bxyz,Bevol,divBsymm,divcurlB,alphaB,iphase,poten,dustfrac) &
!$omp shared(use_ohm,use_hall,use_ambi,ion_rays,ion_thermal,nelements,n_R,n_electronT,ionfrac_eta) &
!$omp shared(ielements,ev_data,np_rho,erot_com,calc_erot,gas_only,track_mass) &
!$omp private(i,j,xi,yi,zi,hi,rhoi,vxi,vyi,vzi,Bxi,Byi,Bzi,epoti,vsigi,v2i) &
!$omp private(ponrhoi,spsoundi,B2i,dumx,dumy,dumz,acci,valfven2i,divBi,hdivBonBi) &
!$omp private(rho1i,shearparam_art,shearparam_phys,ratio_phys_to_av,betai) &
!$omp private(gasfrac,dustfraci,dust_to_gas,n_total,n_ion) &
!$omp private(ierr,temperature,etaart,etaart1,etaohm,etahall,etaambi,vioni,vion,vdrift,data_out) &
!$omp private(erotxi,erotyi,erotzi,fdum) &
!$omp private(ev_data_thread,np_rho_thread) &
!$omp firstprivate(alphai,itype,pmassi) &
#ifdef DUST
!$omp shared(grainsize,graindens,idrag) &
!$omp private(tsi,iregime) &
#endif
#ifdef LIGHTCURVE
!$omp shared(luminosity,track_lum) &
#endif
!$omp reduction(+:np,xmom,ymom,zmom,angx,angy,angz) &
!$omp reduction(+:xmomacc,ymomacc,zmomacc,angaccx,angaccy,angaccz) &
!$omp reduction(+:ekin,etherm,emag,epot)
 call initialise_ev_data(ev_data_thread)
 np_rho_thread = 0
!$omp do
 do i=1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
          if (itype <= 0) call fatal('energies','particle type <= 0')
          pmassi = massoftype(itype)
       endif

       rhoi = rhoh(hi,pmassi)
       call ev_data_update(ev_data_thread,'rho',rhoi)
       if (.not.gas_only) then
          if (npartoftype(igas)        > 0) then
             call ev_data_update(ev_data_thread,'rho gas', rhoi)
             np_rho_thread(igas) =  np_rho_thread(igas) + 1
          endif
          if (npartoftype(idust)       > 0) then
             call ev_data_update(ev_data_thread,'rho dust',rhoi)
             np_rho_thread(idust) =  np_rho_thread(idust) + 1
          endif
          if (npartoftype(iboundary)   > 0) then
             call ev_data_update(ev_data_thread,'rho bdy', rhoi)
             np_rho_thread(iboundary) =  np_rho_thread(iboundary) + 1
          endif
          if (npartoftype(istar)       > 0) then
             call ev_data_update(ev_data_thread,'rho star',rhoi)
             np_rho_thread(istar) =  np_rho_thread(istar) + 1
          endif
          if (npartoftype(idarkmatter) > 0) then
             call ev_data_update(ev_data_thread,'rho dm',  rhoi)
             np_rho_thread(idarkmatter) =  np_rho_thread(idarkmatter) + 1
          endif
          if (npartoftype(ibulge)      > 0) then
             call ev_data_update(ev_data_thread,'rho blg', rhoi)
             np_rho_thread(ibulge) =  np_rho_thread(ibulge) + 1
          endif
       endif

       np   = np + 1

       vxi  = vxyzu(1,i)
       vyi  = vxyzu(2,i)
       vzi  = vxyzu(3,i)

       !  linear momentum
       xmom = xmom + pmassi*vxi
       ymom = ymom + pmassi*vyi
       zmom = zmom + pmassi*vzi

       ! angular momentum
       angx = angx + pmassi*(yi*vzi - zi*vyi)
       angy = angy + pmassi*(zi*vxi - xi*vzi)
       angz = angz + pmassi*(xi*vyi - yi*vxi)

       ! kinetic energy & rms velocity
       v2i  = vxi*vxi + vyi*vyi + vzi*vzi
       ekin = ekin + pmassi*v2i
       call ev_data_update(ev_data_thread,'vrms',v2i)        ! vrms = vrms + v2i

       ! rotational energy around each axis through the Centre of mass
       ! note: centre of mass is updated only when dumpfiles are created
       if (calc_erot) then
          call get_erot(xi,yi,zi,vxi,vyi,vzi,pmassi,erotxi,erotyi,erotzi)
          call ev_data_update(ev_data_thread,'erot_x',erotxi)
          call ev_data_update(ev_data_thread,'erot_y',erotyi)
          call ev_data_update(ev_data_thread,'erot_z',erotzi)
       endif

       if (iexternalforce > 0) then
          dumx = 0.
          dumy = 0.
          dumz = 0.
          call externalforce(iexternalforce,xi,yi,zi,hi,t,dumx,dumy,dumz,epoti,ii=i)
          call externalforce_vdependent(iexternalforce,xyzh(1:3,i),vxyzu(1:3,i),fdum,epoti)
          epot = epot + pmassi*epoti
       endif
       if (nptmass > 0) then
          dumx = 0.
          dumy = 0.
          dumz = 0.
          call get_accel_sink_gas(nptmass,xi,yi,zi,hi,xyzmh_ptmass,dumx,dumy,dumz,epoti)
          epot = epot + pmassi*epoti
       endif
       if (gravity) epot = epot + poten(i)
       !
       ! the following apply ONLY to gas particles
       !
       isgas: if (itype==igas) then

          if (use_dustfrac) then
             dustfraci   = dustfrac(i)
             gasfrac     = 1. - dustfraci
             dust_to_gas = dustfraci/gasfrac
             call ev_data_update(ev_data_thread,'dust/gas',dust_to_gas     )
             call ev_data_update(ev_data_thread,'mgas',    pmassi*gasfrac  ) ! mgas  = mgas  + pmassi*gasfrac
             call ev_data_update(ev_data_thread,'mdust',   pmassi*dustfraci) ! mdust = mdust + pmassi*dustfraci
          else
             dustfraci = 0.
             gasfrac   = 1.
          endif

          ! thermal energy
          if (maxvxyzu >= 4) then
             etherm = etherm + pmassi*utherm(vxyzu(4,i),rhoi)*gasfrac
             call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,vxyzu(4,i))
          else
             call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi)
             if (ieos==2 .and. gamma > 1.001) then
                !--thermal energy using polytropic equation of state
                etherm = etherm + pmassi*ponrhoi/(gamma-1.)*gasfrac
             else if (ieos==9) then
                !--thermal energy using piecewise polytropic equation of state
                etherm = etherm + pmassi*ponrhoi/(gamma_pwp(rhoi)-1.)*gasfrac
             endif
          endif
          vsigi = spsoundi
          ! entropy
          call ev_data_update(ev_data_thread,'totentrop',pmassi*ponrhoi*rhoi**(1.-gamma))

#ifdef DUST
          ! min and mean stopping time
          if (use_dustfrac) then
             call get_ts(idrag,grainsize,graindens,rhoi*gasfrac,rhoi*dustfraci,spsoundi,0.,tsi,iregime)
             call ev_data_update(ev_data_thread,'t_s',tsi)
          endif
#endif

#ifdef LIGHTCURVE
          if (track_lum) call ev_data_update(ev_data_thread,'tot lum',real(luminosity(i)))
#endif

          ! rms mach number (rmsmach = rmsmach + v2i/spsoundi**2)
          if (spsoundi > 0.) call ev_data_update(ev_data_thread,'rmsmach',v2i/spsoundi**2)

          ! max of dissipation parameters
          if (maxalpha==maxp) then
             alphai = alphaind(1,i)
             call ev_data_update(ev_data_thread,'alpha',alphai)
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
             call ev_data_update(ev_data_thread,'visc_rat',ratio_phys_to_av)
          endif

          ! mhd parameters
          if (mhd) then
             if (maxvecp==maxp) then
                Bxi = Bxyz(1,i)
                Byi = Bxyz(2,i)
                Bzi = Bxyz(3,i)
             else
                Bxi = Bevol(1,i)
                Byi = Bevol(2,i)
                Bzi = Bevol(3,i)
             endif
             B2i       = Bxi*Bxi + Byi*Byi + Bzi*Bzi
             rho1i     = 1./rhoi
             valfven2i = B2i*rho1i
             vsigi     = sqrt(valfven2i + spsoundi*spsoundi)
             emag      = emag + pmassi*B2i*rho1i

             divBi     = abs(divcurlB(1,i))
             if (B2i > 0.) then
                hdivBonBi = hi*divBi/sqrt(B2i)
                betai     = 2.0*ponrhoi*rhoi/B2i ! plasma beta
             else
                hdivBonBi = 0.
                betai     = 0.
             endif
             call ev_data_update(ev_data_thread,'divB',   divBi    )
             call ev_data_update(ev_data_thread,'hdivB/B',hdivBonBi)
             call ev_data_update(ev_data_thread,'beta',   betai    )

             if ( mhd_nonideal ) then
                temperature = get_temperature_from_ponrho(ponrhoi)
                call nicil_get_eta(etaohm,etahall,etaambi,sqrt(B2i),rhoi,temperature, &
                                   n_R(:,i),n_electronT(i),ierr,data_out)
                call nicil_get_vion(etaambi,vxi,vyi,vzi,Bxi,Byi,Bzi,real(divcurlB(2:4,i)),vioni,ierr)
                etaart  = 0.5*hi*vsigi*alphaB
                if (etaart > 0.) then
                   etaart1 = 1.0/etaart
                else
                   etaart1 = 0.0
                endif
                call ev_data_update(ev_data_thread,'temp',  temperature)
                call ev_data_update(ev_data_thread,'eta_ar',etaart     )
                if (use_ohm) then
                   call ev_data_update(ev_data_thread,'eta_o',    etaohm              )
                   call ev_data_update(ev_data_thread,'eta_o/art',etaohm*etaart1      )
                endif
                if (use_hall) then
                   call ev_data_update(ev_data_thread,'eta_h',        etahall         )
                   call ev_data_update(ev_data_thread,'|eta_h|',  abs(etahall)        )
                   call ev_data_update(ev_data_thread,'eta_h/art',    etahall*etaart1 )
                   call ev_data_update(ev_data_thread,'|e_h|/art',abs(etahall)*etaart1)
                endif
                if (use_ambi) then
                   vion   = sqrt( dot_product(vioni,vioni) )
                   vdrift = sqrt( (vioni(1)-vxi)**2 + (vioni(2)-vyi)**2 + (vioni(3)-vzi)**2)
                   call ev_data_update(ev_data_thread,'eta_a',    etaambi        )
                   call ev_data_update(ev_data_thread,'eta_a/art',etaambi*etaart1)
                   call ev_data_update(ev_data_thread,'velocity', sqrt(v2i)      )
                   call ev_data_update(ev_data_thread,'v_ion',    vion           )
                   call ev_data_update(ev_data_thread,'v_drift',  vdrift         )
                endif
                n_ion   = data_out(8) + data_out(9) + data_out(10) + data_out(11)
                n_total = n_ion + data_out(7)
                call ev_data_update(ev_data_thread,'ni/n(i+n)',n_ion/n_total)
                call ev_data_update(ev_data_thread,'ne/n(i+n)',data_out(6)/n_total)
                ionfrac_eta(1,i) = real(n_ion/n_total,kind=4)
                ionfrac_eta(2,i) = real(etaohm, kind=4)       ! Save eta_OR for the dump file
                ionfrac_eta(3,i) = real(etahall,kind=4)       ! Save eta_HE for the dump file
                ionfrac_eta(4,i) = real(etaambi,kind=4)       ! Save eta_AD for the dump file
                call ev_data_update(ev_data_thread,   'n_e',      data_out( 6))
                call ev_data_update(ev_data_thread,   'n_n',      data_out( 7))
                if (ion_rays) then
                   call ev_data_update(ev_data_thread,'n_ihR',    data_out( 8))
                   call ev_data_update(ev_data_thread,'n_imR',    data_out( 9))
                   call ev_data_update(ev_data_thread,'n_g(Z=-1)',data_out(12))
                   call ev_data_update(ev_data_thread,'n_g(Z= 0)',data_out(13))
                   call ev_data_update(ev_data_thread,'n_g(Z=+1)',data_out(14))
                endif
                if (ion_thermal) then
                   call ev_data_update(ev_data_thread,'n_isT',    data_out(10))
                   call ev_data_update(ev_data_thread,'n_idT',    data_out(11))
                endif
             endif
          endif
       endif isgas

    elseif (was_accreted(iexternalforce,hi)) then
!
!--count accretion onto fixed potentials (external forces) separately
!
       vxi = vxyzu(1,i)
       vyi = vxyzu(2,i)
       vzi = vxyzu(3,i)
       if (maxphase==maxp) then
          pmassi = massoftype(iamtype(iphase(i)))
       else
          pmassi = massoftype(igas)
       endif
       xmomacc = xmomacc + pmassi*vxi
       ymomacc = ymomacc + pmassi*vyi
       zmomacc = zmomacc + pmassi*vzi

       angaccx = angaccx + pmassi*(yi*vzi - zi*vyi)
       angaccy = angaccy + pmassi*(zi*vxi - xi*vzi)
       angaccz = angaccz + pmassi*(xi*vyi - yi*vxi)

       if (track_mass) call ev_data_update(ev_data_thread,'accretedmas',pmassi)

    endif
 enddo
!$omp enddo
!
!--add contribution from sink particles
!

!$omp do
 do i=1,nptmass
    xi     = xyzmh_ptmass(1,i)
    yi     = xyzmh_ptmass(2,i)
    zi     = xyzmh_ptmass(3,i)
    pmassi = xyzmh_ptmass(4,i)
    !--acci is the accreted mass on the sink
    acci   = xyzmh_ptmass(imacc,i)
    if (track_mass) call ev_data_update(ev_data_thread,'accretedmas',acci)

    vxi    = vxyz_ptmass(1,i)
    vyi    = vxyz_ptmass(2,i)
    vzi    = vxyz_ptmass(3,i)

    !phii   = fxyz_ptmass(4,i)

    xmom   = xmom + pmassi*vxi
    ymom   = ymom + pmassi*vyi
    zmom   = zmom + pmassi*vzi

    angx = angx + pmassi*(yi*vzi - zi*vyi)
    angy = angy + pmassi*(zi*vxi - xi*vzi)
    angz = angz + pmassi*(xi*vyi - yi*vxi)

    angx   = angx + xyzmh_ptmass(ispinx,i)
    angy   = angy + xyzmh_ptmass(ispiny,i)
    angz   = angz + xyzmh_ptmass(ispinz,i)

    v2i    = vxi*vxi + vyi*vyi + vzi*vzi
    ekin   = ekin + pmassi*v2i

    ! rotational energy around each axis through the origin
    if (calc_erot) then
       call get_erot(xi,yi,zi,vxi,vyi,vzi,pmassi,erotxi,erotyi,erotzi)
       call ev_data_update(ev_data_thread,'erot_x',erotxi)
       call ev_data_update(ev_data_thread,'erot_y',erotyi)
       call ev_data_update(ev_data_thread,'erot_z',erotzi)
    endif
 enddo
!$omp enddo
!$omp critical(collatedata)
 call collate_ev_data(ev_data_thread,ev_data)
 if (.not.gas_only) then
    do i = 1,maxtypes
       np_rho(i) = np_rho(i) + np_rho_thread(i)
    enddo
 endif
!$omp end critical(collatedata)
!$omp end parallel

 !--Determing the number of active particles
 nptot     = reduce_fn('+',np)
 if (nptot > 0.) then
    dnptot = 1./real(nptot)
 else
    dnptot = 0.
 endif
 !--Finalise the arrays & correct as necessary
 call finalise_ev_data(ev_data,dnptot)
 ekin = 0.5*ekin
 emag = 0.5*emag
 ekin = reduce_fn('+',ekin)
 if (maxvxyzu >= 4 .or. gamma >= 1.0001) etherm = reduce_fn('+',etherm)
 emag = reduce_fn('+',emag)
 epot = reduce_fn('+',epot)
 if (nptmass > 1) epot = epot + epot_sinksink

 etot = ekin + etherm + emag + epot

 xmom = reduce_fn('+',xmom)
 ymom = reduce_fn('+',ymom)
 zmom = reduce_fn('+',zmom)
 totmom = sqrt(xmom*xmom + ymom*ymom + zmom*zmom)

 angx = reduce_fn('+',angx)
 angy = reduce_fn('+',angy)
 angz = reduce_fn('+',angz)
 angtot = sqrt(angx*angx + angy*angy + angz*angz)

 !--fill in the relevant array elements for energy & momentum
 call ev_data_update(ev_data,'time',  t     )
 call ev_data_update(ev_data,'ekin',  ekin  )
 call ev_data_update(ev_data,'etherm',etherm)
 call ev_data_update(ev_data,'emag',  emag  )
 call ev_data_update(ev_data,'epot',  epot  )
 call ev_data_update(ev_data,'etot',  etot  )
 call ev_data_update(ev_data,'totmom',totmom)
 call ev_data_update(ev_data,'angtot',angtot)

 if (calc_erot) then
    call ev_data_correction(ev_data,'erot_x',0.5)  ! erot = 0.5*erot
    call ev_data_correction(ev_data,'erot_y',0.5)  ! erot = 0.5*erot
    call ev_data_correction(ev_data,'erot_z',0.5)  ! erot = 0.5*erot
    call ev_data_update(ev_data,'erot', sqrt(ev_get_value('erot_x')**2+ev_get_value('erot_y')**2+ev_get_value('erot_z')**2))
 endif

 if (use_dustfrac) then
    mgas  = ev_get_value('mgas')
    mdust = ev_get_value('mdust')
 endif

 if (.not. gas_only) then
    do i = 1,maxtypes
       np_rho(i) = reduce_fn('+',np_rho(i))
    enddo
    ! correct the average densities so that division is by n_p and not n_total
    if (np_rho(igas)        > 0) call ev_data_correction(ev_data,'rho gas', real(nptot)/real(np_rho(igas)),       'a')
    if (np_rho(idust)       > 0) call ev_data_correction(ev_data,'rho dust',real(nptot)/real(np_rho(idust)),      'a')
    if (np_rho(iboundary)   > 0) call ev_data_correction(ev_data,'rho bdy', real(nptot)/real(np_rho(iboundary)),  'a')
    if (np_rho(istar)       > 0) call ev_data_correction(ev_data,'rho star',real(nptot)/real(np_rho(istar)),      'a')
    if (np_rho(idarkmatter) > 0) call ev_data_correction(ev_data,'rho dm',  real(nptot)/real(np_rho(idarkmatter)),'a')
    if (np_rho(ibulge)      > 0) call ev_data_correction(ev_data,'rho blg', real(nptot)/real(np_rho(ibulge)),     'a')
 endif
 call ev_data_correction(ev_data,'vrms',   dnptot)
 call ev_data_correction(ev_data,'rmsmach',dnptot)
 vrms    = ev_get_value('vrms')
 rmsmach = ev_get_value('rmsmach')

 if (iexternalforce > 0) then
    xmomacc   = reduce_fn('+',xmomacc)
    ymomacc   = reduce_fn('+',ymomacc)
    zmomacc   = reduce_fn('+',zmomacc)

    xmomall   = xmom + xmomacc
    ymomall   = ymom + ymomacc
    zmomall   = zmom + zmomacc
    call ev_data_update(ev_data,'totmomall',sqrt(xmomall*xmomall + ymomall*ymomall + zmomall*zmomall))


    angaccx = reduce_fn('+',angaccx)
    angaccy = reduce_fn('+',angaccy)
    angaccz = reduce_fn('+',angaccz)
    angxall = angx + angaccx
    angyall = angy + angaccy
    angzall = angz + angaccz
    call ev_data_update(ev_data,'angall',sqrt(angxall*angxall + angyall*angyall + angzall*angzall))

 endif

 if (track_mass) then
    accretedmass = ev_get_value('accretedmas')
    call ev_data_update(ev_data,'eacc',accretedmass/accradius1) ! total accretion energy
 endif
 if (track_lum) totlum = ev_get_value('tot lum')

 return
end subroutine compute_energies
!----------------------------------------------------------------
!+
!  calculates the centre of mass for use in rotational energy
!+
!----------------------------------------------------------------
subroutine get_erot_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 use centreofmass, only: get_centreofmass
 integer, intent(in) :: npart,nptmass
 real,    intent(in) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real :: xcom(3),vcom(3)

 if (calc_erot_com) then
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
    erot_com(1:3) = real(xcom(1:3))
    erot_com(4:6) = real(vcom(1:3))
 endif

end subroutine get_erot_com
!----------------------------------------------------------------
!+
!  calculates rotational energy
!+
!----------------------------------------------------------------
subroutine get_erot(xi,yi,zi,vxi,vyi,vzi,pmassi,erotxi,erotyi,erotzi)
 real, intent(in)  :: xi,yi,zi,vxi,vyi,vzi,pmassi
 real, intent(out) :: erotxi,erotyi,erotzi
 real              :: dx,dy,dz,dvx,dvy,dvz
 real              :: rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2
 !
 erotxi = 0.0
 erotyi = 0.0
 erotzi = 0.0

 dx  = xi  - erot_com(1)
 dy  = yi  - erot_com(2)
 dz  = zi  - erot_com(3)
 dvx = vxi - erot_com(4)
 dvy = vyi - erot_com(5)
 dvz = vzi - erot_com(6)

 rcrossvx = (dy*dvz - dz*dvy)
 rcrossvy = (dz*dvx - dx*dvz)
 rcrossvz = (dx*dvy - dy*dvx)

 radxy2 = dx*dx + dy*dy
 radyz2 = dy*dy + dz*dz
 radxz2 = dx*dx + dz*dz

 if (radyz2 > 0.) erotxi = pmassi*rcrossvx*rcrossvx/radyz2
 if (radxz2 > 0.) erotyi = pmassi*rcrossvy*rcrossvy/radxz2
 if (radxy2 > 0.) erotzi = pmassi*rcrossvz*rcrossvz/radxy2
 !
end subroutine get_erot
!----------------------------------------------------------------
!+
!  update the ev_data_array
!+
!----------------------------------------------------------------
subroutine ev_data_update(evdata,evtag,val)
 real,              intent(in)    :: val
 real,              intent(inout) :: evdata(0:inumev)
 character(len=*),  intent(in)    :: evtag
 integer                          :: iarray(inumev)
 integer                          :: ival,jval,j_actual
 character(len= 3)                :: cmd
 character(len=11)                :: evtag0

 write(evtag0,'(a)') evtag
 iarray = index(evtag0,ev_tag)
 ival   = maxloc(iarray,1)
 if (iarray(ival) > 0) then
    cmd    = ev_action(ival)
    jval   = ev_istart(ival)
    if (index(cmd,'0') > 0) then
       evdata(jval) =     val
    endif
    if (index(cmd,'s') > 0) then
       evdata(jval) =     evdata(jval)+val
    endif
    if (index(cmd,'x') > 0) then
       j_actual         = jval + index(cmd,'x') - 1
       evdata(j_actual) = max(evdata(j_actual),val)
    endif
    if (index(cmd,'a') > 0) then
       j_actual         = jval + index(cmd,'a') - 1
       evdata(j_actual) =     evdata(j_actual)+val
    endif
    if (index(cmd,'n') > 0) then
       j_actual         = jval + index(cmd,'n') - 1
       evdata(j_actual) = min(evdata(j_actual),val)
    endif
 endif

end subroutine ev_data_update
!----------------------------------------------------------------
!+
!  initiallised the ev_data array to zero or huge (for min)
!+
!----------------------------------------------------------------
subroutine initialise_ev_data(evdata)
 real,    intent(inout) :: evdata(0:inumev)
 integer                :: i,jval,j_actual
 character(len=3)       :: cmd
 !
 evdata = 0.0
 do i = 1,iquantities
    jval = ev_istart(i)
    cmd  = ev_action(i)
    if (index(cmd,'x') > 0) then
       j_actual         = jval + index(cmd,'x') - 1
       evdata(j_actual) = -huge(evdata(j_actual  ))
    endif
    if (index(cmd,'n') > 0) then
       j_actual         = jval + index(cmd,'n') - 1
       evdata(j_actual) =  huge(evdata(j_actual  ))
    endif
 enddo
 !
end subroutine initialise_ev_data
!----------------------------------------------------------------
!+
!  combines the ev_data from the various threads
!+
!----------------------------------------------------------------
subroutine collate_ev_data(evdata_thread,evdata)
 real,            intent(in)    :: evdata_thread(0:inumev)
 real,            intent(inout) :: evdata(0:inumev)
 integer                        :: i,jval,j_actual
 character(len=3)               :: cmd
 !
 do i = 1,iquantities
    jval = ev_istart(i)
    cmd  = ev_action(i)
    if (index(cmd,'s') > 0) then
       evdata(jval    ) =     evdata(jval    )+ evdata_thread(jval    )
    endif
    if (index(cmd,'x') > 0) then
       j_actual         = jval + index(cmd,'x') - 1
       evdata(j_actual) = max(evdata(j_actual), evdata_thread(j_actual))
    endif
    if (index(cmd,'a') > 0) then
       j_actual         = jval + index(cmd,'a') - 1
       evdata(j_actual) =     evdata(j_actual)+ evdata_thread(j_actual)
    endif
    if (index(cmd,'n') > 0) then
       j_actual         = jval + index(cmd,'n') - 1
       evdata(j_actual) = min(evdata(j_actual), evdata_thread(j_actual))
    endif
 enddo
 !
end subroutine collate_ev_data
!----------------------------------------------------------------
!+
!  Performs final generic housekeeping on the ev_data array
!+
!----------------------------------------------------------------
subroutine finalise_ev_data(evdata,dnptot)
 use mpiutils, only:reduceall_mpi
 real,            intent(inout) :: evdata(0:inumev)
 real,            intent(in)    :: dnptot
 integer                        :: i,jval,j_actual
 character(len=3)               :: cmd

 !
 do i = 1,iquantities
    jval = ev_istart(i)
    cmd  = ev_action(i)
    if (index(cmd,'s') > 0) then
       evdata(jval    ) =  reduce_fn('+',evdata(jval  ))
    endif
    if (index(cmd,'x') > 0) then
       j_actual         = jval + index(cmd,'x') - 1
       evdata(j_actual) = reduce_fn('max',evdata(j_actual))
    endif
    if (index(cmd,'a') > 0) then
       j_actual         = jval + index(cmd,'a') - 1
       evdata(j_actual) = reduce_fn('+',  evdata(j_actual))*dnptot
    endif
    if (index(cmd,'n') > 0) then
       j_actual         = jval + index(cmd,'n') - 1
       evdata(j_actual) = reduce_fn('min',evdata(j_actual))
    endif
 enddo
 !
end subroutine finalise_ev_data
!----------------------------------------------------------------
!+
!  Scale value
!+
!----------------------------------------------------------------
subroutine ev_data_correction(evdata,evtag,scalar_in,evcmd)
 real,                        intent(inout) :: evdata(0:inumev)
 real,                        intent(in)    :: scalar_in
 character(len=*),            intent(in)    :: evtag
 character(len= 1), optional, intent(in)    :: evcmd
 character(len=11)                          :: evtag0
 integer                                    :: iarray(inumev)
 integer                                    :: ival,jval

 write(evtag0,'(a)') evtag
 iarray = index(evtag0,ev_tag)
 ival   = maxloc(iarray,1)
 if (iarray(ival) > 0) then
    jval = ev_istart(ival)
    if (present(evcmd)) jval = jval + index(ev_action(jval),evcmd) - 1
    evdata(jval) = evdata(jval)*scalar_in
 endif

end subroutine ev_data_correction
!----------------------------------------------------------------
!+
!  Return a single array value
!+
!----------------------------------------------------------------
real function ev_get_value(evtag,evcmd)
 character(len=*),            intent(in) :: evtag
 character(len= 1), optional, intent(in) :: evcmd
 character(len=11)                       :: evtag0
 integer                                 :: iarray(inumev)
 integer                                 :: ival,jval

 write(evtag0,'(a)') evtag
 iarray = index(evtag0,ev_tag)
 ival   = maxloc(iarray,1)
 if (iarray(ival) > 0) then
    jval = ev_istart(ival)
    if (present(evcmd)) jval = jval + index(ev_action(jval),evcmd) - 1
    ev_get_value = ev_data(jval)
 else
    ev_get_value = 0.0
 endif

end function ev_get_value
!----------------------------------------------------------------
end module energies
