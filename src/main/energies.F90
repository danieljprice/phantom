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
!#define reduce_fn(a,b) b
#define reduce_fn(a,b) reduceall_mpi(a,b)
module energies
 use dim, only: calc_erot, calc_erot_com
 implicit none

 public :: compute_energies,get_erot_com

 integer, public :: itime,iekin,ietherm,iemag,iepot,ietot,itotmom,iangtot,irhoX,irhoA, &
                    idt,ientrop,irms, &
                    idustX,idustA,ibdyX,ibdyA,istarX,istarA,idmX,idmA,iblgX,iblgA,igasX,igasA, &
                    ialphaX,idivBX,idivBA,ihdivBX,ihdivBA,ibetaX,ibetaA,ibetaN, &
                    itX,itA,itN,ietaFX,ietaFA,ietaFN,iohmX,iohmA,iohmN,iohmfX,iohmfA,iohmfN, &
                    ihallX,ihallA,ihallN,iahallX,iahallA,iahallN, &
                    ihallfX,ihallfA,ihallfN,iahallfX,iahallfA,iahallfN, &
                    iambiX,iambiA,iambiN,iambifX,iambifA,iambifN, &
                    ivelX,ivelA,ivelN,ivionX,ivionA,ivionN,ivdftX,ivdftA,ivdftN, &
                    inenX,inenA,inenN,ineX,ineA,innX,innA,inihrX,inihrA,inimrX,inimrA, &
                    ingnX,ingnA,ingX,ingA,ingpX,ingpA,inistX,inistA,inidtX,inidtA, &
                    inhX,inhA,inheX,inheA,innaX,innaA,inmgX,inmgA,inkX,inkA, &
                    inhedX,inhedA,innadX,innadA,inmgdX,inmgdA,inkdX,inkdA, &
                    idtgX,idtgA,idtgN,itsN,itsA,iviscX,iviscA,iviscN, &
                    imomall,iangall,imacc1,imacc2,iamass,ieacc,ilum,ierot,ierotx,ieroty,ierotz
 logical, public :: gas_only,track_mass,track_lum

 real,            public    :: ekin,etherm,emag,epot,etot,totmom,angtot
 real,            public    :: vrms,rmsmach,accretedmass,mdust,mgas
 real,            public    :: xmom,ymom,zmom
 real,            public    :: totlum
 integer,         public    :: ielements
 integer,         parameter :: inumev = 150
 ! Actions to be taken as instructed by ev_action
 integer(kind=1), parameter :: ievU   = -2  ! uninitialised entry to be included in .ev
 integer(kind=1), parameter :: ievI   = -1  ! initial value of ev_action
 integer(kind=1), parameter :: iev0   =  0  ! take no action
 integer(kind=1), parameter :: ievS   =  1  ! sum      the value
 integer(kind=1), parameter :: ievA   =  2  ! average  the value
 integer(kind=1), parameter :: ievN   =  3  ! minimise the value
 integer(kind=1), parameter :: ievX   =  4  ! maximise the value
 ! Additional parameters
 integer(kind=1), public :: ev_action(inumev)
 real,            public :: ev_data(0:inumev),erot_com(6)

contains

!----------------------------------------------------------------
!+
!  Subroutine to compute global quantities on the particles
!+
!----------------------------------------------------------------
subroutine compute_energies(t)
 use dim,  only:maxp,maxvxyzu,maxalpha,maxtypes,use_dustfrac,mhd_nonideal,lightcurve
 use part, only:rhoh,xyzh,vxyzu,massoftype,npart,maxphase,iphase, &
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
 real    :: shearparam_art,shearparam_phys,ratio_phys_to_av
 real    :: gasfrac,dustfraci,dust_to_gas
 real    :: temperature,etaart,etaart1,etaohm,etahall,etaambi,vion,vdrift
 real    :: vioni(3),data_out(17+nelements_max*nlevels-3)
 real    :: erotx,eroty,erotz,erotxi,erotyi,erotzi,fdum(3)
 integer :: i,j,itype,ierr,naccreted
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
 vrms = 0.
 rmsmach = 0.
 np = 0
 naccreted = 0
 accretedmass = 0.
 xmomacc = 0.
 ymomacc = 0.
 zmomacc = 0.
 angaccx = 0.
 angaccy = 0.
 angaccz = 0.
 erotx   = 0.
 eroty   = 0.
 erotz   = 0.
 mgas    = 0.
 mdust   = 0.
 if (maxalpha==maxp) then
    alphai = 0.
 else
    alphai = alpha
 endif
 totlum      = 0.
 ionfrac_eta = 0.
 np_rho      = 0
 call initialise_ev_data(ielements,ev_action,ev_data)
!
!$omp parallel default(none) &
!$omp shared(xyzh,vxyzu,iexternalforce,npart,t,id) &
!$omp shared(alphaind,massoftype,irealvisc) &
!$omp shared(ieos,gamma,nptmass,xyzmh_ptmass,vxyz_ptmass) &
!$omp shared(Bxyz,Bevol,divBsymm,divcurlB,alphaB,iphase,poten,dustfrac) &
!$omp shared(use_ohm,use_hall,use_ambi,ion_rays,ion_thermal,nelements,n_R,n_electronT,ionfrac_eta) &
!$omp shared(ielements,ev_data,ev_action,np_rho,erot_com) &
!$omp shared(calc_erot,gas_only,irhoX,irhoA,ialphaX,iviscX,iviscA,iviscN) &
!$omp shared(idustX,idustA,ibdyX,ibdyA,istarX,istarA,idmX,idmA,iblgX,iblgA,igasX,igasA) &
!$omp shared(idtgX,idtgA,idtgN,itsA,itsN,ientrop) &
!$omp shared(idivBX,idivBA,ihdivBX,ihdivBA,ibetaX,ibetaA,ibetaN) &
!$omp shared(itX,itA,itN,ietaFX,ietaFA,ietaFN,iohmX,iohmA,iohmN,iohmfX,iohmfA,iohmfN) &
!$omp shared(ihallX,ihallA,ihallN,iahallX,iahallA,iahallN) &
!$omp shared(ihallfX,ihallfA,ihallfN,iahallfX,iahallfA,iahallfN) &
!$omp shared(iambiX,iambiA,iambiN,iambifX,iambifA,iambifN) &
!$omp shared(ivelX,ivelA,ivelN,ivionX,ivionA,ivionN,ivdftX,ivdftA,ivdftN) &
!$omp shared(inenX,inenA,inenN,ineX,ineA,innX,innA,inihrX,inihrA,inimrX,inimrA) &
!$omp shared(ingnX,ingnA,ingX,ingA,ingpX,ingpA,inistX,inistA,inidtX,inidtA) &
!$omp shared(inhX,inhA,inheX,inheA,innaX,innaA,inmgX,inmgA,inkX,inkA) &
!$omp shared(inhedX,inhedA,innadX,innadA,inmgdX,inmgdA,inkdX,inkdA) &
!$omp private(i,j,xi,yi,zi,hi,rhoi,vxi,vyi,vzi,Bxi,Byi,Bzi,epoti,vsigi,v2i) &
!$omp private(ponrhoi,spsoundi,B2i,dumx,dumy,dumz,acci,valfven2i,divBi,hdivBonBi) &
!$omp private(rho1i,shearparam_art,shearparam_phys,ratio_phys_to_av,betai) &
!$omp private(gasfrac,dustfraci,dust_to_gas) &
!$omp private(ierr,temperature,etaart,etaart1,etaohm,etahall,etaambi,vioni,vion,vdrift,data_out) &
!$omp private(erotxi,erotyi,erotzi,fdum) &
!$omp private(ev_data_thread,np_rho_thread) &
!$omp firstprivate(alphai,itype,pmassi) &
#ifdef DUST
!$omp shared(grainsize,graindens,idrag) &
!$omp private(tsi,iregime) &
#endif
#ifdef LIGHTCURVE
!$omp shared(luminosity) &
#endif
!$omp reduction(+:np,xmom,ymom,zmom,angx,angy,angz,erotx,eroty,erotz,mgas,mdust) &
!$omp reduction(+:naccreted,xmomacc,ymomacc,zmomacc,angaccx,angaccy,angaccz) &
!$omp reduction(+:ekin,etherm,emag,epot,vrms,rmsmach,accretedmass,totlum)
 call initialise_ev_data(ielements,ev_action,ev_data_thread)
 np_rho_thread = 0
!$omp do
 do i=1,npart
    !OK: skip calculating for particles on other procs
    !Necessary if xyzh not AllReduced, requires sums to be AllReduced
    !if (id  /=  ibelong(xyzh(:,i),id)) cycle

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
       call ev_update(ev_data_thread,rhoi,irhoX,irhoA)
       if (.not.gas_only) then
          select case(itype)
          case(igas)
             call ev_rhoupdate(ev_data_thread,rhoi,igasX, igasA, itype,np_rho_thread)
          case(idust)
             call ev_rhoupdate(ev_data_thread,rhoi,idustX,idustA,itype,np_rho_thread)
          case(iboundary)
             call ev_rhoupdate(ev_data_thread,rhoi,ibdyX, ibdyA, itype,np_rho_thread)
          case(istar)
             call ev_rhoupdate(ev_data_thread,rhoi,istarX,istarA,itype,np_rho_thread)
          case(idarkmatter)
             call ev_rhoupdate(ev_data_thread,rhoi,idmX,  idmA,  itype,np_rho_thread)
          case(ibulge)
             call ev_rhoupdate(ev_data_thread,rhoi,iblgX, iblgA, itype,np_rho_thread)
          end select
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
       vrms = vrms + v2i

       ! rotational energy around each axis through the Centre of mass
       ! note: centre of mass is updated only when dumpfiles are created
       if (calc_erot) then
          call get_erot(xi,yi,zi,vxi,vyi,vzi,pmassi,erotxi,erotyi,erotzi)
          erotx = erotx + erotxi
          eroty = eroty + erotyi
          erotz = erotz + erotzi
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
             call ev_update(ev_data_thread,dust_to_gas,idtgX,idtgA,idtgN)
             mgas        = mgas + pmassi*gasfrac
             mdust       = mdust + pmassi*dustfraci
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
          call ev_update(ev_data_thread,pmassi*ponrhoi*rhoi**(1.-gamma),iA=ientrop)

#ifdef DUST
          ! min and mean stopping time
          if (use_dustfrac) then
             call get_ts(idrag,grainsize,graindens,rhoi*gasfrac,rhoi*dustfraci,spsoundi,0.,tsi,iregime)
             call ev_update(ev_data_thread,tsi,iA=itsA,iN=itsN)
          endif
#endif

#ifdef LIGHTCURVE
          totlum = totlum + luminosity(i)
#endif

          ! rms mach number
          if (spsoundi > 0.) rmsmach = rmsmach + v2i/spsoundi**2

          ! max of dissipation parameters
          if (maxalpha==maxp) then
             alphai = alphaind(1,i)
             call ev_update(ev_data_thread,alphai,ialphaX)
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
             call ev_update(ev_data_thread,ratio_phys_to_av,iviscX,iviscA,iviscN)
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

             divBi     = abs(divBsymm(i))
             if (B2i > 0.) then
                hdivBonBi = hi*divBi/sqrt(B2i)
                betai     = 2.0*ponrhoi*rhoi/B2i ! plasma beta
             else
                hdivBonBi = 0.
                betai     = 0.
             endif
             call ev_update(ev_data_thread,divBi,    idivBX, idivBA )
             call ev_update(ev_data_thread,hdivBonBi,ihdivBX,ihdivBA)
             call ev_update(ev_data_thread,betai,    ibetaX, ibetaA, ibetaN)

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
                call ev_update(ev_data_thread,temperature,      itX,   itA   ,itN   )
                call ev_update(ev_data_thread,etaart,           ietaFX,ietaFA,ietaFN)
                if (use_ohm) then
                   call ev_update(ev_data_thread,etaohm,        iohmX, iohmA, iohmN )
                   call ev_update(ev_data_thread,etaohm*etaart1,iohmfX,iohmfA,iohmfN)
                endif
                if (use_hall) then
                   call ev_update(ev_data_thread,    etahall,         ihallX,  ihallA,  ihallN  )
                   call ev_update(ev_data_thread,abs(etahall),        iahallX, iahallA, iahallN )
                   call ev_update(ev_data_thread,    etahall*etaart1, ihallfX, ihallfA, ihallfN )
                   call ev_update(ev_data_thread,abs(etahall)*etaart1,iahallfX,iahallfA,iahallfN)
                endif
                if (use_ambi) then
                   vion   = sqrt( dot_product(vioni,vioni) )
                   vdrift = sqrt( (vioni(1)-vxi)**2 + (vioni(2)-vyi)**2 + (vioni(3)-vzi)**2)
                   call ev_update(ev_data_thread,etaambi,        iambiX, iambiA, iambiN )
                   call ev_update(ev_data_thread,etaambi*etaart1,iambifX,iambifA,iambifN)
                   call ev_update(ev_data_thread,sqrt(v2i),      ivelX,  ivelA,  ivelN  )
                   call ev_update(ev_data_thread,vion,           ivionX, ivionA, ivionN )
                   call ev_update(ev_data_thread,vdrift,         ivdftX, ivdftA, ivdftN )
                endif
                if (data_out(7) > 0.0) then
                   call ev_update(ev_data_thread,data_out(6)/data_out(7),inenX,inenA,inenN)
                   ionfrac_eta(1,i) = real(data_out(6)/data_out(7),kind=4)
                else
                   call ev_update(ev_data_thread,0.0,                    inenX,inenA,inenN)
                   ionfrac_eta(1,i) = 0.0
                endif
                ionfrac_eta(2,i) = real(etaohm, kind=4)  ! Save eta_OR for the dump file
                ionfrac_eta(3,i) = real(etahall,kind=4)  ! Save eta_HE for the dump file
                ionfrac_eta(4,i) = real(etaambi,kind=4)  ! Save eta_AD for the dump file
                call ev_update(ev_data_thread,      data_out( 6), ineX,   ineA             )
                call ev_update(ev_data_thread,      data_out( 7), innX,   innA             )
                if (ion_rays) then
                   call ev_update(ev_data_thread,   data_out( 8),inihrX,   inihrA          )
                   call ev_update(ev_data_thread,   data_out( 9),inimrX,   inimrA          )
                   call ev_update(ev_data_thread,   data_out(12),ingnX,    ingnA           )
                   call ev_update(ev_data_thread,   data_out(13),ingX,     ingA            )
                   call ev_update(ev_data_thread,   data_out(14),ingpX,    ingpA           )
                endif
                if (ion_thermal) then
                   call ev_update(ev_data_thread,   data_out(10),inistX,   inistA          )
                   call ev_update(ev_data_thread,   data_out(11),inidtX,   inidtA          )
                   j = 17
                   if (nelements>=2) then
                      call ev_update(ev_data_thread,data_out( j),inhX,    inhA             ); j=j+1
                      call ev_update(ev_data_thread,data_out( j),inheX,   inheA            ); j=j+1
                   endif
                   if (nelements>=5) then
                      call ev_update(ev_data_thread,data_out( j),innaX,   innaA            ); j=j+1
                      call ev_update(ev_data_thread,data_out( j),inmgX,   inmgA            ); j=j+1
                      call ev_update(ev_data_thread,data_out( j),inkX,    inkA             ); j=j+1
                   endif
                   if (nlevels>=2) then
                      if (nelements>=2) then
                         call ev_update(ev_data_thread,data_out( j),inhedX,   inhedA       ); j=j+1
                      endif
                      if (nelements>=5) then
                         call ev_update(ev_data_thread,data_out( j),innadX,   innadA       ); j=j+1
                         call ev_update(ev_data_thread,data_out( j),inmgdX,   inmgdA       ); j=j+1
                         call ev_update(ev_data_thread,data_out( j),inkdX,    inkdA        ); j=j+1
                      endif
                   endif
                endif
             endif
          endif
       endif isgas

    elseif (was_accreted(iexternalforce,hi)) then
!
!--count accretion onto fixed potentials (external forces) separately
!
       naccreted = naccreted + 1
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

       accretedmass = accretedmass + pmassi
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
    accretedmass = accretedmass + acci

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
       erotx = erotx + erotxi
       eroty = eroty + erotyi
       erotz = erotz + erotzi
    endif
 enddo
!$omp enddo
!$omp critical(collatedata)
 call collate_ev_data(ielements,ev_action,ev_data_thread,ev_data)
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
 !--Finalise the arrays
 call finalise_ev_data(ielements,ev_data,ev_action,dnptot)

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
 ev_data(itime)   = t
 ev_data(iekin)   = ekin
 ev_data(ietherm) = etherm
 ev_data(iemag)   = emag
 ev_data(iepot)   = epot
 ev_data(ietot)   = etot
 ev_data(itotmom) = totmom
 ev_data(iangtot) = angtot

 if (calc_erot) then
    erotx = 0.5*erotx
    eroty = 0.5*eroty
    erotz = 0.5*erotz
    ev_data(ierot)  = sqrt(erotx*erotx + eroty*eroty + erotz*erotz)
    ev_data(ierotx) = erotx
    ev_data(ieroty) = eroty
    ev_data(ierotz) = erotz
 endif

 if (use_dustfrac) then
    mgas  = reduce_fn('+',mgas)
    mdust = reduce_fn('+',mdust)
 endif

 if (.not. gas_only) then
    do i = 1,maxtypes
       np_rho(i) = reduce_fn('+',np_rho(i))
    enddo
    ! correct the average densities
    if (np_rho(igas)        > 0) ev_data(igasA)  = ev_data(igasA) *nptot/np_rho(igas)
    if (np_rho(idust)       > 0) ev_data(idustA) = ev_data(idustA)*nptot/np_rho(idust)
    if (np_rho(iboundary)   > 0) ev_data(ibdyA)  = ev_data(ibdyA) *nptot/np_rho(iboundary)
    if (np_rho(istar)       > 0) ev_data(istarA) = ev_data(istarA)*nptot/np_rho(istar)
    if (np_rho(idarkmatter) > 0) ev_data(idmA)   = ev_data(idmA)  *nptot/np_rho(idarkmatter)
    if (np_rho(ibulge)      > 0) ev_data(iblgA)  = ev_data(iblgA) *nptot/np_rho(ibulge)
 endif
 vrms      = sqrt(reduce_fn('+',vrms)*dnptot)
 rmsmach   = sqrt(reduce_fn('+',rmsmach)*dnptot)
 ev_data(irms) = rmsmach

 if (iexternalforce > 0) then
    xmomacc   = reduce_fn('+',xmomacc)
    ymomacc   = reduce_fn('+',ymomacc)
    zmomacc   = reduce_fn('+',zmomacc)

    xmomall   = xmom + xmomacc
    ymomall   = ymom + ymomacc
    zmomall   = zmom + zmomacc
    ev_data(itotmom) = sqrt(xmomall*xmomall + ymomall*ymomall + zmomall*zmomall)

    angaccx = reduce_fn('+',angaccx)
    angaccy = reduce_fn('+',angaccy)
    angaccz = reduce_fn('+',angaccz)
    angxall = angx + angaccx
    angyall = angy + angaccy
    angzall = angz + angaccz
    ev_data(iangall) = sqrt(angxall*angxall + angyall*angyall + angzall*angzall)
 endif

 if (track_mass) then
    ev_data(iamass) = reduce_fn('+',accretedmass)
    ev_data(ieacc)  = accretedmass/accradius1 ! total accretion energy
 endif
 totlum = reduce_fn('+',totlum)
 if (track_lum) ev_data(ilum) = totlum

 return
end subroutine compute_energies
!----------------------------------------------------------------
!+
!  calculates the centre of mass for use in rotational energy
!+
!----------------------------------------------------------------
subroutine get_erot_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 use centreofmass, only: get_centreofmass
 implicit none
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
 implicit none
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
!  initiallised the ev_data array to zero or huge (for min)
!+
!----------------------------------------------------------------
subroutine initialise_ev_data(ielements,evaction,evdata)
 implicit none
 integer, intent(in)            :: ielements
 integer(kind=1), intent(in)    :: evaction(:)
 real,            intent(inout) :: evdata(0:inumev)
 integer                        :: i
 !
 evdata = 0.0
 do i = 1,ielements
    if (evaction(i)==ievN) evdata(i) =  huge(evdata(i))
    if (evaction(i)==ievX) evdata(i) = -huge(evdata(i))
 enddo
 !
end subroutine initialise_ev_data
!----------------------------------------------------------------
!+
!  combines the ev_data from the various threads
!+
!----------------------------------------------------------------
subroutine collate_ev_data(ielements,evaction,evdata_thread,evdata)
 implicit none
 integer,         intent(in)    :: ielements
 integer(kind=1), intent(in)    :: evaction(:)
 real,            intent(in)    :: evdata_thread(0:inumev)
 real,            intent(inout) :: evdata(0:inumev)
 integer                        :: i
 !
 do i = 1,ielements
    select case(evaction(i))
    case(ievS,ievA)
       evdata(i) = evdata(i) + evdata_thread(i)
    case(ievX)
       evdata(i) = max(evdata(i),evdata_thread(i))
    case(ievN)
       evdata(i) = min(evdata(i),evdata_thread(i))
    end select
 enddo
 !
end subroutine collate_ev_data
!----------------------------------------------------------------
!+
!  update the ev_data_array
!  first subroutine: generic
!  second subroutine: for density of individual particle types
!+
!----------------------------------------------------------------
subroutine ev_update(evdata,val,iX,iA,iN)
 implicit none
 integer, optional, intent(in)    :: iX,iA,iN
 real,              intent(in)    :: val
 real,              intent(inout) :: evdata(0:inumev)
 !
 if (present(iX)) evdata(iX) = max(evdata(iX),val)
 if (present(iA)) evdata(iA) =     evdata(iA)+val
 if (present(iN)) evdata(iN) = min(evdata(iN),val)
 !
end subroutine ev_update
!
subroutine ev_rhoupdate(evdata,rhoi,iX,iA,itype,nprho)
 implicit none
 integer, intent(in)    :: iX,iA,itype
 real,    intent(in)    :: rhoi
 real,    intent(inout) :: evdata(0:inumev)
 integer(kind=8), intent(inout) :: nprho(:)
 !
 evdata(iX)   = max(evdata(iX),rhoi)
 evdata(iA)   = evdata(iA)+rhoi
 nprho(itype) = nprho(itype) + 1
 !
end subroutine ev_rhoupdate
!----------------------------------------------------------------
!+
!  Performs final generic housekeeping on the ev_data array
!+
!----------------------------------------------------------------
subroutine finalise_ev_data(ielements,evdata,evaction,dnptot)
 use mpiutils, only:reduceall_mpi
 implicit none
 integer,         intent(in)    :: ielements
 integer(kind=1), intent(in)    :: evaction(:)
 real,            intent(inout) :: evdata(0:inumev)
 real,            intent(in)    :: dnptot
 integer                        :: i
 !
 do i = 1,ielements
    select case(evaction(i))
    case(ievS)
       evdata(i) = reduce_fn('+',  evdata(i))
    case(ievA)
       evdata(i) = reduce_fn('+',  evdata(i))*dnptot
    case(ievN)
       evdata(i) = reduce_fn('min',evdata(i))
    case(ievX)
       evdata(i) = reduce_fn('max',evdata(i))
    end select
    if ( evdata(i) >  0.01*huge(evdata(i)) .or. &
         evdata(i) < -0.01*huge(evdata(i))) evdata(i) = 0.0
 enddo
 !
end subroutine finalise_ev_data
!----------------------------------------------------------------
end module energies
