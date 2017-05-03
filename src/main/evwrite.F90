!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: evwrite
!
!  DESCRIPTION:
!  Calculates conserved quantities etc and writes to .ev file
!  Also writes log output
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, energies, extern_binary, externalforces,
!    io, nicil, options, part, ptmass, units, viscosity
!+
!--------------------------------------------------------------------------
module evwrite
 use options,  only: iexternalforce
 use energies, only: inumev,ielements,ievU,ievI,iev0,ievS,ievA,ievX,ievN,ev_data,ev_action
 use energies, only: itime,iekin,ietherm,iemag,iepot,ietot,itotmom,iangtot,irhoX,irhoA,idt,ientrop, &
                     irms,&
                     idustX,idustA,ibdyX,ibdyA,istarX,istarA,idmX,idmA,iblgX,iblgA,igasX,igasA, &
                     ialphaX,idivBX,idivBA,ihdivBX,ihdivBA,ibetaX,ibetaA,ibetaN, &
                     itX,itA,itN,ietaFX,ietaFA,ietaFN,iohmX,iohmA,iohmN,iohmfX,iohmfA,iohmfN, &
                     ihallX,ihallA,ihallN,iahallX,iahallA,iahallN, &
                     ihallfX,ihallfA,ihallfN,iahallfX,iahallfA,iahallfN,&
                     iambiX,iambiA,iambiN,iambifX,iambifA,iambifN, &
                     ivelX,ivelA,ivelN,ivionX,ivionA,ivionN,ivdftX,ivdftA,ivdftN, &
                     inenX,inenA,inenN,ineX,ineA,innX,innA,inihrX,inihrA,inimrX,inimrA, &
                     ingnX,ingnA,ingX,ingA,ingpX,ingpA,inistX,inistA,inidtX,inidtA, &
                     inhX,inhA,inheX,inheA,innaX,innaA,inmgX,inmgA,inkX,inkA, &
                     inhedX,inhedA,innadX,innadA,inmgdX,inmgdA,inkdX,inkdA, &
                     idtgX,idtgA,idtgN,itsN,itsA,iviscX,iviscA,iviscN,&
                     imomall,iangall,imacc1,imacc2,iamass,ieacc,ilum,ierot,ierotx,ieroty,ierotz
 use energies, only: erot_com
 use energies, only: gas_only,track_mass,track_lum
 use externalforces, only: iext_binary

 implicit none
 public                    :: init_evfile, write_evfile, write_evlog
 private                   :: fill_one_label,fill_xan_label,fill_rho_label

 integer,          private :: ievfile
 character(len=19),private :: ev_label(inumev)

 private

contains

!----------------------------------------------------------------
!+
!  opens the .ev file for output
!+
!----------------------------------------------------------------
subroutine init_evfile(iunit,evfile)
 use io,        only: id,master,warning
 use dim,       only: maxalpha,maxp,mhd,mhd_nonideal,use_dustfrac,calc_erot,lightcurve
 use options,   only: ishock_heating,ipdv_heating
 use part,      only: igas,idust,iboundary,istar,idarkmatter,ibulge,nptmass,npartoftype
 use ptmass,    only: icreate_sinks
 use nicil,     only: use_ohm,use_hall,use_ambi,ion_rays,ion_thermal,nelements,nlevels
 use viscosity, only: irealvisc
 integer,            intent(in) :: iunit
 character(len=  *), intent(in) :: evfile
 character(len= 27)             :: ev_fmt
 integer                        :: i
 !
 !--Initialise additional variables
 !
 erot_com  = 0.0
 ev_action = ievI
 ievfile   = iunit

 write(ev_fmt,'(a)') "(1x,'[',i2.2,1x,a11,']',2x)"
 !
 !--Define all the variables to be included in the .ev file.
 !
 i = 1
 call fill_one_label(ev_fmt,'time',     i,itime,  iev0)
 call fill_one_label(ev_fmt,'ekin',     i,iekin,  iev0)
 call fill_one_label(ev_fmt,'etherm',   i,ietherm,iev0)
 call fill_one_label(ev_fmt,'emag',     i,iemag,  iev0)
 call fill_one_label(ev_fmt,'epot',     i,iepot,  iev0)
 call fill_one_label(ev_fmt,'etot',     i,ietot,  iev0)
 call fill_one_label(ev_fmt,'totmom',   i,itotmom,iev0)
 call fill_one_label(ev_fmt,'angtot',   i,iangtot,iev0)
 call fill_one_label(ev_fmt,'rhomax',   i,irhoX,  ievX)
 call fill_one_label(ev_fmt,'rhomean',  i,irhoA,  ievA)
 call fill_one_label(ev_fmt,'dt',       i,idt,    iev0)
 call fill_one_label(ev_fmt,'totentrop',i,ientrop,ievS)
 call fill_one_label(ev_fmt,'rmsmach',  i,irms,   iev0)
 gas_only = .true.
 call fill_rho_label(ev_fmt,'rho dust',i,idustX,idustA,npartoftype(idust),      gas_only)
 call fill_rho_label(ev_fmt,'rho bdy', i,ibdyX, ibdyA, npartoftype(iboundary),  gas_only)
 call fill_rho_label(ev_fmt,'rho star',i,istarX,istarA,npartoftype(istar),      gas_only)
 call fill_rho_label(ev_fmt,'rho dm',  i,idmX,  idmA,  npartoftype(idarkmatter),gas_only)
 call fill_rho_label(ev_fmt,'rho blg', i,iblgX, iblgA, npartoftype(ibulge),     gas_only)
 if (.not.gas_only) then
    call fill_rho_label(ev_fmt,'rho gas',i,igasX,igasA,npartoftype(igas),       gas_only)
 endif
 if (maxalpha==maxp) then
    call fill_one_label(ev_fmt,'alpha max',i,ialphaX,ievX)
 endif
 if ( mhd ) then
    call fill_xan_label(ev_fmt,'divB',   i,idivBX, idivBA       )
    call fill_xan_label(ev_fmt,'hdivB/B',i,ihdivBX,ihdivBA      )
    call fill_xan_label(ev_fmt,'beta',   i,ibetaX, ibetaA,ibetaN)
    if (mhd_nonideal) then
       call fill_xan_label(ev_fmt,'temp',        i,itX,     itA   ,  itN     )
       call fill_xan_label(ev_fmt,'eta_ar',      i,ietaFX,  ietaFA,  ietaFN  )
       if (use_ohm) then
          call fill_xan_label(ev_fmt,'eta_o',    i,iohmX,   iohmA,   iohmN   )
          call fill_xan_label(ev_fmt,'eta_o/art',i,iohmfX,  iohmfA,  iohmfN  )
       endif
       if (use_hall) then
          call fill_xan_label(ev_fmt,'eta_h',    i,ihallX,  ihallA,  ihallN  )
          call fill_xan_label(ev_fmt,'|eta_h|',  i,iahallX, iahallA, iahallN )
          call fill_xan_label(ev_fmt,'eta_h/art',i,ihallfX, ihallfA, ihallfN )
          call fill_xan_label(ev_fmt,'|e_h|/art',i,iahallfX,iahallfA,iahallfN)
       endif
       if (use_ambi) then
          call fill_xan_label(ev_fmt,'eta_a',    i,iambiX,  iambiA,  iambiN  )
          call fill_xan_label(ev_fmt,'eta_a/art',i,iambifX, iambifA, iambifN )
          call fill_xan_label(ev_fmt,'velocity', i,ivelX,   ivelA,   ivelN   )
          call fill_xan_label(ev_fmt,'v_ion',    i,ivionX,  ivionA,  ivionN  )
          call fill_xan_label(ev_fmt,'v_drift',  i,ivdftX,  ivdftA,  ivdftN  )
       endif
       call fill_xan_label(ev_fmt,'n_e/n',       i,inenX,   inenA,   inenN   )
       call fill_xan_label(ev_fmt,'n_e',         i,ineX,    ineA             )
       call fill_xan_label(ev_fmt,'n_n',         i,innX,    innA             )
       if (ion_rays) then
          call fill_xan_label(ev_fmt,'n_ihR',    i,inihrX,  inihrA           )
          call fill_xan_label(ev_fmt,'n_imR',    i,inimrX,  inimrA           )
          call fill_xan_label(ev_fmt,'n_g(Z=-1)',i,ingnX,   ingnA            )
          call fill_xan_label(ev_fmt,'n_g(Z= 0)',i,ingX,    ingA             )
          call fill_xan_label(ev_fmt,'n_g(Z=+1)',i,ingpX,   ingpA            )
       endif
       if (ion_thermal) then
          call fill_xan_label(ev_fmt,   'n_isT', i,inistX,  inistA           )
          if (nlevels>=2) then
             call fill_xan_label(ev_fmt,'n_idT', i,inidtX,  inidtA           )
          endif
          if (nelements>=2) then
             call fill_xan_label(ev_fmt,'n_H+',  i,inhX,    inhA             )
             call fill_xan_label(ev_fmt,'n_He+', i,inheX,   inheA            )
          endif
          if (nelements>=5) then
             call fill_xan_label(ev_fmt,'n_Na+', i,innaX,   innaA            )
             call fill_xan_label(ev_fmt,'n_Mg+', i,inmgX,   inmgA            )
             call fill_xan_label(ev_fmt,'n_K+',  i,inkX,    inkA             )
          endif
          if (nlevels>=2) then
             if (nelements>=2) then
                call fill_xan_label(ev_fmt,'n_He++',i,inhedX,   inhedA       )
             endif
             if (nelements>=5) then
                call fill_xan_label(ev_fmt,'n_Na++',i,innadX,   innadA       )
                call fill_xan_label(ev_fmt,'n_Mg++',i,inmgdX,   inmgdA       )
                call fill_xan_label(ev_fmt,'n_K++', i,inkdX,    inkdA        )
             endif
          endif
       endif
    endif
 endif
 if (use_dustfrac) then
    call fill_xan_label(ev_fmt,'dust/gas',i,idtgX,idtgA,idtgN)
    call fill_xan_label(ev_fmt,'t_s',     i,iA=itsA,iN=itsN)
 endif
 if (iexternalforce > 0) then
    call fill_one_label(ev_fmt,'totmomall',i,imomall,iev0)
    call fill_one_label(ev_fmt,'angall',   i,iangall,iev0)
    if (iexternalforce==iext_binary) then
       call fill_one_label(ev_fmt,'Macc sink 1',i,imacc1,iev0)
       call fill_one_label(ev_fmt,'Macc sink 2',i,imacc2,iev0)
    endif
 endif
 if (iexternalforce>0 .or. nptmass > 0 .or. icreate_sinks > 0) then
    call fill_one_label(ev_fmt,'accretedmas',i,iamass,iev0)
    call fill_one_label(ev_fmt,'eacc',       i,ieacc, iev0)
    track_mass     = .true.
 else
    track_mass     = .false.
 endif
 if (ishock_heating==0 .or. ipdv_heating==0 .or. lightcurve) then
    call fill_one_label(ev_fmt,'tot lum',i,ilum,iev0)
    track_lum      = .true.
 else
    track_lum      = .false.
 endif
 if (calc_erot) then
    call fill_one_label(ev_fmt,'erot',  i,ierot, iev0)
    call fill_one_label(ev_fmt,'erot_x',i,ierotx,iev0)
    call fill_one_label(ev_fmt,'erot_y',i,ieroty,iev0)
    call fill_one_label(ev_fmt,'erot_z',i,ierotz,iev0)
 endif
 if (irealvisc /= 0) then
    call fill_xan_label(ev_fmt,'visc_rat',i,iviscX, iviscA, iviscN )
 endif
 ielements = i - 1 ! The number of columns to be calculates
 !
 !--all threads do above, but only master writes file
 !
 if (id == master) then
    !
    !--open the file for output
    !
    open(unit=ievfile,file=evfile,form='formatted',status='replace')
    !
    !--write a header line
    !
    write(ev_fmt,'(a,I3,a)') '(',ielements+1,'a)'
    write(ievfile,ev_fmt)'#',ev_label(1:ielements)
 endif

end subroutine init_evfile
!
!----------------------------------------------------------------
!+
!  creates a single label and identifies the required action
!+
!----------------------------------------------------------------
subroutine fill_one_label(ev_fmt,label,i,i0,iev)
 implicit none
 integer(kind=1),  intent(in)    :: iev
 integer,          intent(inout) :: i,i0
 character(len=*), intent(in)    :: label,ev_fmt
 !
 i0 = i
 i  = i + 1
 if (i0 > 99) then
    write(ev_label(i0),ev_fmt) 100-i0,trim(label)
 else
    write(ev_label(i0),ev_fmt) i0,trim(label)
 endif
 ev_action(i0) = iev
 !
end subroutine fill_one_label
!----------------------------------------------------------------
!+
!  creates up to three lables per input value: Max, Ave, Min
!+
!----------------------------------------------------------------
subroutine fill_xan_label(ev_fmt,label,i,iX,iA,iN)
 implicit none
 integer,          intent(inout) :: i
 integer,optional, intent(inout) :: iX,iA,iN
 character(len=*), intent(in)    :: label,ev_fmt
 character(len=3)                :: lX,lA,lN
 character(len=11)               :: label0
 !
 if (len(label)<=7) then
    lX = "max"
    lA = "ave"
    lN = "min"
 else
    lX = "X"
    lA = "A"
    lN = "N"
 endif
 !
 if (present(iX)) then
    write(label0,'(a,1x,a)')trim(label),trim(lX); call fill_one_label(ev_fmt,label0,i,iX,ievX)
 endif
 if (present(iA)) then
    write(label0,'(a,1x,a)')trim(label),trim(lA); call fill_one_label(ev_fmt,label0,i,iA,ievA)
 endif
 if (present(iN)) then
    write(label0,'(a,1x,a)')trim(label),trim(lN); call fill_one_label(ev_fmt,label0,i,iN,ievN)
 endif
 !
end subroutine fill_xan_label
!
!----------------------------------------------------------------
!+
!  creates two entries (max, ave) per input density
!+
!----------------------------------------------------------------
subroutine fill_rho_label(ev_fmt,label,i,iX,iA,npartoftypei,gas_only)
 implicit none
 integer,          intent(in)    :: npartoftypei
 integer,          intent(inout) :: i,iX,iA
 character(len=*), intent(in)    :: label,ev_fmt
 logical,          intent(inout) :: gas_only
 character(len=11)               :: label0
 !
 if (npartoftypei > 0) then
    write(label0,'(2a)')trim(label),' X'; call fill_one_label(ev_fmt,label0,i,iX,ievX)
    write(label0,'(2a)')trim(label),' A'; call fill_one_label(ev_fmt,label0,i,iA,ievA)
    gas_only = .false.
 else
    iX = 0
    iA = 0
 endif
 !
end subroutine fill_rho_label
!----------------------------------------------------------------
!+
!  calculates total energy, etc, and writes line to .ev file
!+
!----------------------------------------------------------------
subroutine write_evfile(t,dt)
 use energies,      only:compute_energies,ev_data
 use io,            only:id,master
 use options,       only:iexternalforce
 use extern_binary, only:accretedmass1,accretedmass2
 real, intent(in)  :: t, dt
 character(len=35) :: ev_format

 call compute_energies(t)

 if (id==master) then
    !
    !--fill in additional details that are not calculated in energies.f
    ev_data(idt)   = dt
    if (iexternalforce==iext_binary) then
       ev_data(imacc1) = accretedmass1
       ev_data(imacc2) = accretedmass2
    endif
    !
    !--write line to .ev file (should correspond to header, below)
    !
    write(ev_format,'(a,I3,a)')"(",ielements,"(1pe18.10,1x))"
    write(ievfile,ev_format) ev_data(1:ielements)
    call flush(ievfile)
 endif

 return
end subroutine write_evfile
!----------------------------------------------------------------
!+
!  Writes nicely formatted output to the log file/screen
!  Must be called *after* a call to compute energies has been
!  performed
!+
!----------------------------------------------------------------
subroutine write_evlog(iprint)
 use dim,       only:maxp,maxalpha,mhd,maxvxyzu,periodic,mhd_nonideal,use_dustfrac
 use energies,  only:ev_data,ekin,etherm,emag,epot,etot,rmsmach,vrms,accretedmass,mdust,mgas
 use viscosity, only:irealvisc,shearparam
 use boundary,  only:dxbound,dybound,dzbound
 use units,     only:unit_density
 integer, intent(in) :: iprint
 character(len=120)  :: string

 ! this is currently broken - disabled temporarily
 return

 write(iprint,"(1x,3('E',a,'=',es10.3,', '),('E',a,'=',es10.3))") &
      'tot',etot,'kin',ekin,'therm',etherm,'pot',epot
 if (mhd)        write(iprint,"(1x,('E',a,'=',es10.3))") 'mag',emag
 if (track_mass) write(iprint,"(1x,('E',a,'=',es10.3))") 'acc',ev_data(ieacc)
 write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3))") &
      'Linm',ev_data(itotmom),'Angm',ev_data(iangtot)
 if (iexternalforce > 0) then
    if (abs(ev_data(iangall)-ev_data(iangtot)) > tiny(0.)) then
       write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3),a)") &
       'Linm',ev_data(imomall),'Angm',ev_data(iangall),' [including accreted particles]'
    endif
 endif

 write(iprint,"(1x,a,'(max)=',es10.3,' (mean)=',es10.3,' (max)=',es10.3,a)") &
      'density  ',ev_data(irhoX:irhoA),ev_data(irhoX)*unit_density,' g/cm^3'

 if (use_dustfrac) then
    write(iprint,"(1x,a,'(max)=',es10.3,1x,'(mean)=',es10.3,1x,'(min)=',es10.3)") &
         'dust2gas ',ev_data(idtgX:idtgN)
    write(iprint,"(3x,a,'(mean)=',es10.3,1x,'(min)=',es10.3)") 't_stop ',ev_data(itsA:itsN)
    write(iprint,"(1x,'Mgas = ',es10.3,', Mdust = ',es10.3)") mgas,mdust
 endif

 write(iprint,"(1x,1(a,'=',es10.3))") &
      'Accreted mass',accretedmass

 string = ''
 if (maxalpha==maxp) then
    write(string,"(a,'(max)=',es10.3)") ' alpha',ev_data(ialphaX)
 endif
 if (len_trim(string) > 0) write(iprint,"(a)") trim(string)

 if (irealvisc /= 0) then
    if (periodic) then
       if (irealvisc==1) then
          write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3))") &
             'RMS Mach #',rmsmach,'Reynolds # ',vrms*min(dxbound,dybound,dzbound)/shearparam
       endif
    endif
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),('(mean)=',es10.3),(' (min)=',es10.3))") &
         'Ratio of physical-to-art. visc',ev_data(iviscX:iviscN)
 else
    write(iprint,"(1x,1(a,'=',es10.3))") &
        'RMS Mach #',rmsmach
 endif

 if (mhd) then
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),(a,'(mean)=',es10.3))") &
      'div B ',ev_data(idivBX),'div B ',ev_data(idivBA)
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),(a,'(mean)=',es10.3))") &
      'h|div B|/B ',ev_data(ihdivBX),'h|div B|/B ',ev_data(ihdivBA)
 endif
 write(iprint,"(/)")

 return
end subroutine write_evlog
!----------------------------------------------------------------

end module evwrite
