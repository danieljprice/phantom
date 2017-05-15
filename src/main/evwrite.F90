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
!  Calculates conserved quantities etc and writes to .ev file;
!  Also writes log output
!  To Developer: Adding values to the .ev file is a two step process.
!     In the init_evfile subroutine in evwrite.F90, add the following command:
!        call fill_ev_label(ev_fmt,ev_tag,action,i,j)
!     and in compute_energies subroutine in energies.F90, add the following
!     command:
!        call ev_data_update(ev_data_thread,ev_tag,value)
!     where
!        ev_fmt,ev_data_thread,i,j: pre-defined quantities to included verbatim
!        ev_tag: a string to identify the quantity (e.g. 'c_s' for sound speed)
!        ev_value: the value of the quantity for particle i (e.g., spsoundi
!                  for sound speed)
!        action: a string identifying what action(s) you would like performed
!                on the quantity.  The available options are
!           0: no action taken (e.g. for time)
!           s: sum quantity (e.g. for entropy)
!           x: print the maximum quantity
!           a: print the average (mean) quantity
!           n: print the minimum quantity
!        where any or all of x,a,n can be used as a single action
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
 use io,             only: fatal
 use options,        only: iexternalforce
 use externalforces, only: iext_binary
 use energies,       only: inumev,ielements,iquantities,ev_tag,ev_action,ev_istart,ev_data
 use energies,       only: erot_com,gas_only,track_mass,track_lum

 implicit none
 public                    :: init_evfile, write_evfile, write_evlog
 private                   :: fill_ev_tag, fill_ev_header

 integer,          private :: ievfile
 character(len=19),private :: ev_label(inumev)  ! to make the header for the .ev file

 private

contains

!----------------------------------------------------------------
!+
!  opens the .ev file for output
!+
!----------------------------------------------------------------
subroutine init_evfile(iunit,evfile)
 use io,        only: id,master,warning
 use dim,       only: maxtypes,maxalpha,maxp,mhd,mhd_nonideal,use_dustfrac,calc_erot,lightcurve
 use options,   only: ishock_heating,ipdv_heating
 use part,      only: igas,idust,iboundary,istar,idarkmatter,ibulge,nptmass,npartoftype
 use ptmass,    only: icreate_sinks
 use nicil,     only: use_ohm,use_hall,use_ambi,ion_rays,ion_thermal
 use viscosity, only: irealvisc
 integer,            intent(in) :: iunit
 character(len=  *), intent(in) :: evfile
 character(len= 27)             :: ev_fmt
 integer                        :: i,j
 !
 !--Initialise additional variables
 !
 erot_com  = 0.0
 gas_only  = .true.
 do i = 2,maxtypes
    if (npartoftype(i) > 0) gas_only = .false.
 enddo
 write(ev_fmt,'(a)') "(1x,'[',i2.2,1x,a11,']',2x)"
 !
 !--Define all the variables to be included in the .ev file and their supplementary information
 !
 i = 1
 j = 1
 call fill_ev_tag(ev_fmt,'time',     '0',  i,j)
 call fill_ev_tag(ev_fmt,'ekin',     '0',  i,j)
 call fill_ev_tag(ev_fmt,'etherm',   '0',  i,j)
 call fill_ev_tag(ev_fmt,'emag',     '0',  i,j)
 call fill_ev_tag(ev_fmt,'epot',     '0',  i,j)
 call fill_ev_tag(ev_fmt,'etot',     '0',  i,j)
 call fill_ev_tag(ev_fmt,'totmom',   '0',  i,j)
 call fill_ev_tag(ev_fmt,'angtot',   '0',  i,j)
 call fill_ev_tag(ev_fmt,'rho',      'xa', i,j)
 call fill_ev_tag(ev_fmt,'dt',       '0',  i,j)
 call fill_ev_tag(ev_fmt,'totentrop','s',  i,j)
 call fill_ev_tag(ev_fmt,'rmsmach',  '0',  i,j)
 if (.not. gas_only) then
    if (npartoftype(igas)        > 0) call fill_ev_tag(ev_fmt,'rho gas', 'xa', i,j)
    if (npartoftype(idust)       > 0) call fill_ev_tag(ev_fmt,'rho dust','xa', i,j)
    if (npartoftype(iboundary)   > 0) call fill_ev_tag(ev_fmt,'rho bdy', 'xa', i,j)
    if (npartoftype(istar)       > 0) call fill_ev_tag(ev_fmt,'rho star','xa', i,j)
    if (npartoftype(idarkmatter) > 0) call fill_ev_tag(ev_fmt,'rho dm',  'xa', i,j)
    if (npartoftype(ibulge)      > 0) call fill_ev_tag(ev_fmt,'rho blg', 'xa', i,j)
 endif
 if (maxalpha==maxp)                  call fill_ev_tag(ev_fmt,'alpha',   'x',  i,j)
 if ( mhd ) then
    call fill_ev_tag(ev_fmt,      'divB',     'xa', i,j)
    call fill_ev_tag(ev_fmt,      'hdivB/B',  'xa', i,j)
    call fill_ev_tag(ev_fmt,      'beta',     'xan',i,j)
    if (mhd_nonideal) then
       call fill_ev_tag(ev_fmt,   'temp',     'xan',i,j)
       call fill_ev_tag(ev_fmt,   'eta_ar',   'xan',i,j)
       if (use_ohm) then
          call fill_ev_tag(ev_fmt,'eta_o',    'xan',i,j)
          call fill_ev_tag(ev_fmt,'eta_o/art','xan',i,j)
       endif
       if (use_hall) then
          call fill_ev_tag(ev_fmt,'eta_h',    'xan',i,j)
          call fill_ev_tag(ev_fmt,'|eta_h|',  'xan',i,j)
          call fill_ev_tag(ev_fmt,'eta_h/art','xan',i,j)
          call fill_ev_tag(ev_fmt,'|e_h|/art','xan',i,j)
       endif
       if (use_ambi) then
          call fill_ev_tag(ev_fmt,'eta_a',    'xan',i,j)
          call fill_ev_tag(ev_fmt,'eta_a/art','xan',i,j)
          call fill_ev_tag(ev_fmt,'velocity', 'xan',i,j)
          call fill_ev_tag(ev_fmt,'v_ion',    'xan',i,j)
          call fill_ev_tag(ev_fmt,'v_drift',  'xan',i,j)
       endif
          call fill_ev_tag(ev_fmt,'ni/n(i+n)','xan',i,j)
!         call fill_ev_tag(ev_fmt,'ne/n(i+n)','xan',i,j)
          call fill_ev_tag(ev_fmt,'n_e',      'xa', i,j)
          call fill_ev_tag(ev_fmt,'n_n',      'xa', i,j)
       if (ion_rays) then
          call fill_ev_tag(ev_fmt,'n_ihR',    'xa', i,j)
          call fill_ev_tag(ev_fmt,'n_imR',    'xa', i,j)
          call fill_ev_tag(ev_fmt,'n_g(Z=-1)','xa', i,j)
          call fill_ev_tag(ev_fmt,'n_g(Z= 0)','xa', i,j)
          call fill_ev_tag(ev_fmt,'n_g(Z=+1)','xa', i,j)
       endif
       if (ion_thermal) then
          call fill_ev_tag(ev_fmt,'n_isT',    'xa', i,j)
          call fill_ev_tag(ev_fmt,'n_idT',    'xa', i,j)
       endif
    endif
 endif
 if (use_dustfrac) then
    call fill_ev_tag(ev_fmt,   'dust/gas',   'xan',i,j)
    call fill_ev_tag(ev_fmt,   't_s',        'mn', i,j)
 endif
 if (iexternalforce > 0) then
    call fill_ev_tag(ev_fmt,   'totmomall',  '0',  i,j)
    call fill_ev_tag(ev_fmt,   'angall',     '0',  i,j)
    if (iexternalforce==iext_binary) then
       call fill_ev_tag(ev_fmt,'Macc sink 1','0',  i,j)
       call fill_ev_tag(ev_fmt,'Macc sink 2','0',  i,j)
    endif
 endif
 if (iexternalforce>0 .or. nptmass > 0 .or. icreate_sinks > 0) then
    call fill_ev_tag(ev_fmt,'accretedmas','0',  i,j)
    call fill_ev_tag(ev_fmt,'eacc',       '0',  i,j)
    track_mass     = .true.
 else
    track_mass     = .false.
 endif
 if (ishock_heating==0 .or. ipdv_heating==0 .or. lightcurve) then
    call fill_ev_tag(ev_fmt,'tot lum', '0',  i,j)
    track_lum      = .true.
 else
    track_lum      = .false.
 endif
 if (calc_erot) then
    call fill_ev_tag(ev_fmt,'erot',    '0',  i,j)
    call fill_ev_tag(ev_fmt,'erot_x',  '0',  i,j)
    call fill_ev_tag(ev_fmt,'erot_y',  '0',  i,j)
    call fill_ev_tag(ev_fmt,'erot_z',  '0',  i,j)
 endif
 if (irealvisc /= 0) then
    call fill_ev_tag(ev_fmt,'visc_rat','xan',i,j)
 endif
 iquantities = i - 1 ! The number of different quantities to analyse
 ielements   = j - 1 ! The number of values to be calculated (i.e. the number of columns in .ve)
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
!  creates up to three lables per input value, and fills the required
!  tracking arrays; this includes a check to verify the actions are legal
!+
!----------------------------------------------------------------
subroutine fill_ev_tag(ev_fmt,label,cmd,i,j)
 integer,          intent(inout) :: i,j
 character(len=*), intent(in)    :: ev_fmt,label,cmd
 integer                         :: iindex

 ! set tag
 ev_tag(i) = label
 ! set command
 ev_action(i) = cmd
 ! set first element in the ev_data array
 ev_istart(i) = j
 !
 ! make the headers
 if (index(cmd,'0') > 0) call fill_ev_header(ev_fmt,label,'0',j,index(cmd,'0'))
 if (index(cmd,'s') > 0) call fill_ev_header(ev_fmt,label,'s',j,index(cmd,'s'))
 if (index(cmd,'x') > 0) call fill_ev_header(ev_fmt,label,'x',j,index(cmd,'x'))
 if (index(cmd,'a') > 0) call fill_ev_header(ev_fmt,label,'a',j,index(cmd,'a'))
 if (index(cmd,'n') > 0) call fill_ev_header(ev_fmt,label,'n',j,index(cmd,'n'))
 i = i + 1
 j = j + len(trim(cmd))
 !
 ! verify action command is legal
 if ( (index(cmd,'x') > 0) .or. (index(cmd,'a') > 0) .or. (index(cmd,'n') > 0) ) then
   iindex = 1
 else
   iindex = 0
 endif
 if ( index(cmd,'0') + index(cmd,'s') + iindex > 1) &
    call fatal('fill_ev_tag','using an invalid sequence of actions for element', var=cmd)
 !
end subroutine fill_ev_tag
!----------------------------------------------------------------
!+
!  Fill an array to be used for the header of the .ev file
!+
!----------------------------------------------------------------
subroutine fill_ev_header(ev_fmt,label,cxmn,j,joffset)
 integer,           intent(in) :: j,joffset
 character(len=* ), intent(in) :: ev_fmt,label
 character(len= 1), intent(in) :: cxmn
 character(len=11)             :: label0
 character(len= 3)             :: ext
 integer                       :: j_actual

 if (len(label)>11 .and. (cxmn=='0' .or. cxmn=='s') ) then
    label0 = label(1:11)
 else if (len(label)>9) then
    label0 = label(1:9)
 else
    label0 = label
 endif
 ext = ""
 if (len(label)<=7) then
    if (cxmn=='x') ext = "max"
    if (cxmn=='a') ext = "ave"
    if (cxmn=='n') ext = "min"
 else if (len(label)<=9) then
    if (cxmn=='x') ext = "X"
    if (cxmn=='a') ext = "A"
    if (cxmn=='n') ext = "N"
 endif
 if (ext/="") write(label0,'(a,1x,a)')trim(label0),trim(ext);
 !
 j_actual = j + joffset - 1
 if (j_actual > 99) then
    write(ev_label(j_actual),ev_fmt) 100-j_actual,trim(label0)
 else
    write(ev_label(j_actual),ev_fmt) j_actual,trim(label0)
 endif

end subroutine fill_ev_header
!----------------------------------------------------------------
!+
!  calculates total energy, etc, and writes line to .ev file
!+
!----------------------------------------------------------------
subroutine write_evfile(t,dt)
 use energies,      only:compute_energies,ev_data_update
 use io,            only:id,master
 use options,       only:iexternalforce
 use extern_binary, only:accretedmass1,accretedmass2
 real, intent(in)  :: t, dt
 character(len=35) :: ev_format

 call compute_energies(t)

 if (id==master) then
    !
    !--fill in additional details that are not calculated in energies.f
    call ev_data_update(ev_data,'dt',dt)
    if (iexternalforce==iext_binary) then
       call ev_data_update(ev_data,'Macc sink 1',accretedmass1)
       call ev_data_update(ev_data,'Macc sink 2',accretedmass2)
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
 use energies,  only:ekin,etherm,emag,epot,etot,rmsmach,vrms,accretedmass,mdust,mgas
 use energies,  only:ev_get_value
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
 if (track_mass) write(iprint,"(1x,('E',a,'=',es10.3))") 'acc',ev_get_value('eacc')
 write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3))") &
      'Linm',ev_get_value('totmom'),'Angm',ev_get_value('angtot')
 if (iexternalforce > 0) then
    if (abs(ev_get_value('angall')-ev_get_value('angtot')) > tiny(0.)) then
       write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3),a)") &
       'Linm',ev_get_value('totmomall'),'Angm',ev_get_value('angall'),' [including accreted particles]'
    endif
 endif

 write(iprint,"(1x,a,'(max)=',es10.3,' (mean)=',es10.3,' (max)=',es10.3,a)") &
      'density  ',ev_get_value('rho','x'),ev_get_value('rho','a'),ev_get_value('rho','x')*unit_density,' g/cm^3'

 if (use_dustfrac) then
    write(iprint,"(1x,a,'(max)=',es10.3,1x,'(mean)=',es10.3,1x,'(min)=',es10.3)") &
         'dust2gas ',ev_get_value('dust/gas','x'),ev_get_value('dust/gas','n')
    write(iprint,"(3x,a,'(mean)=',es10.3,1x,'(min)=',es10.3)") 't_stop ',ev_get_value('t_s','a'),ev_get_value('t_s','n')
    write(iprint,"(1x,'Mgas = ',es10.3,', Mdust = ',es10.3)") mgas,mdust
 endif

 write(iprint,"(1x,1(a,'=',es10.3))") &
      'Accreted mass',accretedmass

 string = ''
 if (maxalpha==maxp) then
    write(string,"(a,'(max)=',es10.3)") ' alpha',ev_get_value('alpha','x')
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
         'Ratio of physical-to-art. visc',ev_get_value('visc_rat','x'),ev_get_value('visc_rat','n')
 else
    write(iprint,"(1x,1(a,'=',es10.3))") &
        'RMS Mach #',rmsmach
 endif

 if (mhd) then
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),(a,'(mean)=',es10.3))") &
      'div B ',ev_get_value('divB','x'),'div B ',ev_get_value('divB','a')
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),(a,'(mean)=',es10.3))") &
      'h|div B|/B ',ev_get_value('hdivB/B','x'),'h|div B|/B ',ev_get_value('hdivB/B','a')
 endif
 write(iprint,"(/)")

 return
end subroutine write_evlog
!----------------------------------------------------------------
end module evwrite
