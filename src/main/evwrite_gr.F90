module evwrite
 use io,             only: fatal
 use part,           only: npart
 use options,        only: iexternalforce
 use externalforces, only: iext_binary,was_accreted
 use energies,       only: inumev,iquantities,ev_data
 use energies,       only: ndead
 use energies,       only: erot_com,gas_only,track_mass,track_lum
 use energies,       only: iev_sum,iev_max,iev_min,iev_ave
 use energies,       only: iev_time,iev_ekin,iev_etherm,iev_emag,iev_epot,iev_etot,iev_totmom,&
                           iev_angmom,iev_rho,iev_dt,iev_entrop,iev_rmsmach,iev_vrms,iev_rhop,iev_alpha,&
                           iev_divB,iev_hdivB,iev_beta,iev_temp,iev_etaar,iev_etao,iev_etah,&
                           iev_etaa,iev_vel,iev_vion,iev_vdrift,iev_n,iev_nR,iev_nT,&
                           iev_dtg,iev_ts,iev_momall,iev_angall,iev_angall,iev_maccsink,&
                           iev_macc,iev_eacc,iev_totlum,iev_erot,iev_viscrat

 implicit none
 public                    :: init_evfile, write_evfile, write_evlog
 private                   :: fill_ev_tag, fill_ev_header

 integer,          private :: ievfile,ielements
 integer,          private :: ev_cmd(inumev)    ! array of the actions to be taken
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
 use dim,       only: maxtypes,maxalpha,maxp,mhd,mhd_nonideal,calc_erot,lightcurve
 use options,   only: ishock_heating,ipdv_heating,use_dustfrac
 use part,      only: igas,idust,iboundary,istar,idarkmatter,ibulge,npartoftype
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
 call fill_ev_tag(ev_fmt,iev_time,   'time',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_ekin,   'ekin',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_etherm, 'etherm',   '0', i,j)
 call fill_ev_tag(ev_fmt,iev_emag,   'emag',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_epot,   'epot',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_etot,   'etot',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_totmom, 'totmom',   '0', i,j)
 call fill_ev_tag(ev_fmt,iev_angmom, 'angtot',   '0', i,j)
 call fill_ev_tag(ev_fmt,iev_rho,    'rho',      'xa',i,j)
 call fill_ev_tag(ev_fmt,iev_dt,     'dt',       '0', i,j)
 call fill_ev_tag(ev_fmt,iev_entrop, 'totentrop','s', i,j)
 call fill_ev_tag(ev_fmt,iev_rmsmach,'rmsmach',  's', i,j)
 call fill_ev_tag(ev_fmt,iev_vrms,   'vrms',     's', i,j)
 if (.not. gas_only) then
    if (npartoftype(igas)        > 0) call fill_ev_tag(ev_fmt,iev_rhop(1),'rho gas', 'xa',i,j)
    if (npartoftype(idust)       > 0) call fill_ev_tag(ev_fmt,iev_rhop(2),'rho dust','xa',i,j)
    if (npartoftype(iboundary)   > 0) call fill_ev_tag(ev_fmt,iev_rhop(3),'rho bdy', 'xa',i,j)
    if (npartoftype(istar)       > 0) call fill_ev_tag(ev_fmt,iev_rhop(4),'rho star','xa',i,j)
    if (npartoftype(idarkmatter) > 0) call fill_ev_tag(ev_fmt,iev_rhop(5),'rho dm',  'xa',i,j)
    if (npartoftype(ibulge)      > 0) call fill_ev_tag(ev_fmt,iev_rhop(6),'rho blg', 'xa',i,j)
 endif
 if (maxalpha==maxp)                  call fill_ev_tag(ev_fmt,iev_alpha,  'alpha',   'x' ,i,j)

    track_mass     = .false.
    track_lum      = .false.

 if (calc_erot) then
    call fill_ev_tag(ev_fmt,iev_erot(1),'erot_x',  's',i,j)
    call fill_ev_tag(ev_fmt,iev_erot(2),'erot_y',  's',i,j)
    call fill_ev_tag(ev_fmt,iev_erot(3),'erot_z',  's',i,j)
    call fill_ev_tag(ev_fmt,iev_erot(4),'erot',    '0',i,j)
 endif
 if (irealvisc /= 0) then
    call fill_ev_tag(ev_fmt,iev_viscrat,'visc_rat','xan',i,j)
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
subroutine fill_ev_tag(ev_fmt,itag,label,cmd,i,j)
 integer,          intent(inout) :: i,j
 integer,          intent(out)   :: itag
 character(len=*), intent(in)    :: ev_fmt,label,cmd
 integer                         :: ki,kj,iindex,joffset

 ! initialise command
 itag      = i
 joffset   = 1
 ev_cmd(i) = 0
 !
 ! make the headers & set ev_cmd
 if (index(cmd,'0') > 0) call fill_ev_header(ev_fmt,label,'0',j,joffset)
 if (index(cmd,'s') > 0) call fill_ev_header(ev_fmt,label,'s',j,joffset)
 if (index(cmd,'x') > 0) then
    call fill_ev_header(ev_fmt,label,'x',j,joffset)
    ev_cmd(i) = ev_cmd(i) + 1
    joffset   = joffset   + 1
 endif
 if (index(cmd,'a') > 0) then
    call fill_ev_header(ev_fmt,label,'a',j,joffset)
    ev_cmd(i) = ev_cmd(i) + 2
    joffset   = joffset   + 1
 endif
 if (index(cmd,'n') > 0) then
    call fill_ev_header(ev_fmt,label,'n',j,joffset)
    ev_cmd(i) = ev_cmd(i) + 5
 endif
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
 do ki = 1,len(cmd)-1
    do kj = ki+1,len(cmd)
       if ( cmd(ki:ki)==cmd(kj:kj) ) then
          call fatal('fill_ev_tag','using duplicate actions for the same quantity', var=cmd)
       endif
    enddo
 enddo
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
 real, intent(in)  :: t,dt
 integer           :: i,j
 real              :: ev_data_out(ielements)
 character(len=35) :: ev_format

 call compute_energies(t)

 if (id==master) then
    !--fill in additional details that are not calculated in energies.f
    ev_data(iev_sum,iev_dt) = dt

    ! Fill in the data_out array
    j = 1
    do i = 1,iquantities
       if (ev_cmd(i)==0) then
          ! include the total value
          ev_data_out(j) = ev_data(iev_sum,i)
          j = j + 1
       else
          if (ev_cmd(i)==1 .or. ev_cmd(i)==3 .or. ev_cmd(i)==6 .or. ev_cmd(i)==8) then
             ! include the maximum value
             ev_data_out(j) = ev_data(iev_max,i)
             j = j + 1
          endif
          if (ev_cmd(i)==2 .or. ev_cmd(i)==3 .or. ev_cmd(i)==7 .or. ev_cmd(i)==8) then
             ! include the average value
             ev_data_out(j) = ev_data(iev_ave,i)
             j = j + 1
          endif
          if (ev_cmd(i)==5 .or. ev_cmd(i)==6 .or. ev_cmd(i)==7 .or. ev_cmd(i)==8) then
             ! include the minimum value
             ev_data_out(j) = ev_data(iev_min,i)
             j = j + 1
          endif
       endif
    enddo
    !
    !--write line to .ev file (should correspond to header, below)
    !
    write(ev_format,'(a,I3,a)')"(",ielements,"(1pe18.10,1x))"
    write(ievfile,ev_format) ev_data_out
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
 use dim,       only:maxp,maxalpha,mhd,maxvxyzu,periodic,mhd_nonideal,use_dust
 use energies,  only:ekin,etherm,emag,epot,etot,rmsmach,vrms,accretedmass,mdust,mgas
 use viscosity, only:irealvisc,shearparam
 use boundary,  only:dxbound,dybound,dzbound
 use units,     only:unit_density
 use options,   only:use_dustfrac
 integer, intent(in) :: iprint
 character(len=120)  :: string

 if (ndead > 0) then
    write(iprint,"(1x,a,I10,a,I10)") 'n_alive=',npart-ndead,', n_dead_or_accreted=',ndead
 endif
 write(iprint,"(1x,3('E',a,'=',es10.3,', '),('E',a,'=',es10.3))") &
      'tot',etot,'kin',ekin,'therm',etherm,'pot',epot
 if (mhd)        write(iprint,"(1x,('E',a,'=',es10.3))") 'mag',emag
 if (track_mass) write(iprint,"(1x,('E',a,'=',es10.3))") 'acc',ev_data(iev_sum,iev_eacc)
 write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3))") &
      'Linm',ev_data(iev_sum,iev_totmom),'Angm',ev_data(iev_sum,iev_angmom)
 if (iexternalforce > 0) then
    if (abs(ev_data(iev_sum,iev_angall)-ev_data(iev_sum,iev_angmom)) > tiny(0.)) then
       write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3),a)") &
       'Linm',ev_data(iev_sum,iev_momall),'Angm',ev_data(iev_sum,iev_angall),' [including accreted particles]'
    endif
 endif

 write(iprint,"(1x,a,'(max)=',es10.3,' (mean)=',es10.3,' (max)=',es10.3,a)") &
      'density  ',ev_data(iev_max,iev_rho),ev_data(iev_ave,iev_rho),ev_data(iev_max,iev_rho)*unit_density,' g/cm^3'

 if (use_dustfrac) then
    write(iprint,"(1x,a,'(max)=',es10.3,1x,'(mean)=',es10.3,1x,'(min)=',es10.3)") &
         'dust2gas ',ev_data(iev_max,iev_dtg),ev_data(iev_ave,iev_dtg)
    write(iprint,"(3x,a,'(mean)=',es10.3,1x,'(min)=',es10.3)") 't_stop ',ev_data(iev_ave,iev_ts),ev_data(iev_min,iev_ts)
 endif
 if (use_dust) write(iprint,"(1x,'Mgas = ',es10.3,', Mdust = ',es10.3)") mgas,mdust

 if (track_mass) write(iprint,"(1x,1(a,'=',es10.3))") 'Accreted mass',accretedmass

 string = ''
 if (maxalpha==maxp) then
    write(string,"(a,'(max)=',es10.3)") ' alpha',ev_data(iev_max,iev_alpha)
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
         'Ratio of physical-to-art. visc',ev_data(iev_max,iev_viscrat),ev_data(iev_min,iev_viscrat)
 else
    write(iprint,"(1x,1(a,'=',es10.3))") &
        'RMS Mach #',rmsmach
 endif

 if (mhd) then
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),(a,'(mean)=',es10.3))") &
      'div B ',ev_data(iev_max,iev_divB),'div B ',ev_data(iev_ave,iev_divB)
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),(a,'(mean)=',es10.3))") &
      'h|div B|/B ',ev_data(iev_max,iev_hdivB),'h|div B|/B ',ev_data(iev_ave,iev_hdivB)
 endif
 write(iprint,"(/)")

 return
end subroutine write_evlog
!----------------------------------------------------------------
end module evwrite
