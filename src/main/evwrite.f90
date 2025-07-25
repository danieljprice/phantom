!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module evwrite
!
! Calculates conserved quantities etc and writes to .ev file;
!  Also writes log output
!  To Developer: To add values to the .ev file, follow the following procedure.
!     In the init_evfile subroutine in evwrite.F90, add the following command:
!        call fill_ev_label(ev_fmt,ev_tag_int,ev_tag_char,action,i,j)
!     and in compute_energies subroutine in energies.F90, add the following command:
!        call ev_data_update(ev_data_thread,ev_tag_int,value)
!     where
!        ev_fmt,ev_data_thread,i,j: pre-defined quantities to included verbatim
!        ev_tag_char: a string to identify the quantity for use in the header
!                     (e.g. 'c_s' for sound speed)
!        ev_tag_int:  an integer to identify the quantity (e.g. 'iev_cs' for sound speed);
!                     this integer must be included in energies (as a public variable,
!                     and in the openmp declarations), and passed to evwrite via use energies.
!        ev_value: the value of the quantity for particle i (e.g., spsoundi for sound speed)
!        action: a string identifying what action(s) you would like performed
!                on the quantity.  The available options are
!           0: no action taken (e.g. for time)
!           s: sum quantity (e.g. for entropy)
!           x: print the maximum quantity
!           a: print the average (mean) quantity
!           n: print the minimum quantity
!        where any or all of x,a,n can be used as a single action.  Although 0 & s are treated
!        the same, they are kept separate for clarity without added computational cost
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: boundary, boundary_dyn, dim, energies, eos,
!   externalforces, fileutils, gravwaveutils, io, mpiutils, nicil, options,
!   part, ptmass, timestep, units, viscosity
!
 use io,             only:fatal,iverbose
 use options,        only:iexternalforce
 use timestep,       only:dtmax_dratio
 use externalforces, only:iext_binary,was_accreted
 use energies,       only:inumev,iquantities,ev_data
 use energies,       only:ndead,npartall
 use energies,       only:gas_only,track_mass
 use energies,       only:iev_sum,iev_max,iev_min,iev_ave
 use energies,       only:iev_time,iev_ekin,iev_etherm,iev_emag,iev_epot,iev_etot,iev_totmom,iev_com,&
                          iev_angmom,iev_rho,iev_dt,iev_dtx,iev_entrop,iev_rmsmach,iev_vrms,iev_rhop,iev_alpha,&
                          iev_B,iev_divB,iev_hdivB,iev_beta,iev_temp,iev_etao,iev_etah,&
                          iev_etaa,iev_vel,iev_vhall,iev_vion,iev_n,&
                          iev_dtg,iev_ts,iev_dm,iev_momall,iev_angall,iev_angall,iev_maccsink,&
                          iev_macc,iev_eacc,iev_totlum,iev_erot,iev_viscrat,iev_erad,iev_gws,iev_mass,iev_bdy

 implicit none
 public                    :: init_evfile, write_evfile, write_evlog
 private                   :: fill_ev_tag, fill_ev_header

 integer,          private :: ielements
 integer,          private :: ev_cmd(inumev)    ! array of the actions to be taken
 character(len=19),private :: ev_label(inumev)  ! to make the header for the .ev file

 private

contains

!----------------------------------------------------------------
!+
!  opens the .ev file for output
!+
!----------------------------------------------------------------
subroutine init_evfile(iunit,evfile,open_file)
 use io,        only:id,master,warning
 use dim,       only:maxtypes,maxalpha,maxp,mhd,mhd_nonideal,track_lum
 use options,   only:calc_erot,use_dustfrac
 use units,     only:c_is_unity
 use part,      only:igas,idust,iboundary,istar,idarkmatter,ibulge,npartoftype,ndusttypes,maxtypes
 use nicil,     only:use_ohm,use_hall,use_ambi
 use viscosity, only:irealvisc
 use mpiutils,  only:reduceall_mpi
 use eos,       only:ieos,eos_is_non_ideal,eos_outputs_gasP
 use gravwaveutils, only:calc_gravitwaves
 use boundary_dyn,  only:dynamic_bdy
 integer,            intent(in) :: iunit
 character(len=  *), intent(in) :: evfile
 logical,            intent(in) :: open_file
 character(len= 27)             :: ev_fmt
 character(len= 11)             :: dustname
 integer                        :: i,j,k
 integer(kind=8)                :: npartoftypetot(maxtypes)
 !
 !--Initialise additional variables
 !
 npartoftypetot = reduceall_mpi('+', npartoftype)
 gas_only  = .true.
 do i = 2,maxtypes
    if (npartoftypetot(i) > 0) gas_only = .false.
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
 call fill_ev_tag(ev_fmt,iev_erad,   'erad',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_totmom, 'totmom',   '0', i,j)
 call fill_ev_tag(ev_fmt,iev_angmom, 'angtot',   '0', i,j)
 call fill_ev_tag(ev_fmt,iev_rho,    'rho',      'xa',i,j)
 call fill_ev_tag(ev_fmt,iev_dt,     'dt',       '0', i,j)
 if (dtmax_dratio > 0.) then
    call fill_ev_tag(ev_fmt,iev_dtx, 'dtmax',    '0', i,j)
 endif
 if (dynamic_bdy) then
    call fill_ev_tag(ev_fmt,iev_mass,'mass',     '0', i,j)
 endif
 call fill_ev_tag(ev_fmt,iev_entrop, 'totentrop','s', i,j)
 call fill_ev_tag(ev_fmt,iev_rmsmach,'rmsmach',  '0', i,j)
 call fill_ev_tag(ev_fmt,iev_vrms,   'vrms',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_com(1), 'xcom',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_com(2), 'ycom',     '0', i,j)
 call fill_ev_tag(ev_fmt,iev_com(3), 'zcom',     '0', i,j)
 if (.not. gas_only) then
    if (npartoftypetot(igas)        > 0) call fill_ev_tag(ev_fmt,iev_rhop(1),'rho gas', 'xa',i,j)
    if (npartoftypetot(idust)       > 0) call fill_ev_tag(ev_fmt,iev_rhop(2),'rho dust','xa',i,j)
    if (npartoftypetot(iboundary)   > 0) call fill_ev_tag(ev_fmt,iev_rhop(3),'rho bdy', 'xa',i,j)
    if (npartoftypetot(istar)       > 0) call fill_ev_tag(ev_fmt,iev_rhop(4),'rho star','xa',i,j)
    if (npartoftypetot(idarkmatter) > 0) call fill_ev_tag(ev_fmt,iev_rhop(5),'rho dm',  'xa',i,j)
    if (npartoftypetot(ibulge)      > 0) call fill_ev_tag(ev_fmt,iev_rhop(6),'rho blg', 'xa',i,j)
 endif
 if (maxalpha==maxp) then
    call fill_ev_tag(ev_fmt,      iev_alpha,  'alpha',  'x',  i,j)
 endif
 if (eos_is_non_ideal(ieos) .or. eos_outputs_gasP(ieos)) then
    call fill_ev_tag(ev_fmt,      iev_temp,   'temp',   'xan',i,j)
 endif
 if ( mhd ) then
    call fill_ev_tag(ev_fmt,      iev_B,      'B',      'xan',i,j)
    call fill_ev_tag(ev_fmt,      iev_divB,   'divB',   'xa' ,i,j)
    call fill_ev_tag(ev_fmt,      iev_hdivB,  'hdivB/B','xa' ,i,j)
    call fill_ev_tag(ev_fmt,      iev_beta,   'beta_P', 'xan',i,j)
    if (mhd_nonideal) then
       if (use_ohm) then
          call fill_ev_tag(ev_fmt,iev_etao,   'eta_o',    'xan',i,j)
       endif
       if (use_hall) then
          call fill_ev_tag(ev_fmt,iev_etah(1),'eta_h',    'xan',i,j)
          call fill_ev_tag(ev_fmt,iev_etah(2),'|eta_h|',  'xan',i,j)
          call fill_ev_tag(ev_fmt,iev_vhall,  'v_hall',   'xan',i,j)
       endif
       if (use_ambi) then
          call fill_ev_tag(ev_fmt,iev_etaa,   'eta_a',    'xan',i,j)
          call fill_ev_tag(ev_fmt,iev_vel,    'velocity', 'xan',i,j)
          call fill_ev_tag(ev_fmt,iev_vion,   'v_ion',    'xan',i,j)
       endif
       call fill_ev_tag(ev_fmt,   iev_n(1),   'ni/n(i+n)','xan',i,j)
       call fill_ev_tag(ev_fmt,   iev_n(2),   'ne/n(i+n)','xan',i,j)
       call fill_ev_tag(ev_fmt,   iev_n(3),   'n_e',      'xa', i,j)
       call fill_ev_tag(ev_fmt,   iev_n(4),   'n_n',      'xa', i,j)
       call fill_ev_tag(ev_fmt,   iev_n(5),   'n_g(Z=-1)','xa', i,j)
       call fill_ev_tag(ev_fmt,   iev_n(6),   'n_g(Z= 0)','xa', i,j)
       call fill_ev_tag(ev_fmt,   iev_n(7),   'n_g(Z=+1)','xa', i,j)
    endif
 endif
 if (use_dustfrac) then
    call fill_ev_tag(ev_fmt,   iev_dtg,'dust/gas',     'xan',i,j)
    call fill_ev_tag(ev_fmt,   iev_ts, 't_s',          'xn', i,j)
    do k=1,ndusttypes
       write(dustname,'(a,I3)') 'DustMass',k
       call fill_ev_tag(ev_fmt,iev_dm(k), dustname,    '0',  i,j)
    enddo
 endif
 if (iexternalforce > 0) then
    call fill_ev_tag(ev_fmt,   iev_momall,'totmomall',   '0',i,j)
    call fill_ev_tag(ev_fmt,   iev_angall,'angall',      '0',i,j)
    if (iexternalforce==iext_binary) then
       call fill_ev_tag(ev_fmt,iev_maccsink(1),'Macc sink 1', '0',i,j)
       call fill_ev_tag(ev_fmt,iev_maccsink(2),'Macc sink 2', '0',i,j)
    endif
 endif
 if (was_accreted(iexternalforce,-1.0)) then
    call fill_ev_tag(ev_fmt,iev_macc,     'accretedmas', 's',i,j)
    call fill_ev_tag(ev_fmt,iev_eacc,     'eacc',        '0',i,j)
    track_mass     = .true.
 else
    track_mass     = .false.
 endif
 if (track_lum) then
    call fill_ev_tag(ev_fmt,iev_totlum,'tot lum', '0',i,j)
 endif
 if (calc_erot) then
    call fill_ev_tag(ev_fmt,iev_erot(1),'erot_x',  's',i,j)
    call fill_ev_tag(ev_fmt,iev_erot(2),'erot_y',  's',i,j)
    call fill_ev_tag(ev_fmt,iev_erot(3),'erot_z',  's',i,j)
    call fill_ev_tag(ev_fmt,iev_erot(4),'erot',    '0',i,j)
 endif
 if (irealvisc /= 0) then
    call fill_ev_tag(ev_fmt,iev_viscrat,'visc_rat','xan',i,j)
 endif

 if (calc_gravitwaves) then
    call fill_ev_tag(ev_fmt,iev_gws(1),'hx_0','0',i,j)
    call fill_ev_tag(ev_fmt,iev_gws(2),'hp_0','0',i,j)
    call fill_ev_tag(ev_fmt,iev_gws(3),'hx_{30}','0',i,j)
    call fill_ev_tag(ev_fmt,iev_gws(4),'hp_{30}','0',i,j)
    call fill_ev_tag(ev_fmt,iev_gws(5),'hx_{60}','0',i,j)
    call fill_ev_tag(ev_fmt,iev_gws(6),'hp_{60}','0',i,j)
    call fill_ev_tag(ev_fmt,iev_gws(7),'hx_{90}','0',i,j)
    call fill_ev_tag(ev_fmt,iev_gws(8),'hp_{90}','0',i,j)
 endif
 if (dynamic_bdy) then
    call fill_ev_tag(ev_fmt,iev_bdy(1,1),'min_x','0',i,j)
    call fill_ev_tag(ev_fmt,iev_bdy(1,2),'max_x','0',i,j)
    call fill_ev_tag(ev_fmt,iev_bdy(2,1),'min_y','0',i,j)
    call fill_ev_tag(ev_fmt,iev_bdy(2,2),'max_y','0',i,j)
    call fill_ev_tag(ev_fmt,iev_bdy(3,1),'min_z','0',i,j)
    call fill_ev_tag(ev_fmt,iev_bdy(3,2),'max_z','0',i,j)
 endif
 iquantities = i - 1 ! The number of different quantities to analyse
 ielements   = j - 1 ! The number of values to be calculated (i.e. the number of columns in .ve)
 !
 !--all threads do above, but only master writes file
 !  (the open_file is to prevent an .ev file from being made during the test suite)
 !
 if (open_file .and. id == master) then
    !
    !--open the file for output
    !
    open(unit=iunit,file=evfile,form='formatted',status='replace')
    !
    !--write a header line
    !
    write(ev_fmt,'(a,I3,a)') '(',ielements+1,'a)'
    write(iunit,ev_fmt)'#',ev_label(1:ielements)
 endif

end subroutine init_evfile

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
 elseif (len(label)>9 .and. (cxmn=='x' .or. cxmn=='a' .or. cxmn=='n')) then
    label0 = label(1:9)
 else
    label0 = label
 endif
 ext = ""
 if (len(label)<=7) then
    if (cxmn=='x') ext = "max"
    if (cxmn=='a') ext = "ave"
    if (cxmn=='n') ext = "min"
 elseif (len(label)<=9) then
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
 use io,            only:id,master,ievfile
 use timestep,      only:dtmax_user
 use options,       only:iexternalforce
 use externalforces,only:accretedmass1,accretedmass2
 real, intent(in)  :: t,dt
 integer           :: i,j
 real              :: ev_data_out(ielements)
 character(len=35) :: ev_format

 call compute_energies(t)

 if (id==master) then
    !--fill in additional details that are not calculated in energies.f
    ev_data(iev_sum,iev_dt)  = dt
    ev_data(iev_sum,iev_dtx) = dtmax_user
    if (iexternalforce==iext_binary) then
       ev_data(iev_sum,iev_maccsink(1)) = accretedmass1
       ev_data(iev_sum,iev_maccsink(2)) = accretedmass2
    endif
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

end subroutine write_evfile
!----------------------------------------------------------------
!+
!  Writes nicely formatted output to the log file/screen
!  Must be called *after* a call to compute energies has been
!  performed
!+
!----------------------------------------------------------------
subroutine write_evlog(iprint)
 use dim,           only:maxp,maxalpha,mhd,maxvxyzu,periodic,mhd_nonideal,&
                         use_dust,maxdusttypes,do_radiation,inject_parts
 use energies,      only:ekin,etherm,emag,epot,etot,rmsmach,vrms,accretedmass,mdust,mgas,xyzcom
 use energies,      only:erad
 use part,          only:nptmass,ndusttypes
 use viscosity,     only:irealvisc,shearparam
 use boundary,      only:dxbound,dybound,dzbound
 use units,         only:unit_density
 use options,       only:use_dustfrac
 use fileutils,     only:make_tags_unique
 use ptmass,        only:icreate_sinks
 integer, intent(in) :: iprint
 character(len=120)  :: string,Mdust_label(maxdusttypes)
 integer             :: i

 if (ndead > 0 .or. nptmass > 0 .or. icreate_sinks > 0 .or. inject_parts .or. iverbose > 0) then
    write(iprint,"(1x,4(a,I10))") 'npart=',npartall,', n_alive=',npartall-ndead, &
                                  ', n_dead_or_accreted=',ndead,', nptmass=',nptmass
 endif

 write(iprint,"(1x,3('E',a,'=',es10.3,', '),('E',a,'=',es10.3))") 'tot',etot,'kin',ekin,'therm',etherm,'pot',epot

 if (mhd)          write(iprint,"(1x,('E',a,'=',es10.3))") 'mag',emag
 if (do_radiation) write(iprint,"(1x,('E',a,'=',es10.3))") 'rad',erad
 if (track_mass)   write(iprint,"(1x,('E',a,'=',es10.3))") 'acc',ev_data(iev_sum,iev_eacc)
 write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3))") &
      'Linm',ev_data(iev_sum,iev_totmom),'Angm',ev_data(iev_sum,iev_angmom)
 if (iexternalforce > 0) then
    if (abs(ev_data(iev_sum,iev_angall)-ev_data(iev_sum,iev_angmom)) > tiny(0.)) then
       write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3),a)") &
       'Linm',ev_data(iev_sum,iev_momall),'Angm',ev_data(iev_sum,iev_angall),' [including accreted particles]'
    endif
 endif
 write(iprint,"(1x,3(a,es10.3))") "Centre of Mass = ",xyzcom(1),", ",xyzcom(2),", ",xyzcom(3)

 if (ev_data(iev_max,iev_rho) > 0.) then ! avoid floating point exception if no gas particles
    write(iprint,"(1x,a,'(max)=',es10.3,' (mean)=',es10.3,' (max)=',es10.3,a)") &
      'density  ',ev_data(iev_max,iev_rho),ev_data(iev_ave,iev_rho),ev_data(iev_max,iev_rho)*unit_density,' g/cm^3'
 endif

 if (use_dustfrac) then
    write(iprint,"(1x,a,'(max)=',es10.3,1x,'(mean)=',es10.3,1x,'(min)=',es10.3)") &
         'dust2gas ',ev_data(iev_max,iev_dtg),ev_data(iev_ave,iev_dtg),ev_data(iev_min,iev_dtg)
    write(iprint,"(3x,a,'(mean)=',es10.3,1x,'(min)=',es10.3)") 't_stop ',ev_data(iev_ave,iev_ts),ev_data(iev_min,iev_ts)
 endif
 if (use_dust) then
    write(iprint,"(1x,'Mgas = ',es10.3)") mgas
    Mdust_label = 'Mdust'
    call make_tags_unique(ndusttypes,Mdust_label)
    do i=1,ndusttypes
       write(iprint,"(1x,1(a,' = ',es10.3))") trim(Mdust_label(i)),mdust(i)
    enddo
 endif

 if (track_mass) write(iprint,"(1x,1(a,'=',es10.3))") 'Accreted mass',accretedmass

 string = ''
 if (maxalpha==maxp) then
    if (ev_data(iev_max,iev_alpha) > 0.) write(string,"(a,'(max)=',es10.3)") ' alpha',ev_data(iev_max,iev_alpha)
 endif
 if (len_trim(string) > 0) write(iprint,"(a)") trim(string)

 if (irealvisc /= 0 .and. npartall > 0) then
    if (periodic) then
       if (irealvisc==1) then
          write(iprint,"(1x,1(a,'=',es10.3,', '),(a,'=',es10.3))") &
             'RMS Mach #',rmsmach,'Reynolds # ',vrms*min(dxbound,dybound,dzbound)/shearparam
       endif
    endif
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),('(mean)=',es10.3),(' (min)=',es10.3))") &
         'Ratio of physical-to-art. visc',ev_data(iev_max,iev_viscrat),ev_data(iev_min,iev_viscrat)
 elseif (npartall > 0) then
    write(iprint,"(1x,1(a,'=',es10.3))") 'RMS Mach #',rmsmach
 endif

 if (mhd) then
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),(a,'(mean)=',es10.3))") &
      'div B ',ev_data(iev_max,iev_divB),'div B ',ev_data(iev_ave,iev_divB)
    write(iprint,"(1x,1(a,'(max)=',es10.3,', '),(a,'(mean)=',es10.3))") &
      'h|div B|/B ',ev_data(iev_max,iev_hdivB),'h|div B|/B ',ev_data(iev_ave,iev_hdivB)
    if (ev_data(iev_max,iev_hdivB) > 10.) &
      write(iprint,'(a)') 'WARNING! h|div B|/B is growing!  Recommend increasing hdivbbmax_max for better stability'
 endif
 write(iprint,"(/)")

end subroutine write_evlog

end module evwrite
