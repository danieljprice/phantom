!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module readwrite_dumps_common
!
! readwrite_dumps_common
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, boundary_dyn, checkconserved, dim, dump_utils,
!   dust, dust_formation, eos, externalforces, fileutils, gitinfo, io,
!   options, part, setup_params, sphNGutils, timestep, units
!
 use dump_utils, only:lenid
 implicit none

contains

!--------------------------------------------------------------------
!+
!  contruct header string based on compile-time options
!  these are for information only (ie. not important for restarting)
!+
!--------------------------------------------------------------------
character(len=lenid) function fileident(firstchar,codestring)
 use part,    only:mhd,npartoftype,idust,gravity,lightcurve
 use options, only:use_dustfrac
 use dim,     only:use_dustgrowth,phantom_version_string,use_krome,store_dust_temperature,do_nucleation,h2chemistry
 use gitinfo, only:gitsha
 character(len=2), intent(in) :: firstchar
 character(len=*), intent(in), optional :: codestring
 character(len=10) :: datestring, timestring
 character(len=30) :: string
!
!--print date and time stamp in file header
!
 call date_and_time(datestring,timestring)
 datestring = datestring(7:8)//'/'//datestring(5:6)//'/'//datestring(1:4)
 timestring = timestring(1:2)//':'//timestring(3:4)//':'//timestring(5:)

 string = ' '
 if (gravity) string = trim(string)//'+grav'
 if (npartoftype(idust) > 0) string = trim(string)//'+dust'
 if (use_dustfrac) string = trim(string)//'+1dust'
 if (h2chemistry) string = trim(string)//'+H2chem'
 if (lightcurve) string = trim(string)//'+lightcurve'
 if (use_dustgrowth) string = trim(string)//'+dustgrowth'
 if (use_krome) string = trim(string)//'+krome'
 if (store_dust_temperature) string = trim(string)//'+Tdust'
 if (do_nucleation) string = trim(string)//'+nucleation'
 if (present(codestring)) then
    fileident = firstchar//':'//trim(codestring)//':'//trim(phantom_version_string)//':'//gitsha
 else
    fileident = firstchar//':Phantom'//':'//trim(phantom_version_string)//':'//gitsha
 endif

 if (mhd) then
    fileident = trim(fileident)//' (mhd+clean'//trim(string)//')  : '//trim(datestring)//' '//trim(timestring)
 else
    fileident = trim(fileident)//' (hydro'//trim(string)//'): '//trim(datestring)//' '//trim(timestring)
 endif

end function fileident

!--------------------------------------------------------------------
!+
!  extract various options used in Phantom from the fileid string
!+
!--------------------------------------------------------------------
subroutine get_options_from_fileid(fileid,tagged,phantomdump,smalldump,&
                                   use_onefluiddust,ierr)
 character(len=lenid), intent(in)  :: fileid
 logical,              intent(out) :: tagged,phantomdump,smalldump,use_onefluiddust
 integer,              intent(out) :: ierr
!
!--if file is a small dump, return an error code but still read what
!  can be read from a small dump
!
 ierr = 0
 tagged      = .false.
 smalldump   = .false.
 phantomdump = .false.
 if (fileid(2:2)=='T') tagged = .true.
 if (fileid(1:1) /= 'F') then
    !write(*,*) 'ERROR! file header indicates file is not a full dump'
    ierr = 1
    if (fileid(1:1)=='S') smalldump = .true.
 endif
 if (index(fileid,'Phantom') /= 0) then
    phantomdump = .true.
 elseif (index(fileid,'sphNG') /= 0) then
    phantomdump = .false.
    write(*,*) 'reading dump in sphNG format'
 else
    write(*,*) 'WARNING: could not determine Phantom/sphNG from fileident'
    write(*,*) '(assuming sphNG...)'
    phantomdump = .false.
 endif
 if (index(fileid,'+1dust') /= 0) then
    use_onefluiddust = .true.
 else
    use_onefluiddust = .false.
 endif

end subroutine get_options_from_fileid

!--------------------------------------------------------------------
!+
!  utility to extract header variables to phantom
!+
!-------------------------------------------------------------------
subroutine unfill_header(hdr,phantomdump,got_tags,nparttot, &
                         nblocks,npart,npartoftype, &
                         tfile,hfactfile,alphafile,iprint,id,nprocs,ierr)
 use dim,        only:maxdustlarge,use_dust
 use io,         only:master,iverbose ! check this
 use eos,        only:isink
 use part,       only:maxtypes,igas,idust,ndustsmall,ndustlarge,ndusttypes,&
                      npartoftypetot
 use units,      only:udist,umass,utime,set_units_extra,set_units
 use dump_utils, only:extract,dump_h
 use fileutils,  only:make_tags_unique
 type(dump_h),    intent(in)  :: hdr
 logical,         intent(in)  :: phantomdump,got_tags
 integer(kind=8), intent(out) :: nparttot
 integer,         intent(out) :: nblocks,npart,npartoftype(maxtypes)
 real,            intent(out) :: tfile,hfactfile,alphafile
 integer,         intent(in)  :: iprint,id,nprocs
 integer,         intent(out) :: ierr
 integer         :: nparttoti,npartoftypetoti(maxtypes),ntypesinfile,nptinfile
 integer         :: ierr1,ierrs(3),i,counter
 integer(kind=8) :: ntypesinfile8
 character(len=10) :: dust_label(maxdustlarge)

 ierr = 0
 nparttot = 0
 npartoftypetot(:) = 0
 npart = 0
 npartoftype(:) = 0
 isink = 0
 call extract('ntypes',ntypesinfile,hdr,ierr1)
 if (ierr1 /= 0 .or. ntypesinfile < 1) then
    if (phantomdump .and. got_tags) then
       ierr = 4
       return
    else
       ntypesinfile = 5
    endif
 endif

 ! extract quantities from integer header
 call extract('nparttot',nparttoti,hdr,ierr1)
 if (ierr1 /= 0) then
    ierr = 5
    return
 endif
 if (ntypesinfile > maxtypes) then
    write(*,*) 'WARNING: number of particle types in file exceeds array limits'
    write(*,*) 'READING ONLY FIRST ',maxtypes,' OF ',ntypesinfile,' particle types'
    ntypesinfile = maxtypes
 endif
 call extract('npartoftype',npartoftypetoti(1:ntypesinfile),hdr,ierr1)
 if (ierr1 /= 0) then
    npartoftype(1) = nparttoti  ! assume only gas particles
 endif
 call extract('nblocks',nblocks,hdr,ierr1,default=1)
 if (ierr1 /= 0) write(*,*) 'number of MPI blocks not read: assuming 1'

 nparttot = int(nparttoti,kind=8)
 npartoftypetot = int(npartoftypetoti,kind=8)
 if (nblocks==1) then
    npartoftype(1:ntypesinfile) = int(npartoftypetot(1:ntypesinfile))
    if (npartoftype(idust) > 0) write(*,*) 'n(gas) = ',npartoftype(igas)
    counter = 0
    do i=1,maxdustlarge
       if (npartoftype(idust+i-1) > 0) then
          counter = counter + 1
       endif
    enddo
    dust_label = 'dust'
    call make_tags_unique(counter,dust_label)
    do i=1,counter
       write(*,*) 'n('//trim(dust_label(i))//') = ',npartoftype(idust+i-1)
    enddo
 endif
 call extract('isink',isink,hdr,ierr1)

!--non-MPI dumps
 if (nprocs==1) then
    if (nparttoti > huge(npart)) then
       write (*,*) 'ERROR in readdump: number of particles exceeds 32 bit limit, must use int(kind=8)''s ',nparttoti
       ierr = 4
       return
    endif
 endif
 if (nblocks==1) then
    npart = int(nparttoti)
    nparttot = npart
    if (id==master .and. iverbose >= 0) write (iprint,*) 'npart = ',npart
 endif
 if (got_tags) then
    call extract('ntypes',ntypesinfile8,hdr,ierr1)
    ntypesinfile = int(ntypesinfile8)
 endif
 if (ntypesinfile > maxtypes) then
    write(*,*) 'WARNING: number of particle types in file exceeds array limits'
    write(*,*) 'READING ONLY FIRST ',maxtypes,' OF ',ntypesinfile,' particle types'
    ntypesinfile = maxtypes
 endif
 call extract('nparttot',nparttot,hdr,ierr1)
 if (nblocks > 1) then
    call extract('npartoftype',npartoftype(1:ntypesinfile),hdr,ierr1)
 endif
 if (id==master .and. iverbose >= 0) write(*,*) 'npart(total) = ',nparttot
!
!--number of dust species
!
 if (use_dust) then
    call extract('ndustsmall',ndustsmall,hdr,ierrs(1))
    call extract('ndustlarge',ndustlarge,hdr,ierrs(2))
    if (any(ierrs(1:2) /= 0)) then
       call extract('ndustfluids',ndustsmall,hdr,ierrs(1)) ! for backwards compatibility
       if (ierrs(1) /= 0) write(*,*) 'ERROR reading number of small/large grain types from file header'
    endif
    ndusttypes = ndustsmall + ndustlarge
 endif
!
!--units
!
 call extract('udist',udist,hdr,ierrs(1))
 call extract('umass',umass,hdr,ierrs(2))
 call extract('utime',utime,hdr,ierrs(3))
 if (all(ierrs(1:3)==0)) then
    call set_units_extra()
 else
    write(iprint,*) 'ERROR reading units from dump file, assuming default'
    call set_units()  ! use default units
 endif
 ! get nptmass from header, needed to figure out if gwinspiral info is sensible
 call extract('nptmass',nptinfile,hdr,ierrs(1))
!--default real
 call unfill_rheader(hdr,phantomdump,ntypesinfile,nptinfile,&
                     tfile,hfactfile,alphafile,iprint,ierr)
 if (ierr /= 0) return

 if (id==master .and. iverbose >= 0) write(iprint,*) 'time = ',tfile

end subroutine unfill_header

!--------------------------------------------------------------------
!+
!  subroutine to fill the real header with various things
!+
!-------------------------------------------------------------------
subroutine fill_header(sphNGdump,t,nparttot,npartoftypetot,nblocks,nptmass,hdr,ierr)
 use eos,            only:write_headeropts_eos,polyk2
 use options,        only:tolh,alpha,alphau,alphaB,iexternalforce,ieos
 use part,           only:massoftype,hfact,Bextx,Bexty,Bextz,ndustsmall,ndustlarge,&
                          idust,grainsize,graindens,ndusttypes
 use checkconserved, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in,mtot_in
 use setup_params,   only:rhozero
 use timestep,       only:dtmax_user,idtmax_n_next,idtmax_frac_next,C_cour,C_force
 use externalforces, only:write_headeropts_extern
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax
 use boundary_dyn,   only:dynamic_bdy,dxyz,rho_bkg_ini,irho_bkg_ini
 use dump_utils,     only:reset_header,add_to_rheader,add_to_header,add_to_iheader,num_in_header,dump_h,maxphead
 use dim,            only:use_dust,maxtypes,use_dustgrowth,do_nucleation,use_apr,&
                          phantom_version_major,phantom_version_minor,phantom_version_micro,periodic,idumpfile
 use units,          only:udist,umass,utime,unit_Bfield
 use dust_formation, only:write_headeropts_dust_formation

 logical,         intent(in)    :: sphNGdump
 real,            intent(in)    :: t
 integer(kind=8), intent(in)    :: nparttot,npartoftypetot(:)
 integer,         intent(in)    :: nblocks,nptmass
 type(dump_h),    intent(inout) :: hdr
 integer,         intent(out)   :: ierr
 integer :: number

 ierr = 0
 ! default int
 call add_to_iheader(int(nparttot),'nparttot',hdr,ierr)
 call add_to_iheader(maxtypes,'ntypes',hdr,ierr)
 call add_to_iheader(int(npartoftypetot(1:maxtypes)),'npartoftype',hdr,ierr)
 call add_to_iheader(nblocks,'nblocks',hdr,ierr)
 call add_to_iheader(nptmass,'nptmass',hdr,ierr)
 call add_to_iheader(ndustlarge,'ndustlarge',hdr,ierr)
 call add_to_iheader(ndustsmall,'ndustsmall',hdr,ierr)
 call add_to_iheader(idust,'idust',hdr,ierr)
 call add_to_iheader(idtmax_n_next,'idtmax_n',hdr,ierr)
 call add_to_iheader(idtmax_frac_next,'idtmax_frac',hdr,ierr)
 call add_to_iheader(idumpfile,'idumpfile',hdr,ierr)
 call add_to_iheader(phantom_version_major,'majorv',hdr,ierr)
 call add_to_iheader(phantom_version_minor,'minorv',hdr,ierr)
 call add_to_iheader(phantom_version_micro,'microv',hdr,ierr)

 ! int*8
 call add_to_header(nparttot,'nparttot',hdr,ierr)
 call add_to_header(int(maxtypes,kind=8),'ntypes',hdr,ierr)
 call add_to_header(npartoftypetot(1:maxtypes),'npartoftype',hdr,ierr)

 ! int*4
 call add_to_header(iexternalforce,'iexternalforce',hdr,ierr)
 call add_to_header(ieos,'ieos',hdr,ierr)
 call write_headeropts_eos(ieos,hdr,ierr)

 ! default real variables
 call add_to_rheader(t,'time',hdr,ierr)
 call add_to_rheader(dtmax_user,'dtmax',hdr,ierr)
 call add_to_rheader(rhozero,'rhozero',hdr,ierr)
 if (sphNGdump) then ! number = 23
    call add_to_rheader(0.,'escaptot',hdr,ierr)
    call add_to_rheader(0.,'tkin',hdr,ierr)
    call add_to_rheader(0.,'tgrav',hdr,ierr)
    call add_to_rheader(0.,'tterm',hdr,ierr)
    call add_to_rheader(0.,'anglostx',hdr,ierr)
    call add_to_rheader(0.,'anglosty',hdr,ierr)
    call add_to_rheader(0.,'anglostz',hdr,ierr)
    call add_to_rheader(0.,'specang',hdr,ierr)
    call add_to_rheader(0.,'ptmassin',hdr,ierr)
    call add_to_rheader(0.,'tmag',hdr,ierr)
    call add_to_rheader(Bextx,'Bextx',hdr,ierr)
    call add_to_rheader(Bexty,'Bexty',hdr,ierr)
    call add_to_rheader(Bextz,'Bextz',hdr,ierr)
    call add_to_rheader(0.,'hzero',hdr,ierr)
    call add_to_rheader(1.5*polyk2,'uzero_n2',hdr,ierr)
    call add_to_rheader(0.,'hmass',hdr,ierr)
    call add_to_rheader(0.,'gapfac',hdr,ierr)
    call add_to_rheader(0.,'pmassinitial',hdr,ierr)
 else ! number = 49
    call add_to_rheader(hfact,'hfact',hdr,ierr)
    call add_to_rheader(tolh,'tolh',hdr,ierr)
    call add_to_rheader(C_cour,'C_cour',hdr,ierr)
    call add_to_rheader(C_force,'C_force',hdr,ierr)
    call add_to_rheader(alpha,'alpha',hdr,ierr)
    call add_to_rheader(alphau,'alphau',hdr,ierr)
    call add_to_rheader(alphaB,'alphaB',hdr,ierr)
    call add_to_rheader(massoftype,'massoftype',hdr,ierr) ! array
    if (do_nucleation) call write_headeropts_dust_formation(hdr,ierr)
    call add_to_rheader(Bextx,'Bextx',hdr,ierr)
    call add_to_rheader(Bexty,'Bexty',hdr,ierr)
    call add_to_rheader(Bextz,'Bextz',hdr,ierr)
    call add_to_rheader(0.,'dum',hdr,ierr)
    if (iexternalforce /= 0) call write_headeropts_extern(iexternalforce,hdr,t,ierr)
    if (periodic) then
       call add_to_rheader(xmin,'xmin',hdr,ierr)
       call add_to_rheader(xmax,'xmax',hdr,ierr)
       call add_to_rheader(ymin,'ymin',hdr,ierr)
       call add_to_rheader(ymax,'ymax',hdr,ierr)
       call add_to_rheader(zmin,'zmin',hdr,ierr)
       call add_to_rheader(zmax,'zmax',hdr,ierr)
    endif
    if (dynamic_bdy) then
       call add_to_rheader(dxyz,'dxyz',hdr,ierr)
       call add_to_iheader(irho_bkg_ini,'irho_bkg_ini',hdr,ierr)
       call add_to_rheader(rho_bkg_ini,'rho_bkg_ini',hdr,ierr)
    endif
    call add_to_rheader(get_conserv,'get_conserv',hdr,ierr)
    call add_to_rheader(mtot_in,'mtot_in',hdr,ierr)
    call add_to_rheader(etot_in,'etot_in',hdr,ierr)
    call add_to_rheader(angtot_in,'angtot_in',hdr,ierr)
    call add_to_rheader(totmom_in,'totmom_in',hdr,ierr)
    call add_to_rheader(mdust_in(1:ndusttypes),'mdust_in',hdr,ierr)
    if (use_dust) then
       call add_to_rheader(grainsize(1:ndusttypes),'grainsize',hdr,ierr)
       call add_to_rheader(graindens(1:ndusttypes),'graindens',hdr,ierr)
    endif
 endif

 ! real*8
 call add_to_header(udist,'udist',hdr,ierr)
 call add_to_header(umass,'umass',hdr,ierr)
 call add_to_header(utime,'utime',hdr,ierr)
 call add_to_header(unit_Bfield,'umagfd',hdr,ierr)

 if (ierr /= 0) write(*,*) ' ERROR: arrays too small writing rheader'

 number = num_in_header(hdr%realtags)
 if (number >= maxphead) then
    write(*,*) 'error: header arrays too small for number of items in header: will be truncated'
 endif

end subroutine fill_header

!--------------------------------------------------------------------
!+
!  subroutine to set runtime parameters having read the real header
!+
!-------------------------------------------------------------------
subroutine unfill_rheader(hdr,phantomdump,ntypesinfile,nptmass,&
                          tfile,hfactfile,alphafile,iprint,ierr)
 use io,             only:id,master
 use dim,            only:maxvxyzu,nElements,use_dust,use_dustgrowth,use_krome,do_nucleation,idumpfile
 use eos,            only:extract_eos_from_hdr, read_headeropts_eos
 use options,        only:ieos,iexternalforce
 use part,           only:massoftype,Bextx,Bexty,Bextz,mhd,periodic,&
                          maxtypes,grainsize,graindens,ndusttypes
 use checkconserved, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in,mtot_in
 use setup_params,   only:rhozero
 use externalforces, only:read_headeropts_extern,extract_iextern_from_hdr
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax,set_boundary
 use boundary_dyn,   only:dynamic_bdy,dxyz,irho_bkg_ini,rho_bkg_ini,rho_bkg_ini1
 use dump_utils,     only:extract,dump_h
 use dust,           only:grainsizecgs,graindenscgs
 use units,          only:unit_density,udist
 use timestep,       only:idtmax_n,idtmax_frac
 use dust_formation, only:read_headeropts_dust_formation
 type(dump_h), intent(in)  :: hdr
 logical,      intent(in)  :: phantomdump
 integer,      intent(in)  :: iprint,ntypesinfile,nptmass
 real,         intent(out) :: tfile,hfactfile,alphafile
 integer,      intent(out) :: ierr

 integer, parameter :: lu = 173
 integer            :: ierrs(10),iextern_in_file
 real               :: xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,dtmaxi
 real               :: alphaufile,alphaBfile,C_courfile,C_forcefile,tolhfile
 logical            :: iexist

 ierr  = 0
 call extract('time',tfile,hdr,ierr)
 if (ierr/=0)  call extract('gt',tfile,hdr,ierr)  ! this is sphNG's label for time
 call extract('dtmax',dtmaxi,hdr,ierr)
 call extract('rhozero',rhozero,hdr,ierr)
 Bextx = 0.
 Bexty = 0.
 Bextz = 0.
 if (phantomdump) then
    call extract('hfact',hfactfile,hdr,ierr)
    call extract('tolh',tolhfile,hdr,ierr)
    call extract('C_cour',C_courfile,hdr,ierr)
    call extract('C_force',C_forcefile,hdr,ierr)
    call extract('alpha',alphafile,hdr,ierr)
    if (maxvxyzu >= 4) then
       call extract('alphau',alphaufile,hdr,ierr)
    else
       alphaufile = 0.
    endif
    if (mhd) then
       call extract('alphaB',alphaBfile,hdr,ierr)
    endif

    if (extract_eos_from_hdr) call extract('ieos',ieos,hdr,ierr)

    call extract('massoftype',massoftype(1:ntypesinfile),hdr,ierr)
    if (ierr /= 0) then
       write(*,*) '*** ERROR reading massoftype from dump header ***'
       ierr = 4
    endif
    if (do_nucleation) then
       call read_headeropts_dust_formation(hdr,ierr)
       if (ierr /= 0) ierr = 6
    endif

    call extract('iexternalforce',iextern_in_file,hdr,ierrs(1))
    if (extract_iextern_from_hdr) iexternalforce = iextern_in_file
    if (iexternalforce /= 0) then
       call read_headeropts_extern(iexternalforce,hdr,nptmass,ierrs(1))
       if (ierrs(1) /= 0) ierr = 5
    elseif (iextern_in_file /= 0) then
       call read_headeropts_extern(iextern_in_file,hdr,nptmass,ierrs(1))
       if (ierrs(1) /= 0) ierr = 5
    endif

    call extract('idtmax_n',idtmax_n,hdr,ierr,default=1)
    call extract('idtmax_frac',idtmax_frac,hdr,ierr)
    call extract('idumpfile',idumpfile,hdr,ierr)
 else
    massoftype(1) = 0.
    hfactfile = 0.
 endif

 call read_headeropts_eos(ieos,hdr,ierr)

 if (periodic) then
    call extract('xmin',xmini,hdr,ierrs(1))
    call extract('xmax',xmaxi,hdr,ierrs(2))
    call extract('ymin',ymini,hdr,ierrs(3))
    call extract('ymax',ymaxi,hdr,ierrs(4))
    call extract('zmin',zmini,hdr,ierrs(5))
    call extract('zmax',zmaxi,hdr,ierrs(6))
    if (any(ierrs(1:6) /= 0)) then
       write(*,"(2(/,a))") ' ERROR: dump does not contain boundary positions', &
                           '        but we are using periodic boundaries'
       inquire(file='bound.tmp',exist=iexist)
       if (iexist) then
          open(unit=lu,file='bound.tmp')
          read(lu,*) xmini,xmaxi,ymini,ymaxi,zmini,zmaxi
          close(lu)
          call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
          write(*,"(a,6(es10.3,1x))") ' READ from bound.tmp ',xmin,xmax,ymin,ymax,zmin,zmax
       else
          write(*,"(3(/,a),/,/,a)") ' To silence this error and restart from an older dump file ', &
                           ' create an ascii file called "bound.tmp" in the current directory', &
                           ' with xmin,xmax,ymin,ymax,zmin & zmax in it, e.g.: ', &
                           ' 0. 1. 0. 1. 0. 1.'
          ierr = 5  ! spit fatal error
       endif
    else
       call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
    endif
 endif

 if (dynamic_bdy) then
    call extract('irho_bkg_ini',irho_bkg_ini,hdr,ierrs(1))
    call extract('rho_bkg_ini',rho_bkg_ini,hdr,ierrs(1))
    call extract('dxyz',dxyz,hdr,ierrs(2))
    if (rho_bkg_ini > 0.) then
       rho_bkg_ini1 = 1.0/rho_bkg_ini
    else
       rho_bkg_ini1 = 0.
    endif
 endif

 if (mhd) then
    call extract('Bextx',Bextx,hdr,ierrs(1))
    call extract('Bexty',Bexty,hdr,ierrs(2))
    call extract('Bextz',Bextz,hdr,ierrs(3))
    if (id==master) then
       if (any(ierrs(1:3) /= 0)) then
          write(*,*) 'ERROR reading external field (setting to zero)'
       else
          write(*,*) 'External field found, Bext = ',Bextx,Bexty,Bextz
       endif
    endif
 endif

 ! values to track that conserved values remain conserved
 call extract('get_conserv',get_conserv,hdr,ierrs(1))
 call extract('etot_in',    etot_in,    hdr,ierrs(2))
 call extract('angtot_in',  angtot_in,  hdr,ierrs(3))
 call extract('totmom_in',  totmom_in,  hdr,ierrs(4))
 call extract('mdust_in',   mdust_in(1:ndusttypes), hdr,ierrs(6))
 call extract('mtot_in',    mtot_in,    hdr,ierrs(5))
 if (any(ierrs(1:4) /= 0)) then
    write(*,*) 'ERROR reading values to verify conservation laws.  Resetting initial values.'
    get_conserv = 1.0
 endif


 !--pull grain size and density arrays if they are in the header
 !-- i.e. if dustgrowth is not ON
 if (use_dust .and. .not.use_dustgrowth) then
    call extract('grainsize',grainsize(1:ndusttypes),hdr,ierrs(1))
    call extract('graindens',graindens(1:ndusttypes),hdr,ierrs(2))
    if (any(ierrs(1:2) /= 0)) then
       write(*,*) 'ERROR reading grain size/density from file header'
       grainsize(1) = real(grainsizecgs/udist)
       graindens(1) = real(graindenscgs/unit_density)
    endif
 endif

end subroutine unfill_rheader

!---------------------------------------------------------------
!+
!  make sure required arrays have been read from Phantom file
!  and perform basic sanity checks
!+
!---------------------------------------------------------------
subroutine check_arrays(i1,i2,noffset,npartoftype,npartread,nptmass,nsinkproperties,massoftype,&
                        alphafile,tfile,phantomdump,got_iphase,got_xyzh,got_vxyzu,got_alpha, &
                        got_krome_mols,got_krome_gamma,got_krome_mu,got_krome_T, &
                        got_abund,got_dustfrac,got_sink_data,got_sink_vels,got_sink_sfprop,got_Bxyz,got_psi, &
                        got_dustprop,got_pxyzu,got_VrelVf,got_dustgasprop,got_rad,got_radprop,got_Tdust, &
                        got_eosvars,got_nucleation,got_iorig,got_apr_level,&
                        iphase,xyzh,vxyzu,pxyzu,alphaind,xyzmh_ptmass,Bevol,iorig,iprint,ierr)
 use dim,  only:maxp,maxvxyzu,maxalpha,maxBevol,mhd,h2chemistry,use_dustgrowth,gr,&
                do_radiation,store_dust_temperature,do_nucleation,use_krome,use_apr,store_sf_ptmass
 use eos,  only:ieos,polyk,gamma,eos_is_non_ideal
 use part, only:maxphase,isetphase,set_particle_type,igas,ihacc,ihsoft,imacc,ilum,ikappa,&
                xyzmh_ptmass_label,vxyz_ptmass_label,get_pmass,rhoh,dustfrac,ndusttypes,norig,&
                itemp,iX,iZ,imu,apr_level
 use io,   only:warning,id,master
 use options,        only:alpha,use_dustfrac,use_var_comp
 use sphNGutils,     only:itype_from_sphNG_iphase,isphNG_accreted
 use dust_formation, only:init_nucleation
 integer,         intent(in)    :: i1,i2,noffset,npartoftype(:),npartread,nptmass,nsinkproperties
 real,            intent(in)    :: massoftype(:),alphafile,tfile
 logical,         intent(in)    :: phantomdump,got_iphase,got_xyzh(:),got_vxyzu(:),got_alpha(:),got_dustprop(:)
 logical,         intent(in)    :: got_VrelVf,got_dustgasprop(:)
 logical,         intent(in)    :: got_abund(:),got_dustfrac(:),got_sink_data(:),got_sink_vels(:),got_sink_sfprop(:),got_Bxyz(:)
 logical,         intent(in)    :: got_krome_mols(:),got_krome_gamma,got_krome_mu,got_krome_T
 logical,         intent(in)    :: got_psi,got_Tdust,got_eosvars(:),got_nucleation(:),got_pxyzu(:),got_rad(:)
 logical,         intent(in)    :: got_radprop(:),got_iorig,got_apr_level
 integer(kind=1), intent(inout) :: iphase(:)
 integer(kind=8), intent(inout) :: iorig(:)
 real,            intent(inout) :: vxyzu(:,:),Bevol(:,:),pxyzu(:,:)
 real(kind=4),    intent(inout) :: alphaind(:,:)
 real,            intent(inout) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer,         intent(in)    :: iprint
 integer,         intent(out)   :: ierr
 logical :: use_gas
 integer :: i,itype,nread
 !
 ! particle type information
 !
 if (maxphase==maxp) then
    if (got_iphase) then
       if (phantomdump) then
          do i=i1,i2
             itype = int(iphase(i))
             iphase(i) = isetphase(itype,iactive=.true.)
          enddo
       else
          ! convert from sphNG
          do i=i1,i2
             itype = itype_from_sphNG_iphase(iphase(i))
             iphase(i) = isetphase(itype,iactive=.true.)
             if (itype==isphNG_accreted) then
                !  mark accreted/unknown particle types as dead according
                !  to Phantom (i.e., give negative smoothing length)
                xyzh(4,i) = -abs(xyzh(4,i))
                call set_particle_type(i,igas) ! to give an allowed particle type
             endif
          enddo
       endif
    elseif (any(npartoftype(2:) > 0)) then
       !
       !--iphase is required if there is more than one particle type
       !
       write(*,"(/,a,/)") 'error in rdump: need type information but iamtype not present in dump file'
       ierr = 8
       return
    else
       !
       !--iphase does not need to be read if there is only one particle type
       !  but to start the code it should be set such that all particles are active
       !
       do i=i1,i2
          iphase(i) = isetphase(igas,iactive=.true.)
       enddo
    endif
 endif
 if (maxphase==maxp) then
    use_gas = .false.
 else
    use_gas = .true.
 endif

 !
 ! hydrodynamics arrays
 !
 if (any(.not.got_xyzh)) then
    if (id==master .and. i1==1) write(*,*) 'ERROR: x, y, z or h not found in file'
    ierr = 9
    return
 endif
 if (any(.not.got_vxyzu(1:3))) then
    if (id==master .and. i1==1) write(*,*) 'ERROR: missing velocity information from file'
 endif
 if (maxvxyzu==4 .and. .not.got_vxyzu(4)) then
    if (gamma < 1.01) then
       do i=i1,i2
          vxyzu(4,i) = 1.5*polyk
          !print*,'u = ',vxyzu(4,i)
       enddo
       if (id==master .and. i1==1) write(*,"(/,a,/)") 'WARNING: u not in file but setting u = 3/2 * cs^2'
    else
       do i=i1,i2
          vxyzu(4,i) = (1.0/(gamma-1.0))*polyk*rhoh(xyzh(4,i),get_pmass(i,use_gas))**(gamma - 1.)
          !print*,'u = ',vxyzu(4,i)
       enddo
       if (id==master .and. i1==1) write(*,"(/,a,/)") 'WARNING: u not in file but setting u = (K*rho**(gamma-1))/(gamma-1)'
    endif
 endif
 if (h2chemistry .and. .not.all(got_abund).and. npartread > 0) then
    if (id==master) write(*,*) 'error in rdump: using H2 chemistry, but abundances not found in dump file'
    ierr = 9
    return
 endif
 if (use_krome) then
    if (.not.all(got_krome_mols).and. npartread > 0) then
       if (id==master) write(*,*) 'error in rdump: using KROME chemistry, but abundances not found in dump file'
       !     ierr = 9
       return
    endif
    if (.not.got_krome_gamma .and. npartread > 0) then
       if (id==master) write(*,*) 'error in rdump: using KROME chemistry, but gamma not found in dump file'
       !     ierr = 9
       return
    endif
    if (.not.got_krome_mu .and. npartread > 0) then
       if (id==master) write(*,*) 'error in rdump: using KROME chemistry, but mu not found in dump file'
       !     ierr = 9
       return
    endif
    if (.not.got_krome_T .and. npartread > 0) then
       if (id==master) write(*,*) 'error in rdump: using KROME chemistry, but temperature not found in dump file'
       !     ierr = 9
       return
    endif
 endif
 if (eos_is_non_ideal(ieos) .and. .not.got_eosvars(itemp)) then
    if (id==master .and. i1==1) write(*,"(/,a,/)") 'WARNING: missing temperature information from file'
 endif
 use_var_comp = (got_eosvars(iX) .and. got_eosvars(iZ) .and. got_eosvars(imu))
 if (store_dust_temperature .and. .not.got_Tdust) then
    if (id==master .and. i1==1) write(*,"(/,a,/)") 'WARNING: missing dust temperature information from file'
 endif
 if (maxalpha==maxp) then
    if (got_alpha(1)) then
       if (alphafile < 0.99 .and. tfile > 0.) then
          if (any(alphaind(1,i1:i2) > 1.0 .or. alphaind(1,i1:i2) < 0.)) then
             if (id==master) write(iprint,*) 'ERROR! AV alpha < 0 or alpha > 1 in dump file: using alpha'
             alphaind(1,i1:i2) = real(alpha,kind=4)
          endif
       endif
    else
       if (id==master .and. i1==1) write(*,"(/,a,/)") 'WARNING: alpha not found in file'
       alphaind(1,i1:i2) = real(alpha,kind=4)
    endif
 endif
 if (npartread > 0) then
    do i = 1, size(massoftype)
       if (npartoftype(i) > 0) then
          if (massoftype(i) <= 0.0) then
             if (id==master .and. i1==1) write(*,*) 'ERROR! mass not set in read_dump (Phantom)'
             ierr = 12
             return
          endif
       endif
    enddo
 endif
 if (use_dustfrac .and. .not. all(got_dustfrac(1:ndusttypes))) then
    if (id==master .and. i1==1) write(*,*) 'WARNING! using one-fluid dust, but no dust fraction found in dump file'
    if (id==master .and. i1==1) write(*,*) ' Setting dustfrac = 0'
    dustfrac = 0.
 endif
 if (use_dustgrowth .and. .not.got_dustprop(1)) then
    if (id==master) write(*,*) 'ERROR! using dustgrowth, but no grain mass found in dump file'
    ierr = ierr + 1
 endif
 if (use_dustgrowth .and. .not.got_dustprop(2)) then
    if (id==master) write(*,*) 'ERROR! using dustgrowth, but no grain density found in dump file'
    ierr = ierr + 1
 endif
 if (use_dustgrowth .and. .not.got_VrelVf) then
    if (id==master) write(*,*) 'ERROR! using dustgrowth, but no Vrel/Vfrag found in dump file'
    ierr = ierr + 1
 endif
 if (use_dustgrowth .and. .not.got_dustgasprop(3)) then
    if (id==master) write(*,*) 'ERROR! using dustgrowth, but no St found in dump file'
    ierr = ierr + 1
 endif
 !
 ! sink particle arrays
 !
 nread = 0
 if (nptmass > 0) then
    do i=1,nsinkproperties
       if (.not.got_sink_data(i)) then
          if (i <= 5) then
             if (id==master) write(*,*) 'ERROR! sink particle '//trim(xyzmh_ptmass_label(i))//' not found'
             ierr = 10
             return
          else
             if (id==master) write(*,*) 'WARNING! sink particle '//trim(xyzmh_ptmass_label(i))//' not found'
          endif
       endif
    enddo
    if (.not.all(got_sink_vels(1:3))) then
       if (id==master .and. i1==1) write(*,"(/,a,/)") 'WARNING! sink particle velocities not found'
    endif
    if ( store_sf_ptmass .and. .not.all(got_sink_sfprop(1:2))) then
       if (id==master .and. i1==1) write(*,"(/,a,/)") 'WARNING! sink particle link list not found'
    endif
    if (id==master .and. i1==1) then
       print "(2(a,i4),a)",' got ',nsinkproperties,' sink properties from ',nptmass,' sink particles'
       if (nptmass > 0) print "(1x,58('-'),/,1x,a,'|',5(a9,1x,'|'),/,1x,58('-'))",&
                              'ID',' Mass    ',' Racc    ',' Macc    ',' hsoft   ',' Lsink   '
       do i=1,min(nptmass,999)
          if (xyzmh_ptmass(4,i) >= 0.) print "(i3,'|',5(1pg9.2,1x,'|'))",i,xyzmh_ptmass(4,i),xyzmh_ptmass(ihacc,i),&
                                            xyzmh_ptmass(imacc,i),xyzmh_ptmass(ihsoft,i),xyzmh_ptmass(ilum,i)
       enddo
       if (nptmass > 0) print "(1x,58('-'))"
    endif
 endif
 !
 ! radiation arrays
 !
 if (do_radiation) then
    if (.not.all(got_rad)) then
       if (id==master .and. i1==1) write(*,*) 'ERROR: RADIATION=yes but radiation arrays not found in Phantom dump file'
       ierr = ierr + 1
    endif
    if (.not.got_radprop(ikappa)) then
       if (id==master .and. i1==1) write(*,"(/,1x,a,/)") 'WARNING: RADIATION=yes but opacity not found in Phantom dump file'
    endif
 endif

 !
 ! MHD arrays
 !
 if (mhd) then
    if (.not.all(got_Bxyz(1:3))) then
       if (id==master .and. i1==1) write(*,"(/,1x,a,/)") 'WARNING: MHD but magnetic field arrays not found in Phantom dump file'
    endif
    if (.not.got_psi) then
       if (id==master .and. i1==1) write(*,"(/,1x,a,/)") &
          'WARNING! div B cleaning field (Psi) not found in Phantom dump file: assuming psi=0'
       Bevol(maxBevol,i1:i2) = 0.
    endif
 endif

 !
 ! GR arrays
 !
 if (gr) then
    if (.not.all(got_pxyzu(1:3))) then
       write(*,"(/,1x,a,/)") 'WARNING: GR but momentum arrays not found in Phantom dump file'
       pxyzu(:,i1:i2) = 0.
    endif
 endif
 !
 ! Dust nucleation arrays
 !
 if (do_nucleation) then
    if (.not.all(got_nucleation)) then
       write(*,"(/,1x,a,/)") 'WARNING: DUST_NUCLEATION=yes but nucleation arrays not found in Phantom dump file'
       call init_nucleation()
    endif
 endif

!
! Particle IDs
!
 if (.not.got_iorig) then
    do i=i1,i2
       iorig(i) = i + noffset
    enddo
    norig = i2
    if (id==master .and. i1==1) write(*,"(/,1x,a,/)") 'WARNING: Particle IDs not in dump; resetting IDs'
 else
    norig = 0
    do i=i1,i2
       norig = max(norig,iorig(i))
    enddo
 endif

!
! APR
!
 if (use_apr .and. .not.got_apr_level) then
    do i = i1,i2
       apr_level(i) = 1
    enddo
    if (id==master .and. i1==1) write(*,"(/,1x,a,/)") 'WARNING: APR levels not in dump; setting to default'
 endif

end subroutine check_arrays

end module readwrite_dumps_common
