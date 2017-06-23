!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   This module sets up a disc with embedded planets
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    H_R               -- H/R at R=Rin
!    R_in              -- inner radius
!    R_indust          -- inner radius in dust
!    R_out             -- outer radius in gas
!    R_outdust         -- outer radius in dust
!    accr1             -- star accretion radius
!    alphaSS           -- desired alphaSS
!    disc_m            -- disc mass in gas
!    dust_to_gas_ratio -- dust to gas ratio
!    grainsizeinp      -- grain size (in cm)
!    np                -- number of gas particles
!    np_dust           -- number of dust particles
!    nplanets          -- number of planets
!    pindex            -- p index
!    qindex            -- q index
!    sigma_naught      -- surface density at R=1 in gas
!    sigma_naughtdust  -- surface density at R=1 in dust
!    star_m            -- star mass
!    udist             -- distance unit in cm
!    umass             -- mass unit in g
!    xinc              -- inclination angle
!
!  DEPENDENCIES: centreofmass, dim, dust, infile_utils, io, options, part,
!    physcon, prompting, setdisc, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 use dim,            only:maxp,use_dust,use_dustfrac
 implicit none
 public :: setpart
 real :: R_in,R_out,xinc,star_m,disc_m,pindex,qindex,H_R,sigma_naught,grainsizeinp
 real :: disc_mdust,sigma_naughtdust,dust_to_gas_ratio,R_outdust,R_indust
 real :: accr1,alphaSS
 integer, parameter :: maxplanets = 9
 real    :: mplanet(maxplanets),rplanet(maxplanets),accrplanet(maxplanets)
 integer :: i,nplanets,np,np_dust
 integer :: itype,npart_inter,mass_set
 logical :: ismooth_edge
 character(len=*), dimension(maxplanets), parameter :: planets = &
  (/'1','2','3','4','5','6','7','8','9' /)

 private
 real(kind=8) :: u_dist,u_mass

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,        only:set_disc
 use part,           only:nptmass,xyzmh_ptmass,maxvxyzu,vxyz_ptmass, &
                          ihacc,ihsoft,igas,idust,dustfrac
 use io,             only:master,fatal
 use options,        only:iexternalforce,ieos,alpha
 use units,          only:set_units,umass,udist,utime
 use physcon,        only:au,solarm,jupiterm,pi,years
 use prompting,      only:prompt
 use timestep,       only:tmax,dtmax
 use centreofmass,   only:reset_centreofmass
 use dust,           only:grainsizecgs,set_dustfrac
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=100)               :: filename
 logical :: iexist
 real    :: phi,vphi,sinphi,cosphi,omega,r2,disc_m_within_r,period_longest
 integer :: ierr,j

 !
 !--Set problem parameters
 !

 u_dist = 1.*au
 u_mass = solarm
 call set_units(dist=u_dist,mass=u_mass,G=1.d0)

 np                 = 500000
 if (use_dust .and. .not. use_dustfrac) then
    np_dust = 100000
 else
    np_dust = 0
 endif
 npartoftype(:)     = 0
 npartoftype(igas)  = np
 if (use_dust) npartoftype(idust) = np_dust
 star_m             = 0.55*solarm/umass
 accr1              = 1.*au/udist
 R_in               = accr1
 R_out              = 150.*au/udist
 R_indust           = accr1
 R_outdust          = 150.*au/udist
 mass_set           = 0
 disc_m             = 0.13*solarm/umass
 if (use_dust) then
    dust_to_gas_ratio = 0.01
 else
    dust_to_gas_ratio = 0.
 endif
 sigma_naught       = 1.009E-07
 sigma_naughtdust   = 3.437E-08
 alphaSS            = 0.1
 xinc               = 0.0
 H_R                = 0.1
 pindex             = 0.14
 qindex             = 0.5*0.43

 if (size(vxyzu(:,1)) >= 4) then
    gamma = 5./3.
 else
    gamma = 1.0
 endif
 time   = 0.
 hfact = 1.2

 grainsizeinp = 0.1

 filename=trim(fileprefix)//'.setup'

 print "(/,65('-'),2(/,a),/,65('-'),/)",&
   ' Welcome to the HLTauri Setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_hltauinputfile(filename,ierr)
    call set_units(dist=u_dist,mass=u_mass,G=1.d0)
    if (ierr /= 0) then
       if (id==master) call write_hltauinputfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    !
    !--resolution
    !
    call prompt('Enter total number of gas particles',np,0,maxp)
    if (use_dust .and. .not.use_dustfrac) then
       call prompt('Enter total number of dust particles',np_dust,0,maxp-np)
    else
       np_dust = 0
    endif
    !
    !--units
    !
    call prompt('Enter code units of mass in g (default: 1 solar mass) ',u_mass,0.)
    call prompt('Enter code units of distance in cm (default: 1 AU) ',u_dist,0.)
    call set_units(dist=u_dist,mass=u_mass,G=1.d0)
    !
    !--set star options
    !
    call prompt('Enter mass for central star (code units)',star_m,0.)
    call prompt('Enter accretion radius for central star (code units)',accr1,0.)
    !
    !--set disc options
    !
    call prompt('Enter inner radius (code units) of protoplanetary disc in gas',R_in,0.01)
    call prompt('Enter outer radius (code units) of protoplanetary disc in gas',R_out,R_in)
    if (use_dust .and. .not.use_dustfrac) then
       call prompt('Enter inner radius (code units) of protoplanetary disc in dust',R_indust,0.01)
       call prompt('Enter outer radius (code units) of protoplanetary disc in dust',R_outdust,R_indust)
    endif
    call prompt('How do you want to set disc mass? (0=disc mass, 1=surface density, 2=Toomre Q (not yet implemented))',mass_set,0,2)
    select case(mass_set)
    case(0)
       call prompt('Enter mass (code units) of protoplanetary disc in gas',disc_m,0.0)
       if (use_dust) then
          call prompt('Enter dust to gas ratio',dust_to_gas_ratio,0.)
       endif
    case(1)
       call prompt('Enter suface density (code units) at R=1 (code units) in gas',sigma_naught,0.0)
       if (use_dust .and. .not.use_dustfrac) then
          call prompt('Enter suface density (code units) at R=1 (code units) in dust',sigma_naughtdust,0.0)
       endif
    case(2)
       print "(a)",'Toomre Q not yet implemented'
       stop
    end select
    if (use_dust) call prompt('Enter grain size in cm',grainsizeinp,0.)
    call prompt('Enter desired alpha_SS',alphaSS,0.)
    call prompt('Enter inclination angle to the binary in degrees',xinc)
    call prompt('Do you want Sigma to drop to zero at the inner edge?',ismooth_edge)
    call prompt('Enter H/R at R=Rin',H_R,0.)
    call prompt('Enter p index of surface density profile Sigma = Sigma0*R^-p',pindex,0.)
    call prompt('Enter q index of temperature profile cs = cs0*R^-q',qindex,0.)
    !
    !--add the planets
    !
    nplanets = 3
    call prompt('Enter the number of planets',nplanets,0,maxplanets)
    mplanet(:) = 0.001*star_m*umass/jupiterm
    accrplanet(:) = 0.25*au/udist
    do i=1,nplanets
       rplanet(i) = 10.*i*au/udist
       call prompt('Enter mass (in Jupiter mass) of planet '//trim(planets(i))//' :',mplanet(i),0.,star_m*umass/jupiterm)
       call prompt('Enter distance from the central star (code units) of planet '//trim(planets(i))//' :',rplanet(i),0.,R_out)
       call prompt('Enter accretion radius (code units) of planet '//trim(planets(i))//' :',accrplanet(i),0.)
    enddo

    !
    !--write default input file
    !
    call write_hltauinputfile(filename)

    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'

    stop
 else
    stop
 endif

 npartoftype(igas)  = np
 if (use_dust) npartoftype(idust) = np_dust
 npart              = sum(npartoftype)

 alpha = alphaSS
 xinc  = pi*xinc/180. ! convert inclination to radians

 npart_inter = 1
 itype       = igas
 if (use_dust .and. use_dustfrac) then
    disc_m = (1 + dust_to_gas_ratio)*disc_m
 endif
 select case(mass_set)
 case(0)
    call set_disc(id,master   = master,             &
                npart         = npartoftype(itype), &
                npart_start   = npart_inter,        &
                particle_type = itype,              &
                rmin          = R_in,               &
                rmax          = R_out,              &
                p_index       = pindex,             &
                q_index       = qindex,             &
                HoverR        = H_R,                &
                disc_mass     = disc_m,             &
                star_mass     = star_m,             &
                gamma         = gamma,              &
                particle_mass = massoftype(itype),  &
                hfact         = hfact,              &
                xyzh          = xyzh,               &
                vxyzu         = vxyzu,              &
                polyk         = polyk,              &
                inclination   = xinc,               &
                twist         = .false.,            &
                ismooth       = ismooth_edge,       &
                prefix        = fileprefix)
    if (use_dust) then
       if (use_dustfrac) then
          do i=1,npart
             call set_dustfrac(dust_to_gas_ratio,dustfrac(:,i))
          enddo
       else
          itype       = idust
          npart_inter = npartoftype(igas) + 1
          disc_mdust  = disc_m*dust_to_gas_ratio
          call set_disc(id,master   = master,             &
                      npart         = npartoftype(itype), &
                      npart_start   = npart_inter,        &
                      particle_type = itype,              &
                      rmin          = R_indust,           &
                      rmax          = R_outdust,          &
                      p_index       = pindex,             &
                      q_index       = qindex,             &
                      HoverR        = H_R,                &
                      disc_mass     = disc_mdust,         &
                      star_mass     = star_m,             &
                      gamma         = gamma,              &
                      particle_mass = massoftype(itype),  &
                      hfact         = hfact,              &
                      xyzh          = xyzh,               &
                      vxyzu         = vxyzu,              &
                      polyk         = polyk,              &
                      inclination   = xinc,               &
                      twist         = .false.,            &
                      ismooth       = ismooth_edge,       &
                      prefix        = fileprefix)
       endif
    endif
 case(1)
    call set_disc(id,master   = master,             &
                npart         = npartoftype(itype), &
                npart_start   = npart_inter,        &
                particle_type = itype,              &
                rmin          = R_in,               &
                rmax          = R_out,              &
                p_index       = pindex,             &
                q_index       = qindex,             &
                HoverR        = H_R,                &
                sig_naught    = sigma_naught,       &
                star_mass     = star_m,             &
                gamma         = gamma,              &
                particle_mass = massoftype(itype),  &
                hfact         = hfact,              &
                xyzh          = xyzh,               &
                vxyzu         = vxyzu,              &
                polyk         = polyk,              &
                inclination   = xinc,               &
                twist         = .false.,            &
                ismooth       = ismooth_edge,       &
                prefix        = fileprefix)
    if (use_dust) then
       if (use_dustfrac) then
          do i=1,npart
             call set_dustfrac(dust_to_gas_ratio,dustfrac(:,i))
          enddo
       else
          itype       = idust
          npart_inter = npartoftype(igas) + 1
          call set_disc(id,master   = master,             &
                      npart         = npartoftype(itype), &
                      npart_start   = npart_inter,        &
                      particle_type = itype,              &
                      rmin          = R_indust,           &
                      rmax          = R_outdust,          &
                      p_index       = pindex,             &
                      q_index       = qindex,             &
                      HoverR        = H_R,                &
                      sig_naught    = sigma_naughtdust,   &
                      star_mass     = star_m,             &
                      gamma         = gamma,              &
                      particle_mass = massoftype(itype),  &
                      hfact         = hfact,              &
                      xyzh          = xyzh,               &
                      vxyzu         = vxyzu,              &
                      polyk         = polyk,              &
                      inclination   = xinc,               &
                      twist         = .false.,            &
                      ismooth       = ismooth_edge,       &
                      prefix        = fileprefix)
       endif
    endif
 case(2)
    print "(a)",'Toomre Q not yet implemented'
    stop
 end select

 alpha = 0.1

 grainsizecgs = grainsizeinp

!
!--central star
!
 nptmass = 1
 xyzmh_ptmass(:,:)   = 0.
 xyzmh_ptmass(1:3,1) = 0.
 xyzmh_ptmass(4,1)   = star_m
 xyzmh_ptmass(ihacc,1)  = accr1
 xyzmh_ptmass(ihsoft,1) = accr1
 vxyz_ptmass = 0.

 print "(a,i2,a)",' --------- added ',nplanets,' planets ------------'
 period_longest = 0.
 do i=1,nplanets
    nptmass = nptmass + 1
    phi = 0.*pi/180.
    cosphi = cos(phi)
    sinphi = sin(phi)
    disc_m_within_r = 0.
    do j=1,npart
       r2 = xyzh(1,j)**2 + xyzh(2,j)**2 + xyzh(3,j)**2
       if (r2 < rplanet(i)**2) then
          disc_m_within_r = disc_m_within_r + massoftype(igas)
       endif
    enddo
    xyzmh_ptmass(1:3,nptmass)   = (/rplanet(i)*cosphi,rplanet(i)*sinphi,0./)
    xyzmh_ptmass(4,nptmass)     = mplanet(i)*jupiterm/umass
    xyzmh_ptmass(ihacc,nptmass)  = accrplanet(i) ! 0.25*au/udist
    xyzmh_ptmass(ihsoft,nptmass) = accrplanet(i) ! 0.25*au/udist
    vphi = sqrt((star_m + disc_m_within_r)/rplanet(i))
    vxyz_ptmass(1:3,nptmass)    = (/-vphi*sinphi,vphi*cosphi,0./)
    print "(a,i2,a)",       ' planet ',i,':'
    print "(a,g10.3,a)",    ' radius: ',rplanet(i)*udist/au,' AU'
    print "(a,g10.3,a,2pf7.3,a)",    ' M(<R) : ',(disc_m_within_r + star_m)*umass/solarm, &
          ' MSun, disc mass correction is ',disc_m_within_r/star_m,'%'
    print "(a,2(g10.3,a))", ' mass  : ',mplanet(i),' MJup, or ',mplanet(i)*jupiterm/solarm,' MSun'
    print "(a,2(g10.3,a))", ' period: ',2.*pi*rplanet(i)/vphi*utime/years,' years or ',2*pi*rplanet(i)/vphi,' in code units'
    omega = vphi/rplanet(i)
    period_longest = max(period_longest, 2.*pi/omega)
    print "(a,g10.3,a)",   ' resonances: 3:1: ',(sqrt(star_m)/(3.*omega))**(2./3.),' AU'
    print "(a,g10.3,a)",   '             4:1: ',(sqrt(star_m)/(4.*omega))**(2./3.),' AU'
    print "(a,g10.3,a)",   '             5:1: ',(sqrt(star_m)/(5.*omega))**(2./3.),' AU'
    print "(a,g10.3,a)",   '             9:1: ',(sqrt(star_m)/(9.*omega))**(2./3.),' AU'
 enddo
 print "(1x,45('-'))"

 !
 !--set default options for the input file
 !
 tmax = 10.*period_longest
 dtmax = 0.1*period_longest

 if (maxvxyzu==3 .and. qindex > 0.) then
    ieos = 3  ! use locally isothermal equation of state
 endif

 iexternalforce = 0
 !
 ! reset centre of mass to the origin
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 return
end subroutine setpart

subroutine write_hltauinputfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for HLTauri setup routine'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','number of gas particles',iunit)
 if (use_dust .and. .not.use_dustfrac) then
    call write_inopt(np_dust,'np_dust','number of dust particles',iunit)
 endif
 write(iunit,"(/,a)") '# units'
 call write_inopt(u_dist,'udist','distance unit in cm',iunit)
 call write_inopt(u_mass,'umass','mass unit in g',iunit)
 !--star
 write(iunit,"(/,a)") '# options for central star'
 call write_inopt(star_m,'star_m','star mass',iunit)
 call write_inopt(accr1,'accr1','star accretion radius',iunit)
 !--disc
 write(iunit,"(/,a)") '# options for accretion disc'
 call write_inopt(mass_set,'mass_set',&
     'how to set disc mass (0=disc mass, 1=surface density, 2=Toomre Q (not yet implemented))',iunit)
 select case(mass_set)
 case(0)
    call write_inopt(disc_m,'disc_m','disc mass in gas',iunit)
    if (use_dust) then
       call write_inopt(dust_to_gas_ratio,'dust_to_gas_ratio','dust to gas ratio',iunit)
    endif
 case(1)
    call write_inopt(sigma_naught,'sigma_naught','surface density at R=1 in gas',iunit)
    if (use_dust .and. .not.use_dustfrac) then
       call write_inopt(sigma_naughtdust,'sigma_naughtdust','surface density at R=1 in dust',iunit)
    endif
 case(2)
    print "(a)",'Toomre Q not yet implemented'
    stop
 end select
 call write_inopt(R_in,'R_in','inner radius',iunit)
 call write_inopt(R_out,'R_out','outer radius in gas',iunit)
 if (use_dust .and. .not.use_dustfrac) then
    call write_inopt(R_indust,'R_indust','inner radius in dust',iunit)
    call write_inopt(R_outdust,'R_outdust','outer radius in dust',iunit)
 endif
 call write_inopt(alphaSS,'alphaSS','desired alphaSS',iunit)
 call write_inopt(xinc,'xinc','inclination angle',iunit)
 call write_inopt(H_R,'H_R','H/R at R=Rin',iunit)
 call write_inopt(pindex,'pindex','p index',iunit)
 call write_inopt(qindex,'qindex','q index',iunit)
 if (use_dust) call write_inopt(grainsizeinp,'grainsizeinp','grain size (in cm)',iunit)
 !--planets
 write(iunit,"(/,a)") '# set planets'
 call write_inopt(nplanets,'nplanets','number of planets',iunit)
 do i=1,nplanets
    write(iunit,"(/,a)") '# planet:'//trim(planets(i))
    call write_inopt(mplanet(i),'mplanet'//trim(planets(i)),'planet mass (in Jupiter mass)',iunit)
    call write_inopt(rplanet(i),'rplanet'//trim(planets(i)),'planet distance from star',iunit)
    call write_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),'planet radius',iunit)
 enddo
 close(iunit)

end subroutine write_hltauinputfile

subroutine read_hltauinputfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(np,'np',db,min=0,max=maxp,errcount=nerr)
 if (use_dust .and. .not.use_dustfrac) then
    call read_inopt(np_dust,'np_dust',db,min=0,max=maxp-np,errcount=nerr)
 endif
 call read_inopt(u_dist,'udist',db,min=0.,errcount=nerr)
 call read_inopt(u_mass,'umass',db,min=0.,errcount=nerr)
 call read_inopt(star_m,'star_m',db,min=0.,errcount=nerr)
 call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
 call read_inopt(mass_set,'mass_set',db,min=0,max=2,errcount=nerr)
 select case(mass_set)
 case(0)
    call read_inopt(disc_m,'disc_m',db,min=0.,errcount=nerr)
    if (use_dust) then
       call read_inopt(dust_to_gas_ratio,'dust_to_gas_ratio',db,min=0.,errcount=nerr)
    endif
 case(1)
    call read_inopt(sigma_naught,'sigma_naught',db,min=0.,errcount=nerr)
    if (use_dust .and. .not.use_dustfrac) then
       call read_inopt(sigma_naughtdust,'sigma_naughtdust',db,min=0.,errcount=nerr)
    endif
 case(2)
    print "(a)",'Toomre Q not yet implemented'
    stop
 end select
 call read_inopt(R_in,'R_in',db,min=0.,errcount=nerr)
 call read_inopt(R_out,'R_out',db,min=R_in,errcount=nerr)
 if (use_dust .and. .not.use_dustfrac) then
    call read_inopt(R_indust,'R_indust',db,min=0.,errcount=nerr)
    call read_inopt(R_outdust,'R_outdust',db,min=R_indust,errcount=nerr)
 endif
 call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
 call read_inopt(xinc,'xinc',db,errcount=nerr)
 call read_inopt(H_R,'H_R',db,min=0.,errcount=nerr)
 call read_inopt(pindex,'pindex',db,errcount=nerr)
 call read_inopt(qindex,'qindex',db,errcount=nerr)
 if (use_dust) call read_inopt(grainsizeinp,'grainsizeinp',db,min=0.,errcount=nerr)
 call read_inopt(nplanets,'nplanets',db,min=0,max=maxplanets,errcount=nerr)
 do i=1,nplanets
    call read_inopt(mplanet(i),'mplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
    call read_inopt(rplanet(i),'rplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
    call read_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
!    print*,' read_inopt planet ',rplanet(i), mplanet(i), accrplanet(i)
 enddo
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_hltauinputfile

end module setup
