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
!   This module sets up a dusty disc
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    H_R               -- H/R at R=Rin
!    H_R_dust          -- H/R at R=Rin for dust
!    R_c               -- characteristic radius of the exponential taper
!    R_c_dust          -- characteristic radius of the exponential taper
!    R_in              -- inner radius
!    R_indust          -- inner radius in dust
!    R_out             -- outer radius in gas
!    R_outdust         -- outer radius in dust
!    accradius         -- star accretion radius
!    accretion_radius2 -- accretion radius for secondary
!    alphaSS           -- desired alphaSS
!    binary_set        -- set the binary? (0=no,1=yes)
!    disc_m            -- disc mass in gas
!    dust_to_gas_ratio -- dust to gas ratio
!    eccentricity      -- eccentricity of the binary
!    graindensinp      -- intrinsic grain density (in g/cm^3)
!    grainsizeinp      -- grain size (in cm)
!    icentralforce     -- central mass force choice (0=sink,1=external force)
!    inclination       -- disc inclination [deg] with respect to xy plane
!    iprofiledust      -- set dust surface density profile (0=power-law, 1=tapered power-law)
!    iprofilegas       -- set gas surface density profile (0=power-law, 1=tapered power-law)
!    massratio         -- mass ratio of binary
!    np                -- number of gas particles
!    np_dust           -- number of dust particles
!    nplanets          -- number of planets
!    pindex            -- p index
!    pindex_dust       -- p index
!    profile_set_dust  -- how to set dust density profile (0=equal to gas, 1=custom)
!    qindex            -- q index
!    qindex_dust       -- q index for dust
!    semimajoraxis     -- initial separation of binary (code units)
!    setplanets        -- add planets? (0=no,1=yes)
!    sigma_naught      -- Sigma0 of the gas profile Sigma = Sigma0*(R/Rc)^-p*Exp(-(R/Rc)^(2-p))
!    star_m            -- star mass
!    udist             -- distance unit in cm
!    umass             -- mass unit in g
!    xinc              -- inclination angle
!
!  DEPENDENCIES: centreofmass, dim, dust, externalforces, infile_utils, io,
!    options, part, physcon, prompting, setbinary, setdisc, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 use dim,            only:maxp,use_dust,use_dustfrac,maxalpha
 implicit none
 public :: setpart
 real :: R_in,R_out,xinc,star_m,disc_m,pindex,qindex,pindex_dust,qindex_dust,H_R,H_R_dust,sigma_naught,grainsizeinp,graindensinp
 real :: disc_mdust,dust_to_gas_ratio,R_outdust,R_indust,R_c,R_c_dust
 real :: accradius,alphaSS
 integer, parameter :: maxplanets = 9
 real    :: mplanet(maxplanets),rplanet(maxplanets),accrplanet(maxplanets),inclplan(maxplanets)
 real    :: mprimary,massratio,semimajoraxis,eccentricity,accretion_radius2,inclbin
 integer :: i,nplanets,np,np_dust,icentralforce,setplanets
 integer :: itype,npart_inter,mass_set,profile_set_dust,iprofilegas,iprofiledust,binary
 character(len=*), dimension(maxplanets), parameter :: planets = &
  (/'1','2','3','4','5','6','7','8','9' /)

 private
 real(kind=8) :: u_dist,u_mass

contains

!----------------------------------------------------------------
!
! This subroutine sets up an accretion disc with dust
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use eos,            only:qfacdisc
 use setdisc,        only:set_disc
 use setbinary,      only:set_binary
 use part,           only:nptmass,xyzmh_ptmass,maxvxyzu,vxyz_ptmass,ihacc,ihsoft,igas,idust,dustfrac
 use io,             only:master,fatal
 use options,        only:iexternalforce,ieos,alpha,icooling
 use externalforces, only:accradius1
 use units,          only:set_units,umass,udist,utime
 use physcon,        only:au,solarm,jupiterm,pi,years
 use prompting,      only:prompt
 use timestep,       only:tmax,dtmax
 use centreofmass,   only:reset_centreofmass
 use dust,           only:set_dustfrac,grainsizecgs,graindenscgs
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
 logical :: iexist,isigmapringlegas,isigmapringledust,is_binary,questplanets
 real    :: phi,vphi,sinphi,cosphi,omega,r2,disc_m_within_r,period_longest,idust_to_gas_ratio,ri,sigma_naughtdust,period
 real    :: sini,cosi,polyk_dust
 integer :: ierr,j

 !
 !--Set problem parameters
 !

 u_dist = 1.*au
 u_mass = solarm
 call set_units(dist=u_dist,mass=u_mass,G=1.d0)

 np = 500000
 if (use_dust .and. .not. use_dustfrac) then
    np_dust = 100000
 else
    np_dust = 0
 endif
 npartoftype(:)     = 0
 npartoftype(igas)  = np
 if (use_dust .and. .not. use_dustfrac) npartoftype(idust) = np_dust
 star_m             = 1.*solarm/umass
 accradius          = 1.*au/udist
 R_in               = accradius
 R_out              = 150.*au/udist
 R_indust           = accradius
 R_outdust          = 150.*au/udist
 mass_set           = 0
 profile_set_dust   = 0
 disc_m             = 0.05*solarm/umass
 icentralforce      = 0
 if (use_dust) then
    dust_to_gas_ratio = 0.01
 else
    dust_to_gas_ratio = 0.
 endif
 sigma_naught       = 1.009E-07
 alphaSS            = 0.005
 xinc               = 0.0
 H_R                = 0.05
 H_R_dust           = H_R
 isigmapringlegas   = .false.
 isigmapringledust  = .false.
 is_binary          = .false.
 questplanets       = .false.
 iprofilegas        = 0
 iprofiledust       = 0
 pindex             = 1.
 pindex_dust        = 1.
 R_c                = R_out
 R_c_dust           = R_outdust
 qindex             = 0.25
 qindex_dust        = qindex
 mprimary           = star_m
 massratio          = 1.
 semimajoraxis      = 1.*au/udist
 eccentricity       = 0.
 accretion_radius2  = accradius
 inclbin          = 0.
 setplanets         = 0
 nplanets           = 1
 mplanet(:)         = 1.d0
 accrplanet(:)      = 0.25*au/udist
 do i=1,maxplanets
    rplanet(i) = 10.*i*au/udist
    inclplan(i) = 0.
 enddo

 if (maxvxyzu==3 .and. qindex > 0.) then
    ieos = 3  ! use locally isothermal equation of state
    gamma = 1.0
 else
    ieos = 2
    gamma = 5./3.
    icooling = 1
 endif
 time  = 0.
 hfact = 1.2

 grainsizeinp = 0.1
 graindensinp = 3.

 filename=trim(fileprefix)//'.setup'

 print "(/,65('-'),2(/,a),/,65('-'),/)",&
   ' Welcome to the dustydisc setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_dustydiscinputfile(filename,ierr)
    if (id==master) call write_dustydiscinputfile(filename)
    call set_units(dist=u_dist,mass=u_mass,G=1.d0)
    if (ierr /= 0) then
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
    call prompt('Enter accretion radius for central star (code units)',accradius,0.)
    !
    !--set gas disc options
    !
    call prompt('Enter inner radius (code units) of protoplanetary disc in gas',R_in,0.01)
    call prompt('Enter outer radius (code units) of protoplanetary disc in gas',R_out,R_in)
    call prompt('Do you want to taper the power-law gas surface density profile by an exponential function?',isigmapringlegas)
    call prompt('How do you want to set the gas disc mass? (0=disc mass, 1=surface density, 2=Toomre Q (not yet implemented))'&
                ,mass_set,0,2)
    select case(mass_set)
    case(0)
       call prompt('Enter mass (code units) of protoplanetary disc in gas',disc_m,0.0)
       if (isigmapringlegas) then
          iprofilegas = 1
          call prompt('Enter p index of the gas surface density profile Sigma = Sigma0*(R/Rc)^-p*Exp(-(R/Rc)^(2-p))',pindex)
          call prompt('Enter the characteristic radius Rc of the exponential taper',R_c,R_in,R_out)
       else
          iprofilegas = 0
          call prompt('Enter p index of the power-law gas surface density profile',pindex)
       endif
    case(1)
       if (isigmapringlegas) then
          iprofilegas = 1
          call prompt('Enter the Sigma0 of the gas profile Sigma = Sigma0*(R/Rc)^-p*Exp(-(R/Rc)^(2-p))',sigma_naught,0.0)
          call prompt('Enter p index of the gas surface density profile Sigma = Sigma0*(R/Rc)^-p*Exp(-(R/Rc)^(2-p))',pindex)
          call prompt('Enter the characteristic radius Rc of of the exponential taper',R_c,R_in,R_out)
       else
          iprofilegas = 0
          call prompt('Enter suface density (code units) at R=1 (code units) in gas',sigma_naught,0.0)
          call prompt('Enter p index of the power-law gas surface density profile',pindex)
       endif
    case(2)
       print "(a)",'Toomre Q not yet implemented'
       stop
    end select
    call prompt('Enter q index of temperature profile cs = cs0*R^-q',qindex,0.)
    call prompt('Enter desired alpha_SS',alphaSS,0.)
    call prompt('Enter the disc inclination [deg] with respect to xy plane',xinc,0.,180.)
    call prompt('Enter H/R at R=Rin',H_R,0.)
    call prompt('Do you want to set a binary?',is_binary)
    if (is_binary) then
       binary=1
       call prompt('Enter mass ratio of binary',massratio,0.)
       call prompt('Enter initial separation of binary (code units)',semimajoraxis)
       call prompt('Enter the eccentricity of the binary',eccentricity,0.,1.)
       call prompt('Enter accretion radius for secondary',accretion_radius2,0.,semimajoraxis)
       call prompt('Enter the binary inclination [deg] with respect to xy plane',inclbin,0.,180.)
    else
       binary=0
       call prompt('How is the central object defined? (0=sink particle, 1=external potential)',icentralforce,0,1)
    endif
    !
    !--set dust disc options
    !
    if (use_dust) then
       call prompt('Enter dust to gas ratio',dust_to_gas_ratio,0.)
       call prompt('How do you want to set the dust density profile? (0=equal to gas, 1=custom)',profile_set_dust,0,1)
       select case(profile_set_dust)
       case(0)
          R_indust          = R_in
          R_outdust         = R_out
          pindex_dust       = pindex
          qindex_dust       = qindex
          H_R_dust          = H_R
          isigmapringledust = isigmapringlegas
          iprofiledust      = iprofilegas
          if (isigmapringledust) R_c_dust = R_c
          print "(a)",'dust and gas have the same surface density profile'
       case(1)
          if (use_dustfrac) then
             call prompt('Enter inner radius (code units) of protoplanetary disc in dust',R_indust,R_in)
             call prompt('Enter outer radius (code units) of protoplanetary disc in dust',R_outdust,R_indust,R_out)
          else
             call prompt('Enter inner radius (code units) of protoplanetary disc in dust',R_indust,0.01)
             call prompt('Enter outer radius (code units) of protoplanetary disc in dust',R_outdust,R_indust)
             call prompt('Enter q index of temperature profile cs = cs0*R^-q for dust',qindex_dust,0.)
             call prompt('Enter H/R at R=Rin for dust',H_R_dust,0.)
          endif
          call prompt('Do you want to taper the power-law dust surface density profile by an exponential function?',&
                     isigmapringledust)
          if (isigmapringledust) then
             iprofiledust = 1
             call prompt('Enter p index of the dust surface density profile Sigma = Sigma0*(R/Rc)^-p*Exp(-(R/Rc)^(2-p))',&
                           pindex_dust,0.)
             call prompt('Enter the characteristic radius Rc of the exponential taper',R_c_dust,R_indust,R_outdust)
          else
             iprofiledust = 0
             call prompt('Enter p index of the power-law dust surface density profile',pindex_dust)
          endif
       end select
       call prompt('Enter grain size in cm',grainsizeinp,0.)
       call prompt('Enter intrinsic grain density in g/cm^3',graindensinp,0.)
    endif
    !
    !--add the planets
    !
    call prompt('Do you want to add planets?',questplanets)
    if(questplanets)then
       setplanets = 1
       call prompt('Enter the number of planets',nplanets,0,maxplanets)
       do i=1,nplanets
          call prompt('Enter mass (in Jupiter mass) of planet '//trim(planets(i))//' :',mplanet(i),0.,star_m*umass/jupiterm)
          call prompt('Enter distance from the central star (code units) of planet '//trim(planets(i))//' :',rplanet(i),0.,R_out)
          call prompt('Enter planet inclination [deg] with respect to xy plane '//trim(planets(i))//' :',inclplan(i),0.,180.)
          call prompt('Enter accretion radius (code units) of planet '//trim(planets(i))//' :',accrplanet(i),0.)
       enddo
    else
       setplanets = 0
    endif

    !
    !--write default input file
    !
    call write_dustydiscinputfile(filename)

    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'

    stop
 else
    stop
 endif

 xinc = xinc*(pi/180.0)
 if(maxalpha==0) alpha = alphaSS

 if (profile_set_dust==0)then
    R_indust     = R_in
    R_outdust    = R_out
    pindex_dust  = pindex
    iprofiledust = iprofilegas
    if (iprofiledust==1) R_c_dust = R_c
 endif

 if(binary==1)then
    nptmass  = 0
    mprimary = star_m
    star_m   = mprimary+mprimary*massratio
    print*,'Setup the binary' !--the stars barycentre is put on the origin
    call set_binary(mprimary,massratio,semimajoraxis,eccentricity, &
               accradius,accretion_radius2,xyzmh_ptmass,vxyz_ptmass,nptmass,incl=inclbin)
 else
    print*,'Disc around a single star '
    iexternalforce = icentralforce
    if(iexternalforce==0) then
       print "(a)", 'Central object represented by a sink at the system origin'
       nptmass                      = 1
       xyzmh_ptmass(:,:)            = 0.
       xyzmh_ptmass(1:3,nptmass)    = 0.
       xyzmh_ptmass(4,nptmass)      = star_m
       xyzmh_ptmass(ihacc,nptmass)  = accradius
       xyzmh_ptmass(ihsoft,nptmass) = accradius
       vxyz_ptmass                  = 0.
    else
       print "(a)", 'Central object represented by external force with accretion boundary'
       print*, ' Accretion Radius: ', accradius
       accradius1 = accradius
    endif
 endif

!
! compute the mass of the gas disc if mass_set=1
!
 if(mass_set==1) call sigma_naught2mass(sigma_naught,iprofilegas,pindex,R_in,R_out,R_c,disc_m)

 npartoftype(igas)  = np
 if (use_dust) then
    npartoftype(idust) = np_dust
    disc_mdust = disc_m*dust_to_gas_ratio
 endif
 npart              = sum(npartoftype)


 npart_inter = 1
 itype       = igas
 if(use_dust .and. use_dustfrac) then
    if(mass_set==0) call mass2sigma_naught(sigma_naught,iprofilegas,pindex,R_in,R_out,R_c,disc_m)
    call mass2sigma_naught(sigma_naughtdust,iprofiledust,pindex_dust,R_indust,R_outdust,R_c_dust,disc_m*dust_to_gas_ratio)
    call set_disc(id,master   = master,             &
             mixture          = .true.,             &
             npart            = npartoftype(itype), &
             npart_start      = npart_inter,        &
             rmin             = R_in,               &
             rmax             = R_out,              &
             rmindust         = R_indust,           &
             rmaxdust         = R_outdust,          &
             indexprofile     = iprofilegas,        &
             indexprofiledust = iprofiledust,       &
             rc               = R_c,                &
             rcdust           = R_c_dust,           &
             p_index          = pindex,             &
             p_indexdust      = pindex_dust,        &
             q_index          = qindex,             &
             HoverR           = H_R,                &
             disc_mass        = disc_m,             &
             disc_massdust    = disc_mdust,         &
             star_mass        = star_m,             &
             gamma            = gamma,              &
             particle_mass    = massoftype(itype),  &
             hfact            = hfact,              &
             xyzh             = xyzh,               &
             vxyzu            = vxyzu,              &
             polyk            = polyk,              &
             alpha            = alpha,              &
             ismooth          = .false.,            &
             inclination      = xinc,               &
             twist            = .false.,            &
             prefix           = fileprefix)
    do i=1,npart
       ri = sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       if (ri<R_indust .or. ri>R_outdust) then
          idust_to_gas_ratio = 0.
       else
          call get_d2gratio(idust_to_gas_ratio,ri,iprofilegas,iprofiledust,sigma_naught,&
                            sigma_naughtdust,pindex,pindex_dust,R_c,R_c_dust)
       endif
       call set_dustfrac(idust_to_gas_ratio,dustfrac(i))
    enddo
 else
    call set_disc(id,master  = master,             &
             npart           = npartoftype(itype), &
             npart_start     = npart_inter,        &
             particle_type   = itype,              &
             rmin            = R_in,               &
             rmax            = R_out,              &
             indexprofile    = iprofilegas,        &
             rc              = R_c,                &
             p_index         = pindex,             &
             q_index         = qindex,             &
             HoverR          = H_R,                &
             disc_mass       = disc_m,             &
             star_mass       = star_m,             &
             gamma           = gamma,              &
             particle_mass   = massoftype(itype),  &
             hfact           = hfact,              &
             xyzh            = xyzh,               &
             vxyzu           = vxyzu,              &
             polyk           = polyk,              &
             alpha           = alpha,              &
             ismooth         = .false.,            &
             inclination     = xinc,               &
             twist           = .false.,            &
             prefix          = fileprefix)
    if (use_dust) then
       itype       = idust
       npart_inter = npartoftype(igas) + 1
       call set_disc(id,master   = master,             &
                   npart         = npartoftype(itype), &
                   npart_start   = npart_inter,        &
                   particle_type = itype,              &
                   rmin          = R_indust,           &
                   rmax          = R_outdust,          &
                   indexprofile  = iprofiledust,       &
                   rc            = R_c_dust,           &
                   p_index       = pindex_dust,        &
                   q_index       = qindex_dust,        &
                   HoverR        = H_R_dust,           &
                   disc_mass     = disc_mdust,         &
                   star_mass     = star_m,             &
                   gamma         = gamma,              &
                   particle_mass = massoftype(itype),  &
                   hfact         = hfact,              &
                   xyzh          = xyzh,               &
                   vxyzu         = vxyzu,              &
                   polyk         = polyk_dust,         &
                   ismooth       = .false.,            &
                   inclination   = xinc,               &
                   twist         = .false.,            &
                   prefix        = fileprefix)
    endif
 endif

 qfacdisc = qindex ! reset qfacdisc to gas disc value

 if(use_dust)then
    grainsizecgs = grainsizeinp
    graindenscgs = graindensinp
 endif

 if(setplanets==1)then
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
             if (.not.use_dust .or. (use_dust .and. use_dustfrac))then
                disc_m_within_r = disc_m_within_r + massoftype(igas)
             else
                if(j<=np)then
                   disc_m_within_r = disc_m_within_r + massoftype(igas)
                else
                   disc_m_within_r = disc_m_within_r + massoftype(idust)
                endif
             endif
          endif
       enddo
       if(nplanets>1)then
          do j=1,nplanets
             if(rplanet(j)<rplanet(i)) disc_m_within_r = disc_m_within_r + mplanet(j)*jupiterm/solarm
          enddo
       endif
       xyzmh_ptmass(1:3,nptmass)    = (/rplanet(i)*cosphi,rplanet(i)*sinphi,0./)
       xyzmh_ptmass(4,nptmass)      = mplanet(i)*jupiterm/umass
       xyzmh_ptmass(ihacc,nptmass)  = accrplanet(i) ! 0.25*au/udist
       xyzmh_ptmass(ihsoft,nptmass) = accrplanet(i) ! 0.25*au/udist
       vphi                         = sqrt((star_m + disc_m_within_r)/rplanet(i))
       vxyz_ptmass(1:3,nptmass)     = (/-vphi*sinphi,vphi*cosphi,0./)
       ! Rotate positions and velocities
       if(inclplan(i)  /=  0.) then
          cosi=cos(inclplan(i)*pi/180.)
          sini=sin(inclplan(i)*pi/180.)
          call rotate(xyzmh_ptmass(1:3,nptmass),cosi,sini)
          call rotate(vxyz_ptmass(1:3,nptmass),cosi,sini)
       endif
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
    tmax  = 100.*period_longest
    dtmax = 0.5*period_longest
 else if(binary==1) then
    period = sqrt(4.*pi**2*semimajoraxis**3/(mprimary*(1. + massratio)))
    tmax   = 100.*period
    dtmax  = 0.5*period
 endif

 !
 ! reset centre of mass to the origin
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 return
end subroutine setpart

subroutine write_dustydiscinputfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for dustydisc setup routine'
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
 call write_inopt(accradius,'accradius','star accretion radius',iunit)
 !--disc
 write(iunit,"(/,a)") '# options for gas accretion disc'
 call write_inopt(R_in,'R_in','inner radius',iunit)
 call write_inopt(R_out,'R_out','outer radius in gas',iunit)
 call write_inopt(iprofilegas,'iprofilegas','set gas surface density profile (0=power-law, 1=tapered power-law)',iunit)
 call write_inopt(mass_set,'mass_set',&
     'how to set gas disc mass (0=disc mass, 1=surface density, 2=Toomre Q (not yet implemented))',iunit)
 select case(mass_set)
 case(0)
    call write_inopt(disc_m,'disc_m','disc mass in gas',iunit)
    call write_inopt(pindex,'pindex','p index',iunit)
    if (iprofilegas==1) call write_inopt(R_c,'R_c','characteristic radius of the exponential taper',iunit)
 case(1)
    call write_inopt(pindex,'pindex','p index',iunit)
    if (iprofilegas==0)then
       call write_inopt(sigma_naught,'sigma_naught','Sigma0 of the gas profile Sigma = Sigma0*R^-p',iunit)
    else
       call write_inopt(sigma_naught,'sigma_naught','Sigma0 of the gas profile Sigma = Sigma0*(R/Rc)^-p*Exp(-(R/Rc)^(2-p))',iunit)
       call write_inopt(R_c,'R_c','characteristic radius of the exponential taper',iunit)
    endif
 case(2)
    print "(a)",'Toomre Q not yet implemented'
    stop
 end select
 call write_inopt(qindex,'qindex','q index',iunit)
 if(maxalpha==0) call write_inopt(alphaSS,'alphaSS','desired alphaSS',iunit)
 call write_inopt(xinc,'xinc','inclination angle',iunit)
 call write_inopt(H_R,'H_R','H/R at R=Rin',iunit)
 call write_inopt(binary,'binary_set','set the binary? (0=no,1=yes)',iunit)
 call write_inopt(setplanets,'setplanets','add planets? (0=no,1=yes)',iunit)
 if(binary==1) then
    write(iunit,"(/,a)") '# options for the binary'
    call write_inopt(massratio,'massratio','mass ratio of binary',iunit)
    call write_inopt(semimajoraxis,'semimajoraxis','initial separation of binary (code units)',iunit)
    call write_inopt(eccentricity,'eccentricity','eccentricity of the binary',iunit)
    call write_inopt(accretion_radius2,'accretion_radius2','accretion radius for secondary',iunit)
    call write_inopt(inclbin,'inclination','disc inclination [deg] with respect to xy plane',iunit)
 else
    call write_inopt(icentralforce,'icentralforce', 'central mass force choice (0=sink,1=external force)',iunit)
 endif
 if (use_dust) then
    write(iunit,"(/,a)") '# options for dust accretion disc'
    call write_inopt(dust_to_gas_ratio,'dust_to_gas_ratio','dust to gas ratio',iunit)
    call write_inopt(profile_set_dust,'profile_set_dust','how to set dust density profile (0=equal to gas, 1=custom)',iunit)
    select case(profile_set_dust)
    case(1)
       call write_inopt(iprofiledust,'iprofiledust','set dust surface density profile (0=power-law, 1=tapered power-law)',iunit)
       call write_inopt(pindex_dust,'pindex_dust','p index',iunit)
       if (.not. use_dustfrac) then
          call write_inopt(qindex_dust,'qindex_dust','q index',iunit)
          call write_inopt(H_R_dust,'H_R_dust','H/R at R=Rin',iunit)
       endif
       if (iprofiledust==1) call write_inopt(R_c_dust,'R_c_dust','characteristic radius of the exponential taper',iunit)
       call write_inopt(R_indust,'R_indust','inner radius in dust',iunit)
       call write_inopt(R_outdust,'R_outdust','outer radius in dust',iunit)
    end select
    call write_inopt(grainsizeinp,'grainsizeinp','grain size (in cm)',iunit)
    call write_inopt(graindensinp,'graindensinp','intrinsic grain density (in g/cm^3)',iunit)
 endif
 if(setplanets==1)then
    !--planets
    write(iunit,"(/,a)") '# set planets'
    call write_inopt(nplanets,'nplanets','number of planets',iunit)
    do i=1,nplanets
       write(iunit,"(/,a)") '# planet:'//trim(planets(i))
       call write_inopt(mplanet(i),'mplanet'//trim(planets(i)),'planet mass (in Jupiter mass)',iunit)
       call write_inopt(rplanet(i),'rplanet'//trim(planets(i)),'planet distance from star',iunit)
       call write_inopt(inclplan(i),'inclplanet'//trim(planets(i)),'planet inclination [deg] with respect to xy plane',iunit)
       call write_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),'planet radius',iunit)
    enddo
 endif
 close(iunit)

end subroutine write_dustydiscinputfile

subroutine read_dustydiscinputfile(filename,ierr)
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
 call read_inopt(accradius,'accradius',db,min=0.,errcount=nerr)
 call read_inopt(R_in,'R_in',db,min=0.,errcount=nerr)
 call read_inopt(R_out,'R_out',db,min=R_in,errcount=nerr)
 call read_inopt(iprofilegas,'iprofilegas',db,min=0,max=1,errcount=nerr)
 call read_inopt(mass_set,'mass_set',db,min=0,max=2,errcount=nerr)
 select case(mass_set)
 case(0)
    call read_inopt(disc_m,'disc_m',db,min=0.,errcount=nerr)
    call read_inopt(pindex,'pindex',db,errcount=nerr)
    if (iprofilegas==1) call read_inopt(R_c,'R_c',db,min=0.,errcount=nerr)
 case(1)
    call read_inopt(sigma_naught,'sigma_naught',db,min=0.,errcount=nerr)
    call read_inopt(pindex,'pindex',db,errcount=nerr)
    if (iprofilegas==1) call read_inopt(R_c,'R_c',db,min=0.,errcount=nerr)
 case(2)
    print "(a)",'Toomre Q not yet implemented'
    stop
 end select
 call read_inopt(qindex,'qindex',db,errcount=nerr)
 if(maxalpha==0) call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
 call read_inopt(xinc,'xinc',db,min=0.,max=180.,errcount=nerr)
 call read_inopt(H_R,'H_R',db,min=0.,errcount=nerr)
 call read_inopt(binary,'binary_set',db,min=0,max=1,errcount=nerr)
 call read_inopt(setplanets,'setplanets',db,min=0,max=1,errcount=nerr)
 if (binary==1)then
    call read_inopt(massratio,'massratio',db,min=0.,errcount=nerr)
    call read_inopt(semimajoraxis,'semimajoraxis',db,min=0.,errcount=nerr)
    call read_inopt(eccentricity,'eccentricity',db,min=0.,max=1.,errcount=nerr)
    call read_inopt(accretion_radius2,'accretion_radius2',db,min=0.,max=semimajoraxis,errcount=nerr)
    call read_inopt(inclbin,'inclination',db,min=0.,max=180.,errcount=nerr)
 else
    call read_inopt(icentralforce,'icentralforce',db,min=0,max=1, errcount=nerr)
 endif
 if (use_dust) then
    call read_inopt(dust_to_gas_ratio,'dust_to_gas_ratio',db,min=0.,errcount=nerr)
    call read_inopt(profile_set_dust,'profile_set_dust',db,min=0,max=1,errcount=nerr)
    select case(profile_set_dust)
    case(1)
       call read_inopt(iprofiledust,'iprofiledust',db,min=0,max=1,errcount=nerr)
       call read_inopt(pindex_dust,'pindex_dust',db,errcount=nerr)
       if (.not. use_dustfrac) then
          call read_inopt(qindex_dust,'qindex_dust',db,err=ierr,errcount=nerr)
          if (ierr /= 0) qindex_dust = qindex
          call read_inopt(H_R_dust,'H_R_dust',db,min=0.,err=ierr,errcount=nerr)
          if (ierr /= 0) H_R_dust = H_R
       endif
       if (iprofiledust==1) call read_inopt(R_c_dust,'R_c_dust',db,min=0.,errcount=nerr)
       if (use_dustfrac) then
          call read_inopt(R_indust,'R_indust',db,min=R_in,errcount=nerr)
          call read_inopt(R_outdust,'R_outdust',db,min=R_indust,max=R_out,errcount=nerr)
       else
          call read_inopt(R_indust,'R_indust',db,min=0.,errcount=nerr)
          call read_inopt(R_outdust,'R_outdust',db,min=R_indust,errcount=nerr)
       endif
    end select
    call read_inopt(grainsizeinp,'grainsizeinp',db,min=0.,errcount=nerr)
    call read_inopt(graindensinp,'graindensinp',db,min=0.,errcount=nerr)
 endif
 if(setplanets==1)then
    call read_inopt(nplanets,'nplanets',db,min=0,max=maxplanets,errcount=nerr)
    do i=1,nplanets
       call read_inopt(mplanet(i),'mplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
       call read_inopt(rplanet(i),'rplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
       call read_inopt(inclplan(i),'inclplanet'//trim(planets(i)),db,min=0.,max=180.,errcount=nerr)
       call read_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
!       print*,' read_inopt planet ',rplanet(i), mplanet(i), accrplanet(i)
    enddo
 endif
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_dustydiscinputfile

subroutine mass2sigma_naught(sigma_naught,iprofile,pindex,R_in,R_out,R_c,disc_m)
 use physcon, only:pi
 real,           intent(in)  :: disc_m,pindex,R_in,R_out
 real, optional, intent(in)  :: R_c
 real,           intent(out) :: sigma_naught
 integer,        intent(in)  :: iprofile

 if(iprofile==0) then
    sigma_naught = disc_m/(((2.*pi)/(2.-pindex))*((R_out)**(2.-pindex)-(R_in)**(2.-pindex)))
 else
    sigma_naught = disc_m/(((2.*pi*(R_c**2))/(2.-pindex))*(exp(-(R_in/R_c)**(2.-pindex))-exp(-(R_out/R_c)**(2.-pindex))))
 endif

end subroutine mass2sigma_naught

subroutine sigma_naught2mass(sigma_naught,iprofile,pindex,R_in,R_out,R_c,disc_m)
 use physcon, only:pi
 real,           intent(in)  :: sigma_naught,pindex,R_in,R_out
 real, optional, intent(in)  :: R_c
 real,           intent(out) :: disc_m
 integer,        intent(in)  :: iprofile

 if(iprofile==0) then
    disc_m = sigma_naught*((2.*pi)/(2.-pindex))*((R_out)**(2.-pindex)-(R_in)**(2.-pindex))
 else
    disc_m = sigma_naught*((2.*pi*(R_c**2))/(2.-pindex))*(exp(-(R_in/R_c)**(2.-pindex))-exp(-(R_out/R_c)**(2.-pindex)))
 endif

end subroutine sigma_naught2mass

subroutine get_d2gratio(id2g,ri,iprofilegas,iprofiledust,sigma_naught,sigma_naughtdust,pindex,pindex_dust,R_c,R_c_dust)
 real,           intent(in)  :: sigma_naught,sigma_naughtdust,pindex,pindex_dust,ri
 real, optional, intent(in)  :: R_c,R_c_dust
 real,           intent(out) :: id2g
 integer,        intent(in)  :: iprofilegas,iprofiledust

 if (iprofilegas==0)then
    if (iprofiledust==0)then
       id2g=(sigma_naughtdust/sigma_naught)*(ri**(-pindex_dust))/(ri**(-pindex))
    else
       id2g=(sigma_naughtdust/sigma_naught)*&
             (((ri/R_c_dust)**(-pindex_dust))*exp(-(ri/R_c_dust)**(2.-pindex_dust)))/(ri**(-pindex))
    endif
 else
    if (iprofiledust==0)then
       id2g=(sigma_naughtdust/sigma_naught)*(ri**(-pindex_dust))/(((ri/R_c)**(-pindex))*exp(-(ri/R_c)**(2.-pindex)))
    else
       id2g=(sigma_naughtdust/sigma_naught)*(((ri/R_c_dust)**(-pindex_dust))*exp(-(ri/R_c_dust)**(2.-pindex_dust)))/&
            (((ri/R_c)**(-pindex))*exp(-(ri/R_c)**(2.-pindex)))
    endif
 endif

end subroutine get_d2gratio

pure subroutine rotate(xyz,cosi,sini)
 real, intent(inout) :: xyz(3)
 real, intent(in)    :: cosi,sini
 real :: xi,yi,zi

 xi = xyz(1)
 yi = xyz(2)
 zi = xyz(3)
 xyz(1) =  xi*cosi + zi*sini
 xyz(2) =  yi
 xyz(3) = -xi*sini + zi*cosi

end subroutine rotate

end module setup
