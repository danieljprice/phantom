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
!   This module sets up general binary accretion discs
!   Written by Daniel Price, Giuseppe Lodato and Chris Nixon
!   Cambridge, June 2011
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    accr1             -- primary accretion radius
!    accr2             -- secondary accretion radius
!    alphaSS           -- desired Shakura-Sunyaev alpha viscosity parameter
!    binary_argperi    -- binary angle for argument of periapsis (deg)
!    binary_inc        -- binary inclination in degrees
!    binary_posang     -- binary position angle of ascending node (deg)
!    binary_separation -- binary separation
!    deltat            -- output interval as fraction of binary orbital period
!    dist_unit         -- distance unit (e.g. au)
!    dust_to_gas       -- initial dust-to-gas ratio
!    ecc               -- binary eccentricity
!    m1                -- primary mass
!    m2                -- secondary mass
!    mass_unit         -- mass unit (e.g. solarm)
!    norbits           -- maximum number of binary orbits
!    np                -- number of particles
!
!  DEPENDENCIES: centreofmass, dim, dust, eos, infile_utils, io, options,
!    part, physcon, prompting, setbinary, setdisc, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 use dim, only:use_dustfrac
 implicit none
 public :: setpart
 !--private module variables
 real :: m1,m2,ecc,binary_a,binary_inc,binary_posang,binary_argperi,accr1,accr2,alphaSS,deltat
 logical :: iuse_disc(3),ismooth_edge(3)
 real :: R_in(3),R_out(3),HoverR(3),disc_mass(3),p_index(3),q_index(3),xinc(3)
 real :: dust_to_gas
 integer :: np,norbits
 character(len=*), dimension(3), parameter :: disctype = &
  (/'binary   ',&
    'primary  ', &
    'secondary'/)
 private
 character(len=20) :: dist_unit,mass_unit
 real(kind=8) :: udist,umass

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc, only:set_disc
 use units,   only:set_units,select_unit,umass,utime
 use io,      only:master,fatal,error,warning
 use physcon, only:au,solarm,pi,years
 use options,        only:iexternalforce,alpha,ieos
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,maxvxyzu,igas,dustfrac
 use setbinary,      only:set_binary,Rochelobe_estimate,get_a_from_period
 use prompting,      only:prompt
 use eos,            only:isink
 use timestep,       only:dtmax,tmax
 use centreofmass,   only:reset_centreofmass
 use dust,           only:set_dustfrac
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=100) :: filename
 real :: Rochelobe,totmass,starmass,Rochesizei,minq,period
 integer :: i,ntot,npindisc,j,itest,ierr
 logical :: iexist
 real :: xorigini(3),vorigini(3),alpha_returned(3)
 integer :: number_of_discs
!
!--firstly, setup the binary
!
 alphaSS = 0.1
 alpha_returned(:) = 0.0
 np = size(xyzh(1,:))
 filename=trim(fileprefix)//'.setup'
 dist_unit = '100*au'
 mass_unit = 'solarm'
 norbits = 10
 deltat  = 0.1
 dust_to_gas = 0.01

 print "(/,65('-'),2(/,a),/,65('-'),/)",&
   ' Welcome to the Ultimate Binary Disc Setup Routine^TM', &
   ' Brought to you by Daniel Price, Chris Nixon and Giuseppe Lodato'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_binaryinputfile(filename,ierr)
    if (ierr /= 0) then
       if (id==master) call write_binaryinputfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
       call select_unit(mass_unit,umass,ierr)
       if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    enddo
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo
    !
    !--units
    !
    call set_units(dist=udist,mass=umass,G=1.d0)
    !
    !--set default options
    !
    m1 = 1.
    call prompt('Enter primary mass (code units)',m1,0.)
    m2 = 0.1
    call prompt('Enter secondary mass (code units)',m2,0.,m1)
    binary_a = 1.
    call prompt('Enter the binary semi-major axis',binary_a,0.)
    ecc = 0.
    call prompt('Enter the eccentricity of the binary',ecc,0.,1.)
    binary_inc = 0.
    call prompt('Enter the inclination of the binary in degrees',binary_inc,-180.,180.)
    binary_posang = 0.
    call prompt('Enter position angle of the ascending node in degrees',binary_posang,-180.,180.)
    binary_argperi = 0.
    call prompt('Enter the angle of the argument of periapsis in degrees',binary_argperi,-180.,180.)

    accr1 = 0.25*binary_a
    call prompt('Enter accretion radius for primary (can be adjusted later)',accr1,0.,binary_a)
    accr2 = 0.25*binary_a
    call prompt('Enter accretion radius for secondary (can be adjusted later)',accr2,0.,binary_a)

    xinc(:) = 0.0
    p_index(:) = 1.5
    q_index(:) = 0.75
    iuse_disc = .false.
    ismooth_edge = .false.

    Rochelobe = Rochelobe_estimate(m1,m2,binary_a)
!    R_in  = (/2.**(2./3.)*binary_separation,accr,accr/)
    R_in = (/binary_a + Rochelobe + accr1,accr1,accr2/)
    R_out = (/10.*R_in(1),binary_a - Rochelobe - accr1,Rochelobe - accr2/)
    disc_mass = (/1.e-2*m1,1.e-3*m1,1.e-4*m1/)
    HoverR(:) = 0.1

    do while(.not.(any(iuse_disc)))
       iuse_disc(1) = .true.
       do i=1,3
          call prompt('Do you want a circum'//trim(disctype(i))//' disc?',iuse_disc(i))
          if (iuse_disc(i)) then
             if (i==1) then
                print "(a,es10.3,a)",' Recommended minimum value of Rin = ',R_in(i),' (a + Rochelobe + acc. radius)'
                call prompt('Enter inner radius of circum'//trim(disctype(i))//' disc',R_in(i),binary_a)
                call prompt('Enter outer radius of circum'//trim(disctype(i))//' disc',R_out(i),R_in(i))
             else
                print "(a,es10.3)",' Setting R_in for circum'//trim(disctype(i))//' to ',R_in(i)
                print "(a,es10.3)",' Setting R_out for circum'//trim(disctype(i))//' to ',R_out(i)
             endif
             call prompt('Enter mass of circum'//trim(disctype(i))//' disc',disc_mass(i),0.)
             call prompt('Enter inclination angle to the z=0 plane in degrees',xinc(i))
             call prompt('Do you want Sigma to drop to zero at the inner edge',ismooth_edge(i))
             call prompt('Enter p index of surface density profile Sigma = Sigma0*R^-p',p_index(i),0.)
             call prompt('Enter q index of temperature profile cs = cs0*R^-q',q_index(i),0.)
             call prompt('Enter H/R at R=Rin',HoverR(i),0.)
             if (maxvxyzu==3) then
                do j=1,i-1
                   if (iuse_disc(j) .and. any(q_index(1:j) > 0.)) then
                      call error('setup','q_index differs between discs but '//&
                                 'not storing thermal energy: set maxvxyzu=4 in dim file')
                   endif
                enddo
             endif
          endif
       enddo
       if (.not.any(iuse_disc)) call fatal('setup','need to setup at least one disc!')
    enddo
    if (any(iuse_disc)) call prompt('Enter desired alpha_SS',alphaSS,0.)

    call prompt('Enter total number of gas particles',np,0,size(xyzh(1,:)))

    call prompt('Enter time between dumps as fraction of binary period',deltat,0.)
    call prompt('Enter number of binary periods to simulate',norbits,0)
    if (use_dustfrac) call prompt('Enter initial dust-to-gas ratio',dust_to_gas,0.)
    !
    !--write default input file
    !
    call write_binaryinputfile(filename)

    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
    stop
 else
    stop
 endif
 !
 ! units
 !
 call select_unit(mass_unit,umass,ierr)
 if (ierr /= 0) stop ' ERROR: mass unit not recognised'

 call select_unit(dist_unit,udist,ierr)
 if (ierr /= 0) stop ' ERROR: distance unit not recognised'

 call set_units(dist=udist,mass=umass,G=1.d0)

 Rochelobe = Rochelobe_estimate(m1,m2,binary_a)

 period = sqrt(4.*pi**2*binary_a**3/(m1 + m2))
 print*,' PERIOD = ',period,' or ',period*(utime/years),' YEARS'

! Convert inclination into radians
 xinc(:) = xinc(:)*pi/180.0
!
!--add a binary
!
 nptmass = 0

 call set_binary(m1,massratio=m2/m1,semimajoraxis=binary_a,eccentricity=ecc, &
                 posang_ascnode=binary_posang,arg_peri=binary_argperi,incl=binary_inc,&
                 accretion_radius1=accr1,accretion_radius2=accr2,&
                 xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass)

 gamma = 1.0
 alpha = alphaSS ! this gets reset in call to set_disc
 np = min(np,size(xyzh(1,:)))
 print*,' using npart = ',np
 hfact = 1.2
 time  = 0.

! How many discs do we have
 itest = 0
 do i=1,3
    if (iuse_disc(i)) itest=itest+1
 enddo

! Sanity checks on equation of state for different requests
!
! If no disc then fail.
 if (itest == 0) then
    call fatal('setup_binarydisc','no discs requested...')
! If more than one disc and q_index is not zero than fail
 elseif (itest > 1) then
    minq = abs(minval(q_index(:)))
    if (size(vxyzu) < 4 .and. minq > epsilon(minq)) then
       call fatal('setup_binarydisc','locally isothermal eos for more than one disc requested, no ieos to handle this')
    endif
! If only one disc then set isink for which sink has the disc (0 if circumbinary)
 elseif (itest == 1) then
    do i=1,3
       if (iuse_disc(i)) isink = i-1 ! 1=binary,2=primary,3=secondary...but 1=primary,2=secondary
    enddo
 else
    call fatal('setup_binarydisc','this should never happen...')
 endif

!
!--setup circumbinary accretion disc
!
 totmass = 0.
 do i=1,3
   if (iuse_disc(i)) totmass = totmass + disc_mass(i)
 enddo
 print*,' total mass of system = ',totmass

 ntot = 0
 do i=1,3
    if (iuse_disc(i)) then
       print "(/,a)",'>>> Setting up circum'//trim(disctype(i))//' disc <<<'
       select case(i)
       case(3)
          xorigini(:) = xyzmh_ptmass(1:3,2)
          vorigini(:) = vxyz_ptmass(1:3,2)
          starmass    = m2
          Rochesizei  = Rochelobe
       case(2)
          xorigini(:) = xyzmh_ptmass(1:3,1)
          vorigini(:) = vxyz_ptmass(1:3,1)
          starmass    = m1
          Rochesizei  = binary_a - Rochelobe
       case default
          !
          !--centre of mass of binary defined to be zero
          !  if set_binary is used
          !
          xorigini(:) = 0.
          vorigini(:) = 0.
          starmass    = m1 + m2
          Rochesizei  = huge(0.)
       end select
       if (R_out(i) >  Rochesizei .and. R_out(i) < binary_a) then
          print "(/,a,/)",'*** WARNING: Outer disc radius for circum'//trim(disctype(i))//&
                           ' > Roche lobe of '//trim(disctype(i))//' ***'
       endif

       npindisc = int(disc_mass(i)/totmass*np)

       call set_disc(id,master=master,&
                  npart   = npindisc, &
                  npart_start = ntot + 1, &
                  rmin    = R_in(i),   &
                  rmax    = R_out(i),  &
                  p_index = p_index(i), &
                  q_index = q_index(i), &
                  HoverR  = HoverR(i),  &
                  disc_mass = disc_mass(i),  &
                  star_mass = starmass, &
                  gamma   = gamma,  &
                  particle_mass = massoftype(1), &
                  xyz_origin = xorigini(:), &
                  vxyz_origin = vorigini(:), &
                  hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,polyk=polyk, &
                  inclination=xinc(i), &
                  ismooth=ismooth_edge(i), &
                  alpha=alpha, &
                  isink=isink, &
                  prefix = fileprefix)

! Store requested alpha_AV, and return alpha to original requested alpha_SS for using multiple discs
! if only one disc required, this is corrected below
       alpha_returned(i) = alpha
       alpha = alphaSS

       ntot = ntot + npindisc
    endif
 enddo
 npartoftype(igas) = ntot
 npart = ntot

! set up dust
 if (use_dustfrac) then
    do i=1,npart
       call set_dustfrac(dust_to_gas,dustfrac(i))
    enddo
 endif

! Sort out alpha
 number_of_discs = 0
 do i=1,3
    if (iuse_disc(i)) then
       number_of_discs = number_of_discs + 1
    endif
 enddo
 if (number_of_discs == 1) then
    do i=1,3
       if (iuse_disc(i)) alpha = alpha_returned(i)
    enddo
 else
    call warning('setup_binarydisc','multiple discs: cannot use alpha for alpha_SS, setting equal to 1.0 instead')
    alpha = 1.0
 endif
 !
 ! reset centre of mass to the origin
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 !
 !--set default options for the input file
 !
 iexternalforce = 0
 if (maxvxyzu==3 .and. iuse_disc(1) .and. q_index(1) > 0.) then
    ieos = 3  ! use locally isothermal equation of state, but ONLY for circumbinary disc
 endif
 if (deltat > 0.)  dtmax = deltat*period
 if (norbits >= 0) tmax  = norbits*period

! Sanity check on ieos = 6
 if (ieos==6 .and. isink==0) call fatal('setup_binarydisc','something''s gone wrong with ieos & isink...')

 return
end subroutine setpart

subroutine write_binaryinputfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer :: i
 logical :: donealpha

 donealpha = .false.

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binarydisc setup routines'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','number of particles',iunit)
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 write(iunit,"(/,a)") '# options for binary'
 call write_inopt(m1,'m1','primary mass',iunit)
 call write_inopt(m2,'m2','secondary mass',iunit)
 call write_inopt(ecc,'ecc','binary eccentricity',iunit)
 call write_inopt(binary_a,'binary_a','binary semi-major axis',iunit)
 call write_inopt(binary_inc,'binary_inc','i, inclination (deg)',iunit)
 call write_inopt(binary_posang,'binary_posang','Omega, PA of ascending node (deg)',iunit)
 call write_inopt(binary_argperi,'binary_argperi','w, argument of periapsis (deg)',iunit)

 call write_inopt(accr1,'accr1','primary accretion radius',iunit)
 call write_inopt(accr2,'accr2','secondary accretion radius',iunit)
 do i=1,3
    write(iunit,"(/,a)") '# options for circum'//trim(disctype(i))//' disc'
    call write_inopt(iuse_disc(i),'use_'//trim(disctype(i))//'disc','setup circum'//trim(disctype(i))//' disc',iunit)
    if (iuse_disc(i)) then
       call write_inopt(R_in(i),'R_in'//trim(disctype(i)),'inner radius for circum'//trim(disctype(i))//' disc',iunit)
       call write_inopt(R_out(i),'R_out'//trim(disctype(i)),'outer radius for circum'//trim(disctype(i))//' disc',iunit)
       call write_inopt(HoverR(i),'HoverR'//trim(disctype(i)),'H/R for circum'//trim(disctype(i))//' disc',iunit)
       call write_inopt(xinc(i),'xinc'//trim(disctype(i)),'inclination in degrees for circum'//trim(disctype(i))//' disc',iunit)
       call write_inopt(ismooth_edge(i),'ismooth'//trim(disctype(i)),&
            'smooth density at inner edge for circum'//trim(disctype(i))//' disc',iunit)
       call write_inopt(disc_mass(i),'discmass'//trim(disctype(i)),'disc mass for circum'//trim(disctype(i))//' disc',iunit)
       call write_inopt(p_index(i),'p_index'//trim(disctype(i)),&
                        'surface density profile for circum'//trim(disctype(i))//' disc',iunit)
       call write_inopt(q_index(i),'q_index'//trim(disctype(i)),&
                        'temperature profile for circum'//trim(disctype(i))//' disc',iunit)
       if (.not.donealpha) then
          call write_inopt(alphaSS,'alphaSS','desired Shakura-Sunyaev alpha viscosity parameter',iunit)
          donealpha = .true.
       endif
    endif
 enddo
 write(iunit,"(/,a)") '# timestepping'
 call write_inopt(norbits,'norbits','maximum number of binary orbits',iunit)
 call write_inopt(deltat,'deltat','output interval as fraction of binary orbital period',iunit)

 if (use_dustfrac) then
    write(iunit,"(/,a)") '# dust'
    call write_inopt(dust_to_gas,'dust_to_gas','initial dust-to-gas ratio',iunit)
 endif

 close(iunit)

end subroutine write_binaryinputfile

subroutine read_binaryinputfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 use units,        only:select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: i,nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(np,'np',db,min=0,errcount=nerr)
 call read_inopt(mass_unit,'mass_unit',db,errcount=nerr)
 call read_inopt(dist_unit,'dist_unit',db,errcount=nerr)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,ierr)
 if (ierr /= 0) then
    call error('setup_binarydisc','mass unit not recognised')
    nerr = nerr + 1
 endif

 call select_unit(dist_unit,udist,ierr)
 if (ierr /= 0) then
    call error('setup_binarydisc','length unit not recognised')
    nerr = nerr + 1
 endif

 call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
 call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
 call read_inopt(ecc,'ecc',db,min=0.,errcount=nerr)
 call read_inopt(binary_a,'binary_a',db,errcount=nerr)
 call read_inopt(binary_inc,'binary_inc',db,errcount=nerr)
 call read_inopt(binary_argperi,'binary_argperi',db,errcount=nerr)
 call read_inopt(binary_posang,'binary_posang',db,errcount=nerr)
 call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
 call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)
 do i=1,3
    call read_inopt(iuse_disc(i),'use_'//trim(disctype(i))//'disc',db,errcount=nerr)
    if (iuse_disc(i)) then
       call read_inopt(R_in(i),'R_in'//trim(disctype(i)),db,errcount=nerr)
       call read_inopt(R_out(i),'R_out'//trim(disctype(i)),db,errcount=nerr)
       call read_inopt(HoverR(i),'HoverR'//trim(disctype(i)),db,errcount=nerr)
       call read_inopt(xinc(i),'xinc'//trim(disctype(i)),db,errcount=nerr)
       call read_inopt(ismooth_edge(i),'ismooth'//trim(disctype(i)),db,errcount=nerr)
       call read_inopt(disc_mass(i),'discmass'//trim(disctype(i)),db,errcount=nerr)
       call read_inopt(p_index(i),'p_index'//trim(disctype(i)),db,errcount=nerr)
       call read_inopt(q_index(i),'q_index'//trim(disctype(i)),db,errcount=nerr)
    endif
 enddo
 if (any(iuse_disc(:))) call read_inopt(alphaSS,'alphaSS',db,errcount=nerr)
 call read_inopt(norbits,'norbits',db,errcount=nerr)
 call read_inopt(deltat,'deltat',db,errcount=nerr)
 if (use_dustfrac) call read_inopt(dust_to_gas,'dust_to_gas',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_binaryinputfile

end module setup
