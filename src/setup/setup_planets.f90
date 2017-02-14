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
!   This module sets up general single or binary accretion discs and adds
!   some planets to the disc.
!
!  REFERENCES: None
!
!  OWNER: Chris Nixon
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    accr1             -- primary accretion radius
!    accr2             -- secondary accretion radius
!    alphaSS           -- desired Shakura-Sunyaev alpha viscosity parameter
!    binary_separation -- binary separation
!    ecc               -- binary eccentricity
!    m1                -- primary mass
!    m2                -- secondary mass
!    np                -- number of particles
!    udist             -- distance unit in cm
!    umass             -- mass unit in g
!
!  DEPENDENCIES: eos, infile_utils, io, options, part, physcon, prompting,
!    setbinary, setdisc, setplanets, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart
 !--private module variables
 real :: m1,m2,ecc,binary_separation,accr1,accr2,alphaSS
 logical :: iuse_disc(3),ismooth_edge(3)
 real :: R_in(3),R_out(3),mrat(3),HoverR(3),disc_mass(3),p_index(3),q_index(3),xinc(3)
 integer :: np
 character(len=*), dimension(3), parameter :: disctype = &
  (/'binary   ',&
    'primary  ', &
    'secondary'/)
 private
 real(kind=8) :: udist,umass

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc, only:set_disc
 use units,   only:set_units
 use io,      only:master,fatal,error,warning
 use physcon, only:au,solarm
 use options,        only:iexternalforce,alpha,ieos
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,maxvxyzu,igas
 use setbinary,      only:set_binary,Rochelobe_estimate
 use prompting,      only:prompt
 use eos,            only:isink
 use setplanets,     only:set_planets
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
 real :: Rochelobe,totmass,starmass,Rochesizei,minq
 integer :: i,ntot,npindisc,j,itest
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
 udist = 100.*au
 umass = solarm

 print "(/,65('-'),2(/,a),/,65('-'),/)",&
   ' Welcome to the Ultimate Binary Disc Setup Routine^TM', &
   ' Brought to you by Daniel Price, Chris Nixon and Giuseppe Lodato'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_binaryinputfile(filename)
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    call prompt('Enter code units of mass in g ',umass,0.)
    call prompt('Enter code units of distance in cm ',udist,0.)
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
    binary_separation = 1.
    call prompt('Enter the binary separation',binary_separation,0.)
    ecc = 0.

    call prompt('Enter the eccentricity of the binary',ecc,0.,1.)
    accr1 = 0.25*binary_separation
    call prompt('Enter accretion radius for primary (can be adjusted later)',accr1,0.,binary_separation)
    accr2 = 0.25*binary_separation
    call prompt('Enter accretion radius for secondary (can be adjusted later)',accr2,0.,binary_separation)

    xinc(:) = 0.0
    p_index(:) = 1.5
    q_index(:) = 0.75
    iuse_disc = .false.
    ismooth_edge = .false.

    Rochelobe = Rochelobe_estimate(m1,m2,binary_separation)
!    R_in  = (/2.**(2./3.)*binary_separation,accr,accr/)
    R_in = (/binary_separation + Rochelobe + accr1,accr1,accr2/)
    R_out = (/10.*R_in(1),binary_separation - Rochelobe - accr1,Rochelobe - accr2/)
    disc_mass = (/1.e-2*m1,1.e-3*m1,1.e-4*m1/)
    HoverR(:) = 0.1

    do while(.not.(any(iuse_disc)))
       iuse_disc(1) = .true.
       do i=1,3
          call prompt('Do you want a circum'//trim(disctype(i))//' disc?',iuse_disc(i))
          if (iuse_disc(i)) then
             if (i==1) then
                print "(a,es10.3,a)",' Recommended minimum value of Rin = ',R_in(i),' (a + Rochelobe + acc. radius)'
                call prompt('Enter inner radius of circum'//trim(disctype(i))//' disc',R_in(i),binary_separation)
                call prompt('Enter outer radius of circum'//trim(disctype(i))//' disc',R_out(i),R_in(i))
             else
                print "(a,es10.3)",' Setting R_in for circum'//trim(disctype(i))//' to ',R_in(i)
                print "(a,es10.3)",' Setting R_out for circum'//trim(disctype(i))//' to ',R_out(i)
             endif
             call prompt('Enter mass of circum'//trim(disctype(i))//' disc',disc_mass(i),0.)
             call prompt('Enter inclination angle to the binary in degrees',xinc(i))
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
 !--units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)

 Rochelobe = Rochelobe_estimate(m1,m2,binary_separation)

! Convert inclination into radians
 xinc(:) = xinc(:)*3.1415926535/180.0

!
!--add a binary
!
 nptmass = 0

 print*, 'eccentricity = ',ecc

 call set_binary(m1,massratio=m2/m1,semimajoraxis=binary_separation,eccentricity=ecc,accretion_radius1=accr1,&
                 accretion_radius2=accr2,xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,&
                 nptmass=nptmass)

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
          Rochesizei  = binary_separation - Rochelobe
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
       if (R_out(i) >  Rochesizei .and. R_out(i) < binary_separation) then
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
                  inclination=xinc(i),ismooth=ismooth_edge(i),alpha=alpha,isink=isink, &
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
 !--set default options for the input file
 !
 iexternalforce = 0
 if (maxvxyzu==3 .and. iuse_disc(1) .and. q_index(1) > 0.) then
    ieos = 3  ! use locally isothermal equation of state, but ONLY for circumbinary disc
 endif

! Sanity check on ieos = 6
 if (ieos==6 .and. isink==0) call fatal('setup_binarydisc','something''s gone wrong with ieos & isink...')

 call set_planets(xyzmh_ptmass,vxyz_ptmass,nptmass)

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
 call write_inopt(udist,'udist','distance unit in cm',iunit)
 call write_inopt(umass,'umass','mass unit in g',iunit)
 write(iunit,"(/,a)") '# options for binary'
 call write_inopt(m1,'m1','primary mass',iunit)
 call write_inopt(m2,'m2','secondary mass',iunit)
 call write_inopt(ecc,'ecc','binary eccentricity',iunit)
 call write_inopt(binary_separation,'binary_separation','binary separation',iunit)
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
 close(iunit)

end subroutine write_binaryinputfile

subroutine read_binaryinputfile(filename)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 21
 integer :: ierr,i
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(np,'np',db,ierr)
 call read_inopt(udist,'udist',db,ierr)
 call read_inopt(umass,'umass',db,ierr)
 call read_inopt(m1,'m1',db,ierr)
 call read_inopt(m2,'m2',db,ierr)
 call read_inopt(ecc,'ecc',db,ierr)
 call read_inopt(binary_separation,'binary_separation',db,ierr)
 call read_inopt(accr1,'accr1',db,ierr)
 call read_inopt(accr2,'accr2',db,ierr)
 do i=1,3
    call read_inopt(iuse_disc(i),'use_'//trim(disctype(i))//'disc',db,ierr)
    if (iuse_disc(i)) then
       call read_inopt(R_in(i),'R_in'//trim(disctype(i)),db,ierr)
       call read_inopt(R_out(i),'R_out'//trim(disctype(i)),db,ierr)
       call read_inopt(HoverR(i),'HoverR'//trim(disctype(i)),db,ierr)
       call read_inopt(xinc(i),'xinc'//trim(disctype(i)),db,ierr)
       call read_inopt(ismooth_edge(i),'ismooth'//trim(disctype(i)),db,ierr)
       call read_inopt(disc_mass(i),'discmass'//trim(disctype(i)),db,ierr)
       call read_inopt(p_index(i),'p_index'//trim(disctype(i)),db,ierr)
       call read_inopt(q_index(i),'q_index'//trim(disctype(i)),db,ierr)
    endif
 enddo
 if (any(iuse_disc(:))) call read_inopt(alphaSS,'alphaSS',db,ierr)
 call close_db(db)

end subroutine read_binaryinputfile

end module setup
