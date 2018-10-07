!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  this module sets up the particles
!  (2014 Alice Cerioli)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    HoverRinput  -- H/R at R_in
!    R_in         -- inner radius
!    R_out        -- outer radius
!    a0           -- initial binary separation
!    accradius1   -- primary accretion radius
!    accradius2   -- secondary accretion radius
!    alphaSS      -- desired alpha_SS
!    discm        -- disc mass
!    inc          -- inclination (tilt) in degrees
!    massr        -- mass ratio
!    np           -- number of particles
!    p_indexinput -- surface density profile
!    q_indexinput -- temperature profile
!
!  DEPENDENCIES: extern_binary, externalforces, infile_utils, io, options,
!    physcon, prompting, setdisc, units
!+
!--------------------------------------------------------------------------
module setup
 use extern_binary, only:accradius1,accradius2,massr,a0 !,binary_posvel
 implicit none
 public :: setpart

 integer :: np
 real :: R_in, R_out, HoverRinput, discm, alphaSS
 real :: p_indexinput, q_indexinput, inc

 private

contains

!----------------------------------------------------------------
!
! This subroutine sets up an accretion disc around the center of mass / primary
!  black hole in a binary system that decays due to
!  gravitational waves emission.
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,        only:set_disc
 use units,          only:set_units
 use physcon,        only:pi,solarm
 use io,             only:master
 use options,        only:iexternalforce, alpha
 use externalforces, only:iext_binary
 integer,            intent(in)            :: id
 integer,            intent(out)           :: npart
 integer,            intent(out)           :: npartoftype(:)
 real,               intent(out)           :: xyzh(:,:)
 real,               intent(out)           :: polyk,gamma,hfact
 real,               intent(out)           :: vxyzu(:,:)
 real,               intent(out)           :: massoftype(:)
 real,               intent(inout)         :: time
 character (len=20), intent (in), optional :: fileprefix
 real     :: xinc
 integer  :: ierr
 !real :: xbinary(10),vbinary(6)

 logical :: iexist
 character(len=100) :: filename
 !
 !--set code units
 !
 call set_units(mass=1.*solarm,c=1.)

 filename=trim(fileprefix)//'.setup'

 np = size(xyzh(1,:))
 npart = np
 npartoftype(1) = npart
 gamma = 1.0
 hfact = 1.2
 time  = 0.
 accradius1 = 2.0 ! R_schwarzschild
 a0= 11.5
 R_in  = 6. ! R_isco =3*R_schw
 R_out  = 0.9*a0
 HoverRinput = 0.01
 discm  = 0.1 !1.e-10
 alphaSS = 0.1
 p_indexinput = 1.5
 q_indexinput = 0.75
 inc = 0.

 print "(a,/)",'Phantomsetup: routine to set a shrinking binary. '
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    call setup_interactive()
    call write_setupfile(filename)
    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
    stop
 else
    stop
 endif

 npart = np
 npartoftype(1) = npart
 massoftype(1) = discm/real(npart)
 xinc = inc*(pi/180.0) ! Must be in radians

 alpha=alphaSS

 call set_disc(id,master=master,&
               npart   = np,&
               rmin    = R_in,   &
               rmax    = R_out,  &
               p_index = p_indexinput,    &
               q_index = q_indexinput,   &
               HoverR  = HoverRinput,  &
               disc_mass = discm, &
               star_mass = 1.0,  &
               gamma     = gamma,  &
               particle_mass = massoftype(1), &
               hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,polyk=polyk,alpha=alpha, &
               inclination = xinc, &
               prefix = fileprefix )
 !
 !--set default options for the input file
 !
 iexternalforce = iext_binary

 !--------------------------------------------------
 ! If you want to translate the disc so it is around the primary uncomment the following lines
 !--------------------------------------------------
! call binary_posvel(time,xbinary,vbinary)
! do i=1,npart
!   xyzh(1,i) = xyzh(1,i) + xbinary(1)
!   xyzh(2,i) = xyzh(2,i) + xbinary(2)
!   xyzh(3,i) = xyzh(3,i) + xbinary(3)
!   vxyzu(1,i) = vxyzu(1,i) + vbinary(1)
!   vxyzu(2,i) = vxyzu(2,i) + vbinary(2)
!   vxyzu(3,i) = vxyzu(3,i) + vbinary(3)
! enddo

 return
end subroutine setpart

subroutine setup_interactive()
 use prompting, only:prompt
 use physcon,   only:solarm
 use units,     only:umass
    !
    !--set default options
    !
!    print "(a,f6.3,a,1pe8.2,a)",' Schwarzschild radius is ',2.0*udist/au,' AU for ',umass/solarm,' M_sun black hole'

    call prompt('Enter total number of gas particles ',np,0)
    call prompt('Enter mass ratio of binary',massr,0.,1.)

    accradius2 = massr*accradius1

    call prompt('Enter initial binary separation',a0,0.)
    call prompt('Enter accretion radius of the PRIMARY black hole ',accradius1,accradius1,a0)
    call prompt('Enter accretion radius of the SECONDARY black hole ',accradius2,accradius2,a0)

    call prompt('Enter inner disc edge R_in ',R_in,accradius1)
    call prompt('Enter outer disc edge R_out ',R_out,R_in)
    call prompt('Enter H/R at R=R_in ',HoverRinput,0.)
    call prompt('Enter p index of surface density profile Sigma = Sigma0*R^-p',p_indexinput,0.)
    call prompt('Enter q index of temperature profile cs = cs0*R^-q',q_indexinput,0.)
    print "(a,es10.4,a)",'Enter disc mass in units of ',umass/solarm,' solar masses (Mjup = 8.6 x 10^-4 Msun) '
    call prompt(' ',discm,0.)
    call prompt('Enter desired value of alpha_SS',alphaSS,0.)

    call prompt('Enter disc inclination (deg) ',inc,-180.,180.)

end subroutine setup_interactive

subroutine write_setupfile(filename)
 use infile_utils,  only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)")   '# input file for gwdisc setup routines'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','number of particles',iunit)

 write(iunit,"(/,a)") '# options for binary'
 call write_inopt(massr,'massr','mass ratio',iunit)
 call write_inopt(a0,'a0','initial binary separation',iunit)
 call write_inopt(accradius1,'accradius1','primary accretion radius',iunit)
 call write_inopt(accradius2,'accradius2','secondary accretion radius',iunit)

 write(iunit,"(/,a)") '# options for accretion disc'
 call write_inopt(R_in,'R_in','inner radius',iunit)
 call write_inopt(R_out,'R_out', 'outer radius',iunit)
 call write_inopt(HoverRinput,'HoverRinput','H/R at R_in',iunit)
 call write_inopt(discm,'discm','disc mass',iunit)
 call write_inopt(p_indexinput,'p_indexinput','surface density profile',iunit)
 call write_inopt(q_indexinput,'q_indexinput','temperature profile',iunit)
 call write_inopt(alphaSS,'alphaSS','desired alpha_SS',iunit)
 call write_inopt(inc,'inc','inclination (tilt) in degrees',iunit)

 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), dimension(:), allocatable :: db

 print "(a)",'reading setup options from '//trim(filename)

 call open_db_from_file(db,filename,iunit,ierr)
 nerr = 0
 call read_inopt(np,'np',db,min=1,errcount=nerr)
 call read_inopt(massr,'massr',db,min=0.,errcount=nerr)
 call read_inopt(a0,'a0',db,min=0.,errcount=nerr)
 call read_inopt(accradius1,'accradius1',db,min=0.,errcount=nerr)
 call read_inopt(accradius2,'accradius2',db,min=0.,errcount=nerr)
 call read_inopt(R_in,'R_in',db,min=0.,errcount=nerr)
 call read_inopt(R_out,'R_out',db,min=R_in,errcount=nerr)
 call read_inopt(HoverRinput,'HoverRinput',db,min=0.,errcount=nerr)
 call read_inopt(discm,'discm',db,min=0.,errcount=nerr)
 call read_inopt(p_indexinput,'p_indexinput',db,errcount=nerr)
 call read_inopt(q_indexinput,'q_indexinput',db,errcount=nerr)
 call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
 call read_inopt(inc,'inc',db,errcount=nerr)
 ierr = nerr
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
 endif
 call close_db(db)
 close(iunit)

end subroutine read_setupfile

end module setup
