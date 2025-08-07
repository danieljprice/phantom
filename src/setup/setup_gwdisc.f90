!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! this module sets up the particles
!  (2014 Alice Cerioli)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - HoverRinput         : *H/R at R_in*
!   - R_in                : *inner radius*
!   - R_out               : *outer radius*
!   - a0                  : *initial binary separation*
!   - accradius1          : *primary accretion radius*
!   - accradius2          : *secondary accretion radius*
!   - alphaSS             : *desired alpha_SS*
!   - disc_around_primary : *place disc around primary?*
!   - discm               : *disc mass*
!   - inc                 : *inclination (tilt) in degrees*
!   - mass2               : *mass of secondary*
!   - np                  : *number of particles*
!   - p_indexinput        : *surface density profile*
!   - q_indexinput        : *temperature profile*
!
! :Dependencies: extern_binary, externalforces, infile_utils, io, kernel,
!   options, part, physcon, prompting, setdisc, units
!
 use extern_binary, only:accradius1,accradius2,mass1,mass2,a0 !,binary_posvel
 implicit none
 public :: setpart

 integer :: np
 real :: R_in, R_out, HoverRinput, discm, alphaSS
 real :: p_indexinput, q_indexinput, inc
 logical :: disc_around_primary = .false.

 private

contains

!--------------------------------------------------------------------------------
!+
!  This subroutine sets up an accretion disc around the centre of mass / primary
!  black hole in a binary system that decays due to gravitational wave emission.
!+
!--------------------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,        only:set_disc
 use units,          only:set_units
 use physcon,        only:pi,solarm
 use io,             only:master
 use options,        only:iexternalforce,alpha
 use externalforces, only:iext_binary
 use infile_utils,   only:get_options
 use kernel,         only:hfact_default
 use extern_binary,  only:binary_posvel
 use part,           only:igas
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
 integer  :: ierr,i
 real :: xbinary(10),vbinary(6)
 !
 !--set code units
 !
 call set_units(mass=1.*solarm,c=1.)

 time  = 0.
 hfact = hfact_default

 np = size(xyzh(1,:))
 gamma = 1.0
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
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 npart = np
 npartoftype(igas) = npart
 massoftype(igas) = discm/real(npart)
 xinc = inc*(pi/180.0) ! Must be in radians

 alpha=alphaSS

 call set_disc(id,master=master,npart=np,rmin=R_in,rmax=R_out,p_index=p_indexinput,&
               q_index=q_indexinput,HoverR=HoverRinput,disc_mass=discm,star_mass=1.0,&
               gamma=gamma,particle_mass=massoftype(igas),hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,&
               polyk=polyk,alpha=alpha,inclination=xinc,prefix=fileprefix)
 !
 !--set default options for the input file
 !
 iexternalforce = iext_binary

 if (disc_around_primary) then
    ! translate the disc so it is around the primary
    call binary_posvel(time,xbinary,vbinary)
    do i=1,npart
       xyzh(1,i) = xyzh(1,i) + xbinary(1)
       xyzh(2,i) = xyzh(2,i) + xbinary(2)
       xyzh(3,i) = xyzh(3,i) + xbinary(3)
       vxyzu(1,i) = vxyzu(1,i) + vbinary(1)
       vxyzu(2,i) = vxyzu(2,i) + vbinary(2)
       vxyzu(3,i) = vxyzu(3,i) + vbinary(3)
    enddo
 endif

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Interactive setup routine
!+
!-----------------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 use physcon,   only:solarm
 use units,     only:umass
 !
 !--set default options
 !
!    print "(a,f6.3,a,1pe8.2,a)",' Schwarzschild radius is ',2.0*udist/au,' AU for ',umass/solarm,' M_sun black hole'

 call prompt('Enter total number of gas particles ',np,0)
 call prompt('Enter mass of secondary',mass2,0.,1.)

 accradius2 = mass2/mass1*accradius1

 call prompt('Enter initial binary separation',a0,0.)
 call prompt('Enter accretion radius of the PRIMARY black hole ',accradius1,accradius1,a0)
 call prompt('Enter accretion radius of the SECONDARY black hole ',accradius2,accradius2,a0)

 call prompt('Enter inner disc edge R_in ',R_in,accradius1)
 call prompt('Enter outer disc edge R_out ',R_out,R_in)
 call prompt('Enter H/R at R=R_in ',HoverRinput,0.)
 call prompt('Enter p index of surface density profile Sigma = Sigma0*R^-p',p_indexinput,0.)
 call prompt('Enter q index of temperature profile cs = cs0*R^-q',q_indexinput,0.)
 print "(a,es12.4,a)",'Enter disc mass in units of ',umass/solarm,' solar masses (Mjup = 8.6 x 10^-4 Msun) '
 call prompt(' ',discm,0.)
 call prompt('Enter desired value of alpha_SS',alphaSS,0.)

 call prompt('Enter disc inclination (deg) ',inc,-180.,180.)

end subroutine setup_interactive

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
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
 call write_inopt(mass2,'mass2','mass of secondary',iunit)
 call write_inopt(a0,'a0','initial binary separation',iunit)
 call write_inopt(accradius1,'accradius1','primary accretion radius',iunit)
 call write_inopt(accradius2,'accradius2','secondary accretion radius',iunit)
 call write_inopt(disc_around_primary,'disc_around_primary','place disc around primary?',iunit)

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

!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
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
 call read_inopt(mass2,'mass2',db,min=0.,errcount=nerr)
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
 call read_inopt(disc_around_primary,'disc_around_primary',db,errcount=nerr)
 ierr = nerr
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
 endif
 call close_db(db)
 close(iunit)

end subroutine read_setupfile

end module setup
