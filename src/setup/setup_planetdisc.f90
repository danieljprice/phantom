!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Set up for the de Val Borro et al. planet-disc comparison problem
! this is a simplified disc setup with a planet on a prescribed
! orbit using a time-dependent potential
!
! :References: de Val Borro et al. (2006), MNRAS 370, 529
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - HoverR     : *H/R at R_in*
!   - R_in       : *inner radius*
!   - R_out      : *outer radius*
!   - accradius1 : *primary accretion radius*
!   - accradius2 : *secondary accretion radius*
!   - alphaSS    : *desired alpha_SS*
!   - m2         : *m2*
!   - norbits    : *number of orbits*
!   - np         : *number of particles*
!   - p_index    : *surface density profile*
!   - q_index    : *temperature profile*
!   - sig0       : *disc surface density*
!
! :Dependencies: extern_binary, externalforces, infile_utils, io, options,
!   physcon, prompting, setdisc, timestep, units
!
 implicit none
 public :: setpart

 integer :: np, norbits
 real :: R_in, R_out, HoverR, sig0, sig_in, alphaSS
 real :: p_index, q_index

 private

contains

!----------------------------------------------------------------
!
! This subroutine sets up planet-disc interaction problem
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,       only:set_disc
 use units,         only:set_units
 use physcon,       only:solarm,au,pi
 use io,            only:master
 use options,       only:iexternalforce,alpha
 use timestep,      only:dtmax,tmax
 use prompting,     only:prompt
 use extern_binary, only:accradius1,accradius2 !,binary_posvel
 use extern_binary,  only:mass2,eps_soft1,eps_soft2,ramp
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
 integer :: ierr

 !real :: xbinary(10),vbinary(6)
 real :: a0

 logical :: iexist
 character(len=100) :: filename

 !
 !--set code units
 !
 call set_units(dist=au,mass=solarm,G=1.d0)

 filename=trim(fileprefix)//'.setup'

 np = size(xyzh(1,:))
 npart = np
 npartoftype(1) = npart
 gamma = 1.0
 hfact = 1.2
 time  = 0.
 a0    = 1.
 mass2 = 1.e-3
 HoverR = 0.05
 accradius1 = 0.0
 accradius2 = 0.3
 eps_soft1 = 0.6*HoverR*a0
 R_in  = 0.4
 R_out = 2.5
 sig0  = 0.002/(pi*a0**2)
 alphaSS = 0.01
 p_index = 1.
 q_index = 0.25
 ramp = .true.
 norbits = 100

 print "(a,/)",'Phantomsetup: routine to setup planet-disc interaction with fixed planet orbit '
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       print "(a)",' ERRORS reading '//trim(filename)//', rewriting...'
       call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    !
    !--set default options
    !
    call prompt('Enter total number of gas particles ',np,0,size(xyzh(1,:)))
    call prompt('Enter mplanet',mass2,0.,1.)

    call prompt('Enter accretion radius of the PRIMARY (star)',accradius1,accradius1,1.)
    call prompt('Enter accretion radius of the SECONDARY (planet)',accradius2,accradius2,1.)

    call prompt('Enter softening radius of the PRIMARY (star)',eps_soft1,0.,1.)
    call prompt('Enter softening radius of the SECONDARY (planet)',eps_soft2,0.,1.)

    call prompt('Enter inner disc edge R_in ',R_in,accradius1)
    call prompt('Enter outer disc edge R_out ',R_out,R_in)
    call prompt('Enter H/R at R=R_in ',HoverR,0.)
    call prompt('Enter p index of surface density profile Sigma = Sigma0*R^-p',p_index,0.)
    call prompt('Enter q index of temperature profile cs = cs0*R^-q',q_index,0.)
    call prompt('Enter Sigma0 for disc',sig0)
    call prompt('Enter desired value of alpha_SS',alphaSS,0.)
    !
    !--write default input file
    !
    call write_setupfile(filename)

    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
    stop
 else
    stop
 endif

 npart = np
 npartoftype(1) = npart

 alpha=alphaSS

 !--set sig_in as required for set_disc (sig_norm at R=Rin)
 sig_in = sig0*R_in**p_index

 call set_disc(id,master     = master,        &
               npart         = np,            &
               rmin          = R_in,          &
               rmax          = R_out,         &
               p_index       = p_index,  &
               q_index       = q_index,  &
               HoverR        = HoverR,   &
               sig_norm      = sig_in,        &
               star_mass     = 1.0,           &
               gamma         = gamma,         &
               particle_mass = massoftype(1), &
               hfact         = hfact,         &
               xyzh          = xyzh,          &
               vxyzu         = vxyzu,         &
               polyk         = polyk,         &
               alpha         = alpha,         &
               prefix        = fileprefix)
 !
 !--set default options for the input file
 !
 iexternalforce = iext_binary

 dtmax = 2.*pi
 tmax = norbits*dtmax

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

end subroutine setpart


subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use extern_binary, only:accradius1,accradius2,binary_posvel
 use extern_binary, only:mass2
 implicit none
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for planetdisc setup routine'

 write(iunit,"(/,a)") '# resolution'

 call write_inopt(np,'np','number of particles',iunit)

 write(iunit,"(/,a)") '# options for binary'

 call write_inopt(mass2,'mplanet','m2',iunit)
 call write_inopt(accradius1,'accradius1','primary accretion radius',iunit)
 call write_inopt(accradius2,'accradius2','secondary accretion radius',iunit)
 call write_inopt(norbits, 'norbits', 'number of orbits', iunit)

 write(iunit,"(/,a)") '# options for accretion disc'

 call write_inopt(R_in,'R_in','inner radius',iunit)
 call write_inopt(R_out,'R_out', 'outer radius',iunit)
 call write_inopt(HoverR,'HoverR','H/R at R_in',iunit)
 call write_inopt(sig0,'sig0','disc surface density',iunit)
 call write_inopt(p_index,'p_index','surface density profile',iunit)
 call write_inopt(q_index,'q_index','temperature profile',iunit)
 call write_inopt(alphaSS,'alphaSS','desired alpha_SS',iunit)

 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use extern_binary, only:accradius1,accradius2,binary_posvel
 use extern_binary, only:mass2
 implicit none
 character(len=*), intent(in) :: filename
 integer, intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), dimension(:), allocatable :: db

 print "(a)",'reading setup options from '//trim(filename)

 call open_db_from_file(db,filename,iunit,ierr)

 nerr = 0
 call read_inopt(np,'np',db,errcount=nerr)
 call read_inopt(mass2,'mplanet',db,errcount=nerr)
 call read_inopt(accradius1,'accradius1',db,errcount=nerr)
 call read_inopt(accradius2,'accradius2',db,errcount=nerr)
 call read_inopt(R_in,'R_in',db,errcount=nerr)
 call read_inopt(R_out,'R_out',db,errcount=nerr)
 call read_inopt(HoverR,'HoverR',db,errcount=nerr)
 call read_inopt(sig0,'sig0',db,errcount=nerr)
 call read_inopt(p_index,'p_index',db,errcount=nerr)
 call read_inopt(q_index,'q_index',db,errcount=nerr)
 call read_inopt(alphaSS,'alphaSS',db,errcount=nerr)
 call read_inopt(norbits,'norbits',db,errcount=nerr,min=1)
 ierr = nerr

 call close_db(db)
 close(iunit)

end subroutine read_setupfile

end module setup
