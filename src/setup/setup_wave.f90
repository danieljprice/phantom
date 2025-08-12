!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of a linear sound wave in a box
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - K_drag      : *constant drag coefficient*
!   - ampl        : *perturbation amplitude*
!   - cs          : *sound speed in code units (sets polyk)*
!   - dtg         : *dust to gas ratio*
!   - dust_method : *dust method (1=one fluid,2=two fluid)*
!   - npartx      : *number of gas particles in x direction*
!   - rhozero     : *initial gas density*
!
! :Dependencies: boundary, dim, dust, infile_utils, io, kernel, mpidomain,
!   options, part, physcon, prompting, set_dust, setup_params, unifdis
!
 use dim,          only:use_dust
 use options,      only:use_dustfrac
 use setup_params, only:rhozero
 use dust,         only:K_code
 use part,         only:maxp
 implicit none

 public :: setpart

 private

 ! Module variables for setup parameters
 integer :: npartx,dust_method
 real    :: ampl,kwave,xmin_wave
 real    :: cs,dtg

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis,rho_func
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:set_particle_type,igas,idust,dustfrac,periodic,ndustsmall,ndustlarge,ndusttypes
 use physcon,      only:pi
 use kernel,       only:radkern
 use dim,          only:maxvxyzu
 use dust,         only:idrag
 use set_dust,     only:set_dustfrac
 use mpidomain,    only:i_belong
 use infile_utils, only:get_options
 use kernel,       only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: totmass,fac,deltax,deltay,deltaz
 integer :: i,ierr
 integer :: itype,itypes,ntypes
 integer :: npart_previous
 logical, parameter :: ishift_box =.true.
 real, parameter    :: dust_shift = 0.
 real    :: xmin_dust,xmax_dust,ymin_dust,ymax_dust,zmin_dust,zmax_dust
 real    :: denom,length,uuzero,przero,massfac
 real    :: xmini,xmaxi
 procedure(rho_func), pointer :: density_func
!
! default options
!
 npartx  = 64
 ntypes  = 1
 rhozero = 1.
 massfac = 1.
 cs      = 1.
 ampl    = 1.d-4
 use_dustfrac = .false.
 ndustsmall = 0
 ndustlarge = 0
 dust_method = 2
 dtg = 1.
 idrag = 2
 K_code = 1000.
 hfact = hfact_default
 ! Print setup header
 if (id==master) print "(/,a,/)",'  >>> Setting up particles for linear wave test <<<'

 ! Get setup parameters from file or interactive setup
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 ! Set dust parameters based on method
 if (use_dust) then
    if (dust_method == 1) then
       use_dustfrac = .true.
       massfac = 1. + dtg
       ndustsmall = 1
    else
       use_dustfrac = .false.
       ntypes  = 2
       ndustlarge = 1
    endif
 endif

 ndusttypes = ndustsmall + ndustlarge
!
! boundaries
!
 xmini = -0.5
 xmaxi = 0.5
 length = xmaxi - xmini
 deltax = length/npartx
 ! try to give y boundary that is a multiple of 6 particle spacings in the low density part
 fac = 6.*(int((1.-epsilon(0.))*radkern/6.) + 1)
 deltay = fac*deltax*sqrt(0.75)
 deltaz = fac*deltax*sqrt(6.)/3.
 call set_boundary(xmini,xmaxi,-deltay,deltay,-deltaz,deltaz)
 xmin_wave = xmin
!
! general parameters
!
 time   = 0.
 kwave  = 2.*pi/length
 denom = length - ampl/kwave*(cos(kwave*length)-1.0)
!
! setup particles in the closepacked lattice
!
 if (maxvxyzu >= 4) then
    gamma = 5./3.
 else
    gamma  = 1.
 endif

 if (maxvxyzu < 4) then
    polyk = cs**2
    if (id==master) print*,' polyk = ',polyk
 else
    polyk = 0.
 endif

 npart = 0
 npart_total = 0
 npartoftype(:) = 0
 density_func => dens_func  ! desired density function

 overtypes: do itypes=1,ntypes
    select case (itypes)
    case(1)
       itype = igas
    case(2)
       itype = idust
       rhozero = dtg*rhozero
    end select

    npart_previous = npart

    if (itype == igas) then
       call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                         hfact,npart,xyzh,periodic,nptot=npart_total,rhofunc=density_func,mask=i_belong)
       xmin_dust = xmin + dust_shift*deltax
       xmax_dust = xmax + dust_shift*deltax
       ymin_dust = ymin + dust_shift*deltax
       ymax_dust = ymax + dust_shift*deltax
       zmin_dust = zmin + dust_shift*deltax
       zmax_dust = zmax + dust_shift*deltax
    else
       call set_unifdis('closepacked',id,master,xmin_dust,xmax_dust,ymin_dust, &
                         ymax_dust,zmin_dust,zmax_dust,deltax, &
                         hfact,npart,xyzh,periodic,nptot=npart_total,rhofunc=density_func,mask=i_belong)
    endif

    !--set which type of particle it is
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)

       vxyzu(1,i) = ampl*sin(kwave*(xyzh(1,i)-xmin))
       vxyzu(2:3,i) = 0.
!
!--perturb internal energy if not using a polytropic equation of state
!  (do this before density is perturbed)
!
       if (maxvxyzu >= 4) then
          if (gamma > 1.) then
             uuzero = cs**2/(gamma*(gamma-1.))
             przero = (gamma-1.)*rhozero*uuzero
          else
             uuzero = 3./2.*cs**2
             przero = cs**2*rhozero
          endif
          vxyzu(4,i) = uuzero + przero/rhozero*ampl*sin(kwave*(xyzh(1,i)-xmin))
       endif
!
!--one fluid dust: set dust fraction on gas particles
!
       if (use_dustfrac) then
          if (itype==igas) then
             call set_dustfrac(dtg,dustfrac(:,i))
          else
             dustfrac(:,i) = 0.
          endif
       endif
    enddo

    npartoftype(itype) = npart - npart_previous
    if (id==master) print*,' npart = ',npart,npart_total

    totmass = massfac*rhozero*dxbound*dybound*dzbound
    if (id==master) print*,' box volume = ',dxbound*dybound*dzbound,' rhozero = ',rhozero

    massoftype(itype) = totmass/npartoftype(itype)
    if (id==master) print*,' particle mass = ',massoftype(itype)

 enddo overtypes

end subroutine setpart

!----------------------------------------------------
!+
!  callback function giving desired density profile
!+
!----------------------------------------------------
real function dens_func(x)
 real, intent(in) :: x

 dens_func = 1. + ampl*sin(kwave*(x - xmin_wave))

end function dens_func

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for linear wave setup'
 call write_inopt(npartx,'npartx','number of gas particles in x direction',iunit)
 call write_inopt(rhozero,'rhozero','initial gas density',iunit)
 call write_inopt(cs,'cs','sound speed in code units (sets polyk)',iunit)
 call write_inopt(ampl,'ampl','perturbation amplitude',iunit)
 if (use_dust) then
    call write_inopt(dust_method,'dust_method','dust method (1=one fluid,2=two fluid)',iunit)
    call write_inopt(dtg,'dtg','dust to gas ratio',iunit)
    call write_inopt(K_code(1),'K_drag','constant drag coefficient',iunit)
 endif
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
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartx,'npartx',db,min=8,errcount=nerr)
 call read_inopt(rhozero,'rhozero',db,min=0.,errcount=nerr)
 call read_inopt(cs,'cs',db,min=0.,errcount=nerr)
 call read_inopt(ampl,'ampl',db,errcount=nerr)
 if (use_dust) then
    call read_inopt(dust_method,'dust_method',db,min=1,max=2,errcount=nerr)
    call read_inopt(dtg,'dtg',db,min=0.,errcount=nerr)
    call read_inopt(K_code(1),'K_drag',db,min=0.,errcount=nerr)
 endif
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Interactive setup
!+
!-----------------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt

 print "(/,a,/)",'  >>> Setting up particles for linear wave test <<<'
 call prompt(' enter number of gas particles in x ',npartx,8,int(maxp/144.))
 if (use_dust) then
    call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
    if (dust_method == 1) then
       use_dustfrac = .true.
       K_code = 1000. ! sensible default option for one fluid
    else
       use_dustfrac = .false.
    endif
    call prompt('Enter dust to gas ratio',dtg,0.)
    call prompt('Enter constant drag coefficient',K_code(1),0.)
 endif

 call prompt('enter gas density (gives particle mass)',rhozero,0.)
 call prompt('enter sound speed in code units (sets polyk)',cs,0.)
 call prompt('enter perturbation amplitude',ampl)

end subroutine setup_interactive

end module setup
