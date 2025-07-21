!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of dust settling problem from PL15
!
! :References: Price & Laibe (2015), MNRAS 451, 5332
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - HonR              : *ratio of H/R*
!   - Rdisc             : *radius at which the calculations will be made [au]*
!   - Rmax              : *Complete N revolutions at what radius? [au]*
!   - dust_to_gas_ratio : *dust-to-gas ratio*
!   - graindenscgs      : *grain density [g/cm^3]*
!   - grainsizecgs      : *grain size in [cm]*
!   - ndusttypes        : *number of grain sizes*
!   - norbit            : *Number of orbits at Rmax*
!   - npartx            : *requested number of particles in x-direction*
!   - rhozero           : *midplane density (> 0 for code units; < 0 for cgs)*
!   - sindex            : *power-law index, e.g. MRN*
!   - smaxcgs           : *maximum grain size [cm]*
!   - smincgs           : *minimum grain size [cm]*
!   - stellar_mass      : *mass of the central star [Msun]*
!
! :Dependencies: boundary, dim, dust, externalforces, infile_utils, io,
!   mpidomain, options, part, physcon, prompting, radiation_utils,
!   set_dust, setup_params, table_utils, timestep, unifdis, units
!
 use part,           only:ndusttypes,ndustsmall
 use dust,           only:grainsizecgs,graindenscgs
 use setup_params,   only:rhozero
 use externalforces, only:mass1
 use dim,            only:use_dust,do_radiation
 implicit none
 public  :: setpart

 private
 integer :: npartx,norbit
 real    :: Rdisc_au,Rmax_au
 real    :: H0,HonR,dtg,smincgs,smaxcgs,sindex

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params,   only:npart_total
 use io,             only:master
 use unifdis,        only:set_unifdis,rho_func
 use boundary,       only:set_boundary,xmin,xmax,zmin,zmax,dxbound,dzbound
 use part,           only:labeltype,set_particle_type,igas,dustfrac,&
                          grainsize,graindens,periodic,rad
 use physcon,        only:pi,au,solarm
 use dim,            only:maxvxyzu,maxdustsmall
 use externalforces, only:Rdisc,iext_discgravity
 use options,        only:iexternalforce,use_dustfrac
 use timestep,       only:dtmax,tmax
 use units,          only:set_units,udist,umass,utime
 use dust,           only:init_drag,idrag,get_ts
 use set_dust,       only:set_dustfrac,set_dustbinfrac
 use table_utils,    only:logspace
 use mpidomain,      only:i_belong
 use radiation_utils,only:set_radiation_and_gas_temperature_equal
 use infile_utils,   only:get_options
 use kernel,         only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer            :: i,iregime,ierr
 integer            :: itype
 integer            :: npart_previous
 real               :: totmass,deltax,dz,length
 real               :: omega,ts(maxdustsmall),dustbinfrac(maxdustsmall)
 real               :: xmini,xmaxi,ymaxdisc,cs,t_orb,Rmax
 procedure(rho_func), pointer :: density_func
!
!--default options
!
 hfact        = hfact_default
 npartx       = 32
 rhozero      = 1.e-3
 dtg          = 0.01
 grainsize    = 0.
 graindens    = 0.
 grainsizecgs = 0.1
 graindenscgs = 3.
 ndustsmall   = 1
 smincgs      = 1.e-5
 smaxcgs      = 1.
 sindex       = 3.5
 HonR         = 0.05
 Rdisc        = 5.
 Rdisc_au     = 5.*10.
 Rmax_au      = Rdisc_au
 mass1        = 1.
 norbit       = 15
 time         = 0.
 itype        = igas
 iexternalforce = iext_discgravity
!
!--Read / write .setup file
!
 print "(/,a,/)",'  >>> Setting up particles for dust settling test <<<'
 call set_units(dist=10.*au,mass=solarm,G=1.)

 !--Read values from .setup
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 if (use_dust) then
    use_dustfrac = .true.
    ndusttypes = ndustsmall
    if (ndusttypes > 1) then
       call set_dustbinfrac(smincgs,smaxcgs,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
       grainsize(1:ndusttypes) = grainsize(1:ndusttypes)/udist
       graindens(1:ndusttypes) = graindenscgs/umass*udist**3
    else
       grainsize(1) = grainsizecgs/udist
       graindens(1) = graindenscgs/umass*udist**3
    endif
 endif
!
!--set remaining parameters
!
 if (rhozero < 0.) rhozero = -rhozero*udist**3/umass
 Rdisc = Rdisc_au * 0.1
 Rmax  = Rmax_au  * 0.1
 H0    = HonR*Rdisc
 omega = sqrt(mass1/Rdisc**3)
 t_orb = 2.*pi/sqrt(mass1/Rmax**3)
 cs    = H0*omega
 dtmax = 0.1*t_orb
 tmax  = norbit*t_orb
 if (maxvxyzu >= 4) then
    gamma = 5./3.
 else
    gamma = 1.
    polyk = cs**2
 endif
!
!--get stopping time information
!
 if (use_dust) then
    call init_drag(ierr)
    do i=1,ndustsmall
       call get_ts(idrag,i,grainsize(i),graindens(i),rhozero,0.0*rhozero,cs,0.,ts(i),iregime)
       print*,'s (cm) =',grainsize(i),'   ','St = ts * Omega =',ts(i)*omega
    enddo
 endif
!
!--boundaries
!
 xmini  = -0.25
 xmaxi  =  0.25
 length = xmaxi - xmini
 deltax = length/npartx
 dz     = 2.*sqrt(6.)/npartx
 !deltay = fac*deltax*sqrt(0.75)
 call set_boundary(xmini,xmaxi,-10.*H0,10.*H0,-dz,dz)

 npart = 0
 npart_total = 0
 npartoftype(:) = 0

 !--get total mass from integration of density profile
 ymaxdisc = 3.*H0
 totmass  = 2.*rhozero*sqrt(0.5*pi)*H0*erf(ymaxdisc/(sqrt(2.)*H0))*dxbound*dzbound
 npart_previous = npart
 density_func => rhofunc
 call set_unifdis('closepacked',id,master,xmin,xmax,-ymaxdisc,ymaxdisc,zmin,zmax,deltax, &
                   hfact,npart,xyzh,periodic,nptot=npart_total,&
                   rhofunc=density_func,dir=2,mask=i_belong)

 !--set which type of particle it is
 do i=npart_previous+1,npart
    call set_particle_type(i,itype)
    vxyzu(:,i) = 0.
    !--set internal energy if necessary
    if (maxvxyzu >= 4) then
       if (gamma > 1.) then
          vxyzu(4,i) = cs**2/(gamma-1.)
       else
          vxyzu(4,i) = 1.5*cs**2
       endif
    endif

    !--one fluid dust: set dust fraction on gas particles
    if (use_dustfrac) then
       if (ndusttypes > 1) then
          dustfrac(1:ndusttypes,i) = dustbinfrac(1:ndusttypes)*dtg
       else
          call set_dustfrac(dtg,dustfrac(:,i))
       endif
    elseif (use_dust) then
       dustfrac(:,i) = 0.
    endif
 enddo

 npartoftype(itype) = npart - npart_previous
 massoftype(itype)  = totmass/npartoftype(itype)
 if (use_dust) massoftype(itype) = massoftype(itype)*(1. + dtg)

 if (do_radiation) call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad)

 if (id==master) then
    print*,' npart,npart_total     = ',npart,npart_total
    print*,' particle mass         = ',massoftype(itype),'code units'
    print*,' particle mass         = ',massoftype(itype)*umass,'g'
    print*,' total mass            = ',npart*massoftype(itype)*umass,'g'
    print*,' mid-plane gas density = ',rhozero*umass/udist**3,'g/cm^3'
    print*,' Rdisc                 = ',Rdisc*udist/au,'au'
    print*,' H0                    = ',H0*udist/au,'au'
    print*,' H/R                   = ',HonR
    print*,' cs                    = ',cs*udist/utime,'cm/s'
 endif
 if (HonR > 0.1) then
    print*, ' '
    print*, 'WARNING! This disc is hot, and *may* blow apart rather than settle.  A smaller ratio of H/R may be required.'
    print*, '         The default value of 0.05 yields stable results'
    print*, ' '
 endif

end subroutine setpart
!----------------------------------------------------------------
!+
!  Calculate the vertical density in the disc
!+
!----------------------------------------------------------------
real function rhofunc(x)
 real, intent(in) :: x

 rhofunc = exp(-0.5*(x/H0)**2)

end function rhofunc
!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for dust settling routine'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(npartx,'npartx','requested number of particles in x-direction',iunit)
 write(iunit,"(/,a)") '# Disc options'
 call write_inopt(rhozero,'rhozero','midplane density (> 0 for code units; < 0 for cgs)',iunit)
 call write_inopt(mass1,'stellar_mass','mass of the central star [Msun]',iunit)
 call write_inopt(Rdisc_au,'Rdisc','radius at which the calculations will be made [au]',iunit)
 call write_inopt(HonR,'HonR','ratio of H/R',iunit)
 call write_inopt(Rmax_au,'Rmax','Complete N revolutions at what radius? [au]',iunit)
 call write_inopt(norbit,'norbit','Number of orbits at Rmax',iunit)
 if (use_dust) then
    write(iunit,"(/,a)") '# Dust properties'
    call write_inopt(dtg,'dust_to_gas_ratio','dust-to-gas ratio',iunit)
    call write_inopt(ndusttypes,'ndusttypes','number of grain sizes',iunit)
    if (ndusttypes > 1) then
       call write_inopt(smincgs,'smincgs','minimum grain size [cm]',iunit)
       call write_inopt(smaxcgs,'smaxcgs','maximum grain size [cm]',iunit)
       call write_inopt(sindex, 'sindex', 'power-law index, e.g. MRN',iunit)
       call write_inopt(graindenscgs,'graindenscgs','grain density [g/cm^3]',iunit)
    else
       call write_inopt(grainsizecgs,'grainsizecgs','grain size in [cm]',iunit)
       call write_inopt(graindenscgs,'graindenscgs','grain density [g/cm^3]',iunit)
    endif
 endif
 close(iunit)

end subroutine write_setupfile
!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)

 !--Read values
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartx,'npartx',db,ierr)
 call read_inopt(rhozero,'rhozero',db,ierr)
 call read_inopt(mass1,'stellar_mass',db,ierr)
 call read_inopt(Rdisc_au,'Rdisc',db,ierr)
 call read_inopt(HonR,'HonR',db,ierr)
 call read_inopt(Rmax_au,'Rmax',db,ierr)
 call read_inopt(norbit,'norbit',db,ierr)
 if (use_dust) then
    call read_inopt(dtg,'dust_to_gas_ratio',db,ierr)
    call read_inopt(ndusttypes,'ndusttypes',db,ierr)
    if (ndusttypes > 1) then
       call read_inopt(smincgs,'smincgs',db,ierr)
       call read_inopt(smaxcgs,'smaxcgs',db,ierr)
       call read_inopt(sindex,'cs_sphere',db,ierr)
       call read_inopt(graindenscgs,'graindenscgs',db,ierr)
    else
       call read_inopt(grainsizecgs,'grainsizecgs',db,ierr)
       call read_inopt(graindenscgs,'graindenscgs',db,ierr)
    endif
 endif
 call close_db(db)

 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_dustsettle: ',ierr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

subroutine setup_interactive()
 use prompting, only:prompt
 use part,      only:maxp
 use dim,       only:maxdustsmall

 call prompt('Enter number of gas particles in x ',npartx,8,maxp/144)
 call prompt('Enter gas midplane density (> 0 for code units; < 0 for cgs)',rhozero)
 call prompt('Enter the mass of the central star [Msun]',mass1,0.)
 call prompt('Enter radius at which the calculations will be made, Rdisc [au]',Rdisc_au, 0.)
 call prompt('Enter the ratio of H/R at Rdisc',HonR, 0.)
 Rmax_au = Rdisc_au
 call prompt('Complete N revolutions at what radius, Rmax? [au]',Rmax_au, 0.)
 call prompt('How many orbits at Rmax would you like to complete?',norbit, 1)
 if (use_dust) then
    !--currently assume one fluid dust
    call prompt('Enter total dust to gas ratio',dtg,0.)
    call prompt('How many grain sizes do you want?',ndustsmall,1,maxdustsmall)
    if (ndustsmall > 1) then
       !--grainsizes
       call prompt('Enter minimum grain size in cm',smincgs,0.)
       call prompt('Enter maximum grain size in cm',smaxcgs,0.)
       !--mass distribution
       call prompt('Enter power-law index, e.g. MRN',sindex)
       !--grain density
       call prompt('Enter grain density in g/cm^3',graindenscgs,0.)
    else
       call prompt('Enter grain size in cm',grainsizecgs,0.)
       call prompt('Enter grain density in g/cm^3',graindenscgs,0.)
    endif
 endif

end subroutine setup_interactive

end module setup
