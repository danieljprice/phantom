!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Setup a circumnuclear gas cloud around outflowing TDE
!
! :References: None
!
! :Owner: Fitz) Hu
!
! :Runtime parameters:
!   - ieos             : *equation of state used*
!   - ignore_radius    : *ignore tde particle inside this radius (-ve = ignore all for injection)*
!   - m_target         : *target mass in circumnuclear gas cloud (in Msun) (-ve = ignore and use rho0)*
!   - m_threshold      : *threshold in solving rho0 for m_target (in Msun)*
!   - mu               : *mean molecular density of the cloud*
!   - nbreak           : *number of broken power laws*
!   - nprof            : *number of data points in the cloud profile*
!   - profile_filename : *filename for the cloud profile*
!   - rad_max          : *outer radius of the circumnuclear gas cloud*
!   - rad_min          : *inner radius of the circumnuclear gas cloud*
!   - remove_overlap   : *remove outflow particles overlap with circum particles*
!   - rhof_n_1         : *power law index of the section*
!   - rhof_rho0        : *density at rad_min (in g/cm^3) (-ve = ignore and calc for m_target)*
!   - temperature      : *temperature of the gas cloud (-ve = read from file)*
!   - use_func         : *if use broken power law for density profile*
!
! :Dependencies: eos, infile_utils, io, kernel, mpidomain, part, physcon,
!   setup_params, spherical, stretchmap, timestep, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 public :: modify_dump
 private :: rho,rho_tab,get_temp_r,uerg,calc_rhobreak,calc_rho0,write_setupfile,read_setupfile

 private
 integer           :: ieos_in,nprof,nbreak,nbreak_old
 real              :: temperature,mu,ignore_radius,rad_max,rad_min
 character(len=50) :: profile_filename
 character(len=3)  :: interpolation
 real, allocatable :: rhof_n(:),rhof_rbreak(:),rhof_rhobreak(:)
 real, allocatable :: rhof_n_in(:),rhof_rbreak_in(:)
 real, allocatable :: rad_prof(:),dens_prof(:)
 real              :: rhof_rho0,m_target,m_threshold
 logical           :: use_func,use_func_old,remove_overlap

contains

!----------------------------------------------------------------
!
!  Sets up a circumnuclear gas cloud
!
!----------------------------------------------------------------
subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use physcon,      only:solarm,years,mass_proton_cgs,kb_on_mh,kboltz,radconst
 use setup_params, only:npart_total
 use part,         only:igas,set_particle_type,pxyzu,delete_particles_inside_radius, &
                        delete_particles_outside_sphere,kill_particle,shuffle_part, &
                        eos_vars,itemp,igamma,igasP
 use io,           only:fatal,master,id
 use units,        only:umass,udist,utime,set_units,unit_density
 use timestep,     only:dtmax,tmax
 use eos,          only:ieos,gmw
 use kernel,       only:hfact_default
 use stretchmap,   only:get_mass_r
 use spherical,    only:set_sphere
 use mpidomain,    only:i_belong
 integer,           intent(inout)   :: npart
 integer,           intent(inout)   :: npartoftype(:)
 real,              intent(inout)   :: xyzh(:,:)
 real,              intent(inout)   :: vxyzu(:,:)
 real,              intent(inout)   :: massoftype(:)
 integer                       :: i,ierr,iunit=12,iprof
 integer                       :: np_sphere,npart_old
 real                          :: totmass,delta,r,rhofr,presi
 character(len=120)            :: fileset,fileprefix='radio'
 logical                       :: read_temp,setexists
 real, allocatable             :: masstab(:),temp_prof(:)
 character(len=15), parameter  :: default_name = 'default_profile'
 real, dimension(7), parameter :: dens_prof_default = (/8.9e-21, 5.1e-21, 3.3e-21, 2.6e-21, &
                                                        6.6e-25, 3.4e-25, 8.1e-26/), &
                                  rad_prof_default = (/8.7e16, 1.2e17, 1.4e17, 2.0e17, &
                                                       4.0e17, 4.8e17, 7.1e17/) ! profile from Cendes+2021
 procedure(rho), pointer       :: rhof

 !--Check for existence of the .params files
 fileset=trim(fileprefix)//'.params'
 inquire(file=fileset,exist=setexists)

 !--Set default values
 temperature       = 10.           ! Temperature in Kelvin
 mu                = 1.            ! mean molecular weight
 ieos_in           = 2
 ignore_radius     = 1.e14          ! in cm
 use_func          = .true.
 use_func_old      = use_func
 remove_overlap    = .true.
 !--Power law default setups
 rad_max           = 7.1e16        ! in cm
 rad_min           = 8.7e15        ! in cm
 nbreak            = 1
 nbreak_old        = nbreak
 rhof_rho0         = 1.e4*mu*mass_proton_cgs
 allocate(rhof_n(nbreak),rhof_rbreak(nbreak))
 rhof_n            = -1.7
 rhof_rbreak       = rad_min
 m_target          = dot_product(npartoftype,massoftype)*umass/solarm
 m_threshold       = 1.e-3

 !--Profile default setups
 read_temp         = .false.
 profile_filename  = default_name
 nprof             = 7
 interpolation     = 'log'

 !--Read values from .setup
 if (setexists) call read_setupfile(fileset,ierr)
 if (.not. setexists .or. ierr /= 0) then
    !--Prompt to get inputs and write to file
    call write_setupfile(fileset)
    stop
 elseif (nbreak /= nbreak_old) then
    !--Rewrite setup file
    write(*,'(a)') ' [nbreak] changed. Rewriting setup file ...'
    deallocate(rhof_n,rhof_rbreak)
    allocate(rhof_n(nbreak),rhof_rbreak(nbreak))
    rhof_n = -1.7
    rhof_rbreak = rad_min
    call write_setupfile(fileset)
    stop
 elseif (use_func .neqv. use_func_old) then
    !--Rewrite setup fi.e
    write(*,'(a)') ' [use_func] changed. Rewriting setup file ...'
    call write_setupfile(fileset)
    stop
 endif

 !--allocate memory
 if (use_func) then
    rhof => rho
    deallocate(rhof_n,rhof_rbreak)
    allocate(rhof_n(nbreak),rhof_rbreak(nbreak),rhof_rhobreak(nbreak))
    rhof_n(:) = rhof_n_in(1:nbreak)
    rhof_rbreak(:) = rhof_rbreak_in(1:nbreak)
 else
    if (temperature  <=  0) read_temp = .true.
    rhof => rho_tab

    deallocate(rhof_n,rhof_rbreak)
    allocate(dens_prof(nprof),rad_prof(nprof),masstab(nprof))
    if (read_temp) allocate(temp_prof(nprof))

    !--Read profile from data
    if (profile_filename == default_name) then
       rad_prof = rad_prof_default
       dens_prof = dens_prof_default
    else
       open(iunit,file=profile_filename)
       if (.not. read_temp) then
          do iprof = 1,nprof
             read(iunit,*) rad_prof(iprof), dens_prof(iprof)
          enddo
       else
          do iprof = 1,nprof
             read(iunit,*) rad_prof(iprof), dens_prof(iprof), temp_prof(iprof)
          enddo
       endif
    endif
 endif
 ieos = ieos_in
 gmw = mu
 write(*,'(a,1x,i2)') ' Using eos =', ieos

 !--Everything to code unit
 ignore_radius = ignore_radius/udist
 if (use_func) then
    rad_min = rad_min/udist
    rad_max = rad_max/udist
    rhof_rbreak = rhof_rbreak/udist
    rhof_rhobreak = rhof_rhobreak/unit_density
    m_target = m_target*solarm/umass
    m_threshold = m_threshold*solarm/umass
 else
    rad_prof = rad_prof/udist
    dens_prof = dens_prof/unit_density
    rad_min = rad_prof(1)
    rad_max = rad_prof(nprof)
 endif

 !--Calc rho0 and rhobreak
 if (use_func) then
    if (rhof_rho0 < 0.) then
       call calc_rho0(rhof)
    elseif (m_target < 0.) then
       call calc_rhobreak()
    else
       call fatal('moddump','Must give rho0 or m_target')
    endif
 endif

 !--remove unwanted particles
 if (ignore_radius > 0) then
    npart_old = npart
    call delete_particles_inside_radius((/0.,0.,0./),ignore_radius,npart,npartoftype)
    write(*,'(I10,1X,A23,1X,E8.2,1X,A14)') npart_old - npart, 'particles inside radius', ignore_radius*udist, 'cm are deleted'
    npart_old = npart
    if (remove_overlap) then
       call delete_particles_outside_sphere((/0.,0.,0./),rad_min,npart)
       write(*,'(I10,1X,A24,1X,E8.2,1X,A14)') npart_old - npart, 'particles outside radius', rad_min*udist, 'cm are deleted'
       npart_old = npart
    endif
 else
    write(*,'(a)') ' Ignore all TDE particles'
    do i = 1,npart
       call kill_particle(i,npartoftype)
    enddo
    call shuffle_part(npart)
    npart_old = npart
 endif

 !--setup cloud
 totmass = get_mass_r(rhof,rad_max,rad_min)
 write(*,'(A42,1X,F5.2,1X,A10)') ' Total mass of the circumnuclear gas cloud:', totmass*umass/solarm, 'solar mass'
 np_sphere = nint(totmass/massoftype(igas))
 call set_sphere('random',id,master,rad_min,rad_max,delta,hfact_default,npart,xyzh, &
                 rhofunc=rhof,nptot=npart_total,exactN=.true.,np_requested=np_sphere,mask=i_belong)
 if (ierr /= 0) call fatal('moddump','error setting up the circumnuclear gas cloud')

 npartoftype(igas) = npart
 !--Set particle properties
 do i = npart_old+1,npart
    call set_particle_type(i,igas)
    r = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    rhofr = rhof(r)
    if (read_temp) temperature = get_temp_r(r,rad_prof,temp_prof)
    vxyzu(4,i) = uerg(rhofr,temperature,ieos)
    vxyzu(1:3,i) = 0. ! stationary for now
    pxyzu(4,i) = entropy(rhofr,temperature,ieos)
    pxyzu(1:3,i) = 0.
    eos_vars(itemp,i) = temperature
    presi = pressure(rhofr,temperature,ieos)
    eos_vars(igamma,i) = 1. + presi/(rhofr*vxyzu(4,i))
 enddo
 if (ieos == 12) write(*,'(a,1x,f10.4)') ' Mean gamma =', sum(eos_vars(igamma,npart_old+1:npart))/(npart - npart_old)

 !--Set timesteps
 tmax = 3.*years/utime
 dtmax = tmax/1000.

end subroutine modify_dump

!--Functions

real function rho(r)
 real, intent(in) :: r
 integer          :: i
 logical          :: found_rad

 found_rad = .false.
 do i = 1,nbreak-1
    if (r > rhof_rbreak(i) .and. r < rhof_rbreak(i+1)) then
       rho = rhof_rhobreak(i)*(r/rhof_rbreak(i))**rhof_n(i)
       found_rad = .true.
    endif
 enddo
 if (.not. found_rad) rho = rhof_rhobreak(nbreak)*(r/rhof_rbreak(nbreak))**rhof_n(nbreak)

end function rho

real function rho_tab(r)
 real, intent(in) :: r
 integer          :: i
 real             :: logr1,logr2,logr
 real             :: logrho1,logrho2,logrho_tab
 real             :: gradient

 rho_tab = 0.
 do i = 1,nprof-1
    if (r > rad_prof(i) .and. r < rad_prof(i+1)) then
       select case (interpolation)
       case ('log')
          logr1 = log10(rad_prof(i))
          logr2 = log10(rad_prof(i+1))
          logrho1 = log10(dens_prof(i))
          logrho2 = log10(dens_prof(i+1))
          logr = log10(r)
          gradient = (logrho2-logrho1)/(logr2-logr1)
          logrho_tab = logrho1 + gradient*(logr-logr1)
          rho_tab = 10**logrho_tab
       case ('lin')
          gradient = (dens_prof(i+1)-dens_prof(i))/(rad_prof(i+1)-rad_prof(i))
          rho_tab = dens_prof(i) + gradient*(r-rad_prof(i))
       case default
          write(*,'(a29,1x,a)') 'Unknown interpolation option:', trim(interpolation)
          write(*,'(a53)') "Support only 'lin'ear/'log'arithmic interpolation now"
       end select
    endif
 enddo
end function rho_tab

real function get_temp_r(r,rad_prof,temp_prof)
 real, intent(in) :: r,rad_prof(nprof),temp_prof(nprof)
 integer :: i
 real    :: t1,r1

 get_temp_r = temperature
 do i = 1,nprof
    if (r > rad_prof(i) .and. r < rad_prof(i+1)) then
       t1 = temp_prof(i)
       r1 = rad_prof(i)
       get_temp_r = (temp_prof(i+1)-t1)/(rad_prof(i+1)-r1)*(r-r1) + t1
       exit
    endif
 enddo

end function get_temp_r

real function uerg(rho,T,ieos)
 use physcon, only:kb_on_mh,radconst
 use units,   only:unit_density,unit_ergg
 real, intent(in) :: rho,T
 integer, intent(in) :: ieos
 real :: ucgs_gas,ucgs_rad,rhocgs

 rhocgs = rho*unit_density
 ucgs_gas = 1.5*kb_on_mh*T/mu
 if (ieos == 12) then
    ucgs_rad = radconst*T**4/rhocgs
 else
    ucgs_rad = 0. !radconst*T**4/rhocgs
 endif
 uerg = (ucgs_gas+ucgs_rad)/unit_ergg

end function uerg

real function entropy(rho,T,ieos)
 use physcon, only:kb_on_mh,radconst,kboltz
 use units,   only:unit_density,unit_ergg
 real, intent(in) :: rho,T
 integer, intent(in) :: ieos
 real :: ent_gas,ent_rad,rhocgs

 rhocgs = rho*unit_density
 ent_gas = kb_on_mh/mu*log(T**1.5/rhocgs)
 if (ieos == 12) then
    ent_rad = 4.*radconst*T**3/(3.*rhocgs)
 else
    ent_rad = 0.
 endif
 entropy = (ent_gas+ent_rad)/kboltz/ unit_ergg

end function entropy

real function pressure(rho,T,ieos)
 use physcon, only:kb_on_mh,radconst
 use units,   only:unit_density,unit_pressure
 real, intent(in) :: rho,T
 integer, intent(in) :: ieos
 real :: p_gas,p_rad,rhocgs

 rhocgs = rho*unit_density
 p_gas = rhocgs*kb_on_mh*T/mu
 if (ieos == 12) then
    p_rad = radconst*T**4/3.
 else
    p_rad = 0.
 endif
 pressure = (p_gas+p_rad)/ unit_pressure

end function pressure

subroutine calc_rhobreak()
 integer :: i

 rhof_rhobreak(1) = rhof_rho0
 if (nbreak > 1) then
    do i = 2,nbreak
       rhof_rhobreak(i) = rhof_rhobreak(i-1)*(rhof_rbreak(i)/rhof_rbreak(i-1))**rhof_n(i-1)
    enddo
 endif

end subroutine calc_rhobreak

subroutine calc_rho0(rhof)
 use units,      only:unit_density
 use stretchmap, only:get_mass_r
 procedure(rho), pointer, intent(in) :: rhof
 real    :: rho0_min,rho0_max,totmass
 integer :: iter

 rho0_min = 0.
 rho0_max = 1.
 totmass = -1.
 iter = 0

 do while (abs(totmass - m_target) > m_threshold)
    rhof_rho0 = 0.5*(rho0_min + rho0_max)
    call calc_rhobreak()
    totmass = get_mass_r(rhof,rad_max,rad_min)
    if (totmass > m_target) then
       rho0_max = rhof_rho0
    else
       rho0_min = rhof_rho0
    endif
    iter = iter + 1
 enddo
 write(*,'(a11,1x,es10.2,1x,a12,1x,i3,1x,a10)') ' Get rho0 =', rhof_rho0*unit_density, 'g/cm^-3 with', iter, 'iterations'

end subroutine calc_rho0

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer            :: i
 character(len=20)  :: rstr,nstr

 write(*,"(a)") ' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for setting up a circumnuclear gas cloud'

 write(iunit,"(/,a)") '# geometry'
 call write_inopt(ignore_radius,'ignore_radius','ignore tde particle inside this radius (-ve = ignore all for injection)',iunit)
 call write_inopt(remove_overlap,'remove_overlap','remove outflow particles overlap with circum particles',iunit)
 call write_inopt(use_func,'use_func','if use broken power law for density profile',iunit)
 if (use_func) then
    call write_inopt(rad_min,'rad_min','inner radius of the circumnuclear gas cloud',iunit)
    call write_inopt(rad_max,'rad_max','outer radius of the circumnuclear gas cloud',iunit)
    write(iunit,"(/,a)") '# density broken power law'
    call write_inopt(rhof_rho0,'rhof_rho0','density at rad_min (in g/cm^3) (-ve = ignore and calc for m_target)',iunit)
    call write_inopt(m_target,'m_target','target mass in circumnuclear gas cloud (in Msun) (-ve = ignore and use rho0)',iunit)
    call write_inopt(m_threshold,'m_threshold','threshold in solving rho0 for m_target (in Msun)',iunit)
    call write_inopt(nbreak,'nbreak','number of broken power laws',iunit)
    write(iunit,"(/,a)") '#    section 1 (from rad_min)'
    call write_inopt(rhof_n(1),'rhof_n_1','power law index of the section',iunit)
    if (nbreak > 1) then
       do i=2,nbreak
          write(iunit,"(a,1x,i1)") '#    section',i
          write(rstr,'(a12,i1)') 'rhof_rbreak_',i
          write(nstr,'(a7,i1)') 'rhof_n_',i
          call write_inopt(rhof_rbreak(i),trim(rstr),'inner radius of the section',iunit)
          call write_inopt(rhof_n(i),trim(nstr),'power law index of the section',iunit)
       enddo
    endif
 else
    call write_inopt(profile_filename,'profile_filename','filename for the cloud profile',iunit)
    call write_inopt(nprof,'nprof','number of data points in the cloud profile',iunit)
    call write_inopt(interpolation,'interpolation',"use 'lin'ear/'log'arithmic interpolation between data points",iunit)
 endif

 write(iunit,"(/,a)") '# eos'
 call write_inopt(ieos_in,'ieos','equation of state used',iunit)
 call write_inopt(temperature,'temperature','temperature of the gas cloud (-ve = read from file)',iunit)
 call write_inopt(mu,'mu','mean molecular density of the cloud',iunit)

 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: fatal
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit=21,in_num=50
 integer                       :: i
 type(inopts), allocatable     :: db(:)
 character(len=20)             :: rstr,nstr
 real                          :: use_func_test

 write(*,"(a)")'  reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(ignore_radius,'ignore_radius',db,min=0.,err=ierr)
 call read_inopt(remove_overlap,'remove_overlap',db,err=ierr)
 call read_inopt(use_func,'use_func',db,err=ierr)
 call read_inopt(use_func_test,'nbreak',db,err=ierr)
 if (ierr == -1) use_func_old = .false.
 if (use_func) then
    call read_inopt(rad_min,'rad_min',db,min=ignore_radius,err=ierr)
    call read_inopt(rad_max,'rad_max',db,min=rad_min,err=ierr)
    call read_inopt(rhof_rho0,'rhof_rho0',db,err=ierr)
    call read_inopt(m_target,'m_target',db,err=ierr)
    call read_inopt(m_threshold,'m_threshold',db,err=ierr)
    call read_inopt(nbreak,'nbreak',db,min=1,err=ierr)
    allocate(rhof_rbreak_in(in_num),rhof_n_in(in_num))
    call read_inopt(rhof_n_in(1),'rhof_n_1',db,err=ierr)
    rhof_rbreak_in(1) = rad_min
    do i=2,nbreak+1
       write(rstr,'(a12,i1)') 'rhof_rbreak_',i
       write(nstr,'(a7,i1)') 'rhof_n_',i
       call read_inopt(rhof_rbreak_in(i),trim(rstr),db,min=rhof_rbreak_in(i-1),max=rad_max,err=ierr)
       call read_inopt(rhof_n_in(i),trim(nstr),db,err=ierr)
       if (ierr == 0) nbreak_old = i
    enddo
 else
    call read_inopt(profile_filename,'profile_filename',db,err=ierr)
    call read_inopt(nprof,'nprof',db,min=1,err=ierr)
    call read_inopt(interpolation,'interpolation',db,err=ierr)
 endif

 call read_inopt(ieos_in,'ieos',db,err=ierr)
 call read_inopt(temperature,'temperature',db,err=ierr)
 call read_inopt(mu,'mu',db,err=ierr)

 call close_db(db)

 if (ierr /= 0) then
    call fatal('moddump_radiotde','Error in reading setup file')
 endif

end subroutine read_setupfile

end module moddump

