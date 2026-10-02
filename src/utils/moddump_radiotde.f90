!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! :Dependencies: dynamic_dtmax, eos, infile_utils, io, kernel,
!   moddump_utils, mpidomain, part, physcon, setup_params, spherical,
!   stretchmap, timestep, units
!
 use moddump_utils, only:prompt_for_params,write_moddump_header
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 public :: modify_dump
 private :: rho,rho_tab,get_temp_r,uerg,calc_rhobreak,calc_rho0

 private
 character(len=*), parameter :: default_name = 'default_profile'
 integer           :: ieos_in,nprof,nbreak
 real              :: temperature,mu,ignore_radius,rad_max,rad_min
 character(len=50) :: profile_filename
 character(len=3)  :: interpolation
 real, allocatable :: rhof_n(:),rhof_rbreak(:),rhof_rhobreak(:)
 real, allocatable :: rhof_n_in(:),rhof_rbreak_in(:)
 real, allocatable :: rad_prof(:),dens_prof(:)
 real              :: rhof_rho0,m_target,m_threshold
 logical           :: use_func,remove_overlap

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .false.
 public :: moddump_interactive

contains

!----------------------------------------------------------------
!
!  Sets up a circumnuclear gas cloud
!
!----------------------------------------------------------------
subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use physcon,      only:solarm,years
 use setup_params, only:npart_total
 use part,         only:igas,set_particle_type,pxyzu,delete_particles_inside_radius, &
                        delete_particles_outside_sphere,kill_particle,shuffle_part, &
                        eos_vars,itemp,igamma,igasP
 use io,           only:fatal,master,id
 use units,        only:umass,udist,utime,unit_density
 use timestep,     only:dtmax,tmax
 use dynamic_dtmax,only:idtmax_frac,dtmax_ifactor,idtmax_n
 use eos,          only:ieos,gmw
 use kernel,       only:hfact_default
 use stretchmap,   only:get_mass_r,rho_func
 use spherical,    only:set_sphere
 use mpidomain,    only:i_belong
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(inout) :: massoftype(:)
 integer                       :: i,iunit,iprof
 integer                       :: np_sphere,npart_old
 real                          :: totmass,delta,r,rhofr,presi
 logical                       :: read_temp
 real, allocatable             :: masstab(:),temp_prof(:)
 real, dimension(7), parameter :: dens_prof_default = (/8.9e-21, 5.1e-21, 3.3e-21, 2.6e-21, &
                                                        6.6e-25, 3.4e-25, 8.1e-26/), &
                                  rad_prof_default = (/8.7e16, 1.2e17, 1.4e17, 2.0e17, &
                                                       4.0e17, 4.8e17, 7.1e17/) ! profile from Cendes+2021
 procedure(rho_func), pointer  :: rhof

 read_temp = .false.

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
       open(newunit=iunit,file=profile_filename)
       if (.not. read_temp) then
          do iprof = 1,nprof
             read(iunit,*) rad_prof(iprof), dens_prof(iprof)
          enddo
       else
          do iprof = 1,nprof
             read(iunit,*) rad_prof(iprof), dens_prof(iprof), temp_prof(iprof)
          enddo
       endif
       close(iunit)
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
    write(*,'(i10,1x,a23,1x,es10.2,1x,a14)') npart_old - npart, 'particles inside radius', ignore_radius*udist, 'cm are deleted'
    npart_old = npart
    if (remove_overlap) then
       call delete_particles_outside_sphere((/0.,0.,0./),rad_min,npart)
       write(*,'(i10,1x,a24,1x,es10.2,1x,a14)') npart_old - npart, 'particles outside radius', rad_min*udist, 'cm are deleted'
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
 write(*,'(a42,1x,f5.2,1x,a10)') ' Total mass of the circumnuclear gas cloud:', totmass*umass/solarm, 'solar mass'
 np_sphere = nint(totmass/massoftype(igas))
 call set_sphere('random',id,master,rad_min,rad_max,delta,hfact_default,npart,xyzh, &
                 rhofunc=rhof,nptot=npart_total,exactN=.true.,np_requested=np_sphere,mask=i_belong)

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
 dtmax_ifactor = 0
 idtmax_frac = 0 ! so don't write to .restart
 idtmax_n = 1

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
 real,    intent(in) :: rho,T
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
 real,    intent(in) :: rho,T
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
 real,    intent(in) :: rho,T
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
 use stretchmap, only:get_mass_r,rho_func
 procedure(rho_func), pointer, intent(in) :: rhof
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
!  write parameters to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer            :: i
 character(len=20)  :: rstr,nstr

 !--make sure the broken-power-law arrays match nbreak, which may have been
 !  changed via the file (this is what triggers a top-up + rerun)
 if (use_func) then
    if (.not.allocated(rhof_n) .or. size(rhof_n) /= nbreak) then
       if (allocated(rhof_n)) deallocate(rhof_n)
       if (allocated(rhof_rbreak)) deallocate(rhof_rbreak)
       allocate(rhof_n(nbreak),rhof_rbreak(nbreak))
       rhof_n = -1.7
       rhof_rbreak = rad_min
    endif
 endif

 write(*,"(a)") ' writing moddump options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
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

end subroutine write_moddump

!----------------------------------------------------------------
!+
!  Read parameters from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit=21,in_num=50
 integer                       :: i,nerr
 type(inopts), allocatable     :: db(:)
 character(len=20)             :: rstr,nstr

 nerr = 0
 write(*,"(a)")'  reading moddump options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return

 call read_inopt(ignore_radius,'ignore_radius',db,min=0.,errcount=nerr)
 call read_inopt(remove_overlap,'remove_overlap',db,errcount=nerr)
 call read_inopt(use_func,'use_func',db,errcount=nerr)
 if (use_func) then
    call read_inopt(rad_min,'rad_min',db,min=ignore_radius,errcount=nerr)
    call read_inopt(rad_max,'rad_max',db,min=rad_min,errcount=nerr)
    call read_inopt(rhof_rho0,'rhof_rho0',db,errcount=nerr)
    call read_inopt(m_target,'m_target',db,errcount=nerr)
    call read_inopt(m_threshold,'m_threshold',db,errcount=nerr)
    call read_inopt(nbreak,'nbreak',db,min=1,errcount=nerr)
    if (allocated(rhof_rbreak_in)) deallocate(rhof_rbreak_in)
    if (allocated(rhof_n_in)) deallocate(rhof_n_in)
    allocate(rhof_rbreak_in(in_num),rhof_n_in(in_num))
    call read_inopt(rhof_n_in(1),'rhof_n_1',db,errcount=nerr)
    rhof_rbreak_in(1) = rad_min
    do i=2,nbreak
       write(rstr,'(a12,i1)') 'rhof_rbreak_',i
       write(nstr,'(a7,i1)') 'rhof_n_',i
       call read_inopt(rhof_rbreak_in(i),trim(rstr),db,min=rhof_rbreak_in(i-1),max=rad_max,errcount=nerr)
       call read_inopt(rhof_n_in(i),trim(nstr),db,errcount=nerr)
    enddo
 else
    call read_inopt(profile_filename,'profile_filename',db,errcount=nerr)
    call read_inopt(nprof,'nprof',db,min=1,errcount=nerr)
    call read_inopt(interpolation,'interpolation',db,errcount=nerr)
 endif

 call read_inopt(ieos_in,'ieos',db,errcount=nerr)
 call read_inopt(temperature,'temperature',db,errcount=nerr)
 call read_inopt(mu,'mu',db,errcount=nerr)

 call close_db(db)
 if (nerr > 0) ierr = nerr

 ! Keep the serialised profile in input units before modify_dump converts it.
 if (use_func .and. nerr == 0) then
    rhof_n = rhof_n_in(1:nbreak)
    rhof_rbreak = rhof_rbreak_in(1:nbreak)
 endif

end subroutine read_moddump

subroutine init_moddump()
 use part, only:npartoftype,massoftype
 use units, only:umass
 use physcon, only:solarm,mass_proton_cgs

 !--Set default values
 temperature       = 10.           ! Temperature in Kelvin
 mu                = 1.            ! mean molecular weight
 ieos_in           = 2
 ignore_radius     = 1.e14          ! in cm
 use_func          = .true.
 remove_overlap    = .true.
 !--Power law default setups
 rad_max           = 7.1e16        ! in cm
 rad_min           = 8.7e15        ! in cm
 nbreak            = 1
 rhof_rho0         = 1.e4*mu*mass_proton_cgs
 if (allocated(rhof_n)) deallocate(rhof_n)
 if (allocated(rhof_rbreak)) deallocate(rhof_rbreak)
 allocate(rhof_n(nbreak),rhof_rbreak(nbreak))
 rhof_n            = -1.7
 rhof_rbreak       = rad_min
 m_target          = dot_product(npartoftype,massoftype)*umass/solarm
 m_threshold       = 1.e-3

 !--Profile default setups
 profile_filename  = default_name
 nprof             = 7
 interpolation     = 'log'
end subroutine init_moddump

end module moddump
