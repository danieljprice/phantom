!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module AmusePhantom
!
! AmusePhantom
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, boundary_dyn, cooling, damping, deriv, dim,
!   dust, dust_formation, energies, eos, evolve, gravwaveutils, growth,
!   initial, inject, io, kdtree, memory, metric, mpiutils, nicil_sup,
!   options, part, partinject, physcon, ptmass, ptmass_radiation,
!   radiation_implicit, radiation_utils, step_lf_global, timestep,
!   timestep_ind, units, viscosity
!

 ! Currently, AMUSE only supports up to 32 bit integers for indices.
 ! This may change in the future (see https://github.com/amusecode/amuse/issues/1077)
 ! Until then, we have to ensure that we can safely use 32 bit integers here
 use, intrinsic :: ISO_FORTRAN_ENV, only: INT32, INT64
 implicit none
 integer, parameter :: index_length=INT32
 integer(kind=INT64), parameter :: min_int32 = -2_INT64**31
 integer(kind=INT64), parameter :: max_int32 = 2_INT64**31 - 1

 integer(kind = index_length), allocatable:: amuse_id_lookup(:)
 integer(kind = index_length):: new_particles_since_last_update = 0
 integer(kind = index_length):: particles_added_by_amuse = 0
contains

subroutine construct_id_lookup()
 ! amuse_id_lookup needs to be rebuilt if/when particles are deleted/added/reshuffled
 ! easier (though a bit slower) is to do it at the end of an evolve step, and at recommit_particles
 use dim, only: maxp
 use part, only: iorig, norig
 implicit none
 integer (kind = index_length):: i, j
 integer (kind=index_length):: norig_amuse
 integer (kind=8):: tmp
 write(*,*) "Rebuilding lookup table"
 do i = 1, maxp
    amuse_id_lookup(i) = 0
 enddo
 ! This is a workaround for the fact that norig is an integer(kind=8)
 norig_amuse = int(norig, kind=index_length)
 do i = 1, norig_amuse
    tmp = iorig(i)
    if (tmp < min_int32 .or. tmp > max_int32) then
       error stop "iorig out of range for AMUSE"
    else
       j = int(tmp, kind=index_length)
    endif
    if (j > 0) then
       amuse_id_lookup(j) = i
    endif
 enddo
 write(*,*) size(amuse_id_lookup), "?=", norig, maxp
 write(*,*) "Lookup table rebuilt"
end subroutine construct_id_lookup

subroutine amuse_initialize_code()
 use dim,             only:maxp
 use allocutils,      only:allocate_array
 use mpiutils,        only:init_mpi
 use io,              only:id, nprocs
 implicit none

 id = 0

 call init_mpi(id, nprocs)
 call allocate_array('amuse_id_lookup', amuse_id_lookup, maxp)
end subroutine amuse_initialize_code

subroutine amuse_set_phantom_option(name, valstring, imatch)
 ! This subroutine is meant to be a replacement for read_infile
 ! It should be kept up to date and support all options listed there, or raise an error if something is not caught.
 use inject,           only:read_options_inject
 use dust,             only:read_options_dust
 use growth,           only:read_options_growth
 use metric,           only:read_options_metric
 use nicil_sup,        only:read_options_nicil
 use dust_formation,   only:read_options_dust_formation
 use ptmass_radiation, only:read_options_ptmass_radiation
 use eos,              only:read_options_eos
 use cooling,          only:read_options_cooling
 use ptmass,           only:read_options_ptmass
 use damping,          only:read_options_damping
 use gravwaveutils,    only:read_options_gravitationalwaves
 use boundary_dyn,     only:read_options_boundary
 use timestep, only:tmax, dtmax, C_cour, C_force, C_ent, xtol, tolv, ptol, nout, nmax, &
            dtwallmax, dtmax_min, dtmax_max, dtmax_dratio
 use options,  only:alpha, alphaB, alphamax, alphau, twallmax, tolh, rkill, rhofinal_cgs, &
            psidecayfac, overcleanfac, nmaxdumps, nfulldump, limit_radiation_flux, &
            ishock_heating, iresistive_heating, ireconav, ipdv_heating, iopacity_type, &
            implicit_radiation, exchange_radiation_energy, calc_erot, &
            beta, avdecayconst
 use radiation_implicit, only:tol_rad, itsmax_rad, cv_type
 use radiation_utils,    only:kappa_cgs
 use dim, only:store_dust_temperature
 use viscosity, only:shearparam, irealvisc, bulkvisc
 use io,        only:iverbose
 use part,      only:hfact, ien_type
 implicit none
 character(*), intent(inout):: name, valstring
 logical:: imatch, igotall
 integer:: ierr
 real:: ratio
 logical:: incl_runtime2 = .false.
 imatch = .true.
 select case(trim(name))
 case('logfile')
    ! ignored
 case('dumpfile')
    ! ignored
 case('tmax')
    write(*,*) "Set tmax to ", valstring
    read(valstring, *,iostat = ierr) tmax
 case('dtmax')
    read(valstring, *,iostat = ierr) dtmax
 case('nmax')
    read(valstring, *,iostat = ierr) nmax
 case('nout')
    read(valstring, *,iostat = ierr) nout
 case('nmaxdumps')
    read(valstring, *,iostat = ierr) nmaxdumps
 case('twallmax')
    read(valstring, *,iostat = ierr) twallmax
 case('dtwallmax')
    read(valstring, *,iostat = ierr) dtwallmax
 case('iverbose')
    read(valstring, *,iostat = ierr) iverbose
    write(*,*) "IVERBOSE: ", iverbose
 case('rhofinal_cgs')
    read(valstring, *,iostat = ierr) rhofinal_cgs
    incl_runtime2 = .true.
 case('calc_erot')
    read(valstring, *,iostat = ierr) calc_erot
    incl_runtime2 = .true.
 case('dtmax_dratio')
    read(valstring, *,iostat = ierr) dtmax_dratio
    incl_runtime2 = .true.
 case('dtmax_max')
    read(valstring, *,iostat = ierr) dtmax_max
    if (dtmax_max <= 0.0) dtmax_max = dtmax
    ! to prevent comparison errors from round-off
    ratio = dtmax_max/dtmax
    ratio = int(ratio+0.5)+0.0001
    dtmax_max = dtmax*ratio
 case('dtmax_min')
    read(valstring, *,iostat = ierr) dtmax_min
    ! to prevent comparison errors from round-off
    if (dtmax_min > epsilon(dtmax_min)) then
       ratio = dtmax/dtmax_min
       ratio = int(ratio+0.5)+0.0001
       dtmax_min = dtmax/ratio
    endif
 case('C_cour')
    read(valstring, *,iostat = ierr) C_cour
 case('C_force')
    read(valstring, *,iostat = ierr) C_force
 case('tolv')
    read(valstring, *,iostat = ierr) tolv
 case('C_ent')
    read(valstring, *,iostat = ierr) C_ent
 case('xtol')
    read(valstring, *,iostat = ierr) xtol
 case('ptol')
    read(valstring, *,iostat = ierr) ptol
 case('hfact')
    read(valstring, *,iostat = ierr) hfact
 case('tolh')
    read(valstring, *,iostat = ierr) tolh
 case('rkill')
    read(valstring, *,iostat = ierr) rkill
 case('nfulldump')
    read(valstring, *,iostat = ierr) nfulldump
 case('alpha')
    read(valstring, *,iostat = ierr) alpha
 case('alphamax')
    read(valstring, *,iostat = ierr) alphamax
 case('alphau')
    read(valstring, *,iostat = ierr) alphau
 case('alphaB')
    read(valstring, *,iostat = ierr) alphaB
 case('psidecayfac')
    read(valstring, *,iostat = ierr) psidecayfac
 case('overcleanfac')
    read(valstring, *,iostat = ierr) overcleanfac
 case('beta')
    read(valstring, *,iostat = ierr) beta
 case('ireconav')
    read(valstring, *,iostat = ierr) ireconav
 case('avdecayconst')
    read(valstring, *,iostat = ierr) avdecayconst
 case('ipdv_heating')
    read(valstring, *,iostat = ierr) ipdv_heating
 case('ishock_heating')
    read(valstring, *,iostat = ierr) ishock_heating
 case('iresistive_heating')
    read(valstring, *,iostat = ierr) iresistive_heating
 case('ien_type')
    read(valstring, *,iostat = ierr) ien_type
 case('irealvisc')
    read(valstring, *,iostat = ierr) irealvisc
 case('shearparam')
    read(valstring, *,iostat = ierr) shearparam
 case('bulkvisc')
    read(valstring, *,iostat = ierr) bulkvisc
#ifdef MCFOST
 case('use_mcfost')
    read(valstring, *,iostat = ierr) use_mcfost
 case('Voronoi_limits_file')
    read(valstring, *,iostat = ierr) Voronoi_limits_file
    use_Voronoi_limits_file = .true.
 case('use_mcfost_stars')
    read(valstring, *,iostat = ierr) use_mcfost_stellar_parameters
 case('mcfost_computes_Lacc')
    read(valstring, *,iostat = ierr) mcfost_computes_Lacc
 case('mcfost_uses_PdV')
    read(valstring, *,iostat = ierr) mcfost_uses_PdV
 case('mcfost_keep_part')
    read(valstring, *,iostat = ierr) mcfost_keep_part
 case('ISM')
    read(valstring, *,iostat = ierr) ISM
 case('mcfost_dust_subl')
    read(valstring, *,iostat = ierr) mcfost_dust_subl
#endif
 case('implicit_radiation')
    read(valstring, *,iostat = ierr) implicit_radiation
    if (implicit_radiation) store_dust_temperature = .true.
 case('gas-rad_exchange')
    read(valstring, *,iostat = ierr) exchange_radiation_energy
 case('flux_limiter')
    read(valstring, *,iostat = ierr) limit_radiation_flux
 case('iopacity_type')
    read(valstring, *,iostat = ierr) iopacity_type
 case('cv_type')
    read(valstring, *,iostat = ierr) cv_type
 case('kappa_cgs')
    read(valstring, *,iostat = ierr) kappa_cgs
 case('tol_rad')
    read(valstring, *,iostat = ierr) tol_rad
 case('itsmax_rad')
    read(valstring, *,iostat = ierr) itsmax_rad
 case default
    imatch = .false.
    if (.not.imatch) call read_options_inject(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_dust_formation(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_ptmass_radiation(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_dust(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_growth(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_metric(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_nicil(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_eos(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_cooling(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_damping(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_ptmass(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_gravitationalwaves(name, valstring, imatch, igotall, ierr)
    if (.not.imatch) call read_options_boundary(name, valstring, imatch, igotall, ierr)
 end select
 if (.not.imatch) write(*,*) "Could not set option ", name, ", please check if this is a problem!"
end subroutine amuse_set_phantom_option

subroutine amuse_initialize_wind()
 ! instead of reading a wind setup, set values here
 use inject, only:read_options_inject
 use dust_formation, only: read_options_dust_formation
 use ptmass_radiation, only: read_options_ptmass_radiation
 logical:: imatch, igotall
 integer:: ierr

 call read_options_inject("sonic_type", "0", imatch, igotall, ierr)
 call read_options_inject("wind_velocity", "20.", imatch, igotall, ierr)
 call read_options_inject("wind_inject_radius", "2.000", imatch, igotall, ierr)    ! wind injection radius (au, if 0 takes Rstar)
 call read_options_inject("wind_mass_rate", "1.000E-05", imatch, igotall, ierr)    ! wind mass loss rate (Msun/yr)
 call read_options_inject("wind_temperature", "2500.", imatch, igotall, ierr)    ! wind temperature at injection radius (K, if 0 takes Teff)
 call read_options_inject("iwind_resolution", "5", imatch, igotall, ierr)    ! if<>0 set number of particles on the sphere, reset particle mass
 call read_options_inject("nfill_domain", "0", imatch, igotall, ierr)    ! number of spheres used to set the background density profile
 call read_options_inject("wind_shell_spacing", "1.000", imatch, igotall, ierr)    ! desired ratio of sphere spacing to particle spacing
 call read_options_inject("iboundary_spheres", "5", imatch, igotall, ierr)    ! number of boundary spheres (integer)
 call read_options_inject("outer_boundary", "30.", imatch, igotall, ierr)    ! delete gas particles outside this radius (au)
 call read_options_inject("rkill", "-1.000", imatch, igotall, ierr)    ! deactivate particles outside this radius (<0 is off)

 !# options controlling dust
 call read_options_dust_formation("idust_opacity", "0", imatch, igotall, ierr)    ! compute dust opacity (0 = off, 1 (bowen), 2 (nucleation))

 !# options controlling radiation pressure from sink particles
 call read_options_ptmass_radiation("isink_radiation", "1", imatch, igotall, ierr)    ! sink radiation pressure method (0 = off, 1 = alpha, 2 = dust, 3 = alpha+dust)
 call read_options_ptmass_radiation("alpha_rad", "1.000", imatch, igotall, ierr)

end subroutine amuse_initialize_wind

subroutine amuse_commit_parameters()
 use memory,  only:allocate_memory
 use dim,     only:maxp
 use io,      only:set_io_unit_numbers
 use initial, only:initialise
 use units,   only:set_units, set_units_extra
 use physcon, only:solarm, au

 call allocate_memory(int(maxp, kind = 8))
 call set_io_unit_numbers()
 call initialise()
 call amuse_set_defaults()  ! replaces reading infile

 ! Hard coded defaults, should be set in a different way...
 call set_units(dist = 1*au, mass = 1.*solarm, G = 1.)
 call set_units_extra()
end subroutine amuse_commit_parameters

subroutine amuse_commit_particles()
 use part, only: norig
 use initial, only:startrun
 use timestep, only:nmax, nsteps
 use evolve, only:evol
 implicit none
 character(len = 120):: infile, logfile, evfile, dumpfile
 integer:: nsteps_orig
 integer(kind=index_length):: norig_amuse
 call startrun(infile, logfile, evfile, dumpfile, .true.)

 ! Make sure the evol call only initialises and doesn't do an evolve step
 nsteps_orig = nsteps
 nsteps = nmax
 call evol(infile, logfile, evfile, dumpfile)
 nsteps = nsteps_orig

 call construct_id_lookup()
 ! This is a workaround for the fact that norig is an integer(kind=8)
 if (norig < min_int32 .or. norig > max_int32) then
    error stop "norig out of range for AMUSE"
 endif
 norig_amuse = int(norig, kind=index_length)
 new_particles_since_last_update = norig_amuse-particles_added_by_amuse
end subroutine amuse_commit_particles

subroutine amuse_recommit_particles()
 use deriv, only:get_derivs_global
 !use timestep, only:dtmax
 implicit none
 !double precision:: dt, dtnew_first
 !dtnew_first = dtmax
 !dt = 0

 call get_derivs_global()  ! optional: get_derivs_global(tused, dt_new, dt)
 call construct_id_lookup()
end subroutine amuse_recommit_particles

subroutine amuse_cleanup_code()
 use initial, only:endrun
 implicit none

 call endrun()
 if (allocated(amuse_id_lookup)) deallocate(amuse_id_lookup)
end subroutine amuse_cleanup_code

subroutine amuse_get_norig(norig_out)
 use part, only:norig
 implicit none
 integer(kind = 8), intent(out):: norig_out
 norig_out = norig
end subroutine amuse_get_norig

! Get npart
subroutine amuse_get_npart(npart_out, nodisabled)
 use part, only:npart, xyzh, isdead_or_accreted
 implicit none
 integer, intent(out):: npart_out
 logical, intent(in)  :: nodisabled

 integer:: i

 if (nodisabled) then
    npart_out = 0
    do i = 1, npart
       if (isdead_or_accreted(xyzh(4, i))) then
          npart_out = npart_out+1
       endif
    enddo
 else
    npart_out = npart
 endif
end subroutine amuse_get_npart

! Set default parameters
subroutine amuse_set_defaults()
 use options, only: set_default_options, write_files
 implicit none

 call set_default_options()
 write_files=.false.

end subroutine amuse_set_defaults


! New particles
subroutine amuse_new_sph_particle(i, mass, x, y, z, vx, vy, vz, u, h)
 use part, only:igas, npart, npartoftype, xyzh, vxyzu, massoftype
 !use part, only:abundance, iHI, ih2ratio
 use partinject, only:add_or_update_particle
#ifdef IND_TIMESTEPS
 use part, only:twas, ibin
 use timestep_ind, only:istepfrac, get_dt, nbinmax, change_nbinmax, get_newbin
 use timestep, only:dt, time, dtmax
 use timestep, only:C_cour, rhomaxnow
 use eos, only:gamma
#endif
 use timestep, only:dtextforce
 use units, only:umass, udist, utime
 implicit none
 integer:: i, itype
 double precision:: mass, x, y, z, vx, vy, vz, u, h
 double precision:: position(3), velocity(3)
#ifdef IND_TIMESTEPS
 integer(kind = 1):: nbinmaxprev
 real:: dtinject
 dtinject = huge(dtinject)
 ! dtmax = 0.01  ! TODO This is arbitrarily set. Probably bad.
#endif
 ! Adding a particle of unknown temperature -> use a small cooling step
 dtextforce = 1.e-6

 itype = igas
 i = npart+1
 ! print*, "**************   X position of particle ", i, x, "   ***************"
 position(1) = x
 position(2) = y
 position(3) = z
 velocity(1) = vx
 velocity(2) = vy
 velocity(3) = vz

 if (npartoftype(itype) == 0) then
    massoftype(itype) = mass
 endif
 call add_or_update_particle(itype, position, velocity, h, &
        u, i, npart, npartoftype, xyzh, vxyzu)
 !abundance(:,i) = 0.
 !abundance(ih2ratio, i) = 0.5
 !abundance(iHI, i) = 1.  ! assume all gas is atomic hydrogen initially
 if (i == 1) then
    print*, "xyz vxyz u mass = ", x, y, z, vx, vy, vz, u, mass
    print*, "udist, utime, umass = ", udist, utime, umass
    print*, "x vx u mass in cgs = ", x*udist, vx*udist/utime, u*udist*udist/utime/utime, mass*umass
 endif

#ifdef IND_TIMESTEPS
 dtinject = C_cour*h / (gamma*(gamma-1)*u)**0.5
 nbinmaxprev = nbinmax
 call get_newbin(dtinject, dtmax, nbinmax, allow_decrease=.false.)
 !! not doing clever stuff: all particles go in the shortest possible bin.
 !! FIXME rethink this later...
 nbinmax = 3
 !if (nbinmax > nbinmaxprev) then  ! update number of bins if needed
 !   call change_nbinmax(nbinmax, nbinmaxprev, istepfrac, dtmax, dt)
 !   print*, "nbinmax (prev), time: ", nbinmax, nbinmaxprev, time
 !   print*, "npart:", npart
 !endif
 ! put all injected particles on shortest bin
 ibin(i) = nbinmax
 twas(i) = time+0.5*get_dt(dtmax, ibin(i))
 if (i == 5) then
    print*, "particle ", i, ": "
    print*, x, y, z, vx, vy, vz, mass, ibin(i), twas(i)
 endif
#endif
 particles_added_by_amuse = particles_added_by_amuse+1
end subroutine amuse_new_sph_particle

subroutine amuse_new_dm_particle(i, mass, x, y, z, vx, vy, vz)
 use part, only:idarkmatter, npart, npartoftype, xyzh, vxyzu, massoftype
 use partinject, only:add_or_update_particle
 implicit none
 integer:: i
 integer:: itype
 double precision:: mass, x, y, z, vx, vy, vz, h_smooth, u
 double precision:: position(3), velocity(3)

 u = 0
 itype = idarkmatter
 i = npart+1
 position(1) = x
 position(2) = y
 position(3) = z
 velocity(1) = vx
 velocity(2) = vy
 velocity(3) = vz
 h_smooth = 0.1  ! TODO set this to some default
 if (npartoftype(itype) == 0) then
    massoftype(itype) = mass
 endif

 call add_or_update_particle(itype, position, velocity, h_smooth, &
        u, i, npart, npartoftype, xyzh, vxyzu)
end subroutine amuse_new_dm_particle

subroutine amuse_new_sink_particle(j, mass, x, y, z, vx, vy, vz, &
        radius, accretion_radius, h_smooth)
 use io, only:fatal
 use part, only:nptmass, maxptmass, xyzmh_ptmass, vxyz_ptmass, ihacc, ihsoft, iReff
 implicit none
 integer:: i, j
 double precision:: mass, x, y, z, vx, vy, vz, radius, accretion_radius, h_smooth
 nptmass = nptmass+1
 ! Replace this fatal exception with something AMUSE can handle
 if (nptmass > maxptmass) call fatal('creating new sink', 'nptmass > maxptmass')
 i = nptmass

 ! Sink particles are stored in different arrays than other particles.
 ! To be able to distinguish the particle indices on the AMUSE side, we give sinks
 ! a negative index and other particles a positive index.
 j = -i

 xyzmh_ptmass(:,i) = 0.
 xyzmh_ptmass(1, i) = x
 xyzmh_ptmass(2, i) = y
 xyzmh_ptmass(3, i) = z
 xyzmh_ptmass(4, i) = mass
 xyzmh_ptmass(iReff, i) = radius
 xyzmh_ptmass(ihacc, i) = accretion_radius
 xyzmh_ptmass(ihsoft, i) = h_smooth
 vxyz_ptmass(1, i) = vx
 vxyz_ptmass(2, i) = vy
 vxyz_ptmass(3, i) = vz
end subroutine amuse_new_sink_particle


subroutine amuse_delete_particle(i)
 use part, only:kill_particle, xyzmh_ptmass
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 if (i == abs(i)) then
    write(*,*) "AMUSE killing a gas particle?"
    call amuse_get_index(i, part_index)
    ! call kill_particle(part_index)
 else
    write(*,*) "AMUSE killing a sink particle?"
    ! Sink particles can't be killed-so we just set its mass to negative
    !xyzmh_ptmass(4, -i) = -1
    xyzmh_ptmass(4, -i) = -abs(xyzmh_ptmass(4, i))
 endif
end subroutine amuse_delete_particle

subroutine amuse_get_unit_length(unit_length_out)
 use units, only: udist
 implicit none
 double precision, intent(out):: unit_length_out
 unit_length_out = udist
end subroutine amuse_get_unit_length

subroutine amuse_set_unit_length(unit_length_in)
 use units, only: udist, utime, umass, set_units
 implicit none
 double precision, intent(in):: unit_length_in
 !udist = unit_length_in
 !call set_units(dist = udist, time = utime, G = 1.)
 print*, "set_unit_length called: utime/mass/dist: ", utime, umass, udist
end subroutine amuse_set_unit_length

subroutine amuse_get_unit_mass(unit_mass_out)
 use units, only: umass
 implicit none
 double precision, intent(out):: unit_mass_out
 unit_mass_out = umass
end subroutine amuse_get_unit_mass

subroutine amuse_set_unit_mass(unit_mass_in)
 use units, only: umass, utime, udist, set_units
 implicit none
 double precision, intent(in):: unit_mass_in
 !umass = unit_mass_in
 !call set_units(mass = umass, time = utime, G = 1.)
 print*, "set_unit_mass called: utime/mass/dist: ", utime, umass, udist
end subroutine amuse_set_unit_mass

subroutine amuse_get_unit_time(unit_time_out)
 use units, only: utime
 implicit none
 double precision, intent(out):: unit_time_out
 unit_time_out = utime
end subroutine amuse_get_unit_time

subroutine amuse_set_unit_time(unit_time_in)
 use units, only: utime, umass, udist, set_units
 implicit none
 double precision, intent(out):: unit_time_in
 !utime = unit_time_in
 !call set_units(time = utime, mass = umass, G = 1.)
 print*, "set_unit_time called: utime/mass/dist: ", utime, umass, udist
end subroutine amuse_set_unit_time

subroutine amuse_get_constant_solarm(solarm_out)
 use physcon, only:solarm
 implicit none
 double precision, intent(out):: solarm_out
 solarm_out = solarm
end subroutine amuse_get_constant_solarm

subroutine amuse_get_constant_pc(pc_out)
 use physcon, only:pc
 implicit none
 double precision, intent(out):: pc_out
 pc_out = pc
end subroutine amuse_get_constant_pc

subroutine amuse_get_constant_planckh(planckh_out)
 use physcon, only:planckh
 implicit none
 double precision, intent(out):: planckh_out
 planckh_out = planckh
end subroutine amuse_get_constant_planckh

subroutine amuse_get_potential_energy(epot_out)
 use energies, only:epot
 implicit none
 double precision, intent(out):: epot_out
 epot_out = epot
end subroutine amuse_get_potential_energy

subroutine amuse_get_kinetic_energy(ekin_out)
 use energies, only:ekin
 implicit none
 double precision, intent(out):: ekin_out
 ekin_out = ekin
end subroutine amuse_get_kinetic_energy

subroutine amuse_get_thermal_energy(etherm_out)
 use energies, only:etherm
 implicit none
 double precision, intent(out):: etherm_out
 etherm_out = etherm
end subroutine amuse_get_thermal_energy

subroutine amuse_get_time_step(dt_out)
 use timestep, only:dtmax
 implicit none
 double precision, intent(out):: dt_out
 dt_out = dtmax
end subroutine amuse_get_time_step

subroutine amuse_get_number_of_sph_particles(n)
 use part, only:npartoftype, igas
 implicit none
 integer, intent(out):: n
 logical:: nodisabled
 !nodisabled = .true.
 nodisabled = .false.
 n = npartoftype(igas)
end subroutine amuse_get_number_of_sph_particles

subroutine amuse_get_number_of_particles(n)
 implicit none
 integer, intent(out):: n
 logical:: nodisabled
 nodisabled = .true.
 call amuse_get_npart(n, nodisabled)
end subroutine amuse_get_number_of_particles

subroutine amuse_get_time(time_out)
 use timestep, only:time
 implicit none
 double precision, intent(out):: time_out
 time_out = time
end subroutine amuse_get_time

subroutine amuse_set_time_step(dt_in)
 use timestep, only:dtmax
 implicit none
 double precision, intent(in):: dt_in
 dtmax = dt_in
 print*, "dtmax: ", dtmax
end subroutine amuse_set_time_step

subroutine amuse_get_index(i, part_index)
 use io, only:fatal
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length), intent(out):: part_index
 ! This subroutine maps the unique index (i) to the current index (part_index)
 ! The map is synchronised after each evolve step-but it should also be synchronised after adding particles
 ! Note that i is strictly positive, negative indices are sink particles and they do not use this map!
 if (i /= abs(i)) call fatal('get_index', 'get_index is not for sink particles!')
 part_index = amuse_id_lookup(i)
end subroutine amuse_get_index

subroutine amuse_get_density(i, rho)
 use part, only:rhoh, iphase, massoftype, xyzh
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision:: pmassi
 double precision, intent(out):: rho
 call amuse_get_index(i, part_index)
 if (part_index == 0) then
    rho = 0
 else
    pmassi = massoftype(abs(iphase(part_index)))
    rho = rhoh(xyzh(4, part_index), pmassi)
 endif
end subroutine amuse_get_density

subroutine amuse_get_pressure(i, p)
 use part, only:rhoh, iphase, massoftype, xyzh
 use eos, only:ieos, equationofstate
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 integer:: eos_type
 double precision:: pmassi, ponrho, rho, spsound, x, y, z
 double precision, intent(out):: p
 real:: tempi
 call amuse_get_index(i, part_index)
 if (part_index == 0) then
    p = 0
 else
    eos_type = ieos
    pmassi = massoftype(abs(iphase(part_index)))
    call amuse_get_density(part_index, rho)
    x = xyzh(1, part_index)
    y = xyzh(2, part_index)
    z = xyzh(3, part_index)
    call equationofstate(eos_type, ponrho, spsound, rho, x, y, z, tempi)
    p = ponrho*rho
 endif
end subroutine amuse_get_pressure

subroutine amuse_get_mass(i, part_mass)
 use part, only:iphase, massoftype, xyzmh_ptmass
 implicit none
 double precision, intent(out):: part_mass
 integer(kind=index_length), intent(inout):: i
 integer(kind=index_length):: part_index
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    if (part_index == 0) then
       part_mass = 0
    else
       part_mass = massoftype(abs(iphase(part_index)))
    endif
 else
    part_mass = xyzmh_ptmass(4, -i)
 endif
end subroutine amuse_get_mass

subroutine amuse_get_state_gas(i, mass, x, y, z, vx, vy, vz, u, h)
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision, intent(inout):: mass, x, y, z, vx, vy, vz, u, h
 call amuse_get_mass(i, mass)
 call amuse_get_position(i, x, y, z)
 call amuse_get_velocity(i, vx, vy, vz)
 call amuse_get_internal_energy(i, u)
 call amuse_get_smoothing_length(i, h)
end subroutine amuse_get_state_gas

subroutine amuse_get_state_dm(i, mass, x, y, z, vx, vy, vz)
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision, intent(inout):: mass, x, y, z, vx, vy, vz
 call amuse_get_mass(i, mass)
 call amuse_get_position(i, x, y, z)
 call amuse_get_velocity(i, vx, vy, vz)
 write(*,*) 'getting dm ', i
end subroutine amuse_get_state_dm

subroutine amuse_get_state_sink(i, mass, x, y, z, vx, vy, vz, radius, accretion_radius)
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision, intent(inout):: mass, x, y, z, vx, vy, vz, radius, accretion_radius
 call amuse_get_mass(i, mass)
 call amuse_get_position(i, x, y, z)
 call amuse_get_velocity(i, vx, vy, vz)
 call amuse_get_sink_radius(i, radius)
 call amuse_get_sink_accretion_radius(i, accretion_radius)
end subroutine amuse_get_state_sink

subroutine amuse_set_state_gas(i, mass, x, y, z, vx, vy, vz, u, h)
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision:: mass, x, y, z, vx, vy, vz, u, h
 call amuse_set_mass(i, mass)
 call amuse_set_position(i, x, y, z)
 call amuse_set_velocity(i, vx, vy, vz)
 call amuse_set_internal_energy(i, u)
 call amuse_set_smoothing_length(i, h)
end subroutine amuse_set_state_gas

subroutine amuse_set_state_dm(i, mass, x, y, z, vx, vy, vz)
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision:: mass, x, y, z, vx, vy, vz
 call amuse_set_mass(i, mass)
 call amuse_set_position(i, x, y, z)
 call amuse_set_velocity(i, vx, vy, vz)
end subroutine amuse_set_state_dm

subroutine amuse_set_state_sink(i, mass, x, y, z, vx, vy, vz, radius, accretion_radius)
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision:: mass, x, y, z, vx, vy, vz, radius, accretion_radius
 call amuse_set_mass(i, mass)
 call amuse_set_position(i, x, y, z)
 call amuse_set_velocity(i, vx, vy, vz)
 call amuse_set_sink_radius(i, radius)
 call amuse_set_sink_accretion_radius(i, accretion_radius)
end subroutine amuse_set_state_sink

subroutine amuse_get_sink_radius(i, radius)
 implicit none
 integer(kind=index_length), intent(in):: i
 double precision:: radius
 call amuse_get_sink_effective_radius(i, radius)
end subroutine amuse_get_sink_radius

subroutine amuse_set_sink_radius(i, radius)
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision:: radius
 call amuse_set_sink_effective_radius(i, radius)
 call amuse_set_sink_accretion_radius(i, radius)
end subroutine amuse_set_sink_radius

subroutine amuse_get_sink_effective_radius(i, radius)
 use part, only:xyzmh_ptmass, iReff
 implicit none
 integer(kind=index_length), intent(in):: i
 double precision:: radius
 radius = xyzmh_ptmass(iReff, -i)
end subroutine amuse_get_sink_effective_radius

subroutine amuse_get_sink_accretion_radius(i, radius)
 use part, only:xyzmh_ptmass, ihacc
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision:: radius
 radius = xyzmh_ptmass(ihacc, -i)
end subroutine amuse_get_sink_accretion_radius

subroutine amuse_get_sink_temperature(i, temperature)
 use part, only:xyzmh_ptmass, iTeff
 implicit none
 integer(kind=index_length), intent(in):: i
 double precision:: temperature
 temperature = xyzmh_ptmass(iTeff, -i)
end subroutine amuse_get_sink_temperature

subroutine amuse_get_sink_luminosity(i, luminosity)
 use part, only:xyzmh_ptmass, iLum
 implicit none
 integer(kind=index_length), intent(in):: i
 double precision:: luminosity
 luminosity = xyzmh_ptmass(iLum, -i)
end subroutine amuse_get_sink_luminosity

subroutine amuse_get_position(i, x, y, z)
 use part, only:xyzh, xyzmh_ptmass
 implicit none
 integer(kind=index_length), intent(inout):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: x, y, z
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    if (part_index == 0) then
       x = 0
       y = 0
       z = 0
    else
       x = xyzh(1, part_index)
       y = xyzh(2, part_index)
       z = xyzh(3, part_index)
    endif
 else
    x = xyzmh_ptmass(1, -i)
    y = xyzmh_ptmass(2, -i)
    z = xyzmh_ptmass(3, -i)
 endif
end subroutine amuse_get_position

subroutine amuse_get_velocity(i, vx, vy, vz)
 use part, only:vxyzu, vxyz_ptmass
 implicit none
 integer(kind=index_length), intent(inout):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: vx, vy, vz
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    if (part_index == 0) then
       vx = 0
       vy = 0
       vz = 0
    else
       vx = vxyzu(1, part_index)
       vy = vxyzu(2, part_index)
       vz = vxyzu(3, part_index)
    endif
 else
    vx = vxyz_ptmass(1, -i)
    vy = vxyz_ptmass(2, -i)
    vz = vxyz_ptmass(3, -i)
 endif
end subroutine amuse_get_velocity

subroutine amuse_get_acceleration(i, fx, fy, fz)
 use part, only:fxyzu, fxyz_ptmass
 implicit none
 integer(kind=index_length), intent(inout):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: fx, fy, fz
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    if (part_index == 0) then
       fx = 0
       fy = 0
       fz = 0
    else
       fx = fxyzu(1, part_index)
       fy = fxyzu(2, part_index)
       fz = fxyzu(3, part_index)
    endif
 else
    fx = fxyz_ptmass(1, -i)
    fy = fxyz_ptmass(2, -i)
    fz = fxyz_ptmass(3, -i)
 endif
end subroutine amuse_get_acceleration

subroutine amuse_get_smoothing_length(i, h)
 use part, only:xyzh, xyzmh_ptmass, ihsoft
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: h
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    if (part_index == 0) then
       h = 0
    else
       h = xyzh(4, part_index)
    endif
 else
    h = xyzmh_ptmass(ihsoft, -i)
 endif
end subroutine amuse_get_smoothing_length

subroutine amuse_get_radius(i, radius)
 implicit none
 integer(kind = index_length), intent(in):: i
 double precision, intent(inout):: radius
 if (i == abs(i)) then
    call amuse_get_smoothing_length(i, radius)
 else
    call amuse_get_sink_radius(i, radius)
 endif
end subroutine amuse_get_radius

subroutine amuse_get_internal_energy(i, u)
 use dim, only:maxvxyzu
 use part, only:vxyzu
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: u
 call amuse_get_index(i, part_index)
 if (part_index == 0) then
    u = 0
 elseif (maxvxyzu >= 4) then
    u = vxyzu(4, part_index)
 else
    u = 0
 endif
end subroutine amuse_get_internal_energy

subroutine amuse_set_hi_abundance(i, hi_abundance)
 use part, only:abundance, iHI
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: hi_abundance
 call amuse_get_index(i, part_index)

 abundance(iHI, part_index) = hi_abundance
end subroutine amuse_set_hi_abundance

subroutine amuse_get_hi_abundance(i, hi_abundance)
 use part, only:abundance, iHI
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: hi_abundance
 call amuse_get_index(i, part_index)
 if (part_index == 0) then
    hi_abundance = 0
 else
    hi_abundance = abundance(iHI, part_index)
 endif
end subroutine amuse_get_hi_abundance

subroutine amuse_set_proton_abundance(i, proton_abundance)
 use part, only:abundance, iproton
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: proton_abundance
 call amuse_get_index(i, part_index)

 abundance(iproton, part_index) = proton_abundance
end subroutine amuse_set_proton_abundance

subroutine amuse_get_proton_abundance(i, proton_abundance)
 use part, only:abundance, iproton
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: proton_abundance
 call amuse_get_index(i, part_index)
 if (part_index == 0) then
    proton_abundance = 0
 else
    proton_abundance = abundance(iproton, part_index)
 endif
end subroutine amuse_get_proton_abundance

subroutine amuse_set_electron_abundance(i, electron_abundance)
 use part, only:abundance, ielectron
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: electron_abundance
 call amuse_get_index(i, part_index)

 abundance(ielectron, part_index) = electron_abundance
end subroutine amuse_set_electron_abundance

subroutine amuse_get_electron_abundance(i, electron_abundance)
 use part, only:abundance, ielectron
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: electron_abundance
 call amuse_get_index(i, part_index)
 if (part_index == 0) then
    electron_abundance = 0
 else
    electron_abundance = abundance(ielectron, part_index)
 endif
end subroutine amuse_get_electron_abundance

subroutine amuse_set_co_abundance(i, co_abundance)
 use part, only:abundance, ico
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: co_abundance
 call amuse_get_index(i, part_index)

 abundance(ico, part_index) = co_abundance
end subroutine amuse_set_co_abundance

subroutine amuse_get_co_abundance(i, co_abundance)
 use part, only:abundance, ico
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: co_abundance
 call amuse_get_index(i, part_index)
 if (part_index == 0) then
    co_abundance = 0
 else
    co_abundance = abundance(ico, part_index)
 endif
end subroutine amuse_get_co_abundance

subroutine amuse_set_h2ratio(i, h2ratio)
 use part, only:abundance, iHI, ih2ratio
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: h2ratio
 call amuse_get_index(i, part_index)

 abundance(ih2ratio, part_index) = h2ratio
end subroutine amuse_set_h2ratio

subroutine amuse_get_h2ratio(i, h2ratio)
 use part, only:abundance, iHI, ih2ratio
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(out):: h2ratio
 call amuse_get_index(i, part_index)
 if (part_index == 0) then
    h2ratio = 0
 else
    h2ratio = abundance(ih2ratio, part_index)
 endif
end subroutine amuse_get_h2ratio

subroutine amuse_set_mass(i, part_mass)
 use part, only:iphase, massoftype, xyzmh_ptmass
 implicit none
 double precision, intent(in):: part_mass
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    massoftype(abs(iphase(part_index))) = part_mass
 else
    xyzmh_ptmass(4, -i) = part_mass
 endif
end subroutine amuse_set_mass

subroutine amuse_set_sink_accretion_radius(i, radius)
 use part, only:xyzmh_ptmass, ihacc
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision:: radius
 xyzmh_ptmass(ihacc, i) = radius
end subroutine amuse_set_sink_accretion_radius

subroutine amuse_set_sink_effective_radius(i, radius)
 use part, only:xyzmh_ptmass, iReff
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision:: radius
 xyzmh_ptmass(iReff, i) = radius
end subroutine amuse_set_sink_effective_radius

subroutine amuse_set_position(i, x, y, z)
 use part, only:xyzh, xyzmh_ptmass
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: x, y, z
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    xyzh(1, part_index) = x
    xyzh(2, part_index) = y
    xyzh(3, part_index) = z
 else
    xyzmh_ptmass(1, -i) = x
    xyzmh_ptmass(2, -i) = y
    xyzmh_ptmass(3, -i) = z
 endif
end subroutine amuse_set_position

subroutine amuse_set_velocity(i, vx, vy, vz)
 use part, only:vxyzu, vxyz_ptmass
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: vx, vy, vz
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    vxyzu(1, part_index) = vx
    vxyzu(2, part_index) = vy
    vxyzu(3, part_index) = vz
 else
    vxyz_ptmass(1, -i) = vx
    vxyz_ptmass(2, -i) = vy
    vxyz_ptmass(3, -i) = vz
 endif
end subroutine amuse_set_velocity

subroutine amuse_set_smoothing_length(i, h)
 use part, only:xyzh, xyzmh_ptmass, ihsoft
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: h
 if (i == abs(i)) then
    call amuse_get_index(i, part_index)
    xyzh(4, part_index) = h
 else
    xyzmh_ptmass(ihsoft, -i) = h
 endif
end subroutine amuse_set_smoothing_length

subroutine amuse_set_radius(i, radius)
 implicit none
 integer(kind=index_length), intent(inout):: i
 double precision:: radius
 if (i == abs(i)) then
    call amuse_set_smoothing_length(i, radius)
 else
    call amuse_set_sink_radius(i, radius)
 endif
end subroutine amuse_set_radius

subroutine amuse_set_sink_temperature(i, temperature)
 use part, only:xyzmh_ptmass, iTeff
 implicit none
 integer(kind=index_length), intent(in):: i
 double precision:: temperature
 xyzmh_ptmass(iTeff, -i) = temperature
end subroutine amuse_set_sink_temperature

subroutine amuse_set_sink_luminosity(i, luminosity)
 use part, only:xyzmh_ptmass, iLum
 implicit none
 integer(kind=index_length), intent(in):: i
 double precision:: luminosity
 xyzmh_ptmass(iLum, -i) = luminosity
end subroutine amuse_set_sink_luminosity

subroutine amuse_set_internal_energy(i, u)
 use dim, only:maxvxyzu
 use part, only:vxyzu
 use timestep, only:dtextforce
 implicit none
 integer(kind=index_length), intent(in):: i
 integer(kind=index_length):: part_index
 double precision, intent(in):: u
 call amuse_get_index(i, part_index)
 if (maxvxyzu >= 4) then
    vxyzu(4, part_index) = u
 endif
 ! Changing temperature -> better use a small cooling step
 dtextforce = 1.e-8
end subroutine amuse_set_internal_energy

subroutine amuse_evolve_model(tmax_in)
 use timestep, only:tmax, time, dtmax
 !use timestep, only:dt, dtlast
 !use timestep, only:rhomaxnow
 use evolve, only:evol
#ifdef IND_TIMESTEPS
 use timestep_ind, only:istepfrac
#endif
#ifdef INJECT_PARTICLES
 use inject,           only:inject_particles
 use partinject,       only:update_injected_particles
#endif
 use part, only:norig
 use step_lf_global, only:init_step

 implicit none
 character(len = 120):: infile, logfile, evfile, dumpfile
 integer (kind = index_length):: number_of_particles_at_start
 integer (kind = index_length):: number_of_particles_at_finish
 integer (kind = index_length):: norig_amuse
 logical:: amuse_initialise
 double precision, intent(in):: tmax_in
 real:: tlast
 real:: dtinject
 integer(kind = 1):: nbinmax
#ifndef IND_TIMESTEPS
 integer:: istepfrac
 istepfrac = 0  ! dummy values
#endif
 infile = ""
 logfile = ""
 evfile = ""
 dumpfile = ""
 dtinject  = huge(dtinject)
 ! dtlast = 0
 nbinmax = 0
 ! This is a workaround for the fact that norig is an integer(kind=8)
 if (norig < min_int32 .or. norig > max_int32) then
    error stop "norig out of range for AMUSE"
 endif
 norig_amuse = int(norig, kind=index_length)
 number_of_particles_at_start = norig_amuse

 tmax = tmax_in  ! - epsilon(tmax_in)
 !dtmax = (tmax-time)

 tlast = time
 write(*,*) "TIMESTEPPING: evolve from ", time, " to ", tmax
 ! The reason for doing timestepping here (rather than just in evol) is that we want to be able to use AMUSE stopping conditions,
 ! such as high density detection.
 ! Allowing for a shortage of 1% of dtmax to account for floating point differences
 timestepping: do while (time+0.01*dtmax < tmax)

#ifdef IND_TIMESTEPS
    istepfrac = 0
#endif
    amuse_initialise = .false.
    call evol(infile, logfile, evfile, dumpfile)
    ! Check for stopping conditions here
 enddo timestepping

 call construct_id_lookup()

 ! This is a workaround for the fact that norig is an integer(kind=8)
 if (norig < min_int32 .or. norig > max_int32) then
    error stop "norig out of range for AMUSE"
 endif
 norig_amuse = int(norig, kind=index_length)
 number_of_particles_at_finish = norig_amuse
 new_particles_since_last_update = new_particles_since_last_update+number_of_particles_at_finish-number_of_particles_at_start
end subroutine amuse_evolve_model

!
! Setters and getters for parameters
!

! Setters

subroutine amuse_set_c_courant(C_cour_in)
 use timestep, only:C_cour
 implicit none
 double precision, intent(in):: C_cour_in
 C_cour = C_cour_in
end subroutine amuse_set_c_courant

subroutine amuse_set_c_force(C_force_in)
 use timestep, only:C_force
 implicit none
 double precision, intent(in):: C_force_in
 C_force = C_force_in
end subroutine amuse_set_c_force

subroutine amuse_set_c_cool(C_cool_in)
 use timestep, only:C_cool
 implicit none
 double precision, intent(in):: C_cool_in
 C_cool = C_cool_in
end subroutine amuse_set_c_cool

subroutine amuse_set_tolv(tolv_in)
 use timestep, only:tolv
 implicit none
 double precision, intent(in):: tolv_in
 tolv = tolv_in
end subroutine amuse_set_tolv

subroutine amuse_set_hfact(hfact_in)
 use part, only:hfact
 implicit none
 double precision, intent(in):: hfact_in
 hfact = hfact_in
end subroutine amuse_set_hfact

subroutine amuse_set_tolh(tolh_in)
 use options, only:tolh
 implicit none
 double precision, intent(in):: tolh_in
 tolh = tolh_in
end subroutine amuse_set_tolh

subroutine amuse_set_tree_accuracy(tree_accuracy_in)
 use kdtree, only:tree_accuracy
 implicit none
 double precision, intent(in):: tree_accuracy_in
 tree_accuracy = tree_accuracy_in
end subroutine amuse_set_tree_accuracy

subroutine amuse_set_alpha(alpha_in)
 use options, only:alpha
 implicit none
 double precision, intent(in):: alpha_in
 alpha = alpha_in
end subroutine amuse_set_alpha

subroutine amuse_set_alphamax(alphamax_in)
 use options, only:alphamax
 implicit none
 double precision, intent(in):: alphamax_in
 alphamax = alphamax_in
end subroutine amuse_set_alphamax

subroutine amuse_set_beta(beta_in)
 use options, only:beta
 implicit none
 double precision, intent(in):: beta_in
 beta = beta_in
end subroutine amuse_set_beta

subroutine amuse_set_avdecayconst(avdecayconst_in)
 use options, only:avdecayconst
 implicit none
 double precision, intent(in):: avdecayconst_in
 avdecayconst = avdecayconst_in
end subroutine amuse_set_avdecayconst

subroutine amuse_set_idamp(idamp_in)
 use options, only:idamp
 implicit none
 integer, intent(in):: idamp_in
 idamp = idamp_in
end subroutine amuse_set_idamp

subroutine amuse_set_ieos(ieos_in)
 use eos, only:ieos
 implicit none
 integer, intent(in):: ieos_in
 ieos = ieos_in
end subroutine amuse_set_ieos

subroutine amuse_set_icooling(icooling_in)
 use io, only:id, master, iprint
 use options, only:icooling
 !use options, only:iexternalforce
 !use dim, only:h2chemistry
 !use chem, only:init_chem
 use cooling, only:init_cooling, Tfloor
 !use cooling_ism, only:init_cooling
 implicit none
 integer:: ierr
 integer, intent(in):: icooling_in
 icooling = icooling_in
 if (icooling > 0) then
    Tfloor = 1  ! K
    call init_cooling(id, master, iprint, ierr)
    !if (h2chemistry) then
    !    if (id == master) write(iprint, *) 'initialising cooling function...'
    !    call init_chem()
    !    call init_h2cooling()
    !else
    !    call init_cooling(ierr)
    !    if (ierr /= 0) call fatal('initial','error initialising cooling')
    !endif
 endif
end subroutine amuse_set_icooling

subroutine amuse_set_polyk(polyk_in)
 use eos, only:polyk
 implicit none
 double precision, intent(in):: polyk_in
 polyk = polyk_in
end subroutine amuse_set_polyk

subroutine amuse_set_mu(mu_in)
 use eos, only:gmw
 implicit none
 double precision, intent(in):: mu_in
 gmw = mu_in
end subroutine amuse_set_mu

subroutine amuse_set_rhofinal(rhofinal_in)
 use options, only:rhofinal_cgs, rhofinal1
 use units, only:unit_density
 implicit none
 double precision, intent(in):: rhofinal_in
 rhofinal_cgs = rhofinal_in*unit_density
 if (rhofinal_cgs > 0.) then
    rhofinal1 = unit_density/rhofinal_cgs
 else
    rhofinal1 = 0.0
 endif
end subroutine amuse_set_rhofinal

subroutine amuse_set_rho_crit(rho_crit_in)
 use units, only:unit_density
 use ptmass, only:rho_crit_cgs
 implicit none
 double precision, intent(in):: rho_crit_in
 rho_crit_cgs = rho_crit_in*unit_density
end subroutine amuse_set_rho_crit

subroutine amuse_set_r_crit(r_crit_in)
 use ptmass, only:r_crit
 implicit none
 double precision, intent(in):: r_crit_in
 r_crit = r_crit_in
end subroutine amuse_set_r_crit

subroutine amuse_set_h_acc(h_acc_in)
 use ptmass, only:h_acc
 implicit none
 double precision, intent(in):: h_acc_in
 h_acc = h_acc_in
end subroutine amuse_set_h_acc

subroutine amuse_set_h_soft_sinkgas(h_soft_sinkgas_in)
 use ptmass, only:h_soft_sinkgas
 implicit none
 double precision, intent(in):: h_soft_sinkgas_in
 h_soft_sinkgas = h_soft_sinkgas_in
end subroutine amuse_set_h_soft_sinkgas

subroutine amuse_set_h_soft_sinksink(h_soft_sinksink_in)
 use ptmass, only:h_soft_sinksink
 implicit none
 double precision, intent(in):: h_soft_sinksink_in
 h_soft_sinksink = h_soft_sinksink_in
end subroutine amuse_set_h_soft_sinksink

subroutine amuse_set_f_acc(f_acc_in)
 use ptmass, only:f_acc
 implicit none
 double precision, intent(in):: f_acc_in
 f_acc = f_acc_in
end subroutine amuse_set_f_acc

subroutine amuse_set_iexternalforce(iexternalforce_in)
 use options, only:iexternalforce
 implicit none
 integer, intent(in):: iexternalforce_in
 iexternalforce = iexternalforce_in
end subroutine amuse_set_iexternalforce

subroutine amuse_set_irealvisc(irealvisc_in)
 use viscosity, only:irealvisc
 implicit none
 integer, intent(in):: irealvisc_in
 irealvisc = irealvisc_in
end subroutine amuse_set_irealvisc

subroutine amuse_set_shearparam(shearparam_in)
 use viscosity, only:shearparam
 implicit none
 double precision, intent(in):: shearparam_in
 shearparam = shearparam_in
end subroutine amuse_set_shearparam

subroutine amuse_set_bulkvisc(bulkvisc_in)
 use viscosity, only:bulkvisc
 implicit none
 double precision, intent(in):: bulkvisc_in
 bulkvisc = bulkvisc_in
end subroutine amuse_set_bulkvisc

subroutine amuse_set_gamma(gamma_in)
 use eos, only:gamma
 implicit none
 double precision, intent(in):: gamma_in
 gamma = gamma_in
end subroutine amuse_set_gamma

subroutine amuse_set_umass(umass_in)
 use units, only:umass
 implicit none
 double precision, intent(in):: umass_in
 umass = umass_in
end subroutine amuse_set_umass

subroutine amuse_set_udist(udist_in)
 use units, only:udist
 implicit none
 double precision, intent(in):: udist_in
 udist = udist_in
end subroutine amuse_set_udist

subroutine amuse_set_utime(utime_in)
 use units, only:utime
 implicit none
 double precision, intent(in):: utime_in
 utime = utime_in
end subroutine amuse_set_utime

! End of Setters

! Getters

subroutine amuse_get_c_courant(C_cour_out)
 use timestep, only:C_cour
 implicit none
 double precision, intent(out):: C_cour_out
 C_cour_out = C_cour
end subroutine amuse_get_c_courant

subroutine amuse_get_c_force(C_force_out)
 use timestep, only:C_force
 implicit none
 double precision, intent(out):: C_force_out
 C_force_out = C_force
end subroutine amuse_get_c_force

subroutine amuse_get_c_cool(C_cool_out)
 use timestep, only:C_cool
 implicit none
 double precision, intent(out):: C_cool_out
 C_cool_out = C_cool
end subroutine amuse_get_c_cool

subroutine amuse_get_tolv(tolv_out)
 use timestep, only:tolv
 implicit none
 double precision, intent(out):: tolv_out
 tolv_out = tolv
end subroutine amuse_get_tolv

subroutine amuse_get_hfact(hfact_out)
 use part, only:hfact
 implicit none
 double precision, intent(out):: hfact_out
 hfact_out = hfact
end subroutine amuse_get_hfact

subroutine amuse_get_tolh(tolh_out)
 use options, only:tolh
 implicit none
 double precision, intent(out):: tolh_out
 tolh_out = tolh
end subroutine amuse_get_tolh

subroutine amuse_get_tree_accuracy(tree_accuracy_out)
 use kdtree, only:tree_accuracy
 implicit none
 double precision, intent(out):: tree_accuracy_out
 tree_accuracy_out = tree_accuracy
end subroutine amuse_get_tree_accuracy

subroutine amuse_get_alpha(alpha_out)
 use options, only:alpha
 implicit none
 double precision, intent(out):: alpha_out
 alpha_out = alpha
end subroutine amuse_get_alpha

subroutine amuse_get_alphamax(alphamax_out)
 use options, only:alphamax
 implicit none
 double precision, intent(out):: alphamax_out
 alphamax_out = alphamax
end subroutine amuse_get_alphamax

subroutine amuse_get_beta(beta_out)
 use options, only:beta
 implicit none
 double precision, intent(out):: beta_out
 beta_out = beta
end subroutine amuse_get_beta

subroutine amuse_get_avdecayconst(avdecayconst_out)
 use options, only:avdecayconst
 implicit none
 double precision, intent(out):: avdecayconst_out
 avdecayconst_out = avdecayconst
end subroutine amuse_get_avdecayconst

subroutine amuse_get_idamp(idamp_out)
 use options, only:idamp
 implicit none
 integer, intent(out):: idamp_out
 idamp_out = idamp
end subroutine amuse_get_idamp

subroutine amuse_get_ieos(ieos_out)
 use eos, only:ieos
 implicit none
 integer, intent(out):: ieos_out
 ieos_out = ieos
end subroutine amuse_get_ieos

subroutine amuse_get_icooling(icooling_out)
 use options, only:icooling
 implicit none
 integer, intent(out):: icooling_out
 icooling_out = icooling
end subroutine amuse_get_icooling

subroutine amuse_get_polyk(polyk_out)
 use eos, only:polyk
 implicit none
 double precision, intent(out):: polyk_out
 polyk_out = polyk
end subroutine amuse_get_polyk

subroutine amuse_get_mu(mu_out)
 use eos, only:gmw
 implicit none
 double precision, intent(out):: mu_out
 mu_out = gmw
end subroutine amuse_get_mu

subroutine amuse_get_rhofinal(rhofinal_out)
 use options, only:rhofinal_cgs
 use units, only:unit_density
 implicit none
 double precision, intent(out):: rhofinal_out
 rhofinal_out = rhofinal_cgs/unit_density
end subroutine amuse_get_rhofinal

subroutine amuse_get_rho_crit(rho_crit_out)
 use units, only:unit_density
 use ptmass, only:rho_crit_cgs
 implicit none
 double precision, intent(out):: rho_crit_out
 rho_crit_out = rho_crit_cgs/unit_density
end subroutine amuse_get_rho_crit

subroutine amuse_get_r_crit(r_crit_out)
 use ptmass, only:r_crit
 implicit none
 double precision, intent(out):: r_crit_out
 r_crit_out = r_crit
end subroutine amuse_get_r_crit

subroutine amuse_get_h_acc(h_acc_out)
 use ptmass, only:h_acc
 implicit none
 double precision, intent(out):: h_acc_out
 h_acc_out = h_acc
end subroutine amuse_get_h_acc

subroutine amuse_get_h_soft_sinkgas(h_soft_sinkgas_out)
 use ptmass, only:h_soft_sinkgas
 implicit none
 double precision, intent(out):: h_soft_sinkgas_out
 h_soft_sinkgas_out = h_soft_sinkgas
end subroutine amuse_get_h_soft_sinkgas

subroutine amuse_get_h_soft_sinksink(h_soft_sinksink_out)
 use ptmass, only:h_soft_sinksink
 implicit none
 double precision, intent(out):: h_soft_sinksink_out
 h_soft_sinksink_out = h_soft_sinksink
end subroutine amuse_get_h_soft_sinksink

subroutine amuse_get_f_acc(f_acc_out)
 use ptmass, only:f_acc
 implicit none
 double precision, intent(out):: f_acc_out
 f_acc_out = f_acc
end subroutine amuse_get_f_acc

subroutine amuse_get_iexternalforce(iexternalforce_out)
 use options, only:iexternalforce
 implicit none
 integer, intent(out):: iexternalforce_out
 iexternalforce_out = iexternalforce
end subroutine amuse_get_iexternalforce

subroutine amuse_get_irealvisc(irealvisc_out)
 use viscosity, only:irealvisc
 implicit none
 integer, intent(out):: irealvisc_out
 irealvisc_out = irealvisc
end subroutine amuse_get_irealvisc

subroutine amuse_get_shearparam(shearparam_out)
 use viscosity, only:shearparam
 implicit none
 double precision, intent(out):: shearparam_out
 shearparam_out = shearparam
end subroutine amuse_get_shearparam

subroutine amuse_get_bulkvisc(bulkvisc_out)
 use viscosity, only:bulkvisc
 implicit none
 double precision, intent(out):: bulkvisc_out
 bulkvisc_out = bulkvisc
end subroutine amuse_get_bulkvisc

subroutine amuse_get_gamma(gamma_out)
 use eos, only:gamma
 implicit none
 double precision, intent(out):: gamma_out
 gamma_out = gamma
end subroutine amuse_get_gamma

subroutine amuse_get_umass(umass_out)
 use units, only:umass
 implicit none
 double precision, intent(out):: umass_out
 umass_out = umass
end subroutine amuse_get_umass

subroutine amuse_get_utime(utime_out)
 use units, only:utime
 implicit none
 double precision, intent(out):: utime_out
 utime_out = utime
end subroutine amuse_get_utime

subroutine amuse_get_udist(udist_out)
 use units, only:udist
 implicit none
 double precision, intent(out):: udist_out
 udist_out = udist
end subroutine amuse_get_udist

! End of Getters

end module AmusePhantom
