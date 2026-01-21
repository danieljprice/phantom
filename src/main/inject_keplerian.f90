!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection of material at keplerian speed in an accretion disc
!
! :References:
!
! :Owner: Cristiano Longarini
!
! :Runtime parameters:
!   - HonR_inj    : *aspect ratio to give temperature at rinj*
!   - follow_sink : *injection radius is relative to sink particle 1*
!   - incx        : *inclination about the x-axis (degrees)*
!   - incy        : *inclination about the x-axis (degrees)*
!   - incz        : *inclination about the x-axis (degrees)*
!   - mdot        : *mass injection rate [msun/yr]*
!   - rinj        : *injection radius*
!   - sigma_inj   : *width of gaussian injection profile, =0 is ring injection*
!
! :Dependencies: eos, externalforces, infile_utils, io, options, part,
!   partinject, physcon, random, units, vectorutils
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'keplerian'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject,update_injected_par

 real :: mdot = 0.
 real :: rinj = 25.
 real :: HonR_inj = 0.05
 real :: incx = 0.0
 real :: incy = 0.0
 real :: incz = 0.0
 real :: sigma_inj = 0.0
 logical :: follow_sink = .true.
 integer, private :: iseed=-888

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,   only:warning
 use part, only:nptmass
 integer,  intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0
 if (nptmass > 1) call warning(inject_type,'Using first sink particle to compute Keplerian velocity')

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  set defaults
!+
!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling injection at a given radius rinj
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,             only:fatal,iverbose,warning
 use part,           only:massoftype,igas,nptmass,isdead_or_accreted,maxvxyzu
 use partinject,     only:add_or_update_particle
 use physcon,        only:pi,solarm,years
 use units,          only:umass,utime
 use random,         only:ran2,gauss_random
 use options,        only:iexternalforce,ieos
 use externalforces, only:mass1
 use eos,            only:equationofstate,gamma
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real :: Minject,Mdot_code
 real :: frac_extra,deltat
 real :: x0(3),v0(3),mstar,r2min,dr2,hguess,phi,cosphi,sinphi,r2,xyzi(3),vxyz(3),u,xyzit(3),vxyzt(3)
 real :: vkep,vphi,zi,cs,bigH,gaussian_seed
 real :: dum_ponrho,dum_rho,dum_temp
 integer :: i,k,i_part,ninject
 !
 ! convert mass loss rate from Msun/yr to code units
 !
 Mdot_code = Mdot*(solarm/umass)*(utime/years)

 !
 ! get central mass
 !
 x0 = 0.
 v0 = 0.
 if (iexternalforce > 0) then
    mstar = mass1
 elseif (nptmass >= 1) then
    if (follow_sink) then
       x0 = xyzmh_ptmass(1:3,1)
       v0 = vxyz_ptmass(1:3,1)
    endif
    mstar = xyzmh_ptmass(4,1)
 else
    mstar = 1.
    call fatal(inject_type,'no central object to compute Keplerian velocity')
 endif

 ! for the smoothing length, take it from the closest existing particle to the injection radius
 hguess = 1.
 r2min = huge(r2min)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       r2 = (xyzh(1,i)-x0(1))**2 + (xyzh(2,i)-x0(2))**2
       dr2 = abs(r2 - rinj*rinj)
       if (dr2 < r2min) then
          hguess = xyzh(4,i)
          r2min = dr2
       endif
    endif
 enddo

 vkep = sqrt(mstar/rinj)

 ! for the temperature, call equation of state to get cs at this radius
 if (maxvxyzu >= 4) then
    ! use HonR parameter
    cs = HonR_inj * vkep
 else
    dum_rho = 1.
    dum_temp = 0.
    if (gamma > 1.001) then
       call warning(inject_type,'cannot get temp at r=rinj without knowing density, injecting at z=0')
       cs = 0.
    else
       call equationofstate(ieos,dum_ponrho,cs,dum_rho,rinj,0.,0.,dum_temp)
    endif
 endif

 !
 ! calculate how much mass to inject based on
 ! time interval since last injection
 !
 deltat    = dtlast
 Minject   = Mdot_code*deltat
 !
 ! work out number of particles by divide by mass of gas particles
 !
 ninject = int(Minject/massoftype(igas))
 !
 ! for the residual, roll the dice
 !
 frac_extra = Minject/massoftype(igas) - 2*(ninject/2)
 if (ran2(iseed) < 0.5*frac_extra) ninject = ninject + 2

 if (iverbose >= 2) print*,' injecting ',&
                    ninject,Minject/massoftype(igas),massoftype(igas)

 if (ninject > 0) then
    do k=1,ninject/2
       !
       ! get random position on ring
       !
       phi = 2.*pi*(ran2(iseed) - 0.5)

       cosphi = cos(phi)
       sinphi = sin(phi)

       gaussian_seed = gauss_random(iseed) * sigma_inj

       ! fix sound speed for gaussian injection
       if (sigma_inj > 0.) then
          call equationofstate(ieos,dum_ponrho,cs,dum_rho,abs(rinj+gaussian_seed),0.,0.,dum_temp)
       endif

       bigH = cs*abs(rinj+gaussian_seed)/(vkep*sqrt(rinj/abs(rinj+gaussian_seed)))
       zi = gauss_random(iseed)*bigH

       vphi = vkep*(1. - (zi/rinj)**2)**(-0.75)*sqrt(rinj/abs(rinj+gaussian_seed))  ! see Martire et al. (2024)

       xyzi = (/abs(rinj+gaussian_seed)*cosphi,abs(rinj+gaussian_seed)*sinphi,zi/)
       vxyz = (/-vphi*sinphi, vphi*cosphi, 0./)

       u = 1.5*cs**2

       xyzit = xyzi+x0
       vxyzt = vxyz+v0

       call apply_injection_rotations(xyzit,vxyzt)

       i_part = npart + 1! all particles are new
       call add_or_update_particle(igas, xyzit, vxyzt, hguess, u, i_part, npart, npartoftype, xyzh, vxyzu)

       xyzit = -xyzi+x0
       vxyzt = -vxyz+v0

       call apply_injection_rotations(xyzit,vxyzt)
       i_part = npart + 1! all particles are new
       call add_or_update_particle(igas, xyzit, vxyzt, hguess, u, i_part, npart, npartoftype, xyzh, vxyzu)
    enddo
 endif

 if (iverbose >= 2) then
    print*,'npart = ',npart
 endif
 !
 !-- no constraint on timestep
 !
 dtinject = huge(dtinject)

contains

subroutine apply_injection_rotations(xyz,vxyz)
 use vectorutils,    only:rotatevec
 use physcon,        only:deg_to_rad
 real, intent(inout) :: xyz(3), vxyz(3)

 call rotatevec(xyz,(/1.,0.,0./),incx*deg_to_rad)
 call rotatevec(vxyz,(/1.,0.,0./),incx*deg_to_rad)

 call rotatevec(xyz,(/0.,1.,0./),incy*deg_to_rad)
 call rotatevec(vxyz,(/0.,1.,0./),incy*deg_to_rad)

 call rotatevec(xyz,(/0.,0.,1./),incz*deg_to_rad)
 call rotatevec(vxyz,(/0.,0.,1./),incz*deg_to_rad)

end subroutine apply_injection_rotations

end subroutine inject_particles

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 use part,         only:maxvxyzu,nptmass
 integer, intent(in) :: iunit

 call write_inopt(mdot,'mdot','mass injection rate [msun/yr]',iunit)
 call write_inopt(rinj,'rinj','injection radius',iunit)
 call write_inopt(sigma_inj,'sigma_inj','width of gaussian injection profile, =0 is ring injection',iunit)
 if (maxvxyzu >= 4) call write_inopt(HonR_inj,'HonR_inj','aspect ratio to give temperature at rinj',iunit)
 call write_inopt(incx,'incx','inclination about the x-axis (degrees)',iunit)
 call write_inopt(incy,'incy','inclination about the x-axis (degrees)',iunit)
 call write_inopt(incz,'incz','inclination about the x-axis (degrees)',iunit)
 if (nptmass >= 1) call write_inopt(follow_sink,'follow_sink','injection radius is relative to sink particle 1',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 use part,         only:maxvxyzu,nptmass
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(mdot,'mdot',db,errcount=nerr,min=0.)
 call read_inopt(rinj,'rinj',db,errcount=nerr,min=0.)
 call read_inopt(sigma_inj,'sigma_inj',db,errcount=nerr,min=0.)
 if (maxvxyzu >= 4) call read_inopt(HonR_inj,'HonR_inj',db,errcount=nerr,min=0.)
 call read_inopt(incx,'incx',db,errcount=nerr,min=0.)
 call read_inopt(incy,'incy',db,errcount=nerr,min=0.)
 call read_inopt(incz,'incz',db,errcount=nerr,min=0.)
 if (nptmass >= 1) call read_inopt(follow_sink,'follow_sink',db,errcount=nerr,default=.true.)

end subroutine read_options_inject

end module inject
