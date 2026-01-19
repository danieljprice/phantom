!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection of material on a parabolic orbit like a streamer
!
! :References:
!   Longarini et al., hopefully
!
! :Owner: Cristiano Longarini
!
! :Runtime parameters:
!   - Rimp_streamer : *impact radius on disc*
!   - Rin_streamer  : *injection radius*
!   - Rp_streamer   : *pericentre distance*
!   - Win_streamer  : *streamer cross-section at injection*
!   - incl_streamer : *inclination at impact [deg]*
!   - ingoing       : *TRUE=pre-pericentre*
!   - mdot_streamer : *mass injection rate [Msun/yr]*
!   - phi_streamer  : *node longitude [deg]*
!
! :Dependencies: dim, eos, externalforces, infile_utils, io, options, part,
!   partinject, physcon, random, set_streamer, units
!
 use physcon,        only: pi, solarm, years
 use dim,            only: use_dust
 use units,          only: umass, utime
 use io,             only: fatal, warning, iverbose
 use options,        only: iexternalforce, ieos, use_dustfrac
 use externalforces, only:mass1
 use part,           only: igas, massoftype, nptmass, isdead_or_accreted, dustfrac
 use partinject,     only: add_or_update_particle
 use eos,            only: equationofstate, gamma
 use random,         only: ran2
 use set_streamer,   only: set_streamer_particle

 implicit none

 character(len=*), parameter, public :: inject_type = 'streamer'

 public :: init_inject, inject_particles,&
           write_options_inject, read_options_inject,&
           set_default_options_inject, update_injected_par

 real    :: mdot_streamer = 0.0
 real    :: Rp_streamer   = 1.0
 real    :: Rin_streamer  = 25.0
 real    :: Rimp_streamer = 10.0
 real    :: incl_streamer = 30.0
 real    :: phi_streamer  = 0.0
 real    :: Win_streamer  = 0.5
 logical :: ingoing       = .true.
 integer, private :: iseed = -987654

contains

!-----------------------------------------------------------------------
!  Initialize global variables or arrays needed for injection routine
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,      only: fatal, warning
 use options, only:iexternalforce
 use part,    only: nptmass
 integer, intent(out) :: ierr
 ierr = 0
 if (nptmass < 1 .and. iexternalforce <= 0) then
    call fatal(inject_type,'need a central mass (sink or externalforces) to compute mu')
 endif
end subroutine init_inject

!-----------------------------------------------------------------------
!  Set defaults
!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
! Main routine: inject new particles inside a ring orthogonal to the motion
! at (R_in), centered on the streamline position. Velocity is identical
! for all injected particles and equals the streamline velocity.
!-----------------------------------------------------------------------
subroutine inject_particles(time, dtlast, xyzh, vxyzu, &
                                       xyzmh_ptmass, vxyz_ptmass, &
                                       npart, npart_old, npartoftype, dtinject)
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 real :: mu, mstar
 real :: x0(3), v0(3), x(3), v(3)
 real :: tvec(3), ref(3), tmp(3), n1(3), n2(3)
 real :: Mdot_code, Minject, deltat
 integer :: Nin, k, i_part
 real :: phi, cs, u, R_to_sink
 real :: sink_pos(3)
 real :: dum_ponrho, dum_rho, dum_temp
 real :: hguess, r_random

 ! gravitational parameter mu from sink mass (G=1)
 if (iexternalforce > 0) then
    mstar = mass1
 else
    mstar = xyzmh_ptmass(4,1)
 endif
 mu = mstar

 ! center of streamer: do NOT follow sink position (use origin)
 ! sink may move , but the injection doesn't follow it
 ! possibility to change this - TBD
 ! set_streamer position and velocity at the centre of the streamline
 call set_streamer_particle(mu, Rp_streamer, Rin_streamer, Rimp_streamer, &
                incl_streamer, x0, v0, i_part, phi_streamer, ingoing)
 if (i_part /= 0) then
    call fatal(inject_type,'set_streamer_particle failed (bad geometry or inputs)')
 endif

 ! orthonormal basis in the plane perpendicular to velocity
 ! needed to set up the firehose
 tvec = v0 / max(1e-30, sqrt(dot_product(v0,v0)))
 ref  = (/0.0, 0.0, 1.0/)
 if (abs(dot_product(tvec,ref)) > 0.95) then
    ! security check if vr parallel to z
    ref = x0 / max(1e-30, sqrt(dot_product(x0,x0)))
 endif
 call cross_product(tvec, ref, tmp)
 n1 = tmp / max(1e-30, sqrt(dot_product(tmp,tmp)))
 call cross_product(tvec, n1, tmp)
 n2 = tmp / max(1e-30, sqrt(dot_product(tmp,tmp)))

 Mdot_code = mdot_streamer * (solarm/umass) * (utime/years)
 deltat    = dtlast
 Minject   = Mdot_code * deltat

 Nin = int( Minject / massoftype(igas) )
 ! roll the dice
 if (ran2(iseed) < (Minject/massoftype(igas) - real(Nin))) Nin = Nin + 1

 if (Nin <= 0) then
    dtinject = huge(dtinject)
    return
 endif

 ! smoothing length choice requested: h = 0.25 * W_in
 ! don't know if it is correct, should not matter much
 hguess = 0.25 * Win_streamer

 ! sink position for EOS
 ! now it works only for locally isothermal disc
 if (nptmass >= 1) then
    sink_pos = xyzmh_ptmass(1:3,1)
 else
    sink_pos = 0.0 ! potential
 endif

 do k = 1, Nin
    ! generate random radius with uniform distribution inside circle
    r_random = Win_streamer * sqrt(ran2(iseed))
    phi = 2.0*pi*ran2(iseed)
    x = x0 + r_random*( cos(phi)*n1 + sin(phi)*n2 )
    v = v0

    R_to_sink = sqrt( (x(1)-sink_pos(1))**2 + (x(2)-sink_pos(2))**2 + (x(3)-sink_pos(3))**2 )
    dum_rho = 1.0; dum_temp = 0.0
    call equationofstate(ieos, dum_ponrho, cs, dum_rho, R_to_sink, 0.0, 0.0, dum_temp)

    if (gamma > 1.01) then
       u = cs*cs / (gamma - 1.0)
    else
       u = 1.5 * cs*cs
    endif

    i_part = npart + 1
    call add_or_update_particle(igas, x, v, hguess, u, i_part, npart, npartoftype, xyzh, vxyzu)
 enddo

 if (iverbose >= 2) then
    print '(a,i8,2a,1pg12.4)', ' [streamer] injected N = ', Nin, '  (h=0.4*W_in, W_in=)', Win_streamer
 endif

 dtinject = huge(dtinject)
end subroutine inject_particles

!-----------------------------------------------------------------------
! Write options
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(mdot_streamer,'mdot_streamer','mass injection rate [Msun/yr]',iunit)
 call write_inopt(Rp_streamer,  'Rp_streamer',  'pericentre distance',iunit)
 call write_inopt(Rin_streamer, 'Rin_streamer', 'injection radius',iunit)
 call write_inopt(Rimp_streamer,'Rimp_streamer','impact radius on disc',iunit)
 call write_inopt(incl_streamer,'incl_streamer','inclination at impact [deg]',iunit)
 call write_inopt(phi_streamer, 'phi_streamer', 'node longitude [deg]',iunit)
 call write_inopt(Win_streamer, 'Win_streamer', 'streamer cross-section at injection',iunit)
 call write_inopt(ingoing,      'ingoing',      'TRUE=pre-pericentre',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
! Read options
!-----------------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(mdot_streamer,'mdot_streamer',db,errcount=nerr)
 call read_inopt(Rp_streamer,'Rp_streamer',db,errcount=nerr,min=0.)
 call read_inopt(Rin_streamer,'Rin_streamer',db,errcount=nerr,min=0.)
 call read_inopt(Rimp_streamer,'Rimp_streamer',db,errcount=nerr,min=0.)
 call read_inopt(incl_streamer,'incl_streamer',db,errcount=nerr)
 call read_inopt(phi_streamer,'phi_streamer',db,errcount=nerr)
 call read_inopt(Win_streamer,'Win_streamer',db,errcount=nerr,min=0.)
 call read_inopt(ingoing,'ingoing',db,errcount=nerr,default=.true.)

end subroutine read_options_inject

 !-----------------------------------------------------------------------
 ! Cross product routine
 !-----------------------------------------------------------------------
subroutine cross_product(a,b,c)
 real, intent(in)  :: a(3), b(3)
 real, intent(out) :: c(3)
 c(1) = a(2)*b(3) - a(3)*b(2)
 c(2) = a(3)*b(1) - a(1)*b(3)
 c(3) = a(1)*b(2) - a(2)*b(1)
end subroutine cross_product

 !-----------------------------------------------------------------------
 !+
 !  Updates the injected particles
 !+
 !-----------------------------------------------------------------------
subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

end module inject
