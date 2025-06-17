!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module relaxem
!
! relaxem
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: boundary, deriv, dim, eos, io, kernel, mpidomain, options,
!   part
!
 implicit none

contains

! Subroutine to relax the new set of particles to a reference particle distribution
subroutine relax_particles(npart,n_ref,xyzh_ref,force_ref,nrelax,relaxlist)
 use deriv,     only:get_derivs_global
 use dim,       only:mpi
 use io,        only:error
 integer,           intent(in)    :: npart,n_ref,nrelax
 real,              intent(in)    :: force_ref(3,n_ref),xyzh_ref(4,n_ref)
 integer,           intent(in)    :: relaxlist(1:nrelax)
 real,  allocatable :: a_ref(:,:)
 real :: ke,maxshift,ke_init,shuffle_tol
 logical :: converged
 integer :: ishift,nshifts

 write(*,"(/,70('-'),/,/,2x,a,/,/)") 'APR: time to relax ...'
 if (mpi) then
    call error('APR','relax_particles is not compatible with MPI')
    return
 endif

 write(*,"(1x,1(a,i8,a,i8,a))") 'Relaxing',nrelax,' particles the heavenly way from',n_ref,' references.'

 ! Initialise for the loop
 converged = .false.
 ishift = 0
 nshifts = 50
 shuffle_tol = 0.05

 ! a_ref stores the accelerations at the locations of the new particles as interpolated from the old ones
 allocate(a_ref(3,npart))

 do while (.not.converged)

    ! This gets fxyz of the new particles at their new locations
    call get_derivs_global()

    ! These are the accelerations at the locations of the new particles, interpolated from the parents
    call get_reference_accelerations(npart,a_ref,n_ref,xyzh_ref,force_ref,nrelax,relaxlist)

    ! Shift the particles by minimising the difference between the acceleration at the new particles and
    ! the interpolated values (i.e. what they do have minus what they should have)
    call shift_particles(npart,a_ref,nrelax,relaxlist,ke,maxshift)

    if (ishift == 0) ke_init = ke

    write(*,"(1x,1(a,f5.1,a,i3,a))") 'shuffle decreased to ',ke/ke_init*100.,'% of initial with',ishift,' shifts'

    ! Todo: cut-off criteria
    if (ishift >= nshifts .or. (ke/ke_init < shuffle_tol)) converged = .true.
    ishift = ishift + 1

 enddo

 ! Tidy up
 deallocate(a_ref)

 write(*,"(/,/,2x,a,/,/,70('-'))") 'APR: relaxing finished.'

end subroutine relax_particles

!----------------------------------------------------------------
!+
! Interpolates the accelerations at the locations of the new particles
! from the old set of particles (the reference particles)
!+
!----------------------------------------------------------------

subroutine get_reference_accelerations(npart,a_ref,n_ref,xyzh_ref,&
  force_ref,nrelax,relaxlist)
 use part,         only:xyzh,aprmassoftype,igas,apr_level,rhoh
 use dim,          only:periodic
 use kernel,       only:wkern,grkern,radkern2,cnormk
 use boundary,     only:dxbound,dybound,dzbound
 integer, intent(in) :: npart,n_ref,nrelax
 real,    intent(in) :: force_ref(3,n_ref),xyzh_ref(4,n_ref)
 integer, intent(in) :: relaxlist(nrelax)
 real,    intent(out) :: a_ref(3,npart)
 real :: xi,yi,zi,rij(3),h21,qj2,rij2,rhoj,h31,mass_ref,pmassi
 integer :: i,j,k

 a_ref(:,:) = 0.

 ! Over the new set of particles that are to be shuffled
 !$omp parallel do schedule(guided) default (none) &
 !$omp shared(xyzh,xyzh_ref,npart,n_ref,force_ref,a_ref,relaxlist) &
 !$omp shared(nrelax,apr_level,dxbound,dybound,dzbound) &
 !$omp shared(mass_ref,aprmassoftype) &
 !$omp private(i,j,xi,yi,zi,rij,h21,h31,rhoj,rij2,qj2,pmassi)

 over_new: do k = 1,nrelax
    if (relaxlist(k) == 0) cycle over_new
    i = relaxlist(k)
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    pmassi = aprmassoftype(igas,apr_level(i))

    ! Over the reference set of particles to which we are matching the accelerations
    over_reference: do j = 1,n_ref  ! later this should only be over active particles
       rij(1) = xyzh_ref(1,j) - xi
       rij(2) = xyzh_ref(2,j) - yi
       rij(3) = xyzh_ref(3,j) - zi
       mass_ref = aprmassoftype(igas,apr_level(j))   ! TBD: fix this to allow for dust

       if (periodic) then
          if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
          if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
          if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
       endif

       h21 = 1./(xyzh_ref(4,j))**2
       h31 = 1./(xyzh_ref(4,j))**3
       rhoj = rhoh(xyzh_ref(4,j),mass_ref)

       rij2 = dot_product(rij,rij)
       qj2  = rij2*h21

       if (qj2 < radkern2) then
          ! Interpolate acceleration at the location of the new particle
          a_ref(:,i) = a_ref(:,i) + force_ref(:,j)*wkern(qj2,sqrt(qj2))*cnormk*h31/rhoj
       endif

    enddo over_reference
 enddo over_new
 !$omp end parallel do

end subroutine get_reference_accelerations

!----------------------------------------------------------------
!+
! Calculates the shift for each particle, using the reference
! and interpolated accelerations
! (based off the routine in relax_star)
!+
!----------------------------------------------------------------

subroutine shift_particles(npart,a_ref,nrelax,relaxlist,ke,maxshift)
 use dim,      only:periodic
 use part,     only:xyzh,vxyzu,fxyzu,igas,aprmassoftype,rhoh, &
                     apr_level
 use eos,      only:get_spsound
 use options,  only:ieos
 use boundary, only:cross_boundary
 use mpidomain, only: isperiodic
 integer, intent(in)     :: npart,nrelax
 real,    intent(in)     :: a_ref(3,npart)
 integer, intent(in)     :: relaxlist(nrelax)
 real,    intent(out)    :: ke,maxshift
 real :: hi,rhoi,cs,dti,dx(3),vi(3),err,limit_bound
 real :: pmassi
 integer :: nlargeshift,i,ncross,j

 ke = 0.
 nlargeshift = 0
 ncross = 0
 maxshift = tiny(maxshift)
 limit_bound = 0.4 !! This probably shouldn't be more than 0.5

 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(npart,xyzh,vxyzu,fxyzu,ieos,a_ref,maxshift) &
 !$omp shared(apr_level,aprmassoftype) &
 !$omp shared(isperiodic,ncross,relaxlist,nrelax) &
 !$omp private(i,dx,dti,cs,rhoi,hi,vi,err,pmassi) &
 !$omp reduction(+:nlargeshift,ke)
 do j=1,nrelax
    if (relaxlist(j) == 0) cycle
    i = relaxlist(j)
    hi = xyzh(4,i)
    pmassi = aprmassoftype(igas,apr_level(i))
    rhoi = rhoh(hi,pmassi)
    cs = get_spsound(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
    dti = 0.3*hi/cs   ! h/cs

    dx  = 0.5*dti**2*(fxyzu(1:3,i) - a_ref(1:3,i))
    if (sqrt(dot_product(dx,dx)) > maxshift) maxshift = sqrt(dot_product(dx,dx))
    if (dot_product(dx,dx) > hi**2) then

       dx = dx / sqrt(dot_product(dx,dx)) * hi  ! Avoid large shift in particle position !check with what James has done
       nlargeshift = nlargeshift + 1
    endif

    ! actual shift
    xyzh(1:3,i) = xyzh(1:3,i) + dx(:)

    ! if periodic, move to the other side of the box
    ! (written locally but ideally should call cross_boundary)
    if (periodic) call cross_boundary(isperiodic,xyzh(:,i),ncross)

    ! faux velocities, so we can estimate the magnitude of the shift
    vi(1:3) = dx(:)/dti
    ke = ke + dot_product(vi,vi)

    ! Output for testing purposes
    err = sqrt(dot_product(dx,dx))/hi
    if (err > maxshift) maxshift = err

 enddo
 !$omp end parallel do
 if (nlargeshift > 0) print*,'Warning: Restricted dx for ', nlargeshift, 'particles'


end subroutine shift_particles

!----------------------------------------------------------------
!+
! For each particle that has been shuffled, check the minimum
! distance between to see if we have any pairing issues
! and specifically, is the shuffling making it worse
!+
!----------------------------------------------------------------

subroutine check_for_pairing(nrelax,relaxlist,pair_distance)
 use part, only:xyzh
 integer, intent(in) :: nrelax,relaxlist(nrelax)
 real, intent(out)   :: pair_distance
 real :: dx(3), dx_mag
 integer :: ii,jj

 pair_distance = huge(pair_distance)

 do ii = 1,nrelax
    do jj = 1,nrelax
       if (ii == jj) cycle
       dx = xyzh(1:3,ii) - xyzh(1:3,jj)
       dx_mag = sqrt(dot_product(dx,dx))/xyzh(4,ii) ! scaled by the smoothing length

       if (dx_mag < pair_distance) pair_distance = dx_mag

    enddo
 enddo

end subroutine check_for_pairing

end module relaxem
