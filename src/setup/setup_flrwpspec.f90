!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Cosmological initial conditions for an Einstein-de Sitter universe
! from a matter power spectrum (or a single sine mode): the linear
! growing mode in longitudinal gauge (Macpherson et al. 2019), realised
! by particles moved from a cubic lattice with Zel'dovich displacements
! so that they reproduce the conserved density rho*, and given the
! linear-theory velocities.
!
! The box is 1 code unit across, with comoving length box_Mpc (Mpc/h);
! c = G = 1 and the initial scale factor is 1, so the initial Hubble
! rate is H = (H0 L/c) sqrt(1 + z_init) and the conformal time 2/H. The
! unit of length is the physical size of the box at z_init,
! L/(h (1 + z_init)) Mpc, and the unit of time is that length over c.
! Particles are dark matter by default (collisionless, on geodesics).
!
! A tabulated spectrum at redshift pk_redshift (in the synchronous gauge,
! e.g. from CAMB or CLASS) is scaled to z_init with the Einstein-de Sitter
! growth of the synchronous density contrast, delta_s ~ a on all scales.
!
! :References:
!   Macpherson, Price & Lasky (2019), Phys. Rev. D 99, 063522
!
! :Owner: Not Committed Yet
!
! :Runtime parameters:
!   - box_Mpc      : *comoving box length in Mpc/h*
!   - h_hubble     : *h, for the physical units and the Eisenstein-Hu spectrum*
!   - itype        : *particle type (1=gas, 5=dark matter)*
!   - lambda_min   : *remove modes shorter than this (Mpc/h, 0=none)*
!   - nref         : *grid for the random field (0=nx; same nref and seed give the same field)*
!   - nx           : *number of particles in each direction*
!   - perturbation : *pspec (power spectrum) or single (sine mode in x)*
!   - phi_ampl     : *amplitude of phi for the single mode*
!   - pk_file      : *power spectrum file (k [h/Mpc], P [(Mpc/h)^3]), or none*
!   - pk_redshift  : *redshift of the tabulated power spectrum*
!   - seed         : *random seed*
!   - sigma8       : *sigma8 for the Eisenstein-Hu spectrum*
!   - z_init       : *initial redshift*
!
! :Dependencies: boundary, dim, infile_utils, io, part, set_grf, setunits,
!   setup_params, utils_dft, utils_gr
!
 use setup_params, only:rhozero
 implicit none
 public :: setpart

 integer            :: npartx,itype,iseed,nref
 real               :: box_mpc,z_init,lambda_min,phi_ampl,sigma8,h_hubble,pk_redshift
 character(len=120) :: pk_file
 character(len=20)  :: perturbation

 ! speed of light in km/s, for H0 L/c with H0 = 100 h km/s/Mpc
 real, parameter :: c_kms = 299792.458

 private

contains

!----------------------------------------------------------------
!+
!  setup for cosmological initial conditions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu
 use io,           only:master,fatal
 use boundary,     only:set_boundary,dxbound
 use part,         only:idarkmatter
 use set_grf,      only:zeldovich_displacement
 use setup_params, only:npart_total
 use infile_utils, only:get_options
 use units,        only:set_units
 use physcon,      only:mpc
 use metric_et_utils, only:redshift_init
 integer,          intent(in)    :: id
 integer,          intent(inout) :: npart
 integer,          intent(out)   :: npartoftype(:)
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(out)   :: massoftype(:)
 real,             intent(out)   :: polyk,gamma
 real,             intent(inout) :: hfact
 real,             intent(inout) :: time
 character(len=*), intent(in)    :: fileprefix
 real,             intent(out)   :: vxyzu(:,:)
 real, allocatable :: phi(:,:,:),delta(:,:,:),vel(:,:,:,:),rhostar(:,:,:),psi(:,:,:,:)
 complex, allocatable :: phik(:,:,:)
 real    :: hub,dx,meanrho
 integer :: n,ierr

 ! defaults: the Macpherson et al. (2019) 1 Gpc/h box and spectrum (tabulated
 ! at z = 1000), started at z = 50 as in N-body simulations, so that the
 ! physical density contrast is well above the discreteness of the deposit
 npartx       = 32
 itype        = idarkmatter
 box_mpc      = 1000.
 z_init       = 50.
 perturbation = 'pspec'
 pk_file      = 'EdS_matterpower_synchronous_z1000.dat'
 pk_redshift  = 1000.
 iseed        = -12345
 nref         = 0
 lambda_min   = 0.
 phi_ampl     = 1.e-6
 sigma8       = 0.8
 h_hubble     = 0.45
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 ! c = G = 1, with the unit of length the physical box size at z_init
 call set_units(dist=box_mpc/h_hubble*mpc/(1.+z_init),c=1.d0,G=1.d0)
 ! so that the run reports the redshift (written to the .in file)
 redshift_init = z_init

 ! unit box, initial Hubble rate and density from the Friedmann equation
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.5,0.5)
 n   = npartx
 dx  = dxbound/n
 hub = 100.*box_mpc/c_kms*sqrt(1. + z_init)/dxbound
 rhozero = 3.*hub**2/(8.*acos(-1.))
 time = 2./hub
 if (id==master) print "(a,es12.5,a,es12.5,a,es12.5)",' initial H = ',hub,' rho = ',rhozero,' eta = ',time

 ! linear growing mode on the n^3 grid of lattice sites
 allocate(phik(n,n,n),phi(n,n,n),delta(n,n,n),vel(3,n,n,n),rhostar(n,n,n),psi(3,n,n,n))
 call get_phi_k(n,dx,hub,phik,ierr)
 if (ierr /= 0) call fatal('setup','could not read power spectrum from '//trim(pk_file))
 call get_linear_fields(n,dx,hub,phik,phi,delta,vel,rhostar)
 meanrho = sum(rhostar)/size(rhostar)
 if (id==master) print "(a,3(es10.3,1x))",' max |phi|, rms delta, rms delta* = ',maxval(abs(phi)),&
                       sqrt(sum(delta**2)/size(delta)),sqrt(sum((rhostar/meanrho-1.)**2)/size(rhostar))

 ! particles displaced from the lattice so that they trace rho*
 call zeldovich_displacement(n,dx,rhostar/meanrho - 1.,psi)
 call place_particles(n,dx,psi,vel,rhostar/meanrho,hfact,xyzh,vxyzu)
 npart = n**3
 npart_total = npart
 call set_types(npart,itype)
 npartoftype(:) = 0
 npartoftype(itype) = npart
 massoftype(:) = 0.
 massoftype(itype) = rhozero*meanrho*dxbound**3/npart
 if (id==master) print "(a,i0,a,es12.5)",' npart = ',npart,' particle mass = ',massoftype(itype)

 ! pressureless matter
 polyk = 0.
 gamma = 1.
 if (maxvxyzu >= 4) gamma = 4./3.

end subroutine setpart

!----------------------------------------------------------------
!+
!  Fourier transform of the growing-mode potential phi, either from
!  the power spectrum or for a single sine mode along x
!+
!----------------------------------------------------------------
subroutine get_phi_k(n,dx,hub,phik,ierr)
 use set_grf,   only:read_power_spectrum,pk_eisenstein_hu,white_noise_k,delta_sync_k,phi_k_from_delta_sync
 use utils_dft, only:dft_3d
 use boundary,  only:xmin,dxbound
 integer, intent(in)  :: n
 real,    intent(in)  :: dx,hub
 complex, intent(out) :: phik(n,n,n)
 integer, intent(out) :: ierr
 real, allocatable    :: ktab(:),pktab(:)
 complex, allocatable :: wk(:,:,:),deltak(:,:,:)
 integer :: i,nk

 ierr = 0
 select case(trim(perturbation))
 case('single')
    do i=1,n
       phik(i,:,:) = phi_ampl*sin(2.*acos(-1.)*(xmin + (i - 0.5)*dx)/dxbound)
    enddo
    call dft_3d(phik,-1)
 case default
    if (trim(pk_file)=='none') then
       ! tabulate the Eisenstein-Hu fit for an Einstein-de Sitter universe
       nk = 1000
       allocate(ktab(nk),pktab(nk))
       do i=1,nk
          ktab(i) = 10.**(-5. + 8.*(i-1)/real(nk-1))
          pktab(i) = pk_eisenstein_hu(ktab(i),z_init,h_hubble,1.,0.05,0.96,sigma8)
       enddo
    else
       call read_power_spectrum(pk_file,ktab,pktab,ierr)
       if (ierr /= 0) return
       ! the synchronous (comoving) density contrast of the growing mode
       ! grows as a on all scales in an Einstein-de Sitter universe, so the
       ! tabulated spectrum is scaled from its redshift to z_init
       pktab = pktab*((1. + pk_redshift)/(1. + z_init))**2
    endif
    allocate(wk(n,n,n),deltak(n,n,n))
    if (nref > 0 .and. nref <= n) then
       call white_noise_k(n,nref,iseed,wk)
    else
       call white_noise_k(n,n,iseed,wk)
    endif
    call delta_sync_k(n,box_mpc,wk,ktab,pktab,lambda_min,deltak)
    call phi_k_from_delta_sync(n,dx,hub,deltak,phik)
 end select

end subroutine get_phi_k

!----------------------------------------------------------------
!+
!  real-space potential, density contrast and velocity, and the
!  conserved density rho* = rho u^0 sqrt(-g) in units of rhozero,
!  in the perturbed metric ds^2 = -(1+2phi)dt^2 + (1-2phi)dx^2
!+
!----------------------------------------------------------------
subroutine get_linear_fields(n,dx,hub,phik,phi,delta,vel,rhostar)
 use set_grf,  only:linear_fields_from_phi_k
 use utils_gr, only:perturb_metric,get_u0,get_sqrtg
 integer, intent(in)  :: n
 real,    intent(in)  :: dx,hub
 complex, intent(in)  :: phik(n,n,n)
 real,    intent(out) :: phi(n,n,n),delta(n,n,n),vel(3,n,n,n),rhostar(n,n,n)
 real    :: gcov(0:3,0:3),sqrtg,u0
 integer :: i,j,k,ierr

 call linear_fields_from_phi_k(n,dx,hub,phik,phi,delta,vel)
 !$omp parallel do default(none) shared(n,phi,delta,vel,rhostar) &
 !$omp private(i,j,k,gcov,sqrtg,u0,ierr)
 do k=1,n
    do j=1,n
       do i=1,n
          call perturb_metric(phi(i,j,k),gcov)
          call get_sqrtg(gcov,sqrtg)
          ! coordinate velocity from the Eulerian velocity, v = alpha V
          vel(:,i,j,k) = vel(:,i,j,k)*sqrt(-gcov(0,0))
          call get_u0(gcov,vel(:,i,j,k),u0,ierr)
          rhostar(i,j,k) = (1. + delta(i,j,k))*u0*sqrtg
       enddo
    enddo
 enddo
 !$omp end parallel do

end subroutine get_linear_fields

!----------------------------------------------------------------
!+
!  particles at the lattice sites (cell centres of the grid) moved by
!  psi, wrapped into the periodic box, with velocity vel and smoothing
!  length set from the local conserved density rhos (in mean units)
!+
!----------------------------------------------------------------
subroutine place_particles(n,dx,psi,vel,rhos,hfact,xyzh,vxyzu)
 use boundary, only:xmin,ymin,zmin,dxbound,dybound,dzbound
 integer, intent(in)  :: n
 real,    intent(in)  :: dx,psi(3,n,n,n),vel(3,n,n,n),rhos(n,n,n),hfact
 real,    intent(out) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,j,k,ip

 !$omp parallel do default(none) shared(n,dx,psi,vel,rhos,hfact,xyzh,vxyzu) &
 !$omp shared(xmin,ymin,zmin,dxbound,dybound,dzbound) private(i,j,k,ip)
 do k=1,n
    do j=1,n
       do i=1,n
          ip = i + n*(j-1) + n*n*(k-1)
          xyzh(1,ip) = xmin + modulo((i-0.5)*dx + psi(1,i,j,k),dxbound)
          xyzh(2,ip) = ymin + modulo((j-0.5)*dx + psi(2,i,j,k),dybound)
          xyzh(3,ip) = zmin + modulo((k-0.5)*dx + psi(3,i,j,k),dzbound)
          xyzh(4,ip) = hfact*dx/rhos(i,j,k)**(1./3.)
          vxyzu(1:3,ip) = vel(:,i,j,k)
          if (size(vxyzu,1) >= 4) vxyzu(4,ip) = 0.
       enddo
    enddo
 enddo
 !$omp end parallel do

end subroutine place_particles

!----------------------------------------------------------------
!+
!  set the type of all particles
!+
!----------------------------------------------------------------
subroutine set_types(npart,itypei)
 use part, only:set_particle_type
 integer, intent(in) :: npart,itypei
 integer :: i

 do i=1,npart
    call set_particle_type(i,itypei)
 enddo

end subroutine set_types

!------------------------------------------------------------------------
!+
!  write setup file
!+
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(/,a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for cosmological initial conditions'
 write(iunit,"(/,a)") '# resolution and particles'
 call write_inopt(npartx,'nx','number of particles in each direction',iunit)
 call write_inopt(itype,'itype','particle type (1=gas, 5=dark matter)',iunit)
 write(iunit,"(/,a)") '# cosmology'
 call write_inopt(box_mpc,'box_Mpc','comoving box length in Mpc/h',iunit)
 call write_inopt(z_init,'z_init','initial redshift',iunit)
 write(iunit,"(/,a)") '# perturbations'
 call write_inopt(perturbation,'perturbation','pspec (power spectrum) or single (sine mode in x)',iunit)
 call write_inopt(pk_file,'pk_file','power spectrum file (k [h/Mpc], P [(Mpc/h)^3]), or none',iunit)
 call write_inopt(pk_redshift,'pk_redshift','redshift of the tabulated power spectrum',iunit)
 call write_inopt(iseed,'seed','random seed',iunit)
 call write_inopt(nref,'nref','grid for the random field (0=nx; same nref and seed give the same field)',iunit)
 call write_inopt(lambda_min,'lambda_min','remove modes shorter than this (Mpc/h, 0=none)',iunit)
 call write_inopt(sigma8,'sigma8','sigma8 for the Eisenstein-Hu spectrum',iunit)
 call write_inopt(h_hubble,'h_hubble','h, for the physical units and the Eisenstein-Hu spectrum',iunit)
 call write_inopt(phi_ampl,'phi_ampl','amplitude of phi for the single mode',iunit)
 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!+
!  read setup file
!+
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartx,'nx',db,min=8,errcount=nerr)
 call read_inopt(itype,'itype',db,min=1,errcount=nerr)
 call read_inopt(box_mpc,'box_Mpc',db,min=0.,errcount=nerr)
 call read_inopt(z_init,'z_init',db,min=0.,errcount=nerr)
 call read_inopt(perturbation,'perturbation',db,errcount=nerr)
 call read_inopt(pk_file,'pk_file',db,errcount=nerr)
 call read_inopt(pk_redshift,'pk_redshift',db,min=0.,default=1000.,errcount=nerr)
 call read_inopt(iseed,'seed',db,errcount=nerr)
 call read_inopt(nref,'nref',db,min=0,errcount=nerr)
 call read_inopt(lambda_min,'lambda_min',db,min=0.,errcount=nerr)
 call read_inopt(sigma8,'sigma8',db,min=0.,errcount=nerr)
 call read_inopt(h_hubble,'h_hubble',db,min=0.,errcount=nerr)
 call read_inopt(phi_ampl,'phi_ampl',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
