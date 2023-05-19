!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module relaxstar
!
! Automated relaxation of stellar density profile,
!   iterating towards hydrostatic equilibrium
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - maxits   : *maximum number of relaxation iterations*
!   - tol_dens : *% error in density to stop relaxation*
!   - tol_ekin : *tolerance on ekin/epot to stop relaxation*
!
! :Dependencies: checksetup, damping, deriv, dim, dump_utils, energies,
!   eos, externalforces, fileutils, infile_utils, initial, io, io_summary,
!   memory, options, part, physcon, ptmass, readwrite_dumps, setstar_utils,
!   sortutils, step_lf_global, table_utils, units
!
 implicit none
 public :: relax_star,write_options_relax,read_options_relax

 real,    private :: tol_ekin = 1.e-7 ! criteria for being converged
 real,    private :: tol_dens = 1.   ! allow 1% RMS error in density
 integer, private :: maxits = 1000

 real,    private :: gammaprev,hfactprev,mass1prev
 integer, private :: ieos_prev

 integer, public :: ierr_setup_errors = 1, &
                    ierr_no_pressure  = 2, &
                    ierr_unbound = 3, &
                    ierr_notconverged = 4

 private

contains

!----------------------------------------------------------------
!+
!  relax a star to hydrostatic equilibrium. We run the main
!  code but with a fake equation of state, low neighbour number
!  and fixing the entropy as a function of r
!
!  IN:
!    rho(nt)   - tabulated density as function of r (in code units)
!    pr(nt)    - tabulated pressure as function of r (in code units)
!    r(nt)     - radius for each point in the table
!
!  IN/OUT:
!    xyzh(:,:) - positions and smoothing lengths of all particles
!+
!----------------------------------------------------------------
subroutine relax_star(nt,rho,pr,r,npart,xyzh,use_var_comp,Xfrac,Yfrac,mu,ierr,npin,label)
 use table_utils,     only:yinterp
 use deriv,           only:get_derivs_global
 use dim,             only:maxp,maxvxyzu,gr,gravity
 use part,            only:vxyzu,rad,eos_vars,massoftype,igas
 use step_lf_global,  only:init_step,step
 use initial,         only:initialise
 use memory,          only:allocate_memory
 use energies,        only:compute_energies,ekin,epot,etherm
 use checksetup,      only:check_setup
 use io,              only:error,warning,fatal,id,master
 use fileutils,       only:getnextfilename
 use readwrite_dumps, only:write_fulldump,init_readwrite_dumps
 use eos,             only:gamma,eos_outputs_mu
 use physcon,         only:pi
 use options,         only:iexternalforce
 use io_summary,      only:summary_initialise
 use setstar_utils,   only:set_star_thermalenergy,set_star_composition
 integer, intent(in)    :: nt
 integer, intent(inout) :: npart
 real,    intent(in)    :: rho(nt),pr(nt),r(nt)
 logical, intent(in)    :: use_var_comp
 real,    intent(in), allocatable :: Xfrac(:),Yfrac(:),mu(:)
 real,    intent(inout) :: xyzh(:,:)
 integer, intent(out)   :: ierr
 integer, intent(in), optional :: npin
 character(len=*), intent(in), optional :: label
 integer :: nits,nerr,nwarn,iunit,i1
 real    :: t,dt,dtmax,rmserr,rstar,mstar,tdyn
 real    :: entrop(nt),utherm(nt),mr(nt),rmax,dtext,dtnew
 logical :: converged,use_step,restart
 logical, parameter :: fix_entrop = .true. ! fix entropy instead of thermal energy
 logical, parameter :: write_files = .true.
 character(len=20) :: filename,mylabel

 i1 = 0
 if (present(npin)) i1 = npin  ! starting position in particle array
 !
 ! label for relax_star snapshots
 !
 mylabel = ''
 if (present(label)) mylabel = label
 !
 ! save settings and set a bunch of options
 !
 ierr = 0
 rstar = maxval(r)
 mr = get_mr(rho,r)
 mstar = mr(nt)
 tdyn  = 2.*pi*sqrt(rstar**3/(32.*mstar))
 if (id==master) print*,'rstar = ',rstar,' mstar = ',mstar, ' tdyn = ',tdyn
 !
 ! see if we can restart or skip the relaxation process
 ! based on a previous run
 !
 filename = 'relax'//trim(mylabel)//'_00000'
 if (write_files) call init_readwrite_dumps()
 call check_for_existing_file(filename,npart,massoftype(igas),&
                              xyzh,vxyzu,restart,ierr)
 !
 ! quit with fatal error if non-matching file found, otherwise
 ! this will be overwritten
 !
 if (write_files .and. ierr /= 0) then
    if (id==master) print "(a)",' ERROR: pre-existing relaxed star dump(s) do not match'
    call fatal('relax_star','please delete relax'//trim(mylabel)//'_* and restart...')
 endif

 call set_options_for_relaxation(tdyn)
 call summary_initialise()
 !
 ! check particle setup is sensible
 !
 call check_setup(nerr,nwarn,restart=.true.) ! restart=T allows accreted/masked particles
 if (nerr > 0) then
    call error('relax_star','cannot relax star because particle setup contains errors')
    call restore_original_options(i1,npart)
    ierr = ierr_setup_errors
    return
 endif
 use_step = .false.

 if (iexternalforce > 0 .and. (.not. gr)) then
    call warning('relax_star','asynchronous shifting not implemented with external forces: evolving in time instead')
    use_step = .true.
 endif
 !
 ! define utherm(r) based on P(r) and rho(r)
 ! and use this to set the thermal energy of all particles
 !
 entrop = pr/rho**gamma
 utherm = pr/(rho*(gamma-1.))
 if (any(utherm <= 0.)) then
    call error('relax_star','relax-o-matic needs non-zero pressure array set in order to work')
    call restore_original_options(i1,npart)
    ierr = ierr_no_pressure
    return
 endif
 call reset_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,rho,&
                             utherm,entrop,fix_entrop,rmax,rmserr)
 !
 ! compute derivatives the first time around (needed if using actual step routine)
 !
 t = 0.
 call allocate_memory(int(min(2*npart,maxp),kind=8))
 call get_derivs_global()
 call reset_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,rho,&
                             utherm,entrop,fix_entrop,rmax,rmserr)
 call compute_energies(t)
 !
 ! perform sanity checks
 !
 if (etherm > abs(epot)) then
    call error('relax_star','cannot relax star because it is unbound (etherm > epot)')
    if (id==master) print*,' Etherm = ',etherm,' Epot = ',Epot
    if (maxvxyzu < 4) print "(/,a,/)",' *** Try compiling with ISOTHERMAL=no instead... ***'
    call restore_original_options(i1,npart)
    ierr = ierr_unbound
    return
 endif
 if (id==master) print "(/,3(a,1pg11.3),/,a,0pf6.2,a,es11.3,a,i4)",&
   ' RELAX-A-STAR-O-MATIC: Etherm:',etherm,' Epot:',Epot, ' R*:',maxval(r), &
   '       WILL stop WHEN: dens error < ',tol_dens,'% AND Ekin/Epot < ',tol_ekin,' OR Iter=',maxits

 if (write_files) then
    if (.not.restart) call write_fulldump(t,filename)
    open(newunit=iunit,file='relax'//trim(mylabel)//'.ev',status='replace')
    write(iunit,"(a)") '# nits,rmax,etherm,epot,ekin/epot,L2_{err}'
 endif
 converged = .false.
 dt = epsilon(0.) ! To avoid error in sink-gas substepping
 dtext = huge(dtext)
 if (use_step) then
    dtmax = tdyn
    call init_step(npart,t,dtmax)
 endif
 nits = 0
 do while (.not. converged .and. nits < maxits)
    nits = nits + 1
    !
    ! shift particles by one "timestep"
    !
    t = t + dt
    if (use_step) then
       call step(npart,npart,t,dt,dtext,dtnew)
       dt = dtnew
    else
       call shift_particles(i1,npart,xyzh,vxyzu,dt)
    endif
    !
    ! reset thermal energy and calculate information
    !
    call reset_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,&
         rho,utherm,entrop,fix_entrop,rmax,rmserr)
    !
    ! compute energies and check for convergence
    !
    call compute_energies(t)
    converged = (ekin > 0. .and. ekin/abs(epot) < tol_ekin .and. rmserr < 0.01*tol_dens)
    !
    ! print information to screen
    !
    if (use_step) then
       if (id==master) print "(a,es10.3,a,2pf6.2,2(a,1pg11.3))",&
        ' Relaxing star: t/dyn:',t/tdyn,', dens error:',rmserr,'%, R*:',rmax, &
        ' Ekin/Epot:',ekin/abs(epot)
    else
       if (id==master) print "(a,i4,a,i4,a,2pf6.2,2(a,1pg11.3))",&
        ' Relaxing star: Iter',nits,'/',maxits, &
        ', dens error:',rmserr,'%, R*:',rmax,' Ekin/Epot:',ekin/abs(epot)
    endif
    !
    ! additional diagnostic output, mainly for debugging/checking
    !
    if (write_files) then
       !
       ! write information to the relax.ev file
       !
       write(iunit,*) nits,rmax,etherm,epot,ekin/abs(epot),rmserr
       !
       ! write dump files
       !
       if (mod(nits,100)==0 .or. ((nits==maxits .or. converged).and.nits > 1)) then
          filename = getnextfilename(filename)
          !
          ! before writing a file, set the real thermal energy profile
          ! so the file is useable as a starting file for the main calculation
          !
          if (use_var_comp) call set_star_composition(use_var_comp,&
                                 eos_outputs_mu(ieos_prev),npart,xyzh,&
                                 Xfrac,Yfrac,mu,mr,mstar,eos_vars,npin=i1)

          if (maxvxyzu==4) call set_star_thermalenergy(ieos_prev,rho,pr,&
                                r,nt,npart,xyzh,vxyzu,rad,eos_vars,.true.,&
                                use_var_comp=.false.,initialtemp=1.e3,npin=i1)

          ! write relaxation snapshots
          if (write_files) call write_fulldump(t,filename)

          ! flush the relax.ev file
          call flush(iunit)

          ! restore the fake thermal energy profile
          call reset_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,rho,&
               utherm,entrop,fix_entrop,rmax,rmserr)
       endif
    endif
 enddo
 if (write_files) close(iunit)
 !
 ! warn if relaxation finished due to hitting nits=nitsmax
 !
 if (.not.converged) then
    call warning('relax_star','relaxation did not converge, just reached max iterations')
    ierr = ierr_notconverged
 else
    if (id==master) print "(5(a,/))",&
    "                             _                    _ ",&
    "                    _ __ ___| | __ ___  _____  __| |",&
    "                   | '__/ _ \ |/ _` \ \/ / _ \/ _` |",&
    "                   | | |  __/ | (_| |>  <  __/ (_| |",&
    "          o  o  o  |_|  \___|_|\__,_/_/\_\___|\__,_|  o  o  o"
 endif
 !
 ! unfake some things
 !
 call restore_original_options(i1,npart)

end subroutine relax_star

!----------------------------------------------------------------
!+
!  shift particles: this is like timestepping but done
!  asynchronously. Each particle shifts by dx = 0.5*dt^2*a
!  where dt is the local courant timestep, i.e. h/c_s
!+
!----------------------------------------------------------------
subroutine shift_particles(i1,npart,xyzh,vxyzu,dtmin)
 use deriv, only:get_derivs_global
 use part,  only:fxyzu,fext,xyzmh_ptmass,nptmass,rhoh,massoftype,igas
 use ptmass,only:get_accel_sink_gas
 use eos,   only:get_spsound
 use options, only:ieos
 integer, intent(in) :: i1,npart
 real, intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real, intent(out)   :: dtmin
 real :: dx(3),dti,phi,rhoi,cs,hi
 integer :: i,nlargeshift
!
! shift particles asynchronously
!
 dtmin = huge(dtmin)
 nlargeshift = 0
 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(i1,npart,xyzh,vxyzu,fxyzu,fext,xyzmh_ptmass,nptmass,massoftype,ieos) &
 !$omp private(i,dx,dti,phi,cs,rhoi,hi) &
 !$omp reduction(min:dtmin) &
 !$omp reduction(+:nlargeshift)
 do i=i1+1,npart
    fext(1:3,i) = 0.
    if (nptmass > 0) then
       call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                               xyzmh_ptmass,fext(1,i),fext(2,i),fext(3,i),phi)
    endif
    hi = xyzh(4,i)
    rhoi = rhoh(hi,massoftype(igas))
    cs = get_spsound(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
    dti = 0.3*hi/cs   ! h/cs
    dx  = 0.5*dti**2*(fxyzu(1:3,i) + fext(1:3,i))
    if (dot_product(dx,dx) > hi**2) then
       dx = dx / sqrt(dot_product(dx,dx)) * hi  ! Avoid large shift in particle position
       nlargeshift = nlargeshift + 1
    endif
    xyzh(1:3,i) = xyzh(1:3,i) + dx(:)
    vxyzu(1:3,i) = dx(:)/dti ! fake velocities, so we can measure kinetic energy
    dtmin = min(dtmin,dti)   ! used to print a "time" in the output (but it is fake)
 enddo
 !$omp end parallel do
 if (nlargeshift > 0) print*,'Warning: Restricted dx for ', nlargeshift, 'particles'
 !
 ! get forces on particles
 !
 call get_derivs_global()

end subroutine shift_particles

!----------------------------------------------------------------
!+
!  reset the thermal energy to be exactly p(r)/((gam-1)*rho(r))
!  according to the desired p(r) and rho(r)
!  also compute error between true rho(r) and desired rho(r)
!+
!----------------------------------------------------------------
subroutine reset_u_and_get_errors(i1,npart,xyzh,vxyzu,rad,nt,mr,rho,&
                                  utherm,entrop,fix_entrop,rmax,rmserr)
 use table_utils, only:yinterp
 use sortutils,   only:find_rank,r2func
 use part,        only:rhoh,massoftype,igas,maxvxyzu,iorder=>ll
 use dim,         only:do_radiation
 use eos,         only:gamma
 integer, intent(in) :: i1,npart,nt
 real, intent(in)    :: xyzh(:,:),mr(nt),rho(nt),utherm(nt),entrop(nt)
 real, intent(inout) :: vxyzu(:,:),rad(:,:)
 real, intent(out)   :: rmax,rmserr
 logical, intent(in) :: fix_entrop
 real :: ri,rhor,rhoi,rho1,mstar,massri
 integer :: i

 rho1 = yinterp(rho,mr,0.)
 rmax = 0.
 rmserr = 0.
 call find_rank(npart-i1,r2func,xyzh(1:3,i1+1:npart),iorder)
 mstar = mr(nt)
 do i = i1+1,npart
    ri = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    massri = mstar * real(iorder(i-i1)-1) / real(npart-i1)
    !if (i1 > 0 .and. i-i1 < 10) print*,' r=  ',ri,' massri=',massri,iorder(i-i1),npart-i1
    rhor = yinterp(rho,mr,massri) ! analytic rho(r)
    rhoi = rhoh(xyzh(4,i),massoftype(igas)) ! actual rho
    if (maxvxyzu >= 4) then
       if (fix_entrop) then
          vxyzu(4,i) = (yinterp(entrop,mr,massri)*rhoi**(gamma-1.))/(gamma-1.)
       else
          vxyzu(4,i) = yinterp(utherm,mr,massri)
       endif
    endif
    rmserr = rmserr + (rhor - rhoi)**2
    rmax   = max(rmax,ri)
 enddo
 if (do_radiation) rad = 0.
 rmserr = sqrt(rmserr/npart)/rho1

end subroutine reset_u_and_get_errors

!----------------------------------------------------------------
!+
!  set code options specific to relaxation calculations
!+
!----------------------------------------------------------------
subroutine set_options_for_relaxation(tdyn)
 use eos,  only:ieos,gamma
 use part, only:hfact,maxvxyzu,gr
 use damping, only:damp,tdyn_s
 use options, only:idamp
 use units,          only:utime
 use externalforces, only:mass1
 real, intent(in) :: tdyn

 gammaprev = gamma
 hfactprev = hfact
 ieos_prev = ieos
 mass1prev = mass1
 !
 ! turn on settings appropriate to relaxation
 !
 if (maxvxyzu >= 4) ieos = 2
 if (tdyn > 0.) then
    idamp = 2
    tdyn_s = tdyn*utime
 else
    idamp = 1
    damp = 0.05
 endif
 if (gr) mass1 = 0. ! use Minkowski metric during relaxation

end subroutine set_options_for_relaxation
!----------------------------------------------------------------
!+
!  check if a previous snapshot exists of the relaxed star
!+
!----------------------------------------------------------------
subroutine check_for_existing_file(filename,npart,mgas,xyzh,vxyzu,restart,ierr)
 use dump_utils, only:open_dumpfile_r,read_header,dump_h,lenid,extract
 use fileutils,  only:getnextfilename
 use io,         only:idump,idisk1,id,nprocs,iprint
 use readwrite_dumps, only:read_dump
 character(len=*), intent(inout) :: filename
 integer,          intent(in)    :: npart
 real,             intent(in)    :: mgas
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical,          intent(out)   :: restart
 integer,          intent(out)   :: ierr
 logical :: iexist,tagged
 character(len=len(filename)) :: restart_file,filetmp
 character(len=lenid) :: fileid
 type(dump_h) :: hdr
 integer :: npartfile
 real :: hfactfile,tfile,mfile
 !
 ! check for the last file in the list relax_00000, relax_00001 etc
 !
 ierr = 0
 iexist = .true.
 filetmp = filename
 restart_file = ''
 restart = .false.
 do while (iexist)
    inquire(file=filetmp,exist=iexist)
    if (iexist) then
       restart_file = filetmp
       filetmp = getnextfilename(filetmp)
    endif
 enddo
 if (len_trim(restart_file) <= 0) return

 print "(/,1x,a)",'>> RESTARTING relaxation from '//trim(restart_file)
 call open_dumpfile_r(idump,restart_file,fileid,ierr)
 call read_header(idump,hdr,tagged,ierr)
 close(idump)
 if (ierr /= 0) then
    print "(a)",' ERROR: could not read file header'
    return
 else
    call extract('nparttot',npartfile,hdr,ierr,default=0)
    if (npartfile /= npart) then
       print "(a,i0,a,i0)",' ERROR: np=',npartfile,' in '//&
         trim(restart_file)//' differs from current np=',npart
       ierr = 2
       return
    else
       call extract('massoftype',mfile,hdr,ierr,default=0.)
       if (abs(mfile-mgas) > epsilon(0.)) then
          print "(a,es10.3,a,es10.3)",' ERROR: M=',npart*mfile,' in '//&
            trim(restart_file)//' differs from current M=',npart*mgas
          ierr = 3
          return
       else
          restart = .true.
          call read_dump(restart_file,tfile,hfactfile,&
                         idisk1,iprint,id,nprocs,ierr)
          filename = restart_file
       endif
    endif
 endif

end subroutine check_for_existing_file

!--------------------------------------------------
!+
!  get mass coordinate m(r) = \int 4.*pi*rho*r^2 dr
!+
!--------------------------------------------------
function get_mr(rho,r) result(mr)
 use physcon, only:pi
 real, intent(in)  :: rho(:),r(:)
 real :: mr(size(r))
 integer :: i

 mr(1) = 0.
 do i=2,size(rho)
    mr(i) = mr(i-1) + 4./3.*pi*rho(i) * (r(i) - r(i-1)) * (r(i)**2 + r(i)*r(i-1) + r(i-1)**2)
 enddo

end function get_mr

!----------------------------------------------------------------
!+
!  restore previous settings
!+
!----------------------------------------------------------------
subroutine restore_original_options(i1,npart)
 use eos,     only:ieos,gamma
 use damping, only:damp
 use options, only:idamp
 use part,    only:hfact,vxyzu,gr
 use externalforces, only:mass1
 integer, intent(in) :: i1,npart

 gamma = gammaprev
 hfact = hfactprev
 ieos  = ieos_prev
 idamp = 0
 damp = 0.
 vxyzu(1:3,i1+1:npart) = 0.
 if (gr) mass1 = mass1prev

end subroutine restore_original_options

!----------------------------------------------------------------
!+
!  write relaxation options to .setup file
!+
!----------------------------------------------------------------
subroutine write_options_relax(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(tol_ekin,'tol_ekin','tolerance on ekin/epot to stop relaxation',iunit)
 call write_inopt(tol_dens,'tol_dens','% error in density to stop relaxation',iunit)
 call write_inopt(maxits,'maxits','maximum number of relaxation iterations',iunit)

end subroutine write_options_relax

!----------------------------------------------------------------
!+
!  read relaxation options from .setup file
!+
!----------------------------------------------------------------
subroutine read_options_relax(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(tol_ekin,'tol_ekin',db,errcount=nerr)
 call read_inopt(tol_dens,'tol_dens',db,errcount=nerr)
 call read_inopt(maxits,'maxits',db,errcount=nerr)

end subroutine read_options_relax

end module relaxstar
