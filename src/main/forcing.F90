!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: forcing
!
!  DESCRIPTION:
!  This module is from Christoph Federrath => applies
!  turbulent driving via an Ornstein-Uhlenbeck process.
!
!  Designed to be a drop-in replacement for forcing.F90
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    istir        -- switch to turn stirring on or off at runtime
!    st_amplfac   -- amplitude factor for stirring of turbulence
!    st_decay     -- correlation time for driving
!    st_dtfreq    -- frequency of stirring
!    st_energy    -- energy input/mode
!    st_seed      -- random number generator seed
!    st_solweight -- solenoidal weight
!    st_spectform -- spectral form of stirring
!    st_stirmax   -- maximum stirring wavenumber
!    st_stirmin   -- minimum stirring wavenumber
!
!  DEPENDENCIES: boundary, datafiles, fileutils, infile_utils, io,
!    mpiutils, part
!+
!--------------------------------------------------------------------------

module forcing

 public :: forceit,init_forcing,write_forcingdump,write_options_forcing,read_options_forcing

 integer, parameter :: st_maxmodes = 1000

 !OU variance corresponding to decay time and energy input rate
 real, save :: st_OUvar

 !last time random seeds were updated
 real, save :: tprev

 !Number of modes
 integer, save :: st_nmodes

 real, save :: st_mode(3,st_maxmodes), st_aka(3,st_maxmodes), st_akb(3,st_maxmodes)
 real, save :: st_OUphases(6*st_maxmodes)
 real, save :: st_ampl(  st_maxmodes)

 !Options - give these default values to write to file
 integer, save, public  :: istir = 1
 !logical,save  :: st_computeDt
 real, save, public     :: st_decay = 0.05
 real, save, public     :: st_energy = 2.0
 real, save, public     :: st_stirmin = 6.28
 real, save, public     :: st_stirmax = 18.86
 real, save, public     :: st_solweight = 1.0
 real, save, public     :: st_solweightnorm
 integer, save  :: st_seed = 1
 !integer,save  :: st_freq = 1
 real, save, public    :: st_dtfreq = 0.01
 real, save, public    :: st_amplfac = 1.0
 integer, save, public  :: st_spectform = 1

 !!namelist /force/ istir,st_spectform,st_stirmax,st_stirmin, &
 !                 st_energy,st_decay,st_solweight,st_dtfreq,st_seed
 !private :: force

 integer, save, allocatable :: st_randseed(:)
 integer, save                            :: st_seedLen
 integer, parameter :: ndim = 3

#define N_DIM 3
!!!!#undef CORRECT_MEAN_FORCE
!!!!##define STIR_FROM_FILE

#ifdef STIR_FROM_FILE
 character(len=*), parameter :: forcingfile = 'forcing.dat'
#endif

 private

contains

!! NAME
!!
!!  init_stir
!!
!! SYNOPSIS
!!
!!  init_stir()
!!
!! DESCRIPTION
!!  initialize the stirring module
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!   These are the runtime parameters used in the Stir unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.
!!
!!    istir [INTEGER]
!!        Switch to turn stirring on or off at runtime.
!!    st_decay [REAL]
!!        correlation time for driving
!!    st_energy [REAL]
!!        energy input/mode
!!    st_freq [INTEGER]
!!        frequency of stirring
!!    st_seed [INTEGER]
!!        random number generator seed
!!    st_stirmax [REAL]
!!        maximum stirring *wavenumber*
!!    st_stirmin [REAL]
!!        minimum stirring *wavenumber*
!!    st_solweight [REAL]
!!        solenoidal weight
!!    st_computeDt {BOOLEAN]
!!        whether to restrict timestep based on stirring
!!
!! AUTHOR: Christoph Federrath
!! modified for use in PhantomSPH by Daniel Price 2008
!!
!!***

subroutine init_forcing(dumpfile,infile,time)
 use boundary, only:dxbound,dybound,dzbound, Lx=>dxbound,Ly=>dybound,Lz=>dzbound
 use io, only:id,master,fatal,error,die,iprint,iverbose
 use infile_utils, only:get_inopt
 use part, only:vxyzu
 use fileutils, only:numfromfile


 character(len=*), intent(in) :: dumpfile,infile
 real,             intent(in) :: time
 logical                           :: restart
 integer                           :: ikxmin, ikxmax, ikymin, ikymax, ikzmin, ikzmax
 integer                           :: ikx, iky, ikz
 real                              :: kx, ky, kz, k, kc
#ifdef STIR_FROM_FILE
 real                              :: timeinfile
#endif
 real, parameter                   :: twopi = 6.283185307
 real, parameter                   :: amin = 0.0 ! the amplitude of the modes at kmin and kmax for a parabolic Fourier spectrum wrt 1.0 at the centre kc
 !logical                           :: iexist
 !character(len=len(infile)+5)      :: drvfile
 logical, parameter                :: Debug = .false.

 if (istir == 0) return

!  call get_parm_from_context('st_spectform', st_spectform) ! 0 is peak, 1 is paraboloid
!  call get_parm_from_context('st_decay', st_decay)
!  call get_parm_from_context('st_energy', st_energy)
!  call get_parm_from_context('st_stirmin', st_stirmin)
!  call get_parm_from_context('st_stirmax', st_stirmax)
!  call get_parm_from_context('st_solweight', st_solweight)
!  call get_parm_from_context('st_seed', st_seed)
!  call get_parm_from_context('st_freq', st_freq)
!  call get_parm_from_context('st_computeDt', st_computeDt)

 !initialize some variables, allocate randseed

 st_OUvar         = sqrt(st_energy/st_decay)
 kc               = 0.5*(st_stirmin+st_stirmax)
#if N_DIM == 3
 st_solweightnorm = sqrt(3.0/3.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*st_solweight+3.0*st_solweight**2.0) ! this makes the rms force const irrespective of the solenoidal weight
#endif
#if N_DIM == 2
 st_solweightnorm = sqrt(3.0/2.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*st_solweight+2.0*st_solweight**2.0) ! this makes the rms force const irrespective of the solenoidal weight
#endif
#if N_DIM == 1
 st_solweightnorm = sqrt(3.0/1.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*st_solweight+1.0*st_solweight**2.0) ! this makes the rms force const irrespective of the solenoidal weight
#endif

 ikxmin = 0
 ikymin = 0
 ikzmin = 0

 ikxmax = 20
 ikymax = 0
 ikzmax = 0
#if N_DIM > 1
 ikymax = 8
#if N_DIM > 2
 ikzmax = 8
#endif
#endif

 !Lx = dxbound
 !Ly = dybound
 !Lz = dzbound

 st_nmodes = 0
 !call add_parm_to_context (global_parm_context, "st_nmodes", st_nmodes)

 do ikx = ikxmin, ikxmax
    kx = twopi * ikx / Lx

    do iky = ikymin, ikymax
       ky = twopi * iky / Ly

       do ikz = ikzmin, ikzmax
          kz = twopi * ikz / Lz

          k = sqrt( kx*kx+ky*ky+kz*kz )

          if ((k  >=  st_stirmin) .and. (k  <=  st_stirmax)) then

             if ((st_nmodes + 2**(ndim-1))  >  st_maxmodes) then

                if (id==master) print *,'init_stir:  st_nmodes = ', st_nmodes, ' maxstirmodes = ',st_maxmodes
                call fatal('init_stir','Too many stirring modes')

             endif

             st_nmodes = st_nmodes + 1

             if (st_spectform == 0) st_ampl(st_nmodes) = 1.0
             if (st_spectform == 1) st_ampl(st_nmodes) = 4.0*(amin-1.0)/((st_stirmax-st_stirmin)**2.0)*((k-kc)**2.0)+1.0
             if (st_spectform == 2) st_ampl(st_nmodes) = (st_stirmin/k)**(5./3.)
             if (id==master .and. Debug) print *, "init_stir:  st_ampl(",st_nmodes,") = ", st_ampl(st_nmodes)

             st_mode(1,st_nmodes) = kx
             st_mode(2,st_nmodes) = ky
             st_mode(3,st_nmodes) = kz

#if N_DIM > 1
             st_nmodes = st_nmodes + 1

             if (st_spectform == 0) st_ampl(st_nmodes) = 1.0
             if (st_spectform == 1) st_ampl(st_nmodes) = 4.0*(amin-1.0)/((st_stirmax-st_stirmin)**2.0)*((k-kc)**2.0)+1.0
             if (st_spectform == 2) st_ampl(st_nmodes) = (st_stirmin/k)**(5./3.)
             if (id==master .and. Debug) print *, "init_stir:  st_ampl(",st_nmodes,") = ", st_ampl(st_nmodes)

             st_mode(1,st_nmodes) = kx
             st_mode(2,st_nmodes) =-ky
             st_mode(3,st_nmodes) = kz
#endif

#if N_DIM > 2
             st_nmodes = st_nmodes + 1

             if (st_spectform == 0) st_ampl(st_nmodes) = 1.0
             if (st_spectform == 1) st_ampl(st_nmodes) = 4.0*(amin-1.0)/((st_stirmax-st_stirmin)**2.0)*((k-kc)**2.0)+1.0
             if (st_spectform == 2) st_ampl(st_nmodes) = (st_stirmin/k)**(5./3.)
             if (id==master .and. Debug) print *, "init_stir:  st_ampl(",st_nmodes,") = ", st_ampl(st_nmodes)

             st_mode(1,st_nmodes) = kx
             st_mode(2,st_nmodes) = ky
             st_mode(3,st_nmodes) =-kz

             st_nmodes = st_nmodes + 1

             if (st_spectform == 0) st_ampl(st_nmodes) = 1.0
             if (st_spectform == 1) st_ampl(st_nmodes) = 4.0*(amin-1.0)/((st_stirmax-st_stirmin)**2.0)*((k-kc)**2.0)+1.0
             if (st_spectform == 2) st_ampl(st_nmodes) = (st_stirmin/k)**(5./3.)
             if (id==master .and. Debug) print *, "init_stir:  st_ampl(",st_nmodes,") = ", st_ampl(st_nmodes)

             st_mode(1,st_nmodes) = kx
             st_mode(2,st_nmodes) =-ky
             st_mode(3,st_nmodes) =-kz
#endif

             !call set_parm_in_context(global_parm_context, "st_nmodes", st_nmodes)

          endif

       enddo
    enddo
 enddo

!--restart is defined as starting from dump file > 1
 restart = (numfromfile(dumpfile) > 1)

 !if we are starting from scratch
 if (.not. restart) then

    vxyzu = 0.
    tprev = -1.
    if (id==master .and. iverbose >= 2) write(iprint,*) 'SETTING VELS TO ZERO'

#ifdef STIR_FROM_FILE
    call read_stirring_data_from_file(forcingfile, time, timeinfile)
#else
    ! Everyone, using the same seed, initialize the OU noises for the
    ! nmodes*6 components of the phases. Store seed in randseed
    ! afterward.
    call random_seed(size = st_seedLen)
    if (.not. allocated(st_randseed)) allocate(st_randseed(st_seedLen))
    if (id==master .and. iverbose >= 1) write(iprint,*) 'init_stir:  seed length = ', st_seedLen
    call st_ounoiseinit(st_nmodes*6, st_seed, st_OUvar, st_OUphases)
    call random_seed(get = st_randseed)
    if (id==master .and. iverbose >= 1) write(iprint,*) 'init_stir:  st_randseed = ', st_randseed
#endif

 else ! we are restarting from a checkpoint

#ifdef STIR_FROM_FILE
    call read_stirring_data_from_file(forcingfile, time, timeinfile)
#else
    ! Everyone, using the same seed, initialize the OU noises for the
    ! nmodes*6 components of the phases. Store seed in randseed
    ! afterward.
    call random_seed(size = st_seedLen)
    if (.not. allocated(st_randseed)) allocate(st_randseed(st_seedLen))
    if (id==master .and. iverbose >= 1) write(iprint,*) 'init_stir:  seed length = ', st_seedLen
    call st_ounoiseinit(st_nmodes*6, st_seed, st_OUvar, st_OUphases)
    call random_seed(get = st_randseed)
    if (id==master .and. iverbose >= 1) write(iprint,*) 'init_stir:  st_randseed = ', st_randseed

    call read_forcingdump(dumpfile,ierr)
    if (ierr /= 0) call fatal('init_forcing','error reading forcing file')
#endif

    !call get_parm_from_context(global_parm_context, "st_nmodes", st_nmodes)
    if (id==master) write(iprint,*) 'init_stir:  restarting...  st_nmodes      = ', st_nmodes

 endif ! restart

 if (id==master) then
    write(iprint,*) ' Initialized ',st_nmodes,' modes for stirring.'
    if (st_spectform == 0) write(iprint,*) ' spectral form        = ', st_spectform, ' (Band)'
    if (st_spectform == 1) write(iprint,*) ' spectral form        = ', st_spectform, ' (Parabola)'
    if (st_spectform == 2) write(iprint,*) ' spectral form        = ', st_spectform, ' (k^5/3)'
    write(iprint,*) ' solenoidal weight    = ', st_solweight
    write(iprint,*) ' st_solweightnorm     = ', st_solweightnorm
    write(iprint,*) ' stirring energy      = ', st_energy
    write(iprint,*) ' autocorrelation time = ', st_decay
    write(iprint,*) ' minimum wavenumber   = ', st_stirmin
    write(iprint,*) ' maximum wavenumber   = ', st_stirmax
    write(iprint,*) ' original random seed = ', st_seed
#ifdef CORRECT_MEAN_FORCE
    write(iprint,*) ' bulk motion correction is ON'
#else
    write(iprint,*) ' bulk motion correction is OFF'
#endif
 endif


 ! Then convert those into actual Fourier phases:
 call st_calcPhases()

 return

end subroutine init_forcing

!! NAME
!!
!!  write_options_forcing
!!
!! SYNOPSIS
!!
!!  write_options_forcing(iunit)
!!
!! DESCRIPTION
!!
!!     This routine writes the stirring parameters to the input file
!!
!! ARGUMENTS
!!
!!      iunit :: unit number for the input file which should already be open
!!
!! AUTHOR
!!
!!   Daniel Price, 2009
!!
!!***

subroutine write_options_forcing(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling forcing of turbulence'
 call write_inopt(istir,'istir','switch to turn stirring on or off at runtime',iunit)
#ifndef STIR_FROM_FILE
 call write_inopt(st_spectform,'st_spectform','spectral form of stirring',iunit)
 call write_inopt(st_stirmax,'st_stirmax','maximum stirring wavenumber',iunit)
 call write_inopt(st_stirmin,'st_stirmin','minimum stirring wavenumber',iunit)
 call write_inopt(st_energy,'st_energy','energy input/mode',iunit)
 call write_inopt(st_decay,'st_decay','correlation time for driving',iunit)
 call write_inopt(st_solweight,'st_solweight','solenoidal weight',iunit)
 call write_inopt(st_dtfreq,'st_dtfreq','frequency of stirring',iunit)
 call write_inopt(st_seed,'st_seed','random number generator seed',iunit)
#endif
 call write_inopt(st_amplfac,'st_amplfac','amplitude factor for stirring of turbulence',iunit)

end subroutine write_options_forcing

!! NAME
!!
!!  read_options_forcing
!!
!! SYNOPSIS
!!
!!  read_options_forcing(iunit)
!!
!! DESCRIPTION
!!
!!     This routine is called whilst reading the input file to extract runtime parameters
!!
!! ARGUMENTS
!!
!!      iunit :: unit number for the input file which should already be open
!!
!! AUTHOR
!!
!!   Daniel Price, 2009
!!
!!***

subroutine read_options_forcing(name,valstring,igot,igotall,ierr)
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: igot,igotall
 integer,          intent(out) :: ierr
 integer, save                :: ngot = 0
#ifdef STIR_FROM_FILE
 integer, parameter :: nrequired = 2
#else
 integer, parameter :: nrequired = 10
#endif

 igot = .true.
 select case(trim(name))
 case('istir')
    read(valstring,*,iostat=ierr) istir
#ifndef STIR_FROM_FILE
 case('st_spectform')
    read(valstring,*,iostat=ierr) st_spectform
 case('st_stirmax')
    read(valstring,*,iostat=ierr) st_stirmax
 case('st_stirmin')
    read(valstring,*,iostat=ierr) st_stirmin
 case('st_energy')
    read(valstring,*,iostat=ierr) st_energy
 case('st_decay')
    read(valstring,*,iostat=ierr) st_decay
 case('st_solweight')
    read(valstring,*,iostat=ierr) st_solweight
 case('st_dtfreq')
    read(valstring,*,iostat=ierr) st_dtfreq
 case('st_seed')
    read(valstring,*,iostat=ierr) st_seed
#endif
 case('st_amplfac')
    read(valstring,*,iostat=ierr) st_amplfac
 case default
    igot = .false.
 end select

 if (igot) ngot = ngot + 1

 igotall = .false.
 if (ngot >= nrequired) igotall = .true.

end subroutine read_options_forcing


!! NAME
!!
!!  write_forcingdump
!!
!! SYNOPSIS
!!
!!  write_forcingdump(tdump,dumpfile)
!!
!! DESCRIPTION
!!
!!     This routine writes the stirring information to file
!!     alongside code dumps so the run can be restarted
!!
!! ARGUMENTS
!!
!!      tdump :: time
!!   dumpfile :: name of corresponding Phantom dump file.
!!
!! AUTHOR
!!
!!   Daniel Price, 2008
!!
!!***


subroutine write_forcingdump(tdump,dumpfile)
 use io, only:ifdump,error,iprint
 real,             intent(in) :: tdump
 character(len=*), intent(in) :: dumpfile
 integer :: ierr
 character(len=len(dumpfile)+6) :: filename

 filename = trim(dumpfile)//'.str'

 write(iprint,*) 'writing FORCING INFO to '//trim(filename)
 open(unit=ifdump,file=filename,status='replace',form='formatted',iostat=ierr)

 if (ierr /= 0) then
    call error('write_forcingdump','could not open forcing file for writing')
 else
    write(ifdump,*) tdump,tprev,st_seed
    write(ifdump,*) st_aka
    write(ifdump,*) st_akb
    write(ifdump,*) st_OUphases
    write(ifdump,*) st_ampl
    !write(ifdump,*) tdump,t_turb,tprev,kmax,ifkmax,iseed
    !write(ifdump,*) prank(1:3,1:ifkmax)
    !write(ifdump,*) frank(1:3,1:ifkmax)
 endif

 close(unit=ifdump)

end subroutine write_forcingdump

!! NAME
!!
!!  write_forcingdump
!!
!! SYNOPSIS
!!
!!  write_forcingdump(tdump,dumpfile)
!!
!! DESCRIPTION
!!
!!     Reads forcing arrays from file on restart
!!
!! ARGUMENTS
!!
!!   dumpfile (IN)  :: name of corresponding Phantom dumpfile
!!       ierr (OUT) :: error code
!!
!! AUTHOR
!!
!!   Daniel Price, 2008
!!
!!***

subroutine read_forcingdump(dumpfile,ierr)
 use io, only:ifdumpread,iprint,fatal
 character(len=*), intent(in)  :: dumpfile
 integer,          intent(out) :: ierr
 real :: tdump
 character(len=len(dumpfile)+6) :: filename

 filename = trim(dumpfile)//'.str'
 write(iprint,*) 'reading FORCING INFO from '//trim(filename)

 open(unit=ifdumpread,file=filename,status='old',form='formatted',iostat=ierr)

 if (ierr /= 0) then
    call fatal('read_forcingdump','error opening '//trim(filename)//' for restart')
 else
    read(ifdumpread,*) tdump,tprev,st_seed
    !read(ifdumpread,*) tdump,t_turb,tprev,kmaxtemp,ifkmaxtemp,iseed
    write(iprint,*) ' Time in forcing dump = ',tdump,' tprev = ',tprev, &
                    ' iseed = ',st_seed
    !if (ifkmaxtemp /= ifkmax .or. kmaxtemp /= kmax) then
    !   write(iprint,*) 'ERROR IN RESTART: kmax changed : ', &
    !                   ifkmax,ifkmaxtemp,kmax,kmaxtemp
    !   stop
    !endif
    read(ifdumpread,*) st_aka
    read(ifdumpread,*) st_akb
    read(ifdumpread,*) st_OUphases
    read(ifdumpread,*) st_ampl
 endif

 close(unit=ifdumpread)

end subroutine read_forcingdump


!! NAME
!!
!!  st_calcPhases
!!
!! SYNOPSIS
!!
!!  st_calcPhases()
!!
!! DESCRIPTION
!!
!!     This routine updates the stirring phases from the OU phases.
!!     It copies them over and applies the projection operation.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! AUTHOR
!!
!!   Christoph Federrath, 2008
!!
!!   modified for use in PhantomSPH by Daniel Price 2008
!!***


subroutine st_calcPhases()

 !use dBase, ONLY: ndim
 !use Stir_data, ONLY : st_nmodes, st_mode, st_aka, st_akb, st_OUphases, st_solweight


 real                 :: ka, kb, kk, diva, divb, curla, curlb
 integer              :: i,j
 logical, parameter   :: Debug = .false.

 do i = 1, st_nmodes
    ka = 0.0
    kb = 0.0
    kk = 0.0
    do j=1, ndim
       kk = kk + st_mode(j,i)*st_mode(j,i)
       ka = ka + st_mode(j,i)*st_OUphases(6*(i-1)+2*(j-1)+1+1)
       kb = kb + st_mode(j,i)*st_OUphases(6*(i-1)+2*(j-1)+0+1)
    enddo
    do j=1, ndim

       diva  = st_mode(j,i)*ka/kk
       divb  = st_mode(j,i)*kb/kk
       curla = (st_OUphases(6*(i-1)+2*(j-1) + 0 + 1) - divb)
       curlb = (st_OUphases(6*(i-1)+2*(j-1) + 1 + 1) - diva)

       st_aka(j,i) = st_solweight*curla+(1.0-st_solweight)*divb
       st_akb(j,i) = st_solweight*curlb+(1.0-st_solweight)*diva

! purely compressive
!         st_aka(j,i) = st_mode(j,i)*kb/kk
!         st_akb(j,i) = st_mode(j,i)*ka/kk

! purely solenoidal
!         st_aka(j,i) = bjiR - st_mode(j,i)*kb/kk
!         st_akb(j,i) = bjiI - st_mode(j,i)*ka/kk

       if (Debug) then
          print *, 'st_mode(dim=',j,',mode=',i,') = ', st_mode(j,i)
          print *, 'st_aka (dim=',j,',mode=',i,') = ', st_aka(j,i)
          print *, 'st_akb (dim=',j,',mode=',i,') = ', st_akb(j,i)
       endif

    enddo
 enddo

 return

end subroutine st_calcPhases


!! NAME
!!
!!  st_ounoiseinit
!!
!! SYNOPSIS
!!
!!  st_ounoiseinit(integer,intent (IN)  :: vectorlength,
!!                 integer,intent (IN)  :: iseed,
!!                 real,intent (IN)  :: variance,
!!                 real,intent (INOUT)  :: vector)
!!
!! DESCRIPTION
!!
!!
!! Subroutine initializes a vector of real numbers to be used as a
!!   starting point for the Ornstein-Uhlenbeck, or "colored noise"
!!   generation sequence. Note that this should be invoked once at
!!   the very start of the program; "vector" values and the random
!!   seed should be checkpointed to ensure continuity across restarts.
!!
!! The length of the vector is specified in "vectorlength", and the
!!   random seed to be used is passed in as "seed". The sequence
!!   is initialized using Gaussian values with variance "variance".
!!
!! Please refer to st_ounoiseupdate for further details on algorithm.
!!
!! ARGUMENTS
!!
!!   vectorlength : lenght of the vector
!!
!!   iseed :  input seed for random number generator
!!
!!   variance : variance of distribution
!!
!!   vector : storage for starting vector
!!
!!
!!***


subroutine st_ounoiseinit (vectorlength, iseed, variance, vector)

 !use Stir_data, ONLY : st_randseed


 integer, intent(IN)    :: vectorlength, iseed
 real,    intent(IN)    :: variance
 real,    intent(INOUT) :: vector (vectorlength)
 real                    :: grnval
 integer                 :: i

 !... Initialize pseudorandom sequence with random seed.

 st_randseed = iseed

 call random_seed (put = st_randseed)

 do i = 1, vectorlength
    call st_grn (grnval)
    vector (i) = grnval * variance
 enddo

 return

end subroutine st_ounoiseinit

!! NAME
!!
!!  st_ounoiseupdate
!!
!! SYNOPSIS
!!
!!  st_ounoiseupdate(integer, intent (IN)  :: vectorlength,
!!                   real, intent (INOUT)  :: vector,
!!                   real, intent (IN)     :: variance,
!!                   real, intent (IN)     :: dt,
!!                   real, intent (IN)     :: ts)
!!
!! DESCRIPTION
!!
!! Subroutine updates a vector of real values according to an algorithm
!!   that generates an Ornstein-Uhlenbeck, or "colored noise" sequence.
!!
!! The sequence x_n is a Markov process that takes the previous value,
!!   weights by an exponential damping factor with a given correlation
!!   time "ts", and drives by adding a Gaussian random variable with
!!   variance "variance", weighted by a second damping factor, also
!!   with correlation time "ts". For a timestep of dt, this sequence
!!   can be written as :
!!
!!     x_n+1 = f x_n + sigma * sqrt (1 - f**2) z_n
!!
!! where f = exp (-dt / ts), z_n is a Gaussian random variable drawn
!! from a Gaussian distribution with unit variance, and sigma is the
!! desired variance of the OU sequence. (See Bartosch, 2001).
!!
!! The resulting sequence should satisfy the properties of zero mean,
!!   and stationary (independent of portion of sequence) RMS equal to
!!   "variance". Its power spectrum in the time domain can vary from
!!   white noise to "brown" noise (P (f) = const. to 1 / f^2).
!!
!! References :
!!    Bartosch, 2001
!! http://octopus.th.physik.uni-frankfurt.de/~bartosch/publications/IntJMP01.pdf
!!   Finch, 2004
!! http://pauillac.inria.fr/algo/csolve/ou.pdf
!!         Uhlenbeck & Ornstein, 1930
!! http://prola.aps.org/abstract/PR/v36/i5/p823_1
!!
!!
!! ARGUMENTS
!!
!!   vectorlength : length of vector to be updated
!!
!!   vector :       vector to be updated
!!
!!   variance :     variance of the distributio
!!
!!   dt :           timestep
!!
!!   ts :           correlation time
!!
!!***


subroutine st_ounoiseupdate (vectorlength, vector, variance, dt, ts)


 real,    intent(IN)    :: variance, dt, ts
 integer, intent(IN)    :: vectorlength
 real,    intent(INOUT) :: vector (vectorlength)
 real                              :: grnval, damping_factor
 integer                           :: i

 damping_factor = exp (-dt/ts)

 do i = 1, vectorlength
    call st_grn (grnval)
    vector (i) = vector (i) * damping_factor + variance *   &
          sqrt (1.0 - damping_factor**2) * grnval
 enddo

 return

end subroutine st_ounoiseupdate

!! NAME
!!
!!  st_calcAccel
!!
!! SYNOPSIS
!!
!!  st_calcAccel()
!!
!! DESCRIPTION
!!
!!   Computes components of the zone-averaged forceitational
!!   acceleration.
!!
!! ARGUMENTS
!!
!!   xyzh,fxyzu
!!
!! AUTHOR: Christoph Federrath
!! modified for use in PhantomSPH by Daniel Price 2008
!!
!!***

subroutine st_calcAccel(npart,xyzh,fxyzu)
!  use Stir_data, ONLY : st_maxmodes, st_mode, st_nmodes, st_aka, st_akb, st_solweightnorm, st_ampl
!  use part, only:massoftype
#ifdef IND_TIMESTEPS
 use part, only:iactive,iphase
#endif
#ifdef CORRECT_MEAN_FORCE
 use mpiutils, only:reduceall_mpi
#endif
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: fxyzu(:,:)
 integer                   :: i,m
 real                      :: ampl,kxxi,kyyi,kzzi,kdotx,fxi,fyi,fzi
 real :: xyzi(ndim),fmean(ndim)
 real                      :: realtrigterms, imtrigterms
 logical, parameter        :: Debug = .true.

!==============================================================================

 fmean(:) = 0.

 !! this is the critical loop wrt to computational resources
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(xyzh,fxyzu,npart,st_nmodes,st_mode,st_ampl) &
#ifdef IND_TIMESTEPS
 !$omp shared(iphase) &
#endif
 !$omp shared(st_aka,st_akb,st_solweightnorm,st_amplfac) &
 !$omp private(i,xyzi,fxi,fyi,fzi,kxxi,kyyi,kzzi,kdotx,ampl) &
 !$omp private(realtrigterms,imtrigterms) &
#ifdef CORRECT_MEAN_FORCE
 !$omp reduction(+:fmean) &
#endif
 !$omp private(m)
 do i=1,npart
#ifndef CORRECT_MEAN_FORCE
#ifdef IND_TIMESTEPS
    if (iactive(iphase(i))) then
#endif
#endif
       xyzi(1:ndim) = xyzh(1:ndim,i)
       fxi=0.
       fyi=0.
       fzi=0.

       do m = 1, st_nmodes
          kxxi = st_mode(1,m)*xyzi(1)
          kyyi = st_mode(2,m)*xyzi(2)
          kzzi = st_mode(3,m)*xyzi(3)
          kdotx = kxxi + kyyi + kzzi
          ampl = st_ampl(m) ! was multiplied by 2.0

          !  these are the real and imaginary parts, respectively, of
          !     e^{ i \vec{k} \cdot \vec{x} }
          !          = cos(kx*x + ky*y + kz*z) + i sin(kx*x + ky*y + kz*z)

          realtrigterms = cos(kdotx)
          imtrigterms = sin(kdotx)

          fxi = fxi + ampl*(st_aka(1,m)*realtrigterms - st_akb(1,m)*imtrigterms)
          fyi = fyi + ampl*(st_aka(2,m)*realtrigterms - st_akb(2,m)*imtrigterms)
          fzi = fzi + ampl*(st_aka(3,m)*realtrigterms - st_akb(3,m)*imtrigterms)
       enddo
       !
       !--here we have added an extra parameter st_amplfac to enable
       !  a reduction (or increase) in the stirring amplitude even after
       !  the forcing pattern has been initialised or read from a file
       !
       fxi = 2.*st_amplfac*st_solweightnorm*fxi  ! multiplication by 2.0 moved here
       fyi = 2.*st_amplfac*st_solweightnorm*fyi
       fzi = 2.*st_amplfac*st_solweightnorm*fzi

#ifdef CORRECT_MEAN_FORCE
#ifdef IND_TIMESTEPS
       if (iactive(iphase(i))) then
#endif
#endif
          fxyzu(1,i) = fxi
          fxyzu(2,i) = fyi
          fxyzu(3,i) = fzi
#ifdef IND_TIMESTEPS
       endif
#endif
#ifdef CORRECT_MEAN_FORCE
       fmean(1) = fmean(1) + fxi
       fmean(2) = fmean(2) + fyi
       fmean(3) = fmean(3) + fzi
#endif
    enddo
    !$omp end parallel do

#ifdef CORRECT_MEAN_FORCE
!!  fmean = reduceall_mpi('+',fmean)
    fmean(:) = fmean(:)/real(npart)

!!$omp master
    if (Debug) then
!     print *, 'stir:  mass      = ', massoftype(1)*npart
!     print *, 'stir:  xMomentum = ', xmom(1)
!     print *, 'stir:  yMomentum = ', xmom(2)
!     print *, 'stir:  zMomentum = ', xmom(3)
       print *, 'stir:  xforce_mean = ', fmean(1)
       print *, 'stir:  yforce_mean = ', fmean(2)
       print *, 'stir:  zforce_mean = ', fmean(3)
    endif
!!$omp end master
!
!--correction for bulk motion: note that, unlike in Christoph's code
!  we do not do the time integration here - rather we return the
!  force as part of the total force.
!
    !$omp parallel do private(i) shared(fmean) schedule(static)
    do i=1,npart
#ifdef IND_TIMESTEPS
       if (iactive(iphase(i))) then
#endif
!     vxyzu(1:ndim,i) = vxyzu(1:ndim,i) - xmom(1:ndim)
          fxyzu(1:ndim,i) = fxyzu(1:ndim,i) - fmean(1:ndim)
#ifdef IND_TIMESTEPS
       endif
#endif
    enddo
    !$omp end parallel do
#endif


    return

 end subroutine st_calcAccel

!! NAME
!!
!!  Stir
!!
!! SYNOPSIS
!!
!!  Stir(real(IN) t,integer(IN) npart, real(IN) xyzh,vxyzu,fxyzu)
!!       real(IN) :: dt)
!!
!! DESCRIPTION
!!   Apply the stirring opperator on the list of blocks provided as input
!!
!! ARGUMENTS
!!   t            : the current time
!!   npart        : number of SPH particles
!!   xyzh         : positions of SPH particles
!!   vxyzu        : velocities of SPH particles
!!   fxyzu        : forces for SPH particles
!!   dt           : the current timestep
!!
!! NOTES
!!
!!     modified by Christoph Federrath 2008
!!     modified for use in PhantomSPH by Daniel Price 2008
!!
!!***

 subroutine forceit(t,npart,xyzh,vxyzu,fxyzu)
  real,    intent(in)    :: t
  integer, intent(in)    :: npart
  real,    intent(in)    :: xyzh(:,:)
  real,    intent(inout) :: vxyzu(:,:)
  real,    intent(out)   :: fxyzu(:,:)

  logical   :: update_accel = .true.
#ifdef STIR_FROM_FILE
  real :: tinfile
#endif

  logical, parameter :: Debug = .false.
!!===================================================================

  ! Only update acceleration every st_freq timesteps
  update_accel = .false.
  !nstep = dBasePropertyInteger('CurrentStepNumber')
  !if ( (nstep == 1) .or. (mod (nstep, st_freq) == 0) ) then
  if (t > (tprev + st_dtfreq)) then
     tprev = st_dtfreq*int(t/st_dtfreq)  ! round to last full dtfreq
     update_accel = .true.
#ifdef STIR_FROM_FILE
     call read_stirring_data_from_file(forcingfile,t,tinfile)
     !if (id==master .and. iverbose >= 2) print*,' got new accel, tinfile = ',tinfile
#endif
  endif

  if (Debug) print *, 'stir:  stirring start'

  call st_calcAccel(npart,xyzh,fxyzu)

#ifndef STIR_FROM_FILE
  if (update_accel) then
     if (Debug) print*,'updating accelerations...'
     call st_ounoiseupdate(6*st_nmodes, st_OUphases, st_OUvar, st_dtfreq, st_decay)
     call st_calcPhases()
     !! Store random seed in memory for later checkpoint.
     call random_seed (get = st_randseed)
  endif
#endif

  if (Debug) print *, 'stir:  stirring end'

  return

 end subroutine forceit

!! NAME
!!
!!  st_grn
!!
!! SYNOPSIS
!!
!!  st_grn(real,intent (OUT)  :: grnval)
!!
!! DESCRIPTION
!!
!!  Subroutine draws a number randomly from a Gaussian distribution
!!    with the standard uniform distribution function "random_number"
!!    using the Box-Muller transformation in polar coordinates. The
!!    resulting Gaussian has unit variance.
!!
!!  Reference : Numerical Recipes, section 7.2.
!!
!!
!!
!! ARGUMENTS
!!
!!   grnval : the number drawn from the distribution
!!
!!***
 subroutine st_grn (grnval)


  real, intent(OUT) :: grnval
  real              :: pi, r1, r2, g1

  pi = 4. * atan (1.)

  r1 = 0.; r2 = 0;

  call random_number (r1)
  call random_number (r2)
  g1 = sqrt (2. * log (1. / r1) ) * cos (2. * pi * r2)

  grnval = g1

  return

 end subroutine st_grn

#ifdef STIR_FROM_FILE
!! NAME
!!
!!  read_stirring_data_from_file
!!
!! SYNOPSIS
!!
!!  read_stirring_data_from_file()
!!
!! DESCRIPTION
!!  reads the stirring data necessary to construct
!!  the physical force field from file
!!
!! ARGUMENTS
!!  infile : filename of file containing the forcing data
!!  time   : physical time for which the modes/aka/akb are to be updated
!!
!! AUTHOR: Christoph Federrath
!!
!!***

 subroutine read_stirring_data_from_file(infile, time, timeinfile)
  use io,        only:id,master,fatal
  use datafiles, only:find_phantom_datafile

  character(len=*), intent(in)  :: infile
  real,             intent(in)  :: time
  real,             intent(out) :: timeinfile

  logical, parameter  :: Debug = .false.

  integer             :: nsteps, istep, istepfile, igetstep, ierr
  real                :: end_time
  character(120)      :: my_file

  my_file = find_phantom_datafile(infile,'forcing')

  open (unit=42, file=my_file, iostat=ierr, status='OLD', action='READ', &
        access='SEQUENTIAL', form='UNFORMATTED')
  ! header contains number of times and number of modes, end time, autocorrelation time, ...
  if (ierr==0) then
     if (Debug) write (*,'(A)') 'reading header...'
     read (unit=42) nsteps, st_nmodes, end_time, st_decay, st_energy, st_solweight, &
                    st_solweightnorm, st_stirmin, st_stirmax, st_spectform
     if (Debug) write (*,'(A)') '...finished reading header'
  else
     call fatal('read_stirring_data_from_file','could not open '//trim(my_file)//' for read ')
  endif

  ! these are in the global contex
  st_dtfreq = end_time/nsteps
  igetstep = floor(time/st_dtfreq)

  do istep = 0, nsteps
     if (Debug) write (*,'(A,I6)') 'step = ', istep
     read (unit=42) istepfile, timeinfile, &
                     st_mode    (:, 1:  st_nmodes), &
                     st_aka     (:, 1:  st_nmodes), &
                     st_akb     (:, 1:  st_nmodes), &
                     st_ampl    (   1:  st_nmodes), &
                     st_OUphases(   1:6*st_nmodes)

     if (istep /= istepfile) write(*,'(A,I6)') 'read_stirring_data_from_file: something wrong! step = ', istep
     if (igetstep==istep) then
        if (id==master) then
           write(*,'(A,I6,2(A,E20.6))') 'read_stirring_data_from_file: read new forcing pattern, stepinfile = ',&
                  istep,' , time = ', time, ' , time in stirring table = ', timeinfile
        endif
        close (unit=42)
        exit ! the loop
     endif
  enddo

  return

 end subroutine read_stirring_data_from_file
#endif

end module forcing
