!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!   Analysis routine to call MCFOST code to perform post-processing
!   radiation transport
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, dust, eos, io, mcfost2phantom, part, units
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'mcfost'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use mcfost2phantom, only:init_mcfost_phantom,run_mcfost_phantom
 use part,           only:massoftype,iphase,dustfrac,hfact,npartoftype,&
                          get_ntypes,iamtype,maxphase,maxp,idust,nptmass,&
                          massoftype,xyzmh_ptmass,luminosity
 use units,          only:umass,utime,udist
 use io,             only:fatal
 use dim,            only:use_dust,lightcurve
 use dust,           only:grainsize,graindens
 use eos, only : temperature_coef, gmw, gamma

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:)
 real,             intent(in) :: particlemass,time
 real,           intent(inout):: vxyzu(:,:)
 logical, save :: init_mcfost = .false.
 real    :: mu_gas, factor
 real(kind=4) :: Tdust(npart)
 real    :: grain_size(1)
 integer :: ierr,ntypes,ndusttypes,dustfluidtype,ilen,nlum, i
 integer(kind=1) :: itype(maxp)
 character(len=len(dumpfile) + 20) :: mcfost_para_filename
 logical :: compute_Frad
 real,dimension(6) :: SPH_limits
 real(kind=4),dimension(:,:,:),allocatable :: Frad
 real,dimension(:),allocatable :: dudt

 real, parameter :: Tdefault = 1.
 logical, parameter :: write_T_files = .true. ! ask mcfost to write fits files with temperature structure
 logical, parameter :: ISM = .true. ! turn on ISM heating


 if (.not.init_mcfost) then
    ilen = index(dumpfile,'_')
    mcfost_para_filename = dumpfile(1:ilen-1)//'.para'
    call init_mcfost_phantom(mcfost_para_filename, ierr) !,  np, nptmass, ntypes, ndusttypes, npoftype)
    if (ierr /= 0) call fatal('mcfost-phantom','error in init_mcfost_phantom')
    init_mcfost = .true.
 endif

 ntypes = get_ntypes(npartoftype)
 if (maxphase==maxp) then
    itype = iamtype(iphase)
 else
    itype(:) = 1
 endif
 if (use_dust) then
    ndusttypes = 1
 else
    ndusttypes = 0
 endif
 if (npartoftype(idust) > 0) then
    dustfluidtype = 2
 else
    dustfluidtype = 1
 endif
 grain_size(:) = grainsize
 SPH_limits = 0.
 nlum = 0
 if (lightcurve) nlum = npart
 allocate(dudt(nlum))
 if (lightcurve) dudt(1:nlum) = luminosity(1:nlum)
 allocate(Frad(3,ndusttypes,npart))

 call run_mcfost_phantom(npart,nptmass,ntypes,ndusttypes,dustfluidtype,&
   npartoftype,xyzh,vxyzu,itype,grain_size,graindens,dustfrac,massoftype,&
   xyzmh_ptmass,hfact,umass,utime,udist,nlum,dudt,compute_Frad,SPH_limits,Tdust,&
   Frad,mu_gas,ierr,write_T_files,ISM)
 !print*,' mu_gas = ',mu_gas

 write(*,*) ''
 write(*,*) 'Minimum temperature = ', minval(Tdust, mask=(Tdust > 0.))
 write(*,*) 'Maximum temperature = ', maxval(Tdust)
 write(*,*) ''

 ! set thermal energy
 factor = 1.0/(temperature_coef*gmw*(gamma-1))
 do i=1,npart
    if (Tdust(i) > 0.) then
       vxyzu(4,i) = Tdust(i) * factor
    else
       ! if mcfost doesn't return a temperature set it to Tdefault
       vxyzu(4,i) = Tdefault * factor
    endif
 enddo

 return

end subroutine do_analysis

end module
