!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!  OWNER: Christophe Pinte
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, eos, io, mcfost2phantom, options, part, timestep,
!    units
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
                          massoftype,xyzmh_ptmass,luminosity,igas,&
                          grainsize,graindens,ndusttypes,&
                          do_radiation,radiation,ithick,maxirad,ikappa
 use units,          only:umass,utime,udist
 use io,             only:fatal,warning
 use dim,            only:use_dust,lightcurve,maxdusttypes
 use eos,            only:temperature_coef,gmw,gamma
 use timestep,       only:dtmax
 use options,        only:use_dustfrac,use_mcfost,use_Voronoi_limits_file,Voronoi_limits_file, &
                          use_mcfost_stellar_parameters
 use physcon,        only:cm,gram

 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(in)    :: particlemass,time
 real,             intent(inout) :: vxyzu(:,:)

 logical, save   :: init_mcfost = .false.
 real            :: mu_gas,factor,T_to_u
 real(kind=4), dimension(npart)    :: Tdust, n_packets
 integer         :: ierr,ntypes,dustfluidtype,ilen,nlum,i
 integer(kind=1) :: itype(maxp)
 logical         :: compute_Frad
 real(kind=8), dimension(6), save            :: SPH_limits
 real,         dimension(:),     allocatable :: dudt
 real,    parameter :: Tdefault = 1.
 logical, parameter :: write_T_files = .false. ! ask mcfost to write fits files with temperature structure
 integer, parameter :: ISM = 2 ! ISM heating : 0 -> no ISM radiation field, 1 -> ProDiMo, 2 -> Bate & Keto
 character(len=len(dumpfile) + 20) :: mcfost_para_filename
 integer :: skippedNum = 0

 if (use_mcfost) then
    write(*,*) "Calling mcfost"
    if (.not.init_mcfost) then
       ilen = index(dumpfile,'_',back=.true.) ! last position of the '_' character
       mcfost_para_filename = dumpfile(1:ilen-1)//'.para'
       call init_mcfost_phantom(mcfost_para_filename,ndusttypes,use_Voronoi_limits_file,Voronoi_limits_file,SPH_limits,ierr, &
            fix_star = use_mcfost_stellar_parameters)
       if (ierr /= 0) call fatal('mcfost-phantom','error in init_mcfost_phantom')
       init_mcfost = .true.
    endif

    ntypes = get_ntypes(npartoftype)
    if (maxphase==maxp) then
       itype = iamtype(iphase)
    else
       itype(:) = 1
    endif
    if (.not. use_dust) then
       ndusttypes = 0
    endif
    if (use_dustfrac) then
       dustfluidtype = 1
    else
       dustfluidtype = 2
    endif

    nlum = npart
    allocate(dudt(nlum))
    if (lightcurve) then
       dudt(1:nlum) = luminosity(1:nlum)
    else
       dudt(1:nlum) = vxyzu(4,1:nlum) * massoftype(igas) / dtmax
    endif

    factor = 1.0/(temperature_coef*gmw*(gamma-1))
    ! this this the factor needed to compute u^(n+1)/dtmax from temperature
    T_to_u = factor * massoftype(igas) /dtmax

    call run_mcfost_phantom(&
         npart,nptmass,ntypes,ndusttypes,dustfluidtype,npartoftype,maxirad,&
         xyzh,vxyzu,radiation,ikappa,&
         itype,grainsize,graindens,dustfrac,massoftype,&
         xyzmh_ptmass,hfact,umass,utime,udist,nlum,dudt,compute_Frad,SPH_limits,Tdust,&
         n_packets,mu_gas,ierr,write_T_files,ISM,T_to_u)
    print*, 'just after mcfost'

    write(*,*) ''
    write(*,*) 'Minimum temperature = ', minval(Tdust, mask=(Tdust > 0.))
    write(*,*) 'Maximum temperature = ', maxval(Tdust)
    write(*,*) ''

    ! set thermal energy
    if (do_radiation) then
      do i=1,npart
         if ((Tdust(i) > 1.).and.(radiation(ithick,i) < 0.5)) then
            vxyzu(4,i) = Tdust(i) * factor
         endif
         ! radiation(ikappa,i) = radiation(ikappa,i)*(cm**2/gram)/(udist**2/umass)
         print*, i, ')   [kappa] physical=',radiation(ikappa,i),&
            '[kappa] code=',radiation(ikappa,i)*(cm**2/gram)/(udist**2/umass),&
            '[1] code=', 1*(cm**2/gram)/(udist**2/umass)
         read*
         ! if (radiation(ikappa,i) == 0) read*
         ! if (radiation(ikappa,i) > 1e6) skippedNum = skippedNum + 1
         ! if (i<100) then
         !    print*, vxyzu(4,i), Tdust(i), radiation(ithick,i)
         !    read*
         ! endif
         ! else left the temperature as it was before
      enddo
      ! print*, "number of particles with old kappa: ", skippedNum
    else
      do i=1,npart
         if (Tdust(i) > 1.) then
            vxyzu(4,i) = Tdust(i) * factor
         else
            ! if mcfost doesn't return a temperature set it to Tdefault
            vxyzu(4,i) = Tdefault * factor
         endif
      enddo
    endif

    if (allocated(dudt)) deallocate(dudt)

    write(*,*) "End of analysis mcfost"
 endif ! use_mcfost

 return

end subroutine do_analysis

end module
