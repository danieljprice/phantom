!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to call MCFOST code to perform post-processing
!   radiation transport
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: deriv, dim, energies, eos, growth, io, mcfost2phantom,
!   omp_lib, options, part, physcon, units
!
 use omp_lib

 implicit none
 character(len=20), parameter, public :: analysistype = 'mcfost'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use mcfost2phantom, only:init_mcfost_phantom,&
                             run_mcfost_phantom,&
                             diffusion_opacity,&
                             reset_mcfost_phantom
 use part,           only:massoftype,iphase,dustfrac,hfact,npartoftype,&
                             get_ntypes,iamtype,maxphase,maxp,idust,nptmass,&
                             massoftype,xyzmh_ptmass,vxyz_ptmass,luminosity,igas,&
                             grainsize,graindens,ndusttypes,rad,radprop,&
                             rhoh,ikappa,iradxi,ithick,inumph,drad,ivorcl,eos_vars,itemp
 use units,          only:umass,utime,udist,get_radconst_code
 use io,             only:fatal,iprint
 use dim,            only:use_dust,lightcurve,maxdusttypes,use_dustgrowth,do_radiation
 use eos,            only:temperature_coef,gmw,gamma
 use options,        only:use_dustfrac,use_mcfost,use_Voronoi_limits_file,Voronoi_limits_file, &
                             use_mcfost_stellar_parameters, mcfost_computes_Lacc, mcfost_uses_PdV,&
                             mcfost_keep_part, ISM, mcfost_dust_subl
 use physcon,        only:cm,gram,c,steboltz

 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(in)    :: particlemass,time
 real,             intent(inout) :: vxyzu(:,:)

 logical, save   :: init_mcfost = .false., isinitial = .true.
 real            :: mu_gas,factor
 real(kind=4)    :: Tdust(npart),n_packets(npart)
 integer         :: ierr,ntypes,dustfluidtype,ilen,nlum,i,nerr
 integer(kind=1) :: itype(maxp)
 logical         :: compute_Frad
 logical         :: ISM_heating = .false.
 real(kind=8), dimension(6), save            :: SPH_limits
 real,         dimension(:),     allocatable :: dudt
 real,    parameter :: Tdefault = 1.
 logical, parameter :: write_T_files = .false. ! ask mcfost to write fits files with temperature structure
 character(len=len(dumpfile) + 20) :: mcfost_para_filename
 real :: a_code,rhoi,pmassi,Tmin,Tmax,default_kappa,kappa_diffusion

 if (.not. use_mcfost) return

 if (use_dustgrowth) then
    write(*,*) "Converting to fake multi large grains"
    call growth_to_fake_multi(npart)
 endif

 if (ISM > 0) then
    ISM_heating = .true.
 endif

 if (.not.init_mcfost) then
    ilen = index(dumpfile,'_',back=.true.) ! last position of the '_' character
    mcfost_para_filename = dumpfile(1:ilen-1)//'.para'
    call init_mcfost_phantom(mcfost_para_filename,ndusttypes,use_Voronoi_limits_file,&
         Voronoi_limits_file,SPH_limits,ierr, fix_star = use_mcfost_stellar_parameters, &
         turn_on_Lacc = mcfost_computes_Lacc, keep_particles = mcfost_keep_part, &
         use_ISM_heating = ISM_heating, turn_on_dust_subl = mcfost_dust_subl)
    if (ierr /= 0) call fatal('mcfost-phantom','error in init_mcfost_phantom')
    init_mcfost = .true.
 endif

 ntypes = get_ntypes(npartoftype)
 if (maxphase==maxp) then
    itype = iamtype(iphase)
 else
    itype(:) = 1
 endif
 if (.not. use_dust) ndusttypes = 0
 if (use_dustfrac) then
    dustfluidtype = 1
 else
    dustfluidtype = 2
 endif

 if (lightcurve .and. mcfost_uses_PdV) then
    nlum = npart
 else
    nlum =  0
 endif
 allocate(dudt(nlum))
 if (lightcurve) then
    dudt(1:nlum) = luminosity(1:nlum)
 else
    dudt(1:nlum) = 0.
 endif

 if (gamma <= 1.) then
    write(*,*) 'WARNING: gamma = 1 but should be > 1 for phantom+mcfost'
    write(*,*) 'RESETTING GAMMA = 5/3'
    gamma = 5./3.
 endif
 factor = 1.0/(temperature_coef*gmw*(gamma-1))

 !-- calling mcfost to get Tdust
 call run_mcfost_phantom(npart,nptmass,ntypes,ndusttypes,dustfluidtype,&
         npartoftype,xyzh,vxyzu,itype,grainsize,graindens,dustfrac,massoftype,&
         xyzmh_ptmass,vxyz_ptmass,hfact,umass,utime,udist,nlum,dudt,compute_Frad,SPH_limits,Tdust,&
         n_packets,mu_gas,ierr,write_T_files,ISM,eos_vars(itemp,:))

 Tmin = minval(Tdust, mask=(Tdust > 1.))
 Tmax = maxval(Tdust)
 write(*,*) ''
 write(*,*) 'Minimum temperature = ', Tmin
 write(*,*) 'Maximum temperature = ', Tmax
 write(*,*) ''

 if (use_mcfost .and. use_dustgrowth) then
    write(*,*) "Converting back to normal"
    call back_to_growth(npart)
 endif

 if (do_radiation) then
    a_code = get_radconst_code()
    pmassi = massoftype(igas)

    radprop(inumph,:) = 0.
    if (isinitial) then
       default_kappa = 0.5
    else
       default_kappa = 0.5*(maxval(radprop(ikappa,:))+minval(radprop(ikappa,:)))*(udist**2/umass)
    endif
    write(iprint,"(/,a,f4.2,' cm^2/g')") &
         ' -}+{- RADIATION: cutoff particles kappa = ',&
         default_kappa
    do i=1,npart
       if (maxphase==maxp) then
          if (iamtype(iphase(i)) /= igas) cycle
       endif
       radprop(inumph,i) = n_packets(i)
       if (radprop(inumph,i) > 1e2) then
          radprop(ithick,i) = 0.
       else
          radprop(ithick,i) = 1.
       endif
       if (isinitial.or.(radprop(ithick,i) < 0.5)) then
          ! initial run (t == 0) OR it has got enough info from mcfost
          ! => set new temperature for gas and radiation
          if (Tdust(i) > 1.) then
             vxyzu(4,i) = Tdust(i)*factor
             ! if the temperature is correct and set by mcfost
             ! => suppose we are at equilibrium
             rhoi = rhoh(xyzh(4,i),pmassi)
             rad(iradxi,i) = a_code*Tdust(i)**4.0/rhoi
             drad(iradxi,i) = 0
          else
             ! if I got no temperature from mcfost
             ! => lets try to handle the particle by SPH
             if (isinitial) then
                vxyzu(4,i) = ((Tmax-Tmin)*0.2)*factor
                rhoi = rhoh(xyzh(4,i),pmassi)
                rad(iradxi,i) = a_code*((Tmax-Tmin)*0.2)**4.0/rhoi
                drad(iradxi,i) = 0
             else
                radprop(ithick,i) = 1.
             endif
          endif
          ! else
          ! it is not initial run AND the particle has not got info from mcfost
          ! => temperature is old because of diffusion
       endif
       ! no matter what happend on previous stage, we need to set new
       ! diffusion coefficien with regards to a new/old temperatures
       if (radprop(ivorcl,i) > 0) then
          call diffusion_opacity(vxyzu(4,i)/factor,int(radprop(ivorcl,i)),kappa_diffusion)
          radprop(ikappa,i) = kappa_diffusion*(cm**2/gram)/(udist**2/umass)
       else
          radprop(ikappa,i) = default_kappa*(cm**2/gram)/(udist**2/umass)
       endif
    enddo

    write(iprint,"(/,a,f6.2,'%')") ' -}+{- RADIATION particles done by SPH = ',&
         100.*count(radprop(ithick,:)==1)/real(size(radprop(ithick,:)))
    isinitial = .false.
    nerr = 0
 else ! No diffusion approximation
    nerr = 0
    do i=1,npart
       if (Tdust(i) > 1.) then
          vxyzu(4,i) = Tdust(i) * factor
       else
          ! if mcfost doesn't return a temperature set it to Tdefault
          vxyzu(4,i) = Tdefault * factor
       endif
       if (vxyzu(4,i) <= 0. .or. vxyzu(4,i) > huge(0.)) then
          nerr = nerr + 1
          if (nerr < 10) write(*,*) 'ERROR: part ',i,' u = ',vxyzu(4,i),' Tdust was ',Tdust(i)
       endif
    enddo
 endif
 if (nerr > 0) write(*,*) 'ERROR: ** GOT ',nerr,' particles with invalid u (0 or Infinity), code will crash **'

 call reset_mcfost_phantom()
 if (allocated(dudt)) deallocate(dudt)
 write(*,*) "End of analysis mcfost"

end subroutine do_analysis

subroutine growth_to_fake_multi(npart)
 use growth, only:bin_to_multi,f_smax,size_max,b_per_dex
 use deriv,  only:get_derivs_global

 integer, intent(in)  :: npart

 !- bin sizes
 call bin_to_multi(b_per_dex,f_smax,size_max,verbose=.false.)

 !-- recompute density
 call get_derivs_global()

end subroutine growth_to_fake_multi

subroutine back_to_growth(npart)
 use part,     only:ndusttypes,ndustlarge,idust,massoftype,&
                    npartoftype,iamtype,iphase,idust,&
                    set_particle_type
 use energies, only:mdust
 integer, intent(in)    :: npart
 integer                :: i,j,ndustold,itype

 ndustold = sum(npartoftype(idust:))
 do i=1,npart
    itype = iamtype(iphase(i))
    if (itype > idust) then
       npartoftype(idust) = npartoftype(idust) + 1
       npartoftype(itype) = npartoftype(itype) - 1
       call set_particle_type(i,idust)
    endif
 enddo

 do j=2,ndusttypes
    if (npartoftype(idust+j-1) /= 0) write(*,*) 'ERROR! npartoftype ",idust+j-1 " /= 0'
    massoftype(idust+j-1)      = 0.
    mdust(idust+j-1)           = 0.
 enddo

 ndusttypes                    = 1
 ndustlarge                    = 1
 mdust(idust)                  = npartoftype(idust)*massoftype(idust)

 !- sanity checks for npartoftype
 if (npartoftype(idust) /= ndustold) then
    write(*,*) 'ERROR! npartoftype not conserved'
    write(*,*) npartoftype(idust), " <-- new vs. old --> ",ndustold
 endif

end subroutine back_to_growth

end module analysis
