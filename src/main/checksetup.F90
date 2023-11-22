!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module checksetup
!
! Perform sanity checks of the particle setup
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, boundary_dyn, centreofmass, dim, dust, eos,
!   externalforces, io, metric_tools, nicil, options, part, physcon,
!   ptmass_radiation, sortutils, timestep, units, utils_gr
!
 implicit none
 public :: check_setup

 private

contains

!------------------------------------------------------------------
!+
! This subroutine checks that the particle setup is sensible
!
! OUT:
!    nwarn  -- number of warnings triggered
!    nerror -- number of errors found
! IN:
!    restart -- whether the setup is at t=0, or from a later file
!+
!------------------------------------------------------------------
subroutine check_setup(nerror,nwarn,restart)
 use dim,  only:maxp,maxvxyzu,periodic,use_dust,ndim,mhd,use_dustgrowth, &
                do_radiation,n_nden_phantom,mhd_nonideal,do_nucleation,use_krome
 use part, only:xyzh,massoftype,hfact,vxyzu,npart,npartoftype,nptmass,gravity, &
                iphase,maxphase,isetphase,labeltype,igas,h2chemistry,maxtypes,&
                idust,xyzmh_ptmass,vxyz_ptmass,iboundary,isdeadh,ll,ideadhead,&
                kill_particle,shuffle_part,iamtype,iamdust,Bxyz,rad,radprop, &
                remove_particle_from_npartoftype,ien_type,ien_etotal,gr
 use eos,             only:gamma,polyk,eos_is_non_ideal
 use centreofmass,    only:get_centreofmass
 use options,         only:ieos,icooling,iexternalforce,use_dustfrac,use_hybrid
 use io,              only:id,master
 use externalforces,  only:accrete_particles,update_externalforce,accradius1,iext_star,iext_corotate
 use timestep,        only:time
 use units,           only:G_is_unity,get_G_code
 use boundary,        only:xmin,xmax,ymin,ymax,zmin,zmax
 use boundary_dyn,    only:dynamic_bdy,adjust_particles_dynamic_boundary
 use nicil,           only:n_nden
 use metric_tools,    only:imetric,imet_minkowski
 integer, intent(out) :: nerror,nwarn
 logical, intent(in), optional :: restart
 integer      :: i,nbad,itype,iu,ndead
 integer      :: ncount(maxtypes)
 real         :: xcom(ndim),vcom(ndim)
 real         :: hi,hmin,hmax
 logical      :: accreted,dorestart
 character(len=3) :: string
!
!--check that setup is sensible
!
 nerror = 0
 nwarn = 0
 if (present(restart)) then
    dorestart = restart
 else
    dorestart = .false.
 endif

 if (npart > maxp) then
    print*,'ERROR: npart (',npart,') > maxp (',maxp,')'
    nerror = nerror + 1
 endif
 if (any(npartoftype < 0)) then
    print*,'ERROR: npartoftype -ve: ',npartoftype(:)
    nerror = nerror + 1
 endif
 if (sum(npartoftype) > maxp) then
    print*,'ERROR: sum(npartoftype) > maxp ',sum(npartoftype(:))
    nerror = nerror + 1
 endif
 if (sum(npartoftype) /= npart) then
    print*,'ERROR: sum of npartoftype  /=  npart: np=',npart,' but sum=',sum(npartoftype)
    nerror = nerror + 1
 endif
 if (hfact < 1. .or. isnan(hfact)) then
    print*,'ERROR: hfact = ',hfact,', should be >= 1'
    nerror = nerror + 1
 endif
 if (polyk < 0. .or. isnan(polyk)) then
    print*,'ERROR: polyk = ',polyk,', should be >= 0'
    nerror = nerror + 1
 endif
 if (use_krome) then
    if (ieos /= 19) then
       print*, 'KROME setup. Only eos=19 makes sense.'
       nerror = nerror + 1
    endif
 else
    if (polyk < tiny(0.) .and. ieos /= 2) then
       print*,'WARNING! polyk = ',polyk,' in setup, speed of sound will be zero in equation of state'
       nwarn = nwarn + 1
    endif
    if (gamma <= 0.) then
       print*,'WARNING! gamma not set (should be set > 0 even if not used)'
       nwarn = nwarn + 1
    endif
 endif
 if (npart < 0) then
    print*,'ERROR: npart = ',npart,', should be >= 0'
    nerror = nerror + 1
 elseif (npart==0 .and. nptmass==0) then
    print*,'WARNING! setup: npart = 0 (and no sink particles either)'
    nwarn = nwarn + 1
 elseif (npart==0) then
    print*,'WARNING! setup contains no SPH particles (but has ',nptmass,' point masses)'
    nwarn = nwarn + 1
 endif

 if (maxphase==maxp) then
!
!--particle type should be specified in the dump file for any particles of unusual type
!
    if (all(iphase(1:npart)==0)) then
       if (any(npartoftype(2:) > 0)) then
          print*,'ERROR: npartoftype > 0 for non-gas particles, but types have not been assigned'
          nerror = nerror + 1
       endif
       iphase(1:npart) = isetphase(igas,iactive=.true.)
    elseif (any(iphase(1:npart)==0)) then
       print*,'ERROR: types need to be assigned to all particles (or none)'
       nerror = nerror + 1
    endif
!
!--Check that the numbers of each type add up correctly
!
    ncount(:) = 0
    nbad = 0
    ndead = 0
    do i=1,npart
       itype = iamtype(iphase(i))
       if (itype < 1 .or. itype > maxtypes) then
          nbad = nbad + 1
          if (nbad < 10) print*,'ERROR: unknown particle type ',itype,' on particle ',i
       else
          ncount(itype) = ncount(itype) + 1
       endif
       if (isdeadh(xyzh(4,i))) ndead = ndead + 1
    enddo
    if (nbad > 0) then
       print*,'ERROR: unknown value of particle type on ',nbad,' particles'
       nerror = nerror + 1
    endif
    if (any(ncount /= npartoftype)) then
       print*,'n(via iphase)=',ncount
       print*,'npartoftype  =',npartoftype
       print*,'ERROR: sum of types in iphase is not equal to npartoftype'
       nerror = nerror + 1
    endif
    if (ndead > 0) then
       print*,'ZOMBIE ALERT:',ndead,' DEAD PARTICLES in particle setup'
       if (ideadhead==0) then
          print*,'REBUILDING DEAD PARTICLE LIST...'
          do i=1,npart
             if (isdeadh(xyzh(4,i))) then
                ll(i) = ideadhead
                ideadhead = i
                call remove_particle_from_npartoftype(i,npartoftype)
             endif
          enddo
       endif
       nwarn = nwarn + 1
    endif
!
!--check that no empty particle types have been added
!
    do i=1,size(npartoftype)
       if (labeltype(i)=='empty' .and. npartoftype(i) > 0) then
          print "(/,1x,a,i0,a,i0,a,/)",'ERROR: ',npartoftype(i),' particles of type ',i,&
                ' set up but type=='//trim(labeltype(i))
          if (i==2) print "(a,/)",&
            ' *** This is the old dust particle type: edit setup to use itype=idust ***'
          nerror = nerror + 1
       endif
    enddo
 endif
!
!--should not have negative or zero smoothing lengths in initial setup
!
 nbad = 0
 hmax = 0.
 if (npart > 0) then
    hmin = xyzh(4,1)
 else
    hmin = 0.
 endif
 do i=1,npart
    hi = xyzh(4,i)
    if ((.not.dorestart .and. hi <= 0.) .or. hi > 1.e20) then
       nbad = nbad + 1
       if (nbad <= 10) print*,' particle ',i,' h = ',hi
    endif
    hmin = min(hi,hmin)
    hmax = max(hi,hmax)
 enddo
 if (nbad > 0) then
    print*,'ERROR: negative, zero or ridiculous h on ',nbad,' of ',npart,' particles'
    print*,' hmin = ',hmin,' hmax = ',hmax
    nerror = nerror + 1
 endif
!
!--check for NaNs in arrays
!
 call check_NaN(npart,xyzh,'position/smoothing length (xyzh array)',nerror)
 call check_NaN(npart,vxyzu,'velocity/thermal energy (vxyzu array)',nerror)
 if (mhd) call check_NaN(npart,Bxyz,'magnetic field (Bxyz array)',nerror)
!
!--check for negative thermal energies
!
 if (maxvxyzu==4) then
    nbad = 0
    iu = 4
    do i=1,npart
       if (.not.in_range(vxyzu(iu,i),0.) .and. xyzh(4,i) >= 0.) then !ignore accreted particles that have negative energies
          nbad = nbad + 1
          if (nbad <= 10) print*,' particle ',i,' u = ',vxyzu(iu,i)
       endif
    enddo
    if (nbad > 0) then
       print*,'ERROR: negative thermal energy on ',nbad,' of ',npart,' particles'
       nerror = nerror + 1
    endif
 else
    if (abs(gamma-1.) > tiny(gamma) .and. (ieos /= 2 .and. ieos /=9)) then
       print*,'*** ERROR: using isothermal EOS, but gamma = ',gamma
       gamma = 1.
       print*,'*** Resetting gamma to 1, gamma = ',gamma
       nwarn = nwarn + 1
    endif
 endif
!
!--check that mass of each type has been set
!
 do itype=1,maxtypes
    if (npartoftype(itype) > 0 .and. abs(massoftype(itype)) < tiny(0.)) then
       print*,'WARNING: npartoftype > 0 for '//trim(labeltype(itype))//' particles but massoftype = 0'
       nwarn = nwarn + 1
    endif
    if (npartoftype(itype) > 0 .and. .not.(in_range(massoftype(itype),0.))) then
       print*,'ERROR: massoftype = ',massoftype(itype),' for '//trim(labeltype(itype))// &
              ' particles (n'//trim(labeltype(itype))//' = ',npartoftype(itype),')'
       nerror = nerror + 1
    endif
 enddo
 !
 !  check for particles with identical positions
 !
 call check_for_identical_positions(npart,xyzh,nbad)
 if (nbad > 0) then
    print*,'ERROR: ',nbad,' of ',npart,' particles have identical or near-identical positions'
    nwarn = nwarn + 1
 endif
!
!  check for particles outside boundaries
!
 if (periodic) then
    nbad = 0
    if (dynamic_bdy) then
       call adjust_particles_dynamic_boundary(npart,xyzh)
    else
       do i=1,npart
          if (xyzh(1,i) < xmin .or. xyzh(1,i) > xmax &
          .or.xyzh(2,i) < ymin .or. xyzh(2,i) > ymax &
          .or.xyzh(3,i) < zmin .or. xyzh(3,i) > zmax) then
             nbad = nbad + 1
             if (nbad <= 10) print*,' particle ',i,' xyz = ',xyzh(1:3,i)
          endif
       enddo
       if (nbad > 0) then
          print*,'ERROR: ',nbad,' of ',npart,' particles setup OUTSIDE the periodic box'
          if (.not. dynamic_bdy) nerror = nerror + 1
       endif
    endif
 endif
!
!  warn about external force settings
!
 if (iexternalforce==iext_star .and. nptmass==0) then
    if (id==master) print "(a,/,a)",'WARNING: iexternalforce=1 does not conserve momentum:',&
                                    '         use a sink particle at r=0 if you care about this'
    nwarn = nwarn + 1
 endif
!
!--check for particles placed inside accretion boundaries
!
 if (iexternalforce > 0 .and. .not.dorestart .and. (.not.(gr .and. imetric==imet_minkowski))) then
    call update_externalforce(iexternalforce,time,0.)
    nbad = 0
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,massoftype,time,iexternalforce) &
    !$omp private(i,accreted) &
    !$omp reduction(+:nbad)
    do i=1,npart
       call accrete_particles(iexternalforce,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),massoftype(1),time,accreted)
       if (accreted) nbad = nbad + 1
    enddo
    !$omp end parallel do
    if (nbad > 0) then
       print*,'WARNING: ',nbad,' of ',npart,' particles setup within the accretion boundary'
       nwarn = nwarn + 1
    endif
    !--check if we are using a central accretor
    hi = 0.
    call accrete_particles(iexternalforce,0.,0.,0.,hi,massoftype(1),time,accreted)
    !--if so, check for unresolved accretion radius
    if (accreted .and. accradius1 < 0.5*hmin .and. accradius1 > 0.) then
       print*,'WARNING: accretion radius is unresolved by a factor of hmin/racc = ',hmin/accradius1
       print*,'(this will cause the code to run needlessly slow)'
       nwarn = nwarn + 1
    endif
 endif
!
!--check G=1 in code units where necessary
!
 if (gravity .or. nptmass > 0) then
    if (.not.G_is_unity()) then
       if (gravity) then
          if (id==master) print*,'ERROR: self-gravity ON but G /= 1 in code units, got G=',get_G_code()
       elseif (nptmass > 0) then
          if (id==master) print*,'ERROR: sink particles used but G /= 1 in code units, got G=',get_G_code()
       endif
       nerror = nerror + 1
    endif
 endif
 if (.not. gr .and. (gravity .or. mhd) .and. ien_type == ien_etotal) then
    if (id==master) print*,'Cannot use total energy with self gravity or mhd'
    nerror = nerror + 1
 endif
!
!--sanity checks on magnetic field
!
 if (mhd) then
    if (all(abs(Bxyz(:,1:npart)) < tiny(0.))) then
       if (id==master) print*,'WARNING: MHD is ON but magnetic field is zero everywhere'
       nwarn = nwarn + 1
    endif
    if (mhd_nonideal) then
       if (n_nden /= n_nden_phantom) then
          if (id==master) print*,'ERROR: n_nden in nicil.f90 needs to match n_nden_phantom in config.F90; n_nden = ',n_nden
          nerror = nerror + 1
       endif
    endif
 endif
!
!--check -DDUST is set if dust particles are used
!
 if (npartoftype(idust) > 0) then
    if (.not. use_dust) then
       if (id==master) print*,'ERROR: dust particles present but -DDUST is not set'
       nerror = nerror + 1
    endif
    if (use_dustfrac) then
       call get_environment_variable('PHANTOM_RESTART_ONEFLUID',string)
       if (index(string,'yes') > 0) then
          if (id==master) print "(/,a,/)",' DELETING DUST PARTICLES (from PHANTOM_RESTART_ONEFLUID=yes)'
          if (maxphase==maxp) then
             do i=1,npart
                if (iamdust(iphase(i))) call kill_particle(i)
             enddo
          endif
          call shuffle_part(npart)
          npartoftype(idust) = 0
       else
          if (id==master) then
             if (use_hybrid) then
                print*, "WARNING: HYBRID DUST IMPLEMENTATION IS NOT YET FINISHED (IT'S NOT GONNA RUN ANYWAY)"
             else
                print "(2(/,a),/)",' ** Set PHANTOM_RESTART_ONEFLUID=yes to restart a two fluid', &
                                '    calculation using the one fluid method (dustfrac) **'
             endif
          endif
       endif
    endif
 endif
!
!--check dust grid is sensible
!
 if (use_dustfrac) call check_setup_dustgrid(nerror,nwarn)
!
!--check dust fraction is 0->1 if one fluid dust is used
!
 if (use_dustfrac .and. npart > 0) call check_setup_dustfrac(nerror,nwarn)
!
!--check GR setup
!
 if (gr) call check_gr(npart,nerror,xyzh,vxyzu)
!
!--check radiation setup
!
 if (do_radiation) call check_setup_radiation(npart,nerror,nwarn,radprop,rad)
!
!--check dust growth arrays
!
 if (use_dustgrowth) call check_setup_growth(npart,nerror)
!
!--check dust nucleation arrays
!
 if (do_nucleation) call check_setup_nucleation(npart,nerror)
!
!--check point mass setup
!
 call check_setup_ptmass(nerror,nwarn,hmin)
!
!--check centre of mass
!
 call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (.not.h2chemistry .and. maxvxyzu >= 4 .and. icooling == 3 .and. iexternalforce/=iext_corotate .and. nptmass==0) then
    if (dot_product(xcom,xcom) >  1.e-2) then
       print*,'ERROR: Gammie (2001) cooling (icooling=3) assumes Omega = 1./r^1.5'
       print*,'                but the centre of mass is not at the origin!'
       nerror = nerror + 1
    endif
 endif

 if (nerror==0 .and. nwarn==0) then
!
!--print centre of mass (must be done AFTER types have been checked)
!  ALSO, only print this if there are no warnings or errors to avoid obscuring warnings
!
    if (id==master) then
       write(*,"(a,2(es10.3,', '),es10.3,a)") ' Centre of mass is at (x,y,z) = (',xcom,')'
       write(*,"(1x,a)") 'Particle setup OK'
    endif
 endif

end subroutine check_setup

!----------------------------------------------------
!+
!  function to check for NaNs in particle arrays
!+
!----------------------------------------------------
subroutine check_NaN(npart,array,label,nerror)
 integer,          intent(in)    :: npart
 real,             intent(in)    :: array(:,:)
 character(len=*), intent(in)    :: label
 integer,          intent(inout) :: nerror
 integer :: nbad,i

 nbad = 0
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(npart,array,label) &
 !$omp private(i) &
 !$omp reduction(+:nbad)
 do i=1,npart
    !--check for NaNs in xyzh
    if (any(isnan(array(:,i)))) then
       if (nbad < 10) print*,'NaN in '//trim(label)//' : ', i
       nbad = nbad + 1
    endif
 enddo
 !$omp end parallel do
 if (nbad > 0) then
    print*,'ERROR: NaN in '//trim(label)//' on ',nbad,' of ',npart,' particles'
    nerror = nerror + 1
 endif

end subroutine check_NaN

!----------------------------------------------------
!+
! function to check if a value is
! within the allowed range
! if min/max arguments are given
! Otherwise just checks for NaNs and Infs
!+
!----------------------------------------------------
pure logical function in_range(x,min,max)
 real, intent(in) :: x
 real, intent(in), optional :: min,max

 in_range = .true.
 if ((x /= x) .or. (x > huge(x))) then
    in_range = .false.
 endif
 if (present(min)) then
    if (x < min) then
       in_range = .false.
    endif
 endif
 if (present(max)) then
    if (x > max) then
       in_range = .false.
    endif
 endif

end function in_range

!------------------------------------------------------------------
!+
! check sink particle properties are sensible
!+
!------------------------------------------------------------------
subroutine check_setup_ptmass(nerror,nwarn,hmin)
 use dim,  only:maxptmass
 use part, only:nptmass,xyzmh_ptmass,ihacc,ihsoft,gr,iTeff,sinks_have_luminosity,ilum
 use ptmass_radiation, only:isink_radiation
 integer, intent(inout) :: nerror,nwarn
 real,    intent(in)    :: hmin
 integer :: i,j,n
 real :: dx(3)
 real :: r,hsink

 if (gr .and. nptmass > 0) then
    print*,' ERROR: nptmass = ',nptmass, ' should be = 0 for GR'
    nwarn = nwarn + 1
    return
 endif

 if (nptmass < 0) then
    print*,' ERROR: nptmass = ',nptmass, ' should be >= 0 '
    nerror = nerror + 1
 endif
 if (nptmass > maxptmass) then
    print*,' ERROR: nptmass = ',nptmass,' exceeds ptmass array dimensions of ',maxptmass
    nerror = nerror + 1
    return
 endif

 !
 !  check that sinks have not been placed on top of each other
 !  or within each others accretion radii
 !
 do i=1,nptmass
    if (xyzmh_ptmass(4,i) < 0.) cycle
    do j=i+1,nptmass
       if (xyzmh_ptmass(4,j) < 0.) cycle
       dx = xyzmh_ptmass(1:3,j) - xyzmh_ptmass(1:3,i)
       r  = sqrt(dot_product(dx,dx))
       if (r <= tiny(r)) then
          print*,'ERROR: sink ',j,' on top of sink ',i,' at ',xyzmh_ptmass(1:3,i)
          nerror = nerror + 1
       elseif (r <= max(xyzmh_ptmass(ihacc,i),xyzmh_ptmass(ihacc,j))) then
          print*,'WARNING: sinks ',i,' and ',j,' within each others accretion radii: sep =',&
                  r,' h = ',xyzmh_ptmass(ihacc,i),xyzmh_ptmass(ihacc,j)
          nwarn = nwarn + 1
       endif
    enddo
 enddo

 !
 !  check that sink masses are positive, warn if zero
 !
 n = 0
 do i=1,nptmass
    if (.not.in_range(xyzmh_ptmass(4,i))) then
       nerror = nerror + 1
       print*,' ERROR: sink ',i,' mass = ',xyzmh_ptmass(4,i)
    elseif (xyzmh_ptmass(4,i) < 0.) then
       print*,' Sink ',i,' has previously merged with another sink'
       n = n + 1
    endif
 enddo
 if (n > 0) print*,'The ptmass arrays have ',n,' sinks that have previously merged (i.e. have mass < 0)'
 !
 !  check that accretion radii are positive
 !
 do i=1,nptmass
    if (xyzmh_ptmass(4,i) < 0.) cycle
    hsink = max(xyzmh_ptmass(ihacc,i),xyzmh_ptmass(ihsoft,i))
    if (hsink <= 0.) then
       nerror = nerror + 1
       print*,'ERROR: sink ',i,' has accretion radius ',xyzmh_ptmass(ihacc,i),&
              ' and softening radius ',xyzmh_ptmass(ihsoft,i)
    elseif (hsink <= 0.5*hmin .and. hmin > 0.) then
       nwarn = nwarn + 1
       print*,'WARNING: sink ',i,' has unresolved accretion radius: hmin/racc = ',hmin/hsink
       print*,'         (this makes the code run pointlessly slow)'
    endif
 enddo
 !
 !  check that radiation properties are sensible
 !
 if (isink_radiation > 1 .and. xyzmh_ptmass(ilum,1) < 1e-10) then
    nerror = nerror + 1
    print*,'ERROR: isink_radiation > 1 and sink particle has no luminosity'
    return
 endif
 if (sinks_have_luminosity(nptmass,xyzmh_ptmass)) then
    if (any(xyzmh_ptmass(iTeff,1:nptmass) < 100.)) then
       print*,'WARNING: sink particle temperature less than 100K'
       nwarn = nwarn + 1
    endif
 endif

end subroutine check_setup_ptmass

!------------------------------------------------------------------
!+
! check dust growth arrays are sensible
!+
!------------------------------------------------------------------
subroutine check_setup_growth(npart,nerror)
 use part, only:dustprop,dustprop_label,iamdust,iphase,maxphase,maxp
 integer, intent(in)    :: npart
 integer, intent(inout) :: nerror
 integer :: i,j,nbad(4)

 nbad = 0
 !-- Check that all the parameters are > 0 when needed
 do i=1,npart
    do j=1,2
       if (maxphase==maxp) then
          if (iamdust(iphase(i)) .and. dustprop(j,i) <= 0.) nbad(j) = nbad(j) + 1
       elseif (dustprop(j,i) < 0.) then
          nbad(j) = nbad(j) + 1
       endif
    enddo
 enddo
 do j=1,2
    if (nbad(j) > 0) then
       print*,'ERROR: dustgrowth: ',nbad(j),' of ',npart,' particles with '//trim(dustprop_label(j))//' <= 0'
       nerror = nerror + 1
    endif
 enddo
 !-- check for NaN
 call check_NaN(npart,dustprop,'dust properties (dustprop array)',nerror)

end subroutine check_setup_growth

!------------------------------------------------------------------
!+
! check dust nucleation arrays are sensible
!+
!------------------------------------------------------------------
subroutine check_setup_nucleation(npart,nerror)
 use part, only:nucleation,nucleation_label,n_nucleation,idmu,idgamma
 integer, intent(in)    :: npart
 integer, intent(inout) :: nerror
 integer :: i,j,nbad(n_nucleation)

 nbad = 0
 !-- Check that all the parameters are > 0 when needed
 do i=1,npart
    if (nucleation(idmu,i) < 0.1) nbad(idmu) = nbad(idmu) + 1
    if (nucleation(idgamma,i) < 1.) nbad(idgamma) = nbad(idgamma) + 1
 enddo
 do j=1,n_nucleation
    if (nbad(j) > 0) then
       print*,'ERROR: ',nbad(j),' of ',npart,' particles with '//trim(nucleation_label(j))//' <= 0'
       nerror = nerror + 1
    endif
 enddo
 call check_NaN(npart,nucleation,'nucleation array',nerror)

end subroutine check_setup_nucleation

!------------------------------------------------------------------
!+
! check dust grain properties are sensible
!+
!------------------------------------------------------------------
subroutine check_setup_dustgrid(nerror,nwarn)
 use part,    only:grainsize,graindens,ndustsmall,ndustlarge,ndusttypes
 use options, only:use_dustfrac
 use dim,     only:use_dust,use_dustgrowth
 use units,   only:udist
 use physcon, only:km
 use dust,    only:idrag
 integer, intent(inout) :: nerror,nwarn
 integer :: i

 if (use_dustfrac .and. ndustsmall <= 0) then
    print*,'ERROR: Using one fluid dust but no dust types have been set (ndustsmall=',ndustsmall,')'
    nerror = nerror + 1
 endif
 if (use_dust .and. ndusttypes <= 0) then
    print*,'WARNING: Using dust but no dust species are used (ndusttypes=',ndusttypes,')'
    nwarn = nwarn + 1
 endif
 if (ndusttypes /= ndustsmall + ndustlarge) then
    print*,'ERROR: nsmall + nlarge ',ndustsmall+ndustlarge,&
           ' not equal to ndusttypes: ',ndusttypes
    nerror = nerror + 1
 endif
 !
 ! check that grain size array is sensible and non-zero
 ! except if dustgrowth is switched on, in which case the size defined in
 ! dustprop should be non-zero
 !
 if (use_dustgrowth) then
    ! dust growth not implemented for more than one grain size
    if (ndusttypes > 1) then
       print*,'ERROR: dust growth requires ndusttypes = 1, but ndusttypes = ',ndusttypes
       nerror = nerror + 1
    endif
 else
    do i=1,ndusttypes
       if (idrag == 1 .and. grainsize(i) <= 0.) then
          print*,'ERROR: grainsize = ',grainsize(i),' in dust bin ',i
          nerror = nerror + 1
       endif
    enddo
    do i=1,ndusttypes
       if (idrag == 1 .and. graindens(i) <= 0.) then
          print*,'ERROR: graindens = ',graindens(i),' in dust bin ',i
          nerror = nerror + 1
       endif
    enddo
    do i=1,ndusttypes
       if (grainsize(i) > 10.*km/udist) then
          print*,'WARNING: grainsize is HUGE (>10km) in dust bin ',i,': s = ',grainsize(i)*udist/km,' km'
          nwarn = nwarn + 1
       endif
    enddo
 endif

end subroutine check_setup_dustgrid

!------------------------------------------------------------------
!+
! check dust fractions are sensible
!+
!------------------------------------------------------------------
subroutine check_setup_dustfrac(nerror,nwarn)
 use io,   only:id,master
 use dim,  only:maxdusttypes
 use part, only:npart,ndustsmall,dustfrac
 integer, intent(inout) :: nerror,nwarn
 integer :: i,j,nbad,nunity
 real :: dust_to_gas_mean

 nbad = 0
 nunity = 0
 dust_to_gas_mean = 0.
 do i=1,npart
    do j=1,ndustsmall
       if (dustfrac(j,i) < 0. .or. dustfrac(j,i) > 1.) then
          nbad = nbad + 1
          if (nbad <= 10) print*,' particle ',i,' dustfrac = ',dustfrac(j,i)
       elseif (abs(dustfrac(j,i)-1.) < tiny(1.)) then
          nunity = nunity + 1
       else
          dust_to_gas_mean = dust_to_gas_mean + dustfrac(j,i)/(1. - sum(dustfrac(:,i)))
       endif
    enddo
 enddo
 dust_to_gas_mean = dust_to_gas_mean/real(npart-nbad-nunity)
 if (nbad > 0) then
    print*,'ERROR: ',nbad,' of ',npart,' particles with dustfrac outside [0,1]'
    nerror = nerror + 1
 endif
 if (nunity > 0) then
    print*,'WARNING: ',nunity,' of ',npart,' PARTICLES ARE PURE DUST (dustfrac=1.0)'
    nwarn = nwarn + 1
 endif

 ! check for TOTAL dustfrac exceeding unity
 nbad = 0
 do i=1,npart
    if (sum(dustfrac(1:ndustsmall,i)) > 1.) then
       nbad = nbad + 1
       if (nbad <= 10) print*,' particle ',i,' total dustfrac = ',sum(dustfrac(1:ndustsmall,i))
    endif
 enddo
 if (nbad > 0) then
    print*,'ERROR: ',nbad,' of ',npart,' particles with sum of dust fractions exceeding 1'
    nerror = nerror + 1
 endif

 ! warn if compiled for one-fluid dust but not used
 if (all(dustfrac(:,1:npart) < tiny(dustfrac))) then
    print*,'WARNING: one fluid dust is used but dust fraction is zero everywhere'
    if (maxdusttypes > 1) then
       print*,'WARNING about the previous WARNING: maxdusttypes > 1 so dust arrays are unnecessarily large!'
       print*,'                                    Recompile with maxdusttypes = 1 for better efficiency.'
    endif
    nwarn = nwarn + 1
 endif

 if (id==master) write(*,"(a,es10.3,/)") ' Mean dust-to-gas ratio is ',dust_to_gas_mean

end subroutine check_setup_dustfrac

!------------------------------------------------------------------
!+
! check things necessary for sensible GR evolution
!+
!------------------------------------------------------------------
subroutine check_gr(npart,nerror,xyzh,vxyzu)
 use metric_tools, only:pack_metric,unpack_metric
 use utils_gr,     only:get_u0
 use part,         only:isdead_or_accreted,ien_type,ien_entropy,ien_etotal,ien_entropy_s
 use units,        only:in_geometric_units,get_G_code,get_c_code
 integer, intent(in)    :: npart
 integer, intent(inout) :: nerror
 real,    intent(in)    :: xyzh(:,:),vxyzu(:,:)
 real    :: metrici(0:3,0:3,2),gcov(0:3,0:3),u0
 integer :: ierr,i,nbad

 !
 ! check code units are set for geometric units
 !
 if (.not. in_geometric_units()) then
    print "(/,a)",' ERROR: units are incorrect for GR, need G = c = 1'
    print *,' ...but we have G = ',get_G_code(),' and c = ',get_c_code()
    print*
    nerror = nerror + 1
 endif

 !
 ! check for bad U0, indicating v or u > 1
 !
 nbad = 0
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call pack_metric(xyzh(1:3,i),metrici)
       call unpack_metric(metrici,gcov=gcov)
       call get_u0(gcov,vxyzu(1:3,i),U0,ierr)
       if (ierr/=0) nbad = nbad + 1
    endif
 enddo

 if (nbad > 0) then
    print "(/,a,i10,a,i10,a,/)",' ERROR in setup: ',nbad,' of ',npart,' particles have |v| > 1 or u > 1, giving undefined U^0'
    nerror = nerror + 1
 endif

 if (ien_type /= ien_etotal .and. ien_type /= ien_entropy .and. ien_type /= ien_entropy_s) then
    print "(/,a,i1,a,i1,a,i3,/)",' ERROR: ien_type is incorrect for GR, need ', &
                                 ien_entropy, ', ', ien_etotal, ' or ', ien_entropy_s, &
                                 ' but get ', ien_type
    nerror = nerror + 1
 endif

end subroutine check_gr

!------------------------------------------------------------------
!+
! check for particles with identical positions
! the brute force approach of checking N^2 pairs is WAY too
! slow to be practical here. Instead, we sort the particles by
! radius and check particles with identical radii
!+
!------------------------------------------------------------------
subroutine check_for_identical_positions(npart,xyzh,nbad)
 use sortutils, only:indexxfunc,r2func
 use part,      only:maxphase,maxp,iphase,igas,iamtype,isdead_or_accreted
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(out) :: nbad
 integer :: i,j,itypei,itypej
 real    :: dx(3),dx2
 integer, allocatable :: index(:)
 !
 ! sort particles by radius
 !
 allocate(index(npart))
 call indexxfunc(npart,r2func,xyzh,index)
 !
 ! check for identical positions. Stop checking as soon as non-identical
 ! positions are found.
 !
 nbad = 0
 itypei = igas
 itypej = igas
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,index,maxphase,maxp,iphase) &
 !$omp firstprivate(itypei,itypej) &
 !$omp private(i,j,dx,dx2) &
 !$omp reduction(+:nbad)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,index(i)))) then
       j = i+1
       dx2 = 0.
       if (maxphase==maxp) itypei = iamtype(iphase(index(i)))
       do while (dx2 < epsilon(dx2) .and. j < npart)
          if (isdead_or_accreted(xyzh(4,index(j)))) exit
          dx = xyzh(1:3,index(i)) - xyzh(1:3,index(j))
          if (maxphase==maxp) itypej = iamtype(iphase(index(j)))
          dx2 = dot_product(dx,dx)
          if (dx2 < epsilon(dx2) .and. itypei==itypej) then
             nbad = nbad + 1
             if (nbad <= 10) then
                print*,'WARNING: particles of same type at same position: '
                print*,' ',index(i),':',xyzh(1:3,index(i))
                print*,' ',index(j),':',xyzh(1:3,index(j))
             endif
          endif
          j = j + 1
       enddo
    endif
 enddo
 !$omp end parallel do

 deallocate(index)

end subroutine check_for_identical_positions

!------------------------------------------------------------------
!+
!  1) check for optically thin particles when mcfost is disabled,
!     as the particles will then be overlooked if they are flagged as thin
!  2) check that radiation energy is never negative to begin with
!  3) check for NaNs
!+
!------------------------------------------------------------------
subroutine check_setup_radiation(npart,nerror,nwarn,radprop,rad)
 use part, only:ithick, iradxi, ikappa
 integer, intent(in)    :: npart
 integer, intent(inout) :: nerror,nwarn
 real,    intent(in)    :: rad(:,:), radprop(:,:)
 integer :: i,nthin,nradEn,nkappa,nwarn_en

 nthin = 0
 nradEn = 0
 nkappa = 0
 nwarn_en = 0
 do i=1,npart
    if (radprop(ithick, i) < 0.5) nthin = nthin + 1
    if (rad(iradxi, i) < 0.) nradEn = nradEn + 1
    if (radprop(ikappa, i) <= 0.0 .or. isnan(radprop(ikappa,i))) nkappa = nkappa + 1
    if (rad(iradxi, i) <= 0.) nwarn_en = nwarn_en + 1
 enddo

 if (nthin > 0) then
    print "(/,a,i10,a,i10,a,/)",' ERROR in setup: ',nthin,' of ',npart,&
    ' particles are being treated as optically thin without MCFOST being compiled'
    nerror = nerror + 1
 endif

 if (nradEn > 0) then
    print "(/,a,i10,a,i10,a,/)",' ERROR in setup: ',nradEn,' of ',npart,&
    ' particles have negative radiation energy'
    nerror = nerror + 1
 endif

 if (nwarn_en > 0) then
    print "(/,a,i10,a,i10,a,/)",' WARNING in setup: ',nwarn_en,' of ',npart,&
    ' particles have radiation energy equal to zero'
    nwarn = nwarn + 1
 endif

 if (nkappa > 0) then
    print "(/,a,i10,a,i10,a,/)",' ERROR in setup: ',nkappa,' of ',npart,&
    ' particles have opacity <= 0.0 or NaN'
    nerror = nerror + 1
 endif

 call check_NaN(npart,rad,'radiation_energy',nerror)
 call check_NaN(npart,radprop,'radiation properties',nerror)

end subroutine check_setup_radiation

end module checksetup
