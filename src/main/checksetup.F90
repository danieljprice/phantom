!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: checksetup
!
!  DESCRIPTION:
!   Perform sanity checks of the particle setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, centreofmass, dim, eos, externalforces, io,
!    options, part, physcon, timestep, units
!+
!--------------------------------------------------------------------------
module checksetup
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
 use dim,  only:maxp,maxvxyzu,periodic,use_dust,ndim,mhd,maxdusttypes,use_dustgrowth
 use part, only:xyzh,massoftype,hfact,vxyzu,npart,npartoftype,nptmass,gravity, &
                iphase,maxphase,isetphase,labeltype,igas,h2chemistry,maxtypes,&
                idust,xyzmh_ptmass,vxyz_ptmass,dustfrac,iboundary,&
                kill_particle,shuffle_part,iamdust,Bxyz,ndustsmall
 use eos,             only:gamma,polyk
 use centreofmass,    only:get_centreofmass
 use options,         only:ieos,icooling,iexternalforce,use_dustfrac
 use io,              only:id,master
 use externalforces,  only:accrete_particles,accradius1,iext_star,iext_corotate
 use timestep,        only:time
 use units,           only:umass,udist,utime
 use physcon,         only:gg
 use boundary,        only:xmin,xmax,ymin,ymax,zmin,zmax
 integer, intent(out) :: nerror,nwarn
 logical, intent(in), optional :: restart
 integer      :: i,j,nbad,itype,nunity,iu
 real         :: xcom(ndim),vcom(ndim)
 real(kind=8) :: gcode
 real         :: hi,hmin,hmax,dust_to_gas
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
    print*,'Error in setup: npart (',npart,') > maxp (',maxp,')'
    nerror = nerror + 1
 endif
 if (any(npartoftype < 0)) then
    print*,'Error in setup: npartoftype -ve: ',npartoftype(:)
    nerror = nerror + 1
 endif
 if (sum(npartoftype) > maxp) then
    print*,'Error in setup: sum(npartoftype) > maxp ',sum(npartoftype(:))
    nerror = nerror + 1
 endif
 if (gamma <= 0.) then
    print*,'WARNING! Error in setup: gamma not set (should be set > 0 even if not used)'
    nwarn = nwarn + 1
 endif
 if (hfact < 1. .or. hfact /= hfact) then
    print*,'Error in setup: hfact = ',hfact,', should be >= 1'
    nerror = nerror + 1
 endif
 if (polyk < 0. .or. polyk /= polyk) then
    print*,'Error in setup: polyk = ',polyk,', should be >= 0'
    nerror = nerror + 1
 elseif (polyk < tiny(0.) .and. ieos /= 2) then
    print*,'WARNING! polyk = ',polyk,' in setup, speed of sound will be zero in equation of state'
    nwarn = nwarn + 1
 endif
 if (npart < 0) then
    print*,'Error in setup: npart = ',npart,', should be >= 0'
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
          print*,'Error in setup: npartoftype > 0 for non-gas particles, but types have not been assigned'
          nerror = nerror + 1
       endif
       iphase(1:npart) = isetphase(igas,iactive=.true.)
    elseif (any(iphase(1:npart)==0)) then
       print*,'Error in setup: types need to be assigned to all particles (or none)'
       nerror = nerror + 1
    endif
!
!--If boundary particles are present, then only gas and boundary particles may exist
!
    if (npartoftype(iboundary) > 0) then
       do i = 1,maxtypes
          if (npartoftype(i) > 0 .and. (i/=igas .and. i/=iboundary)) then
             print*, 'Error in setup: boundary particles cannot coexist with non-gas particles'
             nerror = nerror + 1
          endif
       enddo
    endif
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
    !--check for NaNs in xyzh
    if (any(xyzh(:,i) /= xyzh(:,i))) then
       print*,'NaN in position/smoothing length (xyzh array) : ', i
       nerror = nerror + 1
    endif
    !--check for NaNs in velocity
    if (any(vxyzu(:,i) /= vxyzu(:,i))) then
       if (maxvxyzu >= 4) then
          print*,'NaN in velocity/utherm (vxyzu array) : ', i
       else
          print*,'NaN in velocity field (vxyzu array) : ', i
       endif
       nerror = nerror + 1
    endif
    !--check for NaNs in B field
    if (mhd) then
       if (any(Bxyz(:,i) /= Bxyz(:,i))) then
          print*,'NaN in magnetic field (Bxyz array) : ', i
          nerror = nerror + 1
       endif
    endif
    hi = xyzh(4,i)
    if ((.not.dorestart .and. hi <= 0.) .or. hi > 1.e20) then
       nbad = nbad + 1
       if (nbad <= 10) print*,' particle ',i,' h = ',hi
    endif
    hmin = min(hi,hmin)
    hmax = max(hi,hmax)
 enddo
 if (nbad > 0) then
    print*,'Error in setup: negative, zero or ridiculous h on ',nbad,' of ',npart,' particles'
    print*,' hmin = ',hmin,' hmax = ',hmax
    nerror = nerror + 1
 endif
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
       print*,'Error in setup: negative thermal energy on ',nbad,' of ',npart,' particles'
       nerror = nerror + 1
    endif
 else
    if (abs(gamma-1.) > tiny(gamma) .and. (ieos /= 2 .and. ieos /=9)) then
       print*,'*** Error in setup: using isothermal EOS, but gamma = ',gamma
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
       print*,'Error in setup: massoftype = ',massoftype(itype),' for '//trim(labeltype(itype))// &
              ' particles (n'//trim(labeltype(itype))//' = ',npartoftype(itype),')'
       nerror = nerror + 1
    endif
 enddo
!
!  check for particles outside boundaries
!
 if (periodic) then
    nbad = 0
    do i=1,npart
       if (xyzh(1,i) < xmin .or. xyzh(1,i) > xmax &
       .or.xyzh(2,i) < ymin .or. xyzh(2,i) > ymax &
       .or.xyzh(3,i) < zmin .or. xyzh(3,i) > zmax) then
          nbad = nbad + 1
          if (nbad <= 10) print*,' particle ',i,' xyz = ',xyzh(1:3,i)
       endif
    enddo
    if (nbad > 0) then
       print*,'Error in setup: ',nbad,' of ',npart,' particles setup OUTSIDE the periodic box'
       nerror = nerror + 1
    endif
 endif
!
!  warn about external force settings
!
 if (iexternalforce==iext_star .and. nptmass==0) then
    print*,'WARNING: iexternalforce=1 does not conserve momentum - use a sink particle at r=0 if you care about this'
    nwarn = nwarn + 1
 endif
!
!--check for particles placed inside accretion boundaries
!
 if (iexternalforce > 0 .and. .not.dorestart) then
    nbad = 0
    do i=1,npart
       call accrete_particles(iexternalforce,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),massoftype(1),time,accreted)
       if (accreted) nbad = nbad + 1
    enddo
    if (nbad > 0) then
       print*,'Warning: ',nbad,' of ',npart,' particles setup within the accretion boundary'
       nwarn = nwarn + 1
    endif
    !--check if we are using a central accretor
    hi = 0.
    call accrete_particles(iexternalforce,0.,0.,0.,hi,massoftype(1),time,accreted)
    !--if so, check for unresolved accretion radius
    if (accreted .and. accradius1 < 0.5*hmin) then
       print*,'Warning: accretion radius is unresolved by a factor of hmin/racc = ',hmin/accradius1
       print*,'(this will cause the code to run needlessly slow)'
       nwarn = nwarn + 1
    endif
 endif
!
!--check G=1 in code units where necessary
!
 if (gravity .or. nptmass > 0) then
    gcode = gg*umass*utime**2/udist**3
    if (abs(gcode-1.) > max(1.e-15,real(epsilon(gcode)))) then
       if (gravity) then
          print*,'Error in setup: self-gravity ON but G /= 1 in code units'
       elseif (nptmass > 0) then
          print*,'Error in setup: sink particles used but G /= 1 in code units'
          print*,gcode,gcode-1.,epsilon(gcode)
       endif
       nerror = nerror + 1
    endif
 endif
!
!--sanity checks on magnetic field
!
 if (mhd) then
    if (all(abs(Bxyz(:,1:npart)) < tiny(0.))) then
       print*,'WARNING: MHD is ON but magnetic field is zero everywhere'
       nwarn = nwarn + 1
    endif
 endif
!
!--check -DDUST is set if dust particles are used
!
 if (npartoftype(idust) > 0) then
    if (.not. use_dust) then
       if (id==master) print*,'Error in setup: dust particles present but -DDUST is not set'
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
             print*,'ERROR in setup: use of dust particles AND a dust fraction not implemented'
             print*,'                i.e. cannot yet mix two-fluid and one-fluid methods'
             print "(2(/,a),/)",' ** Set PHANTOM_RESTART_ONEFLUID=yes to restart a two fluid', &
                                '    calculation using the one fluid method (dustfrac) **'
          endif
          nerror = nerror + 1
       endif
    endif
 endif
!
!--check dust fraction is 0->1 if one fluid dust is used
!
 if (use_dustfrac) then
    nbad = 0
    nunity = 0
    dust_to_gas = 0.
    do i=1,npart
       do j=1,ndustsmall
          if (dustfrac(j,i) < 0. .or. dustfrac(j,i) > 1.) then
             nbad = nbad + 1
             if (nbad <= 10) print*,' particle ',i,' dustfrac = ',dustfrac(j,i)
          elseif (abs(dustfrac(j,i)-1.) < tiny(1.)) then
             nunity = nunity + 1
          else
             dust_to_gas = dust_to_gas + dustfrac(j,i)/(1. - sum(dustfrac(:,i)))
          endif
       enddo
    enddo
    if (nbad > 0) then
       print*,'ERROR: ',nbad,' of ',npart,' particles with dustfrac outside [0,1]'
       nerror = nerror + 1
    endif
    if (nunity > 0) then
       print*,'WARNING: ',nunity,' of ',npart,' PARTICLES ARE PURE DUST (dustfrac=1.0)'
       nwarn = nwarn + 1
    endif
    ! warn if compiled for one-fluid dust but not used
    if (all(dustfrac(:,1:npart) < tiny(dustfrac))) then
       print*,'WARNING: one fluid dust is used but dust fraction is zero everywhere'
       if (maxdusttypes>1) then
          print*,'WARNING about the previous WARNING: maxdusttypes > 1 so dust arrays are unnecessarily large!'
          print*,'                                    Recompile with maxdusttypes = 1 for better efficiency.'
       endif
       nwarn = nwarn + 1
    endif
    if (id==master) write(*,"(a,es10.3,/)") ' Mean dust-to-gas ratio is ',dust_to_gas/real(npart-nbad-nunity)
 endif

!
!--check dust growth arrays
!
 if (use_dustgrowth) call check_setup_growth(npart,nerror)
!
!--check point mass setup
!
 call check_setup_ptmass(nerror,nwarn,hmin)
!
!--print centre of mass (must be done AFTER types have been checked)
!
 call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 if (id==master) &
    write(*,"(a,2(es10.3,', '),es10.3,a)") ' Centre of mass is at (x,y,z) = (',xcom,')'

 if (.not.h2chemistry .and. maxvxyzu >= 4 .and. icooling >= 1 .and. iexternalforce/=iext_corotate) then
    if (dot_product(xcom,xcom) >  1.e-2) then
       print*,'Error in setup: Gammie (2001) cooling (icooling=1) assumes Omega = 1./r^1.5'
       print*,'                but the centre of mass is not at the origin!'
       nerror = nerror + 1
    endif
 endif

 if (nerror==0 .and. nwarn==0) then
    if (id==master) write(*,"(1x,a)") 'Particle setup OK'
 endif

 return
end subroutine check_setup

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

subroutine check_setup_ptmass(nerror,nwarn,hmin)
 use dim,  only:maxptmass
 use part, only:nptmass,xyzmh_ptmass,ihacc,ihsoft,gr
 integer, intent(inout) :: nerror,nwarn
 real,    intent(in)    :: hmin
 integer :: i,j,n
 real :: dx(3)
 real :: r,hsink

 if (gr .and. nptmass > 0) then
    print*,' Warning! Error in setup: nptmass = ',nptmass, ' should be = 0 for GR'
    nwarn = nwarn + 1
    return
 endif

 if (nptmass < 0) then
    print*,' Error in setup: nptmass = ',nptmass, ' should be >= 0 '
    nerror = nerror + 1
 endif
 if (nptmass > maxptmass) then
    print*,' Error in setup: nptmass = ',nptmass,' exceeds ptmass array dimensions of ',maxptmass
    nerror = nerror + 1
    return
 endif

 !
 !  check that sinks have not been placed on top of each other
 !  or within each others accretion radii
 !
 do i=1,nptmass
    do j=i+1,nptmass
       dx = xyzmh_ptmass(1:3,j) - xyzmh_ptmass(1:3,i)
       r  = sqrt(dot_product(dx,dx))
       if (r <= tiny(r)) then
          print*,'Error in setup: sink ',j,' on top of sink ',i,' at ',xyzmh_ptmass(1:3,i)
          nerror = nerror + 1
       elseif (r <= max(xyzmh_ptmass(ihacc,i),xyzmh_ptmass(ihacc,j))) then
          print*,'Warning: sinks ',i,' and ',j,' within each others accretion radii: sep =',&
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
    if (.not.in_range(xyzmh_ptmass(4,i),0.)) then
       nerror = nerror + 1
       print*,' Error in setup: sink ',i,' mass = ',xyzmh_ptmass(4,i)
    elseif (xyzmh_ptmass(4,i) < tiny(0.)) then
       n = n + 1
    endif
 enddo
 if (n > 0) then
    print*,'WARNING: ',n,' sink particles have zero mass '
    nwarn = nwarn + 1
 endif
 !
 !  check that accretion radii are positive
 !
 do i=1,nptmass
    hsink = max(xyzmh_ptmass(ihacc,i),xyzmh_ptmass(ihsoft,i))
    if (hsink <= 0.) then
       nerror = nerror + 1
       print*,'Error in setup: sink ',i,' has accretion radius ',xyzmh_ptmass(ihacc,i),&
              ' and softening radius ',xyzmh_ptmass(ihsoft,i)
    elseif (hsink <= 0.5*hmin .and. hmin > 0.) then
       nwarn = nwarn + 1
       print*,'Warning: sink ',i,' has unresolved accretion radius: hmin/racc = ',hmin/hsink
       print*,'         (this makes the code run pointlessly slow)'
    endif
 enddo

end subroutine check_setup_ptmass

subroutine check_setup_growth(npart,nerror)
 use part, only:dustprop,dustprop_label
 integer, intent(in)    :: npart
 integer, intent(inout) :: nerror
 integer :: i,j,nbad(4)

 nbad = 0
 !-- Check that all the parameters are > 0 when needed
 do i=1,npart
    do j=1,4
       if (dustprop(j,i) < 0.) nbad(j) = nbad(j) + 1
    enddo
 enddo

 do j=1,4
    if (nbad(j) > 0) then
       print*,'ERROR: ',nbad,' of ',npart,' with '//trim(dustprop_label(j))//' < 0'
       nerror = nerror + 1
    endif
 enddo

end subroutine check_setup_growth

end module checksetup
