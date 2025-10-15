!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! adds dust particles to pre-existing gas particle dump
!
! :References: None
!
! :Owner: Stephane Michoulier
!
! :Runtime parameters: None
!
! :Dependencies: dim, dust, growth, options, part, porosity, prompting,
!   set_dust, table_utils, units
!

 use part,         only:delete_particles_outside_sphere,igas,idust
 use prompting,    only:prompt

 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,          only:use_dust,maxdusttypes,maxdustlarge,maxdustsmall,use_dustgrowth,&
                        update_max_sizes
 use part,         only:igas,idust,set_particle_type,ndusttypes,ndustsmall,ndustlarge,&
                        grainsize,graindens,dustfrac
 use set_dust,     only:set_dustfrac,set_dustbinfrac
 use options,      only:use_dustfrac,use_porosity
 use growth,       only:set_dustprop,convert_to_twofluid
 use porosity,     only:iporosity
 use prompting,    only:prompt
 use dust,         only:grainsizecgs,graindenscgs
 use table_utils,  only:logspace
 use units,   only:umass,udist
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, dimension(3) :: incenter,outcenter
 integer :: i,j,itype,ipart,iloc,dust_method,np_ratio,np_gas,np_dust,maxdust
 real    :: dust_to_gas,smincgs,smaxcgs,sindex,dustbinfrac(maxdusttypes),udens
 integer :: iremoveparttype
 real    :: inradius,outradius,pwl_sizedistrib,R_ref,H_R_ref,q_index
 logical :: icutinside,icutoutside,sizedistrib


 if (.not. use_dust) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes'
    stop
 endif
 udens = umass/(udist**3)
 dust_method = 1
 np_ratio = 5
 dust_to_gas = 0.01
 ndusttypes = 1
 smincgs = 1.e-5
 smaxcgs = 1.
 sindex = 3.5
 dustbinfrac = 0.
 pwl_sizedistrib = -2
 R_ref = 100
 H_R_ref = 0.0895
 q_index = 0.25

 icutinside      = .false.
 icutoutside     = .false.
 iremoveparttype = 0
 incenter(:)     = 0.
 outcenter(:)    = 0.
 inradius        = 10.
 outradius       = 200.

 !- grainsize and graindens already set if convert from one fluid to two fluid with growth
 if (.not. (use_dustfrac .and. use_dustgrowth)) then
    grainsize = 1.
    graindens = 3.
 endif
 grainsizecgs = 1.
 graindenscgs = 3.

 if (use_dustgrowth .and. use_dustfrac) then
    print*,' Detected dustgrowth AND dustfrac: converting from one to two fluid'

    call prompt('Enter ratio between number of gas particles and dust particles',np_ratio,1)
    call prompt('Enter total dust to gas ratio',dust_to_gas,0.)

    call convert_to_twofluid(npart,xyzh,vxyzu,massoftype,npartoftype,np_ratio,dust_to_gas)
 else
    call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)

    if (dust_method==1) then
       maxdust = maxdustsmall
    elseif (dust_method==2) then
       maxdust = maxdustlarge
       call prompt('Enter ratio between number of gas particles and dust particles',np_ratio,1)
       !--We do not care if modulo(npart,np_ratio) is stricly zero, since npart can
       !  be a weird value depdending on the simulation it comes from.
    endif

    call prompt('Enter total dust to gas ratio',dust_to_gas,0.)
    call prompt('How many grain sizes do you want?',ndusttypes,1,maxdust)

    if (ndusttypes > 1) then
       !--grainsizes
       call prompt('Enter minimum grain size in cm',smincgs,0.)
       call prompt('Enter maximum grain size in cm',smaxcgs,0.)
       !--mass distribution
       call prompt('Enter power-law index, e.g. MRN',sindex)
       call set_dustbinfrac(smincgs/udist,smaxcgs/udist,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
       !--grain density
       call prompt('Enter grain density in g/cm^3',graindens(1),0.)
       graindens = graindens(1)/udens
    else
       if (use_dustgrowth) then
          call prompt('Use porosity ? (0=no,1=yes)',iporosity,0,1)
          if (iporosity == 1) then
             use_porosity = .true.
          endif
          call prompt('Set dust size via size distribution ?',sizedistrib)
          if (sizedistrib) then
             call prompt('Enter grain size in cm at Rref',grainsizecgs,0.)
             call prompt('Enter power-law index ',pwl_sizedistrib)
             call prompt('Enter R_ref ',R_ref,0.)
             call prompt('Enter H/R at R_ref',H_R_ref,0.)
             call prompt('Enter q index',q_index)
          else
             call prompt('Enter initial grain size in cm',grainsizecgs,0.)
          endif
       else
          call prompt('Enter grain size in cm',grainsizecgs,0.)
       endif
       call prompt('Enter grain density in g/cm^3',graindenscgs,0.)
       grainsize(1) = grainsizecgs/udist
       graindens(1) = graindenscgs/udens
    endif

    np_gas = npartoftype(igas)

    if (dust_method == 1) then

       use_dustfrac = .true.
       ndustsmall = ndusttypes

       do i=1,np_gas
          if (ndusttypes > 1) then
             dustfrac(1:ndusttypes,i) = dust_to_gas*dustbinfrac(1:ndusttypes)
          else
             call set_dustfrac(dust_to_gas,dustfrac(:,i))
          endif
       enddo

       massoftype(igas) = massoftype(igas)*(1. + dust_to_gas)
       npart = np_gas

    elseif (dust_method == 2) then

       use_dustfrac = .false.
       ndustlarge = ndusttypes
       np_dust = np_gas/np_ratio
       npart = np_gas + np_dust*ndustlarge

       call update_max_sizes(npart)

       do i=1,ndustlarge

          itype = idust + i - 1
          npartoftype(itype) = np_dust
          massoftype(itype)  = massoftype(igas)*dust_to_gas*np_ratio

          do j=1,np_dust
             ipart = np_gas + (i-1)*np_dust + j
             iloc  = np_ratio*j

             xyzh(1,ipart) = xyzh(1,iloc)
             xyzh(2,ipart) = xyzh(2,iloc)
             xyzh(3,ipart) = xyzh(3,iloc)
             xyzh(4,ipart) = xyzh(4,iloc)

             vxyzu(1,ipart) = vxyzu(1,iloc)
             vxyzu(2,ipart) = vxyzu(2,iloc)
             vxyzu(3,ipart) = vxyzu(3,iloc)

             call set_particle_type(ipart,itype)
          enddo

       enddo
    endif
    if (use_dustgrowth) then
       call set_dustprop(npart,xyzh,sizedistrib,pwl_sizedistrib,R_ref,H_R_ref,q_index)
    endif
 endif
 !Delete particles if necessary

 !
 !--set the centers and the radius
 !
 call prompt('Deleting particles inside a given radius ?',icutinside)
 call prompt('Deleting particles outside a given radius ?',icutoutside)
 if (icutinside) then
    call prompt('Enter inward radius in au',inradius,0.)
    call prompt('Enter x coordinate of the center of that sphere',incenter(1))
    call prompt('Enter y coordinate of the center of that sphere',incenter(2))
    call prompt('Enter z coordinate of the center of that sphere',incenter(3))
 endif
 if (icutoutside) then
    call prompt('Enter outward radius in au',outradius,0.)
    call prompt('Enter x coordinate of the center of that sphere',outcenter(1))
    call prompt('Enter y coordinate of the center of that sphere',outcenter(2))
    call prompt('Enter z coordinate of the center of that sphere',outcenter(3))
 endif

 if (icutinside .or. icutoutside) then
    call prompt('Deleting which particles (0=all, 1=gas only, 2=dust only)?', iremoveparttype)
    ! add other types of particles here if needed
    select case (iremoveparttype)
    case (1)
       iremoveparttype = igas
    case (2)
       iremoveparttype = idust
    case default
       iremoveparttype = 0
    end select
 endif

 if (icutinside) then
    print*,'Phantommoddump: Remove particles inside a particular radius'
    print*,'Removing particles inside radius ',inradius
    if (iremoveparttype > 0) then
       print*,'Removing particles type ',iremoveparttype
       call delete_particles_outside_sphere(incenter,inradius,npart,revert=.true.,mytype=iremoveparttype)
    else
       call delete_particles_outside_sphere(incenter,inradius,npart,revert=.true.)
    endif
 endif

 if (icutoutside) then
    print*,'Phantommoddump: Remove particles outside a particular radius'
    print*,'Removing particles outside radius ',outradius
    if (iremoveparttype > 0) then
       print*,'Removing particles type ',iremoveparttype
       call delete_particles_outside_sphere(outcenter,outradius,npart,mytype=iremoveparttype)
    else
       call delete_particles_outside_sphere(outcenter,outradius,npart)
    endif
 endif

end subroutine modify_dump

end module moddump

