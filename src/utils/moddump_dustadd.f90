!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! adds dust particles to pre-existing gas particle dump
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - H_R_ref         : *H/R at R_ref (size distribution)*
!   - R_ref           : *reference radius (size distribution)*
!   - dust_method     : *dust method (1=one fluid, 2=two fluid)*
!   - dust_to_gas     : *total dust to gas ratio*
!   - graindens_cgs   : *grain density [g/cm^3]*
!   - grainsize_cgs   : *grain size [cm]*
!   - icutinside      : *delete particles inside a given radius*
!   - icutoutside     : *delete particles outside a given radius*
!   - incenterx       : *x coordinate of the centre of the inner sphere*
!   - incentery       : *y coordinate of the centre of the inner sphere*
!   - incenterz       : *z coordinate of the centre of the inner sphere*
!   - inradius        : *inward radius [au]*
!   - iporosity       : *use porosity (0=no, 1=yes)*
!   - iremoveparttype : *which particles to delete (0=all, 1=gas only, 2=dust only)*
!   - ngrainsizes     : *number of grain sizes*
!   - np_ratio        : *ratio between number of gas and dust particles (two fluid)*
!   - outcenterx      : *x coordinate of the centre of the outer sphere*
!   - outcentery      : *y coordinate of the centre of the outer sphere*
!   - outcenterz      : *z coordinate of the centre of the outer sphere*
!   - outradius       : *outward radius [au]*
!   - pwl_sizedistrib : *power-law index of the size distribution*
!   - q_index         : *q index (size distribution)*
!   - sindex          : *power-law index of the grain-size mass distribution (e.g. MRN)*
!   - sizedistrib     : *set dust size via a size distribution*
!   - smaxcgs         : *maximum grain size [cm]*
!   - smincgs         : *minimum grain size [cm]*
!
! :Dependencies: dim, dust, growth, infile_utils, io, moddump_utils,
!   options, part, prompting, set_dust, units
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer :: dust_method = 1        ! 1=one fluid, 2=two fluid
 integer :: np_ratio    = 5        ! ratio between number of gas and dust particles
 real    :: dust_to_gas = 0.01     ! total dust to gas ratio
 integer :: ngrainsizes = 1        ! number of grain sizes (-> ndusttypes)
 real    :: smincgs = 1.e-5        ! minimum grain size in cm
 real    :: smaxcgs = 1.           ! maximum grain size in cm
 real    :: sindex  = 3.5          ! power-law index of the grain-size mass distribution
 real    :: grainsize_cgs = 1.     ! grain size in cm (single size / at R_ref)
 real    :: graindens_cgs = 3.     ! grain density in g/cm^3
 integer :: iporosity = 0          ! use porosity (0=no, 1=yes)
 logical :: sizedistrib = .false.  ! set dust size via a size distribution
 real    :: pwl_sizedistrib = -2.  ! power-law index of the size distribution
 real    :: R_ref   = 100.         ! reference radius
 real    :: H_R_ref = 0.0895       ! H/R at R_ref
 real    :: q_index = 0.25         ! q index

 logical :: icutinside  = .false.
 logical :: icutoutside = .false.
 real    :: inradius  = 10.
 real    :: outradius = 200.
 real    :: incenter(3)  = 0.
 real    :: outcenter(3) = 0.
 integer :: iremoveparttype = 0    ! 0=all, 1=gas only, 2=dust only

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,          only:use_dust,maxdusttypes,maxdustlarge,maxdustsmall,use_dustgrowth, &
                        update_max_sizes
 use part,         only:igas,idust,set_particle_type,ndusttypes,ndustsmall,ndustlarge, &
                        grainsize,graindens,dustfrac,delete_particles_outside_sphere
 use set_dust,     only:set_dustfrac,set_dustbinfrac
 use options,      only:use_dustfrac,use_porosity
 use growth,       only:set_dustprop,convert_to_twofluid,iporosity_growth=>iporosity
 use dust,         only:grainsizecgs,graindenscgs
 use units,        only:umass,udist
 use io,           only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,j,itype,ipart,iloc,np_gas,np_dust,maxdust,iremovetype
 real    :: dustbinfrac(maxdusttypes),udens

 if (.not. use_dust) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes'
    stop
 endif
 udens = umass/(udist**3)
 dustbinfrac = 0.

 !- grainsize and graindens already set if convert from one fluid to two fluid with growth
 if (.not. (use_dustfrac .and. use_dustgrowth)) then
    grainsize = 1.
    graindens = 3.
 endif
 grainsizecgs = 1.
 graindenscgs = 3.

 !
 !--prompt only when the driver did not find a parameter file
 !
 if (prompt_for_params) call read_interactive_moddumpfile()

 ndusttypes = ngrainsizes
 iporosity_growth = iporosity
 if (iporosity == 1) use_porosity = .true.

 if (use_dustgrowth .and. use_dustfrac) then
    print*,' Detected dustgrowth AND dustfrac: converting from one to two fluid'
    call convert_to_twofluid(npart,xyzh,vxyzu,massoftype,npartoftype,np_ratio,dust_to_gas)
 else
    if (dust_method==1) then
       maxdust = maxdustsmall
    else
       maxdust = maxdustlarge
    endif
    if (ndusttypes > maxdust) call fatal('moddump_dustadd','ngrainsizes exceeds the maximum for this dust method')

    if (ndusttypes > 1) then
       !--grainsizes + mass distribution
       call set_dustbinfrac(smincgs/udist,smaxcgs/udist,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
       !--grain density
       graindens(:) = graindens_cgs/udens
    else
       grainsizecgs = grainsize_cgs
       graindenscgs = graindens_cgs
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

 !
 !--delete particles if necessary
 !
 iremovetype = 0
 if (icutinside .or. icutoutside) then
    select case (iremoveparttype)
    case (1)
       iremovetype = igas
    case (2)
       iremovetype = idust
    case default
       iremovetype = 0
    end select
 endif

 if (icutinside) then
    print*,'Phantommoddump: Remove particles inside a particular radius'
    print*,'Removing particles inside radius ',inradius
    if (iremovetype > 0) then
       print*,'Removing particles type ',iremovetype
       call delete_particles_outside_sphere(incenter,inradius,npart,revert=.true.,mytype=iremovetype)
    else
       call delete_particles_outside_sphere(incenter,inradius,npart,revert=.true.)
    endif
 endif

 if (icutoutside) then
    print*,'Phantommoddump: Remove particles outside a particular radius'
    print*,'Removing particles outside radius ',outradius
    if (iremovetype > 0) then
       print*,'Removing particles type ',iremovetype
       call delete_particles_outside_sphere(outcenter,outradius,npart,mytype=iremovetype)
    else
       call delete_particles_outside_sphere(outcenter,outradius,npart)
    endif
 endif

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactively set the moddump parameters
!+
!----------------------------------------------------------------
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt
 use dim,       only:use_dustgrowth,maxdustsmall,maxdustlarge
 use options,   only:use_dustfrac
 integer :: maxdust

 if (use_dustgrowth .and. use_dustfrac) then
    call prompt('Enter ratio between number of gas particles and dust particles',np_ratio,1)
    call prompt('Enter total dust to gas ratio',dust_to_gas,0.)
 else
    call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
    if (dust_method==1) then
       maxdust = maxdustsmall
    else
       maxdust = maxdustlarge
       call prompt('Enter ratio between number of gas particles and dust particles',np_ratio,1)
    endif

    call prompt('Enter total dust to gas ratio',dust_to_gas,0.)
    call prompt('How many grain sizes do you want?',ngrainsizes,1,maxdust)

    if (ngrainsizes > 1) then
       call prompt('Enter minimum grain size in cm',smincgs,0.)
       call prompt('Enter maximum grain size in cm',smaxcgs,0.)
       call prompt('Enter power-law index, e.g. MRN',sindex)
       call prompt('Enter grain density in g/cm^3',graindens_cgs,0.)
    else
       if (use_dustgrowth) then
          call prompt('Use porosity ? (0=no,1=yes)',iporosity,0,1)
          call prompt('Set dust size via size distribution ?',sizedistrib)
          if (sizedistrib) then
             call prompt('Enter grain size in cm at Rref',grainsize_cgs,0.)
             call prompt('Enter power-law index ',pwl_sizedistrib)
             call prompt('Enter R_ref ',R_ref,0.)
             call prompt('Enter H/R at R_ref',H_R_ref,0.)
             call prompt('Enter q index',q_index)
          else
             call prompt('Enter initial grain size in cm',grainsize_cgs,0.)
          endif
       else
          call prompt('Enter grain size in cm',grainsize_cgs,0.)
       endif
       call prompt('Enter grain density in g/cm^3',graindens_cgs,0.)
    endif
 endif

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
 endif

end subroutine read_interactive_moddumpfile

!----------------------------------------------------------------
!+
!  write options to .moddump file
!+
!----------------------------------------------------------------
subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 use dim,          only:use_dustgrowth
 use options,      only:use_dustfrac
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)

 if (use_dustgrowth .and. use_dustfrac) then
    write(iunit,"(/,a)") '# converting from one fluid to two fluid (dustgrowth + dustfrac)'
    call write_inopt(np_ratio,'np_ratio','ratio between number of gas and dust particles (two fluid)',iunit)
    call write_inopt(dust_to_gas,'dust_to_gas','total dust to gas ratio',iunit)
 else
    write(iunit,"(/,a)") '# dust setup'
    call write_inopt(dust_method,'dust_method','dust method (1=one fluid, 2=two fluid)',iunit)
    if (dust_method == 2) call write_inopt(np_ratio,'np_ratio','ratio between number of gas and dust particles (two fluid)',iunit)
    call write_inopt(dust_to_gas,'dust_to_gas','total dust to gas ratio',iunit)
    call write_inopt(ngrainsizes,'ngrainsizes','number of grain sizes',iunit)
    if (ngrainsizes > 1) then
       call write_inopt(smincgs,'smincgs','minimum grain size [cm]',iunit)
       call write_inopt(smaxcgs,'smaxcgs','maximum grain size [cm]',iunit)
       call write_inopt(sindex,'sindex','power-law index of the grain-size mass distribution (e.g. MRN)',iunit)
       call write_inopt(graindens_cgs,'graindens_cgs','grain density [g/cm^3]',iunit)
    else
       if (use_dustgrowth) then
          call write_inopt(iporosity,'iporosity','use porosity (0=no, 1=yes)',iunit)
          call write_inopt(sizedistrib,'sizedistrib','set dust size via a size distribution',iunit)
          call write_inopt(grainsize_cgs,'grainsize_cgs','grain size [cm] (at R_ref if using a size distribution)',iunit)
          if (sizedistrib) then
             call write_inopt(pwl_sizedistrib,'pwl_sizedistrib','power-law index of the size distribution',iunit)
             call write_inopt(R_ref,'R_ref','reference radius (size distribution)',iunit)
             call write_inopt(H_R_ref,'H_R_ref','H/R at R_ref (size distribution)',iunit)
             call write_inopt(q_index,'q_index','q index (size distribution)',iunit)
          endif
       else
          call write_inopt(grainsize_cgs,'grainsize_cgs','grain size [cm]',iunit)
       endif
       call write_inopt(graindens_cgs,'graindens_cgs','grain density [g/cm^3]',iunit)
    endif
 endif

 write(iunit,"(/,a)") '# particle removal'
 call write_inopt(icutinside,'icutinside','delete particles inside a given radius',iunit)
 call write_inopt(icutoutside,'icutoutside','delete particles outside a given radius',iunit)
 call write_inopt(inradius,'inradius','inward radius [au]',iunit)
 call write_inopt(incenter(1),'incenterx','x coordinate of the centre of the inner sphere',iunit)
 call write_inopt(incenter(2),'incentery','y coordinate of the centre of the inner sphere',iunit)
 call write_inopt(incenter(3),'incenterz','z coordinate of the centre of the inner sphere',iunit)
 call write_inopt(outradius,'outradius','outward radius [au]',iunit)
 call write_inopt(outcenter(1),'outcenterx','x coordinate of the centre of the outer sphere',iunit)
 call write_inopt(outcenter(2),'outcentery','y coordinate of the centre of the outer sphere',iunit)
 call write_inopt(outcenter(3),'outcenterz','z coordinate of the centre of the outer sphere',iunit)
 call write_inopt(iremoveparttype,'iremoveparttype','which particles to delete (0=all, 1=gas only, 2=dust only)',iunit)

 close(iunit)

end subroutine write_moddump

!----------------------------------------------------------------
!+
!  read options from .moddump file
!+
!----------------------------------------------------------------
subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use dim,          only:use_dustgrowth
 use options,      only:use_dustfrac
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return

 if (use_dustgrowth .and. use_dustfrac) then
    call read_inopt(np_ratio,'np_ratio',db,min=1,errcount=nerr)
    call read_inopt(dust_to_gas,'dust_to_gas',db,min=0.,errcount=nerr)
 else
    call read_inopt(dust_method,'dust_method',db,min=1,max=2,errcount=nerr)
    if (dust_method == 2) call read_inopt(np_ratio,'np_ratio',db,min=1,errcount=nerr)
    call read_inopt(dust_to_gas,'dust_to_gas',db,min=0.,errcount=nerr)
    call read_inopt(ngrainsizes,'ngrainsizes',db,min=1,errcount=nerr)
    if (ngrainsizes > 1) then
       call read_inopt(smincgs,'smincgs',db,min=0.,errcount=nerr)
       call read_inopt(smaxcgs,'smaxcgs',db,min=0.,errcount=nerr)
       call read_inopt(sindex,'sindex',db,errcount=nerr)
       call read_inopt(graindens_cgs,'graindens_cgs',db,min=0.,errcount=nerr)
    else
       if (use_dustgrowth) then
          call read_inopt(iporosity,'iporosity',db,min=0,max=1,errcount=nerr)
          call read_inopt(sizedistrib,'sizedistrib',db,errcount=nerr)
          call read_inopt(grainsize_cgs,'grainsize_cgs',db,min=0.,errcount=nerr)
          if (sizedistrib) then
             call read_inopt(pwl_sizedistrib,'pwl_sizedistrib',db,errcount=nerr)
             call read_inopt(R_ref,'R_ref',db,min=0.,errcount=nerr)
             call read_inopt(H_R_ref,'H_R_ref',db,errcount=nerr)
             call read_inopt(q_index,'q_index',db,errcount=nerr)
          endif
       else
          call read_inopt(grainsize_cgs,'grainsize_cgs',db,min=0.,errcount=nerr)
       endif
       call read_inopt(graindens_cgs,'graindens_cgs',db,min=0.,errcount=nerr)
    endif
 endif

 call read_inopt(icutinside,'icutinside',db,errcount=nerr)
 call read_inopt(icutoutside,'icutoutside',db,errcount=nerr)
 call read_inopt(inradius,'inradius',db,min=0.,errcount=nerr)
 call read_inopt(incenter(1),'incenterx',db,errcount=nerr)
 call read_inopt(incenter(2),'incentery',db,errcount=nerr)
 call read_inopt(incenter(3),'incenterz',db,errcount=nerr)
 call read_inopt(outradius,'outradius',db,min=0.,errcount=nerr)
 call read_inopt(outcenter(1),'outcenterx',db,errcount=nerr)
 call read_inopt(outcenter(2),'outcentery',db,errcount=nerr)
 call read_inopt(outcenter(3),'outcenterz',db,errcount=nerr)
 call read_inopt(iremoveparttype,'iremoveparttype',db,min=0,max=2,errcount=nerr)

 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

end module moddump
