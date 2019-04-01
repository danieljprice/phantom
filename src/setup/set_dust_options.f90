!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: set_dust_options
!
!  DESCRIPTION:
!  Contains interactive set up for dust
!
!  REFERENCES:
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    dust_method   -- dust method (1=one fluid,2=two fluid)
!    dust_to_gas   -- dust to gas ratio
!    graindensinp  -- intrinsic grain density (in g/cm^3)
!    grainsizeinp  -- grain size (in cm)
!    igraindens    -- grain density input (0=equal,1=manually)
!    ndusttypesinp -- number of grain sizes
!    sindex        -- grain size power-law index (e.g. MRN = 3.5)
!    smaxcgs       -- max grain size (in cm)
!    smincgs       -- min grain size (in cm)
!
!  DEPENDENCIES: dim, dust, eos, fileutils, growth, infile_utils, io,
!    options, part, prompting
!+
!--------------------------------------------------------------------------
module set_dust_options
 use dim,       only:maxdusttypes,maxdustsmall,maxdustlarge,use_dustgrowth
 use prompting, only:prompt
 implicit none

 integer, public :: dust_method
 real,    public :: dust_to_gas
 integer, public :: ndusttypesinp
 real,    public :: grainsizeinp(maxdusttypes)
 real,    public :: graindensinp(maxdusttypes)
 integer, public :: igrainsize
 integer, public :: igraindens
 integer, public :: isetdust
 real,    public :: smincgs
 real,    public :: smaxcgs
 real,    public :: sindex
 real,    public :: dustbinfrac(maxdusttypes)
 real,    public :: Kdrag
 logical, public :: ilimitdustfluxinp

 public :: set_dust_default_options
 public :: set_dust_interactively
 public :: read_dust_setup_options
 public :: write_dust_setup_options
 public :: check_dust_method

 private

contains

!--------------------------------------------------------------------------
!+
!  Subroutine for setting default dust setup properties
!+
!--------------------------------------------------------------------------
subroutine set_dust_default_options()

 dust_method = 2
 dust_to_gas = 0.01
 ndusttypesinp = 1
 grainsizeinp(:) = 1.
 graindensinp(:) = 3.
 igrainsize = 0
 igraindens = 0
 isetdust = 0
 smincgs = 1.e-4
 smaxcgs = 1.
 sindex = 3.5
 dustbinfrac(:) = 0.
 dustbinfrac(1) = 1.
 Kdrag = 1000.
 ilimitdustfluxinp = .false.

end subroutine set_dust_default_options

!--------------------------------------------------------------------------
!+
!  Subroutine for setting dust properties interactively
!+
!--------------------------------------------------------------------------
subroutine set_dust_interactively()

 !--can only use one dust method
 call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
 call prompt('Enter total dust to gas ratio',dust_to_gas,0.)

 if (dust_method==1) then
    call prompt('How many grain sizes do you want?',ndusttypesinp,1,maxdustsmall)
    call prompt('Do you want to limit the dust flux?',ilimitdustfluxinp)
 elseif (dust_method==2) then
    call prompt('How many grain sizes do you want?',ndusttypesinp,1,maxdustlarge)
 endif

 if (ndusttypesinp > 1 .and. .not.use_dustgrowth) then
    call prompt('How do you want to set the grain sizes?'//new_line('A')// &
               ' 0=log-spaced'//new_line('A')// &
               ' 1=manually'//new_line('A'),igrainsize,0,1)
    call prompt('How do you want to set the (intrinsic) grain density?'//new_line('A')// &
               ' 0=equal'//new_line('A')// &
               ' 1=manually'//new_line('A'),igraindens,0,1)
 endif

 call prompt('How do you want to set the dust density profile?'//new_line('A')// &
            ' 0=equal to the gas'//new_line('A')// &
            ' 1=custom'//new_line('A')// &
            ' 2=equal to the gas, but with unique cutoffs'//new_line('A'),isetdust,0,2)

end subroutine set_dust_interactively

!--------------------------------------------------------------------------
!+
!  Subroutine for reading dust properties from a setup file
!+
!--------------------------------------------------------------------------
subroutine read_dust_setup_options(db,nerr)
 use growth,        only:read_growth_setup_options
 use infile_utils,  only:inopts,read_inopt
 use io,            only:error
 use fileutils,     only:make_tags_unique

 type(inopts), allocatable, intent(inout) :: db(:)
 integer, intent(inout) :: nerr

 character(len=120) :: varlabel(maxdusttypes)
 integer :: i,ierr

 !--can only use one method currently
 call read_inopt(dust_method,'dust_method',db,min=1,max=2,errcount=nerr)
 call read_inopt(dust_to_gas,'dust_to_gas',db,min=0.,errcount=nerr)
 if (dust_method == 1) then
    call read_inopt(ilimitdustfluxinp,'ilimitdustfluxinp',db,err=ierr,errcount=nerr)
 endif

 !--options for setting up the dust grid
 if (dust_method == 1) then
    call read_inopt(ndusttypesinp,'ndusttypesinp',db,min=1,max=maxdustsmall,errcount=nerr)
 elseif (dust_method == 2) then
    call read_inopt(ndusttypesinp,'ndusttypesinp',db,min=1,max=maxdustlarge,errcount=nerr)
 endif

 if (ndusttypesinp > 1) then
    call read_inopt(igrainsize,'igrainsize',db,min=0,max=1,errcount=nerr)
    select case(igrainsize)
    case(0)
       call read_inopt(smincgs,'smincgs',db,min=0.,errcount=nerr)
       call read_inopt(smaxcgs,'smaxcgs',db,min=smincgs,errcount=nerr)
       call read_inopt(sindex ,'sindex' ,db,errcount=nerr)
    case(1)
       varlabel = 'grainsizeinp'
       call make_tags_unique(ndusttypesinp,varlabel)
       do i=1,ndusttypesinp
          call read_inopt(grainsizeinp(i),trim(varlabel(i)),db,min=0.,err=ierr,errcount=nerr)
       enddo
       varlabel = 'dustbinfrac'
       call make_tags_unique(ndusttypesinp,varlabel)
       do i=1,ndusttypesinp
          call read_inopt(dustbinfrac(i),trim(varlabel(i)),db,min=0.,max=1.,err=ierr,errcount=nerr)
       enddo
       if (abs(sum(dustbinfrac(:)) - 1.) > epsilon(1.)) then
          call error('set_dust','dust bin fraction needs to add up to 1!')
          nerr = nerr+1
       endif
    end select
    call read_inopt(igraindens,'igraindens',db,min=0,errcount=nerr)
    select case(igraindens)
    case(0)
       call read_inopt(graindensinp(1),'graindensinp',db,min=0.,err=ierr,errcount=nerr)
    case(1)
       varlabel = 'graindensinp'
       call make_tags_unique(ndusttypesinp,varlabel)
       do i=1,ndusttypesinp
          call read_inopt(graindensinp(i),trim(varlabel(i)),db,min=0.,err=ierr,errcount=nerr)
       enddo
    end select
 else
    call read_inopt(grainsizeinp(1),'grainsizeinp',db,min=0.,err=ierr,errcount=nerr)
    call read_inopt(graindensinp(1),'graindensinp',db,min=0.,errcount=nerr)
 endif

 if (use_dustgrowth) call read_growth_setup_options(db,nerr)

end subroutine read_dust_setup_options

!--------------------------------------------------------------------------
!+
!  Subroutine for writing dust properties to a setup file
!+
!--------------------------------------------------------------------------
subroutine write_dust_setup_options(iunit)
 use growth,        only:write_growth_setup_options
 use infile_utils,  only:write_inopt
 use fileutils,     only:make_tags_unique

 integer, intent(in) :: iunit

 character(len=120) :: varlabel(maxdusttypes)
 character(len=120) :: varstring(maxdusttypes)
 integer :: i

 write(iunit,"(/,a)") '# options for dust'

 !--can only use one method per calculation currently
 call write_inopt(dust_method,'dust_method','dust method (1=one fluid,2=two fluid)',iunit)
 call write_inopt(dust_to_gas,'dust_to_gas','dust to gas ratio',iunit)

 call write_inopt(ndusttypesinp,'ndusttypesinp','number of grain sizes',iunit)

 if (dust_method==1) then
    call write_inopt(ilimitdustfluxinp,'ilimitdustfluxinp',&
       'limit dust diffusion using Ballabio et al. (2018)',iunit)
 endif

 if (ndusttypesinp > 1) then

    call write_inopt(igrainsize,'igrainsize', &
       'grain size distribution (0=log-space,1=manually)',iunit)
    select case(igrainsize)
    case(0)
       call write_inopt(smincgs,'smincgs','min grain size (in cm)',iunit)
       call write_inopt(smaxcgs,'smaxcgs','max grain size (in cm)',iunit)
       call write_inopt(sindex ,'sindex' ,'grain size power-law index (e.g. MRN = 3.5)',iunit)
    case(1)
       varlabel = 'grainsizeinp'
       varstring = 'grain size'
       call make_tags_unique(ndusttypesinp,varlabel)
       call make_tags_unique(ndusttypesinp,varstring)
       do i=1,ndusttypesinp
          varstring(i) = trim(varstring(i))//' (in cm)'
          call write_inopt(grainsizeinp(i),trim(varlabel(i)),trim(varstring(i)),iunit)
       enddo
       varlabel = 'dustbinfrac'
       varstring = 'dust bin fraction'
       call make_tags_unique(ndusttypesinp,varlabel)
       call make_tags_unique(ndusttypesinp,varstring)
       do i=1,ndusttypesinp
          varstring(i) = trim(varstring(i))//' (frac. of total dust)'
          call write_inopt(dustbinfrac(i),trim(varlabel(i)),trim(varstring(i)),iunit)
       enddo
    end select
    call write_inopt(igraindens,'igraindens','grain density input (0=equal,1=manually)',iunit)
    select case(igraindens)
    case(0)
       call write_inopt(graindensinp(1),'graindensinp','intrinsic grain density (in g/cm^3)',iunit)
    case(1)
       varlabel = 'graindensinp'
       varstring = 'grain density'
       call make_tags_unique(ndusttypesinp,varlabel)
       call make_tags_unique(ndusttypesinp,varstring)
       do i=1,ndusttypesinp
          varstring(i) = trim(varstring(i))//' (in g/cm^3)'
          call write_inopt(graindensinp(i),trim(varlabel(i)),trim(varstring(i)),iunit)
       enddo
    end select

 else

    if (use_dustgrowth) then
       call write_inopt(grainsizeinp(1),'grainsizeinp','initial grain size (in cm)',iunit)
       call write_inopt(graindensinp(1),'graindensinp','intrinsic grain density (in g/cm^3)',iunit) ! Modify this is graindens becomes variable
    else
       call write_inopt(grainsizeinp(1),'grainsizeinp','grain size (in cm)',iunit)
       call write_inopt(graindensinp(1),'graindensinp','intrinsic grain density (in g/cm^3)',iunit)
    endif

 endif

 call write_inopt(isetdust,'isetdust', &
    'how to set dust density profile (0=equal to gas,1=custom,2=equal to gas with cutoffs)',iunit)

 if (use_dustgrowth) call write_growth_setup_options(iunit)

end subroutine write_dust_setup_options

!--------------------------------------------------------------------------
!+
!  Subroutine for deciding wheather to use one-fluid or two-fluid dust
!+
!--------------------------------------------------------------------------
subroutine check_dust_method(dust_method,ichange_method)
 use dust,    only:init_drag,get_ts,idrag
 use eos,     only:ieos,get_spsound
 use io,      only:master
 use options, only:use_dustfrac
 use part,    only:npart,massoftype,xyzh,vxyzu,rhoh,igas,dustfrac,&
                   grainsize,graindens,ndusttypes
 integer,          intent(inout) :: dust_method
 logical,          intent(out)   :: ichange_method
 integer :: i,l,iregime,ierr,icheckdust
 real    :: r,rhogasi,rhodusti,rhoi,dustfracisum,spsoundi
 real    :: dustfraci(maxdusttypes),tsi(maxdusttypes)
 character(len=120) :: string
 logical :: iforce_dust_method

 iforce_dust_method = .false.

 call init_drag(ierr)

 dustfraci(:) = 0.
 icheckdust = 0
 do i=1,npart
    r = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
    if (use_dustfrac) then
       rhoi = rhoh(xyzh(4,i),massoftype(igas))
       dustfraci(1:ndusttypes) = dustfrac(1:ndusttypes,i)
       dustfracisum = sum(dustfraci(1:ndusttypes))
       rhogasi      = rhoi*(1.-dustfracisum)
       spsoundi     = get_spsound(ieos,xyzh(:,i),rhogasi,vxyzu(:,i))
       do l=1,ndusttypesinp
          rhodusti = rhoi*dustfraci(l)
          call get_ts(idrag,grainsize(l),graindens(l),rhogasi,rhodusti,spsoundi,0.,tsi(l),iregime)
       enddo
       if (any(tsi(1:ndusttypes) > xyzh(4,i)/spsoundi)) icheckdust = icheckdust + 1
    endif
 enddo

 call get_environment_variable('IFORCE_DUST_METHOD',string)
 if (trim(string)=='yes') iforce_dust_method = .true.

 ichange_method = .false.
 if (real(icheckdust)/real(npart) > 0.1 .and. .not.iforce_dust_method) then
    if (dust_method == 1) then
      ! use_dustfrac = .false.
       ichange_method = .true.
       dust_method = 2
    endif

    print*,''
    print*,'-------------------------------------------------------------------------------'
    print*,''
    print*,'    WARNING! More than 10% of particles have Stokes number greater than'
    print*,'    the threshold under which the terminal velocity approximation is valid.'
    print*,'    We suggest you switch to the "two-fluid" method. You can set the'
    print*,'    environment variable IFORCE_DUST_METHOD=yes to not see this message'
    print*,'    again.'
    print*,''
    print*,"Particles not satisfying the condition:",real(icheckdust)/real(npart)*100,"%"
    print*,'-------------------------------------------------------------------------------'
    print*,''
 endif

end subroutine check_dust_method

end module set_dust_options
