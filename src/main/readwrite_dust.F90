!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: readwrite_dust
!
!  DESCRIPTION:
!  Contains most I/O routines for dust
!
!  REFERENCES:
!
!  OWNER: Mark Hutchison
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    K_code            -- drag constant when constant drag is used
!    dust_method       -- dust method (1=one fluid,2=two fluid)
!    dust_to_gas_ratio -- dust to gas ratio
!    graindens         -- Intrinsic grain density in g/cm^3
!    graindensinp      -- intrinsic grain density (in g/cm^3)
!    grainsize         -- Grain size in cm
!    grainsizeinp      -- grain size (in cm)
!    icut_backreaction -- cut the drag on the gas phase (0=no, 1=yes)
!    idrag             -- gas/dust drag (0=off,1=Epstein/Stokes,2=const K,3=const ts)
!    ilimitdustflux    -- limit the dust flux using Ballabio et al. (2018)
!    io_graindens      -- grain density input (0=equal,1=manually)
!    sindex            -- grain size power-law index (e.g. MRN = 3.5)
!    smaxcgs           -- max grain size (in cm)
!    smincgs           -- min grain size (in cm)
!
!  DEPENDENCIES: dim, dump_utils, dust, eos, growth, infile_utils, io,
!    options, part, prompting, units
!+
!--------------------------------------------------------------------------

module readwrite_dust
 use dim,  only:ndusttypes
 use dust, only:smincgs,smaxcgs,sindex,grainsizecgs,graindenscgs,K_code,idrag, &
                ilimitdustflux
 implicit none
 !--Default values for the dust in the infile
 integer, public :: io_grainsize = 0
 integer, public :: io_graindens = 0
 integer, public :: dust_method = 1
 logical, public :: ichange_method = .false.

 public :: get_onefluiddust
 public :: set_dustfrac_from_inopts
 public :: nduststrings
 ! generic interface to interactively_set_dust
 interface interactively_set_dust
  module procedure interactively_set_dust_simple, interactively_set_dust_full
 end interface interactively_set_dust
 public :: interactively_set_dust
 public :: read_dust_setup_options
 public :: read_dust_infile_options
 public :: extract_dust_header
 public :: write_dust_to_header
 public :: write_dust_setup_options
 public :: write_dust_infile_options
 public :: write_temp_grains_file
 public :: check_dust_method

 private

contains

!--------------------------------------------------------------------
!+
!  extract whether one-fluid dust used in Phantom from the fileid
!+
!--------------------------------------------------------------------
subroutine get_onefluiddust(dumpfile,use_onefluiddust,fileid)
 use io,         only:idisk1
 use dump_utils, only:open_dumpfile_r,lenid
 character(len=lenid), optional, intent(inout)  :: fileid
 character(len=*), intent(in)  :: dumpfile
 logical,          intent(out) :: use_onefluiddust
 integer :: ierr
 character(len=lenid) :: fileident

 if (.not.present(fileid)) then
    call open_dumpfile_r(idisk1,dumpfile,fileident,ierr)
    close(idisk1)
 else
    fileident = fileid
 endif

 if (index(fileident,'+1dust') /= 0) then
    use_onefluiddust = .true.
 else
    use_onefluiddust = .false.
 endif

end subroutine get_onefluiddust


!----------------------------------------------------------------
!+
!  utility wrapper function to initialise the dust fraction
!  after reading the temparary dust file
!+
!----------------------------------------------------------------
subroutine set_dustfrac_from_inopts(dust_to_gas,percent,ipart)
 use part, only:dustfrac
 use dust, only:set_dustfrac
 integer, intent(in), optional :: ipart
 real,    intent(in), optional :: percent(:)
 real,    intent(in) :: dust_to_gas
 integer :: i
 real    :: dustfrac_temp(ndusttypes)
 real    :: dustfrac_multiplier(ndusttypes) = 1.

 select case(io_grainsize)
 case(0)
    call set_dustfrac(dust_to_gas,dustfrac_temp,smincgs,smaxcgs,sindex)
 case(1)
    call set_dustfrac(dust_to_gas,dustfrac_temp)
    dustfrac_multiplier = 1./real(ndusttypes)
 case(2)
    call set_dustfrac(dust_to_gas,dustfrac_temp)
    if (present(percent)) then
       dustfrac_multiplier = 1.e-2*percent
    else
       dustfrac_multiplier = 1./real(ndusttypes)
       print*,'Warning! percent of dust fraction not passed to set_dustfrac_from_inopts!'
       print*,'   ...guessing dustfrac_multiplier = ',dustfrac_multiplier
    endif
 end select

 if (present(ipart)) then
    dustfrac(:,ipart) = dustfrac_multiplier(:)*dustfrac_temp(:)
 else
    do i = 1,ndusttypes
       dustfrac(i,:) = dustfrac_multiplier(i)*dustfrac_temp(i)
    enddo
 endif

end subroutine set_dustfrac_from_inopts


!-----------------------------------------------------------------------------
!+
!  utility function to number strings for N dust species
!+
!-----------------------------------------------------------------------------
subroutine nduststrings(pre_string,post_string,complete_string)
 use io,  only:fatal
 character(len=*), intent(in) :: pre_string,post_string
 character(len=120), intent(out) :: complete_string(ndusttypes)

 integer :: i,total_len,int_len
 character(len=20) :: num_string

 int_len = 0
 if (ndusttypes > 1) int_len = floor(log10(real(ndusttypes) + tiny(0.))) + 1
 total_len = len(pre_string) + int_len  + len(post_string)
 if (len(complete_string) < total_len) then
    call fatal('dust','N dust string is not long enough!')
 endif

 if (ndusttypes > 1) then
    do i = 1,ndusttypes
       write(num_string,'(I0)') i
       write(complete_string(i),'(A)') pre_string//trim(adjustl(num_string))//post_string
    enddo
 else
    write(complete_string,'(A)') pre_string//post_string
 endif

 return
end subroutine nduststrings


!-----------------------------------------------------------------------
!+
!  Subroutine for setting simple dust properties interactively
!+
!-----------------------------------------------------------------------
subroutine interactively_set_dust_simple(dust_to_gas,imethod,Kdrag,units)
 use dim,       only:use_dustgrowth
 use options,   only:use_dustfrac
 use prompting, only:prompt
 use io,        only:fatal
 use dust,      only:set_grainsize,idrag,K_code,grainsizecgs,graindenscgs
 real,    intent(out) :: dust_to_gas
 logical, intent(in),  optional :: Kdrag
 integer, intent(out), optional :: imethod
 character(len=*), intent(in), optional :: units

 real    :: units_to_cgs = 1.
 real    :: grainsizeinp(ndusttypes),graindensinp(ndusttypes)
 character(len=120) :: message
 character(len=10)  :: abbrev = 'cm'

 if (present(imethod)) dust_method = imethod
 if (present(Kdrag)) then
    if (Kdrag) idrag = 2
 endif
 if (present(units)) then
    call get_units_factor(units,units_to_cgs)
    abbrev = units
 endif

 !--initialise the size/density to the default values
 grainsizeinp = grainsizecgs
 graindensinp = graindenscgs

 if (use_dustfrac) ilimitdustflux = .false.

 if (dust_to_gas <= 0.) dust_to_gas = 1. ! for a more sensible better option
 call prompt('Enter dust to gas ratio',dust_to_gas,0.)

 io_grainsize = 1
 io_graindens = 0

 select case(idrag)
 case(1)
    if (use_dustgrowth) then
       message = 'Enter initial grain size in'
    else
       message = 'Enter grain size in'
    endif
    message = trim(message)//' '//trim(abbrev)
    grainsizeinp(1) = grainsizeinp(1)/units_to_cgs
    call prompt(trim(message),grainsizeinp(1),0.)
    grainsizecgs(:) = units_to_cgs*grainsizeinp(1)

    call prompt('Enter the intrinsic grain density in g/cm^3',graindensinp(1),0.)
    graindenscgs(:) = graindensinp(1)
 case(2,3)
    if (use_dustfrac) K_code = 1000. ! for a more sensible better option
    call prompt('Enter constant drag coefficient',K_code,0.)
 case default
    stop 'ERROR! Invalid idrag option received in interactively_set_dust_simple'
 end select

end subroutine interactively_set_dust_simple


!-----------------------------------------------------------------------
!+
!  Subroutine for setting dust properties interactively
!  (originally written for setup_disc, but also needed for moddump)
!+
!-----------------------------------------------------------------------
subroutine interactively_set_dust_full(dust_to_gas,dustfrac_percent,grainsizeinp,graindensinp, &
                                  imethod,iprofile,Kdrag,units)
 use dim,       only:use_dustgrowth
 use options,   only:use_dustfrac
 use prompting, only:prompt
 use io,        only:fatal
 use dust,      only:set_grainsize,idrag,grainsizecgs,graindenscgs,K_code
 real,    intent(out) :: dust_to_gas
 real,    intent(out) :: dustfrac_percent(:)
 real,    intent(out) :: grainsizeinp(:),graindensinp(:)
 logical, intent(in),  optional :: Kdrag
 integer, intent(out), optional :: imethod,iprofile
 character(len=*), intent(in), optional :: units

 integer :: i
 integer :: dust_method = -1
 integer :: profile_set_dust = -1
 real    :: units_to_cgs = 1.
 logical :: simple_grainsize = .false.
 logical :: simple_graindens = .false.
 character(len=120) :: varstring(ndusttypes),varstringalt(ndusttypes)
 character(len=120) :: message
 character(len=10)  :: abbrev = 'cm'

 if (present(imethod)) dust_method = imethod
 if (present(iprofile)) profile_set_dust = iprofile
 if (present(Kdrag)) then
    idrag = 2
 else
    idrag = 1
 endif
 if (present(units)) then
    call get_units_factor(units,units_to_cgs)
    abbrev = units
 endif

 !--initialise the size/density to the default values
 grainsizeinp = grainsizecgs
 graindensinp = graindenscgs
 
 if (dust_method /= -1) then
    !
    !--dust method
    !
    if (ndusttypes > 1) then
       dust_method  = 1
    else
       dust_method  = 2
    endif

    call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
    if (dust_method == 1) then
       use_dustfrac = .true.
    else
       use_dustfrac = .false.
       if (ndusttypes > 1) call fatal('setup','dust_method=2 is currently only compatible with ndusttypes=1!')
    endif
    if (present(imethod)) imethod = dust_method
 endif
 
 if (use_dustfrac) call prompt('Do you want to limit the dust flux?',ilimitdustflux)
 select case(idrag)
 case(1)
    if (profile_set_dust /= -1) then
       call prompt('How do you want to set the dust density profile?'//new_line('A')// &
                  ' 0=equal to the gas'//new_line('A')// &
                  ' 1=custom'//new_line('A')// &
                  ' 2=equal to the gas, but with unique cutoffs'//new_line('A'),profile_set_dust,0,2)
       if (present(iprofile)) iprofile = profile_set_dust
    endif
    if (ndusttypes > 1 .and. .not.use_dustgrowth) then
       call prompt('Enter total dust to gas ratio',dust_to_gas,0.)
       call prompt('How do you want to set the grain sizes?'//new_line('A')// &
                  ' 0=power-law'//new_line('A')// &
                  ' 1=equal'//new_line('A')// &
                  ' 2=manually'//new_line('A'),io_grainsize,0,2)
       select case(io_grainsize)
       case(0)

          message = 'Enter minimum grain size in '//trim(abbrev)
          smincgs = smincgs/units_to_cgs
          call prompt(trim(message),smincgs,0.)
          smincgs = units_to_cgs*smincgs

          message = 'Enter maximum grain size in '//trim(abbrev)
          smaxcgs = smaxcgs/units_to_cgs
          call prompt(trim(message),smaxcgs,smincgs)
          smaxcgs = units_to_cgs*smaxcgs

          call prompt('Enter power-law index, e.g. MRN',sindex)
          call set_grainsize(smincgs,smaxcgs)
          grainsizeinp(:) = grainsizecgs(:)
       case(1)
          message = 'Enter grain size for all grains in '//trim(abbrev)
          grainsizeinp(1) = grainsizeinp(1)/units_to_cgs
          call prompt(trim(message),grainsizeinp(1),0.)
          grainsizeinp(:) = units_to_cgs*grainsizeinp(1)
       case(2)
          message = ' (in '//trim(abbrev)//')'
          call nduststrings('Enter grain size ',trim(message),varstring)
          call nduststrings('Enter dust fraction ',' (in % of total dust)',varstringalt)
          dustfrac_percent(:) = 100./ndusttypes
          simple_grainsize = .true.
          do i = 1,ndusttypes
             grainsizeinp(i) = grainsizeinp(i)/units_to_cgs
             call prompt(trim(varstring(i)),grainsizeinp(i),0.)
             grainsizeinp(i) = units_to_cgs*grainsizeinp(i)
             if (i == 1 .or. ndusttypes == 1) then
                call prompt(trim(varstringalt(i)),dustfrac_percent(i),0.,100.)
             elseif (i == ndusttypes .and. ndusttypes > 1) then
                dustfrac_percent(i) = 100. - sum(dustfrac_percent(1:i-1))
                print*,trim(varstringalt(i))//'...based on previous choices :', &
                      dustfrac_percent(i)
             else
                call prompt(trim(varstringalt(i)),dustfrac_percent(i),0., &
                           100. - sum(dustfrac_percent(1:i-1)))
             endif
          enddo
          if (any(grainsizeinp /= grainsizeinp(1)) .or. &
              any(dustfrac_percent /= dustfrac_percent(1))) simple_grainsize = .false.
       end select
       call prompt('How do you want to set the (intrinsic) grain density?'//new_line('A')// &
                  ' 0=equal'//new_line('A')// &
                  ' 1=manually'//new_line('A'),io_graindens,0,1)
       select case(io_graindens)
       case(0)
          call prompt('Enter the intrinsic grain density for all grain sizes in g/cm^3',graindensinp(1),0.)
          graindensinp(:) = graindensinp(1)
       case(1)
          call nduststrings('Enter grain density ',' (in g/cm^3)',varstring)
          simple_graindens = .true.
          do i = 1,ndusttypes
             call prompt(trim(varstring(i)),graindensinp(i),0.)
          enddo
          if (any(graindensinp /= graindensinp(1))) simple_graindens = .false.
       end select
    else
       call prompt('Enter dust to gas ratio',dust_to_gas,0.)
       if (use_dustgrowth) then
          message = 'Enter initial grain size in'
       else
          message = 'Enter grain size in'
       endif
       message = trim(message)//' '//trim(abbrev)
       grainsizeinp(1) = grainsizeinp(1)/units_to_cgs
       call prompt(trim(message),grainsizeinp(1),0.)
       grainsizeinp(1) = units_to_cgs*grainsizeinp(1)
       call prompt('Enter the intrinsic grain density in g/cm^3',graindensinp(1),0.)
    endif
    if (simple_grainsize) then
       print*,'Grainsize values chosen are equivalent to the ''equal'' case'
       print*,'   ...switching for a simpler setup file'
       io_grainsize = 1
    endif
    if (simple_graindens) then
       print*,'Grain density values chosen are equivalent to the ''equal'' case'
       print*,'   ...switching for a simpler setup file'
       io_graindens = 0
    endif
 case(2)
    call prompt('Enter dust to gas ratio',dust_to_gas,0.)
    if (use_dustgrowth) then
       message = 'Enter initial grain size in '//trim(abbrev)
       grainsizeinp(1) = grainsizeinp(1)/units_to_cgs
       call prompt(trim(message),grainsizeinp(1),0.)
       grainsizeinp(1) = units_to_cgs*grainsizeinp(1)
    else
       call prompt('Enter constant drag coefficient',K_code,0.)
    endif
 case default
 end select

end subroutine interactively_set_dust_full


!-----------------------------------------------------------------------
!+
!  Subroutine to select factor used in converting distance units
!+
!-----------------------------------------------------------------------
subroutine get_units_factor(units,factor)
 real,             intent(out) :: factor
 character(len=*), intent(in)  :: units
 select case(units)
 case('nm')
    factor = 1.e-7
 case('microns')
    factor = 1.e-4
 case('mm')
    factor = 1.e-1
 case('cm')
    factor = 1.
 case('m')
    factor = 1.e2
 case('km')
    factor = 1.e5
 case default
    stop 'Error! unrecognised unit'
 end select

end subroutine get_units_factor


!-----------------------------------------------------------------------
!+
!  Subroutine for reading dust properties from an input file
!  (originally written for setup_disc, but also needed for moddump/old setups)
!+
!-----------------------------------------------------------------------
subroutine read_dust_setup_options(db,nerr,dust_to_gas,df,gs,gd,isimple,imethod)
 use options,      only:use_dustfrac
 use dim,          only:use_dustgrowth
 use growth,       only:read_growth_setup_options
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use dust,         only:set_grainsize,grainsize,graindens
 use units,        only:udist,umass
 type(inopts), allocatable, intent(inout) :: db(:)
 integer, intent(out), optional :: imethod
 real,    intent(out), optional :: df(:),gs(:),gd(:)
 logical, intent(in),  optional :: isimple
 integer, intent(inout)         :: nerr
 real,    intent(out)           :: dust_to_gas
 integer            :: i,ierr
 real               :: grainsizeinp(ndusttypes)
 real               :: graindensinp(ndusttypes)
 real               :: dustfrac_percent(ndusttypes)
 logical            :: simple_grainsize = .false.
 logical            :: simple_graindens = .false.
 logical            :: simple_output    = .false.
 character(len=120) :: varlabel(ndusttypes)

 call read_inopt(dust_method,'dust_method',db,min=1,max=2,errcount=nerr)
 if (present(imethod)) imethod = dust_method

 if (use_dustfrac) then
    if (dust_method == 2) then
       print*,'Warning! use_dustfrac = .true. AND two-fluid dust are incompatible'
       print*,'   ...resetting use_dustfrac = .false.'
       use_dustfrac = .false.
    else
       call read_inopt(ilimitdustflux,'ilimitdustflux',db,err=ierr,errcount=nerr)
       if (ierr /= 0) nerr = nerr + 1
    endif
 endif

 call read_inopt(dust_to_gas,'dust_to_gas_ratio',db,min=0.,errcount=nerr)

 if (present(isimple)) simple_output = isimple
 if (.not.simple_output) then
    if (use_dustfrac .and. ndusttypes > 1) then
       call read_inopt(io_grainsize,'io_grainsize',db,min=0,max=2,errcount=nerr)
       select case(io_grainsize)
       case(0)
          call read_inopt(smincgs,'smincgs',db,min=0.,errcount=nerr)
          call read_inopt(smaxcgs,'smaxcgs',db,min=smincgs,errcount=nerr)
          call read_inopt(sindex ,'sindex' ,db,errcount=nerr)
       case(1)
          call read_inopt(grainsizeinp(1),'grainsizeinp',db,min=0.,err=ierr,errcount=nerr)
          if (ierr /= 0) then
             grainsizeinp(:) = 0.1
          else
             grainsizeinp(:) = grainsizeinp(1)
          endif
          grainsizecgs = grainsizeinp
       case(2)
          !--Make N grain size labels
          call nduststrings('grainsizeinp','',varlabel)
          do i = 1,ndusttypes
             call read_inopt(grainsizeinp(i),trim(varlabel(i)),db,min=0.,err=ierr,errcount=nerr)
             if (ierr /= 0) then
                call set_grainsize(smincgs,smaxcgs)
                grainsizeinp(:) = grainsizecgs(:)
             endif
          enddo
          !--Make N dust fraction labels
          call nduststrings('dustfrac','',varlabel)
          do i = 1,ndusttypes
             call read_inopt(dustfrac_percent(i),trim(varlabel(i)),db,min=0.,max=100.,err=ierr,errcount=nerr)
          enddo
          if (sum(dustfrac_percent(:)) /= 100.) then
             print*,'ERROR: dust fraction percentages need to add up to 100!'
             nerr = nerr+1
          endif
          simple_grainsize = .true.
          if (any(grainsizeinp /= grainsizeinp(1)) .or. &
              any(dustfrac_percent /= dustfrac_percent(1))) simple_grainsize = .false.
          if (simple_grainsize) then
             print*,'Grainsize values chosen are equivalent to the ''equal'' case'
             print*,'   ...switching for a simpler setup file'
             io_grainsize = 1
          endif
       end select
       call read_inopt(io_graindens,'io_graindens',db,min=0,errcount=nerr)
       select case(io_graindens)
       case(0)
          call read_inopt(graindensinp(1),'graindensinp',db,min=0.,err=ierr,errcount=nerr)
          if (ierr /= 0) then
             graindensinp(:) = 3.
          else
             graindensinp(:) = graindensinp(1)
          endif
          graindenscgs = graindensinp
       case(1)
          !--Make N grain size labels
          call nduststrings('graindensinp','',varlabel)
          do i = 1,ndusttypes
             call read_inopt(graindensinp(i),trim(varlabel(i)),db,min=0.,err=ierr,errcount=nerr)
             if (ierr /= 0) then
                graindensinp(i) = 2. + (i-1)*2./real(ndusttypes-1)
             endif
          enddo
          simple_graindens = .true.
          if (any(graindensinp /= graindensinp(1))) simple_graindens = .false.
          if (simple_graindens) then
             print*,'Grain density values chosen are equivalent to the ''equal'' case'
             print*,'   ...switching for a simpler setup file'
             io_graindens = 0
          endif
       end select
    else
       call read_inopt(grainsizeinp(1),'grainsizeinp',db,min=0.,err=ierr,errcount=nerr)
       if (ierr /= 0) then
          grainsizeinp = 0.1
       endif
       grainsizecgs = grainsizeinp
       call read_inopt(graindensinp(1),'graindensinp',db,min=0.,errcount=nerr)
       if (ierr /= 0) then
          graindensinp = 3.
       endif
       graindenscgs = graindensinp
    endif
    if (use_dustgrowth) call read_growth_setup_options(db,nerr)
 endif

 if (io_grainsize == 0) then
    grainsizeinp(:) = grainsizecgs(:)
 else
    grainsizecgs(:) = grainsizeinp(:)
 endif
 graindenscgs(:)    = graindensinp(:)

 grainsize(:) = grainsizecgs(:)/udist
 graindens(:) = graindenscgs(:)*udist**3/umass

 if (present(df)) df = dustfrac_percent
 if (present(gs)) gs = grainsizecgs
 if (present(gd)) gd = graindenscgs

end subroutine read_dust_setup_options


!-----------------------------------------------------------------------
!+
!  reads input dust options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_dust_infile_options(name,valstring,imatch,igotall,ierr)
 use units,          only:udist,umass
 use options,        only:use_dustfrac
 use dust,           only:icut_backreaction,grainsize,graindens
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 real(kind=8)  :: udens
 integer, save :: ngrainsize = 0
 integer, save :: ngraindens = 0
 integer, parameter :: nvars = 9
 integer, parameter :: narrs = 2
 integer, parameter :: nvalues = nvars + narrs*(ndusttypes-1)
 integer, parameter :: iidrag        = 1, &
                       ismin         = 2, &
                       ismax         = 3, &
                       isindex       = 4, &
                       iKcode        = 5, &
                       ibackreact    = 6, &
                       ilimitflux    = 7, &
                       !--dust arrays initial index
                       igrainsize    = 8, &
                       igraindens    = 9 + (ndusttypes-1), &
                       !--dust arrays final index
                       igrainsizeend = igraindens-1, &
                       igraindensend = nvars
 integer, save :: igot(nvalues)  = 0
 integer       :: ineed(nvalues)

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('idrag')
    read(valstring,*,iostat=ierr) idrag
    igot(iidrag) = 1
 case('grainsize')
    ngrainsize = ngrainsize + 1
    read(valstring,*,iostat=ierr) grainsizecgs(ngrainsize)
    grainsize(ngrainsize) = grainsizecgs(ngrainsize)/udist
    igot(igrainsize + ngrainsize - 1) = 1
 case('smin')
    read(valstring,*,iostat=ierr) smincgs
    igot(ismin) = 1
 case('smax')
    read(valstring,*,iostat=ierr) smaxcgs
    igot(ismax) = 1
 case('p')
    read(valstring,*,iostat=ierr) sindex
    igot(isindex) = 1
 case('graindens')
    ngraindens = ngraindens + 1
    read(valstring,*,iostat=ierr) graindenscgs(ngraindens)
    udens = umass/udist**3
    graindens(ngraindens) = graindenscgs(ngraindens)/udens
    igot(igraindens + ngraindens - 1) = 1
 case('K_code')
    read(valstring,*,iostat=ierr) K_code
    igot(iKcode) = 1
 case('icut_backreaction')
    read(valstring,*,iostat=ierr) icut_backreaction
    igot(ibackreact) = 1
 case('ilimitdustflux')
    read(valstring,*,iostat=ierr) ilimitdustflux
    igot(ilimitflux) = 1
 case default
    imatch = .false.
 end select

 ineed = 0

 !--Parameters needed by all combinations
 ineed(iidrag)     = 1
 ineed(ibackreact) = 1

 !--Parameters specific to particular setups
 select case(idrag)
 case(1)
    if ((ndusttypes == 1) .or. (io_grainsize == 1 .and. io_graindens == 0)) then
       ineed(igrainsize) = 1
       ineed(igraindens) = 1
    endif
 case(2,3)
    ineed(iKcode) = 1
 case default
    stop 'Error! Invalid idrag option passed to read_dust_infile_options'
 end select

 if (use_dustfrac) ineed(ilimitflux) = 1

 !--Check that we have just the *necessary* parameters
 if (all(ineed == igot)) igotall = .true.

end subroutine read_dust_infile_options


!-----------------------------------------------------------------------
!+
!  Wrapper for writing dust properties to grains.tmp file
!  (needed for moddump and old setups)
!+
!-----------------------------------------------------------------------
subroutine extract_dust_header(multidustdump,hdr,ierr)
 use dim,             only:ndusttypes
 use dump_utils,      only:extract,dump_h
 use units,           only:udist,umass
 use options,         only:use_dustfrac
 integer,      intent(inout) :: ierr
 logical,      intent(in)    :: multidustdump
 type(dump_h), intent(in)    :: hdr

 integer :: i,ierr1,ierrs(max(10,ndusttypes)),nerr
 integer :: io_old
 real    :: smin,smax,sind,grainsizeinp(ndusttypes),graindensinp(ndusttypes)
 real    :: dust_to_gas
 real    :: dustfrac_percent(ndusttypes) = 0.
 logical :: missingdust = .false.
 character(len=120) :: varlabel(ndusttypes)

 if (use_dustfrac .and. multidustdump) then
    call extract('io_grainsize',io_grainsize,hdr,ierr1)
    if (ierr1 /= 0) then
       io_grainsize = 0
       missingdust  = .true.
    endif
    call extract('io_graindens',io_graindens,hdr,ierr1)
    if (ierr1 /= 0) then
       io_graindens = 0
       missingdust  = .true.
    endif
 endif

 if (missingdust) then
    call write_temp_grains_file(dust_to_gas,dustfrac_percent,imethod=dust_method,ireadwrite=1)

    !--save the input values from the dust file
    smin = smincgs
    smax = smaxcgs
    sind = sindex
    grainsizeinp = grainsizecgs
    graindensinp = graindenscgs

    !--extract all possible dust options from the dump file and compare them with
    !  new values calculated from the temporary dust file
    nerr = 0

    print*,'Checking grain sizes...'
    call extract('smincgs',smincgs,hdr,ierrs(1))
    call check_dust_value(ierrs(1),nerr,'smincgs',smin,smincgs)

    call extract('smaxcgs',smaxcgs,hdr,ierrs(2))
    call check_dust_value(ierrs(2),nerr,'smaxcgs',smax,smaxcgs)

    call extract('sindex' ,sindex,hdr,ierrs(3))
    call check_dust_value(ierrs(3),nerr,'sindex',sind,sindex)

    io_old = io_grainsize
    if (all(ierrs(1:3) == 0)) then
       print*,'   ...smincgs, smaxcgs, and sindex found'
       io_grainsize = 0
    else
       call extract('grainsize',grainsizecgs(1),hdr,ierrs(1))
       if (ierrs(1) == 0) then
          print*,'   ...grain size found'
          grainsizecgs = grainsizecgs(1)*udist
          io_grainsize = 1
          call check_dust_value(ierrs(1),nerr,'grainsize',grainsizeinp(1),grainsizecgs(1))
       else
          print*,'Warning! No grain size information found in dump file!'
          print*,'         Carefully check your grain sizes by hand....'
       endif
    endif
    if (io_old /= io_grainsize) nerr = nerr + 1

    print*,'Checking grain densities...'
    call extract('graindens',graindenscgs(1),hdr,ierrs(1))
    io_old = io_graindens
    if (ierrs(1) == 0) then
       print*,'   ...grain density found'
       graindenscgs = graindenscgs(1)*umass/udist**3
       io_graindens = 0
       call check_dust_value(ierrs(1),nerr,'graindens',graindensinp(1),graindenscgs(1))
    else
       print*,'Warning! No grain density information found in dump file!'
       print*,'         Carefully check your grain density by hand....'
    endif

    if (io_old /= io_graindens) nerr = nerr + 1

    call write_temp_grains_file(dust_to_gas,dustfrac_percent,imethod=dust_method,ireadwrite=2)
    if (nerr /= 0) then
       print*,'ERROR: inconsistent dust options! Rewriting dust file...'
       print "(/,a)",' >>> Try re-running the executable <<<'
       ierr = 5
    endif
 else
    if (io_grainsize == 0) then
       !--grainsize is reconstructed from these three quantities
       call extract('smincgs',smin,hdr,ierrs(1))
       call extract('smaxcgs',smax,hdr,ierrs(2))
       call extract('sindex' ,sind,hdr,ierrs(3))
       if (any(ierrs(1:3) /= 0) .and. idrag == 1) then
       else
          smincgs = smin
          smaxcgs = smax
          sindex  = sind
       endif
       write(*,*) 'Grain-size distribution properties, [smin,smax,sindex] = ',smincgs,smaxcgs,sindex
    elseif (io_grainsize == 1) then
       !--grainsize is the same for each species
       call extract('grainsizecgs',grainsizeinp(1),hdr,ierrs(1))
       if (ierrs(1) /= 0) then
          print*,'Warning: constant grainsize not found...using default values!'
       else
          grainsizecgs = grainsizeinp(1)
       endif
       write(*,*) 'Grain size (cm; same for all species) =',grainsizecgs(1)
    elseif (io_grainsize == 2) then
       write(*,*) 'Using a custom grain-size distribution:'
       call nduststrings('grainsizecgs','',varlabel)
       do i = 1,ndusttypes
          call extract(trim(varlabel(i)),grainsizeinp(i),hdr,ierrs(1))
          if (ierrs(1) /= 0) then
             print*,'Warning: user defined grainsize not found...using default value!'
          else
             grainsizecgs(i) = grainsizeinp(i)
          endif
          write(*,*) '   grain size (cm) =',grainsizecgs(i)
       enddo
    endif
    if (io_graindens == 0) then
       !--graindens is the same for each species
       call extract('graindenscgs',graindensinp(1),hdr,ierr)
       if (ierrs(1) /= 0) then
          print*,'Warning: constant graindens not found...using default values!'
       else
          graindenscgs = graindensinp(1)
       endif
       write(*,*) 'Intrinsic grain density (g/cm^3; same for all species) =',graindenscgs(1)
    elseif (io_graindens == 1) then
       write(*,*) 'Using a custom intrinsic grain density for the dust:'
       call nduststrings('graindenscgs','',varlabel)
       do i = 1,ndusttypes
          call extract(trim(varlabel(i)),graindensinp(i),hdr,ierrs(1))
          if (ierrs(1) /= 0) then
             print*,'Warning: user defined graindens not found...using default value!'
          else
             graindenscgs(i) = graindensinp(i)
          endif
          write(*,*) '   grain density (g/cm^3) =',graindenscgs(i)
       enddo
    endif
 endif

end subroutine extract_dust_header


!-----------------------------------------------------------------
!+
!  Check if new dust options are compatible with old dump files
!+
!-----------------------------------------------------------------
subroutine check_dust_value(ierr,nerr,tag,newval,oldval)
 integer, intent(inout) :: nerr
 integer, intent(in)    :: ierr
 real,    intent(in)    :: newval,oldval
 character(len=*), intent(in) :: tag

 real :: tol = 1.e-5

 if (ierr == 0) then
    if (abs((newval-oldval)/oldval) > tol) then
       print*,'ERROR: '//tag//' is not the same as in dump file'
       print*,'Change ',newval,'to ',oldval,'in temparary dust file'
       nerr = nerr + 1
    endif
 endif

end subroutine check_dust_value



!-----------------------------------------------------------------------
!+
!  Writing dust properties to dump file header
!+
!-----------------------------------------------------------------------
subroutine write_dust_to_header(multidustdump,hdr,ierr)
 use dim,             only:ndusttypes,use_dustgrowth
 use dump_utils,      only:add_to_iheader,add_to_rheader,dump_h
 integer,      intent(inout) :: ierr
 type(dump_h), intent(inout) :: hdr
 logical,      intent(in)    :: multidustdump

 integer :: i
 character(len=120) :: varlabel(ndusttypes)

 ! write dust information
 if (use_dustgrowth) then
    write(*,*) 'writing dust growth properties to header'
 elseif (multidustdump) then
    write(*,*) 'writing multi-grain properties to header'
 else
    write(*,*) 'writing graindens and grainsize to header'
 endif
 if (multidustdump) then
    call add_to_iheader(ndusttypes,  'ndusttypes',  hdr,ierr)
    call add_to_iheader(io_grainsize,'io_grainsize',hdr,ierr)
    call add_to_iheader(io_graindens,'io_graindens',hdr,ierr)

    !--grain size
    if (io_grainsize == 0) then
       !--grainsize is reconstructed from these three quantities
       call add_to_rheader(smincgs,'smincgs',hdr,ierr)
       call add_to_rheader(smaxcgs,'smaxcgs',hdr,ierr)
       call add_to_rheader(sindex, 'sindex', hdr,ierr)
    elseif (io_grainsize == 1) then
       !--grainsize is the same for each species
       call add_to_rheader(grainsizecgs(1),'grainsizecgs',hdr,ierr)
    elseif (io_grainsize == 2) then
       !--grainsize is defined by the user at setup
       call nduststrings('grainsizecgs','',varlabel)
       do i = 1,ndusttypes
          call add_to_rheader(grainsizecgs(i),trim(varlabel(i)),hdr,ierr)
       enddo
    endif
    !--grain density
    if (io_graindens == 0) then
       !--graindens is the same for each species
       call add_to_rheader(graindenscgs(1),'graindenscgs',hdr,ierr)
    elseif (io_graindens == 1) then
       !--graindens is defined by the user at setup
       call nduststrings('graindenscgs','',varlabel)
       do i = 1,ndusttypes
          call add_to_rheader(graindenscgs(i),trim(varlabel(i)),hdr,ierr)
       enddo
    endif
 else
    call add_to_rheader(grainsizecgs(1),'grainsizecgs',hdr,ierr)
    call add_to_rheader(graindenscgs(1),'graindenscgs',hdr,ierr)
 endif

end subroutine write_dust_to_header

!-----------------------------------------------------------------------
!+
!  Subroutine for deciding wheather to use one-fluid or two-fluid dust
!+
!-----------------------------------------------------------------------
subroutine check_dust_method(id,filename,ichange_method)
 use options, only:use_dustfrac
 use dust,    only:init_drag,get_ts,grainsize,graindens,idrag
 use part,    only:npart,massoftype,xyzh,vxyzu,rhoh,igas,dustfrac
 use eos,     only:ieos,get_spsound
 use io,      only:master
 character(len=*), intent(in)  :: filename
 integer,          intent(in)  :: id
 logical,          intent(out) :: ichange_method
 integer :: i,l,iregime,ierr,icheckdust
 real    :: r,rhogasi,rhodusti,rhoi,dustfracisum,spsoundi
 real    :: dustfraci(ndusttypes),tsi(ndusttypes)
 character(len=120) :: string
 logical :: iforce_dust_method = .false.

 call init_drag(ierr)

 icheckdust = 0
 do i = 1,npart
    r = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
    if (use_dustfrac) then
       rhoi = rhoh(xyzh(4,i),massoftype(igas))
       dustfraci(:) = dustfrac(:,i)
       dustfracisum = sum(dustfraci(:))
       rhogasi      = rhoi*(1.-dustfracisum)
       spsoundi = get_spsound(ieos,xyzh(:,i),rhogasi,vxyzu(:,i))
       do l = 1,ndusttypes
          rhodusti = rhoi*dustfraci(l)
          call get_ts(idrag,grainsize(l),graindens(l),rhogasi,rhodusti,spsoundi,0.,tsi(l),iregime)
       enddo
       if (any(tsi(:) > xyzh(4,i)/spsoundi)) icheckdust = icheckdust + 1
    endif
 enddo

 call get_environment_variable('IFORCE_DUST_METHOD',string)
 if (trim(string)=='yes') iforce_dust_method = .true.

 ichange_method = .false.
 if (real(icheckdust)/real(npart) > 0.1 .and. .not.iforce_dust_method) then
    if (dust_method == 1) then
       use_dustfrac = .false.
       ichange_method = .true.
       dust_method = 2
    endif

    print*,''
    print*,'*******************************************************************************'
    print*,'WARNING! More than 10% of the particles have a Stokes Number larger than the'
    print*,'threshold under which the terminal velocity approximation is valid. We suggest'
    print*,'you switch to using two-fluid. If you absolutely insist on using the one-fluid,'
    print*,'you can set the environment variable IFORCE_DUST_METHOD=yes and rerun the setup.'
    print*,'*******************************************************************************'
    print*,''
 elseif (iforce_dust_method) then
    print*,''
    print*,'*******************************************************************************'
    print*,'WARNING! You have chosen to manually select the dust method. Care should be taken'
    print*,'to ensure that you do not violate the terminal velocity approximation.'
    print*,'*******************************************************************************'
    print*,''
 endif

end subroutine check_dust_method



!-----------------------------------------------------------------------
!+
!  Subroutine for writing dust properties to an input file
!  (originally written for setup_disc, but also needed for moddump/old setups)
!+
!-----------------------------------------------------------------------
subroutine write_dust_setup_options(iunit,dust_to_gas,df,gs,gd,imethod,iprofile,isimple)
 use dim,          only:use_dustgrowth
 use growth,       only:write_growth_setup_options
 use options,      only:use_dustfrac
 use infile_utils, only:write_inopt
 use dust,         only:grainsizecgs,graindenscgs
 integer, intent(in), optional :: imethod,iprofile
 real,    intent(in), optional :: df(:),gs(:),gd(:)
 logical, intent(in), optional :: isimple
 integer, intent(in)    :: iunit
 real,    intent(in)    :: dust_to_gas

 integer :: i
 integer :: profile_set_dust = -1
 real    :: grainsizeinp(ndusttypes)
 real    :: graindensinp(ndusttypes)
 real    :: dustfrac_percent(ndusttypes) = 0.
 logical :: simple_output = .false.
 character(len=120) :: varlabel(ndusttypes),varstring(ndusttypes)

 grainsizeinp = grainsizecgs
 graindensinp = graindenscgs

 if (present(df)) dustfrac_percent = df
 if (present(gs)) grainsizeinp = gs
 if (present(gd)) graindensinp = gd
 if (present(imethod)) dust_method = imethod
 if (present(iprofile)) profile_set_dust = iprofile
 if (present(isimple)) simple_output = isimple

 write(iunit,"(/,a)") '# options for dust'
 call write_inopt(dust_method,'dust_method','dust method (1=one fluid,2=two fluid)',iunit)
 if (use_dustfrac) call write_inopt(ilimitdustflux,'ilimitdustflux', &
                                    'limit dust diffusion using Ballabio et al. (2018)',iunit)
 call write_inopt(dust_to_gas,'dust_to_gas_ratio','dust to gas ratio',iunit)
 if (profile_set_dust /= -1) then
    call write_inopt(profile_set_dust,'profile_set_dust', &
       'how to set dust density profile (0=equal to gas,1=custom,2=equal to gas with cutoffs)',iunit)
 endif
 if (.not.simple_output) then
    if (use_dustgrowth) then
       call write_inopt(grainsizeinp(1),'grainsizeinp','initial grain size (in cm)',iunit)
       call write_inopt(graindensinp(1),'graindensinp','intrinsic grain density (in g/cm^3)',iunit) ! Modify this is graindens becomes variable
    elseif (use_dustfrac .and. ndusttypes > 1) then
       call write_inopt(io_grainsize,'io_grainsize', &
          'grain size distribution (0=power-law,1=equal,2=manually)',iunit)
       select case(io_grainsize)
       case(0)
          call write_inopt(smincgs,'smincgs','min grain size (in cm)',iunit)
          call write_inopt(smaxcgs,'smaxcgs','max grain size (in cm)',iunit)
          call write_inopt(sindex ,'sindex' ,'grain size power-law index (e.g. MRN = 3.5)',iunit)
       case(1)
          call write_inopt(grainsizeinp(1),'grainsizeinp','grain size (in cm)',iunit)
       case(2)
          !--Make N grain size labels
          call nduststrings('grainsizeinp','',varlabel)
          call nduststrings('grain size ',' (in cm)',varstring)
          do i = 1,ndusttypes
             call write_inopt(grainsizeinp(i),trim(varlabel(i)),trim(varstring(i)),iunit)
          enddo
          !--Make N dust fraction labels
          call nduststrings('dustfrac','',varlabel)
          call nduststrings('dust fraction ',' (in % of total dust)',varstring)
          if (all(dustfrac_percent == 0.) .or. sum(dustfrac_percent) /= 100.) dustfrac_percent(:) = 100./ndusttypes
          do i = 1,ndusttypes
             call write_inopt(dustfrac_percent(i),trim(varlabel(i)),trim(varstring(i)),iunit)
          enddo
       end select
       call write_inopt(io_graindens,'io_graindens','grain density input (0=equal,1=manually)',iunit)
       select case(io_graindens)
       case(0)
          call write_inopt(graindensinp(1),'graindensinp','intrinsic grain density (in g/cm^3)',iunit)
       case(1)
          !--Make N grain density labels
          call nduststrings('graindensinp','',varlabel)
          call nduststrings('grain density ',' (in g/cm^3)',varstring)
          do i = 1,ndusttypes
             call write_inopt(graindensinp(i),trim(varlabel(i)),trim(varstring(i)),iunit)
          enddo
       end select
    else
       call write_inopt(grainsizeinp(1),'grainsizeinp','grain size (in cm)',iunit)
       call write_inopt(graindensinp(1),'graindensinp','intrinsic grain density (in g/cm^3)',iunit) ! Modify this is graindens becomes variable
    endif
    if (use_dustgrowth) call write_growth_setup_options(iunit)
 endif

end subroutine write_dust_setup_options


!-----------------------------------------------------------------------
!+
!  writes input dust options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_dust_infile_options(iunit)
 use dim,          only:use_dustgrowth
 use infile_utils, only:write_inopt
 use options,      only:use_dustfrac
 use dust,         only:icut_backreaction
 integer, intent(in) :: iunit
 character(len=10) :: numdust

 write(numdust,'(I10)') ndusttypes
 write(iunit,"(/,a)") '# options controlling dust ('//trim(adjustl(numdust))//' dust species)'

 call write_inopt(idrag,'idrag','gas/dust drag (0=off,1=Epstein/Stokes,2=const K,3=const ts)',iunit)

 if (ndusttypes > 1) then
    !--the grainsize (and powerlaw index) should be set in the setup file
    !--the intrinsic grain density should be set in the setup file
    select case(idrag)
    case(1)
       print "(/,a)",'*************************************************************************'
       print "(a)",  '*************************************************************************'
       print "(/,a)",'Warning! Grain size/density are now set during setup when ndusttypes > 1 '
       print "(a)",  '         and only limited setups (e.g. dustydisc) support this ability.  '
       print "(a)",  '         If not using one of these setups, switch to using idrag = [2,3].'
       print "(/,a)",'*************************************************************************'
       print "(a)",  '*************************************************************************'
    case(2,3)
       call write_inopt(K_code,'K_code','drag constant when constant drag is used',iunit)
    end select

 else
    select case(idrag)
    case(1)
       if (use_dustgrowth) then
          call write_inopt(grainsizecgs(ndusttypes),'grainsize','Initial grain size in cm',iunit)
       else
          call write_inopt(grainsizecgs(ndusttypes),'grainsize','Grain size in cm',iunit)
       endif
       call write_inopt(graindenscgs(ndusttypes),'graindens','Intrinsic grain density in g/cm^3',iunit)
    case(2,3)
       call write_inopt(K_code,'K_code','drag constant when constant drag is used',iunit)
    end select
 endif

 call write_inopt(icut_backreaction,'icut_backreaction','cut the drag on the gas phase (0=no, 1=yes)',iunit)

 if (use_dustfrac) then
    call write_inopt(ilimitdustflux,'ilimitdustflux','limit the dust flux using Ballabio et al. (2018)',iunit)
 endif

end subroutine write_dust_infile_options


!-----------------------------------------------------------------------
!+
!  Wrapper for writing dust properties to grains.tmp file
!  (needed for moddump and old setups)
!+
!-----------------------------------------------------------------------
subroutine write_temp_grains_file(dust_to_gas,dustfrac_percent,imethod,iprofile,ireadwrite)
 use dim,          only:ndusttypes,use_dustgrowth
 use io,           only:id,master
 use infile_utils, only:open_db_from_file,close_db,inopts
 use options,      only:use_dustfrac
 real,    intent(out) :: dust_to_gas
 real,    intent(out) :: dustfrac_percent(:)
 integer, intent(in),  optional :: iprofile
 integer, intent(in),  optional :: ireadwrite
 integer, intent(out), optional :: imethod

 integer, parameter :: iunit = 20
 integer            :: nerr = 0
 !integer            :: dust_method = -1
 integer            :: profile_set_dust = -1
 integer            :: ioselect = 0
 real               :: grainsizeinp(ndusttypes),graindensinp(ndusttypes)
 logical            :: iexist
 character(len=10)  :: grainsfile = 'grains.tmp'
 type(inopts), allocatable :: db(:)

 print "(/,a)",' Warning: missing essential dust options. Checking if '//trim(grainsfile)//' exists...'

 inquire(file=grainsfile,exist=iexist)
 if (present(imethod)) dust_method = imethod
 if (present(iprofile)) profile_set_dust = iprofile
 if (present(ireadwrite) .and. iexist) ioselect = ireadwrite

 select case(ioselect)
 case(0)
    !
    !--Read and write the dust file and check for errors
    !
    dustfrac_percent = 0.
    inquire(file=grainsfile,exist=iexist)
    if (iexist) then
       print*,trim(grainsfile)//' found...opening file and reading dust parameters...'
       !--read from grains.tmp file
       call open_db_from_file(db,grainsfile,iunit,nerr)
       call read_dust_setup_options(db,nerr,dust_to_gas,df=dustfrac_percent,gs=grainsizeinp, &
                                    gd=graindensinp,imethod=dust_method)
       call close_db(db)
       if (dust_method == 2 .and. use_dustfrac) use_dustfrac = .false.
       if (id==master) then
          open(unit=iunit,file=grainsfile,status='replace',form='formatted')
          call write_dust_setup_options(iunit,dust_to_gas,df=dustfrac_percent,gs=grainsizeinp, &
                                        gd=graindensinp,imethod=dust_method,iprofile=profile_set_dust)
          close(iunit)
       endif
       if (nerr /= 0) then
          print*,'ERROR: inconsistent dust options! Will try to rewrite '//trim(grainsfile)//'...'
          print "(/,a)",' >>> Try re-running the executable <<<'
          stop
       else
          print "(a,/)",' Finished reading '//trim(grainsfile)
       endif
    elseif (id == master) then
       ! interactive setup
       print "(a,/)",' '//trim(grainsfile)//' not found: using interactive setup'
       dust_to_gas = 0.01
       grainsizeinp = grainsizecgs
       graindensinp = graindenscgs
       call interactively_set_dust(dust_to_gas,dustfrac_percent,grainsizeinp,graindensinp, &
                                   imethod=dust_method,iprofile=profile_set_dust)
       !
       !--write default input file
       !
       open(unit=iunit,file=grainsfile,status='replace',form='formatted')
       call write_dust_setup_options(iunit,dust_to_gas,df=dustfrac_percent,gs=grainsizeinp, &
                                     gd=graindensinp,imethod=dust_method,iprofile=profile_set_dust)
       if (use_dustgrowth) print*,'--> Setting default growth values'
       close(iunit)
       print "(/,a)",' >>> please edit '//trim(grainsfile)//' to set parameters for your problem then rerun <<<'
       stop
    else
       stop 'Something is not right...'
    endif
 case(1)
    !
    !--Read the dust file and check for errors
    !
    call open_db_from_file(db,grainsfile,iunit,nerr)
    call read_dust_setup_options(db,nerr,dust_to_gas,df=dustfrac_percent,gs=grainsizeinp, &
                                 gd=graindensinp)
    call close_db(db)
 case(2)
    !
    !--Write the dust file
    !
    if (id==master) then
       open(unit=iunit,file=grainsfile,status='replace',form='formatted')
       call write_dust_setup_options(iunit,dust_to_gas,df=dustfrac_percent,gs=grainsizecgs, &
                                     gd=graindenscgs,imethod=dust_method,iprofile=profile_set_dust)
       close(iunit)
    endif
 case default
    stop 'Invalid ireadwrite option in write_temp_grains_file'
 end select
 if (present(imethod)) imethod = dust_method

end subroutine write_temp_grains_file

end module readwrite_dust
