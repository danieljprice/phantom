!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: set_dust
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
!    dust_method       -- dust method (1=one fluid,2=two fluid)
!    dust_to_gas_ratio -- dust to gas ratio
!    graindensinp      -- intrinsic grain density (in g/cm^3)
!    grainsizeinp      -- grain size (in cm)
!    io_graindens      -- grain density input (0=equal,1=manually)
!    sindex            -- grain size power-law index (e.g. MRN = 3.5)
!    smaxcgs           -- max grain size (in cm)
!    smincgs           -- min grain size (in cm)
!
!  DEPENDENCIES: dim, dust, eos, growth, infile_utils, io, options, part,
!    prompting, units
!+
!--------------------------------------------------------------------------

module set_dust
 use dim,  only:ndusttypes
 use dust, only:smincgs,smaxcgs,sindex,grainsizecgs,graindenscgs,K_code,idrag, &
                ilimitdustflux
 implicit none
 !--Default values for the dust in the infile
 integer, public :: io_grainsize = 0
 integer, public :: io_graindens = 0
 logical, public :: ichange_method = .false.

 public :: set_dustfrac_from_inopts
 public :: nduststrings
 ! generic interface to interactively_set_dust
 interface interactively_set_dust
  module procedure interactively_set_dust_simple, interactively_set_dust_full
 end interface interactively_set_dust
 public :: interactively_set_dust
 public :: read_dust_setup_options
 public :: write_dust_setup_options
 public :: write_temp_grains_file
 public :: check_dust_method

 private

contains

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
 real,    intent(out)           :: dust_to_gas
 logical, intent(in),  optional :: Kdrag
 integer, intent(out), optional :: imethod
 character(len=*), intent(in), optional :: units

 integer :: dust_method = -1
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

 if (dust_method /= -1) then
    !
    !--dust method
    !
    if (ndusttypes > 1) dust_method  = 1

    call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
    if (dust_method == 1) then
       use_dustfrac = .true.
    else
       use_dustfrac = .false.
       if (ndusttypes > 1) call fatal('setup','dust_method=2 is currently only compatible with ndusttypes=1!')
    endif
    if (present(imethod)) imethod = dust_method
 endif

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
    if (ndusttypes > 1) dust_method  = 1

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
          call grainsize_dist_assist(units_to_cgs,abbrev)

          !message = 'Enter minimum grain size in '//trim(abbrev)
          !smincgs = smincgs/units_to_cgs
          !call prompt(trim(message),smincgs,0.)
          !smincgs = units_to_cgs*smincgs

          !message = 'Enter maximum grain size in '//trim(abbrev)
          !smaxcgs = smaxcgs/units_to_cgs
          !call prompt(trim(message),smaxcgs,smincgs)
          !smaxcgs = units_to_cgs*smaxcgs

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
 case('micron')
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
 integer            :: dust_method = -1
 real               :: tol = 1.e-5
 real               :: percentsum
 real               :: grainsizeinp(ndusttypes)
 real               :: graindensinp(ndusttypes)
 real               :: dustfrac_percent(ndusttypes)
 logical            :: simple_grainsize = .false.
 logical            :: simple_graindens = .false.
 logical            :: simple_output    = .false.
 character(len=120) :: varlabel(ndusttypes),string(3)

 call read_inopt(dust_method,'dust_method',db,min=1,max=2,errcount=nerr)
 if (present(imethod)) imethod = dust_method

 string(1) = 'Warning! dust_method and use_dustfrac are currently incompatible. This is'
 string(2) = '         normal behaviour when using the moddump utility, but otherwise'
 string(3) = '         may indicate an initialisation problem for one or both variables.'
 if (dust_method == 1 .and. .not.use_dustfrac) then
    write(*,'(/a,/a,/a/)') (trim(string(i)), i = 1,3)
    print*,'...setting use_dustfrac = .true.'
    use_dustfrac = .true.
 elseif (dust_method == 2 .and. use_dustfrac) then
    write(*,'(/a,/a,/a/)') (trim(string(i)), i = 1,3)
    print*,'...setting use_dustfrac = .false.'
    use_dustfrac = .false.
 endif

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
          percentsum = sum(dustfrac_percent(:))
          if (0.01*abs(percentsum - 100.) > tol) then
             print*,'ERROR: dust fraction percentages only add up to:',percentsum
             print*,'sum(dustfrac_percent) = ',sum(dustfrac_percent(:))
             nerr = nerr+1
          elseif (0.01*abs(percentsum - 100.) < tol) then
             print*,'Warning: dust fraction percentages only add up to:',percentsum
             print*,'...adding the difference to the largest dust fraction'
             where (dustfrac_percent == maxval(dustfrac_percent))
                dustfrac_percent = dustfrac_percent + (100.-percentsum)
             endwhere
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
!  Subroutine for deciding wheather to use one-fluid or two-fluid dust
!+
!-----------------------------------------------------------------------
subroutine check_dust_method(id,filename,dust_method,ichange_method)
 use options, only:use_dustfrac
 use dust,    only:init_drag,get_ts,grainsize,graindens,idrag
 use part,    only:npart,massoftype,xyzh,vxyzu,rhoh,igas,dustfrac
 use eos,     only:ieos,get_spsound
 use io,      only:master
 integer,          intent(in)    :: id
 integer,          intent(inout) :: dust_method
 logical,          intent(out)   :: ichange_method
 character(len=*), intent(in)    :: filename
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
 integer :: dust_method = -1
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
 if (dust_method /= -1) then
    call write_inopt(dust_method,'dust_method','dust method (1=one fluid,2=two fluid)',iunit)
 endif
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
 integer            :: dust_method = -1
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



!------------------------------------------------------
!+
!  Functions used to determine the ideal min/max sizes
!  and spacing betweeen grain sizes for a simulation
!+
!------------------------------------------------------
real function smin_func(N,s1,sN)
 integer, intent(in)  :: N
 real,    intent(in)  :: s1,sN
 smin_func = s1*(s1/sN)**(1./(2.*real(N-1)))
end function smin_func

real function smax_func(N,smin,s1,sN)
 integer, intent(in)  :: N
 real,    intent(in)  :: smin,s1,sN
 smax_func = smin*(sN/s1)**(real(N)/real(N-1))
end function smax_func

real function s1_func(N,sN,logds)
 integer, intent(in)  :: N
 real,    intent(in)  :: sN,logds
 s1_func = sN*10.**(-real(N-1)*logds)
end function s1_func

real function sN_func(N,s1,logds)
 integer, intent(in)  :: N
 real,    intent(in)  :: s1,logds
 sN_func = s1*10.**(real(N-1)*logds)
end function sN_func

real function logds_func(N,smin,smax)
 integer, intent(in)  :: N
 real,    intent(in)  :: smin,smax
 logds_func = log10(smax/smin)/real(N)
end function logds_func
!------------------------------------------------------



!------------------------------------------------------
!+
!  Utility to help users choose more intuitive grain
!  size distribution parameters for a simulation
!+
!------------------------------------------------------
subroutine grainsize_dist_assist(units_to_cgs,abbrev)
 use draw_ascii_art, only:draw_vesta,draw_grainsize_dist_assist
 use prompting,      only:prompt
 real,              intent(in) :: units_to_cgs
 character(len=10), intent(in) :: abbrev
 character(len=120) :: message
 real,  allocatable :: loggrid(:),grid(:),s(:)
 real      :: s1,sN,smin,smax,logds
 logical   :: iloop = .true.
 logical   :: itryagain
 integer   :: i,N,ichoosedist

 call draw_vesta
 call draw_grainsize_dist_assist

 print*,''
 print*,'The grain size is selected based on a few important quantites:'
 print*,'   N    :: number of discrete dust phases'
 print*,'   s1   :: smallest simulated grain size'
 print*,'   sN   :: largest simulated grain size'
 print*,'   dS   :: spacing between grain sizes in log10 space'
 print*,'   smin :: minimum grain size in the grain-size distribution'
 print*,'           (used for calculating the mass, i.e. not simulated)'
 print*,'   smax :: maximum grain size in the grain-size distribution'
 print*,'           (used for calculating the mass, i.e. not simulated)'
 print*,''
 print*,'The following is an illustration of the grid in log10 space:'
 print*,''
 print*,'          dS                                                                   '
 print*,'  |<-------------->|                                                           '
 print*,' G(1)             G(2)                                   G(N)            G(N+1)'
 print*,'  |------S(1)------|------S(2)------ // ------S(N-1)------|------S(N)------|   '
 print*,' Smin     |<-------------->|                                             Smax  '
 print*,'                  dS                                                           '
 print*,''
 print*,'where capital letters denote log10 quantities and lowercase actual CGS values'
 print*,'   for example :: S(i) = log10(s(i))'
 print*,'G(i) are the cell edges (important for calculating the mass of simulated grains'
 print*,'   for example :: m(1) is the integrated mass between g(1) and g(2)'
 print*,'Note that S(i) is the midpoint of each log10 grid cell, i.e. s1 = sqrt(g(1)*g(2)),'
 print*,'   which weights the effective grain size to smaller values, reflecting the fact'
 print*,'   that the number density of a physical distribution decreases with grain size'
 print*,''

 do while(iloop)
    !--initialise/reset the values for prompts
    itryagain   = .false.
    ichoosedist = 1
    N           = ndusttypes
    smin        = smincgs/units_to_cgs
    smax        = smaxcgs/units_to_cgs
    logds       = logds_func(N,smin,smax)
    s1          = smin
    sN          = smax

    call prompt('Which quantities do you want to define?'// &
                ' (values <= 0 allow you to trial different N)'//new_line('A')// &
                ' 0 = s(1), s(N), and dS'//new_line('A')// &
                ' 1 = s(1) and s(N)'//new_line('A')// &
                ' 2 = s(1) and dS'//new_line('A')// &
                ' 3 = s(N) and dS'//new_line('A')// &
                ' 4 = smin and smax'//new_line('A'),ichoosedist,-4,4)

    select case(abs(ichoosedist))
    case(0)
       message = 'Choose the smallest simulated grain size s(1) in '//trim(abbrev)
       call prompt(trim(message),s1,0.)
       message = 'Choose the largest simulated grain size s(N) in '//trim(abbrev)
       call prompt(trim(message),sN,s1,1000.)
       message = 'Choose the cell width in log10 space dS assuming units of '//trim(abbrev)
       call prompt(trim(message),logds,0.01,100.)
       N    = int(1./logds*log10(sN/s1)) + 1
       sN   = sN_func(N,s1,logds)
       smin = smin_func(N,s1,sN)
       smax = smax_func(N,smin,s1,sN)
    case(1)
       if (ichoosedist < 0) then
          call prompt('Choose the number of discrete dust phases (N)',N,0)
       endif
       message = 'Choose the smallest simulated grain size s(1) in '//trim(abbrev)
       call prompt(trim(message),s1,0.)
       message = 'Choose the largest simulated grain size s(N) in '//trim(abbrev)
       call prompt(trim(message),sN,s1,1000.)
       smin  = smin_func(N,s1,sN)
       smax  = smax_func(N,smin,s1,sN)
       logds = logds_func(N,smin,smax)
    case(2)
       if (ichoosedist < 0) then
          call prompt('Choose the number of discrete dust phases (N)',N,0)
       endif
       message = 'Choose the smallest simulated grain size s(1) in '//trim(abbrev)
       call prompt(trim(message),s1,0.)
       message = 'Choose the cell width in log10 space dS assuming units of '//trim(abbrev)
       call prompt(trim(message),logds,0.01,100.)
       sN   = sN_func(N,s1,logds)
       smin = smin_func(N,s1,sN)
       smax = smax_func(N,smin,s1,sN)
    case(3)
       if (ichoosedist < 0) then
          call prompt('Choose the number of discrete dust phases (N)',N,0)
       endif
       message = 'Choose the largest simulated grain size s(N) in '//trim(abbrev)
       call prompt(trim(message),sN,s1,1000.)
       message = 'Choose the cell width in log10 space dS assuming units of '//trim(abbrev)
       call prompt(trim(message),logds,0.01,100.)
       s1   = s1_func(N,sN,logds)
       smin = smin_func(N,s1,sN)
       smax = smax_func(N,smin,s1,sN)
    case(4)
       if (ichoosedist < 0) then
          call prompt('Choose the number of discrete dust phases (N)',N,0)
       endif
       message = 'Enter smin for the distribution in '//trim(abbrev)
       call prompt(trim(message),smin,0.)
       message = 'Enter smax for the distribution in '//trim(abbrev)
       call prompt(trim(message),smax,smin)
       logds = logds_func(N,smin,smax)
    end select

    print*,''
    write(*,'(A)') '------------------------------------------'
    write(*,'(A)') '   Summary of the selected grain sizes:   '
    write(*,'(A)') '------------------------------------------'
    if (ichoosedist <= 0 .and. N /= ndusttypes) then
       write(*,'(A,I2.0)') 'WARNING! If you like these values, you will '//new_line('A')// &
                            '         have to recompile with ndusttypes = ',N
       write(*,'(A,I2.0,A/)') '         (currently ndusttypes = ',ndusttypes,')'
    endif
    write(*,'(A,E21.14E2,1X,A )') 'smin = ',smin,trim(abbrev)
    write(*,'(A,E21.14E2,1X,A/)') 'smax = ',smax,trim(abbrev)

    if (allocated(loggrid)) deallocate(loggrid)
    allocate(loggrid(N+1))
    if (allocated(grid)) deallocate(grid)
    allocate(grid(N+1))
    if (allocated(s)) deallocate(s)
    allocate(s(N))

    do i = 1,N+1
       loggrid(i) = log10(smin) + (i-1)*logds
    enddo

    grid = 10.**loggrid

    s = 0.
    do i = 1,N
       s(i) = sqrt(grid(i)*grid(i+1))
       write(*,'(A,I0.2,A,E21.14E2,1X,A)') 's',i,' = ',s(i),trim(abbrev)
    enddo
    write(*,'(A)') '------------------------------------------'
    print*,''

    if (N == ndusttypes) then
       call prompt('Would you like to try different values?',itryagain)

       if (.not. itryagain) then
          smincgs = units_to_cgs*smin
          smaxcgs = units_to_cgs*smax
          iloop = .false.
       endif
    endif
 enddo

end subroutine grainsize_dist_assist


end module set_dust
