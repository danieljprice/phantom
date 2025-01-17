!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module growth
!
! Contains routine for dust growth and fragmentation
!
! :References:
!  Stepinski & Valageas (1997), A&A 319, 1007
!  Kobayashi & Tanaka (2010), Icarus 206, 735
!  Rozner, Grishin & Perets (2020), MNRAS 496, 4827
!
! :Owner: Arnaud Vericel
!
! :Runtime parameters:
!   - Tsnow         : *snow line condensation temperature in K*
!   - bin_per_dex   : *(mcfost) number of bins of sizes per dex*
!   - cohacc        : *strength of the cohesive acceleration in g/s^2*
!   - dsize         : *size of ejected grain during erosion in cm*
!   - flyby         : *use primary for keplerian freq. calculation*
!   - force_smax    : *(mcfost) set manually maximum size for binning*
!   - grainsizemin  : *minimum allowed grain size in cm*
!   - ieros         : *erosion of dust (0=off,1=on)*
!   - ifrag         : *fragmentation of dust (0=off,1=on,2=Kobayashi)*
!   - isnow         : *snow line (0=off,1=position based,2=temperature based)*
!   - rsnow         : *snow line position in AU*
!   - size_max_user : *(mcfost) maximum size for binning in cm*
!   - tsmincgs      : *minimum allowed stopping time*
!   - vfrag         : *uniform fragmentation threshold in m/s*
!   - vfragin       : *inward fragmentation threshold in m/s*
!   - vfragout      : *inward fragmentation threshold in m/s*
!
! :Dependencies: checkconserved, dim, dust, eos, infile_utils, io, options,
!   part, physcon, table_utils, units, viscosity
!
 use units,        only:udist,umass,utime,unit_density,unit_velocity
 use physcon,      only:au,Ro
 use part,         only:xyzmh_ptmass,nptmass,this_is_a_flyby
 use options,      only:use_mcfost
 implicit none

 !--Default values for the growth and fragmentation of dust in the input file
 integer, public        :: ifrag        = 1
 integer, public        :: isnow        = 0
 integer, public        :: ieros        = 0

 real, public           :: gsizemincgs  = 5.e-3
 real, public           :: tsmincgs     = 1.e5
 real, public           :: rsnow        = 100.
 real, public           :: Tsnow        = 150.
 real, public           :: vfragSI      = 15.
 real, public           :: vfraginSI    = 5.
 real, public           :: vfragoutSI   = 15.
 real, public           :: cohacccgs    = 100
 real, public           :: dsizecgs     = 1.0e-3

 real, public           :: vfrag
 real, public           :: vref
 real, public           :: vfragin
 real, public           :: vfragout
 real, public           :: grainsizemin
 real, public           :: tsmin
 real, public           :: cohacc
 real, public           :: dsize

 logical, public        :: f_smax    = .false.
 real,    public        :: size_max  = 0.2 !- cm
 integer, public        :: b_per_dex = 5

 public                 :: get_growth_rate,get_vrelonvfrag,check_dustprop
 public                 :: write_options_growth,read_options_growth,print_growthinfo,init_growth
 public                 :: vrelative,read_growth_setup_options,write_growth_setup_options
 public                 :: comp_snow_line,bin_to_multi,convert_to_twofluid,get_size

contains

!------------------------------------------------
!+
!  Initialise variables for computing growth rate
!+
!------------------------------------------------
subroutine init_growth(ierr)
 use io,        only:error
 use physcon,   only:fourpi
 use viscosity, only:irealvisc,shearparam
 use dust,      only:grainsizecgs
 use options,   only:use_porosity
 integer, intent(out) :: ierr

 ierr = 0

 !--initialise variables in code units
 vref           = 100 / unit_velocity
 vfrag          = vfragSI * 100 / unit_velocity
 vfragin        = vfraginSI * 100 / unit_velocity
 vfragout       = vfragoutSI * 100 / unit_velocity
 rsnow          = rsnow * au / udist
 grainsizemin   = gsizemincgs / udist
 cohacc         = cohacccgs * utime * utime / umass
 dsize          = dsizecgs / udist
 tsmin          = tsmincgs / utime

 if (ifrag > 0) then
    if (grainsizemin < 0.) then
       call error('init_growth','grainsizemin < 0',var='grainsizemin',val=grainsizemin)
       ierr = 1
    endif
    if (gsizemincgs > grainsizecgs .and. .not. use_porosity) then
       call error('init_growth','grainsizemin > grainsize',var='grainsizemin',val=grainsizemin)
       ierr = 1
    endif
    select case(isnow)
    case(0) !-- uniform vfrag
       if (vfrag <= 0.) then
          call error('init_growth','vfrag <= 0',var='vfrag',val=vfrag)
          ierr = 2
       endif
    case(1) !--position based snow line
       if (rsnow <= 0.) then
          call error('init_growth','rsnow <= 0',var='rsnow',val=rsnow)
          ierr = 2
       endif
    case(2) !-- temperature based snow line
       if (Tsnow <= 0.) then
          call error('init_growth','Tsnow <= 0',var='Tsnow',val=Tsnow)
          ierr = 2
       endif
    case default
       ierr = 0
    end select
 endif

 if (isnow > 0) then
    if (vfragin <= 0) then
       call error('init_growth','vfragin <= 0',var='vfragin',val=vfragin)
       ierr = 3
    endif
    if (vfragout <= 0.) then
       call error('init_growth','vfragout <= 0',var='vfragout',val=vfragout)
       ierr = 3
    endif
 endif

 if (ifrag > -1) then
    if (irealvisc == 1) then
       call error('init_growth','shearparam should be used for growth when irealvisc /= 1',var='shearparam',val=shearparam)
       ierr = 4
    endif
 endif

 if (ieros == 1) then
    if (cohacc < 0) then
       call error('init_growth','cohacc < 0',var='cohacc',val=cohacc)
       ierr = 5
    endif
    if (dsize < 0) then
       call error('init_growth','dsize < 0',var='dsize',val=dsize)
       ierr = 5
    endif
 endif

end subroutine init_growth

!----------------------------------------------------------
!+
!  print information about growth and fragmentation of dust
!+
!----------------------------------------------------------
subroutine print_growthinfo(iprint)
 use viscosity, only:shearparam

 integer, intent(in) :: iprint

 if (ifrag == 0) write(iprint,"(a)")    ' Using pure growth model where dm/dt = + 4pi*rhod*s**2*vrel*dt    '
 if (ifrag == 1) write(iprint,"(a)")    ' Using growth/frag where dm/dt = (+ or -) 4pi*rhod*s**2*vrel*dt   '
 if (ifrag == 2) write(iprint,"(a)")    ' Using growth/frag with Kobayashi fragmentation model '
 if (ifrag > -1) write(iprint,"((a,1pg10.3))")' Computing Vrel with alphaSS = ',shearparam
 if (ifrag > 0) then
    write(iprint,"(2(a,1pg10.3),a)")' grainsizemin = ',gsizemincgs,' cm = ',grainsizemin,' (code units)'
    if (isnow == 1) then
       write(iprint,"(a)")              ' ===> Using position based snow line <=== '
       write(iprint,"(2(a,1pg10.3),a)") ' rsnow = ',rsnow*udist/au,'    AU = ',rsnow, ' (code units)'
    endif
    if (isnow == 2) then
       write(iprint,"(a)")              ' ===> Using temperature based snow line <=== '
       write(iprint,"(2(a,1pg10.3),a)") ' Tsnow = ',Tsnow,' K = ',Tsnow,' (code units)'
    endif
    if (isnow == 0) then
       write(iprint,"(2(a,1pg10.3),a)") ' vfrag = ',vfragSI,' m/s = ',vfrag ,' (code units)'
    else
       write(iprint,"(2(a,1pg10.3),a)") ' vfragin = ',vfraginSI,' m/s = ',vfragin,' (code units)'
       write(iprint,"(2(a,1pg10.3),a)") ' vfragin = ',vfragoutSI,' m/s = ',vfragout,' (code units)'
    endif
 endif
 if (ieros == 1) then
    write(iprint,"(a)")    ' Using aeolian-erosion model where ds = -fourpi*rhos*rhog*s*(deltav**3)*(dsize**2)/(3*cohacc)*dt    '
    write(iprint,"(2(a,1pg10.3),a)")' dsize = ',dsizecgs,' cm = ',dsize,' (code units)'
 endif

end subroutine print_growthinfo

!-----------------------------------------------------------------------
!+
!  Main routine that returns dm/dt and calculate Vrel/Vfrag.
!  This growth model is currently only available for the
!  two-fluid dust method.
!+
!-----------------------------------------------------------------------
subroutine get_growth_rate(npart,xyzh,vxyzu,dustgasprop,VrelVf,dustprop,filfac,dmdt)
 use part,            only:rhoh,idust,igas,iamtype,iphase,isdead_or_accreted,&
                           massoftype,Omega_k,dustfrac,tstop,deltav
 use options,         only:use_dustfrac,use_porosity
 use physcon,         only:fourpi
 use eos,             only:ieos,get_spsound
 real, intent(in)     :: dustprop(:,:)
 real, intent(inout)  :: dustgasprop(:,:)
 real, intent(in)     :: xyzh(:,:)
 real, intent(in)     :: filfac(:)
 real, intent(inout)  :: VrelVf(:),vxyzu(:,:)
 real, intent(out)    :: dmdt(:)
 integer, intent(in)  :: npart
 !
 real                 :: rhog,rhod,vrel,rho,sdust
 integer              :: i,iam

 vrel = 0.
 rhod = 0.
 rho  = 0.

 !--get dm/dt over all particles

 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,ieos,massoftype,use_dustfrac,dustfrac,use_porosity) &
 !$omp shared(ifrag,ieros,utime,umass,dsize,cohacc) &
 !$omp shared(xyzh,vxyzu,dustprop,dustgasprop,dmdt,filfac,VrelVf,tstop,deltav) &
 !$omp private(i,iam,rho,rhog,rhod,vrel,sdust)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       iam = iamtype(iphase(i))

       if (iam == idust .or. (iam == igas .and. use_dustfrac)) then

          if (use_dustfrac .and. iam == igas) then
             !- no need for interpolations
             rho              = rhoh(xyzh(4,i),massoftype(igas))
             rhog             = rho*(1-dustfrac(1,i))
             rhod             = rho*dustfrac(1,i)
             dustgasprop(1,i) = get_spsound(ieos,xyzh(:,i),rhog,vxyzu(:,i))
             dustgasprop(2,i) = rhog
             dustgasprop(3,i) = tstop(1,i) * Omega_k(i)
             dustgasprop(4,i) = sqrt(deltav(1,1,i)**2 + deltav(2,1,i)**2 + deltav(3,1,i)**2)
          else
             rhod = rhoh(xyzh(4,i),massoftype(idust))
          endif

          !--dust size from mass and filling factor
          if (use_porosity) then
             sdust = get_size(dustprop(1,i),dustprop(2,i),filfac(i))
          else
             sdust = get_size(dustprop(1,i),dustprop(2,i))
          endif

          call get_vrelonvfrag(xyzh(:,i),vxyzu(:,i),vrel,VrelVf(i),dustgasprop(:,i))

          !
          !--dustprop(1) = mass, dustprop(2) = intrinsic density,
          !
          !--if statements to compute dm/dt
          !
          if (ifrag == -1) dmdt(i) = 0.
          if ((VrelVf(i) < 1. .or. ifrag == 0) .and. ifrag /= -1) then ! vrel/vfrag < 1 or pure growth --> growth
             dmdt(i) = fourpi*sdust**2*rhod*vrel
          elseif (VrelVf(i) >= 1. .and. ifrag > 0) then ! vrel/vfrag > 1 --> fragmentation
             select case(ifrag)
             case(1)
                dmdt(i) = -fourpi*sdust**2*rhod*vrel ! Symmetrical of Stepinski & Valageas
             case(2)
                dmdt(i) = -fourpi*sdust**2*rhod*vrel*(VrelVf(i)**2)/(1+VrelVf(i)**2) ! Kobayashi model
             end select
          endif
          if (ieros == 1) then  !sqrt(0.0123)=0.110905    !1.65 -> surface energy in cgs
             ! Erosion model of Rozner, Grishin & Perets (2020)
             if (dustgasprop(4,i) >= 0.110905*sqrt(1.65*utime*utime/umass/dustprop(2,i)/dsize)) then
                dmdt(i) = dmdt(i) - fourpi*sdust*dustprop(2,i)*dustgasprop(2,i)*(dustgasprop(4,i)**3)*(dsize**2)/(3.*cohacc)
             endif
          endif
       endif
    else
       dmdt(i) = 0.
    endif
 enddo
 !$omp end parallel do


end subroutine get_growth_rate

!-----------------------------------------------------------------------
!+
!  Compute the local ratio vrel/vfrag and vrel
!+
!-----------------------------------------------------------------------
subroutine get_vrelonvfrag(xyzh,vxyzu,vrel,VrelVf,dustgasprop)
 use viscosity,       only:shearparam
 use physcon,         only:Ro,roottwo
 real, intent(in)     :: xyzh(:)
 real, intent(in)     :: dustgasprop(:)
 real, intent(inout)  :: vrel,vxyzu(:)
 real, intent(out)    :: VrelVf
 real                 :: Vt
 integer              :: izone

 !--compute turbulent velocity
 Vt   = sqrt(roottwo*Ro*shearparam)*dustgasprop(1)
 !--compute vrel
 vrel = vrelative(dustgasprop,Vt)
 !
 !--If statements to compute local ratio vrel/vfrag
 !
 VrelVf = 0. ! default value
 if (ifrag == 0) then
    if (vref > 0.) VrelVf = vrel/vref ! for pure growth, vrel/vfrag gives vrel in m/s
 elseif (ifrag > 0) then
    call comp_snow_line(xyzh,vxyzu,dustgasprop(2),izone)
    select case(izone)
    case(2)
       if (vfragout > 0.) VrelVf = vrel/vfragout
    case(1)
       if (vfragin > 0.) VrelVf = vrel/vfragin
    case default
       if (vfrag > 0.) VrelVf = vrel/vfrag
    end select
 endif

end subroutine get_vrelonvfrag

!----------------------------------------------------------------------------
!+
! Get the location of a given particle (dust or gas or mixture) with respect
! to a snow line. (You can add more with this structure)
!+
!----------------------------------------------------------------------------
subroutine comp_snow_line(xyzh,vxyzu,rhogas,izone)
 use eos,           only:ieos,get_temperature
 integer, intent(out) :: izone
 real, intent(inout)  :: vxyzu(:)
 real, intent(in)     :: xyzh(:)
 real, intent(in)     :: rhogas
 real                 :: r,Tgas

 select case(isnow)
 case(0)
    izone = 0
 case(1)
    r = sqrt(xyzh(1)**2 + xyzh(2)**2 + xyzh(3)**2)
    if (r<=rsnow) izone = 1
    if (r>rsnow) izone = 2
 case(2)
    Tgas = get_temperature(ieos,xyzh,rhogas,vxyzu)
    if (Tgas >= Tsnow) izone = 1
    if (Tgas < Tsnow) izone = 2
 case default
    izone = 0
 end select

end subroutine comp_snow_line

!-----------------------------------------------------------------------
!+
!  Write growth options in the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_growth(iunit)
 use infile_utils, only:write_inopt
 use options,      only:use_porosity
 integer, intent(in)        :: iunit

 write(iunit,"(/,a)") '# options controlling growth'
 if (nptmass > 1) call write_inopt(this_is_a_flyby,'flyby','use primary for keplerian freq. calculation',iunit)
 call write_inopt(ifrag,'ifrag','dust fragmentation (0=off,1=on,2=Kobayashi)',iunit)
 call write_inopt(ieros,'ieros','erosion of dust (0=off,1=on)',iunit)
 if (ifrag /= 0) then
    if (use_porosity) then
       call write_inopt(tsmincgs,'tsmincgs','minimum allowed stopping time',iunit)
    else
       call write_inopt(gsizemincgs,'grainsizemin','minimum grain size in cm',iunit)
    endif
    call write_inopt(isnow,'isnow','snow line (0=off,1=position based,2=temperature based)',iunit)
    if (isnow == 1) call write_inopt(rsnow,'rsnow','position of the snow line in AU',iunit)
    if (isnow == 2) call write_inopt(Tsnow,'Tsnow','snow line condensation temperature in K',iunit)
    if (isnow == 0) call write_inopt(vfragSI,'vfrag','uniform fragmentation threshold in m/s',iunit)
    if (isnow > 0) then
       call write_inopt(vfraginSI,'vfragin','inward fragmentation threshold in m/s',iunit)
       call write_inopt(vfragoutSI,'vfragout','outward fragmentation threshold in m/s',iunit)
    endif
 endif
 if (ieros == 1) then
    call write_inopt(cohacccgs,'cohacc','strength of the cohesive acceleration in g/s^2',iunit)
    call write_inopt(dsizecgs,'dsize','size of ejected grain during erosion in cm',iunit)
 endif

 if (use_mcfost) then
    call write_inopt(f_smax,'force_smax','(mcfost) set manually maximum size for binning',iunit)
    call write_inopt(size_max,'size_max_user','(mcfost) maximum size for binning in cm',iunit)
    call write_inopt(b_per_dex,'bin_per_dex','(mcfost) number of bins of sizes per dex',iunit)
 endif

end subroutine write_options_growth

!-----------------------------------------------------------------------
!+
!  Read growth options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_growth(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)        :: name,valstring
 logical,intent(out)                 :: imatch,igotall
 integer,intent(out)                 :: ierr

 integer,save                        :: ngot = 0
 integer                             :: imcf = 0
 integer                             :: goteros = 1
 logical                             :: tmp = .false.

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('ifrag')
    read(valstring,*,iostat=ierr) ifrag
    ngot = ngot + 1
 case('ieros')
    read(valstring,*,iostat=ierr) ieros
    ngot = ngot + 1
 case('grainsizemin')
    read(valstring,*,iostat=ierr) gsizemincgs
    ngot = ngot + 1
 case('tsmincgs')
    read(valstring,*,iostat=ierr) tsmincgs
    ngot = ngot + 1
 case('isnow')
    read(valstring,*,iostat=ierr) isnow
    ngot = ngot + 1
 case('rsnow')
    read(valstring,*,iostat=ierr) rsnow
    ngot = ngot + 1
 case('Tsnow')
    read(valstring,*,iostat=ierr) Tsnow
    ngot = ngot + 1
 case('vfrag')
    read(valstring,*,iostat=ierr) vfragSI
    ngot = ngot + 1
 case('vfragin')
    read(valstring,*,iostat=ierr) vfraginSI
    ngot = ngot + 1
 case('vfragout')
    read(valstring,*,iostat=ierr) vfragoutSI
    ngot = ngot + 1
 case('cohacc')
    read(valstring,*,iostat=ierr) cohacccgs
    ngot = ngot + 1
 case('dsize')
    read(valstring,*,iostat=ierr) dsizecgs
    ngot = ngot + 1
 case('flyby')
    read(valstring,*,iostat=ierr) this_is_a_flyby
    ngot = ngot + 1
    if (nptmass < 2) tmp = .true.
 case('force_smax')
    read(valstring,*,iostat=ierr) f_smax
    ngot = ngot + 1
 case('size_max_user')
    read(valstring,*,iostat=ierr) size_max
    ngot = ngot + 1
 case('bin_per_dex')
    read(valstring,*,iostat=ierr) b_per_dex
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 if (use_mcfost) imcf = 3

 if (ieros == 1) goteros = 3

 if (nptmass > 1 .or. tmp) then
    if ((ifrag <= 0) .and. ngot == 2+imcf+goteros) igotall = .true.
    if (isnow == 0) then
       if (ngot == 5+imcf+goteros) igotall = .true.
    elseif (isnow > 0) then
       if (ngot == 7+imcf+goteros) igotall = .true.
    else
       igotall = .false.
    endif
 else
    if ((ifrag <= 0) .and. ngot == 1+imcf+goteros) igotall = .true.
    if (isnow == 0) then
       if (ngot == 4+imcf+goteros) igotall = .true.
    elseif (isnow > 0) then
       if (ngot == 6+imcf+goteros) igotall = .true.
    else
       igotall = .false.
    endif
 endif

end subroutine read_options_growth

!-----------------------------------------------------------------------
!+
!  Write growth options to the .setup file
!+
!-----------------------------------------------------------------------
subroutine write_growth_setup_options(iunit)
 use infile_utils,      only:write_inopt
 use options,           only:use_porosity
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for growth and fragmentation of dust'

 call write_inopt(ifrag,'ifrag','fragmentation of dust (0=off,1=on,2=Kobayashi)',iunit)
 call write_inopt(ieros,'ieros','erosion of dust (0=off,1=on)',iunit)
 call write_inopt(isnow,'isnow','snow line (0=off,1=position based,2=temperature based)',iunit)
 call write_inopt(rsnow,'rsnow','snow line position in AU',iunit)
 call write_inopt(Tsnow,'Tsnow','snow line condensation temperature in K',iunit)
 call write_inopt(vfragSI,'vfrag','uniform fragmentation threshold in m/s',iunit)
 call write_inopt(vfraginSI,'vfragin','inward fragmentation threshold in m/s',iunit)
 call write_inopt(vfragoutSI,'vfragout','inward fragmentation threshold in m/s',iunit)
 if (use_porosity) then
    call write_inopt(tsmincgs,'tsmincgs','minimum allowed stopping time',iunit)
 else
    call write_inopt(gsizemincgs,'grainsizemin','minimum allowed grain size in cm',iunit)
 endif

end subroutine write_growth_setup_options

!-----------------------------------------------------------------------
!+
!  Read growth options from the .setup file
!+
!-----------------------------------------------------------------------
subroutine read_growth_setup_options(db,nerr)
 use infile_utils,    only:read_inopt,inopts
 use options,         only:use_porosity
 type(inopts), allocatable, intent(inout) :: db(:)
 integer, intent(inout)                   :: nerr

 call read_inopt(ifrag,'ifrag',db,min=-1,max=2,errcount=nerr)
 call read_inopt(ieros,'ieros',db,min=0,max=1,errcount=nerr)
 if (ifrag > 0) then
    call read_inopt(isnow,'isnow',db,min=0,max=2,errcount=nerr)
    if (use_porosity) then
       call read_inopt(tsmincgs,'tsmincgs',db,min=1.e-5,errcount=nerr)
    else
       call read_inopt(gsizemincgs,'grainsizemin',db,min=1.e-5,errcount=nerr)
    endif
    select case(isnow)
    case(0)
       call read_inopt(vfragSI,'vfrag',db,min=0.,errcount=nerr)
    case(1)
       call read_inopt(rsnow,'rsnow',db,min=0.,errcount=nerr)
       call read_inopt(vfraginSI,'vfragin',db,min=0.,errcount=nerr)
       call read_inopt(vfragoutSI,'vfragout',db,min=0.,errcount=nerr)
    case(2)
       call read_inopt(Tsnow,'Tsnow',db,min=0.,errcount=nerr)
       call read_inopt(vfraginSI,'vfragin',db,min=0.,errcount=nerr)
       call read_inopt(vfragoutSI,'vfragout',db,min=0.,errcount=nerr)
    end select
 endif

end subroutine read_growth_setup_options

!-----------------------------------------------------------------------
!+
!  In case of fragmentation, limit masses to a minimum value
!+
!-----------------------------------------------------------------------
subroutine check_dustprop(npart,dustprop,filfac,mprev,filfacprev)
 use part,                 only:iamtype,iphase,idust,igas,dustgasprop,Omega_k
 use options,              only:use_dustfrac,use_porosity
 real,intent(inout)        :: dustprop(:,:)
 integer,intent(in)        :: npart
 real, intent(in)          :: filfac(:),mprev(:),filfacprev(:)
 integer                   :: i,iam
 real                      :: tsnew,sdustprev,sdustmin,sdust

 !$omp parallel do default(none) &
 !$omp shared(iphase,dustgasprop,use_dustfrac,use_porosity) &
 !$omp shared(npart,ifrag,dustprop,filfac,mprev,filfacprev) &
 !$omp shared(tsmin,grainsizemin) &
 !$omp private(i,iam,tsnew,sdustprev,sdustmin,sdust)
 do i=1,npart
    iam = iamtype(iphase(i))
    if ((iam == idust .or. (use_dustfrac .and. iam == igas))  .and. ifrag > 0 .and. dustprop(1,i) <= mprev(i)) then
       if (use_porosity) then
          sdustprev = get_size(mprev(i),dustprop(2,i),filfacprev(i))
          sdust = get_size(dustprop(1,i),dustprop(2,i),filfac(i))
          tsnew = dustgasprop(3,i)*sdustprev*filfacprev(i)/sdust/filfac(i)/Omega_k(i)
          if (tsnew < tsmin) then
             sdustmin = tsmin*sdustprev*filfacprev(i)*Omega_k(i)/filfac(i)/dustgasprop(3,i)
             dustprop(1,i) = dustprop(1,i) * (sdustmin/sdust)**3.
          endif
       else
          sdust = get_size(dustprop(1,i),dustprop(2,i))
          if (sdust < grainsizemin) then
             dustprop(1,i) = dustprop(1,i) * (grainsizemin/sdust)**3.   ! fragmentation at constant density and filling factor
          endif
       endif
    endif
 enddo
 !$omp end parallel do

end subroutine check_dustprop

!-----------------------------------------------------------------------
!+
!  Set dustprop (used by moddump)
!+
!-----------------------------------------------------------------------
subroutine set_dustprop(npart,xyzh,sizedistrib,pwl_sizedistrib,R_ref,H_R_ref,q_index)
 use dust,    only:grainsizecgs,graindenscgs
 use part,    only:iamtype,iphase,idust,igas,dustprop,filfac,probastick
 use physcon, only:fourpi
 use options, only:use_dustfrac
 integer, intent(in)           :: npart
 real, intent(in)              :: xyzh(:,:)
 integer                       :: i,iam
 real                          :: r,h
 logical, optional, intent(in) :: sizedistrib
 real, optional, intent(in)    :: pwl_sizedistrib,R_ref,H_R_ref,q_index

 do i=1,npart
    iam = iamtype(iphase(i))
    if (iam == idust .or. (iam == igas .and. use_dustfrac)) then
       dustprop(2,i) = graindenscgs / unit_density
       if (sizedistrib) then
          r = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
          h = H_R_ref * R_ref * au / udist * (r * udist / au / R_ref)**(1.5-q_index)
          dustprop(1,i) = grainsizecgs/udist * (r * udist / au / R_ref)**pwl_sizedistrib * exp(-0.5*xyzh(3,i)**2/h**2)
          if (dustprop(1,i) < 2.e-5/udist) then
             dustprop(1,i) = 2.e-5/udist
          endif
          dustprop(1,i) = fourpi/3. * dustprop(2,i) * (dustprop(1,i))**3
       else
          dustprop(1,i) = fourpi/3. * dustprop(2,i) * (grainsizecgs / udist)**3
       endif
    else
       dustprop(:,i) = 0.
    endif
    filfac(i) = 0.
    probastick(i) = 1.
 enddo

end subroutine set_dustprop

!-----------------------------------------------------------------------
!+
!  Bin sizes to fake multi large grains (used by moddump and live mcfost)
!+
!-----------------------------------------------------------------------
subroutine bin_to_multi(bins_per_dex,force_smax,smax_user,verbose)
 use part,           only:npart,npartoftype,massoftype,ndusttypes,&
                          ndustlarge,grainsize,dustprop,&
                          iamtype,iphase,set_particle_type,idust,filfac
 use options,        only:use_porosity
 use units,          only:udist
 use table_utils,    only:logspace
 use io,             only:fatal
 use checkconserved, only:mdust_in
 integer, intent(in)    :: bins_per_dex
 real,    intent(in)    :: smax_user
 logical, intent(inout) :: force_smax
 logical, intent(in)    :: verbose
 real                   :: smaxtmp,smintmp,smax,smin,tolm
 real                   :: mdustold,mdustnew,code_to_mum
 logical                :: init
 integer                :: nbinsize,nbinsizemax,i,j,itype,ndustold,ndustnew,npartmin,imerge,iu
 integer                :: nbinfilfacmax,ndustsizetypes
 real, allocatable, dimension(: )  :: grid
 real, allocatable, dimension(:,:) :: dustpropmcfost  !dustpropmcfost(1=size,2=filfac)
 character(len=20)                 :: outfile = "bin_distrib.dat"

 !- initialise
 code_to_mum   = udist*1.e4
 tolm          = 1.e-5
 smaxtmp       = 0.
 smintmp       = 1.e26
 ndustold      = 0
 ndustnew      = 0
 mdustold      = 0.
 mdustnew      = 0.
 nbinsizemax   = 25
 nbinfilfacmax = 10
 npartmin      = 50 !- limit to find neighbours
 init          = .false.
 allocate(dustpropmcfost(2,npart))
 !graindens     = maxval(dustprop(2,:))

 !- loop over particles, find min and max on non-accreted dust particles
 do i = 1,npart
    itype = iamtype(iphase(i))
    if (itype==idust) then
       if (use_porosity) then
          dustpropmcfost(1,i) = get_size(dustprop(1,i),dustprop(2,i),filfac(i))
          dustpropmcfost(2,i) = filfac(i)
       else
          dustpropmcfost(1,i) = get_size(dustprop(1,i),dustprop(2,i))
          dustpropmcfost(2,i) = 1
       endif
       if (dustpropmcfost(1,i) < smintmp) smintmp = dustpropmcfost(1,i)
       if (dustpropmcfost(1,i) > smaxtmp) smaxtmp = dustpropmcfost(1,i)
       !if (dustpropmcfost(2,i) < fmintmp) fmintmp = dustpropmcfost(2,i)
       !if (dustpropmcfost(2,i) > fmaxtmp) fmaxtmp = dustpropmcfost(2,i)
    endif
 enddo

 !- overrule force_smax if particles are small, avoid empty bins
 if ((maxval(dustpropmcfost(1,1:npart))*udist < smax_user) .and. force_smax) then
    force_smax = .false.
    write(*,*) "Overruled force_smax from T to F"
 endif

 !- force smax if needed, check for flat size distribution
 if (force_smax .and. (smintmp /= smaxtmp)) then
    smax = smax_user/udist
 elseif (smintmp /= smaxtmp) then
    smax = smaxtmp
 else
    init = .true.
    write(*,*) "Detected initial condition, restraining nbinsize = 1"
 endif

 if (.not. init) then
    smin = smintmp

    !- set ndusttypes based on desired log size spacing
    nbinsize   = int((log10(smax)-log10(smin))*bins_per_dex + 1.)
    !nbinfilfac = int((log10(smax)-log10(smin))*bins_per_dex + 1.)
    ndustsizetypes = min(nbinsize, nbinsizemax) !- prevent memory allocation errors
    !ndustfilfactypes = min(nbinfilfac,nbinfilfacmax)
    ndustlarge = ndustsizetypes!*ndustfilfactypes !- this is written to the header

    !- allocate memory for a grid of ndusttypes+1 elements
    allocate(grid(ndustsizetypes+1))!,ndustfilfactypes))

    !- bin sizes in ndusttypes bins
    write(*,"(a,f10.1,a,f10.1,a,i3,a)") "Binning sizes between ",smin*code_to_mum, " (µm) and ",&
                                     smax*code_to_mum," (µm) in ",ndusttypes, " bins"

    call logspace(grid(1:ndustsizetypes+1),smin,smax) !- bad for live mcfost, need to compile it before growth.F90
    !call logspace(grid(2,1:ndustfilfactypes),fmin,fmax)

    !- find representative size for each bin
    do i = 1,ndustsizetypes
       grainsize(i) = sqrt(grid(i)*grid(i+1))
    enddo

    !- transfer particles from bin to bin depending on their size
    ndustold = npartoftype(idust)
    mdustold = massoftype(idust)*npartoftype(idust) !- initial total mass
    do i=1,npart
       itype = iamtype(iphase(i))
       if (itype==idust) then
          !- figure out which bin
          do j=1,ndusttypes
             if ((dustpropmcfost(1,i) >= grid(j)) .and. (dustpropmcfost(1,i) < grid(j+1))) then
                if (j > 1) then
                   npartoftype(idust+j-1) = npartoftype(idust+j-1) + 1
                   npartoftype(idust)     = npartoftype(idust) - 1
                   call set_particle_type(i,idust+j-1)
                endif
             endif
             !- if smax has been forced, put larger grains inside last bin
             if ((j==ndusttypes) .and. force_smax .and. (dustpropmcfost(1,i) >= grid(j+1))) then
                npartoftype(idust+j-1) = npartoftype(idust+j-1) + 1
                npartoftype(idust)     = npartoftype(idust) - 1
                call set_particle_type(i,idust+j-1)
             endif
          enddo
       endif
    enddo

    !- check npart inside each bin, merge bins if necessary (needed for mcfost live)
    do while (any(npartoftype(idust:idust+ndusttypes-1) < npartmin))
       call merge_bins(npart,grid,npartmin)
       imerge = imerge + 1
       if (imerge > 50) call fatal('bin merging','merging number of iterations exceeded limit',var="imerge",ival=imerge)
    enddo

    !- set massoftype for each bin and print info
    if (verbose) open(newunit=iu,file=outfile,status="replace")
    write(*,"(a3,a1,a10,a1,a10,a1,a10,a5,a6)") "Bin #","|","s_min","|","s","|","s_max","|==>","npart"

    do itype=idust,idust+ndusttypes-1
       write(*,"(i3,a1,f10.1,a1,f10.1,a1,f10.1,a5,i6)") itype-idust+1,"|",grid(itype-idust+1)*code_to_mum,"|", &
                                                                     grainsize(itype-idust+1)*code_to_mum, &
                                                                     "|",grid(itype-idust+2)*code_to_mum,"|==>",npartoftype(itype)

       if (itype > idust) massoftype(itype) = massoftype(idust)
       mdust_in(itype-idust+1) = massoftype(itype)*npartoftype(itype)
       mdustnew        = mdustnew + mdust_in(itype-idust+1)
       ndustnew        = ndustnew + npartoftype(itype)

       if (verbose) write(iu,*) itype-idust+1,grid(itype-idust+1)*code_to_mum,grainsize(itype-idust+1)*code_to_mum,&
                     grid(itype-idust+2)*code_to_mum,npartoftype(itype)
    enddo

    if (verbose) close(unit=iu)

    !- sanity check for total number of dust particles
    if (sum(npartoftype(idust:)) /= ndustold) then
       write(*,*) 'ERROR! npartoftype not conserved'
       write(*,*) sum(npartoftype(idust:)), " <-- new vs. old --> ",ndustold
    endif

    !- sanity check for total dust mass
    if (abs(mdustold-mdustnew)/mdustold > tolm) then
       write(*,*) 'ERROR! total dust mass not conserved'
       write(*,*) mdustnew, " <-- new vs. old --> ",mdustold
    endif
 else !- init
    grainsize(ndusttypes) = smaxtmp !- only 1 bin, all particles have same size
 endif

 ! clean up
 if (allocated(dustpropmcfost)) deallocate(dustpropmcfost)
 if (allocated(grid)) deallocate(grid)

end subroutine bin_to_multi

!-----------------------------------------------------------------------
!+
!  Merge bins if too many particles are in them
!+
!-----------------------------------------------------------------------
subroutine merge_bins(npart,grid,npartmin)
 use part,           only:ndusttypes,ndustlarge,set_particle_type,&
                         npartoftype,massoftype,idust,iphase,iamtype,&
                         grainsize,graindens
 use checkconserved, only:mdust_in
 integer, intent(in) :: npart,npartmin
 real, intent(inout) :: grid(:)
 integer             :: i,iculprit,itype,idusttype,nculprit,nother,iother
 logical             :: backward = .true.

!- initialise
 i         = 0
 iculprit  = 0
 itype     = 0
 idusttype = 0
 nculprit  = 0
 nother    = 0
 iother    = 0

!- scan npartoftype, get index of most upper bin
 do i=ndusttypes+idust-1,idust,-1
    if (npartoftype(i) < npartmin) then
       iculprit  = i
       idusttype = iculprit - idust + 1
       nculprit  = npartoftype(iculprit)
       if (iculprit == idust) then
          iother = iculprit + 1
          write(*,*) "Merging bin number ",idusttype," forward"
          backward = .false.
       else
          iother = iculprit - 1
          write(*,*) "Merging bin number ",idusttype," backward"
       endif
       nother = npartoftype(iother)
       exit
    endif
 enddo

!- transfer particles of that type to bin n-1, set particle type to n-1
 do i=1,npart
    itype = iamtype(iphase(i))
    if (backward) then
       if (itype == iculprit) then
          npartoftype(iculprit) = npartoftype(iculprit) - 1
          npartoftype(iother)   = npartoftype(iother) + 1
          call set_particle_type(i,iother)
       endif
    else !- translate everyone to the bin to their left, except culprit bin
       if (itype /= iculprit) then
          npartoftype(itype)   = npartoftype(itype) - 1
          npartoftype(itype-1) = npartoftype(itype-1) + 1
          call set_particle_type(i,itype-1)
       endif
    endif
 enddo

!- recompute grainsize and grid borders
 if (backward) then
    massoftype(iculprit) = 0.
    mdust_in(iculprit)   = 0.
    graindens(idusttype) = 0.

    !- recompute size with weigthed sum
    grainsize(idusttype-1) = (grainsize(idusttype-1)*nother + grainsize(idusttype)*nculprit) / (nother + nculprit)
    grid(idusttype)        = grid(idusttype+1)
 else
    do i=1,ndusttypes
       if (i==1) then
          grainsize(i) = (grainsize(i)*npartoftype(i+idust-1) + grainsize(i+1)*npartoftype(i+idust)) &
                        / (npartoftype(i+idust-1)+npartoftype(i+idust))
       else
          grainsize(i) = grainsize(i+1)
          grid(i)      = grid(i+1)
       endif
    enddo
    massoftype(ndusttypes+idust-1) = 0.
    mdust_in(ndusttypes+idust-1)   = 0.
    graindens(ndusttypes+idust-1)  = 0.
    grid(ndusttypes+1)             = 0.
 endif

!- reduce ndusttypes
 ndusttypes           = ndusttypes - 1
 ndustlarge           = ndusttypes

end subroutine merge_bins

!-----------------------------------------------------------------------
!+
!  Convert a one-fluid dustgrowth sim into a two-fluid one (used by moddump_dustadd)
!+
!-----------------------------------------------------------------------
subroutine convert_to_twofluid(npart,xyzh,vxyzu,massoftype,npartoftype,np_ratio,dust_to_gas)
 use part,            only: dustprop,dustgasprop,ndustlarge,ndustsmall,igas,idust,VrelVf,&
                             dustfrac,iamtype,iphase,deltav,set_particle_type,filfac
 use options,         only: use_dustfrac,use_porosity
 use dim,             only: update_max_sizes
 integer, intent(inout)  :: npart,npartoftype(:)
 real, intent(inout)     :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 integer, intent(in)     :: np_ratio
 real, intent(in)        :: dust_to_gas
 integer                 :: np_gas,np_dust,j,ipart,iloc,iam

 !- add number of dust particles
 np_gas = npartoftype(igas)
 ndustlarge = 1
 np_dust = np_gas/np_ratio
 npart = np_gas + np_dust

 !- update memore allocation
 call update_max_sizes(npart)

 !- set dust quantities
 do j=1,np_dust
    ipart = np_gas + j
    iloc  = np_ratio*j

    xyzh(1,ipart) = xyzh(1,iloc)
    xyzh(2,ipart) = xyzh(2,iloc)
    xyzh(3,ipart) = xyzh(3,iloc)
    xyzh(4,ipart) = xyzh(4,iloc) * (np_ratio*dust_to_gas/dustfrac(1,iloc))**(1./3.) !- smoothing lengths

    !- dust velocities out of the barycentric frame
    vxyzu(1,ipart) = vxyzu(1,iloc) + (1 - dustfrac(1,iloc)) * deltav(1,1,iloc)
    vxyzu(2,ipart) = vxyzu(2,iloc) + (1 - dustfrac(1,iloc)) * deltav(2,1,iloc)
    vxyzu(3,ipart) = vxyzu(3,iloc) + (1 - dustfrac(1,iloc)) * deltav(3,1,iloc)

    dustprop(1,ipart)    = dustprop(1,iloc)
    dustprop(2,ipart)    = dustprop(2,iloc)
    dustgasprop(1,ipart) = dustgasprop(1,iloc)
    dustgasprop(2,ipart) = dustgasprop(2,iloc)
    dustgasprop(3,ipart) = dustgasprop(3,iloc)
    dustgasprop(4,ipart) = dustgasprop(4,iloc)
    VrelVf(ipart)        = VrelVf(iloc)
    if (use_porosity) then
       filfac(ipart)     = filfac(iloc)
    endif

    call set_particle_type(ipart,idust)
 enddo

 !- gas quantities out of barycentric frame
 do j=1,npart
    iam = iamtype(iphase(j))
    if (iam == igas) then
       !- smoothing lenghts
       xyzh(4,j) =  xyzh(4,j) * (1-dustfrac(1,j))**(-1./3.)
       !- velocities
       vxyzu(1,j) = vxyzu(1,j) - dustfrac(1,j) * deltav(1,1,j)
       vxyzu(2,j) = vxyzu(2,j) - dustfrac(1,j) * deltav(2,1,j)
       vxyzu(3,j) = vxyzu(3,j) - dustfrac(1,j) * deltav(3,1,j)

       !- reset dust quantities of the mixture
       dustprop(1,j)    = 0.
       dustprop(2,j)    = 0.
       dustgasprop(1,j) = 0.
       dustgasprop(2,j) = 0.
       dustgasprop(3,j) = 0.
       dustgasprop(4,j) = 0.
       VrelVf(j)        = 0.
       if (use_porosity) then
          filfac(j)     = 0.
       endif

    endif
 enddo

 !- unset onefluid, add mass of type dust
 massoftype(idust)  = massoftype(igas)*dust_to_gas*np_ratio
 npartoftype(idust) = np_dust
 ndustsmall         = 0
 use_dustfrac       = .false.

end subroutine convert_to_twofluid

!-----------------------------------------------------------------------
!+
!  Compute the relative velocity following Stepinski & Valageas (1997)
!+
!-----------------------------------------------------------------------
real function vrelative(dustgasprop,Vt)
 use physcon,     only:roottwo
 real, intent(in) :: dustgasprop(:),Vt
 real             :: Sc

 vrelative = 0.
 Sc        = 0.

 !--compute Schmidt number Sc
 if (Vt > 0. .and. dustgasprop(4) > 0.) then
    Sc = (1. + dustgasprop(3)) * sqrt(1. + (dustgasprop(4)/Vt)**2)
 else
    Sc = 1. + dustgasprop(3)
 endif
 if (Sc > 0.) vrelative = roottwo*Vt*sqrt(Sc-1.)/(Sc)

end function vrelative

!--Compute size from mass and filling factor
real function get_size(mass,dens,filfac)
 use physcon,            only:fourpi
 real, intent(in)           :: mass,dens
 real, optional, intent(in) :: filfac
 real                       :: f

 if (present(filfac)) then
    f = filfac
 else
    f = 1.0
 endif

 if (dens > 0. .and. f > 0.) then
    get_size = ( 3.*mass / (fourpi*dens*f) )**(1./3.)
 else
    get_size = 0.
 endif

end function get_size

end module growth
