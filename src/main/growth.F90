!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: growth
!
!  DESCRIPTION:
!  Contains routine for dust growth and fragmentation
!
!  REFERENCES:
!  Stepinski & Valageas (1997)
!  Kobayashi & Tanaka (2009)
!
!  OWNER: Arnaud Vericel
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    Tsnow        -- snow line condensation temperature in K
!    grainsizemin -- minimum allowed grain size in cm
!    ifrag        -- fragmentation of dust (0=off,1=on,2=Kobayashi)
!    isnow        -- snow line (0=off,1=position based,2=temperature based)
!    rsnow        -- snow line position in AU
!    vfrag        -- uniform fragmentation threshold in m/s
!    vfragin      -- inward fragmentation threshold in m/s
!    vfragout     -- inward fragmentation threshold in m/s
!
!  DEPENDENCIES: dust, eos, infile_utils, io, options, part, physcon, units
!+
!--------------------------------------------------------------------------
module growth
 use units,        only:udist,unit_density,unit_velocity
 use physcon,      only:au,Ro
 use part,         only:xyzmh_ptmass
 implicit none

 !--Default values for the growth and fragmentation of dust in the input file
 integer, public        :: ifrag        = 1
 integer, public        :: isnow        = 0
 logical, public        :: iinterpol    = .true.

 real, public           :: gsizemincgs  = 1.e-3
 real, public           :: rsnow        = 100.
 real, public           :: Tsnow        = 20.
 real, public           :: vfragSI      = 15.
 real, public           :: vfraginSI    = 5.
 real, public           :: vfragoutSI   = 15.

 real, public           :: vfrag
 real, public           :: vref
 real, public           :: vfragin
 real, public           :: vfragout
 real, public           :: grainsizemin

 public                 :: get_growth_rate,get_vrelonvfrag,check_dustprop
 public                 :: write_options_growth,read_options_growth,print_growthinfo,init_growth
 public                 :: vrelative,read_growth_setup_options,write_growth_setup_options
 public                 :: comp_snow_line

contains

!------------------------------------------------
!+
!  Initialise variables for computing growth rate
!+
!------------------------------------------------
subroutine init_growth(ierr)
 use io,        only:error
 use part,        only:dustprop,npart
 integer, intent(out) :: ierr

 integer              :: i

 i = 0
 ierr = 0

 !--initialise variables in code units
 vref           = 100 / unit_velocity
 vfrag          = vfragSI * 100 / unit_velocity
 vfragin        = vfraginSI * 100 / unit_velocity
 vfragout       = vfragoutSI * 100 / unit_velocity
 rsnow          = rsnow * au / udist
 grainsizemin   = gsizemincgs / udist

 !-- Check that all the parameters are > 0 when needed
 do i=1,npart
    if (dustprop(1,i) < 0) then
       call error('init_growth','grainsize < 0',var='dustprop',val=dustprop(1,i))
       ierr = 1
    endif
    if (dustprop(2,i) < 0) then
       call error('init_growth','graindens < 0',var='dustprop',val=dustprop(2,i))
       ierr = 1
    endif
    if (dustprop(3,i) < 0) then
       call error('init_growth','vrel < 0',var='dustprop',val=dustprop(3,i))
       ierr = 1
    endif
    if (dustprop(4,i) < 0) then
       call error('init_growth','vrel/vfrag < 0',var='dustprop',val=dustprop(4,i))
       ierr = 1
    endif
 enddo

 if (ifrag > 0) then
    if (grainsizemin < 0) then
       call error('init_growth','grainsizemin < 0',var='grainsizemin',val=grainsizemin)
       ierr = 1
    endif
    select case(isnow)
    case(0) !-- uniform vfrag
       if (vfrag <= 0) then
          call error('init_growth','vfrag <= 0',var='vfrag',val=vfrag)
          ierr = 2
       endif
    case(1) !--position based snow line
       if (rsnow <= 0) then
          call error('init_growth','rsnow <= 0',var='rsnow',val=rsnow)
          ierr = 2
       endif
    case(2) !-- temperature based snow line
       if (Tsnow <= 0) then
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
    if (vfragout <= 0) then
       call error('init_growth','vfragout <= 0',var='vfragout',val=vfragout)
       ierr = 3
    endif
 endif

end subroutine init_growth

!----------------------------------------------------------
!+
!  print information about growth and fragmentation of dust
!+
!----------------------------------------------------------
subroutine print_growthinfo(iprint)
 integer, intent(in) :: iprint

 if (ifrag == 0) write(iprint,"(a)")    ' Using pure growth model where ds = + vrel*rhod/graindens*dt    '
 if (ifrag == 1) write(iprint,"(a)")    ' Using growth/frag where ds = (+ or -) vrel*rhod/graindens*dt   '
 if (ifrag == 2) write(iprint,"(a)")    ' Using growth with Kobayashi fragmentation model '
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

end subroutine print_growthinfo

!-----------------------------------------------------------------------
!+
!  Main routine that make the dust grow and shatter.
!  It is currently available only for the
!  two-fluid dust method.
!+
!-----------------------------------------------------------------------
subroutine get_growth_rate(npart,xyzh,vxyzu,dustprop,dsdt)
 use part,            only:get_pmass,rhoh,idust,iamtype,iphase,St,maxvxyzu
 use eos,             only:get_spsound,ieos
 real, intent(inout)  :: dustprop(:,:),vxyzu(:,:)
 real, intent(in)     :: xyzh(:,:)
 real, intent(out)    :: dsdt(:)
 integer, intent(in)  :: npart
 !
 real                 :: rhod,cs,vrel
 integer              :: i,iam

 !--get ds/dt over all dust particles
 do i=1,npart

    iam = iamtype(iphase(i))

    if (iam==idust) then

       rhod = rhoh(xyzh(4,i),get_pmass(idust,.false.))
       cs   = get_spsound(ieos,xyzh(:,i),rhod,vxyzu(:,i))
       call get_vrelonvfrag(xyzh(:,i),vrel,dustprop(:,i),cs,St(i))
       !
       !--dustprop(1)= size, dustprop(2) = intrinsic density,
       !  dustprop(3) = vrel/vfrag, dustprop(4) = vd - vg
       !
       !--if statements to compute ds/dt
       !
       if (ifrag == -1) dsdt(i) = 0.
       if ((dustprop(3,i) < 1. .or. ifrag == 0) .and. ifrag /= -1) then ! vrel/vfrag < 1 or pure growth --> growth
          dsdt(i) = rhod/dustprop(2,i)*vrel
       elseif (dustprop(3,i) >= 1. .and. ifrag > 0) then ! vrel/vfrag > 1 --> fragmentation
          select case(ifrag)
          case(1)
             dsdt(i) = -rhod/dustprop(2,i)*vrel ! Symmetrical of Stepinski & Valageas
          case(2)
             dsdt(i) = -rhod/dustprop(2,i)*vrel*(dustprop(3,i)**2)/(1+dustprop(3,i)**2) ! Kobayashi model
          end select
       endif
    endif
 enddo
end subroutine get_growth_rate

!-----------------------------------------------------------------------
!+
!  Compute the local ratio vrel/vfrag and vrel
!+
!-----------------------------------------------------------------------
subroutine get_vrelonvfrag(xyzh,vrel,dustprop,cs,St)
 use options,         only:alpha
 use physcon,         only:Ro
 real, intent(in)     :: xyzh(:)
 real, intent(in)     :: cs,St
 real, intent(inout)  :: dustprop(:)
 real, intent(out)    :: vrel
 real                 :: Vt
 integer              :: izone

 !--compute terminal velocity
 Vt = sqrt(sqrt(2)*Ro*alpha)*cs

 !--compute vrel
 vrel = vrelative(St,dustprop(4),Vt)
 !
 !--If statements to compute local ratio vrel/vfrag
 !
 if (ifrag == 0) then
    dustprop(3) = vrel/vref ! for pure growth, vrel/vfrag gives vrel in m/s
 elseif (ifrag > 0) then
    call comp_snow_line(xyzh,cs,izone)
    select case(izone)
    case(0)
       dustprop(3) = vrel/vfrag
    case(1)
       dustprop(3) = vrel/vfragin
    case(2)
       dustprop(3) = vrel/vfragout
    end select
 endif
end subroutine get_vrelonvfrag

!----------------------------------------------------------------------------
!+
! Get the location of a given particle (dust or gas or mixture) with respect
! to a (or more) snow line(s)
!+
!----------------------------------------------------------------------------
subroutine comp_snow_line(xyzh,cs,izone)
 use eos,    only:temperature_coef,gmw
 integer, intent(out) :: izone
 real, intent(in)     :: xyzh(:),cs
 real                 :: cs_snow,r

 select case(isnow)
 case(0)
    izone = 0
 case(1)
    r = sqrt(xyzh(1)**2 + xyzh(2)**2 + xyzh(3)**2)
    if (r<=rsnow) izone = 1
    if (r>rsnow) izone = 2
 case(2)
    cs_snow = sqrt(Tsnow/(temperature_coef*gmw))
    if (cs > cs_snow) izone = 1
    if (cs < cs_snow) izone = 2
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
 use infile_utils,        only:write_inopt
 integer, intent(in)        :: iunit

 write(iunit,"(/,a)") '# options controlling growth'
 call write_inopt(ifrag,'ifrag','dust fragmentation (0=off,1=on,2=Kobayashi)',iunit)
 if (ifrag /= 0) then
    call write_inopt(gsizemincgs,'grainsizemin','minimum grain size in cm',iunit)
    call write_inopt(isnow,'isnow','snow line (0=off,1=position based,2=temperature based)',iunit)
    if (isnow == 1) call write_inopt(rsnow,'rsnow','position of the snow line in AU',iunit)
    if (isnow == 2) call write_inopt(Tsnow,'Tsnow','snow line condensation temperature in K',iunit)
    if (isnow == 0) call write_inopt(vfragSI,'vfrag','uniform fragmentation threshold in m/s',iunit)
    if (isnow > 0) then
       call write_inopt(vfraginSI,'vfragin','inward fragmentation threshold in m/s',iunit)
       call write_inopt(vfragoutSI,'vfragout','outward fragmentation threshold in m/s',iunit)
    endif
 endif

end subroutine write_options_growth

!-----------------------------------------------------------------------
!+
!  Read growth options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_growth(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)        :: name,valstring
 logical,intent(out)                        :: imatch,igotall
 integer,intent(out)                        :: ierr

 integer,save                                        :: ngot = 0

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('ifrag')
    read(valstring,*,iostat=ierr) ifrag
    ngot = ngot + 1
 case('grainsizemin')
    read(valstring,*,iostat=ierr) gsizemincgs
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
 case default
    imatch = .false.
 end select

 if ((ifrag <= 0) .and. ngot == 1) igotall = .true.
 if (isnow == 0) then
    if (ngot == 4) igotall = .true.
 elseif (isnow > 0) then
    if (ngot == 6) igotall = .true.
 else
    igotall = .false.
 endif
end subroutine read_options_growth

!-----------------------------------------------------------------------
!+
!  Write growth options to the .setup file
!+
!-----------------------------------------------------------------------
subroutine write_growth_setup_options(iunit)
 use infile_utils,      only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for growth and fragmentation of dust'

 call write_inopt(ifrag,'ifrag','fragmentation of dust (0=off,1=on,2=Kobayashi)',iunit)
 call write_inopt(isnow,'isnow','snow line (0=off,1=position based,2=temperature based)',iunit)
 call write_inopt(rsnow,'rsnow','snow line position in AU',iunit)
 call write_inopt(Tsnow,'Tsnow','snow line condensation temperature in K',iunit)
 call write_inopt(vfragSI,'vfrag','uniform fragmentation threshold in m/s',iunit)
 call write_inopt(vfraginSI,'vfragin','inward fragmentation threshold in m/s',iunit)
 call write_inopt(vfragoutSI,'vfragout','inward fragmentation threshold in m/s',iunit)
 call write_inopt(gsizemincgs,'grainsizemin','minimum allowed grain size in cm',iunit)

end subroutine write_growth_setup_options

!-----------------------------------------------------------------------
!+
!  Read growth options from the .setup file
!+
!-----------------------------------------------------------------------
subroutine read_growth_setup_options(db,nerr)
 use infile_utils,    only:read_inopt,inopts
 type(inopts), allocatable, intent(inout) :: db(:)
 integer, intent(inout)                   :: nerr

 call read_inopt(ifrag,'ifrag',db,min=-1,max=2,errcount=nerr)
 if (ifrag > 0) then
    call read_inopt(isnow,'isnow',db,min=0,max=2,errcount=nerr)
    call read_inopt(grainsizemin,'grainsizemin',db,min=1.e-5,errcount=nerr)
    select case(isnow)
    case(0)
       call read_inopt(vfrag,'vfrag',db,min=0.,errcount=nerr)
    case(1)
       call read_inopt(rsnow,'rsnow',db,min=0.,errcount=nerr)
       call read_inopt(vfragin,'vfragin',db,min=0.,errcount=nerr)
       call read_inopt(vfragout,'vfragout',db,min=0.,errcount=nerr)
    case(2)
       call read_inopt(Tsnow,'Tsnow',db,min=0.,errcount=nerr)
       call read_inopt(vfragin,'vfragin',db,min=0.,errcount=nerr)
       call read_inopt(vfragout,'vfragout',db,min=0.,errcount=nerr)
    end select
 endif

end subroutine read_growth_setup_options

!-----------------------------------------------------------------------
!+
!  Update dustprop and make sure dustprop(1,:) isn't too small
!+
!-----------------------------------------------------------------------
subroutine check_dustprop(npart,size)
 use part,                only:iamtype,iphase,idust
 real,intent(inout)        :: size(:)
 integer,intent(in)        :: npart
 integer                   :: i

 do i=1,npart
    if (iamtype(iphase(i))==idust) then
       if (ifrag > 0 .and. size(i) < grainsizemin) size(i) = grainsizemin
    endif
 enddo
end subroutine check_dustprop

!-----------------------------------------------------------------------
!+
!  Set dustprop (used by moddump)
!+
!-----------------------------------------------------------------------
subroutine set_dustprop(npart)
 use dust, only:grainsizecgs,graindenscgs
 use part, only:dustprop
 integer,intent(in) :: npart
 integer            :: i

 do i=1,npart
    dustprop(1,i) = grainsizecgs / udist
    dustprop(2,i) = graindenscgs / unit_density
 enddo

end subroutine set_dustprop

!--Compute the relative velocity following Stepinski & Valageas (1997)
real function vrelative(St,dv,Vt)
 real, intent(in) :: St,dv,Vt
 real             :: Sc

 !--compute Schmidt number Sc
 Sc = (1.+St)*sqrt(1.+dv**2/Vt**2)
 !--then compute vrel
 vrelative = sqrt(2.)*Vt*sqrt(Sc-1.)/(Sc)

 return
end function vrelative

end module growth
