!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
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
!  Kobayashi & Tanaka (2008)
!
!  OWNER: Arnaud Vericel
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    ifrag             -- fragmentation of dust (0=off, 1=on, 2=Kobayashi)
!    grainsizemin      -- minimum grain size in cm
!    isnow             -- snow line (0=off, 1=position based, 2=temperature based)
!    rsnow             -- snow line position in AU
!	 Tsnow			   -- snow line condensation temperature in K
!    vfragin           -- inward fragmentation threshold in m/s
!    vfragout          -- outward fragmentation threshold in m/s
!    vfrag			   -- uniform fragmentation threshold in m/s
!
!  DEPENDENCIES: part,physcon,units,eos,infile_utiles,io,dim,options
!+
!--------------------------------------------------------------------------

module growth
 use units,	only:udist,unit_density,unit_velocity
 use physcon,	only:au
 implicit none
 
 !--Default values for the growth and fragmentation of dust in the input file
 integer, public	:: ifrag	= 1         
 integer, public	:: isnow	= 0         
 
 real, public		:: grainsizemin	= 1.e-3
 real, public		:: rsnow	= 100.
 real, public		:: Tsnow	= 20.
 real, public		:: vfrag	= 15.   
 real, public		:: vfragin	= 5.       
 real, public		:: vfragout	= 15.    
 
 public			:: evolve_grainsize,get_vrelonvfrag
 public			:: write_options_growth,read_options_growth,print_growthinfo,init_growth
 public			:: vrelative
 
 integer, public	:: i
 real, public		:: r,T
 real, public		:: vrel
 
contains
	
!------------------------------------------------
!+
!  Initialise variables for computing growth rate
!+
!------------------------------------------------
subroutine init_growth(ierr)
 use io,	only:error
 use part,	only:dustprop
 use dim,	only:maxp_growth
 integer, intent(out) :: ierr
 ierr = 0
 
 !--Convert in code units
 dustprop(1,:)	= dustprop(1,:) / udist
 dustprop(2,:)	= dustprop(2,:) / unit_density
 grainsizemin	= grainsizemin / udist
 vfrag			= vfrag * 100 / unit_velocity
 vfragin		= vfragin * 100 / unit_velocity
 vfragout		= vfragout * 100 / unit_velocity
 rsnow			= rsnow * au / udist
 
 !-- Check that all the parameters are > 0 when needed
 do i=1,maxp_growth
 	if (dustprop(1,i) < 0) then
		call error('init_growth','grainsize < 0',var='dustprop',val=dustprop(1,i))
 		ierr = 1 
 	endif
 	if (dustprop(2,i) < 0) then
		call error('init_growth','graindens < 0',var='dustprop',val=dustprop(2,i))
    	ierr = 1 
    endif
 	if (dustprop(3,i) < 0) then
		call error('init_growth','vrel/vfrag < 0',var='dustprop',val=dustprop(3,i))
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
 if (ifrag == 2) write(iprint,"(a)")    ' Using growth with Kobayashi fragmentation model       		 '
 if (ifrag > 0) then
 	write(iprint,"(2(a,1pg10.3),a)")' grainsizemin = ',grainsizemin*udist,' cm = ',grainsizemin,' (code units)'
 	if (isnow == 1) then
 		write(iprint,"(a)")              ' ===> Using position based snow line <===	           					    '
		write(iprint,"(2(a,1pg10.3),a)") '   rsnow = ',rsnow*udist/au            ,'    AU = ',rsnow   ,' (code units)'    
    endif
	if (isnow == 2) then
  		write(iprint,"(a)")              ' ===> Using temperature based snow line <===							    '
 		write(iprint,"(2(a,1pg10.3),a)") '   Tsnow = ',Tsnow				        ,'     K = ',Tsnow   ,' (code units)'
 	endif
	if (isnow == 0) then
		write(iprint,"(2(a,1pg10.3),a)") '   vfrag = ',vfrag*unit_velocity/100   ,'   m/s = ',vfrag   ,' (code units)'
    else
		write(iprint,"(2(a,1pg10.3),a)") ' vfragin = ',vfragin*unit_velocity/100 ,'   m/s = ',vfragin ,' (code units)'
		write(iprint,"(2(a,1pg10.3),a)") ' vfragin = ',vfragout*unit_velocity/100,'   m/s = ',vfragout,' (code units)'
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
subroutine evolve_grainsize(xyzh,dustprop,rhod,spsound,T,St,dt) 
 real, intent(inout)	:: dustprop(3)
 real, intent(in)		:: rhod,spsound,St,dt,T
 real, intent(in)		:: xyzh(4)
 
 !--compute vrel and vrel/vfrag from get_vrelonvfrag subroutine
 call get_vrelonvfrag(xyzh,dustprop(3),spsound,T,St)
 !
 !--If statements to evolve grainsize
 !--dustprop(1)= size, dustprop(2) = intrinsic density, dustprop(3) = local vrel/vfrag
 !
 if (dustprop(3) >= 1.) then ! vrel/vfrag < 1 --> growth
 	dustprop(1) = dustprop(1) + rhod/dustprop(2)*vrel*dt
 elseif (dustprop(3) < 1. .and. ifrag > 0) then ! vrel/vfrag > 1 --> fragmentation 
	select case(ifrag)
	case(1)
		dustprop(1) = dustprop(1) - rhod/dustprop(2)*vrel*dt ! Symmetrical of Stepinski & Valageas
	case(2)
	    dustprop(1) = dustprop(1) - rhod/dustprop(2)*vrel*(dustprop(3)**2)/(1+dustprop(3)**2) ! Kobayashi model
	case default
    end select
	
	if (dustprop(1) < grainsizemin) then
		dustprop(1) = grainsizemin ! Prevent dust from becoming too small
	endif
 endif
end subroutine evolve_grainsize
 
!-----------------------------------------------------------------------
!+
!  Compute the local ratio vrel/vfrag and vrel
!+
!-----------------------------------------------------------------------
subroutine get_vrelonvfrag(xyzh,vrelonvfrag,spsound,T,St) 
 use eos,			only:ieos	
 use options,		only:alpha
 real, intent(in)	:: xyzh(4)
 real, intent(in)	:: spsound,St,T
 real, intent(out)	:: vrelonvfrag
 
 !--Compute relative velocity of the dust particle
 vrel = vrelative(spsound,St,alpha)
 !
 !--If statements to compute local ratio vrel/vfrag 
 !
 select case(isnow)
 case(0) !--uniform vfrag
	vrelonvfrag = vrel / vfrag
 case(1) !--position based snow line in cylindrical geometry
 	r = sqrt(xyzh(1)**2+xyzh(2)**2)
	if (r < rsnow) vrelonvfrag = vrel / vfragin
	if (r > rsnow) vrelonvfrag = vrel / vfragout
 case(2) !--temperature based snow line wrt eos
	if (T > Tsnow) vrelonvfrag = vrel / vfragin
	if (T < Tsnow) vrelonvfrag = vrel / vfragout
 case default
 	vrelonvfrag = 0.
	vrel = 0.
 end select
 
end subroutine get_vrelonvfrag

!-----------------------------------------------------------------------
!+
!  Write growth options in the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_growth(iunit)
 use infile_utils,	only:write_inopt
 integer, intent(in)	:: iunit

 write(iunit,"(/,a)") '# options controlling growth'
 call write_inopt(ifrag,'ifrag','dust fragmentation (0=off,1=on,2=Kobayashi)',iunit)
 if (ifrag /= 0) then
 	call write_inopt(grainsizemin,'grainsizemin','minimum grain size in cm',iunit)
 	call write_inopt(isnow,'isnow','snow line (0=off,1=position based,2=temperature based)',iunit)
 	if (isnow == 1) call write_inopt(rsnow,'rsnow','position of the snow line in AU',iunit)
 	if (isnow == 2) call write_inopt(rsnow,'Tsnow','snow line condensation temperature in K',iunit)
 	if (isnow == 0) call write_inopt(vfrag,'vfrag','uniform fragmentation threshold in m/s',iunit)
 	if (isnow > 0) then
 		call write_inopt(vfragin,'vfragin','inward fragmentation threshold in m/s',iunit)
 		call write_inopt(vfragout,'vfragout','outward fragmentation threshold in m/s',iunit)
	endif
 endif

end subroutine write_options_growth

!-----------------------------------------------------------------------
!+
!  Read growth options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_growth(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)	:: name,valstring
 logical,intent(out)			:: imatch,igotall
 integer,intent(out)			:: ierr
 
 integer,save					:: ngot = 0
  
 imatch  = .true.
 igotall = .false.
 
 select case(trim(name))
 case('ifrag')
    read(valstring,*,iostat=ierr) ifrag
    ngot = ngot + 1
 case('grainsizemin')
    read(valstring,*,iostat=ierr) grainsizemin
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
    read(valstring,*,iostat=ierr) vfrag
	ngot = ngot + 1
 case('vfragin')
    read(valstring,*,iostat=ierr) vfragin
	ngot = ngot + 1
 case('vfragout')
    read(valstring,*,iostat=ierr) vfragout
	ngot = ngot + 1
    case default
    imatch = .false.
 end select
 
 if (ifrag == 0 .and. ngot == 1) igotall = .true.
 if (isnow == 0) then
 	if (ngot == 4) igotall = .true.
 elseif (isnow > 0) then
 	if (ngot == 6) igotall = .true.
 else
 	igotall = .false.
 endif
end subroutine read_options_growth

!--Compute the relative velocity following Stepinski & Valageas (1997)
real function vrelative(spsound,St,alpha)
 real, intent(in) :: spsound,St,alpha
 real, parameter  :: Ro = 3.
 vrel = sqrt(2**(1.5)*Ro*alpha)*spsound*sqrt(St)/(1+St)
 return
end function vrelative

end module growth

