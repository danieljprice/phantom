!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  test common envelope - put point source star next to gas sphere
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, dim, externalforces, infile_utils, io,
!    options, part, physcon, prompting, rho_profile, setbinary, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,&
                             delete_dead_or_accreted_particles,mhd,rhoh,shuffle_part,kill_particle
 use setbinary,         only:set_binary
 use units,             only:umass,udist,utime
 use physcon,           only:au,solarm,solarr,gg
 use centreofmass,      only:reset_centreofmass,get_centreofmass
 use prompting,         only:prompt
 use options,           only:iexternalforce
 use externalforces,    only:omega_corotate,iext_corotate
 use infile_utils,      only:open_db_from_file,inopts,read_inopt,close_db
 use rho_profile,       only:read_red_giant_file
 use dim,               only:maxptmass
 use io,                only:fatal

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i, ierr, setup_case, two_sink_case = 1, npts, rhomaxi, n
 real                   :: primary_mass, companion_mass, mass_ratio, a, e, omega_vec(3), omegacrossr(3), vr = 0.0, hsoft_default = 3
 real                   :: hacc1,hacc2,mcore,comp_shift=100, sink_dist, vel_shift
 real                   :: mcut,rcut,Mstar,radi,rhopart,rhomax = 0.0
 integer, parameter     :: ng_max = 5000
 real                   :: r(ng_max),den(ng_max),pres(ng_max),temp(ng_max),enitab(ng_max)
 logical                :: corotate_answer
 character(len=20)      :: filename = 'binary.in'
 character(len=100)     :: densityfile
 type(inopts), allocatable :: db(:)

!control: more than one sink particle already in the code could cause problems
 if (nptmass > 2) then
    stop 'ERROR: Number of sink particles > 1'
 elseif (nptmass == 2) then
    print*, 'Two sinks present. If this is intentional, then choose option below.'

    print "(2(/,a))",'1) Switch from corotating frame to normal frame', &
                  '2) Change the position of the companion'
    call prompt('Select option above : ',two_sink_case)
    select case(two_sink_case)

    case(1)
       call prompt('Please write the name of the input file : ',filename)
       call open_db_from_file(db,filename,20,ierr)
       call read_inopt(omega_corotate,'omega_corotate',db)
       call close_db(db)
       omega_vec = (/ 0.,0.,omega_corotate /)
       do i=1,npart
          call cross(omega_vec,xyzh(:3,i),omegacrossr)
          vxyzu(1,i) = vxyzu(1,i) + omegacrossr(1)
          vxyzu(2,i) = vxyzu(2,i) + omegacrossr(2)
          vxyzu(3,i) = vxyzu(3,i) + omegacrossr(3)
       enddo
       do i=1,nptmass
          call cross(omega_vec,xyzmh_ptmass(:3,i),omegacrossr)
          vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + omegacrossr(1)
          vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + omegacrossr(2)
          vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + omegacrossr(3)
       enddo

    case(2)
       call prompt('How many code units to shift companion (+ve is towards primary)?',comp_shift)
       sink_dist = sqrt((xyzmh_ptmass(1,1)-xyzmh_ptmass(1,2))**2 &
                 + (xyzmh_ptmass(2,1)-xyzmh_ptmass(2,2))**2 &
                 + (xyzmh_ptmass(3,1)-xyzmh_ptmass(3,2))**2)

       xyzmh_ptmass(1,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1)) - xyzmh_ptmass(1,2))
       xyzmh_ptmass(2,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(2,2)-xyzmh_ptmass(2,1)) - xyzmh_ptmass(2,2))
       xyzmh_ptmass(3,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(3,2)-xyzmh_ptmass(3,1)) - xyzmh_ptmass(3,2))

       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
       iexternalforce = iext_corotate
       omega_corotate = sqrt((sink_dist-comp_shift)*(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2)))/(sink_dist-comp_shift)**2

       do i=1,npart
          vxyzu(1,i) = 0.0
          vxyzu(2,i) = 0.0
          vxyzu(3,i) = 0.0
       enddo

       do i=1,nptmass
          vxyz_ptmass(1,i) = 0.0
          vxyz_ptmass(1,i) = 0.0
          vxyz_ptmass(1,i) = 0.0
       enddo
    case(3)
       sink_dist = sqrt((xyzmh_ptmass(1,1)-xyzmh_ptmass(1,2))**2 &
                 + (xyzmh_ptmass(2,1)-xyzmh_ptmass(2,2))**2 &
                 + (xyzmh_ptmass(3,1)-xyzmh_ptmass(3,2))**2)

       print*, utime, umass, udist
       vel_shift = sink_dist/(5.0*60*60*24*365/utime)

       call prompt('Give velocity to add in direction of the primary : ',vel_shift, 0.0)

       vxyz_ptmass(1,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1)) - vxyz_ptmass(1,2))
       vxyz_ptmass(2,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(2,2)-xyzmh_ptmass(2,1)) - vxyz_ptmass(2,2))
       vxyz_ptmass(3,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(3,2)-xyzmh_ptmass(3,1)) - vxyz_ptmass(3,2))

       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
    end select
 else

!choose what to do with the star: set a binary or setup a magnetic field
    print "(3(/,a))",'1) setup a binary system', &
                  '2) setup a magnetic field in the star', &
                  '3) manually create sink in core'

    setup_case = 1
    call prompt('Choose a setup option ',setup_case,1,3)

    select case(setup_case)

    case(1)

       !takes necessary inputs from user 1
       print*, 'Current mass unit is ', umass,'g):'
       companion_mass = 0.6
       call prompt('Enter companion mass in code units',companion_mass,0.)

       print*, 'Current length unit is ', udist ,'cm):'
       a = 100.
       call prompt('Enter orbit semi-major axis in code units', a, 0.0)

       e = 0.0
       call prompt('Enter orbit eccentricity', e, 0.0, 1.0)
       call prompt('Enter companion radial velocity', vr)

       print*, 'Current length unit is ', udist ,'cm):'
       hacc1 = 0.0
       hacc2 = 0.0
       call prompt('Enter accretion radius for the primary in code units', hacc1, 0.0)
       call prompt('Enter accretion radius for the companion in code units', hacc2, 0.0)

       corotate_answer = .false.
       call prompt('Do you want a corotating frame with a corotating binary?', corotate_answer)

       !resets to (0,0,0) position and velocity of centre of mass for whole system before creating the binary
       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

       !removes the dead or accreted particles for a correct total mass computation
       call delete_dead_or_accreted_particles(npart,npartoftype)
       print*,' Got ',npart,npartoftype(igas),' after deleting accreted particles'

       !sets up the binary system orbital parameters
       if (nptmass > 0) then
          mcore = xyzmh_ptmass(4,1)
       else
          mcore = 0.
       endif

       primary_mass = npartoftype(igas) * massoftype(igas) + mcore
       mass_ratio = companion_mass / primary_mass

       !sets the binary
       !corotating frame
       if (corotate_answer) then
          !turns on corotation
          iexternalforce = iext_corotate

          call set_binary(primary_mass,mass_ratio,a,e,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,omega_corotate)

          print "(/,a,es18.10,/)", ' The angular velocity in the corotating frame is: ', omega_corotate

          !sets all the gas velocities in corotating frame to 0, implying that the binary is corotating
          !for the moment only a corotating binary can be built in the corotating frame
          do i=1,npart
             vxyzu(1,i) = 0.0
             vxyzu(2,i) = 0.0
             vxyzu(3,i) = 0.0
          enddo
          !non corotating frame
       else
          call set_binary(primary_mass,mass_ratio,a,e,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass)
       endif

       !"set_binary" newly created sinks shifted to the first and second element of the xyzmh_ptmass array (original sink overwritten)
       if (nptmass > 2) then
          !new primary from pos 2 to 1
          !positions
          xyzmh_ptmass(1,1) = xyzmh_ptmass(1,2)
          xyzmh_ptmass(2,1) = xyzmh_ptmass(2,2)
          xyzmh_ptmass(3,1) = xyzmh_ptmass(3,2)

          !velocities
          vxyz_ptmass(1,1) = vxyz_ptmass(1,2)
          vxyz_ptmass(2,1) = vxyz_ptmass(2,2)
          vxyz_ptmass(3,1) = vxyz_ptmass(3,2)

          !new companion from pos 3 to 2
          !positions
          xyzmh_ptmass(:,2) = xyzmh_ptmass(:,3)

          !velocities
          vxyz_ptmass(1,2) = vxyz_ptmass(1,3) + vr
          vxyz_ptmass(2,2) = vxyz_ptmass(2,3)
          vxyz_ptmass(3,2) = vxyz_ptmass(3,3)
       endif

       !takes necessary inputs from user 2 (the softening lengths for the sinks have to be taken in input after using the "set_binary" function since it resets them)
       xyzmh_ptmass(ihsoft,1) = xyzmh_ptmass(ihsoft,1)
       xyzmh_ptmass(ihsoft,2) = xyzmh_ptmass(ihsoft,1)
       call prompt('Enter softening length for primary',xyzmh_ptmass(ihsoft,1),0.)
       call prompt('Enter softening length for secondary',xyzmh_ptmass(ihsoft,2),0.)

       !shifts gas to the primary point mass created in 'set_binary'
       do i=1,npart
          !positions
          xyzh(1,i) = xyzh(1,i) + xyzmh_ptmass(1,1)
          xyzh(2,i) = xyzh(2,i) + xyzmh_ptmass(2,1)
          xyzh(3,i) = xyzh(3,i) + xyzmh_ptmass(3,1)

          !velocities
          vxyzu(1,i) = vxyzu(1,i) + vxyz_ptmass(1,1)
          vxyzu(2,i) = vxyzu(2,i) + vxyz_ptmass(2,1)
          vxyzu(3,i) = vxyzu(3,i) + vxyz_ptmass(3,1)
       enddo

       !deletes third point mass
       nptmass = 2

       !resets to (0,0,0) position and velocity of centre of mass for whole system after creating the binary
       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    case(2)

       if (mhd) then
          print "(/,a,/)", 'Automatic insertion of the magnetic field through the setBfield module'
       else
          print "(/,a,/)", 'Code not compiled with MHD=yes, no changes to the dump have been made'
       endif

    case(3)

       densityfile = 'P12_Phantom_Profile.data'
       call prompt('Enter filename of the input stellar profile', densityfile)
       call prompt('Enter mass of the created point mass core', mcut)
       call prompt('Enter softening length of the point mass', hsoft_default)
       call read_red_giant_file(trim(densityfile),ng_max,npts,r,den,pres,temp,enitab,Mstar,ierr,mcut,rcut)

       do i=1,npart
          rhopart = rhoh(xyzh(4,i), massoftype(igas))
          if (rhopart > rhomax) then
             rhomax = rhopart
             rhomaxi = i
          endif
       enddo

       nptmass = nptmass + 1
       if (nptmass > maxptmass) call fatal('ptmass_create','nptmass > maxptmass')
       n = nptmass
       xyzmh_ptmass(:,n)   = 0.  ! zero all quantities by default
       xyzmh_ptmass(1:3,n) = xyzh(1:3,rhomaxi)
       xyzmh_ptmass(4,n)   = 0.  ! zero mass
       xyzmh_ptmass(ihsoft,n) = hsoft_default
       vxyz_ptmass(:,n) = 0.     ! zero velocity, get this by accreting


       do i=1,npart
          radi = sqrt((xyzh(1,i)-xyzh(1,rhomaxi))**2 + &
                   (xyzh(2,i)-xyzh(2,rhomaxi))**2 + &
                   (xyzh(3,i)-xyzh(3,rhomaxi))**2)
          if (radi < rcut) then
             xyzmh_ptmass(4,n) = xyzmh_ptmass(4,n) + massoftype(igas)
             npartoftype(igas) = npartoftype(igas) - 1
             call kill_particle(i)
          endif
       enddo

       call shuffle_part(npart)

    end select

 endif

 return
end subroutine modify_dump

subroutine cross(a,b,c)

 ! Return the vector cross product of two 3d vectors
 implicit none
 real,intent(in),dimension(3)  :: a,b
 real,intent(out),dimension(3) :: c

 c(1) = a(2)*b(3)-b(2)*a(3)
 c(2) = a(3)*b(1)-b(3)*a(1)
 c(3) = a(1)*b(2)-b(1)*a(2)

end subroutine cross

end module moddump

