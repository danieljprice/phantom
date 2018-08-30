!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
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
 use rho_profile,       only:read_mesa_file
 use dim,               only:maxptmass
 use io,                only:fatal

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i, ierr, setup_case, two_sink_case = 1, three_sink_case = 1, npts, irhomax, n
 integer                :: iremove = 2
 real                   :: primary_mass, companion_mass_1, companion_mass_2, mass_ratio
 real                   :: a1, a2, e, omega_vec(3), omegacrossr(3), vr = 0.0, hsoft_default = 3
 real                   :: hacc1,hacc2,hacc3,mcore,comp_shift=100, sink_dist, vel_shift
 real                   :: mcut,rcut,Mstar,radi,rhopart,rhomax = 0.0
 integer, parameter     :: ng_max = 5000
 real                   :: r(ng_max),den(ng_max),pres(ng_max),temp(ng_max),enitab(ng_max)
 logical                :: corotate_answer
 character(len=20)      :: filename = 'binary.in'
 character(len=100)     :: densityfile
 type(inopts), allocatable :: db(:)

!control: more than one sink particle already in the code could cause problems
 if (nptmass > 3) then
    stop 'ERROR: Number of sink particles > 1'
 elseif (nptmass == 3) then
    print "(1(/,a))",'1) Remove a sink from the simulation'
    call prompt('Select option above : ',three_sink_case)
    select case(three_sink_case)

    case(1)
       do i=1,nptmass
          write(*,'(A,I2,A,ES10.3,A,ES10.3)') 'Point mass ',i,': M = ',xyzmh_ptmass(4,i),&
                                              ' and radial position = ',sqrt(dot_product(xyzmh_ptmass(1:3,i),xyzmh_ptmass(1:3,i)))
       enddo
       call prompt('Which sink would you like to remove : ',iremove)
       if (iremove == 3) then
          xyzmh_ptmass(:,iremove) = 0.
          vxyz_ptmass(:,iremove) = 0.
          nptmass = 2
       elseif (iremove == 2) then
          xyzmh_ptmass(:,2) = xyzmh_ptmass(:,3)
          vxyz_ptmass(:,2) = vxyz_ptmass(:,3)
          nptmass = 2
       endif
    end select

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
    print "(4(/,a))",'1) setup a binary system', &
                  '2) setup a magnetic field in the star', &
                  '3) manually create sink in core', &
                  '4) setup trinary system'

    setup_case = 1
    call prompt('Choose a setup option ',setup_case,1,4)

    select case(setup_case)

    case(1)

       !takes necessary inputs from user 1
       print*, 'Current mass unit is ', umass,'g):'
       companion_mass_1 = 0.6
       call prompt('Enter companion mass in code units',companion_mass_1,0.)

       print*, 'Current length unit is ', udist ,'cm):'
       a1 = 100.
       call prompt('Enter orbit semi-major axis in code units', a1, 0.0)

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
       mass_ratio = companion_mass_1 / primary_mass

       !sets the binary
       !corotating frame
       if (corotate_answer) then
          !turns on corotation
          iexternalforce = iext_corotate

          call set_binary(primary_mass,mass_ratio,a1,e,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,omega_corotate)

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
          call set_binary(primary_mass,mass_ratio,a1,e,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass)
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
       call read_mesa_file(trim(densityfile),ng_max,npts,r,den,pres,temp,enitab,Mstar,ierr,mcut,rcut)

       irhomax = 1
       do i=1,npart
          rhopart = rhoh(xyzh(4,i), massoftype(igas))
          if (rhopart > rhomax) then
             rhomax = rhopart
             irhomax = i
          endif
       enddo

       nptmass = nptmass + 1
       if (nptmass > maxptmass) call fatal('ptmass_create','nptmass > maxptmass')
       n = nptmass
       xyzmh_ptmass(:,n)   = 0.  ! zero all quantities by default
       xyzmh_ptmass(1:3,n) = xyzh(1:3,irhomax)
       xyzmh_ptmass(4,n)   = 0.  ! zero mass
       xyzmh_ptmass(ihsoft,n) = hsoft_default
       vxyz_ptmass(:,n) = 0.     ! zero velocity, get this by accreting


       do i=1,npart
          radi = sqrt((xyzh(1,i)-xyzh(1,irhomax))**2 + &
                   (xyzh(2,i)-xyzh(2,irhomax))**2 + &
                   (xyzh(3,i)-xyzh(3,irhomax))**2)
          if (radi < rcut) then
             xyzmh_ptmass(4,n) = xyzmh_ptmass(4,n) + massoftype(igas)
             npartoftype(igas) = npartoftype(igas) - 1
             call kill_particle(i)
          endif
       enddo

       call shuffle_part(npart)

    case(4)

       !takes necessary inputs from user 1
       print*, 'Current mass unit is ', umass,'g):'
       companion_mass_1 = 0.0095
       call prompt('Enter 1st companion mass in code units',companion_mass_1,0.)
       companion_mass_2 = 0.0095
       call prompt('Enter 2nd companion mass in code units',companion_mass_2,0.)

       print*, 'Current length unit is ', udist ,'cm):'
       a1 = 166.5
       call prompt('Enter 1st companion orbit semi-major axis in code units', a1, 0.0)
       a2 = 336.8
       call prompt('Enter 2nd companion orbit semi-major axis in code units', a2, 0.0)

       print*, 'Current length unit is ', udist ,'cm):'
       hacc1 = 0.0
       hacc2 = 0.0
       hacc3 = 0.0
       call prompt('Enter accretion radius for the primary in code units', hacc1, 0.0)
       call prompt('Enter accretion radius for the 1st companion in code units', hacc2, 0.0)
       call prompt('Enter accretion radius for the 2nd companion in code units', hacc3, 0.0)

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

       call set_trinary(primary_mass,companion_mass_1,companion_mass_2,&
                        a1,a2,hacc1,hacc2,hacc3,&
                        xyzmh_ptmass,vxyz_ptmass,nptmass)


       if (nptmass > 3) then
          xyzmh_ptmass(1:3,1) = xyzmh_ptmass(1:3,2)
          xyzmh_ptmass(:,2) = xyzmh_ptmass(:,3)
          xyzmh_ptmass(:,3) = xyzmh_ptmass(:,4)

          vxyz_ptmass(:,1) = vxyz_ptmass(:,2)
          vxyz_ptmass(:,2) = vxyz_ptmass(:,3)
          vxyz_ptmass(:,3) = vxyz_ptmass(:,4)
       endif

       xyzmh_ptmass(ihsoft,1) = xyzmh_ptmass(ihsoft,1)
       xyzmh_ptmass(ihsoft,2) = xyzmh_ptmass(ihsoft,1)
       xyzmh_ptmass(ihsoft,3) = xyzmh_ptmass(ihsoft,1)
       call prompt('Enter softening length for primary',xyzmh_ptmass(ihsoft,1),0.)
       call prompt('Enter softening length for secondary',xyzmh_ptmass(ihsoft,2),0.)
       call prompt('Enter softening length for tertiary',xyzmh_ptmass(ihsoft,3),0.)


       !shifts gas to the primary point mass created in 'set_binary'
       do i=1,npart
          !positions
          xyzh(1:3,i) = xyzh(1:3,i) + xyzmh_ptmass(1:3,1)

          !velocities
          vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz_ptmass(1:3,1)
       enddo

       !deletes third point mass
       nptmass = 3

       !resets to (0,0,0) position and velocity of centre of mass for whole system after creating the binary
       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

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

subroutine set_trinary(mprimary,msecondary,mtertiary,semimajoraxis12,semimajoraxis13,&
                      accretion_radius1,accretion_radius2,accretion_radius3,&
                      xyzmh_ptmass,vxyz_ptmass,nptmass)
 use part,    only:ihacc,ihsoft
 real,    intent(in)    :: mprimary,msecondary,mtertiary
 real,    intent(in)    :: semimajoraxis12,semimajoraxis13
 real,    intent(in)    :: accretion_radius1,accretion_radius2,accretion_radius3
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass

 integer :: i1,i2,i3
 real    :: m1,m2,m3,mtot,dx12(3),dx13(3),dv12(3),dv13(3)
 real    :: x1(3),x2(3),x3(3),v1(3),v2(3),v3(3)

 i1 = nptmass + 1
 i2 = nptmass + 2
 i3 = nptmass + 3
 nptmass = nptmass + 3

 ! masses
 m1 = mprimary
 m2 = msecondary
 m3 = mtertiary
 mtot = m1 + m2 + m3

!
!--check for stupid parameter choices
!
! if (mprimary <= 0.)      stop 'ERROR: primary mass <= 0'
! if (massratio < 0.)      stop 'ERROR: binary mass ratio < 0'
! if (semimajoraxis <= 0.) stop 'ERROR: semi-major axis <= 0'
! if (eccentricity > 1. .or. eccentricity < 0.) &
!    stop 'ERROR: eccentricity must be between 0 and 1'

 dx12 = (/semimajoraxis12,0.,0./)
 dv12 = (/0.,sqrt((m1+m2)/dx12(1)),0./)

 dx13 = (/semimajoraxis13,0.,0./)
 dv13 = (/0.,sqrt(mtot/dx13(1)),0./)

 ! positions of each star so centre of mass is at zero
 x1 = -(dx12*m2 + dx13*m3)/mtot
 x2 = (dx12*m1 + dx12*m3 - dx13*m3)/mtot
 x3 = (dx13*m1 + dx13*m2 - dx12*m2)/mtot

 ! velocities
 v1 = -(dv12*m2 + dv13*m3)/mtot
 v2 = (dv12*m1 + dv12*m3 - dv13*m3)/mtot
 v3 = (dv13*m1 + dv13*m2 - dv12*m2)/mtot

!
!--positions and accretion radii
!
 xyzmh_ptmass(:,i1:i3) = 0.
 xyzmh_ptmass(1:3,i1) = x1
 xyzmh_ptmass(1:3,i2) = x2
 xyzmh_ptmass(1:3,i3) = x3
 xyzmh_ptmass(4,i1) = m1
 xyzmh_ptmass(4,i2) = m2
 xyzmh_ptmass(4,i3) = m3
 xyzmh_ptmass(ihacc,i1) = accretion_radius1
 xyzmh_ptmass(ihacc,i2) = accretion_radius2
 xyzmh_ptmass(ihacc,i3) = accretion_radius3
 xyzmh_ptmass(ihsoft,i1) = 0.0
 xyzmh_ptmass(ihsoft,i2) = 0.0
 xyzmh_ptmass(ihsoft,i3) = 0.0
!
!--velocities
!
 vxyz_ptmass(:,i1) = v1
 vxyz_ptmass(:,i2) = v2
 vxyz_ptmass(:,i3) = v3

end subroutine set_trinary


end module moddump
