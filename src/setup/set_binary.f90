!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setbinary
!
! This module is contains utilities for setting up binaries
!   Our conventions for binary orbital parameters are consistent with
!   those produced by the imorbel code (Pearce, Wyatt & Kennedy 2015)
!   which can be used to produce orbits matching observed orbital
!   arcs of binary companions on the sky
!
! :References:
!   Eggleton (1983) ApJ 268, 368-369 (ref:eggleton83)
!   Lucy (2014), A&A 563, A126
!   Pearce, Wyatt & Kennedy (2015), MNRAS 448, 3679
!   https://en.wikipedia.org/wiki/Orbital_elements
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: binaryutils
!
 implicit none
 public :: set_binary,set_multiple,Rochelobe_estimate,L1_point,get_a_from_period
 public :: get_mean_angmom_vector,get_eccentricity_vector

 private
 interface get_eccentricity_vector
  module procedure get_eccentricity_vector,get_eccentricity_vector_sinks
 end interface get_eccentricity_vector

 real, parameter :: pi = 4.*atan(1.)
 integer, parameter :: &
   ierr_m1   = 1, &
   ierr_m2   = 2, &
   ierr_ecc  = 3, &
   ierr_semi = 4, &
   ierr_HIER1 = 5, &
   ierr_HIER2 = 6, &
   ierr_subststar = 7, &
   ierr_Omegasubst = 8, &
   ierr_missstar = 9

contains

!----------------------------------------------------------------
!+
!  setup for a binary
!+
!----------------------------------------------------------------
subroutine set_binary(m1,m2,semimajoraxis,eccentricity, &
                      accretion_radius1,accretion_radius2, &
                      xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                      posang_ascnode,arg_peri,incl,f,verbose)
 use binaryutils, only:get_E
 real,    intent(in)    :: m1,m2
 real,    intent(in)    :: semimajoraxis,eccentricity
 real,    intent(in)    :: accretion_radius1,accretion_radius2
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 integer, intent(out)   :: ierr
 real,    intent(in),  optional :: posang_ascnode,arg_peri,incl,f
 real,    intent(out), optional :: omega_corotate
 logical, intent(in),  optional :: verbose
 integer :: i1,i2,i
 real    :: mtot,dx(3),dv(3),Rochelobe1,Rochelobe2,period
 real    :: x1(3),x2(3),v1(3),v2(3),omega0,cosi,sini,xangle,reducedmass,angmbin
 real    :: a,E,E_dot,P(3),Q(3),omega,big_omega,inc,ecc,tperi,term1,term2,theta
 logical :: do_verbose

 ierr = 0
 do_verbose = .true.
 if (present(verbose)) do_verbose = verbose

 i1 = nptmass + 1
 i2 = nptmass + 2
 nptmass = nptmass + 2

 ! masses
 mtot = m1 + m2

 Rochelobe1 = Rochelobe_estimate(m2,m1,semimajoraxis)
 Rochelobe2 = Rochelobe_estimate(m1,m2,semimajoraxis)
 period = sqrt(4.*pi**2*semimajoraxis**3/mtot)
 reducedmass = m1*m2/mtot
 angmbin = reducedmass*sqrt(mtot*semimajoraxis*(1. - eccentricity**2))

 if (do_verbose) then
    print "(/,2x,a)",'---------- binary parameters ----------- '
    print "(8(2x,a,g12.3,/),2x,a,g12.3)", &
        'primary mass     :',m1, &
        'secondary mass   :',m2, &
        'mass ratio m2/m1 :',m2/m1, &
        'reduced mass     :',reducedmass, &
        'semi-major axis  :',semimajoraxis, &
        'period           :',period, &
        'eccentricity     :',eccentricity, &
        'pericentre       :',semimajoraxis*(1. - eccentricity), &
        'apocentre        :',semimajoraxis*(1. + eccentricity)
 endif
 if (accretion_radius1 > Rochelobe1) then
    print "(1x,a)",'WARNING: set_binary: accretion radius of primary > Roche lobe'
 endif
 if (accretion_radius2 > Rochelobe2) then
    print "(1x,a)",'WARNING: set_binary: accretion radius of secondary > Roche lobe'
 endif
!
!--check for stupid parameter choices
!
 if (m1 <= 0.) then
    print "(1x,a)",'ERROR: set_binary: primary mass <= 0'
    ierr = ierr_m1
 endif
 if (m2 < 0.) then
    print "(1x,a)",'ERROR: set_binary: secondary mass <= 0'
    ierr = ierr_m2
 endif
 if (semimajoraxis <= 0.) then
    print "(1x,a)",'ERROR: set_binary: semi-major axis <= 0'
    ierr = ierr_semi
 endif
 if (eccentricity > 1. .or. eccentricity < 0.) then
    print "(1x,a)",'ERROR: set_binary: eccentricity must be between 0 and 1'
    ierr = ierr_ecc
 endif
 ! exit routine if cannot continue
 if (ierr /= 0) return

 dx = 0.
 dv = 0.
 if (present(posang_ascnode) .and. present(arg_peri) .and. present(incl)) then
    ! Campbell elements
    a = semimajoraxis
    ecc = eccentricity
    omega     = arg_peri*pi/180.
    ! our conventions here are Omega is measured East of North
    big_omega = posang_ascnode*pi/180. + 0.5*pi
    inc       = incl*pi/180.

    if (present(f)) then
       ! get eccentric anomaly from true anomaly
       ! (https://en.wikipedia.org/wiki/Eccentric_anomaly#From_the_true_anomaly)
       theta = f*pi/180.
       E = atan2(sqrt(1. - ecc**2)*sin(theta),(ecc + cos(theta)))
    else
       ! set binary at apastron
       tperi = 0.5*period ! time since periastron: half period = apastron

       ! Solve Kepler equation for eccentric anomaly
       call get_E(period,eccentricity,tperi,E)
    endif

    ! Positions in plane (Thiele-Innes elements)
    P(1) = cos(omega)*cos(big_omega) - sin(omega)*cos(inc)*sin(big_omega)
    P(2) = cos(omega)*sin(big_omega) + sin(omega)*cos(inc)*cos(big_omega)
    P(3) = sin(omega)*sin(inc)
    Q(1) = -sin(omega)*cos(big_omega) - cos(omega)*cos(inc)*sin(big_omega)
    Q(2) = -sin(omega)*sin(big_omega) + cos(omega)*cos(inc)*cos(big_omega)
    Q(3) = sin(inc)*cos(omega)

    term1 = cos(E)-eccentricity
    term2 = sqrt(1.-(eccentricity*eccentricity))*sin(E)
    E_dot = sqrt((m1 + m2)/(a**3))/(1.-eccentricity*cos(E))

    if (do_verbose) then
       print "(4(2x,a,g12.4,/),2x,a,g12.4)", &
             'Eccentric anomaly:',E, &
             'E_dot            :',E_dot, &
             'inclination     (i, deg):',incl, &
             'angle asc. node (O, deg):',posang_ascnode, &
             'arg. pericentre (w, deg):',arg_peri
       if (present(f)) print "(2x,a,g12.4)", &
             'true anomaly    (f, deg):',f
    endif

    ! Rotating everything
    ! Set the positions for the primary and the central secondary
    dx(:) = a*(term1*P(:) + term2*Q(:)) ! + xyzmh_ptmass(1,1)

    ! Set the velocities
    dv(:) = -a*sin(E)*E_dot*P(:) + a*sqrt(1.-(ecc*ecc))*cos(E)*E_dot*Q(:)

 else
    ! set binary at apastron
    dx = (/semimajoraxis*(1. + eccentricity),0.,0./)
    dv = (/0.,sqrt(semimajoraxis*(1.-eccentricity**2)*mtot)/dx(1),0./)
 endif

 ! positions of each star so centre of mass is at zero
 x1 = -dx*m2/mtot
 x2 =  dx*m1/mtot

 ! velocities
 v1 = -dv*m2/mtot !(/0.,-m2/mtot*vmag,0./)
 v2 =  dv*m1/mtot !(/0.,m1/mtot*vmag,0./)

 omega0 = dv(2)/semimajoraxis

 ! print info about positions and velocities
 if (do_verbose) then
    print "(7(2x,a,g12.4,/),2x,a,g12.4)", &
        'angular momentum :',angmbin, &
        'mean ang. speed  :',omega0, &
        'Omega_0 (prim)   :',v1(2)/x1(1), &
        'Omega_0 (second) :',v1(2)/x1(1), &
        'R_accretion (1)  :',accretion_radius1, &
        'R_accretion (2)  :',accretion_radius2, &
        'Roche lobe  (1)  :',Rochelobe1, &
        'Roche lobe  (2)  :',Rochelobe2
 endif

 if (present(omega_corotate)) then
    if (do_verbose) print "(a)",' SETTING VELOCITIES FOR COROTATING FRAME: '
    omega_corotate = omega0
    v1(2) = v1(2) - omega0*x1(1)
    v2(2) = v2(2) - omega0*x2(1)
    if (do_verbose) print "(2(2x,a,g12.4,/))", &
     'Omega_0 (primary)     :',v1(2)/x1(1), &
     'Omega_0 (secondary)   :',v2(2)/x2(1)
 endif

 ! conclude printout
 if (do_verbose) print "(2x,40('-'),/)"

!
!--positions and accretion radii
!
 xyzmh_ptmass(:,i1:i2) = 0.
 xyzmh_ptmass(1:3,i1) = x1
 xyzmh_ptmass(1:3,i2) = x2
 xyzmh_ptmass(4,i1) = m1
 xyzmh_ptmass(4,i2) = m2
 xyzmh_ptmass(5,i1) = accretion_radius1
 xyzmh_ptmass(5,i2) = accretion_radius2
 xyzmh_ptmass(6,i1) = 0.0
 xyzmh_ptmass(6,i2) = 0.0
!
!--velocities
!
 vxyz_ptmass(:,i1) = v1
 vxyz_ptmass(:,i2) = v2
!
! rotate if inclination is non-zero
!
 if (present(incl) .and. .not.(present(arg_peri) .and. present(posang_ascnode))) then
    xangle = incl*pi/180.
    cosi = cos(xangle)
    sini = sin(xangle)
    do i=i1,i2
       call rotate(xyzmh_ptmass(1:3,i),cosi,sini)
       call rotate(vxyz_ptmass(1:3,i),cosi,sini)
    enddo
 endif

end subroutine set_binary

!----------------------------------------------------------------
!+
!  setup for a multiple, using set_binary
!+
!----------------------------------------------------------------
subroutine set_multiple(m1,m2,semimajoraxis,eccentricity, &
                      accretion_radius1,accretion_radius2, &
                      xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                      posang_ascnode,arg_peri,incl,f,verbose,subst)
 real,    intent(in)    :: m1,m2
 real,    intent(in)    :: semimajoraxis,eccentricity
 real,    intent(in)    :: accretion_radius1,accretion_radius2
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 integer, intent(out)   :: ierr
 real,    intent(in),  optional :: posang_ascnode,arg_peri,incl,f
 integer, intent(in),  optional :: subst
 real,    intent(out), optional :: omega_corotate
 logical, intent(in),  optional :: verbose
 integer :: i1,i2,i,subst_index
 real    :: mtot,period
 real    :: x_subst(3),v_subst(3)
 real    :: omega,inc
 !logical :: do_verbose

 real, dimension(24,10) :: data
 character(len=20)      :: hier_prefix
 logical                :: iexist
 integer                :: io, lines
 real                   :: period_ratio,criterion,q_comp,a_comp,e_comp,m_comp
 real                   :: rel_posang_ascnode=0.,rel_arg_peri=0.,rel_incl=0.
 real                   :: q2,mprimary,msecondary
 real                   :: alpha_y, beta_y, gamma_y, alpha_z, beta_z, gamma_z, sign_alpha, sign_gamma

 ierr = 0
 !do_verbose = .true.
 !if (present(verbose)) do_verbose = verbose

 !--- Load/Create HIERARCHY file: xyzmh_ptmass index | hierarchical index | star mass | companion star mass | semi-major axis | eccentricity | period | inclination | argument of pericenter | ascending node longitude
 inquire(file='HIERARCHY', exist=iexist)
 if (present(subst)) then
    if (iexist) then
       open(1, file = 'HIERARCHY', status = 'old')
       lines=0
       do
          read(1, *, iostat=io) data(lines+1,:)
          if (io/=0) exit
          lines = lines + 1
       enddo
       close(1)
    else
       print "(1x,a)",'ERROR: set_multiple: there is no HIERARCHY file, cannot perform subtitution.'
       ierr = ierr_HIER2
    endif
 else
    if (iexist) then
       print "(1x,a)",'WARNING: set_multiple: deleting an existing HIERARCHY file.'
       open(1, file='HIERARCHY', status='old')
       close(1, status='delete')
    endif

    mtot = m1 + m2
    period = sqrt(4.*pi**2*semimajoraxis**3/mtot)

    open(1, file = 'HIERARCHY', status = 'new')
    if (present(incl)) then
       if (present(posang_ascnode) .and. present(arg_peri)) then
          write(1,*) 1, 11, m1, m2, semimajoraxis, eccentricity, period, incl, arg_peri, posang_ascnode
          write(1,*) 2, 12, m2, m1, semimajoraxis, eccentricity, period, incl, arg_peri, posang_ascnode
       else ! set binary at apastron with inclination
          write(1,*) 1, 11, m1, m2, semimajoraxis, eccentricity, period, incl, 0, 0
          write(1,*) 2, 12, m2, m1, semimajoraxis, eccentricity, period, incl, 0, 0
       endif
    else ! set binary at apastron without inclination
       write(1,*) 1, 11, m1, m2, semimajoraxis, eccentricity, period, 0, 0, 0
       write(1,*) 2, 12, m2, m1, semimajoraxis, eccentricity, period, 0, 0, 0
    endif
    close(1)
 endif

 !--- Checks to avoid bad substitutions
 if (present(subst)) then
    write(hier_prefix, *) subst
    io=0
    subst_index = 0
    mtot = 0.
    do i=1,lines
       if (data(i,2)==abs(subst)) then ! Check that star to be substituted exists in HIERARCHY file
          if (data(i,1)==0) then ! Check that star to be substituted has not already been substituted
             print "(1x,a)",'ERROR: set_multiple: star '//trim(hier_prefix)//' substituted yet.'
             ierr = ierr_subststar
          endif
          subst_index = int(data(i,1))
          data(i,1) = 0

          if (subst>0) then
             rel_posang_ascnode = data(i, 10)

             if (rel_posang_ascnode /= 0) then
                print "(1x,a)",'ERROR: set_multiple: at the moment phantom can subst only Omega=0 binaries.'
                ierr = ierr_Omegasubst
             endif

             rel_arg_peri= data(i, 9)
             rel_incl = data(i, 8)
          else
             rel_posang_ascnode = posang_ascnode
             rel_arg_peri = arg_peri
             rel_incl = incl
          endif

          mtot = data(i, 3)
          m_comp = data(i, 4)
          a_comp = data(i, 5)
          e_comp = data(i, 6)

          q_comp = mtot/m_comp
          if (q_comp>1) q_comp=q_comp**(-1)

          ! Mardling&Aarseth (2001) criterion check

          period_ratio = sqrt((a_comp*a_comp*a_comp)/(m_comp+mtot)/(semimajoraxis*semimajoraxis*semimajoraxis)*(mtot)) ! Po/Pi
          criterion = 4.7*(1-e_comp)**(-1.8)*(1+e_comp)**(0.6)*(1+q_comp)**(0.1)

          if (criterion > period_ratio) then
             print "(1x,a)",'WARNING: set_multiple: orbital parameters does not satisfy Mardling and Aarseth stability criterion.'
          endif

          q2=m2/m1
          mprimary = mtot/(1+q2)
          msecondary = mtot*q2/(1+q2)

          io=1
          exit
       endif
    enddo

    if (io == 0) then
       print "(1x,a)",'ERROR: set_multiple: star '//trim(hier_prefix)//' not present in HIERARCHY file.'
       ierr = ierr_missstar
    endif

    if (subst_index > 0 .and. subst_index <= size(xyzmh_ptmass(1,:))) then ! check for seg fault
       x_subst(:)=xyzmh_ptmass(1:3,subst_index)
       v_subst(:)=vxyz_ptmass(:,subst_index)
    endif
    !i1 = subst_index
    !i2 = nptmass + 1
    !nptmass = nptmass + 1

    period = sqrt(4.*pi**2*semimajoraxis**3/mtot)
 else
    mprimary = m1
    msecondary = m2

    if (present(posang_ascnode)) rel_posang_ascnode = posang_ascnode
    if (present(arg_peri)) rel_arg_peri= arg_peri
    if (present(incl)) rel_incl = incl

 endif

 !--- Create the binary
 call set_binary(mprimary,msecondary,semimajoraxis=semimajoraxis,eccentricity=eccentricity, &
            posang_ascnode=rel_posang_ascnode,arg_peri=rel_arg_peri,incl=rel_incl, &
            f=f,accretion_radius1=accretion_radius1,accretion_radius2=accretion_radius2, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass, ierr=ierr)

 if (present(subst)) then
    !--- lower nptmass, copy one of the new sinks to the subst star
    nptmass = nptmass-1
    i1 = subst_index
    i2 = nptmass

    ! positions and accretion radii
    xyzmh_ptmass(1:6,i1) = xyzmh_ptmass(1:6,nptmass+1)

    ! test Jolien
!    print "(5(2x,a,g12.3,/),2x,a,g12.3)", &
!    'i1     :',i1, &
!     'mass i1:',xyzmh_ptmass(4,i1), &
!     'i2     :',i2, &
!     'mass i2:',xyzmh_ptmass(4,i2)

    ! velocities
    vxyz_ptmass(:,i1) = vxyz_ptmass(:,nptmass+1)

    !---
    ! Rotate the substituting binary with orientational parameters
    ! referring to the substituted star's orbital plane
    if (subst>0) then

       omega     = rel_arg_peri*pi/180.
       !big_omega = rel_posang_ascnode*pi/180.! + 0.5*pi
       inc       = rel_incl*pi/180.

       ! Retrieve eulerian angles of the substituted star orbit's semi-major axis (y axis)
       if (omega <= pi/2) then
          beta_y = omega
          sign_alpha=-1
          if (inc <= pi) then
             sign_gamma=1
          else
             sign_gamma=-1
          endif
       else
          beta_y = 2*pi-omega
          sign_alpha=1
          if (inc <= pi) then
             sign_gamma=-1
          else
             sign_gamma=1
          endif
       endif
       gamma_y=acos(sign_gamma*sin(beta_y)*sin(inc))
       alpha_y=acos(sign_alpha*sqrt(abs(sin(beta_y)**2-cos(gamma_y)**2))) ! Needs abs cause float approx for cos

       ! Retrieve eulerian angles of the axis perpendicular to the substituted star orbital plane (z axis)
       beta_z = pi/2.
       gamma_z = inc
       alpha_z = pi/2. - inc
       if (inc <= pi) then
          gamma_z=inc
          if (inc <= pi/2.) then
             alpha_z = pi/2.-inc
          elseif (inc > pi/2.) then
             alpha_z = inc-pi/2.
          endif
       elseif (inc < 2.*pi .and. inc > pi) then
          gamma_z = 2.*pi-inc
          if (inc <= 3.*pi/2.) then
             alpha_z = inc-pi/2
          elseif (inc > 3.*pi/2.) then
             alpha_z = 5.*pi/2.-inc
          endif
       endif

       ! Rotate substituting sinks by argument of pericenter around the z axis
       call gen_rotate(xyzmh_ptmass(1:3,i1),alpha_z,beta_z,gamma_z, arg_peri*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i1),alpha_z,beta_z,gamma_z, arg_peri*pi/180)
       call gen_rotate(xyzmh_ptmass(1:3,i2),alpha_z,beta_z,gamma_z, arg_peri*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i2),alpha_z,beta_z,gamma_z, arg_peri*pi/180)

       ! Rotate substituting sinks by inclination around the y axis
       call gen_rotate(xyzmh_ptmass(1:3,i1),alpha_y,beta_y,gamma_y, incl*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i1),alpha_y,beta_y,gamma_y, incl*pi/180)
       call gen_rotate(xyzmh_ptmass(1:3,i2),alpha_y,beta_y,gamma_y, incl*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i2),alpha_y,beta_y,gamma_y, incl*pi/180)

       ! Rotate substituting sinks by ascending node longitude around the z axis
       call gen_rotate(xyzmh_ptmass(1:3,i1),alpha_z,beta_z,gamma_z, posang_ascnode*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i1),alpha_z,beta_z,gamma_z, posang_ascnode*pi/180)
       call gen_rotate(xyzmh_ptmass(1:3,i2),alpha_z,beta_z,gamma_z, posang_ascnode*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i2),alpha_z,beta_z,gamma_z, posang_ascnode*pi/180)
    endif

    ! Move the substituting binary's center of mass in the substituted star position
    xyzmh_ptmass(1:3,i1) = xyzmh_ptmass(1:3,i1)+x_subst
    xyzmh_ptmass(1:3,i2) = xyzmh_ptmass(1:3,i2)+x_subst
    ! Set the substituting binary's center of mass velocity
    vxyz_ptmass(:,i1) = vxyz_ptmass(:,i1)+v_subst
    vxyz_ptmass(:,i2) = vxyz_ptmass(:,i2)+v_subst

    ! Write updated HIERARCHY file with the two new stars and the substituted one
    open(1, file = 'HIERARCHY', status = 'old')
    do i=1,lines
       write(1,*) int(data(i,1)), int(data(i,2)), data(i,3:)
    enddo
    write(1,*) i1, trim(hier_prefix)//"1", mprimary, msecondary, semimajoraxis, eccentricity, &
         period, incl, arg_peri, posang_ascnode
    write(1,*) i2, trim(hier_prefix)//"2", msecondary, mprimary, semimajoraxis, eccentricity, &
         period, incl, arg_peri, posang_ascnode
    close(1)
 endif

end subroutine set_multiple

!------------------------------------
! Rotate an (x,y,z) point by theta
! radiants around an axis with alpha,
! beta and gamma eulerian angles
!------------------------------------
pure subroutine gen_rotate(xyz,alpha,beta,gamma,theta)
 real, intent(inout) :: xyz(3)
 real, intent(in)    :: alpha, beta, gamma,theta
 real :: xi,yi,zi,A,B,C,D,E,F,G,H,I,nx,ny,nz

 nx=cos(alpha)
 ny=cos(beta)
 nz=cos(gamma)

 A=cos(theta)+nx**2*(1-cos(theta))
 B=nx*ny*(1-cos(theta))-nz*sin(theta)
 C=nx*nz*(1-cos(theta))+ny*sin(theta)
 D=nx*ny*(1-cos(theta))+nz*sin(theta)
 E=cos(theta)+ny**2*(1-cos(theta))
 F=ny*nz*(1-cos(theta))-nx*sin(theta)
 G=nx*nz*(1-cos(theta))-ny*sin(theta)
 H=ny*nz*(1-cos(theta))+nx*sin(theta)
 I=cos(theta)+nz**2*(1-cos(theta))

 xi = xyz(1)
 yi = xyz(2)
 zi = xyz(3)
 xyz(1) = A*xi+B*yi+C*zi
 xyz(2) = D*xi+E*yi+F*zi
 xyz(3) = G*xi+H*yi+I*zi

end subroutine gen_rotate

pure subroutine rotate(xyz,cosi,sini)
 real, intent(inout) :: xyz(3)
 real, intent(in)    :: cosi,sini
 real :: xi,yi,zi

 xi = xyz(1)
 yi = xyz(2)
 zi = xyz(3)
 xyz(1) =  xi*cosi + zi*sini
 xyz(2) =  yi
 xyz(3) = -xi*sini + zi*cosi

end subroutine rotate

!------------------------------------
! Compute estimate of the Roche Lobe
! Eggleton (1983) ApJ 268, 368-369
!------------------------------------
real function Rochelobe_estimate(m1,m2,sep)
 real, intent(in) :: m1,m2,sep
 real :: q,q13,q23

 q = m2/m1
 q13 = q**(1./3.)
 q23 = q13*q13
 Rochelobe_estimate = sep * 0.49*q23/(0.6*q23 + log(1. + q13))

end function Rochelobe_estimate

!---------------------------------------------
! Find first Lagrange point (L1)
! via Newton-Raphson solution of quintic
!
! INPUT: mass ratio of binary
! OUTPUT: L1 point, as distance from primary
!---------------------------------------------
real function L1_point(qinv)
 real, intent(in) :: qinv
 real :: fL, dfL, dL, L, q11

 q11 = 1./(1.+qinv)
 L = 0.5 + 0.2222222*log10(qinv)

 dL = 1.e7
 do while (abs(dL)>1.e-6)
    fL = qinv/L**2- 1./(1.-L)**2 - (1.+qinv)*L + 1.
    dfL=-2*qinv/L**3 - 2./(1.-L)**3 - (1.+qinv)
    dL = -fL/(dfL*L)
    L = L*(1.+dL)
 enddo

 L1_point = L

end function L1_point

!-------------------------------------------------------------
! Function to determine the semi-major axis given the period
!-------------------------------------------------------------
function get_a_from_period(m1,m2,period) result(a)
 real, intent(in) :: m1,m2,period
 real :: a

 a = ((m1 + m2)*(period/(2.*pi))**2)**(1./3.)

end function get_a_from_period

!----------------------------------------------------
! Eccentricity vector, for second body
! https://en.wikipedia.org/wiki/Eccentricity_vector
!----------------------------------------------------
function get_eccentricity_vector(m1,m2,x1,x2,v1,v2)
 real, intent(in) :: m1,m2
 real, intent(in) :: x1(3),x2(3),v1(3),v2(3)
 real :: x0(3),v0(3),r(3),v(3),dr,mu
 real :: get_eccentricity_vector(3)

 ! centre of mass position and velocity
 x0 = (m1*x1 + m2*x2)/(m1 + m2)
 v0 = (m1*v1 + m2*v2)/(m1 + m2)

 ! position and velocity vectors relative to each other
 r = x2 - x1
 v = v2 - v1

 ! intermediate quantities
 dr = 1./sqrt(dot_product(r,r))
 mu = m1 + m2  ! "standard gravitational parameter"

 ! formula for eccentricity vector
 get_eccentricity_vector = (dot_product(v,v)/mu - dr)*r - dot_product(r,v)/mu*v

end function get_eccentricity_vector

!----------------------------------------------------
! interface to above assuming two sink particles
!----------------------------------------------------
function get_eccentricity_vector_sinks(xyzmh_ptmass,vxyz_ptmass,i1,i2)
 real,    intent(in) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: i1, i2
 real :: get_eccentricity_vector_sinks(3)

 if (i1 > 0 .and. i2 > 0) then
    get_eccentricity_vector_sinks = get_eccentricity_vector(&
        xyzmh_ptmass(4,i1),xyzmh_ptmass(4,i2),&
        xyzmh_ptmass(1:3,i1),xyzmh_ptmass(1:3,i2),&
        vxyz_ptmass(1:3,i1),vxyz_ptmass(1:3,i2))
 else
    get_eccentricity_vector_sinks = 0.
 endif

end function get_eccentricity_vector_sinks

!-------------------------------------------------------------
! Function to find mean angular momentum vector from a list
! of positions and velocities
!-------------------------------------------------------------
function get_mean_angmom_vector(n,xyz,vxyz) result(l)
 integer, intent(in) :: n
 real,    intent(in) :: xyz(:,:),vxyz(:,:)
 real    :: l(3),li(3)
 integer :: i

 l = 0.
 do i=1,n
    call get_cross_product(xyz(:,i),vxyz(:,i),li)
    l = l + li
 enddo
 l = l/real(n)

end function get_mean_angmom_vector

!-------------------------------------------------------------
!
! cross product routine
!
!-------------------------------------------------------------
pure subroutine get_cross_product(veca,vecb,vecc)
 real, intent(in)  :: veca(3),vecb(3)
 real, intent(out) :: vecc(3)

 vecc(1) = veca(2)*vecb(3) - veca(3)*vecb(2)
 vecc(2) = veca(3)*vecb(1) - veca(1)*vecb(3)
 vecc(3) = veca(1)*vecb(2) - veca(2)*vecb(1)

end subroutine get_cross_product

end module setbinary
