!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testorbits
!
! Unit tests for routines in the orbits module (utils_orbits.f90)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: orbits, part, setbinary, sethierarchical, setorbit,
!   testutils
!
 use testutils, only:checkval,update_test_scores
 use orbits,    only:get_eccentricity,get_semimajor_axis, get_orbital_period, &
                     get_E_from_true_anomaly,get_E_from_mean_anomaly, &
                     get_true_anomaly_from_separation,get_time_to_separation, &
                     convert_flyby_to_elements,escape,pi,rad_to_deg, &
                     get_orbital_elements,get_pericentre_distance,&
                     get_inclination,get_true_anomaly,get_dx_dv_ptmass,orbit_is_parabolic
 use setbinary, only:set_binary
 implicit none
 real, parameter :: tol = 2.e-14

 public :: test_orbits

 private

contains

!--------------------------------------------
!+
!  Run all orbital utility tests
!+
!--------------------------------------------
subroutine test_orbits(ntests,npass)
 integer, intent(inout) :: ntests,npass

 call test_circular_orbit(ntests,npass)
 call test_true_anomaly_from_separation(ntests,npass)
 call test_time_to_separation_ellipse(ntests,npass)
 call test_convert_flyby(ntests,npass)
 call test_escape_velocity(ntests,npass)
 call test_set_binary_simple(ntests,npass)
 call test_set_binary_full_elements(ntests,npass)
 call test_set_orbits_nstar(ntests,npass)

end subroutine test_orbits

!--------------------------------------------
!+
!  Test e, a, P on a circular Keplerian orbit
!+
!--------------------------------------------
subroutine test_circular_orbit(ntests,npass)
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(5)
 real :: mu,a,P
 real :: dx(3),dv(3)

 if (npass < 0) return
 nfailed = 0

 mu = 1.0
 a  = 2.5
 P  = 2.*pi*sqrt(a**3/mu)

 dx = 0.; dv = 0.
 dx(1) = a
 dv(2) = sqrt(mu/a)

 call checkval(get_eccentricity(mu,dx,dv),0.0,tol,nfailed(1),'eccentricity (circular)')
 call checkval(get_semimajor_axis(mu,dx,dv),a,tol,nfailed(2),'semimajor axis (circular)')
 call checkval(get_orbital_period(mu,a),P,tol,nfailed(3),'orbital period (circular)')
 call checkval(get_pericentre_distance(mu,dx,dv),a,tol,nfailed(4),'pericentre (circular)')
 call checkval(get_true_anomaly(mu,dx,dv),0.,tol,nfailed(5),'true anomaly (circular)')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_circular_orbit

!--------------------------------------------
!+
!  Test true anomaly from radius for an ellipse
!+
!--------------------------------------------
subroutine test_true_anomaly_from_separation(ntests,npass)
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(5)
 real :: a,e,rp,ra,f0,fp,fa,fp1,fp2

 nfailed = 0
 a  = 1.0
 e  = 0.5
 rp = a*(1.0 - e)
 ra = a*(1.0 + e)

 f0 = get_true_anomaly_from_separation(a,e,a)  ! r=a occurs at f = acos((a(1-e^2)/a -1)/e)
 fp = get_true_anomaly_from_separation(a,e,rp) ! pericentre -> f=0
 fa = get_true_anomaly_from_separation(a,e,ra) ! apocentre -> f=180.0
 fp1 = get_true_anomaly_from_separation(a,1.0,a) ! pericentre of parabolic orbit, f=0
 fp2 = get_true_anomaly_from_separation(a,0.0,a) ! circular orbit, ambiguous but should return f=180.

 call checkval(fp,0.0,tol,nfailed(1),'true anomaly at pericentre')
 call checkval(fa,180.0,tol,nfailed(2),'true anomaly at apocentre')
 call checkval(f0,acos((a*(1.-e**2)/a - 1.)/e)*rad_to_deg,tol,nfailed(3),'true anomaly at r=a')
 call checkval(fp1,0.0,tol,nfailed(4),'true anomaly at rp (e=1)')
 call checkval(fp2,180.0,tol,nfailed(5),'true anomaly at rp (e=0)')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_true_anomaly_from_separation

!--------------------------------------------
!+
!  Test time to go from pericentre to apocentre = P/2 for ellipse
!+
!--------------------------------------------
subroutine test_time_to_separation_ellipse(ntests,npass)
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(1)
 real :: mu,a,e,rp,ra,time,P

 nfailed = 0
 mu = 1.0
 a  = 1.0
 e  = 0.5
 rp = a*(1.0 - e)
 ra = a*(1.0 + e)
 P  = 2.*pi*sqrt(a**3/mu)

 time = get_time_to_separation(mu,a,e,rp,ra)
 call checkval(time,P*0.5,tol,nfailed(1),'time peri->apo equals P/2')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_time_to_separation_ellipse

!--------------------------------------------
!+
!  Test flyby conversion helper (basic sanity)
!+
!--------------------------------------------
subroutine test_convert_flyby(ntests,npass)
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(6)
 real :: rp,e,d,a,f

 nfailed = 0

 ! elliptical case (here distance > apocentre so should get f = 180.0)
 rp = 2.0; e = 0.6; d = 10.0
 call convert_flyby_to_elements(rp,d,e,a,f)
 call checkval(a,rp/(1.-e),tol,nfailed(1),'a from flyby (e<1)')
 call checkval(f,-180.0,tol,nfailed(2),'f from flyby (e<1)')

 ! parabolic case
 e = 1.0
 call convert_flyby_to_elements(rp,d,e,a,f)
 call checkval(a,rp,tol,nfailed(3),'a from flyby (e=1)')

 ! hyperbolic case
 e = 1.5
 call convert_flyby_to_elements(rp,d,e,a,f)
 call checkval(a,rp/(1.-e),tol,nfailed(4),'a from flyby (e>1)')

 ! hyperbolic with requested distance < pericentre
 e = 10.0; d = 1.5
 call convert_flyby_to_elements(rp,d,e,a,f)
 call checkval(a,rp/(1.-e),tol,nfailed(5),'a from flyby (e>1)')
 call checkval(f,0.0,tol,nfailed(6),'f from flyby (e>1)')

 call update_test_scores(ntests,nfailed,npass)

end subroutine test_convert_flyby

!--------------------------------------------
!+
!  Test escape velocity condition
!+
!--------------------------------------------
subroutine test_escape_velocity(ntests,npass)
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2)
 real :: mu,r,vesc

 nfailed = 0

 mu = 1.0
 r  = 2.0
 vesc = sqrt(2.*mu/r)

 call checkval(escape(1.1*vesc,mu,r),.true., nfailed(1),'escape when v>vesc')
 call checkval(escape(0.9*vesc,mu,r),.false.,nfailed(2),'bound when v<vesc')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_escape_velocity

!--------------------------------------------
!+
!  Test set_binary with minimal elements (a,e,i only)
!+
!--------------------------------------------
subroutine test_set_binary_simple(ntests,npass)
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(3),nptmass,ierr
 real :: m1,m2,mu,a,e,inc
 real :: xyzmh(6,2),vxyz(3,2),dx(3),dv(3)

 nfailed = 0
 m1 = 2.0; m2 = 3.0; mu = m1+m2
 a  = 2.0
 e  = 0.3
 inc = 20.0
 xyzmh = 0.; vxyz = 0.; nptmass = 0

 call set_binary(m1,m2,a,e,0.,0.,xyzmh,vxyz,nptmass,ierr,verbose=.false.,incl=inc)
 call get_dx_dv_ptmass(xyzmh,vxyz,dx,dv)

 call checkval(get_eccentricity(mu,dx,dv),e,tol,nfailed(1),'ecc (simple)')
 call checkval(get_semimajor_axis(mu,dx,dv),a,tol,nfailed(2),'a (simple)')
 call checkval(get_inclination(dx,dv),inc,tol,nfailed(3),'inc (simple)')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_set_binary_simple

!--------------------------------------------
!+
!  Test set_binary with full elements (elliptic, hyperbolic, parabolic)
!+
!--------------------------------------------
subroutine test_set_binary_full_elements(ntests,npass)
 integer, intent(inout) :: ntests,npass
 ! Test parameters for different orbit types
 character(len=*), parameter :: orbit_types(3) = (/'elliptic  ', 'hyperbolic', 'parabolic '/)
 real, parameter :: a_vals(3) = (/3.0, -5.0, 3.0/)  ! semi-major axis (parabolic uses rp)
 real, parameter :: e_vals(3) = (/0.5, 1.4, 1.0/)   ! eccentricity
 real, parameter :: inc_vals(3) = (/30.0, 15.0, 5.0/)  ! inclination
 real, parameter :: O_vals(3) = (/-10.0, -25.0, 20.0/) ! longitude of ascending node
 real, parameter :: w_vals(3) = (/-40.0, 60.0, -31.0/) ! argument of periapsis
 real, parameter :: f_vals(3) = (/120.0, 30.0, 120.0/)  ! true anomaly
 integer :: nfailed(10),nptmass,ierr,i
 real :: m1,m2,mu,a_rec,e_rec,i_rec,O_rec,w_rec,f_rec
 real :: xyzmh(6,2),vxyz(3,2),dx(3),dv(3),rp_vals(3)
 real :: dx0(3),dv0(3)
 real :: a0,e0,i0,O0,w0,f0

 m1 = 3.0; m2 = 2.0; mu = m1+m2

 do i = 1, 3
    nfailed = 0
    xyzmh = 0.; vxyz = 0.; nptmass = 0

    ! check for pericentre distance
    if (orbit_is_parabolic(e_vals(i))) then
       rp_vals(i) = a_vals(i) ! for parabolic orbits, use rp = a
    else
       rp_vals(i) = a_vals(i) * (1.0 - e_vals(i))
    endif

    call set_binary(m1,m2,a_vals(i),e_vals(i),0.,0.,xyzmh,vxyz,nptmass,ierr,&
                    posang_ascnode=O_vals(i),arg_peri=w_vals(i),incl=inc_vals(i),f=f_vals(i),verbose=.false.)
    call get_dx_dv_ptmass(xyzmh,vxyz,dx,dv)

    call get_orbital_elements(mu,dx,dv,a_rec,e_rec,i_rec,O_rec,w_rec,f_rec)
    call checkval(a_rec,a_vals(i),tol,nfailed(1),'a ('//trim(orbit_types(i))//')')
    call checkval(e_rec,e_vals(i),tol,nfailed(2),'e ('//trim(orbit_types(i))//')')
    call checkval(i_rec,inc_vals(i),tol,nfailed(3),'inc ('//trim(orbit_types(i))//')')
    call checkval(O_rec,O_vals(i),tol,nfailed(4),'Omega ('//trim(orbit_types(i))//')')
    call checkval(w_rec,w_vals(i),tol,nfailed(5),'w ('//trim(orbit_types(i))//')')
    call checkval(f_rec,f_vals(i),tol,nfailed(6),'f ('//trim(orbit_types(i))//')')

    ! also check individual extraction functions
    call checkval(get_semimajor_axis(mu,dx,dv),a_vals(i),tol,nfailed(7),'get_a ('//trim(orbit_types(i))//')')
    call checkval(get_eccentricity(mu,dx,dv),e_vals(i),tol,nfailed(8),'get_ecc ('//trim(orbit_types(i))//')')
    call checkval(get_inclination(dx,dv),inc_vals(i),tol,nfailed(9),'get_inc ('//trim(orbit_types(i))//')')
    call checkval(get_pericentre_distance(mu,dx,dv),rp_vals(i),tol,nfailed(10),'get_peri ('//trim(orbit_types(i))//')')

    call update_test_scores(ntests,nfailed,npass)
 enddo

 ! also try setting up a binary manually
 nptmass = 0
 xyzmh(1:3,1) = (/-1.2,1.0,0.0/)
 xyzmh(1:3,2) = (/0.8,0.0,0.0/)
 vxyz(1:3,1) = (/1.0,0.0,0.0/)
 vxyz(1:3,2) = (/0.0,1.0,0.0/)
 ! query the orbital elements
 call get_dx_dv_ptmass(xyzmh,vxyz,dx0,dv0)
 call get_orbital_elements(mu,dx0,dv0,a_rec,e_rec,i_rec,O_rec,w_rec,f_rec)

 ! set up a binary with these elements
 call set_binary(m1,m2,a_rec,e_rec,0.,0.,xyzmh,vxyz,nptmass,ierr,&
                 posang_ascnode=O_rec,arg_peri=w_rec,incl=i_rec,f=f_rec,verbose=.true.)

 ! check that the original separation and velocity difference are recovered
 call get_dx_dv_ptmass(xyzmh,vxyz,dx,dv)
 nfailed = 0
 do i=1,3
    call checkval(dx(i),dx0(i),tol,nfailed(i),'dx (manual)')
 enddo
 do i=1,3
    call checkval(dv(i),dv0(i),tol,nfailed(3+i),'dv (manual)')
 enddo
 call update_test_scores(ntests,nfailed,npass)

 call get_orbital_elements(mu,dx,dv,a0,e0,i0,O0,w0,f0)
 call checkval(a0,a_rec,tol,nfailed(1),'a (manual)')
 call checkval(e0,e_rec,tol,nfailed(2),'e (manual)')
 call checkval(i0,i_rec,tol,nfailed(3),'inc (manual)')
 call checkval(O0,O_rec,tol,nfailed(4),'Omega (manual)')
 call checkval(w0,w_rec,tol,nfailed(5),'w (manual)')
 call checkval(f0,f_rec,tol,nfailed(6),'f (manual)')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_set_binary_full_elements

!--------------------------------------------
!+
!  Test set_orbits for 1, 2 and several stars: star i must end
!  up on sink i with its own mass and accretion radius, and the
!  system must have zero total momentum
!+
!--------------------------------------------
subroutine test_set_orbits_nstar(ntests,npass)
 use setorbit,        only:orbit_t,set_defaults_orbits,set_orbits
 use sethierarchical, only:hierarchy,set_hierarchical_default_options
 use sethier_utils,   only:load_hierarchy_file,hier_db_size,hier_db_prop
 use part,            only:nsinkproperties
 use io,              only:id
 integer, intent(inout) :: ntests,npass
 integer, parameter :: maxstar = 4
 character(len=20) :: prefix
 real, parameter :: m(maxstar)    = (/1.0, 0.5, 0.3, 0.2/)
 real, parameter :: hacc(maxstar) = (/0.01, 0.02, 0.03, 0.04/)
 type(orbit_t) :: orbit
 integer :: nfailed(5),nstar,nptmass,ierr,i,iu,lines
 real :: xyzmh(nsinkproperties,maxstar+1),vxyz(3,maxstar+1),mom(3)
 real :: data(hier_db_size,hier_db_prop)
 integer :: isink(hier_db_size)
 real :: msink(maxstar),hsink(maxstar)
 logical :: iexist

 ! one hierarchy file per MPI task, since every task runs this test
 write(prefix,"(a,i0)") 'test_orbits',id
 do nstar=1,maxstar
    nfailed = 0
    xyzmh = 0.; vxyz = 0.; nptmass = 0
    call set_defaults_orbits(nstar,orbit)
    call set_orbits(nstar,m,hacc,orbit,prefix,xyzmh,vxyz,nptmass,.false.,ierr)
    call checkval(ierr,0,0,nfailed(1),'set_orbits ierr')
    call checkval(nptmass,nstar,0,nfailed(2),'number of sinks')
    ! copy to contiguous arrays to avoid array temporaries in the calls below
    msink = xyzmh(4,1:maxstar)
    hsink = xyzmh(5,1:maxstar)
    call checkval(nstar,msink,m,0.,nfailed(3),'star i has the mass of star i')
    call checkval(nstar,hsink,hacc,0.,nfailed(4),'star i has the accretion radius of star i')
    do i=1,3
       mom(i) = sum(xyzmh(4,1:nstar)*vxyz(i,1:nstar))
    enddo
    call checkval(sqrt(dot_product(mom,mom)),0.,tol,nfailed(5),'total momentum is zero')
    call update_test_scores(ntests,nfailed,npass)
 enddo

 ! a binary with a planet around the first star, given so that the stars are sinks 1 and 2
 ! and the planet is sink 3: the sinks must follow the order of the hierarchy string
 nfailed = 0
 nstar = 3
 xyzmh = 0.; vxyz = 0.; nptmass = 0
 call set_defaults_orbits(nstar,orbit)
 hierarchy = '111,12,112'
 call set_hierarchical_default_options()
 call set_orbits(nstar,m,hacc,orbit,prefix,xyzmh,vxyz,nptmass,.false.,ierr)
 call checkval(ierr,0,0,nfailed(1),'set_orbits ierr (111,12,112)')
 call checkval(nptmass,nstar,0,nfailed(2),'number of sinks (111,12,112)')
 msink = xyzmh(4,1:maxstar)
 hsink = xyzmh(5,1:maxstar)
 call checkval(nstar,msink,m,0.,nfailed(3),'sinks in order of hierarchy string (mass)')
 call checkval(nstar,hsink,hacc,0.,nfailed(4),'sinks in order of hierarchy string (racc)')
 do i=1,3
    mom(i) = sum(xyzmh(4,1:nstar)*vxyz(i,1:nstar))
 enddo
 call checkval(sqrt(dot_product(mom,mom)),0.,tol,nfailed(5),'total momentum is zero (111,12,112)')
 call update_test_scores(ntests,nfailed,npass)

 ! a triple appended after an existing sink: the sinks of the triple are sinks 2 to 4,
 ! and the .hierarchy file must give these indices
 nfailed = 0
 nstar = 3
 xyzmh = 0.; vxyz = 0.; nptmass = 1
 xyzmh(4,1) = 2.
 call set_defaults_orbits(nstar,orbit)
 call set_orbits(nstar,m,hacc,orbit,prefix,xyzmh,vxyz,nptmass,.false.,ierr)
 call checkval(ierr,0,0,nfailed(1),'set_orbits ierr (appended)')
 call checkval(nptmass,nstar+1,0,nfailed(2),'number of sinks (appended)')
 msink = xyzmh(4,2:maxstar+1)
 call checkval(nstar,msink,m,0.,nfailed(3),'appended sinks after the existing sink')
 call load_hierarchy_file(trim(prefix)//'.hierarchy',data,lines,ierr)
 isink = 0
 isink(1:lines) = int(data(1:lines,1))
 call checkval(minval(isink(1:lines),mask=isink(1:lines) > 0),2,0,nfailed(4),'first sink index in .hierarchy file')
 call checkval(maxval(isink(1:lines)),nstar+1,0,nfailed(5),'last sink index in .hierarchy file')
 call update_test_scores(ntests,nfailed,npass)

 ! remove the file written by the hierarchical setup
 inquire(file=trim(prefix)//'.hierarchy',exist=iexist)
 if (iexist) then
    open(newunit=iu,file=trim(prefix)//'.hierarchy',status='old')
    close(iu,status='delete')
 endif

end subroutine test_set_orbits_nstar

end module testorbits

