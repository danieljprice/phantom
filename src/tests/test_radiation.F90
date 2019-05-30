module testradiation
 use part, only:ithick,iradxi,ifluxx,ifluxy,ifluxz,idflux,ikappa

 implicit none

 public :: test_radiation
 private

contains

subroutine test_radiation(ntests,npass)
 use dim,          only:maxp,exchange_radiation_energy
 use io,           only:id,master,iverbose
 use part,         only:npart,xyzh,fxyzu,vxyzu,massoftype,igas,&
                        iphase,maxphase,isetphase,rhoh,&
                        npartoftype,&
                        radiation,iradxi,ikappa,idflux,ifluxx,ifluxy,ifluxz,&
                        maxvxyzu
 use kernel,       only:hfact_default
 use unifdis,      only:set_unifdis
 use units,        only:set_units,unit_density
 use eos,          only:gmw,gamma,polyk
 use linklist,     only:set_linklist
 use densityforce, only:densityiterate
 use forces,       only:force
 use boundary,     only:dxbound,dybound,dzbound
 use physcon,      only:au,solarm

 integer, intent(inout) :: ntests,npass

 real :: psep,hfact
 real :: pmassi,rhozero,totmass

 if (id==master) write(*,"(/,a,/)") '--> TESTING RADIATION MODULE'

 iverbose = 1
 exchange_radiation_energy = .false.

 radiation(:,:) = 0.
 radiation(ithick,:) = 1.

 psep = 1./16.
 hfact = hfact_default
 call set_units(dist=au,mass=solarm,G=1.d0)
 npart = 0

 call set_unifdis('cubic',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,psep,hfact,npart,xyzh)
 rhozero = 1.e-7/unit_density  ! 1e-7 g/cm^3
 totmass = rhozero*(dxbound*dybound*dzbound)
 massoftype(igas) = totmass/npart
 gamma = 5./3.
 gmw = 2.0
 polyk = 0.
 if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(1) = npart
 pmassi = massoftype(igas)

 call test_exchange_terms(npart,pmassi,radiation,xyzh,vxyzu,fxyzu,ntests,npass)

#ifndef PERIODIC
  if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF RADIATION MODULE (need -DPERIODIC)'
#else

 call test_uniform_derivs(ntests,npass)
#endif

 if (id==master) write(*,"(/,a)") '<-- RADIATION TEST COMPLETE'
end subroutine test_radiation

subroutine test_exchange_terms(npart,pmassi,radiation,xyzh,vxyzu,fxyzu,ntests,npass)
 use radiation_utils, only:update_radenergy
 use units,      only:unit_ergg,unit_density,umass,udist,utime
 use part,       only:rhoh
 use testutils,  only:checkval
 use physcon,    only:seconds

 implicit none

 real,intent(inout) :: &
   pmassi
 real,intent(inout) :: &
   radiation(:,:),xyzh(:,:),vxyzu(:,:),fxyzu(:,:)
 integer,intent(inout) ::&
   npart,ntests,npass

 real :: dt,t,physrho,rhoi,maxt,laste
 integer :: i,ierr

  do i=1,npart
     rhoi         = rhoh(xyzh(4,i),pmassi)
     radiation(iradxi,i) = 1e12/(unit_ergg*unit_density)/rhoi
     radiation(ikappa,i)  = 0.4/(udist**2/umass)
     vxyzu(4,i)   = 1e10/(unit_ergg*unit_density)
     vxyzu(4,i)   = vxyzu(4,i)/rhoi
     fxyzu(4,i)  = 0
  enddo

  maxt = 5e-7*seconds
  t = 0.
  rhoi    = rhoh(xyzh(4,1),pmassi)
  physrho = rhoi*unit_density
  i = 0
  do while(t < maxt/utime)
     dt = max(1e-18*seconds/utime,0.05*t)
     ! dt = maxt/utime
     call update_radenergy(1,xyzh,fxyzu,vxyzu,radiation,dt)
     ! call solve_internal_energy_implicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
     ! call solve_internal_energy_explicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
     t = t + dt
     if (mod(i,10)==0) then
        laste = (vxyzu(4,1)*unit_ergg)*physrho
        write(24,*) t*utime, laste,(radiation(iradxi,1)*unit_ergg)*physrho
      end if
     i = i + 1
  enddo
  call checkval(laste,21195027.055207778,1e-10,ierr,'energy exchange for gas cooling')
  ntests = ntests + 1
  if (ierr == 0) npass = npass + 1

  do i=1,npart
     rhoi         = rhoh(xyzh(4,i),pmassi)
     radiation(iradxi,i) = 1e12/(unit_ergg*unit_density)/rhoi
     radiation(ikappa,i) = 0.4/(udist**2/umass)
     vxyzu(4,i)   = 1e2/(unit_ergg*unit_density)
     vxyzu(4,i)   = vxyzu(4,i)/rhoi
     fxyzu(4,i)  = 0
  enddo

  dt = 1e-11*seconds/utime
  t = 0.
  physrho = rhoi*unit_density
  i = 0
  do while(t < maxt/utime)
     dt = max(1e-18*seconds/utime,0.05*t)
     ! dt = maxt/utime
     call update_radenergy(1,xyzh,fxyzu,vxyzu,radiation,dt)
     t = t + dt
     if (mod(i,10)==0) then
        laste = (vxyzu(4,1)*unit_ergg)*physrho
        write(25,*) t*utime, laste,(radiation(iradxi,1)*unit_ergg)*physrho
      end if
     i = i + 1
  enddo
  call checkval(laste,21142367.365743987,1e-10,ierr,'energy exchange for gas heating')
  ntests = ntests + 1
  if (ierr == 0) npass = npass + 1
end subroutine

subroutine test_uniform_derivs(ntests,npass)
 use dim,             only:maxp
 use io,              only:id,master
 use part,            only:npart,xyzh,fxyzu,vxyzu,massoftype,igas,divcurlB,&
                           iphase,maxphase,isetphase,rhoh,&
                           bevol,divcurlv,fext,npartoftype,&
                           radiation,ifluxx,&
                           Bextx,Bexty,Bextz,maxvxyzu,dBevol,ddustevol,ddustprop,&
                           dustfrac,temperature,Bpred,vpred,dustproppred
 use kernel,          only:hfact_default
 use unifdis,         only:set_unifdis
 use units,           only:set_units,udist,utime,umass,unit_energ,unit_velocity,unit_ergg
 use linklist,        only:set_linklist
 use forces,          only:force
 use physcon,         only:c,Rg,pi,steboltz
 use eos,             only:gamma,gmw
 use readwrite_dumps, only:write_fulldump
 use boundary,        only:set_boundary
 use testutils,       only:checkvalbuf,checkvalbuf_end
 use deriv,           only:derivs
 use step_lf_global,  only:init_step,step

 integer,intent(inout) ::&
    ntests,npass

 real :: psep,hfact,a,c_code,cv1,rhoi,steboltz_code
 real :: dtmax,dtext,pmassi, dt,t,kappa_code
 real :: xmin,xmax,ymin,ymax,zmin,zmax,Tref,xi0,D0,rho0,l0
 real :: dtsph,dtnew,timei
 real :: exact_grE,exact_DgrF,exact_xi
 real :: errmax_e,errmax_f,tol_e,tol_f,errmax_xi,tol_xi

 integer :: i,j
 integer :: nactive,nerr_e,ncheck_e,nerr_f,ncheck_f,nerr_xi,ncheck_xi

 character(len=1024) :: filename

 psep = 1./32.
 hfact = hfact_default
 npart = 0
 xmin = -0.5
 xmax =  0.5
 ymin = -0.1
 ymax =  0.1
 zmin = -0.1
 zmax =  0.1
 call set_boundary(xmin,xmax,ymin,ymax,zmin,zmax)
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,hfact,npart,xyzh)
 massoftype(igas) = 1./npart*1e-25
 pmassi = massoftype(igas)
 if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(1) = npart

 vxyzu(4,:) = 0.
 fxyzu(:,:) = 0.
 fext(:,:)  = 0.
 Bevol(:,:) = 0.
 Bextx = 0.
 Bexty = 0.
 Bextz = 0.
 dBevol(:,:) = 0.
 divcurlv(:,:) = 0.
 nactive = npart

 call set_linklist(npart,nactive,xyzh,vxyzu)

 c_code = c/unit_velocity
 steboltz_code = steboltz/(unit_energ/(udist**2*utime))
 cv1 = (gamma-1.)*gmw/Rg*unit_velocity**2
 a   = 4.*steboltz_code/c_code
 pmassi = massoftype(igas)
 radiation(ithick,:) = 1.
 kappa_code = 1.0/(udist**2/umass)
 Tref = 100

 rho0 = rhoh(xyzh(4,1),pmassi)
 xi0 = a*Tref**4.0/rho0
 do i=1,npart
    vxyzu(4,i) = (Tref/cv1)/(unit_ergg)
    radiation(ikappa,i) = kappa_code
    radiation(iradxi,i) = xi0*(1 + 1e-1*sin(xyzh(1,i)*2*pi/(xmax-xmin)))
    ! etot = vxyzu(4,i) + radiation(iradxi,i)
    ! Tgas = vxyzu(4,i)*unit_ergg*cv1
    ! print*, vxyzu(4,i),radiation(iradxi,i),etot
    ! Trad = (rhoi*(etot-vxyzu(4,i))/a)**(1./4.)
    ! print*, Tref, Trad, Tgas
 enddo

do i = 1,2
 call derivs(1,npart,nactive,xyzh,vpred,fxyzu,fext,divcurlv,&
             divcurlB,Bpred,dBevol,dustproppred,ddustprop,dustfrac,&
             ddustevol,temperature,timei,dtsph,dtnew)
enddo
 nerr_e = 0
 ncheck_e = 0
 errmax_e = 0.
 tol_e = 1e-10

 nerr_f = 0
 ncheck_f = 0
 errmax_f = 0.
 tol_f = 2e-2

 l0 = 2*pi/(xmax-xmin)
 do i=1,npart
    rhoi = rhoh(xyzh(4,i),pmassi)
    D0  = c_code*(1./3)/kappa_code/rhoi
    exact_grE  =  xi0*rho0*0.1*l0   *cos(xyzh(1,i)*l0)
    exact_DgrF = -xi0*D0  *0.1*l0*l0*sin(xyzh(1,i)*l0)

    call checkvalbuf(radiation(ifluxx,i),exact_grE,tol_e, '  grad{E}',nerr_e,ncheck_e,errmax_e)
    call checkvalbuf(radiation(idflux,i),exact_DgrF,tol_f,'D*grad{F}',nerr_f,ncheck_f,errmax_f)
 enddo
 call checkvalbuf_end('  grad{E}',ncheck_e,nerr_e,errmax_e,tol_e)
 call checkvalbuf_end('D*grad{F}',ncheck_f,nerr_f,errmax_f,tol_f)

 ntests = ntests + 1
 if (nerr_e == 0) npass = npass + 1
 ntests = ntests + 1
 if (nerr_f == 0) npass = npass + 1

 dt = 1e-23
 t  = 0
 dtmax = dt
 dtext = dt
 call init_step(npart,t,dtmax)
 do i = 1,50
    t = t + dt
    dtext = dt
    call step(npart,nactive,t,dt,dtext,dtnew)

    if (mod(i,10) == 0) then
       nerr_xi = 0
       ncheck_xi = 0
       errmax_xi = 0.
       tol_xi = 2e-2
       do j = 1,npart
          rhoi = rhoh(xyzh(4,i),pmassi)
          D0  = c_code*(1./3)/kappa_code/rhoi
          exact_xi = xi0*(1.+0.1*sin(xyzh(1,i)*l0)*exp(-l0*l0*t*D0))
          write (filename,'(A,I2.2,A)') ' xi(t_', i, ')'
          call checkvalbuf(&
             radiation(iradxi,i),exact_xi,tol_xi,&
             trim(filename),nerr_xi,ncheck_xi,errmax_xi)
       enddo
       call checkvalbuf_end(trim(filename),ncheck_xi,nerr_xi,errmax_xi,tol_xi)
       ntests = ntests + 1
       if (nerr_e == 0) npass = npass + 1
    endif
    ! write (filename,'(A5,I2.2)') 'rad_test_', i
    ! call write_fulldump(t,filename)
 enddo
end subroutine

end module testradiation
