#ifdef RADIATION
module testradiation
 implicit none
 public :: test_radiation

 private

contains

subroutine test_radiation(ntests,npass)
 use dim,          only:maxp
 use io,           only:id,master,iverbose
 use part,         only:npart,xyzh,fxyzu,vxyzu,massoftype,igas,divcurlB,&
                        iphase,maxphase,isetphase,radenergy,radkappa,rhoh,&
                        alphaind,bevol,gradh,divcurlv,fext,npartoftype,&
                        radenevol,radkappa,radenflux,radthick,dradenevol,&
                        maxvxyzu,dBevol,ddustevol,ddustprop,&
                        dustfrac,dustprop,temperature
 use kernel,       only:hfact_default
 use unifdis,      only:set_unifdis
 use units,        only:set_units
 use eos,          only:gmw,gamma,polyk
 use linklist,     only:set_linklist
 use ptmass,       only:ipart_rhomax
 use densityforce, only:densityiterate
 use forces,       only:force

 integer, intent(inout) :: ntests,npass

 real :: psep,hfact,rr
 real :: stressmax,pmassi,dt
 integer :: i,j,nactive

 if (id==master) write(*,"(/,a,/)") '--> TESTING RADIATION MODULE'

 iverbose = 1

 radenergy(:)    = 0
 radenevol(:)    = 0
 radkappa(:)     = 0
 radenflux(:,:)  = 0
 dradenevol(:)   = 0

 psep = 1./16.
 hfact = hfact_default
 call set_units(G=1.d0)
 npart = 0
 call set_unifdis('cubic',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,psep,hfact,npart,xyzh)
 massoftype(igas) = 1e-7/npart
 gamma = 5./3.
 gmw = 2.0
 polyk = 0.
 if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(1) = npart
 pmassi = massoftype(igas)

 call test_exchange_terms(npart,pmassi,radenergy,radenevol,radkappa,xyzh,vxyzu,fxyzu)

#ifndef PERIODIC
  if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF RADIATION MODULE (need -DPERIODIC)'
#else

  call test_uniform_fluxes_and_derivs()
  call test_random_fluxes_and_derivs()

#endif

 if (id==master) write(*,"(/,a)") '<-- RADIATION TEST COMPLETE'
end subroutine test_radiation

subroutine test_exchange_terms(npart,pmassi,radenergy,radenevol,radkappa,xyzh,vxyzu,fxyzu)
  use units,      only:unit_ergg,unit_density,umass,udist,utime
  use part,       only:rhoh
  use radiation,  only:update_radenergy
  use testutils,  only:checkval

  implicit none

  real,intent(inout) :: &
    pmassi
  real,intent(inout) :: &
    radenergy(:),radenevol(:),radkappa(:),xyzh(:,:),vxyzu(:,:),fxyzu(:,:)
  integer,intent(inout) ::&
    npart

  real :: dt,t,physrho,rhoi,maxt,laste
  integer :: i,ierr

  do i=1,npart
     radenergy(i) = 1e12/(unit_ergg*unit_density)
     rhoi         = rhoh(xyzh(4,i),pmassi)
     radenevol(i) = radenergy(i)/rhoi
     radkappa(i)  = 0.4/(udist**2/umass)
     vxyzu(4,i)   = 1e10/(unit_ergg*unit_density)
     vxyzu(4,i)   = vxyzu(4,i)/rhoi
     fxyzu(4,i)  = 0
  enddo

  maxt = 9e-7
  t = 0.
  rhoi    = rhoh(xyzh(4,1),pmassi)
  physrho = rhoi*unit_density
  i = 0
  do while(t < maxt/utime)
     ! dt = max(1e-16/utime,0.05*t)
     dt = maxt/utime/4
     call update_radenergy(1,xyzh,fxyzu,vxyzu,radenevol,radkappa,dt)
     t = t + dt
     ! if (mod(i,10)==0) then
        laste = (vxyzu(4,1)*unit_ergg)*physrho
        write(26,*) t*utime, laste,(radenevol(1)*unit_ergg)*physrho
      ! end if
     i = i + 1
  enddo
  call checkval(laste,21195027.055207778,3.e-4,ierr,'gas energy')

  do i=1,npart
     radenergy(i) = 1e12/(unit_ergg*unit_density)
     rhoi         = rhoh(xyzh(4,i),pmassi)
     radenevol(i) = radenergy(i)/rhoi
     radkappa(i)  = 0.4/(udist**2/umass)
     vxyzu(4,i)   = 1e2/(unit_ergg*unit_density)
     vxyzu(4,i)   = vxyzu(4,i)/rhoi
     fxyzu(4,i)  = 0
  enddo

  dt = 1e-11/utime
  t = 0.
  physrho = rhoi*unit_density
  i = 0
  do while(t < maxt/utime)
     ! dt = max(1e-16/utime,0.05*t)
     dt = maxt/utime/4
     call update_radenergy(1,xyzh,fxyzu,vxyzu,radenevol,radkappa,dt)
     t = t + dt
     ! if (mod(i,10)==0) then
        laste = (vxyzu(4,1)*unit_ergg)*physrho
        write(26,*) t*utime, laste,(radenevol(1)*unit_ergg)*physrho
      ! end if
     i = i + 1
  enddo
  call checkval(laste,21142736.646201313,3.e-4,ierr,'gas energy')
end subroutine

subroutine test_uniform_fluxes_and_derivs()
  use dim,          only:maxp
  use io,           only:id,master
  use part,         only:npart,xyzh,fxyzu,vxyzu,massoftype,igas,divcurlB,&
                         iphase,maxphase,isetphase,radenergy,radkappa,rhoh,&
                         alphaind,bevol,gradh,divcurlv,fext,npartoftype,&
                         radenevol,radkappa,radenflux,radthick,dradenevol,&
                         Bextx,Bexty,Bextz,maxvxyzu,dBevol,ddustevol,ddustprop,&
                         dustfrac,dustprop,temperature
  use kernel,       only:hfact_default
  use unifdis,      only:set_unifdis
  use units,        only:set_units
  use linklist,     only:set_linklist
  use ptmass,       only:ipart_rhomax
  use densityforce, only:densityiterate
  use forces,       only:force

  real :: psep,hfact
  real :: stressmax,dtmax,dtext,pmassi, dt,t
  integer :: i,j,nactive

  psep = 1./16.
  hfact = hfact_default
  npart = 0
  call set_unifdis('cubic',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,psep,hfact,npart,xyzh)
  massoftype(igas) = 1./npart
  pmassi = massoftype(igas)
  if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
  npartoftype(:) = 0
  npartoftype(1) = npart

  dt = 1e-15
  t  = 0
  dtmax = dt
  dtext = dt
  if (maxvxyzu >= 4) vxyzu(4,:) = 0.
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
  do i=1,20
     radenevol(:) = 1.0
     radenflux(:,:) = 0.
     radthick(:) = .true.
     radkappa(:) = 4e1
     dradenevol(:) = 0.

    call densityiterate(1,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,&
                        Bevol,stressmax,fxyzu,fext,alphaind,gradh,&
                        radenevol,radenflux,radenergy,radthick)
    call force(1,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,&
               dBevol,dustprop,ddustprop,dustfrac,ddustevol,&
               ipart_rhomax,dt,stressmax,temperature,&
               radenevol,radenflux,radkappa,dradenevol)

     do j=2,npart
       radenflux(:,1) = radenflux(:,1) + radenflux(:,j)
       dradenevol(1)  = dradenevol(1)  + dradenevol(j)
     enddo
     print'(a,g5.2,a,g9.2,2x,g9.2,2x,g9.2,a,g9.2)',&
           "Uniform | iter = ",i," | d{rad}/dx_i = ", radenflux(:,1)," | d^2{rad}/dt^2 = ", dradenevol(1)
  enddo
end subroutine

subroutine test_random_fluxes_and_derivs()
  use dim,          only:maxp
  use io,           only:id,master
  use part,         only:npart,xyzh,fxyzu,vxyzu,massoftype,igas,divcurlB,&
                         iphase,maxphase,isetphase,radenergy,radkappa,rhoh,&
                         alphaind,bevol,gradh,divcurlv,fext,npartoftype,&
                         radenevol,radkappa,radenflux,radthick,dradenevol,&
                         Bextx,Bexty,Bextz,maxvxyzu,dBevol,ddustevol,ddustprop,&
                         dustfrac,dustprop,temperature
  use kernel,       only:hfact_default
  use unifdis,      only:set_unifdis
  use units,        only:set_units
  use linklist,     only:set_linklist
  use ptmass,       only:ipart_rhomax
  use densityforce, only:densityiterate
  use forces,       only:force

  real :: psep,hfact,rr
  real :: stressmax,dtmax,dtext,pmassi,dt,t
  integer :: i,j,nactive

  psep = 1./16.
  hfact = hfact_default
  npart = 0
  call set_unifdis('random',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,psep,hfact,npart,xyzh)
  massoftype(igas) = 1./npart
  pmassi = massoftype(igas)
  if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
  npartoftype(:) = 0
  npartoftype(1) = npart

  dt = 1e-15
  t  = 0
  dtmax = dt
  dtext = dt
  if (maxvxyzu >= 4) vxyzu(4,:) = 0.
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
  do i=1,20
    do j=1,npart
       call random_number(rr)
       radenevol(j) = 1.0 + 0.1*(rr - 0.5)
       radenflux(:,j) = 0.
       radthick(j) = .true.
       radkappa(j) = 4e1
       dradenevol(j) = 0.
    enddo

    call densityiterate(1,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,&
                        Bevol,stressmax,fxyzu,fext,alphaind,gradh,&
                        radenevol,radenflux,radenergy,radthick)
    call force(1,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,&
               dBevol,dustprop,ddustprop,dustfrac,ddustevol,&
               ipart_rhomax,dt,stressmax,temperature,&
               radenevol,radenflux,radkappa,dradenevol)

     do j=2,npart
       radenflux(:,1) = radenflux(:,1) + radenflux(:,j)
       dradenevol(1)  = dradenevol(1)  + dradenevol(j)
     enddo
     print'(a,g5.2,a,g9.2,2x,g9.2,2x,g9.2,a,g9.2)',&
           "Random  | iter = ",i," | d{rad}/dx_i = ", radenflux(:,1)," | d^2{rad}/dt^2 = ", dradenevol(1)
  enddo
end subroutine

end module testradiation
#endif
