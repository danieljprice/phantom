!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine computing particle tracks for a box
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'tracks'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use io,   only:fatal
 use part, only:iphase,iamtype,igas,idust

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: time,particlemass

 real               :: min_xgas_ini,min_xdust_ini,vinitrgas,vinitrdust,dragcste
 real               :: vgas_pred,vdust_pred,xgas_pred,xdust_pred,v_cdm
 real               :: expo, inv, inv_pow, alpha1,eps_deltav, sqrt_a3, func_a3, ch, sh, sqrt_a2, deltavx_PM
 integer            :: epsilon

 real               :: get_tstop
 integer            :: i,itrgas,itrdust
 logical            :: initialised = .false.
 real, parameter                  :: hrhogas    = 1.0008253226
 real, parameter                  :: hrhodust   = 0.5*1.0008253226!2.0016506443

 !--Select the drag structure:
 !--1: linear  -  2: powerlaw with exponent alpha  - 3: quadratic
 !--4: cubic expansion          -  5: Paardekooper and Mellema drag
 integer, parameter               :: idrag_structure = 1
 real, parameter                  :: dragcoeff  = 0.1
 real, parameter                  :: alpha      = 0.4
 real, parameter                  :: a2 = 0.2
 real, parameter                  :: a3 = 0.2

 integer, parameter               :: iepstein   = 0!
 real, parameter                  :: size       = 9.50523997326203253d-13 !1421.98 cm
 real, parameter                  :: massgrain  = 6.04924796450298316d-24
 real, parameter                  :: fracvolgas = 1.
 real, parameter                  :: spsoundgas = 1.

50 format(i8.8,i8.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8)
60 format(1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8,1pe16.8)

 !--Select a particle with the smallest x at the initial timestep
 if (.not.initialised) then
    min_xgas_ini  = huge(min_xgas_ini)
    min_xdust_ini = huge(min_xdust_ini)

    do i=1,npart
       if (iamtype(iphase(i))==igas) then
          if (xyzh(1,i) < min_xgas_ini) then
             min_xgas_ini = xyzh(1,i)
             itrgas = i
             vinitrgas = vxyzu(1,i)
          endif
       else
          if (xyzh(1,i) < min_xdust_ini) then
             min_xdust_ini = xyzh(1,i)
             itrdust = i
             vinitrdust = vxyzu(1,i)
          endif
       endif
    enddo
    open(81,file='output_tracked')
    write(81,50) itrgas,itrdust,vinitrgas,vinitrdust,min_xgas_ini,min_xdust_ini
    close(81)
    initialised = .true.
 else
    open(81,file='output_tracked')
    read(81,50) itrgas,itrdust,vinitrgas,vinitrdust,min_xgas_ini,min_xdust_ini
    close(81)
 endif

 if (itrgas==0 .or. itrdust==0) then
    call fatal('analysis_trackbox','gas or dust particle initially not tracked')
 endif

 !--calculate the theoretical evolution of the positions and the velocities
 !--initialisation of some general quantities
 dragcste   = dragcoeff*(1./hrhogas + 1./hrhodust)
 v_cdm      = (hrhogas*vinitrgas + hrhodust*vinitrdust)/(hrhogas + hrhodust)

 if (vinitrgas  >  vinitrdust) then
    epsilon = 1
 else
    epsilon = -1
 endif

 if (dragcoeff  <=   0. .or. a2  <=   0. .or. a3  <=   0.) then
    print*,'incorrect negative coeficient'
    stop
 endif

 !--calculation of the velocities
 select case(idrag_structure)
 case(1)
    if (iepstein == 1) then
       dragcste   = 4./3.*sqrt(8*3.141592653)*size*size/massgrain*spsoundgas*hrhogas/fracvolgas &
                      *hrhodust*(1./hrhogas + 1./hrhodust)
    endif
    expo       = exp(-dragcste*time)
    vgas_pred  = v_cdm + hrhodust/(hrhogas + hrhodust)*(vinitrgas - vinitrdust)*expo
    vdust_pred = v_cdm -  hrhogas/(hrhogas + hrhodust)*(vinitrgas - vinitrdust)*expo

    xgas_pred  = min_xgas_ini + ((hrhogas*vinitrgas + hrhodust*vinitrdust)*time - &
              hrhodust/dragcste*(vinitrgas - vinitrdust)*(expo-1.))/(hrhogas + hrhodust)
    xdust_pred = min_xdust_ini + ((hrhogas*vinitrgas + hrhodust*vinitrdust)*time + &
              hrhogas/dragcste*(vinitrgas - vinitrdust)*(expo-1.))/(hrhogas + hrhodust)

 case(2)
    alpha1     = 1./alpha
    eps_deltav = epsilon*(vinitrgas - vinitrdust)
    if (eps_deltav  <  0) then
       print*,'negative value of eps_deltav'
       stop
    endif
    inv_pow    = 1./(1. + alpha*dragcste*time*(eps_deltav)**alpha)**alpha1
    vgas_pred  = v_cdm + hrhodust/(hrhogas + hrhodust)*(vinitrgas - vinitrdust)*inv_pow
    vdust_pred = v_cdm -  hrhogas/(hrhogas + hrhodust)*(vinitrgas - vinitrdust)*inv_pow

 case(3)
    inv        = 1./(1. + epsilon*dragcste*(vinitrgas - vinitrdust)*time)
    vgas_pred  = v_cdm + hrhodust/(hrhogas + hrhodust)*(vinitrgas - vinitrdust)*inv
    vdust_pred = v_cdm -  hrhogas/(hrhogas + hrhodust)*(vinitrgas - vinitrdust)*inv

 case(4)
    expo    = exp(-dragcste*time)
    sqrt_a3 = sqrt(1. + a3*(vinitrgas - vinitrdust)**2)
    func_a3 = sqrt(1./(1. - a3*(vinitrgas - vinitrdust)**2*expo*expo/sqrt_a3**2  ))
    vgas_pred  = v_cdm + hrhodust/(hrhogas + hrhodust)*(vinitrgas - vinitrdust)*expo*func_a3/sqrt_a3
    vdust_pred = v_cdm -  hrhogas/(hrhogas + hrhodust)*(vinitrgas - vinitrdust)*expo*func_a3/sqrt_a3

 case(5)
    ch = cosh(dragcste*time)
    sh = sinh(dragcste*time)
    sqrt_a2 = sqrt(1. + a2*(vinitrgas - vinitrdust)**2)
    deltavx_PM = epsilon*sqrt( ( (sh + sqrt_a2*ch)/(ch + sqrt_a2*sh) )**2 - 1.)/sqrt(a2)
    vgas_pred  = v_cdm + hrhodust/(hrhogas + hrhodust)*deltavx_PM
    vdust_pred = v_cdm -  hrhogas/(hrhogas + hrhodust)*deltavx_PM

 end select

 !--track the particles and output the positions and the velocities
 open(224,file='outputgas',position='append')
 print *,'Writing to file... ','outputgas'
 open(225,file='outputdust',position='append')
 print *,'Writing to file... ','outputdust'
 open(226,file='outputpred',position='append')
 print *,'Writing to file... ','outputpred'
 open(227,file='outputsimu',position='append')
 print *,'Writing to file... ','outputsimu'

 if (idrag_structure ==1 ) then
    write(224,60)   time,vxyzu(1,itrgas),vgas_pred,vgas_pred -  vxyzu(1,itrgas), &
                       xyzh(1,itrgas),xgas_pred,xgas_pred -  xyzh(1,itrgas)
    write(225,60) time,vxyzu(1,itrdust),vdust_pred, vdust_pred -  vxyzu(1,itrdust), &
                    xyzh(1,itrdust),xdust_pred,xdust_pred - xyzh(1,itrdust)

    write(226,60) time,vdust_pred,vgas_pred,xdust_pred,xgas_pred
    write(227,60) time,vxyzu(1,itrdust),vxyzu(1,itrgas),xyzh(1,itrdust),xyzh(1,itrgas)
 else
    write(224,60)   time,vxyzu(1,itrgas),vgas_pred,vgas_pred -  vxyzu(1,itrgas)
    write(225,60) time,vxyzu(1,itrdust),vdust_pred, vdust_pred -  vxyzu(1,itrdust)

    write(226,60) time,vdust_pred,vgas_pred
    write(227,60) time,vxyzu(1,itrdust),vxyzu(1,itrgas)
 endif

 close(224)
 close(225)
 close(226)
 close(227)

 get_tstop = time/log((vinitrgas - vinitrdust)/(vxyzu(1,itrgas) - vxyzu(1,itrdust)))
 print*,'tsop from the simulation is...',get_tstop

end subroutine do_analysis



end module
