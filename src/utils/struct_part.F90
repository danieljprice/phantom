!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: structurefn_part
!
!  DESCRIPTION:
!  module for obtaining structure functions
!  direct from SPH particles
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: random, timing
!+
!--------------------------------------------------------------------------
module structurefn_part
 implicit none

contains

subroutine get_structure_fn(sf,nbins,norder,distmin,distmax,xbins,ncount,npart,xyz,vel,&
                            rho,dxbox,dybox,dzbox,massweighted,ierr)
 !use fastmath, only:finvsqrt
 use timing,   only:get_timings,print_time
 use random,   only:ran2
 integer,         intent(in)  :: npart,nbins,norder
 real,            intent(in)  :: xyz(:,:)
 real,            intent(in)  :: vel(:,:)
 real,            intent(in)  :: rho(:)
 real(kind=8),    intent(out) :: sf(2,norder,nbins)
 real,            intent(out) :: xbins(nbins)
 integer(kind=8), intent(out) :: ncount(nbins)
 real,            intent(in)  :: distmax,distmin
 real,            intent(in)  :: dxbox,dybox,dzbox
 logical,         intent(in)  :: massweighted
 integer,         intent(out) :: ierr

 real(kind=8) :: sfprev(2,norder,nbins)
 integer, allocatable :: list(:)
 integer                   :: i,iran,ipart,ipt,iorder,ibin,iseed,npts,isf,nptstot,its
 real :: err(norder),sfmax(norder)
 real :: xpt(3),velpt(3)
 real                      :: dxbin,dvx,dvy,dvz,dx,dy,dz,rij1,rij
 real(kind=4)              :: t1,t2,tcpu1,tcpu2
 real                      :: rij2,distmin2,ddxbin,minusdistminddxbin
 real                      :: dvdotr,dvterm,dvtrans,rhomax,errtot,temp
 real(kind=8)              :: dvdotrterm,dvtransterm
 !$ integer                   :: omp_get_num_threads
 logical                   :: converged
!
!--set up the distance bins (linear)
!
 dxbin = (distmax-distmin)/float(nbins-1)
 do ibin=1,nbins
    xbins(ibin) = distmin + (ibin-0.5)*dxbin
 enddo
 distmin2 = distmin*distmin
 ddxbin = 1./dxbin
 minusdistminddxbin = -distmin*ddxbin
 ierr = 0
!
!--set structure functions to zero
!
 sf(:,:,:) = 0.
 sfprev(:,:,:) = 0.
 ncount(:) = 0
 iseed = -128
 npts = min(128,npart)
 nptstot = 0
 its = 0
!
!--start with a low number of points, and we keep adding more
!  points until the structure function calculation is converged
!
 converged = .false.
 !$omp parallel
 !$omp master
 !$ print*,' Using ',omp_get_num_threads(),' cpus'
 !$omp end master
 !$omp end parallel

 iterations: do while(nptstot  <=  npart .and. .not.converged)

    its = its + 1
    nptstot = nptstot + npts
    print "(a,i2,2(a,i10),a)",' Iteration ',its,': adding ',npts,' sample particles (',nptstot,' in total)'
    if (allocated(list)) deallocate(list)
    allocate(list(npts),stat=ierr)
    if (ierr /= 0) then
       print*,' error: cannot allocate memory for ',npts,' sample particles '
       sf = sfprev
       return
    endif
    print*,' iseed = ',iseed,' ncount(1:10) = ',ncount(1:10)

    !
    !--choose a random selection of npts particles
    !
    if (massweighted) then
       !
       !--select particles randomly according to particle id
       !  (this preferentially selects particles in dense regions)
       !
       do ipt=1,npts
          iran = int(ran2(iseed)*npart) + 1
          list(ipt) = iran
       enddo
    else
       !
       !--alternatively, select particles but weight selection by
       !  the volume element m/rho, i.e., inversely proportional to rho
       !
       rhomax = 0.
       !$omp parallel do schedule(static) private(i) reduction(max:rhomax)
       do i=1,npart
          rhomax = max(rho(i),rhomax)
       enddo
       if (rhomax <= 0.) then
          print*,' ERROR: max density on particles <= 0'
          print*,' cannot use volume element weighting for structure fns'
          return
       endif
       ipt = 0
       write(*,"(2x,a,i8,a)",ADVANCE='NO') 'choosing ',npts,' volume-weighted points...'
       do while(ipt < npts)
!--first random number chooses the particle
          iran = int(ran2(iseed)*npart) + 1
!--then select particle if rho/rhomax (0..1) is less than
!  a second random number
          if (rho(iran)/rhomax  <  ran2(iseed)) then
             ipt = ipt + 1
             list(ipt) = iran
          endif
       enddo
       print*,' done'
    endif

    call get_timings(t1,tcpu1)
    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(npts,xyz,vel,list,npart) &
    !$omp firstprivate(distmin2,dxbox,dybox,dzbox,ddxbin,norder,minusdistminddxbin) &
    !$omp private(ipt,xpt,velpt,dx,dy,dz,rij2,rij1,rij,dvdotr) &
    !$omp private(i,dvx,dvy,dvz) &
    !$omp private(dvterm,dvtrans,dvdotrterm,dvtransterm,ibin) &
    !$omp reduction(+:ncount) &
    !$omp reduction(+:sf)
    do ipt=1,npts
#ifndef _OPENMP
       if (mod(ipt,100)==0) then
          call cpu_time(tcpu2)
          print*,' ipt = ',ipt,tcpu2-tcpu1
       endif
#endif
       i = list(ipt)
       xpt(1) = xyz(1,i)
       xpt(2) = xyz(2,i)
       xpt(3) = xyz(3,i)
       velpt(1) = vel(1,i)
       velpt(2) = vel(2,i)
       velpt(3) = vel(3,i)

       do ipart=1,npart
          dx = xyz(1,ipart) - xpt(1)
          dy = xyz(2,ipart) - xpt(2)
          dz = xyz(3,ipart) - xpt(3)
          !--mod distances with periodic boundary
          if (abs(dx) > 0.5*dxbox) dx = dx - dxbox*sign(1.0,dx)
          if (abs(dy) > 0.5*dybox) dy = dy - dybox*sign(1.0,dy)
          if (abs(dz) > 0.5*dzbox) dz = dz - dzbox*sign(1.0,dz)

          rij2 = dx*dx + dy*dy + dz*dz
!
!--work out which distance bin this pair lies in
!  exclude pairs which lie closer than the minimum
!  separation bin
!
          if (rij2 > distmin2) then
             dvx = vel(1,ipart) - velpt(1)
             dvy = vel(2,ipart) - velpt(2)
             dvz = vel(3,ipart) - velpt(3)

             !       rij1 = finvsqrt(rij2)
             rij1 = 1./sqrt(rij2)

             dvdotr = abs((dvx*dx + dvy*dy + dvz*dz)*rij1)
             dvterm = (dvx*dvx + dvy*dvy + dvz*dvz) - dvdotr*dvdotr
             if (dvterm < 0.) dvterm = 0.
             dvtrans = sqrt(dvterm)

             rij = 1./rij1
             ibin = int(rij*ddxbin + minusdistminddxbin) + 1
             !if (ibin < 1 .or. ibin > nbins) stop 'ibin out of range'

             dvdotrterm = 1.0d0
             dvtransterm = 1.0d0
             do iorder=1,norder
                dvdotrterm = dvdotrterm*dvdotr     ! dvdotrterm = dvdotr**iorder
                dvtransterm = dvtransterm*dvtrans  ! dvtransterm = dvtrans**iorder

                sf(1,iorder,ibin) = sf(1,iorder,ibin) + dvdotrterm
                sf(2,iorder,ibin) = sf(2,iorder,ibin) + dvtransterm
             enddo
             ncount(ibin) = ncount(ibin) + 1_8
          endif
       enddo
    enddo
    !$omp end parallel do
    call get_timings(t2,tcpu2)
    call print_time(t2-t1,' wall time :')
    call print_time(tcpu2-tcpu1,' cpu time  :')

    err(:) = 0.
    sfmax(:) = 0.
    !$omp parallel do schedule(runtime) private(ibin) &
    !$omp reduction(+:err) &
    !$omp reduction(max:sfmax)
    do ibin=1,nbins
       if (ncount(ibin) > 0) then
          do iorder=1,norder
             do isf=1,2
                temp = sf(isf,iorder,ibin)/real(ncount(ibin))
                err(iorder) = err(iorder) + (temp - sfprev(isf,iorder,ibin))**2
                sfmax(iorder) = max(sfmax(iorder),temp)
                sfprev(isf,iorder,ibin) = temp
             enddo
          enddo
       else
          sfprev(:,:,ibin) = 0.
       endif
    enddo
    !$omp end parallel do

    errtot = 0.
    do iorder=1,norder
       if (sfmax(iorder) > 0.) then
          err(iorder) = err(iorder)/sfmax(iorder)**2/real(nbins*2)
       endif
       errtot = errtot + err(iorder)
       print*,' Error in structure function of order ',iorder,' = ',sqrt(err(iorder))
    enddo
    errtot = sqrt(errtot/real(norder))
    print*,' mean square error = ',errtot
    converged = maxval(sqrt(err(1:norder))) < 1.e-2 .and. errtot < 1.e-2
    npts = min(nptstot,npart-nptstot)

    !
    !--write the iterations to file (debugging only)
    !
    !do i=1,nbins
    !   write(10+its,*) xbins(i),(sfprev(1,iorder,i),iorder=1,norder)
    !enddo

 enddo iterations

 print*,' Converged!'

 !$omp parallel do schedule(static) private(ibin)
 do ibin=1,nbins
    sf(:,:,ibin) = sfprev(:,:,ibin)
 enddo

 if (allocated(list)) deallocate(list)

end subroutine get_structure_fn

end module structurefn_part
