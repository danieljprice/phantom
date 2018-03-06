!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: powerspec
!
!  DESCRIPTION: An attempt to calculate powerspectra directly
!   from the SPH particle data. NOT RECOMMENDED.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: icosahedron
!+
!--------------------------------------------------------------------------
module powerspec
 implicit none
 real, parameter, private :: pi = 3.1415926536, twopi = 2.*pi

! resolution = 1        ! 0 pixels per face,  so  12 pixels in total.
! resolution = 2        ! 4 pixels per face,  so  92 pixels in total.
! resolution = 3        ! 12 pixels per face, so 252 pixels in total.
! resolution = 4        ! 24 pixels per face, so 492 pixels in total.

contains

!
! computes a grid of k values directly from the particle data
!
subroutine sphfft3D(dat,xyzh,rho,pmass,npart,datft,nkx,nky,nkz,dxmax)
 integer, intent(in)  :: npart,nkx,nky,nkz
 real,    intent(in)  :: dat(:),rho(:)
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(in)  :: pmass,dxmax
 real,    intent(out) :: datft(:,:,:)

 integer            :: k,j,i,nsubsample
 real :: xkx(3)
 real               :: fac,p2

 fac = twopi/dxmax
 nsubsample = 1

!--loop over k values
!$omp parallel do  default(none) schedule(static)  &
!$omp shared       (nkx,nky,nkz,npart)             &
!$omp shared       (xyzh,dat,rho,pmass,datft)      &
!$omp private      (k,j,i,xkx,p2)                  &
!$omp firstprivate (nsubsample,fac)
 do k=1,nkz
    xkx(3) = fac*real(k/2)
    print*,' k = ',k,xkx(3)
    do j=1,nky
       xkx(2) = fac*real(j/2)
       print*,' j = ',j,xkx(2)
       do i=1,nkx/2
          xkx(1) = fac*real(i)
          !print*,' i = ',1,xkx(1)
          call powerk_fourier3D(npart,nsubsample,xyzh,dat,rho,pmass,xkx,datft(2*i-1,j,k),datft(2*i,j,k),p2)
          !print*,' power = ',datft(i,j,k)**2
       enddo
    enddo
 enddo
!$omp end parallel do

end subroutine sphfft3D

!
! computes power on a grid of k values directly from the particle data
! (does not store the fourier coefficients - just returns the power summed and binned into k shells)
!
subroutine sphpow3D(dat,xyzh,rho,pmass,npart,pk,xk,nk,numk,ptot,nkx,nky,nkz,dxmax,do_average)
 integer, intent(in)  :: npart,nkx,nky,nkz,numk
 real,    intent(in)  :: dat(:),rho(:)
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(in)  :: pmass,dxmax
 real,    intent(out) :: pk(numk),xk(numk)
 integer, intent(out) :: nk(numk)
 real,    intent(out) :: ptot
 logical, intent(in)  :: do_average

 integer            :: k,j,i,nsubsample,ik
 real :: xkx(3)
 real               :: fac,p2,kx,ky,kz,cx,cy,cz,k2,kk,an,bn

 fac = twopi/dxmax
 nsubsample = 1
 xk(:) = 0.
 pk(:) = 0.
 ptot = 0.
 nk(:) = 0

!--loop over k values
!!$omp parallel do  default(none) schedule(static)  &
!!$omp shared       (nkx,nky,nkz,npart)             &
!!$omp shared       (xyzh,dat,rho,pmass,pk,xk)      &
!!$omp private      (k,j,i,xkx)                     &
!!$omp firstprivate (nsubsample,fac)
 do k=1,nkz
    kz = k/2
    cz = 2.
    if (k==1 .or. k==nkz) cz = 1.
    xkx(3) = fac*kz
    !print*,' k = ',k,xkx(3)
    do j=1,nky
       ky = j/2
       cy = 2.
       if (j==1 .or. j==nky) cy = 1.
       xkx(2) = fac*ky
       !print*,' j = ',j,xkx(2)
       do i=1,nkx
          kx = i/2
          cx = 2.
          if (i==1 .or. i==nkx) cx = 1.
          xkx(1) = fac*kx

          !print*,' i = ',1,xkx(1)
          !if (mod(i,2)==0) then
          call powerk_fourier3D(npart,nsubsample,xyzh,dat,rho,pmass,xkx,an,bn,p2)
          !else
          !   p2 = 0.
          !endif

          k2 = kx**2 + ky**2 + kz**2
          kk = sqrt(k2)
          ik = 1.5 + kk
          if (ik  <=  numk) then
             !print*,'bin ',ik,'i,j,k = ',i,j,k
             xk(ik) = xk(ik) + kk
             if (ik <= 2) print*,'ik=',ik,i,j,k,' power  = ',pk(ik),' + ',(cx*cy*cz)*p2,' an,bn=',an,bn
             if (do_average) then
                pk(ik) = pk(ik) + (cx*cy*cz)*p2*k2*4.*pi
             else
                pk(ik) = pk(ik) + (cx*cy*cz)*p2
             endif
             ptot = ptot + pk(ik)
             nk(ik) = nk(ik) + 1
          endif
          !print*,' power = ',datft(i,j,k)**2
       enddo
    enddo
 enddo
!!$omp end parallel do

 xk(:) = xk(:)/(nk(:)+epsilon(0.))
 xk(1) = 0.
 if (do_average) pk(:) = pk(:)/(nk(:)+epsilon(0.))

end subroutine sphpow3D


!
! computes a grid of k values directly from the particle data
!
! The output has the same format as fftpack, but normalised
! (i.e., this routine produces the same output as the fft3df
!  routine that we use to call fftpack for grids).
!
subroutine sphfft3D_fast(dat,xyzh,rho,pmass,npart,datft,nkx,nky,nkz,dxmax)
 integer, intent(in)  :: npart,nkx,nky,nkz
 real,    intent(in)  :: dat(:),rho(:)
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(in)  :: pmass,dxmax
 real,    intent(out) :: datft(:,:,:)

 integer            :: k,j,i,ip,lx,ly,lz,nsubsample
 real :: xkx(3),xyzpart(3)
 real               :: fac,hi2,rhoi,pmassonrhoi,dati,q2i,wk
 real               :: xkdotx,cos_term,sin_term,t1,t2

 fac = twopi/dxmax
 call cpu_time(t1)

 datft(:,:,:) = 0.
 nsubsample = 5

!!$omp parallel do  default(none) schedule(static)  &
!!$omp shared       (nkx,nky,nkz,npart)             &
!!$omp shared       (xyzh,dat,rho,pmass,datft)      &
!!$omp private      (k,j,i,xkx)                     &
!!$omp firstprivate (nsubsample,fac)
 do ip=1,npart,nsubsample
    if (mod(i-1,npart/100)==0) then
       call cpu_time(t2)
       print*,i,' timing = ',t2-t1,'s'
    endif
    xyzpart(1) = xyzh(1,ip)
    xyzpart(2) = xyzh(2,ip)
    xyzpart(3) = xyzh(3,ip)
    hi2 = xyzh(4,ip)**2
    rhoi = rho(ip)
    pmassonrhoi = pmass/rhoi
    dati = dat(ip)
    if (dati == rhoi) then
       dati = pmass
    else
       dati = pmassonrhoi*dati
    endif
    lz = nkz/2
    ly = nky/2
    lx = nkx/2

    !
    !--loop over k values
    !
    datft(:,:,1) = datft(:,:,1) + dati
    do k=2,lz
       xkx(3) = fac*(k-1)
       do j=2,ly
          xkx(2) = fac*(j-1)
          do i=2,lx
             xkx(1) = fac*(i-1)
             xkdotx = xkx(1)*xyzpart(1) + xkx(2)*xyzpart(2) + xkx(3)*xyzpart(3)

             cos_term = cos(xkdotx)
             sin_term = sin(xkdotx)
             !
             !--kernel smoothing
             !
             q2i = (xkx(1)*xkx(1) + xkx(2)*xkx(2) + xkx(3)*xkx(3))*hi2
             wk = real(nsubsample)*exp(-0.25*q2i)

             datft(2*i-2,2*j-2,2*k-2) = datft(2*i-2,2*j-2,2*k-2) + dati*cos_term*wk
             datft(2*i-1,2*j-2,2*k-2) = datft(2*i-1,2*j-2,2*k-2) - dati*sin_term*wk
             datft(2*i-2,2*j-1,2*k-2) = datft(2*i-2,2*j-1,2*k-2) + dati*cos_term*wk
             datft(2*i-1,2*j-1,2*k-2) = datft(2*i-1,2*j-1,2*k-2) - dati*sin_term*wk
             datft(2*i-2,2*j-2,2*k-1) = datft(2*i-2,2*j-2,2*k-1) + dati*cos_term*wk
             datft(2*i-1,2*j-2,2*k-1) = datft(2*i-1,2*j-2,2*k-1) - dati*sin_term*wk
             datft(2*i-2,2*j-1,2*k-1) = datft(2*i-2,2*j-1,2*k-1) + dati*cos_term*wk
             datft(2*i-1,2*j-1,2*k-1) = datft(2*i-1,2*j-1,2*k-1) - dati*sin_term*wk
          enddo
       enddo
    enddo
 enddo
 call cpu_time(t2)
 print*,'timing= ',t2-t1,'s'

end subroutine sphfft3D_fast

subroutine power_part(dxmax,dxmin,numk,xkgrid,pow,xyzh,dat,rho,massoftype,npart,iunit)
 use icosahedron, only:compute_matrices,compute_corners,pixel2vector
 real,    intent(in)  :: dxmin,dxmax
 integer, intent(in)  :: numk,npart,iunit
 real,    intent(out) :: xkgrid(numk),pow(numk)
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(in)  :: dat(:),rho(:),massoftype(:)
 integer, parameter :: iresolution = 3
 integer, parameter :: nsubk = 20*(2*iresolution*(iresolution-1)) + 12
 integer :: ik,isubk,iktot,i
 real :: xkmin,xkmax,xk,deltalogk,xkx(3),datmean
 real :: Rmat(0:19,3,3), vmat(0:11,3)
 real :: powerk,t2,t1,dati,xkdotx,dk,xkmaxbin,xkminbin!,dvolk
 real(kind=8) :: ptot
 real :: xkstore(3,numk*nsubk)
 real :: sums(4,numk*nsubk)
 real :: xyzpart(3)
 integer, parameter :: nkx = 128, nky=nkx, nkz=nkx
 real :: rhoi,pmassonrhoi,darea,sum1,sum2,sum3,sum4,ck2,cos_term,sin_term
 real :: q2i,wk,hi2,ak,bk
 logical, parameter :: loopoverk = .true.
 integer :: nsubsample

!---------------------------------------------
! divide up k-space between kmax and kmin

 xkmin = twopi/dxmax
 xkmax = twopi/dxmin
 print*,' max wavelength = ',dxmax,' min wavelength = ',dxmin,' frac = ',dxmin/dxmax
 print*,' min k = ',twopi/dxmax,' max k = ',twopi/dxmin
 deltalogk = (log10(xkmax) - log10(xkmin))/real(numk-1)
 print*,' numk = ',numk,' nk per shell = ',nsubk

!--initialise matrices for icosahedron package
 call compute_matrices(Rmat)
 call compute_corners(vmat)

 pow(:) = 0.
 datmean = 0. !sum(dat(1:npart))/real(npart)
 print*,' mean data value = ',datmean
 ptot = 0.
 nsubsample = 1

!--loop over k values
 iktot = 0
 do ik = 1,numk
    xk = xkmin*10**((ik-0.5)*deltalogk)
    xkmaxbin = xkmin*10**(ik*deltalogk)
    xkminbin = xkmin*10**((ik-1)*deltalogk)
    dk = xkmaxbin - xkminbin

    xkgrid(ik) = xk/(2.*pi)
    call cpu_time(t1)
    !--get the power at each k from a slow fourier transform
    do isubk = 1, nsubk
       iktot = iktot + 1
       !--choose kx, ky and kz evenly distributed on the sphere of radius k
       call pixel2vector(isubk-1,iresolution,Rmat,vmat,xkx)
       xkx(:) = xkx(:)*xk
       xkstore(:,iktot) = xkx(:)

       !print*,'xkx = ',xkx,' k = ',sqrt(dot_product(xkx,xkx)),xk
       if (loopoverk) then
          call powerk_fourier3D(npart,nsubsample,xyzh,dat,rho,massoftype(1),xkx,ak,bk,powerk)
          !call powerk_scargle3D(npart,nsubsample,xyzh(:,:),dat,datmean,xkx,powerk)
          !--total power in the shell is power at k
          pow(ik) = pow(ik) + powerk
       endif
    enddo
    if (loopoverk) then
       !dvolk = 4./3.*pi*(xkmaxbin**3 - xkminbin**3)
       !pow(ik) = pow(ik)*dvolk/real(xk*xk*nsubk)   ! using volume element

       darea = (4.*pi*xk**2)/real(nsubk)
       pow(ik) = pow(ik)*darea       ! assuming thin shell

       ptot = ptot + pow(ik)*dk

       call cpu_time(t2)
       print*,ik,': k = ',xkgrid(ik),' pow = ',pow(ik),' timing= ',t2-t1,'s'
       !print*,ik,xk,dk,dvolk,4.*pi*xk*xk*dk,pow(ik)
       if (iunit > 0) write(iunit,*) xkgrid(ik),pow(ik)
    endif
 enddo
 print*,'ptot = ',ptot
!
!--Alternative method: outer loop is over particles (slightly faster)
!
 if (.not.loopoverk) then
    call cpu_time(t1)
    !--loop over particles
    sums(:,:) = 0.
    nsubsample = 5
    do i=1,npart,nsubsample
       if (mod(i-1,npart/100)==0) then
          call cpu_time(t2)
          print*,i,' timing = ',t2-t1,'s'
       endif
       xyzpart(1) = xyzh(1,i)
       xyzpart(2) = xyzh(2,i)
       xyzpart(3) = xyzh(3,i)
       hi2 = xyzh(4,i)**2
       rhoi = rho(i)
       pmassonrhoi = massoftype(1)/rhoi
       dati = dat(i)
       if (dati == rhoi) then
          dati = massoftype(1)
       else
          dati = pmassonrhoi*dati
       endif
       iktot = 0
       do ik=1,numk
          !
          !--kernel smoothing
          !
          q2i = (2.*pi*xkgrid(ik))**2*hi2
          wk = exp(-0.25*q2i)
          !
          !--integrate over area
          !
          do isubk=1,nsubk
             iktot = iktot + 1
             xkdotx = xkstore(1,iktot)*xyzpart(1) + xkstore(2,iktot)*xyzpart(2) + xkstore(3,iktot)*xyzpart(3)

             cos_term = cos(xkdotx)
             sin_term = sin(xkdotx)
             sums(1,iktot) = sums(1,iktot) + dati*cos_term*wk
             sums(2,iktot) = sums(2,iktot) + dati*sin_term*wk
             sums(3,iktot) = sums(3,iktot) + pmassonrhoi*cos_term*cos_term
             sums(4,iktot) = sums(4,iktot) + pmassonrhoi*sin_term*sin_term
          enddo
       enddo
    enddo
    call cpu_time(t2)
    print*,'timing= ',t2-t1,'s'

    iktot = 0
    do ik=1,numk
       pow(ik) = 0.
       !--sum up the power in each shell by summing up the fourier coefficients
       !  over the shell area
       do isubk=1,nsubk
          iktot = iktot + 1
          sum1 = nsubsample*sums(1,iktot)
          sum2 = nsubsample*sums(2,iktot)
          sum3 = nsubsample*sums(3,iktot)
          sum4 = nsubsample*sums(4,iktot)
          ck2 = sum1**2/sum3**2 + sum2**2/sum4**2
          print*,ik,isubk,sum3,sum4
          pow(ik) = pow(ik) + 0.5*ck2
       enddo
       xk = xkgrid(ik)
       dk = xkmin*(10**(ik*deltalogk) - 10**((ik-1)*deltalogk))

       darea = (4.*pi*xk**2)/real(nsubk)
       pow(ik) = pow(ik)*darea
       !print*,ik,xk,dk,pow(ik)
       !write(8,*) ik,xk,pow(ik)
    enddo
 endif

 return
end subroutine power_part

!
!--calculate power at a given kx, ky and kz
!  using a slow fourier transform (looping over all particles)
!
subroutine powerk_fourier3D(npts,nsubsample,xyzh,dat,rho,pmassi,xk,an,bn,power)
 integer, intent(in)  :: npts,nsubsample
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(in)  :: dat(:),rho(:)
 real,    intent(in)  :: xk(3)
 real,    intent(in)  :: pmassi
 real,    intent(out) :: an,bn,power
 integer :: i
 real    :: dvoli,dati,xk2,xkdotx,q2i,wk
 real    :: sum1,sum2,sum3,sum4,cos_term,sin_term

 sum1 = 0.
 sum2 = 0.
 sum3 = 0.
 sum4 = 0.
!--get k^2, pre-multiply by 0.25 here (outside loop)
 xk2 = 0.25*(xk(1)**2 + xk(2)**2 + xk(3)**2)
 do i=1,npts,nsubsample
    xkdotx = xk(1)*xyzh(1,i) + xk(2)*xyzh(2,i) + xk(3)*xyzh(3,i)
    q2i = xk2*xyzh(4,i)**2
    wk = exp(-q2i)
    cos_term = cos(xkdotx)
    sin_term = sin(xkdotx)
    dvoli = pmassi/rho(i)
    dati = dat(i)
    sum1 = sum1 + dvoli*dati*cos_term*wk
    sum2 = sum2 + dvoli*dati*sin_term*wk
    sum3 = sum3 + dvoli*cos_term*cos_term
    sum4 = sum4 + dvoli*sin_term*sin_term
 enddo
 sum1 = 2.*sum1
 sum2 = 2.*sum2
 sum3 = 2. !*sum3
 sum4 = 2. !*sum4
 !print*,'sum3,sum4 = ',sum3,sum4
 an = real(nsubsample)*sum1/sum3
 if (abs(sum4) > epsilon(sum4)) then
    bn = real(nsubsample)*sum2/sum4
 else
    bn = 0.   ! gracefully handle the k=0 case
 endif
 power = 0.5*(an**2 + bn**2)

 return
end subroutine powerk_fourier3D

!
!--calculate power at a given kx, ky and kz
!  using a slow fourier transform (looping over all particles)
!  normalised as in Scargle (1982) ApJ 263, 835
!
subroutine powerk_scargle3D(npts,nsubsample,xyz,dat,datmean,xk,power)
 integer, intent(in)  :: npts,nsubsample
 real,    intent(in)  :: xyz(:,:)
 real,    intent(in)  :: dat(:)
 real,    intent(in)  :: datmean
 real,    intent(in)  :: xk(3)
 real,    intent(out) :: power
 integer :: i
 real    :: dati,sum1,sum2,xkdotx
 real    :: sum1denom,sum2denom,sin_term,cos_term

 power = 0.
 sum1 = 0.
 sum2 = 0.
 sum1denom = 0.
 sum2denom = 0.
 do i=1,npts,nsubsample
    xkdotx = xk(1)*xyz(1,i) + xk(2)*xyz(2,i) + xk(3)*xyz(3,i)
    !xkdotx = dot_product(xk(1:3),xyz(1:3,i))
    dati = dat(i) - datmean
    cos_term = cos(xkdotx)
    sin_term = sin(xkdotx)
    sum1 = sum1 + dati*cos_term
    sum2 = sum2 + dati*sin_term
    sum1denom = sum1denom + cos_term*cos_term
    sum2denom = sum2denom + sin_term*sin_term
 enddo
 power= real(nsubsample)*((sum1**2)/sum1denom + (sum2**2)/sum2denom)/real(npts) !*0.5

 return
end subroutine powerk_scargle3D

!
! opens the output file and writes a header
!
subroutine open_power_ascii(iunit,basename,variable,nk)
 integer,          intent(in) :: iunit
 character(len=*), intent(in) :: basename,variable
 integer,          intent(in) :: nk

 print "(a)",' writing to '//trim(basename)//'_'//trim(variable)//'.pow'
 open(unit=iunit,file=trim(basename)//'_'//trim(variable)//'.pow',form='formatted',status='replace')
 write(iunit,"('#',a,1x,i5,1x,1pe10.3)") trim(variable),nk

end subroutine open_power_ascii

!
! opens the output file and writes it
!
subroutine write_power_ascii(iunit,basename,variable,nk,xk,pk,ptot,comp)
 integer,          intent(in) :: iunit
 character(len=*), intent(in) :: basename,variable
 integer,          intent(in) :: nk
 real,             intent(in) :: xk(nk),pk(nk)
 real,             intent(in) :: ptot
 real,             intent(in),           optional :: comp(:)
 integer :: i,j

 print "(a)",' writing to '//trim(basename)//'_'//trim(variable)//'.pow'
 open(unit=iunit,file=trim(basename)//'_'//trim(variable)//'.pow',form='formatted',status='replace')
 write(iunit,"('# variable, nk, ptot  # power spectrum file written by phantom2power')")
 write(iunit,"('# ',a,1x,i5,1x,es10.3)") trim(variable),nk,ptot

 if (present(comp)) then
    write(iunit,"('#',2('[',i2.2,a12,']'),8('[',i2.2,a8,f4.2,']'))") &
          1,'k',2,'pk',(j+2,'pk*k^',comp(j),j=1,size(comp))
    do i=1,nk
       write(iunit,"(24(es10.4,1x))") xk(i),pk(i),(pk(i)*xk(i)**comp(j),j=1,size(comp))
    enddo
 else
    write(iunit,"('#',2('[',i2.2,a12,']'))") 1,'k',2,'pk'
    do i=1,nk
       write(iunit,"(2(es10.4,1x))") xk(i),pk(i)
    enddo
 endif
 close(iunit)

end subroutine write_power_ascii

end module powerspec
