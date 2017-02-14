! $Id$

!*******************************************************************************
REAL FUNCTION interp(mx,my,mz,f,x,y,z)
!
! Tri-linear interpolation in f(:,:,:)
!
!-------------------------------------------------------------------------------
  implicit none
  integer mx,my,mz
  real :: f(mx,my,mz)
  real x,y,z
  real px,py,pz,qx,qy,qz
  integer ix,iy,iz,ixp1,iyp1,izp1

  px=x; ix=px; px=px-ix; qx=1.-px                                       ! interpolation weights
  py=y; iy=py; py=py-iy; qy=1.-py
  pz=z; iz=pz; pz=pz-iz; qz=1.-pz

  ix=1+mod(ix,mx); ixp1=1+mod(ix+1,mx)                                  ! periodic folding
  iy=1+mod(iy,my); iyp1=1+mod(iy+1,my)
  iz=1+mod(iz,mz); izp1=1+mod(iz+1,mz)

  interp = qz*(qy*(qx*f(ix,iy  ,iz  )+px*f(ixp1,iy  ,iz  ))  &          ! tri-linear interpolation
             + py*(qx*f(ix,iyp1,iz  )+px*f(ixp1,iyp1,iz  ))) &
         + pz*(qy*(qx*f(ix,iy  ,izp1)+px*f(ixp1,iy  ,izp1))  &
             + py*(qx*f(ix,iyp1,izp1)+px*f(ixp1,iyp1,izp1)))
  !print '(3f8.2,3i5,f10.3)',x,y,z,ix,iy,iz,interp
END

!*******************************************************************************
subroutine massflux2vel(mx,my,mz,lnrho,ui,i,off)                                ! convert mass flux to velocity
 use params
  use omputils, only:limits_omp
  implicit none
  integer mx,my,mz,i,iz,izs,ize
  real :: lnrho(mx,my,mz),ui(mx,my,mz)
  real, allocatable :: lnrho1(:,:,:)
  real off,off_in(3),off_out(3)
!
  allocate(lnrho1(mx,my,mz))                                                    ! shifted lnrho

  off_in=(/0.,0.,0./)
  if      (i==1) then
    off_out=(/off,0.,0./)
  else if (i==2) then
    off_out=(/0.,off,0./)
  else
    off_out=(/0.,0.,off/)
  endif

  !$omp parallel private(iz,izs,ize)
  call interpolate(lnrho,mx,my,mz,lnrho1,mx,my,mz,off_in,off_out)               ! interpolate to shifted position
  call limits_omp(1,mz,izs,ize)
  do iz=izs,ize
    ui(:,:,iz)=ui(:,:,iz)/exp(lnrho1(:,:,iz))
  enddo
  !$omp end parallel

  deallocate(lnrho1)
END

!*******************************************************************************
subroutine interpolate(vin,mxin,myin,mzin,vout,mxout,myout,mzout,off_in,off_out)
!
! Tri-linear interpolation, with offsets (e.g. staggering) on input and output.
! Sign check: If input data is Stagger Code data, the index position needs to
! be +0.5 higher to pick up the right velocity for centered data.  The sign on
! off_in should thus be negative.
!
!-------------------------------------------------------------------------------
 use params
  use omputils, only:limits_omp
  implicit none
  real vin(mxin,myin,mzin), vout(mxout,myout,mzout)
  integer mxin,myin,mzin,mxout,myout,mzout,ix,iy,iz,jx,jy,jz,jxp1,jyp1,jzp1,izs,ize
  real px,py,pz,qx,qy,qz,f,off_in(3),off_out(3)
  !logical omp_master
  !character(len=80),save:: myid="$Id$"
!-------------------------------------------------------------------------------
  !if (myid /= ' ' .and. omp_master()) print'(1x,a)',hl,myid,hl; myid=' '

  call limits_omp(1,mzout,izs,ize)
  do iz=izs,ize
    f=mzin/real(mzout)
    pz=(iz+off_out(3))*f-off_in(3)+mzin
    jz=pz
    pz=pz-jz
    qz=1.-pz
    jz=mod(jz-1,mzin)+1
    jzp1=mod(jz,mzin)+1
    do iy=1,myout
      f=myin/real(myout)
      py=(iy+off_out(2))*f-off_in(2)+myin
      jy=py
      py=py-jy
      qy=1.-py
      jy=mod(jy-1,myin)+1
      jyp1=mod(jy,myin)+1
      do ix=1,mxout
        f=mxin/real(mxout)
        px=(ix+off_out(1))*f-off_in(1)+mxin
        jx=px
        px=px-jx
        qx=1.-px
        jx=mod(jx-1,mxin)+1
        jxp1=mod(jx,mxin)+1
        vout(ix,iy,iz) = &
          qz*(qy*(qx*vin(jx,jy  ,jz  ) + px*vin(jxp1,jy  ,jz  )) + &
              py*(qx*vin(jx,jyp1,jz  ) + px*vin(jxp1,jyp1,jz  ))) + &
          pz*(qy*(qx*vin(jx,jy  ,jzp1) + px*vin(jxp1,jy  ,jzp1)) + &
              py*(qx*vin(jx,jyp1,jzp1) + px*vin(jxp1,jyp1,jzp1)))
      enddo
    enddo
    !print '(3i4,6f8.3)',iz,jz,jzp1,pz,qz,off_in(3),off_out(3),vin(1,1,jz),vout(1,1,iz)
  enddo
END

!***********************************************************************
subroutine interpolate3(mx,my,mz,ux,uy,uz,np,x,y,z,vx,vy,vz,offset)     ! tri-linear interpolation
!
! Since we assume the first cell center to be at 0.5/mx, the indices
! should change from 1 to 2 at 1/mx.  This routine assumes coordinates
! to lie between 0.0 and mx-epsilon.
!
!-----------------------------------------------------------------------
  use omputils, only:limits_omp
  implicit none
  integer mx,my,mz,np,ip,ix,iy,iz,ixp1,iyp1,izp1,nps,npe
  real :: ux(mx,my,mz),uy(mx,my,mz),uz(mx,my,mz)
  real :: x(np),y(np),z(np),vx(np),vy(np),vz(np)
  real px,py,pz,qx,qy,qz,offset

  call limits_omp(1,np,nps,npe)                                         ! OpenMP loop limits
  do ip=nps,npe                                                         ! interpolation loop
    px=x(ip)+mx-offset
              ix=px; px=px-ix; qx=1.-px                                 ! ix=0 -> 1 when x=1
    py=y(ip); iy=py; py=py-iy; qy=1.-py
    pz=z(ip); iz=pz; pz=pz-iz; qz=1.-pz
    ix=1+mod(ix-1+mx,mx); ixp1=1+mod(ix+mx,mx)                          ! ix=1 -> 2 when x=1
    iy=1+mod(iy-1+my,my); iyp1=1+mod(iy+my,my)
    iz=1+mod(iz-1+mz,mz); izp1=1+mod(iz+mz,mz)
    vx(ip) = qz*(qy*(qx*ux(ix,iy  ,iz  )+px*ux(ixp1,iy  ,iz  ))  &      ! tri-linear interpolation
               + py*(qx*ux(ix,iyp1,iz  )+px*ux(ixp1,iyp1,iz  ))) &
           + pz*(qy*(qx*ux(ix,iy  ,izp1)+px*ux(ixp1,iy  ,izp1))  &
               + py*(qx*ux(ix,iyp1,izp1)+px*ux(ixp1,iyp1,izp1)))
    px=x(ip); ix=px; px=px-ix; qx=1.-px
    py=y(ip)+my-offset
              iy=py; py=py-iy; qy=1.-py
    ix=1+mod(ix-1+mx,mx); ixp1=1+mod(ix+mx,mx)
    iy=1+mod(iy-1+my,my); iyp1=1+mod(iy+my,my)
    vy(ip) = qz*(qy*(qx*uy(ix,iy  ,iz  )+px*uy(ixp1,iy  ,iz  ))  &
               + py*(qx*uy(ix,iyp1,iz  )+px*uy(ixp1,iyp1,iz  ))) &
           + pz*(qy*(qx*uy(ix,iy  ,izp1)+px*uy(ixp1,iy  ,izp1))  &
               + py*(qx*uy(ix,iyp1,izp1)+px*uy(ixp1,iyp1,izp1)))
    py=y(ip); iy=py; py=py-iy; qy=1.-py
    pz=z(ip)+mz-offset
              iz=pz; pz=pz-iz; qz=1.-pz
    iy=1+mod(iy-1+my,my); iyp1=1+mod(iy+my,my)
    iz=1+mod(iz-1+mz,mz); izp1=1+mod(iz+mz,mz)
    vz(ip) = qz*(qy*(qx*uz(ix,iy  ,iz  )+px*uz(ixp1,iy  ,iz  ))  &
               + py*(qx*uz(ix,iyp1,iz  )+px*uz(ixp1,iyp1,iz  ))) &
           + pz*(qy*(qx*uz(ix,iy  ,izp1)+px*uz(ixp1,iy  ,izp1))  &
               + py*(qx*uz(ix,iyp1,izp1)+px*uz(ixp1,iyp1,izp1)))
  enddo
END

!-------------------------------------------------------------------------------
subroutine smooth3t (mx,my,mz,vin)
!
! Triple point smoothing, with weights w (1-2*w) w
!
  implicit none
  integer mx,my,mz,ix,iy,iz,ixp1,ixm1,iyp1,iym1,izp1,izm1
  real, parameter :: w=0.2
  real :: vin(mx,my,mz)
  real, allocatable :: vtmp(:,:,:)

  !print*,'smoothing ...',mx,my,mz
  allocate(vtmp(mx,my,mz))

  !$omp parallel do private(ix,iy,iz,izm1,izp1)
  do iz=1,mz
    izp1=mod(iz+mz,mz)+1
    izm1=mod(iz-2+mz,mz)+1
    do iy=1,my
      do ix=1,mx
        vtmp(ix,iy,iz)= w*vin(ix,iy,izm1) + (1.-2.*w)*vin(ix,iy,iz  ) + w*vin(ix,iy,izp1)
      enddo
    enddo
  enddo

  !print*,'loop2'
  !$omp parallel do private(ix,iy,iz,iym1,iyp1)
  do iz=1,mz
    do iy=1,my
      iyp1=mod(iy+my,my)+1
      iym1=mod(iy-2+my,my)+1
      do ix=1,mx
        vin(ix,iy,iz)= w*vtmp(ix,iym1,iz) + (1.-2.*w)*vtmp(ix,iy  ,iz) + w*vtmp(ix,iyp1,iz)
      enddo
    enddo
  enddo

  !print*,'loop3'
  !$omp parallel do private(ix,iy,iz,ixp1,ixm1)
  do iz=1,mz
    do iy=1,my
      !print*,iy,iz
      do ix=1,mx
        ixp1=mod(ix+mx,mx)+1
        ixm1=mod(ix-2+mx,mx)+1
        !print*,ixm1,ixp1,ix,iy,iz
        vtmp(ix,iy,iz)= w*vin(ixm1,iy,iz) + (1.-2.*w)*vin(ix  ,iy,iz) + w*vin(ixp1,iy,iz)
      enddo
      !print*,iy,iz
      do ix=1,mx
        vin(ix,iy,iz)=vtmp(ix,iy,iz)
      enddo
    enddo
  enddo

  !print*,'deall'
  deallocate(vtmp)

END
