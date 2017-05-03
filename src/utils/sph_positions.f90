!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: sph_positions
!
!  DESCRIPTION:
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, params
!+
!--------------------------------------------------------------------------

module sph_positions

contains

!***********************************************************************
subroutine positions(mx,my,mz,rho,npc,nptot,np,x,y,z,d,seed,id)
!
! This routine places particles in the ranges ! 0 < x < mx, 0 < y < my,
! 0 < z < mz, in such a way that they represent the mass density of a
! cube of size mx,my,mz, with periodic repetition in all three directions.
! The center of the 1st cell is at {0.5/mx,0.5/my,0.5/mz}.
!
!-----------------------------------------------------------------------
 use params,only:hl
 !use domain, only:ibelong
 use io, only:fatal
 integer, intent(in)    :: id
 integer, intent(inout) :: np
 integer, intent(out)   :: nptot
 integer mx,my,mz,ip,seed,ix,iy,iz,kpc,ipc,jx,jy,jz,np1,iptot
 real :: rho(mx,my,mz)
 real :: x(np),y(np),z(np),d(np)
 real fpc,npc,xx,yy,zz,interp,rhop,rho_max,facx,facy,facz
 character(len=80), save :: myid="$Id$"
!-------------------------------------------------------------------------------
 print'(1x,a)',hl,myid,hl

 ip = 0
 iptot = 0
 np1 = 0
 facx = 1./mx
 facy = 1./my
 facz = 1./mz

 print*,'positioning particles:'
 do iz=1,mz
    !if (mod(iz,100)==0) print'(i6,i12)',iz,ip
    do iy=1,my
       do ix=1,mx
          fpc = npc*rho(ix,iy,iz)                                             ! <number of particles>
          kpc = fpc                                                           ! integer number
          if (ran1s(seed) < (fpc-kpc)) kpc = kpc+1                            ! adjust expectation

          if (kpc > 9) then                                                   ! enough to make sense
             rho_max=0.                                                        ! max rho in cell
             do jz=0,1; do jy=0,1; do jx=0,1
                      xx=ix+jx-0.5; yy=iy+jy-0.5; zz=iz+jz-0.5                        ! corner coordinates
                      rho_max = max(rho_max,interp(mx,my,mz,rho,xx,yy,zz))            ! corner values
                   enddo; enddo; enddo
             ipc = iptot                                                          ! current particle index
             do while (iptot < ipc+kpc)                                           ! add kpc particles
                xx = (ix-1)+ran1s(seed)                                         ! proposed cell x
                yy = (iy-1)+ran1s(seed)                                         ! proposed cell y
                zz = (iz-1)+ran1s(seed)                                         ! proposed cell z
                rhop = interp(mx,my,mz,rho,xx,yy,zz)                            ! density there
                if (ran1s(seed) < rhop/rho_max) then                            ! check probability
                   iptot = iptot + 1
                   xx = xx*facx
                   yy = yy*facy
                   zz = zz*facz
                   !if (ibelong(xx,yy,zz,id)==id) then
                   ip = ip+1                                                     ! next particle
                   if (ip > size(x)) call fatal('position','ip > maxp')
                   x(ip)=xx; y(ip)=yy; z(ip)=zz                                  ! position OK
                   d(ip)=rhop
                   !endif
                   if(mod(iptot,10000000)==0 .and. ip > 0) print'(3i6,i6,2i12,4f9.2)',ix,iy,iz,kpc,iptot,ip,x(ip),y(ip),z(ip),d(ip)
                endif
             enddo
             np1 = np1+kpc
             !print*,kpc,np1,rho(ix,iy,iz)

          else                                                                ! at most a few
             do ipc=1,kpc                                                      ! add them uniformly
                xx = (ix-1)+ran1s(seed)                                      ! distribute evenly
                yy = (iy-1)+ran1s(seed)                                      ! distribute evenly
                zz = (iz-1)+ran1s(seed)                                      ! distribute evenly
                iptot = iptot + 1
                xx = xx*facx
                yy = yy*facy
                zz = zz*facz
                !if (ibelong(xx,yy,zz,id)==id) then
                ip = ip+1                                                       ! next particle
                if (ip > size(x)) call fatal('position','ip > maxp')
                x(ip) = xx; y(ip) = yy; z(ip) = zz
                d(ip) = rho(ix,iy,iz)                                           ! density there
                !endif
                if(mod(iptot,10000000)==0 .and. ip > 0) print'(3i6,i6,2i12,4f9.2)',ix,iy,iz,kpc,iptot,ip,x(ip),y(ip),z(ip),d(ip)
             enddo
          endif
       enddo
    enddo
 enddo

 np = ip                                                               ! actual number
 print*,'thread ',id,' got ',np,' of ',iptot
 nptot = iptot
 print*,'precision placement fraction =',np1/real(np)
end subroutine positions

!***********************************************************************
real function ran1s(idum)
!
! Simple pseudo-random number generator, after Press et al
!
!-----------------------------------------------------------------------
 integer idum,ia,im,iq,ir,ntab,ndiv,k,iy
 real am,eps,rnmx
 parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836, &
  ntab=32,ndiv=1+(im-1)/ntab,eps=1.2E-7,rnmx=1.-eps)

 if (idum <= 0) then
    idum=max(-idum,1)
 endif
 k=idum/iq
 idum=ia*(idum-k*iq)-ir*k
 if (idum < 0) idum=idum+im
 iy=idum
 ran1s=min(am*iy,rnmx)
end function ran1s

end module sph_positions
