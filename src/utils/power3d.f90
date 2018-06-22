!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: power
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: omputils
!+
!--------------------------------------------------------------------------

!***********************************************************************
module power
 integer, parameter :: mfile=60,mlabel=60
 integer, parameter :: power_unit=14

contains

subroutine FAKE_NAMELIST(unit,name,var,string)
 implicit none
 integer unit
 character(len=*) :: name, var, string
 write (unit,'(1x,a)') "&"//trim(name)
 write (unit,'(1x,a)') trim(var)//" = '"//trim(string)//"'"
 write (unit,'(1x,a)') "/"
end subroutine

!***********************************************************************
subroutine open_power (file,origin,mk,xk,nk,n)
 implicit none
 character(len=mfile) file,origin
 character(len=mlabel) label
 integer mk,n
 real :: xk(mk)
 integer :: nk(mk)
 namelist /size/mk,n
 namelist /describe/label
!
 open(power_unit,file=file,status='unknown',form='formatted')
 write (power_unit,size)
 call fake_namelist(power_unit,"FILES","ORIGIN",origin)
 label='xk: average bin wave numbers'
 call fake_namelist(power_unit,"DESCRIBE","LABEL",label)
 write (power_unit,2) xk
 label='nk: number of wavenumbers in bin'
 call fake_namelist(power_unit,"DESCRIBE","LABEL",label)
 write (power_unit,3) nk
1 format(a)
2 format(8g15.7)
3 format(10i12)
end subroutine

!***********************************************************************
subroutine write_power (mk,rho_power,pk,ptot,case)
 implicit none
 integer mk
 real rho_power,ptot
 character(len=mlabel) case
 real :: pk(mk)
 namelist /component/rho_power,ptot
!
 call fake_namelist(power_unit,"DESCRIBE","LABEL",case)
 write (power_unit,component)
 write (power_unit,1) pk
1 format(8g15.7)
end subroutine

!***********************************************************************
subroutine power3d (ft,mx,my,mz,pk,xk,nk,mk,ptot,compensate,do_average,fty,ftz)
!
!  3D power spectrum, based on FFT output from FFTPACK (cos/sin real coeffs).
!  The calculations of kx,ky,kz would need to be modified if the transform
!  values are ordered differently.
!
!  The power spectrum is defined as the average power in each interval
!  times the k-space area of the interval, rather than the sum of the power
!  in each interval.  Empirically (and this may be deep), this gives
!  smoother power spectra.  So, nature cares more about the actual amplitudes
!  as a faction of wave number than about the total power in each interval.
!
 use omputils, only:limits_omp
 implicit none
 integer, intent(in)  :: mx, my, mz, mk
 integer, intent(out) :: nk(mk)
 real,    intent(in)  :: ft(mx,my,mz)
 real,    intent(out) :: pk(mk), xk(mk)
 real,    intent(out) :: ptot
 real,    intent(in),                            optional :: fty(mx,my,mz),ftz(mx,my,mz)
 logical do_average

 real, parameter :: pi=3.1415926
 real pk_thread(mk), xk_thread(mk), compensate, fk
 real(kind=8) ptotx, ptoty, ptot_thread
 real kx, ky, kz, k2, kk, cx, cy, cz, p2
 integer i, j, k, ik, nk_thread(mk)
 integer izs, ize
 character(len=72) :: id='$Id$'
!-----------------------------------------------------------------------
!$omp master
 if (id  /=  ' ') then; print *,id; id = ' '; endif
!$omp end master

!$omp single
 nk = 0
 xk = 0.
 pk = 0.
 ptot = 0.
!$omp end single
!$omp barrier

 nk_thread = 0
 xk_thread = 0.
 pk_thread = 0.
 ptot_thread = 0.

 call limits_omp(1,mz,izs,ize)
 do k=izs,ize
    kz = k/2
    cz = 2.
    if (k==1 .or. k==mz) cz = 1.
    ptoty = 0.
    do j=1,my
       ky = j/2
       cy = 2.
       if (j==1 .or. j==my) cy = 1.
       ptotx = 0.
       do i=1,mx
          kx = i/2
          cx = 2.
          if (i==1 .or. i==mx) cx = 1.
          k2 = kx**2 + ky**2 + kz**2
          kk = sqrt(k2)
          if (compensate == 0.) then
             fk = 1.
          else if (compensate == 2.) then
             fk = k2
          else
             fk = kk**compensate
          endif
          ik = 1.5 + kk                                                           ! wave number bin to add to
          p2 = (cx*cy*cz)*ft(i,j,k)**2                                            ! power contribution from this wavenumber
          if (present(fty)) then
             p2 = p2 + (cx*cy*cz)*fty(i,j,k)**2
          endif
          if (present(ftz)) then
             p2 = p2 + (cx*cy*cz)*ftz(i,j,k)**2
          endif
!       print '(1x,4i4,6g14.7)',i,j,k,ik,kx,ky,kz,k2,kk,p2                      ! debugging
          if (ik  <=  mk) then
             xk_thread(ik) = xk_thread(ik) + kk                                    ! wave number magnitude
             nk_thread(ik) = nk_thread(ik) + 1                                     ! number of wave numbers in this bin
             if (do_average) then
                pk_thread(ik) = pk_thread(ik) + fk*p2*k2*4.*pi                      ! power times 4*pi*k^2 (area factor)
             else
                pk_thread(ik) = pk_thread(ik) + fk*p2                               ! power
             endif
          endif
          ptotx = ptotx + p2                                                      ! Parseval
       enddo
       ptoty = ptoty + ptotx                                                     ! Parseval
    enddo
    ptot_thread = ptot_thread + ptoty                                           ! Parseval
 enddo

!$omp critical
 ptot = ptot + ptot_thread                                                     ! total power (Parseval)
 xk = xk + xk_thread                                                           ! wave number sum
 pk = pk + pk_thread                                                           ! power sum
 nk = nk + nk_thread                                                           ! number of wave numbers
!$omp end critical
!$omp barrier

!$omp single
 xk = xk/(nk+1e-10)                                                            ! average wave number
 xk(1) = 0.                                                                    ! first one always zero
 if (do_average) pk = pk/(nk+1e-10)                                            ! average power
!$omp end single
!$omp barrier

end subroutine
end module
