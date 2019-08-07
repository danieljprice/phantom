!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for measuring alpha
!  from an original code by Giuseppe Lodato
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'angmom'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 integer, parameter :: nr = 350
 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: xyzh(4,npart),vxyz(3,npart)
 real,             intent(in) :: pmass,time
 integer,          intent(in) :: npart,iunit

 real vr(nr)
 real rad(nr), sigma(nr), h_smooth(nr)
 real hr, r, hoverr, h_pressure
 real pi ,rmin,rmax,rad_i,cos_th,sin_th,vr_i
 real Star_M,G,disc_m,a,alpha_AV,cs0,R_in,R_out,b,C
 integer i, rad_int,no_u_pts(nr),numfile
 real, parameter :: rtracked = 5.
 real, parameter :: rtracked_tol = 0.01
 real retraced, rtracked_min, rtracked_max
 logical :: found_rtracked = .false.
 integer :: itracked

 logical :: vr_flag = .true.
 integer, parameter :: itsmax = 20
 real, parameter :: tol = 1.e-8
 real :: k1(nr), k2(nr), k3, ci(nr), f_fit(nr), df_fit(nr), ddf_fit(nr)
 real :: s1_fit, s2_fit, s3_fit
 real :: nu_i,nu_prev, tau, func, fderiv
 integer :: its

 character*8 output

 write(output,"(a5,i3.3)") 'alpha',numfile

 rmin=.5
 rmax=10.
 hr=(rmax-rmin)/(nr-1)
 r=rmin
 do i=1,nr
    rad(i)=r
    r=r+hr
 enddo

 pi = ACos(-1.0)

 do i=1,nr
    no_u_pts(i)=0
    sigma(i)=0.
    vr(i)=0.
    h_smooth(i) = 0.
 enddo


 hoverr = 0.02
 print*,' Assuming H/R = 0.02 at R=1: check that this is consistent with the setup used!!'
 G = 1.
 Star_M = 1.
 print*,' Assuming Star_M = 1: check that this is consistent with the setup used!!'
 disc_m = 0.01
 print*,' Assuming disc_m = 0.01: check that this is consistent with the setup used!!'
 a       = 2.
 print*,' Assuming a = 1: check that this is consistent with the setup used!!'
 R_in = 0.5
 print*,' Assuming R_in = 0.5: check that this is consistent with the setup used!!'
 R_out = 10.
 print*,' Assuming R_out = 10.: check that this is consistent with the setup used!!'

 rad_i = 0.
 cos_th = 0.
 sin_th = 0.
 vr_i = 0.

 rtrackedmin = rtracked*(1. - rtracked_tol)
 rtrackedmax = rtracked*(1. + rtracked_tol)

 do i = 1, npart

    rad_i=sqrt((xyzh(1,i))**2+(xyzh(2,i))**2)
    if (.not.found_rtracked .and. rad_i >= rtracked_min .and. rad_i <= rtracked_max) then
       itracked = i
       found_rtracked = .true.
    endif

    rad_int = INT((rad_i-rad(1))/hr + 1)

    If (rad_int  >  nr) rad_int = nr
    If (rad_int  <  1) rad_int = 1

    sigma(rad_int) = sigma(rad_int) + pmass/ &
            (pi*((rad(rad_int)+hr/2.)**2-(rad(rad_int)- hr/2.)**2))

    cos_th = xyzh(1,i)/rad_i
    sin_th = xyzh(2,i)/rad_i
    vr_i = cos_th*vxyz(1,i) + sin_th*vxyz(2,i)

    vr(rad_int) = vr(rad_int) + vr_i

    h_smooth(rad_int) = h_smooth(rad_int) + xyzh(4,i)

    no_u_pts(rad_int) = no_u_pts(rad_int) + 1

 enddo

 do i = 1, nr
    vr(i)=vr(i)/no_u_pts(i)

    h_pressure = hoverr * rad(i)**(.75)

    h_smooth(i) = h_smooth(i)/no_u_pts(i)/h_pressure

    if (no_u_pts(i)==0) then
       vr(i)=0.
       sigma(i)=0.
    endif
 enddo


! --- Determine the alphaSS of the flow
 if (vr_flag) then

! --- Fit the r*vr expression of LP74
    cs0 = (hoverr)*Sqrt(G*Star_M)
    alpha_AV = 0.01
    print *,'Guessed value of alphaAV: 0.01'
    nu_i = alpha_AV*2.1*(50000./real(npart))**(1./3.)
    b = sqrt(a)*(G*Star_M)*(R_out-R_in)
    C = 3.d0*nu_i*cs0*hoverr*a*disc_m*(G*Star_M)**2/&
       (1.d0 - exp(-(b)**2) - R_in*sqrt(pi)*sqrt(a)*(G*Star_M)*erf(-b))

    print*,'fit the r*vr profile at each timestep'

    do i=1,nr
       k1(i) = 1.5*cs0*hoverr
       k2(i) = 4.*a*(rad(i)- R_in)**2.
       ci(i) = vr(i)*(rad(i)- R_in)
       print*,'love',(rad(i)- R_in), vr(i),ci(i)
       f_fit(i) = 0.
       df_fit(nr) = 0.
       ddf_fit(nr) = 0.
    enddo
    k3 = 12.*a*(G*Star_M)**2.*cs0*hoverr*time

! --- Newton-Raphson for the err' (=func) function

    nu_prev = 10.*nu_i
    its = 0

    do while ((abs(nu_i-nu_prev) > tol).and.(its < itsmax))
       tau = k3 * nu_i + 1.
       nu_prev = nu_i

       s1_fit = 0.
       s2_fit = 0.
       s3_fit = 0.

       do i=1,nr
          if (no_u_pts(i) /= 0) then
             f_fit(i) = -k1(i)*nu_i*(1. - k2(i)/tau )
             df_fit(i) = -k1(i)*(1. - k2(i)/tau) - nu_i*k1(i)*k2(i)*k3/tau**2.
             ddf_fit(i) = -2.*k1(i)*k2(i)*k3/tau**2. + nu_i*2*k1(i)*k2(i)*k3**2./tau**3.

             s1_fit =  s1_fit + (f_fit(i) - ci(i))*df_fit(i)
             s2_fit =  s2_fit + (f_fit(i) - ci(i))**2.
             s3_fit =  s3_fit + (df_fit(i))**2. + (f_fit(i) - ci(i))*ddf_fit(i)
          endif
       enddo

       print*,'bite',nu_i,det,k1(i),k2(i),s1_fit,s2_fit,s3_fit
       func =  s1_fit/sqrt(nr*s2_fit)
       fderiv = (-(s1_fit)**2. + s2_fit*s3_fit)/(sqrt(real(nr))*(s2_fit)**1.5)
       print*,'bite2',func,fderiv

       nu_i = nu_i - func/fderiv     ! newton-raphson iteration
       its = its + 1
    enddo


! --- Fit the Sigma expression of LP74
    print*,'fit the Sigma profile at each timestep'

    do i=1,nr
       k1(i) = a*(rad(i)-R_in)**2.
       ci(i) = sigma(i)
       print*,'love',ci(i)
       f_fit(i) = 0.
       df_fit(nr) = 0.
       ddf_fit(nr) = 0.
    enddo
    k3 = 12.*a*(G*Star_M)**2.*cs0*hoverr*time

! --- Newton-Raphson for the err' (=func) function

    nu_prev = 10.*nu_i
    its = 0

    do while ((abs(nu_i-nu_prev) > tol).and.(its < itsmax))
       tau = k3 * nu_i + 1.
       nu_prev = nu_i

       s1_fit = 0.
       s2_fit = 0.
       s3_fit = 0.

       do i=1,nr
          if (no_u_pts(i) /= 0) then
             f_fit(i) = C*(tau)**(-1.25)/(3*pi*nu_i)*exp(-k1(i)/tau)
             df_fit(i) = -C*exp(-k1(i)/tau)*(9*k3**2.*nu_i**2.+13.*k3*nu_i+4.-4.*k1(i)*k3*nu_i)/&
                    (12.*pi*nu_i**2.*tau**(3.25))
             ddf_fit(i) = C*exp(-k1(i)/tau)*(nu_i**4.*(117.*k3**4.) + nu_i**3.*((338.-104.*k1(i))*k3**3.) +&
                                           nu_i**2.*(k3**2*(357. - 136.*k1(i))) + nu_i*(k3*(168. - 32.*k1(i))) + 32.)/&
                    (12.*pi*nu_i**3.*tau**(5.25))

             s1_fit =  s1_fit + (f_fit(i) - ci(i))*df_fit(i)
             s2_fit =  s2_fit + (f_fit(i) - ci(i))**2.
             s3_fit =  s3_fit + (df_fit(i))**2. + (f_fit(i) - ci(i))*ddf_fit(i)
          endif
       enddo

       print*,'bite',nu_i,det,k1(i),k2(i),s1_fit,s2_fit,s3_fit
       func =  s1_fit/sqrt(nr*s2_fit)
       fderiv = (-(s1_fit)**2. + s2_fit*s3_fit)/(sqrt(real(nr))*(s2_fit)**1.5)
       print*,'bite2',func,fderiv

       nu_i = nu_i - func/fderiv     ! newton-raphson iteration
       its = its + 1
    enddo

 endif

 if (its >= itsmax) then
    print*,'warning: bump max - too many iterations'
 endif
 print*,'iteration',its,'nu_i =',nu_i

 open(iunit,file=output)
 print *,'Writing to file... ',output

 do i=1,nr
    write(iunit,50) rad(i),sigma(i),vr(i),h_smooth(i),nu_i,h_smooth(itracked)
50  format(10(1pe16.8,1x))
 enddo

 write(iunit,*) time

 close(iunit)

 return
end subroutine do_analysis

end module analysis

