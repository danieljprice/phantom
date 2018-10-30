!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: mesa_microphysics
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
!  DEPENDENCIES: datafiles
!+
!--------------------------------------------------------------------------
module mesa_microphysics
 use datafiles, only:find_phantom_datafile

 implicit none
 public

!For mesa_microphysics
 integer            :: mesa_eos_nz=3, mesa_eos_nh=5   !Number of Z values in EOS tables, number of hydrogen values in EOS tables
 real               :: mesa_eos_z1=0.d0, mesa_eos_h1=0.d0, mesa_eos_dz=0.02d0, mesa_eos_dh=0.2d0 !lowest Z, lowest H, delta Z, delta H in EOS tables
 character(len=200) :: mesa_eos_dir="mesa_data",mesa_opacs_dir="mesa_dir",mesa_eos_prefix="eos_",mesa_opacs_suffix="" ! path to EOS tables, path to opacity tables, prefix for EOS table files, prefix for oapcity table files
 logical            :: mesa_eos_full_output=.false.,mesa_binary_data=.true.  ! EOS table format (if ASCII), EOS table format (ASCII=t,binary=f)

!Opacity
 integer                           :: mesa_opacs_nz,mesa_opacs_nx,mesa_opacs_nr,mesa_opacs_nt
 real                              :: mesa_opacs_dr,mesa_opacs_dt
 real, dimension(:), allocatable   :: mesa_opacs_zs,mesa_opacs_xs, mesa_opacs_rs,mesa_opacs_ts
 real, dimension(:,:), allocatable :: mesa_opacs_k,mesa_opacs_kd,mesa_opacs_kt

!EOS
 integer                                 :: mesa_eos_ne, mesa_eos_nv, mesa_eos_nvar2
 real                                    :: mesa_eos_v1, mesa_eos_e1, mesa_eos_de, mesa_eos_dv
 real, dimension(:), allocatable         :: mesa_eos_z, mesa_eos_h, mesa_eos_logEs, mesa_eos_logVs
 integer, dimension(:,:), allocatable    :: mesa_eos_data_exists
 real, dimension(:,:,:), allocatable     :: mesa_de_data
 real, dimension(:,:,:,:), allocatable   :: mesa_eos0
 real, dimension(:,:,:,:,:), allocatable :: mesa_de_data0

 public :: get_opacity_constants_mesa
 public :: read_opacity_mesa
 public :: get_kappa_mesa
 public :: get_eos_constants_mesa
 public :: read_eos_mesa
 !public :: fasteossc_mesa
 !public :: sloweossc_mesa
 !public :: getgamma_mesa
 !public :: getgamma1_mesa
 !public :: getpressure_mesa
 public :: getvalue_mesa

 private :: eos_cubic_spline_mesa
 private :: cubic_spline_mesa

contains

! Get the constants to be used in the MESA opacities
subroutine get_opacity_constants_mesa
 character(len=500) :: opacs_file,filename
 integer            :: fnum

 ! Find the opacity table
 filename   = 'opacs'//trim(mesa_opacs_suffix)//'.bindata'
 opacs_file = find_phantom_datafile(filename,'eos/mesa')

 ! Read the constants from the header of the opacity file
 open(newunit=fnum, file=trim(opacs_file), status='old', action='read', form='unformatted')
 read(fnum) mesa_opacs_nz,mesa_opacs_nx,mesa_opacs_nr,mesa_opacs_nt
 close(fnum)

 ! Allocate ALL the arrays!
 allocate(mesa_opacs_zs(mesa_opacs_nz),mesa_opacs_xs(mesa_opacs_nx),mesa_opacs_rs(mesa_opacs_nr),mesa_opacs_ts(mesa_opacs_nt))
 allocate(mesa_opacs_k(mesa_opacs_nr,mesa_opacs_nt),mesa_opacs_kd(mesa_opacs_nr,mesa_opacs_nt))
 allocate(mesa_opacs_kt(mesa_opacs_nr,mesa_opacs_nt))

 return

end subroutine get_opacity_constants_mesa

! Read opacity values for a single composition
subroutine read_opacity_mesa(x,z)
 real, intent(in) :: x,z
 real :: dz, dx
 real, dimension(:,:,:,:), allocatable :: kappas
 character(len=500) :: opacs_file,filename
 integer :: zz, xx, k, i
 character(len=7) :: empty
 integer :: nz2,nx2
 integer :: fnum

 fnum=121

 empty='       '


 !Find the opacity table
 filename = trim(mesa_opacs_dir)//'/'//'opacs'//trim(mesa_opacs_suffix)//'.bindata'
 opacs_file = find_phantom_datafile(filename,'eos/mesa')
 open(unit=fnum, file=trim(opacs_file), status='old', action='read', form='unformatted')
 read(fnum) mesa_opacs_nz,mesa_opacs_nx,mesa_opacs_nr,mesa_opacs_nt

 ! Read in the size of the table and the data
 ! zz and xx hold the Z and X values respectively
 ! i and k hold the values of r and temperature respectively

 allocate(kappas(mesa_opacs_nz,mesa_opacs_nx,mesa_opacs_nr,mesa_opacs_nt))

 read(fnum) (mesa_opacs_zs(zz),zz=1,mesa_opacs_nz)
 read(fnum) (mesa_opacs_xs(xx),xx=1,mesa_opacs_nx)
 read(fnum) (mesa_opacs_ts(k),k=1,mesa_opacs_nt)
 read(fnum) (mesa_opacs_rs(i),i=1,mesa_opacs_nr)
 do zz=1,mesa_opacs_nz
    do xx=1,mesa_opacs_nx
       do k=1,mesa_opacs_nt
          read(fnum) (kappas(zz,xx,i,k),i=1,mesa_opacs_nr)
       enddo
    enddo
 enddo
 close(fnum)

 ! Chooses the Z values from the available values which is greater than the input Z
 nz2=mesa_opacs_nz
 do zz=2,mesa_opacs_nz-1
    if(mesa_opacs_zs(zz) >= z) then
       nz2=zz
       exit
    endif
 enddo

 ! Finds fraction of input Z distance between two enclosing Z (sorry if confusing, just read the code)
 dz=(z-mesa_opacs_zs(nz2-1))/(mesa_opacs_zs(nz2)-mesa_opacs_zs(nz2-1))

 ! Same for X
 nx2=mesa_opacs_nx
 do xx=2,mesa_opacs_nx-1
    if(mesa_opacs_xs(xx) >= x) then
       nx2=xx
       exit
    endif
 enddo
 dx=(x-mesa_opacs_xs(nx2-1))/(mesa_opacs_xs(nx2)-mesa_opacs_xs(nx2-1))

 ! Interpolate kappas between X and Z values
 do k=1,mesa_opacs_nt
    do i=1,mesa_opacs_nr
       mesa_opacs_k(i,k) =      dz  *       dx  * kappas(nz2,nx2,i,k)   + &
                                dz  * (1.d0-dx) * kappas(nz2,nx2-1,i,k) + &
                          (1.d0-dz) *       dx  * kappas(nz2-1,nx2,i,k) + &
                          (1.d0-dz) * (1.d0-dx) * kappas(nz2-1,nx2-1,i,k)
    enddo
 enddo

 ! Calculate the partial derivative of kappa with respect to r
 do k=1,mesa_opacs_nt
    do i=2,mesa_opacs_nr-1
       mesa_opacs_kd(i,k) = (mesa_opacs_k(i+1,k) - mesa_opacs_k(i-1,k)) / &
                            (mesa_opacs_rs(i+1)  - mesa_opacs_rs(i-1))
    enddo

    mesa_opacs_kd(1,k) = mesa_opacs_kd(3,k) + (mesa_opacs_rs(1)   - mesa_opacs_rs(3))   * &
                                              (mesa_opacs_kd(2,k) - mesa_opacs_kd(3,k)) / &
                                              (mesa_opacs_rs(2)   - mesa_opacs_rs(3))

    mesa_opacs_kd(mesa_opacs_nr,k) =  mesa_opacs_kd(mesa_opacs_nr+1-3,k) + &
                                     (mesa_opacs_rs(mesa_opacs_nr+1-1)   - mesa_opacs_rs(mesa_opacs_nr+1-3))   * &
                                     (mesa_opacs_kd(mesa_opacs_nr+1-2,k) - mesa_opacs_kd(mesa_opacs_nr+1-3,k)) / &
                                     (mesa_opacs_rs(mesa_opacs_nr+1-2)   - mesa_opacs_rs(mesa_opacs_nr+1-3))

 enddo

 ! Calculate the partial derivative of kappa with respect to temperature
 do i=1,mesa_opacs_nr
    do k=2,mesa_opacs_nt-1
       mesa_opacs_kt(i,k)= (mesa_opacs_k(i,k+1) - mesa_opacs_k(i,k-1)) / &
                           (mesa_opacs_ts(k+1)  - mesa_opacs_ts(k-1)) - &
                            3.d0 * mesa_opacs_kd(i,k) ! I don't know why this part is here (or below)
    enddo

    mesa_opacs_kt(i,1) = mesa_opacs_kt(i,3) + (mesa_opacs_ts(1)   - mesa_opacs_ts(3))   * &
                                              (mesa_opacs_kt(i,2) - mesa_opacs_kt(i,3)) / &
                                              (mesa_opacs_ts(2)   - mesa_opacs_ts(3))   - &
                                               3.d0 * mesa_opacs_kd(i,1)


    mesa_opacs_kt(i,mesa_opacs_nt+1-1) =  mesa_opacs_kt(i,mesa_opacs_nt+1-3) + &
                                         (mesa_opacs_ts(mesa_opacs_nt+1-1)   - mesa_opacs_ts(mesa_opacs_nt+1-3))   * &
                                         (mesa_opacs_kt(i,mesa_opacs_nt+1-2) - mesa_opacs_kt(i,mesa_opacs_nt+1-3)) / &
                                         (mesa_opacs_ts(mesa_opacs_nt+1-2)   - mesa_opacs_ts(mesa_opacs_nt+1-3))   - &
                                          3.d0 * mesa_opacs_kd(i,mesa_opacs_nt+1-1)

 enddo

 ! Save some values for faster interpolation
 mesa_opacs_dr = mesa_opacs_rs(2) - mesa_opacs_rs(1)
 mesa_opacs_dt = mesa_opacs_ts(2) - mesa_opacs_ts(1)

 deallocate(kappas)

 return
end subroutine



! Return value of kappa for a value of density and temperature
subroutine get_kappa_mesa(rho,temp,kap,kapt,kapr)

 real, intent(in) :: rho,temp
 real, intent(out) :: kap,kapt,kapr
 real :: opac_k,opac_kd,opac_kt
 real :: dnr, dnt
 integer :: nr,nt
 real :: logrho,logt,logr,dr,dt

 logrho = log10(rho)
 logt   = log10(temp)
 logr   = logrho + 18.d0 - 3.d0 * logt

 ! Get the r and T indices for looking up the table
 dnr = 1.d0 + ((logr - mesa_opacs_rs(1)) / mesa_opacs_dr)
 dnt = 1.d0 + ((logt - mesa_opacs_ts(1)) / mesa_opacs_dt)

 nr = int(dnr)
 nt = int(dnt)

 ! Allow extrapolation
 if (nr > (mesa_opacs_nr-1)) nr = mesa_opacs_nr-1
 if (nr < 1) nr = 1

 if (nt > (mesa_opacs_nt-1)) nt = mesa_opacs_nt-1
 if (nt < 1) nt = 1

 dr = (logr - mesa_opacs_rs(nr)) / mesa_opacs_dr
 dt = (logt - mesa_opacs_ts(nt)) / mesa_opacs_dt

 ! Interpolate kappa (and derivatives?)
 opac_k  =       dr  * (dt * mesa_opacs_k(nr+1,nt+1)  + (1.d0-dt) * mesa_opacs_k(nr+1,nt)) + &
           (1.d0-dr) * (dt * mesa_opacs_k(nr,nt+1)    + (1.d0-dt) * mesa_opacs_k(nr,nt))

 opac_kt =       dr  * (dt * mesa_opacs_kt(nr+1,nt+1) + (1.d0-dt) * mesa_opacs_kt(nr+1,nt)) + &
           (1.d0-dr) * (dt * mesa_opacs_kt(nr,nt+1)   + (1.d0-dt) * mesa_opacs_kt(nr,nt))

 opac_kd =       dr  * (dt * mesa_opacs_kd(nr+1,nt+1) + (1.d0-dt) * mesa_opacs_kd(nr+1,nt)) + &
           (1.d0-dr) * (dt * mesa_opacs_kd(nr,nt+1)   + (1.d0-dt) * mesa_opacs_kd(nr,nt))


 kap  = 10.d0**opac_k
 kapt = opac_kt * kap / temp
 kapr = opac_kd * kap / rho

 return

end subroutine get_kappa_mesa



! Get the constants to be used in the MESA EoS
subroutine get_eos_constants_mesa(ierr)
 integer, intent(out) :: ierr
 character (len=20)   :: zz, hh
 character (len=300)  :: filename
 character (len=8)    :: fmt1
 integer              :: fnum, i

 ! Allocate the arrays to carry X and Z values.
 allocate(mesa_eos_z(mesa_eos_nz), mesa_eos_h(mesa_eos_nh))

 ! Assigning the X and Z values of the EoS tables to be read in (based on parameters at the beginning of the file.
 mesa_eos_z(1)=mesa_eos_z1
 mesa_eos_h(1)=mesa_eos_h1

 do i=2,mesa_eos_nz
    mesa_eos_z(i)=mesa_eos_z(i-1)+mesa_eos_dz
 enddo

 do i=2,mesa_eos_nh
    mesa_eos_h(i)=mesa_eos_h(i-1)+mesa_eos_dh
 enddo

 fmt1 = '(F4.2)' !

 fnum=122

 ! Convert first X and Z to string
 write (zz,fmt1) mesa_eos_z1
 write (hh,fmt1) mesa_eos_h1

 ! Find the first EoS tables
 filename = trim(mesa_eos_prefix)//'z'//trim(zz)//'x'//trim(hh)//'.bindata'
 filename = find_phantom_datafile(filename,'eos/mesa')

 ! Read constants from the header of first EoS tables
 open(unit=fnum, file=trim(filename), status='old', action='read', form='unformatted',iostat=ierr)
 if (ierr /= 0) return
 read(fnum) mesa_eos_ne, mesa_eos_nv, mesa_eos_nvar2
 close(fnum)

 ! Allocate ALL the arrays!
 allocate(mesa_eos_logEs(mesa_eos_ne), mesa_eos_logVs(mesa_eos_nv))
 allocate(mesa_eos_data_exists(mesa_eos_nz,mesa_eos_nh))
 allocate(mesa_de_data(mesa_eos_ne,mesa_eos_nv,mesa_eos_nvar2))
 allocate(mesa_eos0(mesa_eos_nz,mesa_eos_nh,mesa_eos_ne,mesa_eos_nv))
 allocate(mesa_de_data0(mesa_eos_nz,mesa_eos_nh,mesa_eos_ne,mesa_eos_nv,mesa_eos_nvar2))

 return
end subroutine get_eos_constants_mesa


! Read MESA EoS tables, and then construct a new array for the specific values of X and Z
subroutine read_eos_mesa(x,z,ierr)
 real,    intent(in)  :: x, z
 integer, intent(out) :: ierr
 real, parameter      :: arad=7.5657d-15
 integer              :: i,j,k,l,m
 character (len=300)  :: filename
 character (len=8)    :: fmt1
 character (len=20)   :: zz, hh
 integer              :: nz1,nz2,nx1,nx2
 real                 :: dz,dx
 integer              :: fnum

 fmt1 = '(F4.2)'

 fnum = 120
 ierr = 0

 ! following lines to prevent compiler warnings
 nx1 = 0
 nx2 = 0
 dx = 0.

 !Loop over files
 do i=1,mesa_eos_nz
    write (zz,fmt1) mesa_eos_z(i) ! Z value as string

    do j=1,mesa_eos_nh
       write (hh,fmt1) mesa_eos_h(j) ! X value as string

       ! Find the EoS tables
       filename = trim(mesa_eos_prefix)//'z'//trim(zz)//'x'//trim(hh)//'.bindata'
       filename = find_phantom_datafile(filename,'eos/mesa')

       ! Note that the data exists
       mesa_eos_data_exists(i,j)=1

       ! Read in the size of the tables and the data
       ! i and j hold the Z and X values respectively
       ! k, l and m hold the values of V, Eint and the data respectively
       open(unit=fnum, file=trim(filename), status='old', action='read', form='unformatted')
       read(fnum) mesa_eos_ne, mesa_eos_nv, mesa_eos_nvar2
       read(fnum)(mesa_eos_logVs(k),k=1,mesa_eos_nv)
       read(fnum)(mesa_eos_logEs(l),l=1,mesa_eos_ne)
       do k=1,mesa_eos_nv
          do l=1,mesa_eos_ne
             read(fnum) (mesa_de_data0(i,j,l,k,m),m=1,mesa_eos_nvar2)
          enddo
       enddo

       close(fnum)
    enddo
 enddo

 ! Chooses the Z values from the available values which is greater than the input Z
 nz2=mesa_eos_nz
 do i=2,mesa_eos_nz-1
    if(mesa_eos_z(i) >= z) then
       nz2=i
       exit
    endif
 enddo

 ! Chooses Z values below input Z, and then finds fraction of input Z distance between two enclosing Z (sorry if confusing, just read the code)
 nz1 = nz2 - 1
 dz=(z-mesa_eos_z(nz1))/(mesa_eos_z(nz2)-mesa_eos_z(nz1))

 ! Seems to do the same as above, except for X. I can't quite figure out the weird if statements.
 do j=2,mesa_eos_nh
    if(mesa_eos_h(j) > x.or.j==mesa_eos_nh) then
       nx2=j
       if(j==2.and.mesa_eos_data_exists(i,j-2)==0) then
          nx2=3
       else if(j==mesa_eos_nh.and.mesa_eos_data_exists(i,mesa_eos_nh-1)==0) then
          nx2=mesa_eos_nh-1
       endif
       exit
    endif
 enddo
 nx1 = nx2 - 1
 dx = (x-mesa_eos_h(nx1))/(mesa_eos_h(nx2)-mesa_eos_h(nx1))

 ! Bilinear interpolation of values from the four tables (two values each of X and Z)
 do l=1,mesa_eos_ne
    do k=1,mesa_eos_nv
       do m=1,mesa_eos_nvar2
          mesa_de_data(l,k,m)=(1.d0-dx)*(1.d0-dz)*mesa_de_data0(nz1,nx1,l,k,m)+dz*(1.d0-dx)*mesa_de_data0(nz2,nx1,l,k,m) &
               +dx*(1.d0-dz)*mesa_de_data0(nz1,nx2,l,k,m)+dx*dz*mesa_de_data0(nz2,nx2,l,k,m)
       enddo
    enddo
 enddo


 ! Save some things to make the interpolation fast
 mesa_eos_dv=mesa_eos_logVs(2)-mesa_eos_logVs(1)
 mesa_eos_de=mesa_eos_logEs(2)-mesa_eos_logEs(1)
 mesa_eos_v1=mesa_eos_logVs(1)
 mesa_eos_e1=mesa_eos_logEs(1)

 !deallocate(mesa_eos0)
 !deallocate(mesa_de_data0)

 return

end subroutine read_eos_mesa



! Get value from the MESA tables given a specific value of density and internal energy
! The columns in the data are:
! 1. logRho        2. logP          3. logPgas       4. logT
! 5. dlnP/dlnrho|e 6. dlnP/dlne|rho 7. dlnT/dlnrho|e 8. dlnT/dlne|rho
! 9. logS         10. dlnT/dlnP|S  11. Gamma1       12. gamma
subroutine getvalue_mesa(rho,eint,ivout,vout,ierr)
 real, intent(in) :: rho, eint
 real, intent(out) :: vout
 integer, intent(in) :: ivout
 integer, intent(out), optional :: ierr
 real :: loge, logv, de, dv
 integer :: ne, nv
 real :: dx
 integer :: nx


 ! logRho = logV + 0.7*logE - 20

 loge = log10(eint)
 logv = 20.d0 + log10(rho) - 0.7d0 * loge

 ! Get the Eint and V indices for looking up the tables
 ne = 1 + int((loge - mesa_eos_e1) / mesa_eos_de)
 nv = 1 + int((logv - mesa_eos_v1) / mesa_eos_dv)

 ! Allow extrapolation
 if (ne < 1) ne=1
 if (ne > mesa_eos_ne-1) ne=mesa_eos_ne-1
 if (nv < 1) nv=1
 if (nv > mesa_eos_nv-1) nv=mesa_eos_nv-1

 de = (loge - mesa_eos_logEs(ne)) / mesa_eos_de
 dv = (logv - mesa_eos_logVs(nv)) / mesa_eos_dv

 ! If the given Eint and V fall within the limits of the table, then use cubic spline interpolation to find value
 ! Else use linear extrapolation beyond the limits of the table (I think? I can't tell if this is the correct way to do this)
 if (ne > 1 .and. nv > 1 .and. ne < mesa_eos_ne-1 .and. nv < mesa_eos_nv-1) then
    call eos_cubic_spline_mesa(ne,nv,loge,logv,ivout,vout,nx,dx)
    if (ivout < 5) vout=10.d0**vout
    if (present(ierr)) ierr = 0
 else
    vout = 10.d0**((1.d0-de) * (1.d0-dv) * mesa_de_data(ne,nv,ivout)   + &
                         de  * (1.d0-dv) * mesa_de_data(ne+1,nv,ivout) + &
                   (1.d0-de) *       dv  * mesa_de_data(ne,nv+1,ivout) + &
                         de  *       dv  * mesa_de_data(ne+1,nv+1,ivout))
    if (present(ierr)) ierr = 1  ! warn if extrapolating
 endif

 return

end subroutine getvalue_mesa


!only use if between e(2) < e < e(n_e-1) and v(2) < v < v(n_v-1)

subroutine  eos_cubic_spline_mesa(e1,v1,e,v,n_var,z,h1,dh)

! use mesa_eos_logEs, mesa_eos_logVs
 implicit none

 integer, intent(in) :: e1, v1, h1, n_var
 real, intent(in) :: e, v, dh

 real, intent(out) :: z

 real :: x0,x1,x2,x3,y0,y1,y2,y3
 real :: as,bs,tt

 real :: yy(4)
 integer :: j

 x0=mesa_eos_logEs(e1-1)
 x1=mesa_eos_logEs(e1+0)
 x2=mesa_eos_logEs(e1+1)
 x3=mesa_eos_logEs(e1+2)

 tt=(e-x1)/mesa_eos_de

 do j=1,4

    y0 = mesa_de_data(e1-1,v1+j-2,n_var)
    y1 = mesa_de_data(e1+0,v1+j-2,n_var)
    y2 = mesa_de_data(e1+1,v1+j-2,n_var)
    y3 = mesa_de_data(e1+2,v1+j-2,n_var)

    call cubic_spline_mesa(x0,x1,x2,x3,y0,y1,y2,y3,as,bs)
    yy(j)=(1.d0-tt)*y1+tt*y2+tt*(1.d0-tt)*(as*(1.d0-tt)+bs*tt)

 enddo

 x0=mesa_eos_logVs(v1-1)
 x1=mesa_eos_logVs(v1+0)
 x2=mesa_eos_logVs(v1+1)
 x3=mesa_eos_logVs(v1+2)

 y0=yy(1)
 y1=yy(2)
 y2=yy(3)
 y3=yy(4)

 tt=(v-x1)/mesa_eos_dv

 call cubic_spline_mesa(x0,x1,x2,x3,y0,y1,y2,y3,as,bs)
 z=(1.d0-tt)*y1+tt*y2+tt*(1.d0-tt)*(as*(1.d0-tt)+bs*tt)

 return

end subroutine eos_cubic_spline_mesa



subroutine cubic_spline_mesa(x0,x1,x2,x3,y0,y1,y2,y3,as,bs)

 implicit none

 real, intent(in) :: x0,x1,x2,x3,y0,y1,y2,y3
 real :: k1,k2  !k1 and k2 are derivatives
 real, intent(out) :: as,bs

 k2=(y3-y1)/(x3-x1)
 k1=(y2-y0)/(x2-x0)

 as=k1*(x2-x1)-(y2-y1)
 bs=-k2*(x2-x1)+(y2-y1)

 return
end subroutine cubic_spline_mesa

subroutine deallocate_arrays_mesa

 if (allocated(mesa_opacs_zs)) deallocate(mesa_opacs_zs)
 if (allocated(mesa_opacs_xs)) deallocate(mesa_opacs_xs)
 if (allocated(mesa_opacs_rs)) deallocate(mesa_opacs_rs)
 if (allocated(mesa_opacs_ts)) deallocate(mesa_opacs_ts)

 if (allocated(mesa_opacs_k)) deallocate(mesa_opacs_k)
 if (allocated(mesa_opacs_kd)) deallocate(mesa_opacs_kd)
 if (allocated(mesa_opacs_kt)) deallocate(mesa_opacs_kt)

 if (allocated(mesa_eos_z)) deallocate(mesa_eos_z)
 if (allocated(mesa_eos_h)) deallocate(mesa_eos_h)
 if (allocated(mesa_eos_logEs)) deallocate(mesa_eos_logEs)
 if (allocated(mesa_eos_logVs)) deallocate(mesa_eos_logVs)
 if (allocated(mesa_eos_data_exists)) deallocate(mesa_eos_data_exists)

 if (allocated(mesa_eos0)) deallocate(mesa_eos0)
 if (allocated(mesa_de_data)) deallocate(mesa_de_data)
 if (allocated(mesa_de_data0)) deallocate(mesa_de_data0)

end subroutine deallocate_arrays_mesa

end module mesa_microphysics
