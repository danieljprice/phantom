!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric
!
! Superposed Kerr-Schild metric for a binary black hole system
!
! :References:
!   Combi & Ressler (2024) arXiv:2403.13308
!   Original file from: https://zenodo.org/records/10841021
!   translated from C into Fortran by Daniel Price (2026)
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - a1x             : *spin in x dir for black hole 1*
!   - a1y             : *spin in y dir for black hole 1*
!   - a1z             : *spin in z dir for black hole 1*
!   - a2x             : *spin in x dir for black hole 2*
!   - a2y             : *spin in y dir for black hole 2*
!   - a2z             : *spin in z dir for black hole 2*
!   - mass1           : *mass of black hole 1*
!   - mass2           : *mass of black hole 2*
!   - trajectory_file : *file containing binary black hole trajectory*
!
! :Dependencies: datafiles, dump_utils, infile_utils, inverse4x4, io
!
 implicit none
 character(len=*), parameter :: metric_type = 'binarybh'
 integer,          parameter :: imetric = 5

 ! extern void SuperposedBbh(const double *xx, double gcov[][NDIM], const double *traj_array)
 integer, parameter :: nparams = 12
 real, public :: metric_params(nparams)

 ! parameters from the original SuperposedBbh.c file, could in principle be user options
 real, parameter :: tiny_param = 1e-40
 real, parameter :: AST_adjust_mass1 = 1.
 real, parameter :: AST_adjust_mass2 = 1.
 real, parameter :: AST_a1_buffer = 1.e-4
 real, parameter :: AST_a2_buffer = 1.e-4
 real, parameter :: AST_cutoff_floor = 0.1

 real, private :: a1x = 0.,a1y = 0.,a1z = 0.,a2x = 0.,a2y = 0.,a2z = 0.
 real, public :: mass1 = 0.5
 real, public :: mass2 = 0.5
 real, public  :: a = 0.       ! black hole 1 spin
 real, private :: a_bh2 = 0.   ! black hole 2 spin
 character(len=128), public :: trajectory_file = 'cbwaves.txt'
 logical, private :: metric_initialised = .false.

contains

!-------------------------------------------------------------------------------
!+
!  Subroutine to update the metric inputs if time dependent
!+
!-------------------------------------------------------------------------------
subroutine update_metric(time)
 use io,        only:id,master,fatal
 use datafiles, only:find_phantom_datafile
 real, intent(in) :: time
 real :: x1(3),x2(3),v1(3),v2(3)
 integer :: ierr

 ! defaults for a single black hole at the origin, for testing
 x1 = [0.,0.,0.]
 x2 = [10000000.,0.,0.]
 v1 = [0.,0.,0.]
 v2 = [0.,0.,0.]
 !print*,' updating metric at time ',time

 if (.not.metric_initialised) then
    metric_initialised = .true.
    trajectory_file = find_phantom_datafile(trajectory_file,'binarybh')
    if (id==master) print "(a)",' Reading black hole trajectory from '//trim(trajectory_file)
 endif
 call get_trajectory_from_file(time,x1,x2,v1,v2,ierr)
 if (ierr /= 0) call fatal('metric_binarybh','could not open trajectory file '//trim(trajectory_file))
 !if (id==master) print*,' time = ',time,' binary separation = ',sqrt(dot_product(x1 - x2,x1 - x2))

 metric_params(1:3) = x1
 metric_params(4:6) = x2
 metric_params(7:9) = v1
 metric_params(10:12) = v2

end subroutine update_metric

!----------------------------------------------------------------
!+
!  Read the binary black hole trajectory from file with linear
!  interpolation in time to match the specified time.
!+
!----------------------------------------------------------------
subroutine get_trajectory_from_file(time,x1,x2,v1,v2,ierr)
 use io, only:error
 real,    intent(in)  :: time
 real,    intent(out) :: x1(3),x2(3),v1(3),v2(3)
 integer, intent(out) :: ierr
 integer :: iu,nlines
 real :: t_prev,t_next
 real :: x1_prev(3),x2_prev(3),v1_prev(3),v2_prev(3)
 real :: x1_next(3),x2_next(3),v1_next(3),v2_next(3)
 real :: frac

 t_prev = 0.
 x1 = [0.,0.,0.]
 x2 = [0.,0.,0.]
 v1 = [0.,0.,0.]
 v2 = [0.,0.,0.]

 open(newunit=iu,file=trajectory_file,status='old',action='read',iostat=ierr)
 if (ierr /= 0) then
    call error('metric_binarybh','could not open trajectory file '//trim(trajectory_file))
    return
 endif
 read(iu,*,iostat=ierr) nlines
 if (ierr /= 0) then
    call error('metric_binarybh','could not read nlines from '//trim(trajectory_file))
 endif
 read(iu,*,iostat=ierr) t_prev,x1_prev,x2_prev,v1_prev,v2_prev
 if (ierr /= 0) then
    call error('metric_binarybh','could not read first trajectory line from '//trim(trajectory_file))
    return
 endif
 if (time <= t_prev) then
    x1 = x1_prev
    x2 = x2_prev
    v1 = v1_prev
    v2 = v2_prev
    close(iu)
    return
 endif
 do
    read(iu,*,iostat=ierr) t_next,x1_next,x2_next,v1_next,v2_next
    if (ierr /= 0) then
       call error('metric_binarybh','time beyond range of trajectory file: '//trim(trajectory_file))
    endif
    if (time <= t_next) then
       frac = (time - t_prev) / (t_next - t_prev)
       x1 = x1_prev + frac * (x1_next - x1_prev)
       x2 = x2_prev + frac * (x2_next - x2_prev)
       v1 = v1_prev + frac * (v1_next - v1_prev)
       v2 = v2_prev + frac * (v2_next - v2_prev)
       close(iu)
       return
    endif
    t_prev = t_next
    x1_prev = x1_next
    x2_prev = x2_next
    v1_prev = v1_next
    v2_prev = v2_next
 enddo

end subroutine get_trajectory_from_file

!----------------------------------------------------------------
!+
!  Compute the metric tensor in both covariant (gcov) and
!  contravariant (gcon) form
!+
!----------------------------------------------------------------
pure subroutine get_metric_cartesian(xx,gcov,gcon,sqrtg)
 use inverse4x4, only:inv4x4
 real, intent(in)  :: xx(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: det
 real :: x,y,z
 real :: ks1(0:3,0:3), ks2(0:3,0:3), j1(0:3,0:3), j2(0:3,0:3)
 real :: xi1x, xi1y, xi1z, xi2x, xi2y, xi2z
 real :: v1x, v1y, v1z, v2x, v2y, v2z
 real :: m1, m2, m1_t, m2_t
 real :: a1, a2, a1_t, a2_t
 real :: o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12
 real :: o13, o14, o15, o16, o17, o18, o19, o20, o21, o22, o23, o24
 real :: o25, o26, o27, o28, o29, o30, o31, o32, o33, o34, o35, o36
 real :: o37, o38, o39, o40, o41, o42, o43, o44, o45, o46, o47, o48
 real :: o49, o50, o51, o52, o53, o54, o55, o56, o57, o58, o59, o60
 real :: o61, o62, o63, o64, o65, o66, o67, o68, o69, o70, o71, o72
 real :: o73, o74, o75, o76, o77, o78, o79, o80, o81, o82, o83, o84
 real :: o85, o86, o87, o88, o89, o90, o91, o92, o93, o94, o95, o96
 real :: o97, o98, o99, o100, o101, o102, o103, o104, o105, o106, o107, o108
 real :: o109, o110, o111, o112, o113, o114, o115, o116, o117, o118, o119, o120
 real :: o121, o122, o123, o124, o125, o126, o127, o128, o129, o130, o131, o132
 real :: oo1, oo2, oo3, oo4, oo5, oo6, oo7, oo8, oo9, oo10, oo11, oo12
 real :: oo13, oo14, oo15, oo16, oo17, oo18, oo19, oo20, oo21, oo22, oo23, oo24
 real :: oo25, oo26, oo27, oo28, oo29, oo30, oo31, oo32, oo33, oo34, oo35, oo36, oo37
 real :: rbh1, rbh2, rbh1_Cutoff, rbh2_Cutoff, sum, term1, term2, term3, term4
 real :: v1, v2, x0bh1, x0bh2, x1bh1, x1bh2, x2bh1, x2bh2, x3bh1, x3bh2
 integer :: i, j, m, n

 ! code below is translated from the original SuperposedBBH.c file

 ! mask
 x = xx(1)
 y = xx(2)
 z = xx(3)

 ! Load trajectories
 xi1x = metric_params(1)
 xi1y = metric_params(2)
 xi1z = metric_params(3)
 xi2x = metric_params(4)
 xi2y = metric_params(5)
 xi2z = metric_params(6)
 v1x = metric_params(7) + tiny_param
 v1y = metric_params(8) + tiny_param
 v1z = metric_params(9) + tiny_param
 v2x = metric_params(10) + tiny_param
 v2y = metric_params(11) + tiny_param
 v2z = metric_params(12) + tiny_param

 v2 = sqrt( v2x * v2x + v2y * v2y + v2z * v2z )
 v1 = sqrt( v1x * v1x + v1y * v1y + v1z * v1z )

 m1_t = mass1
 m2_t = mass2
 a1_t = a     !sqrt(a1x*a1x + a1y*a1y + a1z*a1z + tiny_param)
 a2_t = a_bh2 !sqrt(a2x*a2x + a2y*a2y + a2z*a2z + tiny_param)

 ! Load coordinates
 oo1 = v1 * v1
 oo2 = -oo1
 oo3 = 1. + oo2
 oo4 = sqrt(oo3)
 oo5 = 1. / oo4
 oo6 = -x
 oo7 = oo6 + xi1x
 oo8 = v1x * oo7
 oo9 = -y
 oo10 = -z
 oo11 = v2 * v2
 oo12 = -oo11
 oo13 = 1. + oo12
 oo14 = sqrt(oo13)
 oo15 = 1. / oo14
 oo16 = oo6 + xi2x
 oo17 = v2x * oo16
 oo18 = -xi1x
 oo19 = 1. / oo1
 oo20 = -1. + oo4
 oo21 = -xi1y
 oo22 = -xi1z
 oo23 = -xi2x
 oo24 = 1. / oo11
 oo25 = -1. + oo14
 oo26 = -xi2y
 oo27 = -xi2z
 oo28 = xi1y * v1y
 oo29 = xi1z * v1z
 oo30 = v1y * (-y)
 oo31 = v1z * (-z)
 oo32 = oo28 + (oo29 + (oo30 + (oo31 + oo8)))
 oo33 = xi2y * v2y
 oo34 = xi2z * v2z
 oo35 = v2y * (-y)
 oo36 = v2z * (-z)
 oo37 = oo17 + (oo33 + (oo34 + (oo35 + oo36)))
 x0bh1 = (oo8 + ((oo9 + xi1y) * v1y + (oo10 + xi1z) * v1z)) * oo5
 x0bh2 = (oo17 + ((oo9 + xi2y) * v2y + (oo10 + xi2z) * v2z)) * oo15
 x1bh1 = (oo18 + x) - oo20 * (oo5 * (v1x * (((oo18 + x) * v1x + ((oo21 + y) * v1y + (oo22 + z) * v1z)) * oo19)))
 x1bh2 = (oo23 + x) - oo24 * (oo25 * (v2x * (((oo23 + x) * v2x + ((oo26 + y) * v2y + (oo27 + z) * v2z)) * oo15)))
 x2bh1 = oo21 + (oo20 * (oo32 * (oo5 * (v1y * oo19))) + y)
 x2bh2 = oo26 + (oo24 * (oo25 * (oo37 * (v2y * oo15))) + y)
 x3bh1 = oo22 + (oo20 * (oo32 * (oo5 * (v1z * oo19))) + z)
 x3bh2 = oo27 + (oo24 * (oo25 * (oo37 * (v2z * oo15))) + z)

 ! Adjust mass
 ! This is useful for reducing the effective mass of each bh
 ! Adjust by hand to get the correct irreducible mass of the bh
 a1 = a1_t * AST_adjust_mass1
 m1 = m1_t * AST_adjust_mass1
 a2 = a2_t * AST_adjust_mass2
 m2 = m2_t * AST_adjust_mass2

 !============================================
 ! Regularize horizon and apply excision mask
 !============================================

 ! Define radius with respect to bh frame
 rbh1 = sqrt(x1bh1*x1bh1 + x2bh1*x2bh1 + x3bh1*x3bh1)
 rbh2 = sqrt(x1bh2*x1bh2 + x2bh2*x2bh2 + x3bh2*x3bh2)

 ! Define radius cutoff
 rbh1_Cutoff = abs(a1) * (1.0 + AST_a1_buffer) + AST_cutoff_floor
 rbh2_Cutoff = abs(a2) * (1.0 + AST_a2_buffer) + AST_cutoff_floor

 ! Apply excision
 if ((rbh1) < rbh1_Cutoff) then
    if (x3bh1>0) then
       x3bh1 = rbh1_Cutoff
    else
       x3bh1 = -rbh1_Cutoff
    endif
 endif
 if ((rbh2) < rbh2_Cutoff) then
    if (x3bh2>0) then
       x3bh2 = rbh2_Cutoff
    else
       x3bh2 = -rbh2_Cutoff
    endif
 endif

 !=================
 !     Metric
 !=================
 o1 = 1.4142135623730951
 o2 = 1. / o1
 o3 = a1x * a1x
 o4 = -o3
 o5 = a1z * a1z
 o6 = -o5
 o7 = a2x * a2x
 o8 = -o7
 o9 = x1bh1 * x1bh1
 o10 = x2bh1 * x2bh1
 o11 = x3bh1 * x3bh1
 o12 = x1bh1 * a1x
 o13 = x2bh1 * a2x
 o14 = x3bh1 * a1z
 o15 = o12 + (o13 + o14)
 o16 = o15 * o15
 o17 = o16 * 4.
 o18 = o10 + (o11 + (o4 + (o6 + (o8 + o9))))
 o19 = o18 * o18
 o20 = o17 + o19
 o21 = sqrt(o20)
 o22 = o10 + (o11 + (o21 + (o4 + (o6 + (o8 + o9)))))
 o23 = o22**1.5
 o24 = o22 * o22
 o25 = o24 * 0.25
 o26 = o16 + o25
 o27 = 1. / o26
 o28 = x2bh1 * a1z
 o29 = a2x * (-x3bh1)
 o30 = sqrt(o22)
 o31 = 1. / o30
 o32 = o1 * (o15 * (o31 * a1x))
 o33 = o30 * (x1bh1 * o2)
 o34 = o28 + (o29 + (o32 + o33))
 o35 = o22 * 0.5
 o36 = o3 + (o35 + (o5 + o7))
 o37 = 1. / o36
 o38 = o2 * (o23 * (o27 * (o34 * (o37 * m1))))
 o39 = a1z * (-x1bh1)
 o40 = x3bh1 * a1x
 o41 = o1 * (o15 * (o31 * a2x))
 o42 = o30 * (x2bh1 * o2)
 o43 = o39 + (o40 + (o41 + o42))
 o44 = o2 * (o23 * (o27 * (o37 * (o43 * m1))))
 o45 = x1bh1 * a2x
 o46 = a1x * (-x2bh1)
 o47 = o1 * (o15 * (o31 * a1z))
 o48 = o30 * (x3bh1 * o2)
 o49 = o45 + (o46 + (o47 + o48))
 o50 = o2 * (o23 * (o27 * (o37 * (o49 * m1))))
 o51 = o36 * o36
 o52 = 1. / o51
 o53 = o2 * (o23 * (o27 * (o34 * (o43 * (o52 * m1)))))
 o54 = o2 * (o23 * (o27 * (o34 * (o49 * (o52 * m1)))))
 o55 = o2 * (o23 * (o27 * (o43 * (o49 * (o52 * m1)))))
 o56 = a2y * a2y
 o57 = -o56
 o58 = a2z * a2z
 o59 = -o58
 o60 = x1bh2 * x1bh2
 o61 = x2bh2 * x2bh2
 o62 = x3bh2 * x3bh2
 o63 = x1bh2 * a2x
 o64 = x2bh2 * a2y
 o65 = x3bh2 * a2z
 o66 = o63 + (o64 + o65)
 o67 = o66 * o66
 o68 = o67 * 4.
 o69 = o57 + (o59 + (o60 + (o61 + (o62 + o8))))
 o70 = o69 * o69
 o71 = o68 + o70
 o72 = sqrt(o71)
 o73 = o57 + (o59 + (o60 + (o61 + (o62 + (o72 + o8)))))
 o74 = o73**1.5
 o75 = o73 * o73
 o76 = o75 * 0.25
 o77 = o67 + o76
 o78 = 1. / o77
 o79 = x2bh2 * a2z
 o80 = a2y * (-x3bh2)
 o81 = sqrt(o73)
 o82 = 1. / o81
 o83 = o1 * (o66 * (o82 * a2x))
 o84 = o81 * (x1bh2 * o2)
 o85 = o79 + (o80 + (o83 + o84))
 o86 = o73 * 0.5
 o87 = o56 + (o58 + (o7 + o86))
 o88 = 1. / o87
 o89 = o2 * (o74 * (o78 * (o85 * (o88 * m2))))
 o90 = a2z * (-x1bh2)
 o91 = x3bh2 * a2x
 o92 = o1 * (o66 * (o82 * a2y))
 o93 = o81 * (x2bh2 * o2)
 o94 = o90 + (o91 + (o92 + o93))
 o95 = o2 * (o74 * (o78 * (o88 * (o94 * m2))))
 o96 = x1bh2 * a2y
 o97 = a2x * (-x2bh2)
 o98 = o1 * (o66 * (o82 * a2z))
 o99 = o81 * (x3bh2 * o2)
 o100 = o96 + (o97 + (o98 + o99))
 o101 = o100 * (o2 * (o74 * (o78 * (o88 * m2))))
 o102 = o87 * o87
 o103 = 1. / o102
 o104 = o103 * (o2 * (o74 * (o78 * (o85 * (o94 * m2)))))
 o105 = o100 * (o103 * (o2 * (o74 * (o78 * (o85 * m2)))))
 o106 = o100 * (o103 * (o2 * (o74 * (o78 * (o94 * m2)))))
 o107 = v1 * v1
 o108 = -o107
 o109 = 1. + o108
 o110 = sqrt(o109)
 o111 = 1. / o110
 o112 = o111 * (-v1x)
 o113 = o111 * (-v1y)
 o114 = o111 * (-v1z)
 o115 = 1. / o107
 o116 = -1. + o111
 o117 = o116 * (v1x * (v1y * o115))
 o118 = o116 * (v1x * (v1z * o115))
 o119 = o116 * (v1y * (v1z * o115))
 o120 = v2 * v2
 o121 = -o120
 o122 = 1. + o121
 o123 = sqrt(o122)
 o124 = 1. / o123
 o125 = o124 * (-v2x)
 o126 = o124 * (-v2y)
 o127 = o124 * (-v2z)
 o128 = 1. / o120
 o129 = -1. + o124
 o130 = o129 * (v2x * (v2y * o128))
 o131 = o129 * (v2x * (v2z * o128))
 o132 = o129 * (v2y * (v2z * o128))
 ks1(0,0) = o2 * (o23 * (o27 * m1))
 ks1(1,0) = o38
 ks1(2,0) = o44
 ks1(3,0) = o50
 ks1(0,1) = o38
 ks1(1,1) = o2 * (o23 * (o27 * ((o34 * o34) * (o52 * m1))))
 ks1(2,1) = o53
 ks1(3,1) = o54
 ks1(0,2) = o44
 ks1(1,2) = o53
 ks1(2,2) = o2 * (o23 * (o27 * ((o43 * o43) * (o52 * m1))))
 ks1(3,2) = o55
 ks1(0,3) = o50
 ks1(1,3) = o54
 ks1(2,3) = o55
 ks1(3,3) = o2 * (o23 * (o27 * ((o49 * o49) * (o52 * m1))))
 ks2(0,0) = o2 * (o74 * (o78 * m2))
 ks2(1,0) = o89
 ks2(2,0) = o95
 ks2(3,0) = o101
 ks2(0,1) = o89
 ks2(1,1) = o103 * (o2 * (o74 * (o78 * ((o85 * o85) * m2))))
 ks2(2,1) = o104
 ks2(3,1) = o105
 ks2(0,2) = o95
 ks2(1,2) = o104
 ks2(2,2) = o103 * (o2 * (o74 * (o78 * ((o94 * o94) * m2))))
 ks2(3,2) = o106
 ks2(0,3) = o101
 ks2(1,3) = o105
 ks2(2,3) = o106
 ks2(3,3) = (o100 * o100) * (o103 * (o2 * (o74 * (o78 * m2))))
 j1(0,0) = o111
 j1(1,0) = o112
 j1(2,0) = o113
 j1(3,0) = o114
 j1(0,1) = o112
 j1(1,1) = 1. + o116 * ((v1x * v1x) * o115)
 j1(2,1) = o117
 j1(3,1) = o118
 j1(0,2) = o113
 j1(1,2) = o117
 j1(2,2) = 1. + o116 * ((v1y * v1y) * o115)
 j1(3,2) = o119
 j1(0,3) = o114
 j1(1,3) = o118
 j1(2,3) = o119
 j1(3,3) = 1. + o116 * ((v1z * v1z) * o115)
 j2(0,0) = o124
 j2(1,0) = o125
 j2(2,0) = o126
 j2(3,0) = o127
 j2(0,1) = o125
 j2(1,1) = 1. + o129 * ((v2x * v2x) * o128)
 j2(2,1) = o130
 j2(3,1) = o131
 j2(0,2) = o126
 j2(1,2) = o130
 j2(2,2) = 1. + o129 * ((v2y * v2y) * o128)
 j2(3,2) = o132
 j2(0,3) = o127
 j2(1,3) = o131
 j2(2,3) = o132
 j2(3,3) = 1. + o129 * ((v2z * v2z) * o128)

 ! Initialize the flat part
 gcov = reshape([-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],[4,4])

 ! Load symmetric gcov
 do i=0,3
    do j=i,3
       sum = 0.
       do m=0,3
          term1 = j2(i,m)
          term2 = j1(i,m)
          do n=0,3
             term3 = j2(j,n)
             term4 = j1(j,n)
             sum = sum + (term1 * term3 * ks2(n,m) + term2 * term4 * ks1(n,m))
          enddo
       enddo
       gcov(i,j) = gcov(i,j) + sum
       gcov(j,i) = gcov(i,j)
    enddo
 enddo

 ! end of original SuperposedBbh.c file contents

 if (present(gcon)) then
    gcon = 0.
    call inv4x4(gcov,gcon,det)
    if (present(sqrtg)) sqrtg = sqrt(-det)
 endif

end subroutine get_metric_cartesian

pure subroutine metric_cartesian_derivatives(xx,dgcovdx,dgcovdy,dgcovdz)
 real,    intent(in)  :: xx(3)
 real,    intent(out) :: dgcovdx(0:3,0:3),dgcovdy(0:3,0:3),dgcovdz(0:3,0:3)
 real, parameter :: eps = 1.e-10
 real, parameter :: two_eps = 2.*eps
 real :: gcov1(0:3,0:3), gcov2(0:3,0:3)

 dgcovdx = 0.
 dgcovdy = 0.
 dgcovdz = 0.

 ! metric has no time dependence centred differencing: (f(x+h)-f(x-h))/(2h)
 ! x direction
 call get_metric_cartesian(xx+[eps,0.,0.],gcov1)
 call get_metric_cartesian(xx+[-eps,0.,0.],gcov2)
 dgcovdx = (gcov1 - gcov2) / two_eps

 ! y direction
 call get_metric_cartesian(xx+[0.,eps,0.],gcov1)
 call get_metric_cartesian(xx+[0.,-eps,0.],gcov2)
 dgcovdy = (gcov1 - gcov2) / two_eps

 ! z direction
 call get_metric_cartesian(xx+[0.,0.,eps],gcov1)
 call get_metric_cartesian(xx+[0.,0.,-eps],gcov2)
 dgcovdz = (gcov1 - gcov2) / two_eps

end subroutine metric_cartesian_derivatives

!----------------------------------------------------------------
!+
!  The metric tensor in SPHERICAL-like form
!  (these are dummy routines for compatibility with other metrics)
!+
!----------------------------------------------------------------
pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg

 gcov = 0.
 if (present(gcon)) gcon = 0.
 if (present(sqrtg)) sqrtg = 0.

end subroutine get_metric_spherical

pure subroutine cartesian2spherical(xcart,xspher)
 real, intent(in)  :: xcart(3)
 real, intent(out) :: xspher(3)

 xspher = xcart

end subroutine cartesian2spherical

pure subroutine metric_spherical_derivatives(position,dgcovdr,dgcovdtheta,dgcovdphi)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi

 dgcovdr     = 0.
 dgcovdtheta = 0.
 dgcovdphi   = 0.

end subroutine metric_spherical_derivatives

!----------------------------------------------------------------
!+
!  Check if a particle should be accreted by either black hole
!+
!----------------------------------------------------------------
subroutine accrete_particles_metric(xi,yi,zi,mi,ti,accradius1,accreted)
 real,    intent(in)  :: xi,yi,zi,mi,ti,accradius1
 logical, intent(out) :: accreted
 real :: r1,r2,x1(3),x2(3)

 x1 = metric_params(1:3)
 x2 = metric_params(4:6)

 r1 = ((xi-x1(1))**2 + (yi-x1(2))**2 + (zi-x1(3))**2)/mass1**2
 r2 = ((xi-x2(1))**2 + (yi-x2(2))**2 + (zi-x2(3))**2)/mass2**2
 if (r1 < accradius1**2 .or. r2 < accradius1**2) then
    accreted = .true.
 else
    accreted = .false.
 endif

end subroutine accrete_particles_metric

!-----------------------------------------------------------------------
!+
!  writes relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_metric(hdr,time,accradius,ierr)
 use dump_utils, only:lentag,dump_h,add_to_rheader
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time,accradius
 integer,      intent(out)   :: ierr
 character(len=lentag)       :: tags(16)
 real    :: rheader(16)
 integer :: i

 ierr = 0
 rheader(1:3) = metric_params(1:3)
 rheader(4) = mass1
 rheader(5) = accradius*mass1
 rheader(6:8) = metric_params(4:6)
 rheader(9) = mass2
 rheader(10) = accradius*mass2
 rheader(11:13) = metric_params(7:9)
 rheader(14:16) = metric_params(10:12)

 !  rheader(17) = accretedmass1
 !  rheader(18) = accretedmass2

 tags(1:16) = (/'x1 ','y1 ','z1 ','m1 ','h1 ','x2 ','y2 ','z2 ','m2 ', &
                'h2 ','vx1','vy1','vz1','vx2','vy2','vz2'/)

 do i=1,16
    call add_to_rheader(rheader(i),tags(i),hdr,ierr)
 enddo

end subroutine write_headeropts_metric

!-----------------------------------------------------------------------
!+
!  reads relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_metric(hdr,ierr)
 use dump_utils, only:dump_h,extract
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr

 ierr  = 0
! call extract('accretedmass1',accretedmass1,hdr,ierr1)
! call extract('accretedmass2',accretedmass2,hdr,ierr2)

! if (ierr1 /= 0 .or. ierr2 /= 0) then
 !   write(*,*) ' ERROR extracting accretedmass1 and accretedmass2 from file'
 !   ierr = 1
 !endif

end subroutine read_headeropts_metric

!-----------------------------------------------------------------------
!+
!  writes metric options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_metric(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# Options relating to the '//trim(metric_type)//' metric'
 call write_inopt(a1x,'a1x','spin in x dir for black hole 1',iunit)
 call write_inopt(a1y,'a1y','spin in y dir for black hole 1',iunit)
 call write_inopt(a1z,'a1z','spin in z dir for black hole 1',iunit)
 call write_inopt(a2x,'a2x','spin in x dir for black hole 2',iunit)
 call write_inopt(a2y,'a2y','spin in y dir for black hole 2',iunit)
 call write_inopt(a2z,'a2z','spin in z dir for black hole 2',iunit)
 call write_inopt(mass1,'mass1','mass of black hole 1',iunit)
 call write_inopt(mass2,'mass2','mass of black hole 2',iunit)
 call write_inopt(trajectory_file,'trajectory_file','file containing binary black hole trajectory',iunit)

end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(db,nerr)
 use infile_utils, only:inopts,read_inopt
 use io,           only:error
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(a1x,'a1x',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(a1y,'a1y',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(a1z,'a1z',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(a2x,'a2x',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(a2y,'a2y',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(a2z,'a2z',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(mass1,'mass1',db,errcount=nerr,min=0.,max=1.e12)
 call read_inopt(mass2,'mass2',db,errcount=nerr,min=0.,max=1.e12)
 call read_inopt(trajectory_file,'trajectory_file',db,errcount=nerr)

 a = sqrt(a1x**2 + a1y**2 + a1z**2 + tiny_param)
 a_bh2 = sqrt(a2x**2 + a2y**2 + a2z**2 + tiny_param)
 if (a > 1.) then
    call error('metric','black hole spin: a > 1 for black hole 1')
    nerr = nerr + 1
 endif
 if (a_bh2 > 1.) then
    call error('metric','black hole spin: a > 1 for black hole 2')
    nerr = nerr + 1
 endif

end subroutine read_options_metric

end module metric
