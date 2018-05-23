!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  Input is a relaxed star, output is two relaxed stars in binary orbit
!
!  REFERENCES: None
!
!  OWNER: Terrence Tricco
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, dim, externalforces, initial_params, io,
!    options, part, prompting, readwrite_dumps
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas
 use prompting,      only: prompt
 use centreofmass,   only: reset_centreofmass,get_centreofmass
 use initial_params, only: get_conserv
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: opt, synchro, Nstar1, Nstar2, nx
 real :: sep,mtot,angvel,vel1,vel2
 real :: xcom(3), vcom(3), x1com(3), v1com(3), x2com(3), v2com(3)
 real :: pmassi,m1,m2

 !
 ! Option selection
 !
 print *, 'Running moddump_binarystar:'
 print *, ''
 print *, 'This utility sets two stars in binary orbit around each other, or modifies an existing binary.'
 print *, ''
 print *, 'Options:'
 print *, '   1) Duplicate a relaxed star'
 print *, '   2) Add a star from another dumpfile'
 print *, '   3) Adjust separation of existing binary'
! print *, '   4) Add low-density background fluid'

 opt = 1
 call prompt('Choice',opt, 1, 3)

 sep = 1.0
 print *, ''
 call prompt('Enter radial separation between stars (code unit)', sep, 0.)

 synchro = 1
 call prompt('Synchronise binaries? [0 false; 1 true]',synchro, 0, 1)
 print *, ''

 !
 ! Create binary
 !
 if (opt == 1) then
    call duplicate_star(npart, npartoftype, xyzh, vxyzu, Nstar1, Nstar2)
 endif

 if (opt == 2) then
    call add_star(npart, npartoftype, xyzh, vxyzu, Nstar1, Nstar2)
 endif

 if (opt == 3) then
    call determine_Nstar(npart,xyzh,Nstar1,Nstar2)
 endif

 if (opt == 4) then
    nx = 128
    call prompt('Specify number of particles in x-direction for low-density background',nx,0)
    print *, ''
    call determine_Nstar(npart,xyzh,Nstar1,Nstar2)
    call get_centreofmass(xcom,  vcom,  npart,  xyzh,                   vxyzu)
    call get_centreofmass(x1com, v1com, Nstar1, xyzh(:,1:Nstar1),       vxyzu(:,1:Nstar1))
    call get_centreofmass(x2com, v2com, Nstar2, xyzh(:,Nstar1+1:npart), vxyzu(:,Nstar1+1:npart))
!    call add_background(npart,npartoftype,massoftype,xyzh,vxyzu,Nstar1,Nstar2,x1com,x2com,nx)
!     call add_magneticfield(npart,xyzh,Nstar1,Nstar2,x1com,x2com)
 endif


 !
 ! Set binary at specified separation with corresponding orbital velocity
 !
 call get_centreofmass(xcom,  vcom,  npart,  xyzh,                   vxyzu)
 call get_centreofmass(x1com, v1com, Nstar1, xyzh(:,1:Nstar1),       vxyzu(:,1:Nstar1))
 call get_centreofmass(x2com, v2com, Nstar2, xyzh(:,Nstar1+1:npart), vxyzu(:,Nstar1+1:npart))

 pmassi = massoftype(igas)
 mtot   = npart  * pmassi
 m1     = Nstar1 * pmassi
 m2     = Nstar2 * pmassi

 print *, '   Mass of first star:  ', m1
 print *, '   Mass of second star: ', m2
 print *, ''

 call adjust_sep(npart,xyzh,vxyzu,Nstar1,Nstar2,sep,x1com,v1com,x2com,v2com)

 angvel = sqrt(1.0 * mtot / sep**3) ! angular velocity
 vel1   = m1 * sep / mtot * angvel
 vel2   = m2 * sep / mtot * angvel

 call get_centreofmass(x1com, v1com, Nstar1, xyzh(:,1:Nstar1),       vxyzu(:,1:Nstar1))
 call get_centreofmass(x2com, v2com, Nstar2, xyzh(:,Nstar1+1:npart), vxyzu(:,Nstar1+1:npart))

 call reset_velocity(npart,vxyzu)

 call set_velocity(npart,vxyzu,Nstar1,Nstar2,angvel,vel1,vel2)
! call set_corotate_velocity(angvel)


 !
 ! synchronise rotation
 !
 if (synchro == 1) then
    call synchronise(npart,xyzh,vxyzu,Nstar1,Nstar2,angvel,x1com,x2com)
 endif


 !
 ! reset centre of mass to origin
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 get_conserv = 1.

end subroutine modify_dump


!
! Take the star from the input file and duplicate it some distance apart.
! This assumes the dump file only has one star.
!
subroutine duplicate_star(npart,npartoftype,xyzh,vxyzu,Nstar1,Nstar2)
 use part,         only: igas,set_particle_type,temperature
 use dim,          only: store_temperature
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(out)   :: Nstar1, Nstar2
 integer :: i
 real :: sep

 print *, 'Duplicating star'
 print *, ''

 npart = npartoftype(igas)

 sep = 10.0

 ! duplicate relaxed star
 do i = npart+1, 2*npart
    ! place star a distance rad away
    xyzh(1,i) = xyzh(1,i-npart) + sep
    xyzh(2,i) = xyzh(2,i-npart)
    xyzh(3,i) = xyzh(3,i-npart)
    xyzh(4,i) = xyzh(4,i-npart)
    vxyzu(1,i) = vxyzu(1,i-npart)
    vxyzu(2,i) = vxyzu(2,i-npart)
    vxyzu(3,i) = vxyzu(3,i-npart)
    vxyzu(4,i) = vxyzu(4,i-npart)
    if (store_temperature) then
       temperature(i) = temperature(i-npart)
    endif
    call set_particle_type(i,igas)
 enddo

 Nstar1 = npart
 Nstar2 = npart

 npart = 2 * npart
 npartoftype(igas) = npart

end subroutine duplicate_star


!
! Place a star that is read from another dumpfile
!
subroutine add_star(npart,npartoftype,xyzh,vxyzu,Nstar1,Nstar2)
 use part,            only: igas,set_particle_type,temperature,alphaind
 use prompting,       only: prompt
 use dim,             only: maxp,maxvxyzu,nalpha,maxalpha,store_temperature
 use readwrite_dumps, only: read_dump
 use io,              only: idisk1,iprint
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(out)   :: Nstar1, Nstar2
 character(len=120) :: fn
 real, allocatable :: xyzh2(:,:)
 real, allocatable :: vxyzu2(:,:)
 real, allocatable :: temperature2(:)
 real, allocatable :: alphaind2(:,:)
 integer :: i,ierr
 real    :: time2,hfact2,sep

 print *, 'Adding a new star read from another dumpfile'
 print *, ''

 fn = ''
 call prompt('Name of second dumpfile',fn)

 ! read_dump will overwrite the current particles, so store them in a temporary array
 allocate(xyzh2(4,maxp),stat=ierr)  ! positions + smoothing length
 if (ierr /= 0) stop ' error allocating memory to store positions'
 allocate(vxyzu2(maxvxyzu,maxp),stat=ierr)  ! velocity + thermal energy
 if (ierr /= 0) stop ' error allocating memory to store velocity'
 if (store_temperature) then        ! temperature
    allocate(temperature2(maxp),stat=ierr)
    if (ierr /= 0) stop ' error allocating memory to store temperature'
 endif
 if (maxalpha == maxp) then         ! artificial viscosity alpha
    allocate(alphaind2(nalpha,maxp),stat=ierr)
    if (ierr /= 0) stop ' error allocating memory to store alphaind'
 endif

 Nstar2 = npart
 xyzh2  = xyzh
 vxyzu2 = vxyzu
 if (store_temperature) then
    temperature2 = temperature
 endif
 if (maxalpha == maxp) then
    alphaind2 = alphaind
 endif


 ! read second dump file
 call read_dump(trim(fn),time2,hfact2,idisk1+1,iprint,0,1,ierr)
 if (ierr /= 0) stop 'error reading second dumpfile'


 Nstar1 = npart
 sep = 10.0

 ! insert saved star (from original dump file)
 do i = npart+1, npart+Nstar2
    ! place star a distance rad away
    xyzh(1,i) = xyzh2(1,i-npart) + sep
    xyzh(2,i) = xyzh2(2,i-npart)
    xyzh(3,i) = xyzh2(3,i-npart)
    xyzh(4,i) = xyzh2(4,i-npart)
    vxyzu(1,i) = vxyzu2(1,i-npart)
    vxyzu(2,i) = vxyzu2(2,i-npart)
    vxyzu(3,i) = vxyzu2(3,i-npart)
    vxyzu(4,i) = vxyzu2(4,i-npart)
    if (store_temperature) then
       temperature(i) = temperature2(i-npart)
    endif
    if (maxalpha == maxp) then
       alphaind(1,i) = real(alphaind2(1,i-npart),kind=4)
       alphaind(2,i) = real(alphaind2(2,i-npart),kind=4)
    endif
    call set_particle_type(i,igas)
 enddo

 npart = npart + Nstar2
 npartoftype(igas) = npart

end subroutine add_star


!
! If editing an existing binary, determine number of particles in each star based on their position
!
subroutine determine_Nstar(npart,xyzh,Nstar1,Nstar2)
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 integer, intent(out)   :: Nstar1, Nstar2
 integer :: Nstar1xneg, Nstar1xpos, Nstar1yneg, Nstar1ypos
 integer :: Nstar2xneg, Nstar2xpos, Nstar2yneg, Nstar2ypos
 integer :: i
 logical :: done

 print *, 'Auto-Determining which particles belong in each star'
 print *, '  (warning: assumes stars are separated, in x-y plane, and com is at x=y=0)'
 print *, '        (and does not work if ambient background already present)'
 print *, ''

 !
 ! We do this in a brute force simple way.
 !
 ! There exists a cut along x=0 or y=0 that separates the two stars.
 ! So try both and take the one that works.
 !
 ! We also need to know which star is ordered first in memory.
 ! Means we try 4 cuts in total (left/right and right/left for x=0, etc)
 !
 ! Strategy: Try looping from i=1 onwards for 4 cuts to find Nstar1
 !    then loop from i=npart downwards to find Nstar2
 !
 ! Solution exists when Nstar1 + Nstar2 = npart
 ! (If they don't, then cut is going through the stars or is backwards.)
 !

 Nstar1xneg = 0
 Nstar1xpos = 0
 Nstar1yneg = 0
 Nstar1ypos = 0

 i = 1
 done = .false.
 do while (.not.done)
    if (xyzh(1,i) < 0.0) then
       Nstar1xneg = NStar1xneg + 1
    else
       done = .true.
    endif
    i = i + 1
 enddo
 i = 1
 done = .false.
 do while (.not.done)
    if (xyzh(1,i) > 0.0) then
       Nstar1xpos = NStar1xpos + 1
    else
       done = .true.
    endif
    i = i + 1
 enddo
 i = 1
 done = .false.
 do while (.not.done)
    if (xyzh(2,i) < 0.0) then
       Nstar1yneg = NStar1yneg + 1
    else
       done = .true.
    endif
    i = i + 1
 enddo
 i = 1
 done = .false.
 do while (.not.done)
    if (xyzh(2,i) > 0.0) then
       Nstar1ypos = NStar1ypos + 1
    else
       done = .true.
    endif
    i = i + 1
 enddo


 Nstar2xneg = 0
 Nstar2xpos = 0
 Nstar2yneg = 0
 Nstar2ypos = 0

 i = npart
 done = .false.
 do while (.not.done)
    if (xyzh(1,i) < 0.0) then
       Nstar2xneg = NStar2xneg + 1
    else
       done = .true.
    endif
    i = i - 1
 enddo
 i = npart
 done = .false.
 do while (.not.done)
    if (xyzh(1,i) > 0.0) then
       Nstar2xpos = NStar2xpos + 1
    else
       done = .true.
    endif
    i = i - 1
 enddo
 i = npart
 done = .false.
 do while (.not.done)
    if (xyzh(2,i) < 0.0) then
       Nstar2yneg = NStar2yneg + 1
    else
       done = .true.
    endif
    i = i - 1
 enddo
 i = npart
 done = .false.
 do while (.not.done)
    if (xyzh(2,i) > 0.0) then
       Nstar2ypos = NStar2ypos + 1
    else
       done = .true.
    endif
    i = i - 1
 enddo

 if (Nstar1xneg + Nstar2xpos == npart) then
    Nstar1 = Nstar1xneg
    Nstar2 = Nstar2xpos
 endif
 if (Nstar1xpos + Nstar2xneg == npart) then
    Nstar1 = Nstar1xpos
    Nstar2 = Nstar2xneg
 endif
 if (Nstar1yneg + Nstar2ypos == npart) then
    Nstar1 = Nstar1yneg
    Nstar2 = Nstar2ypos
 endif
 if (Nstar1ypos + Nstar2yneg == npart) then
    Nstar1 = Nstar1ypos
    Nstar2 = Nstar2yneg
 endif
 print *, '   Nstar1 = ', Nstar1
 print *, '   Nstar2 = ', Nstar2
 print *, ''

end subroutine determine_Nstar


!
! Add an ambient background fluid
!
!subroutine add_background(npart,npartoftype,massoftype,xyzh,vxyzu,Nstar1,Nstar2,x1com,x2com,nx)
! use part,     only: hfact,igas,set_particle_type
! use unifdis,  only: set_unifdis
! use io,       only: master
! use boundary, only: set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,totvol
! use units,    only: unit_density
! integer, intent(inout) :: npart
! integer, intent(inout) :: npartoftype(:)
! real,    intent(in)    :: massoftype(:)
! real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
! integer, intent(in)    :: Nstar1, Nstar2
! real,    intent(in)    :: x1com(:),x2com(:)
! integer, intent(in)    :: nx
! integer :: id
! integer :: i
! real :: dx, dy, dz, dr
! real :: rad1, rad2
! real :: xlen, deltax, bgdens
! integer(kind=8) :: npart_total
!
! print *, 'Adding uniform low-density background fluid'
! print *, ''
!
! rad1 = 0.0
! do i = 1, Nstar1
!    dx   = xyzh(1,i) - x1com(1)
!    dy   = xyzh(2,i) - x1com(2)
!    dz   = xyzh(3,i) - x1com(3)
!    dr   = sqrt(dx*dx + dy*dy + dz*dz)
!    rad1 = max(dr, rad1)
! enddo
!
! rad2 = 0.0
! do i = Nstar1+1, npart
!    dx   = xyzh(1,i) - x2com(1)
!    dy   = xyzh(2,i) - x2com(2)
!    dz   = xyzh(3,i) - x2com(3)
!    dr   = sqrt(dx*dx + dy*dy + dz*dz)
!    rad2 = max(dr, rad2)
! enddo
!
! print *, "radius of star 1: ", rad1
! print *, "radius of star 2: ", rad2
!
! xlen = 0.437 * 0.5  ! 10^10 cm
! call set_boundary(-xlen,xlen,-xlen,xlen,-xlen,xlen)
!
! deltax = (2.0 * xlen) / nx
!
! if (deltax > rad1 .or. deltax > rad2) then
!    print *, 'WARNING: particle separation greater than star radius'
!    print *, '        suggest using nx=', int(2.0 * xlen / rad1)
!    print *, ''
! endif
!
! id = 1
! npart_total = int(npart,kind=4)
! call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total) 
!
! npart = npart_total
! do i = Nstar1 + Nstar2 + 1, npart
!    call set_particle_type(i,igas)
! enddo
! 
! npartoftype(igas) = npart
!
! bgdens = (npart-Nstar1-Nstar2) * massoftype(igas) / totvol 
! print *, '   Added background fluid of density: ', (bgdens*unit_density), ' g/cm^3'
! print *, ''
!
!end subroutine add_background


!
! Adjust the separation of the two stars.
! First star is placed at the origin, second star is placed sep away in x
!
subroutine adjust_sep(npart,xyzh,vxyzu,Nstar1,Nstar2,sep,x1com,v1com,x2com,v2com)
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in)    :: Nstar1, Nstar2
 real,    intent(in)    :: x1com(:),v1com(:),x2com(:),v2com(:)
 real,    intent(in)    :: sep
 integer :: i

 print *, 'Placing stars apart at a distance: ', sep
 print *, ''

 do i = 1, Nstar1
    xyzh(1,i) = xyzh(1,i) - x1com(1)
    xyzh(2,i) = xyzh(2,i) - x1com(2)
    xyzh(3,i) = xyzh(3,i) - x1com(3)
    vxyzu(1,i) = vxyzu(1,i) - v1com(1)
    vxyzu(2,i) = vxyzu(2,i) - v1com(2)
    vxyzu(3,i) = vxyzu(3,i) - v1com(3)
 enddo

 do i = Nstar1+1, npart
    xyzh(1,i) = xyzh(1,i) - x2com(1) + sep
    xyzh(2,i) = xyzh(2,i) - x2com(2)
    xyzh(3,i) = xyzh(3,i) - x2com(3)
    vxyzu(1,i) = vxyzu(1,i) - v2com(1)
    vxyzu(2,i) = vxyzu(2,i) - v2com(2)
    vxyzu(3,i) = vxyzu(3,i) - v2com(3)
 enddo

end subroutine adjust_sep


!
! reset velocities
!
subroutine reset_velocity(npart,vxyzu)
 integer, intent(in)    :: npart
 real,    intent(inout) :: vxyzu(:,:)

 vxyzu(1:3,:) = 0.0

end subroutine reset_velocity


!
! Set corotation external force on using angular velocity
!
subroutine set_corotate_velocity(angvel)
 use options,        only:iexternalforce
 use externalforces, only: omega_corotate,iext_corotate
 real,    intent(in)    :: angvel

 print "(/,a,es18.10,/)", ' The angular velocity in the corotating frame is: ', angvel

 ! Turns on corotation
 iexternalforce = iext_corotate
 omega_corotate = angvel

end subroutine


!
! Set orbital velocity in normal space
!
subroutine set_velocity(npart,vxyzu,Nstar1,Nstar2,angvel,vel1,vel2)
 integer, intent(in)    :: npart
 real,    intent(inout) :: vxyzu(:,:)
 integer, intent(in)    :: Nstar1, Nstar2
 real,    intent(in)    :: angvel
 real,    intent(in)    :: vel1,vel2
 integer :: i

 print *, "Setting stars in mutual orbit with angular velocity: ", angvel
 print *, "  Adding bulk velocity |v| = ", vel1, " to first star"
 print *, "                       |v| = ", vel2, " to second star"
 print *, ""

 ! Adjust bulk velocity of relaxed star towards second star
 do i = 1, Nstar1
    vxyzu(2,i) = vxyzu(2,i) + vel1
 enddo

 do i = Nstar1+1, npart
    vxyzu(2,i) = vxyzu(2,i) - vel2
 enddo

end subroutine set_velocity


!
! Set binaries in synchronised orbit
!
subroutine synchronise(npart,xyzh,vxyzu,Nstar1,Nstar2,angvel,x1com,x2com)
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 integer, intent(in)    :: Nstar1, Nstar2
 real,    intent(in)    :: angvel
 real,    intent(in)    :: x1com(:), x2com(:)
 integer :: i

 print *, "Synchronising rotation to orbital period"
 print *, ""

 do i = 1, Nstar1
    vxyzu(1,i) = vxyzu(1,i) + angvel * (xyzh(2,i) - x1com(2))
    vxyzu(2,i) = vxyzu(2,i) - angvel * (xyzh(1,i) - x1com(1))
 enddo

 do i = Nstar1+1, npart
    vxyzu(1,i) = vxyzu(1,i) + angvel * (xyzh(2,i) - x2com(2))
    vxyzu(2,i) = vxyzu(2,i) - angvel * (xyzh(1,i) - x2com(1))
 enddo

end subroutine synchronise

end module moddump

