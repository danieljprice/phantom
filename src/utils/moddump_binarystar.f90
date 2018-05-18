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
!    options, part, prompting, readwrite_dumps, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas
 use units,          only: set_units
 use prompting,      only: prompt
 use centreofmass,   only: reset_centreofmass,get_centreofmass
 use initial_params, only: get_conserv
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: opt, Nstar1, Nstar2
 real :: sep,mtot,angvel,vel1,vel2
 real :: xcom(3), vcom(3), x1com(3), v1com(3), x2com(3), v2com(3)
 real :: pmassi,m1,m2

 print *, 'Running moddump_binarystar:'
 print *, ''
 print *, 'This utility sets two stars in binary orbit around each other, or modifies an existing binary.'
 print *, ''
 print *, 'Options:'
 print *, '   1) Duplicate a relaxed star'
 print *, '   2) Add a star from another dumpfile'
 print *, '   3) Adjust separation of existing binary'

 opt = 1
 call prompt('Choice',opt, 1, 3)

 if (opt  /=  1 .and. opt  /=  2 .and. opt /= 3) then
    print *, 'Incorrect option selected. Doing nothing.'
    return
 endif

 sep = 10.0
 print *, ''
 call prompt('Enter radial separation between stars (in code unit)', sep, 0.)
 print *, ''

 ! duplicate star if chosen
 if (opt == 1) then
    call duplicate_star(npart, npartoftype, massoftype, xyzh, vxyzu, Nstar1, Nstar2)
 endif

 ! add a new star from another dumpfile
 if (opt == 2) then
    call add_star(npart, npartoftype, massoftype, xyzh, vxyzu, Nstar1, Nstar2)
 endif

 ! add a uniform low density background fluid
! if (opt == 3) then
!    call add_background(npart, npartoftype, massoftype, xyzh, vxyzu)
! endif

 ! get centre of mass of total system, and each star individually
 call get_centreofmass(xcom,  vcom,  npart,  xyzh,                   vxyzu)
 call get_centreofmass(x1com, v1com, Nstar1, xyzh(:,1:Nstar1),       vxyzu(:,1:Nstar1))
 call get_centreofmass(x2com, v2com, Nstar2, xyzh(:,Nstar1+1:npart), vxyzu(:,Nstar1+1:npart))

 ! we only work with gas, so this is ok
 pmassi = massoftype(igas)
 mtot   = npart  * pmassi
 m1     = Nstar1 * pmassi
 m2     = Nstar2 * pmassi

 ! adjust separation of binary
 call adjust_sep(npart,npartoftype,massoftype,xyzh,vxyzu,Nstar1,Nstar2,sep,x1com,v1com,x2com,v2com)


! mtot = npart*massoftype(igas)
 angvel = sqrt(1.0 * mtot / sep**3)   ! angular velocity
 vel1   = m1 * sep / mtot * angvel
 vel2   = m2 * sep / mtot * angvel

 ! find the centre of mass position and velocity for each star
 call get_centreofmass(x1com, v1com, Nstar1, xyzh(:,1:Nstar1),       vxyzu(:,1:Nstar1))
 call get_centreofmass(x2com, v2com, Nstar2, xyzh(:,Nstar1+1:npart), vxyzu(:,Nstar1+1:npart))

 ! set orbital velocity
 call set_velocity(npart,npartoftype,massoftype,xyzh,vxyzu,Nstar1,Nstar2,angvel,vel1,vel2)
! call set_corotate_velocity(angvel)




 ! reset centre of mass of the binary system
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 get_conserv = 1.

end subroutine modify_dump


!
! Take the star from the input file and duplicate it some distance apart.
! This assumes the dump file only has one star.
!
subroutine duplicate_star(npart,npartoftype,massoftype,xyzh,vxyzu,Nstar1,Nstar2)
 use part,         only: igas,set_particle_type,igas,temperature
 use prompting,    only: prompt
 use centreofmass, only: reset_centreofmass
 use dim,          only: store_temperature
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(out)   :: Nstar1, Nstar2
 integer :: i
 real :: sep

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
subroutine add_star(npart,npartoftype,massoftype,xyzh,vxyzu,Nstar1,Nstar2)
 use part,            only: igas,set_particle_type,igas,temperature,alphaind
 use prompting,       only: prompt
 use centreofmass,    only: reset_centreofmass
 use dim,             only: maxp,maxvxyzu,nalpha,maxalpha,store_temperature
 use readwrite_dumps, only: read_dump
 use io,              only: idisk1,iprint
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(out)   :: Nstar1, Nstar2
 character(len=120) :: fn
 real, allocatable :: xyzh2(:,:)
 real, allocatable :: vxyzu2(:,:)
 real, allocatable :: temperature2(:)
 real, allocatable :: alphaind2(:,:)
 integer :: i,ierr
 real    :: time2,hfact2,sep


 print *, ''
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

 print *, npart
end subroutine add_star


!
! Adjust the separation of the two stars.
! First star is placed at the origin, second star is placed sep away in x
!
subroutine adjust_sep(npart,npartoftype,massoftype,xyzh,vxyzu,Nstar1,Nstar2,sep,x1com,v1com,x2com,v2com)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in)    :: Nstar1, Nstar2
 real,    intent(in)    :: x1com(:),v1com(:),x2com(:),v2com(:)
 real,    intent(in)    :: sep
 integer :: i

 print *, sep
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
subroutine set_velocity(npart,npartoftype,massoftype,xyzh,vxyzu,Nstar1,Nstar2,angvel,vel1,vel2)
 use units,        only: unit_velocity
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in)    :: Nstar1, Nstar2
 real,    intent(in)    :: angvel
 real,    intent(in)    :: vel1,vel2
 integer :: i

 print *, "Setting stars in mutual orbit with angular velocity ", angvel
 print *, "  Adding bulk velocity |v| = ", vel1, "( = ", (vel1*unit_velocity), &
                  " physical units) to first star"
 print *, "                       |v| = ", vel2, "( = ", (vel2*unit_velocity), &
                  " physical units) to second star"
 print *, ""

 ! Adjust bulk velocity of relaxed star towards second star
 vxyzu(1,:) = 0.
 vxyzu(2,:) = 0.
 vxyzu(3,:) = 0.
 do i = 1, Nstar1
    vxyzu(2,i) = vxyzu(2,i) + vel1
 enddo

 do i = Nstar1+1, npart
    vxyzu(2,i) = vxyzu(2,i) - vel2
 enddo

end subroutine set_velocity


end module moddump

