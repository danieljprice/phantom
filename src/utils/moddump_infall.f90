!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Adds either a sphere or cylinder of material to infall
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, eos, infile_utils, io, part, partinject,
!   physcon, prompting, setdisc, vectorutils
!
 implicit none

 integer,parameter :: nr = 200

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use partinject, only:add_or_update_particle
 use part,       only:igas,isdead_or_accreted,xyzmh_ptmass,vxyz_ptmass,nptmass
 use eos,        only:gamma,polyk
 use io,         only:id,master,fatal
 use spherical,  only:set_sphere
 use physcon,    only:pi
 use vectorutils,only:rotatevec
 use prompting,      only:prompt
 use centreofmass,   only:reset_centreofmass,get_total_angular_momentum
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, dimension(:,:), allocatable :: xyzh_add,vxyzu_add(:,:)
 integer :: in_shape,in_orbit,ii,ipart,n_killed,i,sigmaprofile,n_add,ii_match,ierr,jj,np
 integer(kind=8) :: nptot
 integer, parameter :: iunit = 23
 real    :: r_close,in_mass,hfact,pmass,ecc,delta,r_init,r_in
 real    :: vp(3), xp(3)
 real    :: dma,n0,pf,m0,x0,y0,z0,r0,vx0,vy0,vz0,mtot
! real    :: R_in,R_out,R_ref,q_value,p_value,HonR
! real    :: radius,R_c,add_mass_disc,R_ext
! real    :: sigma_norm,star_m,toomre_min,enc_m(4096),dummy(4096)
! real    :: Ltot(3,nr),Lunit(3,nr),rad(nr),dr,area,sigma_match
! real    :: sigma(nr),tilt(nr),twist(nr),rotate_about_z,rotate_about_y
! real    :: L_match(3),L_mag,term(3),termmag
 logical :: iexist
 type(inopts), allocatable :: db(:)

 r_close = 100.
 in_mass = 0.01
 r_in = 100.0
 r_init = 1000.0
 in_orbit = 1
 in_shape = 0

 ! Prompt user for infall material shape
 call prompt('Enter the infall material shape (0=sphere, 1=cylinder)',in_shape,0,1)
 call prompt('Enter radius of shape:', r_in, 0.1)

 ! Prompt user for the infall material orbit
 call prompt('Enter orbit type (0=bound, 1=parabolic, 2=hyperbolic)', in_orbit,0) 
 if (in_orbit == 1) then
   ecc = 1.0
   call prompt('Enter closest approach in au:', r_close, 0.) 
   call prompt('Enter infall mass in Msun:', in_mass, 0.01) 
   call prompt('Enter initial radial distance in au:', r_init, 0.0) 
 endif


 dma = r_close
 n0  = r_init/r_close

 mtot=sum(xyzmh_ptmass(4,:))

 !--focal parameter dma = pf/2
 pf = 2*dma

 !--define m0 = -x0/dma such that r0 = n0*dma
 !  companion starts at negative x and y
 !  positive root of 1/8*m**4 + m**2 + 2(1-n0**2) = 0
 !  for n0 > 1
 m0 = 2*sqrt(n0-1.0)

 !--perturber initial position
 x0 = -m0*dma
 y0 = dma*(1.0-(x0/pf)**2)
 z0 = 0.0
 xp = (/x0,y0,z0/)
 print*, 'Initial cloud centre is: ', xp

 !--perturber initial velocity
 r0  = sqrt(x0**2+y0**2+z0**2)
 vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
 vy0 = -(x0/r0)*sqrt(mtot/pf)
 vz0 = 0.0
 vp  = (/vx0,vy0,vz0/)


 ! Number of injected particles is given by existing particle mass and total
 ! added disc mass
 pmass = massoftype(igas)
 n_add = int(in_mass/pmass)
 print* , 'Number of particles to add ', n_add
 allocate(xyzh_add(4,n_add),vxyzu_add(4,n_add))
 hfact = 1.2
 delta = 1.0 ! no idea what this is
 nptot= n_add + npartoftype(igas) 
 np = 0
 call set_sphere('random',id,master,0.,r_in,delta,hfact,np,xyzh_add,xyz_origin=(/x0, y0, z0/),&
      np_requested=n_add, nptot=nptot) 
 print*, 'We have set the sphere'
 ! Initiate initial velocity of the partiles in the shape
 vxyzu_add(1, :) = vx0
 vxyzu_add(2, :) = vy0
 vxyzu_add(3, :) = vz0
 vxyzu_add(4, :) = vxyzu(4, 1)

!subroutine set_sphere(lattice,id,master,rmin,rmax,delta,hfact,np,xyzh, &
!                      rhofunc,rhotab,rtab,xyz_origin,nptot,dir,exactN,np_requested,mask)
! use stretchmap, only:set_density_profile,rho_func
! character(len=*), intent(in)    :: lattice
! integer,          intent(in)    :: id,master
! integer,          intent(inout) :: np
! real,             intent(in)    :: rmin,rmax,hfact
! real,             intent(out)   :: xyzh(:,:)
! real,             intent(inout) :: delta
! procedure(rho_func), pointer, optional :: rhofunc
! real,             intent(in),    optional :: rhotab(:), rtab(:)
! integer,          intent(in),    optional :: dir
! integer,          intent(in),    optional :: np_requested
! integer(kind=8),  intent(inout), optional :: nptot
! real,             intent(in),    optional :: xyz_origin(3)
! logical,          intent(in),    optional :: exactN


 ! Run a basic disc analysis to determine tilt and twist as a function of radius
 !sigma = 0.
 !tilt = 0.
 !twist = 0.
 !Ltot = 0.
 !n_killed = 0
 !Lunit = 0.
 !ii_match = -1

 !dr = (R_out - R_in)/real(nr-1)
 !do i=1,nr
 !   rad(i)=R_in + real(i-1)*dr
 !   if (ii_match < 0 .and. rad(i) > R_match) ii_match = i
 !enddo

 !if (ii_match < 0) then
 !   call fatal('moddump_extenddisc','cannot match the discs')
 !endif

 !R_match = rad(ii_match)

 !do i = 1,npart
 !   ! i for particle number, ii for radial bin
 !   radius = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
 !   ii = int((radius-rad(1))/dr + 1)
 !   if (xyzh(4,i) > 0. .and. ii <=nr .and. ii > 0) then
 !      area = (pi*((rad(ii)+dr/2.)**2-(rad(ii)- dr/2.)**2))
 !      sigma(ii) = sigma(ii) + pmass/area
 !      Ltot(1,ii) = Ltot(1,ii) + pmass*(xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i))
 !      Ltot(2,ii) = Ltot(2,ii) + pmass*(xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i))
 !      Ltot(3,ii) = Ltot(3,ii) + pmass*(xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i))
 !   endif
 !   if (radius > R_match) then
 !      xyzh(4,i) = -1.0
!       call kill_particle(i)
  !     n_killed = n_killed + 1
 !   endif
 !enddo

! call shuffle_part(npart)

 ! Average a little
 !sigma_match = (sigma(ii_match - 1) + sigma(ii_match) + sigma(ii_match + 1))/3.
 !sigma_norm = scaled_sigma(R_match,sigmaprofile,p_value,R_ref,R_in,R_out,R_c)
 !sigma_norm = sigma_match/sigma_norm


 ! Tilt and twist the particles
 ! The tilt could be set by set_warp instead but this way gives tilt AND twist
 ! We also have to rotate velocities (for accurate L)
 ! NB: these rotations are (correctly) in the opposite order than
 ! when these rotations are used elsewhere in the code

 !L_match = (Ltot(:,ii_match - 1) + Ltot(:,ii_match) + Ltot(:,ii_match + 1))/3.
 !L_mag = sqrt(dot_product(L_match,L_match))
 !term = (/L_match(1),L_match(2),0./)
 !termmag = sqrt(dot_product(term,term))

 !rotate_about_z = acos(dot_product((/1.,0.,0./),term/termmag))
 !rotate_about_y = acos(dot_product((/0.,0.,1./),L_match/L_mag))

 ! Now add those new particles to existing disc
 ipart = npart ! The initial particle number (post shuffle)
 do ii = 1,n_add
    ! Rotate particle to correct position and velocity
    !call rotatevec(xyzh_add(1:3,ii),(/0.,1.0,0./),rotate_about_y)
    !call rotatevec(xyzh_add(1:3,ii),(/0.,0.,1.0/),rotate_about_z)
    !call rotatevec(vxyzu_add(1:3,ii),(/0.,1.0,0./),rotate_about_y)
    !call rotatevec(vxyzu_add(1:3,ii),(/0.,0.,1.0/),rotate_about_z)
    !radius = sqrt(xyzh_add(1,ii)**2 + xyzh_add(2,ii)**2 + xyzh_add(3,ii)**2)
    !jj = int((radius-rad(1))/dr + 1)
    !if (jj > ii_match-1) then
       ! Add the particle
       ipart = ipart + 1
       call  add_or_update_particle(igas, xyzh_add(1:3,ii), vxyzu_add(1:3,ii), xyzh_add(4,ii), &
            vxyzu_add(4,ii), ipart, npart, npartoftype, xyzh, vxyzu)
    !endif
 enddo

 print "(a,f5.1,/)",' Added sphere successfully '

 deallocate(xyzh_add,vxyzu_add)

 return
end subroutine modify_dump

end module moddump
