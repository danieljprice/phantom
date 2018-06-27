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
!  radially extends an accretion disc, matching surface density, tilt
!  and twist (but assumes a certain sigma profile)
!
!  REFERENCES: None
!
!  OWNER: Bec Nealon
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, eos, infile_utils, io, part, partinject,
!    physcon, prompting, setdisc, vectorutils
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

 integer,parameter :: nr = 200

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use partinject, only:add_or_update_particle
 use part,       only:igas,isdead_or_accreted,xyzmh_ptmass,vxyz_ptmass,nptmass
 use eos,        only:gamma,polyk
 use io,         only:id,master,fatal
 use setdisc,    only:set_disc,get_disc_mass,scaled_sigma
 use physcon,    only:pi
 use vectorutils,only:rotatevec
 use prompting,      only:prompt
 use centreofmass,   only:reset_centreofmass
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=40)      :: filename,infile
 real, dimension(:,:), allocatable :: xyzh_add,vxyzu_add(:,:)
 integer :: ii,ipart,n_killed,i,sigmaprofile,n_add,ii_match,n_accreted,ierr
 integer, parameter :: iunit = 23
 real    :: hfact,pmass,R_match
 real    :: R_in,R_out,R_ref,q_value,p_value,HonR
 real    :: radius,R_c,add_mass_disc,R_ext
 real    :: sigma_norm,star_m,toomre_min,enc_m(4096),dummy(4096)
 real    :: Ltot(3,nr),Lunit(3,nr),rad(nr),dr,area,sigma_match
 real    :: sigma(nr),tilt(nr),twist(nr),rotate_about_z,rotate_about_y
 real    :: L_match(3),L_mag,term(3),termmag
 logical :: iexist
 type(inopts), allocatable :: db(:)

 ! Prompt user for old disc parameters and new outer radius
 print "(/,2a,/)",'Most of the disc parameters are read from the *.discparams file.'
 call prompt('Enter the extension (e.g. [].discparams):',filename)

 infile = trim(filename)//'.discparams'
 inquire(file=trim(infile),exist=iexist)
 if (iexist) then
    print "(/,2a,/)",' Reading default values from ', trim(infile)
    call open_db_from_file(db,infile,iunit,ierr)
    call read_inopt(R_ref,'R_ref',db,ierr)
    call read_inopt(R_in,'R_in',db,ierr)
    call read_inopt(R_out,'R_out',db,ierr)
    call read_inopt(p_value,'p_index',db,ierr)
    call read_inopt(q_value,'q_index',db,ierr)
    call read_inopt(HonR,'H/R_ref',db,ierr)
    call read_inopt(star_m,'M_star',db,ierr)
    call close_db(db)
 else
    call fatal('moddump_extenddisc','Cannot find the *.discparams file, try again.')
 endif

 call prompt('Outer radius of new disc?',R_ext)

 ! Run a couple of checks
 if (R_out > R_ext) then
    call fatal('moddump_extenddisc','the extended radius is inside the current disc, try again')
 endif

 R_c = 1.0
 R_match = 2./3.*R_out

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 ! Old disc parameters
 pmass = massoftype(igas)
 sigmaprofile = 2
 sigma_norm = 1.0
 hfact = 1.2

 ! Run a basic disc analysis to determine tilt and twist as a function of radius
 sigma = 0.
 tilt = 0.
 twist = 0.
 Ltot = 0.
 n_killed = 0
 n_accreted = 0
 Lunit = 0.
 ii_match = -1

 dr = (R_out - R_in)/real(nr-1)
 do i=1,nr
    rad(i)=R_in + real(i-1)*dr
    if (ii_match < 0 .and. rad(i) > R_match) ii_match = i
 enddo

 R_match = rad(ii_match)

 do i = 1,npart
    ! i for particle number, ii for radial bin
    radius = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
    if (xyzh(4,i) > 0.) then
       ii = int((radius-rad(1))/dr + 1)
       if (ii > nr .or. ii < 1) cycle
       area = (pi*((rad(ii)+dr/2.)**2-(rad(ii)- dr/2.)**2))
       sigma(ii) = sigma(ii) + pmass/area
       Ltot(1,ii) = Ltot(1,ii) + pmass*(xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i))
       Ltot(2,ii) = Ltot(2,ii) + pmass*(xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i))
       Ltot(3,ii) = Ltot(3,ii) + pmass*(xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i))
    else
       n_accreted = n_accreted + 1
    endif
    if (radius > R_match) then
       xyzh(4,i) = -1.0
       n_killed = n_killed + 1
    endif
 enddo

 ! Average a little
 sigma_match = (sigma(ii_match - 1) + sigma(ii_match) + sigma(ii_match + 1))/3.

 sigma_norm = scaled_sigma(R_match,sigmaprofile,p_value,R_ref,R_in,R_c)
 sigma_norm = sigma_match/sigma_norm

 ! Guess a better p value, assuming the new sigma function is not smoothed
 p_value = log(sigma_match/sigma_norm)/log(R_c/R_match)

 call get_disc_mass(add_mass_disc,enc_m,dummy,toomre_min,0,sigma_norm, &
                         star_m,p_value,q_value,R_match,R_ext,R_ref,R_c,HonR)

 print*,n_killed,' particles removed outside of R=',R_match

 ! Number of injected particles is given by existing particle mass and total added disc mass
 n_add = int(add_mass_disc/pmass)
 allocate(xyzh_add(4,n_add),vxyzu_add(4,n_add))

 ! Construct the new disc
 call set_disc(id,master        = master,             &
               mixture          = .false.,            &
               npart            = n_add,              &
               rref             = R_ref,              &
               rmin             = R_match,            &
               rmax             = R_ext,              &
               p_index          = p_value,            &
               q_index          = q_value,            &
               HoverR           = HonR,               &
               disc_mass        = add_mass_disc,      &
               gamma            = gamma,              &
               particle_mass    = pmass,              &
               hfact            = hfact,              &
               xyzh             = xyzh_add,           &
               vxyzu            = vxyzu_add,          &
               polyk            = polyk,              &
               ismooth          = .false.,            &
               inclination      = 0.,                 &
               rwarp            = R_in,               &
               warp_smoothl     = 0.0)

 ! Tilt        and twist the particles
 ! The tilt could be set by set_warp instead but this way gives tilt AND twist
 ! We also have to rotate velocities (for accurate L)

 L_match = (Ltot(:,ii_match - 1) + Ltot(:,ii_match) + Ltot(:,ii_match + 1))/3.
 L_mag = sqrt(dot_product(L_match,L_match))
 term = (/L_match(1),L_match(2),0./)
 termmag = sqrt(dot_product(term,term))

 rotate_about_z = -acos(dot_product((/1.,0.,0./),term/termmag))
 rotate_about_y = acos(dot_product((/0.,0.,1./),L_match/L_mag))

 ! Now add those new particles to existing disc
 ipart = npart ! The initial particle number
 do ii = 1,n_add
    call rotatevec(xyzh_add(1:3,ii),(/0.,0.,1.0/),rotate_about_z)
    call rotatevec(xyzh_add(1:3,ii),(/0.,1.0,0./),rotate_about_y)
    call rotatevec(vxyzu_add(1:3,ii),(/0.,0.,1.0/),rotate_about_z)
    call rotatevec(vxyzu_add(1:3,ii),(/0.,1.0,0./),rotate_about_y)

    ipart = ipart + 1
    call  add_or_update_particle(igas, xyzh_add(1:3,ii), vxyzu_add(1:3,ii), xyzh_add(4,ii), &
             vxyzu_add(4,ii), ipart, npart, npartoftype, xyzh, vxyzu)
 enddo

 print "(a,f5.1,/)",' Disc successfully extended to ',R_ext

 deallocate(xyzh_add,vxyzu_add)

 return
end subroutine modify_dump

end module moddump

