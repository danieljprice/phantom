!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! radially extends an accretion disc, matching surface density, tilt
!  and twist (but assumes a certain sigma profile)
!
! :References: None
!
! :Owner: Rebecca Nealon
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
 use eos,        only:gamma,polyk,get_spsound,ieos
 use io,         only:id,master,fatal
 use setdisc,    only:set_disc,get_disc_mass,scaled_sigma
 use physcon,    only:pi
 use vectorutils,only:rotatevec
 use prompting,      only:prompt
 use centreofmass,   only:reset_centreofmass,get_total_angular_momentum
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=40)      :: filename,infile
 real, dimension(:,:), allocatable :: xyzh_add,vxyzu_add(:,:)
 integer :: ii,ipart,n_killed,i,sigmaprofile,n_add,ii_match,ierr,jj
 integer, parameter :: iunit = 23
 real    :: xp(3),va(3)
 real    :: hfact,pmass,mass_in,rhoi,HonR_r_cyl,spsound,rad_vel,vel_x,vel_y,theta
 real    :: R_in,R_out,R_ref,q_value,p_value,HonR
 real    :: cyl_radius,cyl_width,cyl_in,cyl_out
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
    call close_db(db)
 else
    call fatal('moddump_extenddisc','Cannot find the *.discparams file, try again.')
 endif
 cyl_radius = 10.0
 cyl_width = .5
 mass_in = 0.001
 rad_vel = 5.
 call prompt('Radius of the ring',cyl_radius)
 call prompt('Width of the ring',cyl_width)
 call prompt('Mass of the ring',mass_in)
 call prompt('Radial velocity (spsound units)',rad_vel)

 cyl_in = cyl_radius-cyl_width/2
 cyl_out = cyl_radius+cyl_width/2
 hfact = 1.2

 pmass = massoftype(igas)
 n_add = int(mass_in/pmass)
 allocate(xyzh_add(4,n_add),vxyzu_add(4,n_add))
 ! Construct the new disc
 call set_disc(id,master        = master,             &
               mixture          = .false.,            &
               npart            = n_add,              &
               rref             = R_ref,              &
               rmin             = cyl_in,             &
               rmax             = cyl_out,            &
               p_index          = 0.0,                &
               q_index          = q_value,            &
               HoverR           = HonR,               &
               disc_mass        = mass_in,            &
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



 xp = (/cyl_radius, 0., 0./)
 HonR_r_cyl = HonR*(cyl_radius/R_ref)**(1/2+q_value)
 print*,HonR_r_cyl
 rhoi = 1/sqrt(2*pi)*mass_in/(2*pi*(cyl_out - cyl_in))*1/HonR_r_cyl
 rad_vel = get_spsound(ieos,xp,rhoi,vxyzu_add(:,1))*rad_vel
 ipart = npart
 do ii = 1,n_add
       ipart = ipart + 1
       theta = atan(xyzh_add(2, ii), xyzh_add(1, ii))
       vel_x = rad_vel*cos(theta)
       vel_y = rad_vel*sin(theta)
       va = (/vel_x, vel_y, 0.0/)
       vxyzu_add(1:3,ii) = vxyzu_add(1:3,ii) + va
       call  add_or_update_particle(igas, xyzh_add(1:3,ii), vxyzu_add(1:3,ii), xyzh_add(4,ii), &
            vxyzu_add(4,ii), ipart, npart, npartoftype, xyzh, vxyzu)
 enddo


 print "(a,f5.1,/)",' Disc successfully added to ',cyl_radius

 deallocate(xyzh_add,vxyzu_add)

 return
end subroutine modify_dump

end module moddump
