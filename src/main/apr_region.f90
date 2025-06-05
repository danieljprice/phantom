!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module apr_region
!
! Everything for setting the adaptive particle refinement regions
!
! Current options include:
! 
! 1: A static sphere; absolute position set in *.in file
!
! 2: A sphere around a sink particle; sink is identified as track_part in *.in
!
! 3: A sphere around fragments; designed for use in an SG disc
!
! 4: A sphere centred on two sinks; this is centred between two sequential sinks
!    (i.e. modelling the region around a binary star)
! 
! 5: A sphere centered on the CoM of the system; this is intended for use with sinks
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part
!
 implicit none

 public :: identify_clumps,read_options_apr,write_options_apr

 ! default values for runtime parameters are stored here
 integer :: apr_max_in = 3, ref_dir = 1, apr_type = 1, apr_max = 4
 integer :: top_level = 1, ntrack = 0, ntrack_max = 10,read_track_part
 integer, allocatable :: npart_regions(:), track_part(:)
 real :: apr_rad = 1.0, apr_drad = 0.1, apr_centre_in(3)
 real, allocatable :: apr_regions(:), apr_centre(:,:)
 real, save :: apr_H(2,100)  ! we enforce this to be 100


 logical :: dynamic_apr = .false., apr_region_is_circle = .false.
contains
   
!-----------------------------------------------------------------------
!+
!  Setting/updating the centre of the apr region (as it may move)
!+
!-----------------------------------------------------------------------
subroutine set_apr_centre(apr_type,apr_centre,ntrack,track_part)
 use part, only: xyzmh_ptmass,xyzh,npart,vxyzu,nptmass,vxyz_ptmass
 use centreofmass, only:get_centreofmass
 integer, intent(in)  :: apr_type
 real,    intent(out) :: apr_centre(3,ntrack_max)
 integer, optional, intent(in) :: ntrack,track_part(:)
 integer :: ii
 real :: xcom(3), vcom(3)

 select case (apr_type)

 case(1) ! a static circle
    ! do nothing here

 case(2) ! around sink particle named track_part
    dynamic_apr = .true.
    apr_centre(1,1) = xyzmh_ptmass(1,track_part(1))
    apr_centre(2,1) = xyzmh_ptmass(2,track_part(1))
    apr_centre(3,1) = xyzmh_ptmass(3,track_part(1))

 case(3) ! to derefine a clump - only activated when the centre of the clump
    ! has been found
    dynamic_apr = .true.
    do ii = 1,ntrack
      !print*,ii,apr_centre(:,ii)
     ! apr_centre(1,ii) = xyzh(1,track_part(ii))
     ! apr_centre(2,ii) = xyzh(2,track_part(ii))
     ! apr_centre(3,ii) = xyzh(3,track_part(ii))
      !if (ii > 10) cycle
   enddo
   !if (ntrack > 0) read*
   
 case(4) ! averaging two sequential sinks
    dynamic_apr = .true.
    apr_centre(1,1) = 0.5*(xyzmh_ptmass(1,track_part(1)) + xyzmh_ptmass(1,track_part(1) + 1))
    apr_centre(2,1) = 0.5*(xyzmh_ptmass(2,track_part(1)) + xyzmh_ptmass(2,track_part(1) + 1))
    apr_centre(3,1) = 0.5*(xyzmh_ptmass(3,track_part(1)) + xyzmh_ptmass(3,track_part(1) + 1))

 case(5) ! centre of mass
    dynamic_apr = .true.
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
    apr_centre(:,1) = xcom

 case(6) ! vertically uniformly resolved disc
    dynamic_apr = .true.
    !call run_apr_disc_analysis(100,xyzh,vxyzu,apr_H)

 case default ! used for the test suite
    apr_centre(:,1) = 0.

 end select

end subroutine set_apr_centre

!-----------------------------------------------------------------------
!+
!  Initialising all the apr region arrays that decide
!  the spatial arrangement of the regions
!+
!-----------------------------------------------------------------------
subroutine set_apr_regions(ref_dir,apr_max,apr_regions,apr_rad,apr_drad)
 integer, intent(in) :: ref_dir,apr_max
 real, intent(in)    :: apr_rad,apr_drad
 real, intent(inout) :: apr_regions(apr_max)
 integer :: ii,kk

 if (ref_dir == 1) then
    apr_regions(1) = huge(apr_regions(1)) ! this needs to be a number that encompasses the whole domain
    do ii = 2,apr_max
       kk = apr_max - ii + 2
       apr_regions(kk) = apr_rad + (ii-2)*apr_drad
    enddo
 else
    apr_regions(apr_max) = huge(apr_regions(apr_max)) ! again this just needs to encompass the whole domain
    do ii = 1,apr_max-1
       apr_regions(ii) = apr_rad + (ii-1)*apr_drad
    enddo
 endif

end subroutine set_apr_regions

!-----------------------------------------------------------------------
!+
!  Analysis to see where the s-g clumps are - this function identifies
!  *potential* track_part, these are refined and created in 
!  create_or_update_clumps
!+
!-----------------------------------------------------------------------
subroutine identify_clumps(npart,xyzh,vxyzu,poten,apr_level,xyzmh_ptmass,aprmassoftype, &
                           ntrack_temp,track_part_temp)
 use part, only:igas,rhoh
 use ptmass, only:rho_crit_cgs
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: apr_level(:)
 real, intent(in) :: xyzh(:,:), vxyzu(:,:), aprmassoftype(:,:),xyzmh_ptmass(:,:)
 real(kind=4), intent(in) :: poten(:)
 integer, intent(out) :: ntrack_temp, track_part_temp(ntrack_max)
 integer :: nbins, ii, ibin, nmins, jj, kk
 integer, allocatable :: counter(:), minima(:), min_particle(:)
 real, allocatable :: radius(:), ave_poten(:)
 real :: rin, rout, dbin, dx, dy, dz, rad, gradleft, gradright
 real :: minpoten, pmassi, rhoi

 ! set up arrays
 nbins = 50
 allocate(counter(nbins),radius(nbins),ave_poten(nbins),&
    minima(nbins),min_particle(nbins))

 ! Currently hardwired but this is problematic
 call find_inner_and_outer_radius(npart,xyzh,rin,rout)
 rin = 2.0
 rout = 30.
 dbin = (rout-rin)/real(nbins-1)
 do ii = 1,nbins
    radius(ii) = rin + real(ii-1)*dbin
 enddo

 ave_poten = 0.
 counter = 0
 ! Create an azimuthally averaged potential energy vs. radius profile
 do ii = 1,npart
    dx = xyzh(1,ii) - xyzmh_ptmass(1,1)
    dy = xyzh(2,ii) - xyzmh_ptmass(2,1)
    dz = xyzh(3,ii) - xyzmh_ptmass(3,1)
    rad = sqrt(dx**2 + dy**2 + dz**2)
    pmassi = aprmassoftype(igas,apr_level(ii))

    ibin = int((rad - radius(1))/dbin + 1)
    if ((ibin > nbins) .or. (ibin < 1)) cycle

    ave_poten(ibin) = ave_poten(ibin) + poten(ii)/pmassi
    counter(ibin) = counter(ibin) + 1
 enddo

 ! average with the number of particles in the bin
 do ii = 1,nbins
    if (counter(ii) > 0) then
       ave_poten(ii) = ave_poten(ii)/counter(ii)
    else
       ave_poten(ii) = 0.
    endif
 enddo

 ! Identify what radius the local minima are at
 minima = 0
 nmins = 0
 do ii = 2, nbins-1
    gradleft = (ave_poten(ii) - ave_poten(ii-1))/(radius(ii) - radius(ii-1))
    gradright = (ave_poten(ii+1) - ave_poten(ii))/(radius(ii+1) - radius(ii))
    if (gradleft * gradright < 0.) then
       nmins = nmins + 1
       minima(nmins) = ii
    endif
 enddo
 if (nmins == 0) return

 ! Identify the particles in these minima that have the lowest potential energy
 ! this is quite inefficient, in future should save these above into the bins so
 ! you just need to cycle through the subset? Don't know if this is faster
 do jj = 1,nmins
    minpoten = 1.0
    do ii = 1,npart
       dx = xyzh(1,ii) - xyzmh_ptmass(1,1)
       dy = xyzh(2,ii) - xyzmh_ptmass(2,1)
       dz = xyzh(3,ii) - xyzmh_ptmass(3,1)
       rad = sqrt(dx**2 + dy**2 + dz**2)
       pmassi = aprmassoftype(igas,apr_level(ii))

       ibin = int((rad - radius(1))/dbin + 1)
       if ((ibin == (minima(jj))) .or. &
        (ibin - 1 == (minima(jj))) .or. &
        (ibin + 1 == (minima(jj)))) then ! if it is in the minima bin or adjacents
          if ((poten(ii)/pmassi) < minpoten) then
             minpoten = poten(ii)/pmassi ! the one in this bin with the lowest minimum potential
             min_particle(jj) = ii
          endif
       endif
    enddo
 enddo

 ! Check they are not already within a region of low potential energy
 ! If they are, replace the existing particle as the one to be tracked
 ntrack_temp = 0
 track_part_temp(:) = 0
 over_mins: do jj = 1,nmins
    ii = min_particle(jj)
    ! check that the particle at the lowest potential energy has also met the
    ! density criteria
    pmassi = aprmassoftype(igas,apr_level(ii))
    rhoi = rhoh(xyzh(4,ii),pmassi)
    if (rhoi < rho_crit_cgs) cycle over_mins

    ! check we haven't seen this particle already and that we're not already tracking it
    do kk = 1,ntrack_temp
      if (track_part_temp(kk) == ii) cycle over_mins
    enddo
    do kk = 1,ntrack
      if (track_part(kk) == ii) cycle over_mins
    enddo

    ! otherwise it is a new particle or particle in existing region to track
    ntrack_temp = ntrack_temp + 1
    track_part_temp(ntrack_temp) = ii
 enddo over_mins

 ! tidy up
 deallocate(counter,ave_poten,radius,minima,min_particle)

end subroutine identify_clumps


!-----------------------------------------------------------------------
!+
!  Find inner and outer radius - assuming that the disc is around the
!  central star for now
!+
!-----------------------------------------------------------------------
subroutine find_inner_and_outer_radius(npart,xyzh,rmin,rmax)
 use part, only: xyzmh_ptmass
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:)
 real, intent(out)   :: rmin,rmax
 integer :: ii
 real    :: rmin_test, rmax_test, xi, yi, zi, r2_test

 rmin_test = huge(rmin_test)
 rmax_test = huge(rmax_test) ! just big and small initial guesses

 do ii = 1,npart
   xi = xyzh(1,ii) - xyzmh_ptmass(1,1)
   yi = xyzh(2,ii) - xyzmh_ptmass(2,1)
   zi = xyzh(3,ii) - xyzmh_ptmass(3,1)
   r2_test = xi**2 + yi**2 + zi**2

   if (r2_test < rmin_test) rmin_test = r2_test
   if (r2_test > rmax_test) rmax_test = r2_test

 enddo

 rmin = sqrt(rmin_test)
 rmax = sqrt(rmax_test)

end subroutine find_inner_and_outer_radius


!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_apr(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_apr'
 logical :: igotall1,igotall2

 imatch  = .true.
 igotall1 = .true.
 igotall2 = .true.
 select case(trim(name))
 case('apr_max')
    read(valstring,*,iostat=ierr) apr_max_in
    ngot = ngot + 1
    if (apr_max_in  <  0) call fatal(label,'apr_max < 0 in input options')
 case('ref_dir')
    read(valstring,*,iostat=ierr) ref_dir
    ngot = ngot + 1
 case('apr_type')
    read(valstring,*,iostat=ierr) apr_type
    ngot = ngot + 1
 case('apr_rad')
    read(valstring,*,iostat=ierr) apr_rad
    ngot = ngot + 1
    if (apr_rad  <  tiny(apr_rad)) call fatal(label,'apr_rad too small in input options')
 case('apr_drad')
    read(valstring,*,iostat=ierr) apr_drad
    ngot = ngot + 1
    if (apr_drad  <  tiny(apr_drad)) call fatal(label,'apr_drad too small in input options')
 case default
    imatch = .false.
    select case(apr_type)
    case(1)
       call read_options_apr1(name,valstring,imatch,igotall1,ierr)
    case(2,4)
       call read_options_apr2(name,valstring,imatch,igotall2,ierr)
    end select
 end select
 igotall = (ngot >= 5) .and. igotall1 .and. igotall2

end subroutine read_options_apr

!-----------------------------------------------------------------------
!+
!  extra subroutines for reading in different styles of apr zones
!+
!-----------------------------------------------------------------------
subroutine read_options_apr1(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_apr1'

 imatch  = .true.
 select case(trim(name))
 case('apr_centre(1)')
    read(valstring,*,iostat=ierr) apr_centre_in(1)
    ngot = ngot + 1
 case('apr_centre(2)')
    read(valstring,*,iostat=ierr) apr_centre_in(2)
    ngot = ngot + 1
 case('apr_centre(3)')
    read(valstring,*,iostat=ierr) apr_centre_in(3)
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 3)

end subroutine read_options_apr1

subroutine read_options_apr2(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_apr2'

 imatch  = .true.
 select case(trim(name))
 case('track_part')
    read(valstring,*,iostat=ierr) read_track_part
    ngot = ngot + 1
    if (read_track_part  <  1) call fatal(label,'track_part not chosen in input options')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)

end subroutine read_options_apr2

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_apr(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for adaptive particle refinement'
 call write_inopt(apr_max_in,'apr_max','number of additional refinement levels (3 -> 2x resolution)',iunit)
 call write_inopt(ref_dir,'ref_dir','increase (1) or decrease (-1) resolution',iunit)
 call write_inopt(apr_type,'apr_type','1: static, 2: sink, 3: clumps, 4: sequential sinks, 5: com, 6: vertical',iunit)

 select case (apr_type)
 case (1)
    call write_inopt(apr_centre_in(1),'apr_centre(1)','centre of region x position',iunit)
    call write_inopt(apr_centre_in(2),'apr_centre(2)','centre of region y position',iunit)
    call write_inopt(apr_centre_in(3),'apr_centre(3)','centre of region z position',iunit)
 case(2,4)
    call write_inopt(read_track_part,'track_part','number of sink to track',iunit)
 case default
   ! write nothing
 end select

 call write_inopt(apr_rad,'apr_rad','radius of innermost region',iunit)
 call write_inopt(apr_drad,'apr_drad','size of step to next region',iunit)

end subroutine write_options_apr

end module apr_region