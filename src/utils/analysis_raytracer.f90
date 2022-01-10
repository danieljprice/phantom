!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which computes neighbour lists for all particles
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: getneighbours
!
   use raytracer,        only:get_all_tau_inwards, get_all_tau_outwards
   use part,             only:rhoh
   use dump_utils,       only:read_array_from_file
   use getneighbours,    only:generate_neighbour_lists, read_neighbours, write_neighbours, &
                              neighcount,neighb,neighmax
   implicit none
   character(len=20), parameter, public :: analysistype = 'raytracer'
   public :: do_analysis

   private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
   character(len=*), intent(in) :: dumpfile
   integer,          intent(in) :: num,npart,iunit
   real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
   real,             intent(in) :: particlemass,time
   
   logical :: existneigh
   character(100) :: neighbourfile
   real :: rho(npart), kappa(npart), tau(npart)
   integer :: i,ierr,iu1,iu2
   real :: start, finish

   !***************************************
   ! Assign neighbour lists to particles by searching shared list of host cell
   !***************************************

   call read_array_from_file(123,dumpfile,'kappa',kappa(:),ierr)
   if (ierr/=0) then
      write(*,*)
      write(*,'("WARNING: could not read kappa from file. It will be set to zero")')
      write(*,*)
      kappa = 0.
   endif

   ! Construct neighbour lists for derivative calculations
   ! write points
   open(newunit=iu1, file='points.txt', status='replace', action='write')
   do i=1, size(xyzh(1,:))
      write(iu1, *) xyzh(1:3,i)
   enddo

   ! get neighbours
   neighbourfile = 'neigh_'//TRIM(dumpfile)
   inquire(file=neighbourfile,exist = existneigh)

   if (existneigh.eqv..true.) then
      print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
      call read_neighbours(neighbourfile,npart)
   else
      ! If there is no neighbour file, generate the list
      print*, 'No neighbour file found: generating'
      call generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile,.true.)
   endif

   do i=1,npart
      rho(i) = rhoh(xyzh(4,i), particlemass)
   enddo
   print*,maxval(kappa)
   print*,''
   print*, 'Start calculating optical depth'
   call cpu_time(start)
   call get_all_tau_outwards(npart-1, npart, xyzh(1:3,:), neighb, rho*kappa, 5, tau)
   call cpu_time(finish)
   print*,'Time = ',finish-start,' seconds.'
   print*, 'Calculating optical depth complete'

   open(newunit=iu2, file='taus.txt', status='replace', action='write')
   do i=1, size(tau)
      write(iu2, *) tau(i)
   enddo
   print*,maxval(tau)
 
end subroutine do_analysis
end module analysis