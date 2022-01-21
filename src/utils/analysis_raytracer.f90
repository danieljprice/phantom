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
   use part,             only:rhoh,isdead_or_accreted
   use dump_utils,       only:read_array_from_file
   use getneighbours,    only:generate_neighbour_lists, read_neighbours, write_neighbours, &
                              neighcount,neighb,neighmax
   use dust_formation,   only:kappa_dust_bowen
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
   real           :: primsec(4,2), rho(npart), kappa(npart), temp(npart), xyzh2(4,npart)
   real(kind=8), dimension(:), allocatable :: tau
   integer :: i,j,ierr,iu1,iu2,iu3, npart2
   integer            :: nums(8,4)
   integer :: start, finish
   real :: timeNeigh, timeTau

   call read_array_from_file(123,dumpfile,'kappa',kappa(:),ierr, 1)
   print*,(ierr)
   if (ierr/=0) then
      print*,''
      print*,'("WARNING: could not read kappa from file. It will be set to zero")'
      print*,''
      kappa = 0.
   endif

   if (kappa(1) <= 0. .and. kappa(2) <= 0. .and. kappa(2) <= 0.) then
      call read_array_from_file(123,dumpfile,'temperature',temp(:),ierr, 1)
      do i=1,npart
         kappa(i)=kappa_dust_bowen(temp(i))
      enddo
   endif
   
   j=1
   do i=1,npart
      if (.not.isdead_or_accreted(xyzh(4,i))) then
         xyzh2(:,j) = xyzh(:,i)
         kappa(j) = kappa(i)
         j=j+1
      endif
   enddo
   npart2 = j-1
   print*,'npart = ',npart2
   allocate(tau(npart2))

   open(newunit=iu3, file='rho_'//dumpfile//'.txt', status='replace', action='write')
   do i=1,npart2
      rho(i) = rhoh(xyzh2(4,i), particlemass)
      write(iu3, *) rho(i)
   enddo

   call read_array_from_file(123,dumpfile,'x',primsec(1,:),ierr, 2)
   print*,(ierr)
   call read_array_from_file(123,dumpfile,'y',primsec(2,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'z',primsec(3,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'h',primsec(4,:),ierr, 2)
   xyzh2(:,npart2+1) = primsec(:,1)
   xyzh2(:,npart2+2) = primsec(:,2)

   ! Construct neighbour lists for derivative calculations
   ! write points
   open(newunit=iu1, file='points_'//dumpfile//'.txt', status='replace', action='write')
   do i=1, npart2+2
      write(iu1, *) xyzh2(1:3,i)
   enddo

   ! get neighbours
   neighbourfile = 'neigh_'//TRIM(dumpfile)
   inquire(file=neighbourfile,exist = existneigh)

   if (existneigh.eqv..true.) then
      print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
      call read_neighbours(neighbourfile,npart2+2)
   else
      ! If there is no neighbour file, generate the list
      print*, 'No neighbour file found: generating'
      call system_clock(start)
      call generate_neighbour_lists(xyzh2,vxyzu,npart2+2,dumpfile,.false.)
      call system_clock(finish)
      timeNeigh = (finish-start)/1000.
      print*,'Time = ',timeNeigh,' seconds.'
      call write_neighbours(neighbourfile, npart2+2)
      print*, 'Neighbour finding complete for file ', TRIM(dumpfile)
   endif

   print*,''
   print*, 'Start calculating optical depth'
   call system_clock(start)
   call get_all_tau_inwards(npart2+1, xyzh2(1:3,:), neighb, rho*kappa, 2.37686663, tau, npart2+2,0.1)
   call system_clock(finish)
   timeTau = (finish-start)/1000.
   print*,'Time = ',timeTau,' seconds.'

   open(newunit=iu2, file='taus_'//dumpfile//'_inwards.txt', status='replace', action='write')
   do i=1, size(tau)
      write(iu2, *) tau(i)
   enddo
   print*,''
   print*,'Total time of the calculation = ',timeNeigh+timeTau,' seconds.'
 
end subroutine do_analysis
end module analysis
