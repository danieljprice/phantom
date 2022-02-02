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
! :Owner: Mats Esseldeurs
!
! :Runtime parameters: None
!
! :Dependencies: getneighbours
!
   use raytracer,        only:get_all_tau_inwards, get_all_tau_outwards, get_all_tau_optimised
   use part,             only:rhoh,isdead_or_accreted
   use dump_utils,       only:read_array_from_file
   use getneighbours,    only:generate_neighbour_lists, read_neighbours, write_neighbours, &
                              neighcount,neighb,neighmax
   use dust_formation,   only:kappa_dust_bowen
   use linklist, only:set_linklist,allocate_linklist,deallocate_linklist
   implicit none
   character(len=20), parameter, public :: analysistype = 'raytracer'
   public :: do_analysis

   private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
   Use omp_lib
   character(len=*), intent(in) :: dumpfile
   integer,          intent(in) :: num,npart,iunit
   real(kind=8),     intent(in) :: xyzh(:,:),vxyzu(:,:)
   real(kind=8),     intent(in) :: particlemass,time
   
   logical :: existneigh
   character(100) :: neighbourfile
   character(100)   :: jstring
   real(kind=8)   :: primsec(4,2), rho(npart), kappa(npart), temp(npart), u(npart), xyzh2(4,npart), vxyzu2(4,npart)
   real(kind=8), dimension(:), allocatable :: tau
   integer :: i,j,ierr,iu1,iu2,iu3,iu4, npart2
   integer :: start, finish
   real :: totalTime, timeTau
   logical :: tess = .false.

   print*,'("Reading kappa from file")'
   call read_array_from_file(123,dumpfile,'kappa',kappa(:),ierr, 1)
   if (ierr/=0) then
      print*,''
      print*,'("WARNING: could not read kappa from file. It will be set to zero")'
      print*,''
      kappa = 0.
   endif

   if (kappa(1) <= 0. .and. kappa(2) <= 0. .and. kappa(2) <= 0.) then
      print*,'("Reading temperature from file")'
      call read_array_from_file(123,dumpfile,'temperature',temp(:),ierr, 1)
      if (temp(1) <= 0. .and. temp(2) <= 0. .and. temp(2) <= 0.) then
         print*,'("Reading internal energy from file")'
         call read_array_from_file(123,dumpfile,'u',u(:),ierr, 1)
         do i=1,npart
            temp(i)=(1.2-1)*2.381*u(i)*1.6735337254999998e-24*1.380649e-16
         enddo
      endif
      do i=1,npart
         kappa(i)=kappa_dust_bowen(temp(i))
      enddo
   endif
   
   j=1
   do i=1,npart
      if (.not.isdead_or_accreted(xyzh(4,i))) then
         xyzh2(:,j) = xyzh(:,i)
         vxyzu2(:,j) = vxyzu(:,i)
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
   close(iu3)

   call read_array_from_file(123,dumpfile,'x',primsec(1,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'y',primsec(2,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'z',primsec(3,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'h',primsec(4,:),ierr, 2)
   xyzh2(:,npart2+1) = primsec(:,1)
   xyzh2(:,npart2+2) = primsec(:,2)

   ! write points
   open(newunit=iu1, file='points_'//dumpfile//'.txt', status='replace', action='write')
   do i=1, npart2+2
      write(iu1, *) xyzh2(1:3,i)
   enddo
   close(iu1)

   ! get neighbours
   if (.not.tess) then
      neighbourfile = 'neigh_'//TRIM(dumpfile)
      inquire(file=neighbourfile,exist = existneigh)
      if (existneigh.eqv..true.) then
         print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
         call read_neighbours(neighbourfile,npart2+2)
      else
         ! If there is no neighbour file, generate the list
         print*, 'No neighbour file found: generating'
         call system_clock(start)
         call generate_neighbour_lists(xyzh2,vxyzu2,npart2+2,dumpfile, .false.)
         call system_clock(finish)
         totalTime = (finish-start)/1000.
         print*,'Time = ',totalTime,' seconds.'
         call write_neighbours(neighbourfile, npart2+2)
         print*, 'Neighbour finding complete for file ', TRIM(dumpfile)
      endif
   else
      allocate(neighb(npart2+2,100))
      neighb = 0
      inquire(file='neighborb_tess.txt',exist = existneigh)
      if (existneigh.eqv..true.) then
         print*, 'Neighbour file neighbors.txt found'
      else
         call execute_command_line('python3 getNeigh.py -f '//'points_'//dumpfile//'.txt')
      endif
      open(newunit=iu4, file='neighbors_tess.txt', status='old', action='read')
      do i=1, npart2+2
            read(iu4,*) neighb(i,:)
      enddo
      close(iu4)
   endif
   
   call set_linklist(npart2,npart2,xyzh2,vxyzu)

   if (.false..and..not.tess) then
      print*,''
      print*, 'Start calculating optical depth inwards'
      call system_clock(start)
      call get_all_tau_inwards(npart2+1, xyzh2(1:3,:), xyzh2, neighb, rho*kappa*1.496e+13, &
                                 real(2.37686663,8), tau, npart2+2,real(0.1,8))
      call system_clock(finish)
      timeTau = (finish-start)/1000.
      print*,'Time = ',timeTau,' seconds.'
      open(newunit=iu4, file='times_'//dumpfile//'.txt', status='replace', action='write')
      write(iu4, *) timeTau
      close(iu4)
      totalTime = totalTime + timeTau
      open(newunit=iu2, file='taus_'//dumpfile//'_inwards.txt', status='replace', action='write')
      do i=1, size(tau)
         write(iu2, *) tau(i)
      enddo
      close(iu2)

      do j = 0, 9
         write(jstring,'(i0)') j
         print*,''
         print*, 'Start calculating optical depth outwards: ', trim(jstring)
         call system_clock(start)
         call get_all_tau_outwards(npart2+1, xyzh2(1:3,:), xyzh2, neighb, rho*kappa*1.496e+13, &
                                    real(2.37686663,8), j, tau, npart2+2,real(0.1,8))
         call system_clock(finish)
         timeTau = (finish-start)/1000.
         print*,'Time = ',timeTau,' seconds.'
         open(newunit=iu4, file='times_'//dumpfile//'.txt',position='append', status='old', action='write')
         write(iu4, *) timeTau
         close(iu4)
         totalTime = totalTime + timeTau
         open(newunit=iu2, file='taus_'//dumpfile//'_'//trim(jstring)//'.txt', status='replace', action='write')
         do i=1, size(tau)
            write(iu2, *) tau(i)
         enddo
         close(iu2)
      enddo
      print*,''
      print*,'Total time of the calculation = ',totalTime,' seconds.'
   else if (tess) then
      print*,''
      print*, 'Start calculating optical depth inwards'
      call system_clock(start)
      call get_all_tau_inwards(npart2+1, xyzh2(1:3,:), xyzh2, neighb, rho*kappa*1.496e+13, &
                                 real(2.37686663,8), tau, npart2+2,real(0.1,8))
      call system_clock(finish)
      timeTau = (finish-start)/1000.
      print*,'Time = ',timeTau,' seconds.'
      open(newunit=iu4, file='times_'//dumpfile//'_tess.txt', status='replace', action='write')
      write(iu4, *) timeTau
      close(iu4)
      totalTime = totalTime + timeTau
      open(newunit=iu2, file='taus_'//dumpfile//'_tess_inwards.txt', status='replace', action='write')
      do i=1, size(tau)
         write(iu2, *) tau(i)
      enddo
      close(iu2)

      do j = 0, 9
         write(jstring,'(i0)') j
         print*,''
         print*, 'Start calculating optical depth outwards: ', trim(jstring)
         call system_clock(start)
         call get_all_tau_outwards(npart2+1, xyzh2(1:3,:), xyzh2, neighb, rho*kappa*1.496e+13, &
                                    real(2.37686663,8), j, tau, npart2+2,real(0.1,8))
         call system_clock(finish)
         timeTau = (finish-start)/1000.
         print*,'Time = ',timeTau,' seconds.'
         open(newunit=iu4, file='times_'//dumpfile//'_tess.txt',position='append', status='old', action='write')
         write(iu4, *) timeTau
         close(iu4)
         totalTime = totalTime + timeTau
         open(newunit=iu2, file='taus_'//dumpfile//'_tess_'//trim(jstring)//'.txt', status='replace', action='write')
         do i=1, size(tau)
            write(iu2, *) tau(i)
         enddo
         close(iu2)
      enddo
      print*,''
      print*,'Total time of the calculation = ',totalTime,' seconds.'
   endif

   if (.false.) then
      print*,''
      print*, 'Start calculating optical depth optimised'
      call system_clock(start)
      call get_all_tau_optimised(npart2+1, xyzh2(1:3,:), xyzh2, neighb, rho*kappa*1.496e+13, &
                                 real(2.37686663,8),0 , 7, tau, npart2+2,real(0.1,8))
      call system_clock(finish)
      timeTau = (finish-start)/1000.
      print*,'Time = ',timeTau,' seconds.'
      totalTime = totalTime + timeTau
      open(newunit=iu2, file='taus_'//dumpfile//'_optimised.txt', status='replace', action='write')
      do i=1, size(tau)
         write(iu2, *) tau(i)
      enddo
      close(iu2)
   endif

   if (.false.) then
      open(newunit=iu4, file='times_'//dumpfile//'_scaling.txt', status='replace', action='write')
      close(iu4)
      do i=0, omp_get_num_procs()-1!int(log(real(omp_get_num_procs()))/log(2.))
         call omp_set_num_threads(omp_get_num_procs()-i)
         print*,omp_get_max_threads()
         call deallocate_linklist
         call allocate_linklist
         call set_linklist(npart2,npart2,xyzh2,vxyzu)
         call system_clock(start)
         call get_all_tau_outwards(npart2+1, xyzh2(1:3,:), xyzh2, neighb, rho*kappa*1.496e+13, &
                                    real(2.37686663,8), 6, tau, npart2+2,real(0.1,8))
         call system_clock(finish)
         timeTau = (finish-start)/1000.
         print*,'nthread = ',omp_get_max_threads(),': Time = ',timeTau,' seconds.'
         open(newunit=iu4, file='times_'//dumpfile//'_scaling.txt',position='append', status='old', action='write')
         write(iu4, *) omp_get_max_threads(), timeTau
         close(iu4)
      enddo
   endif

end subroutine do_analysis
end module analysis