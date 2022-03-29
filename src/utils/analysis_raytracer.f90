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
! :Dependencies: raytracer, part, dump_utils, dust_formation, linklist
!
   use raytracer_copy,   only:get_all_tau_inwards, get_all_tau_adaptive
   use raytracer,        only:get_all_tau_outwards
   use part,             only:rhoh,isdead_or_accreted
   use dump_utils,       only:read_array_from_file
   use getneighbours,    only:generate_neighbour_lists, read_neighbours, write_neighbours, &
                              neighcount,neighb,neighmax
   use dust_formation,   only:calc_kappa_bowen
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
   character(100)   :: jstring, kstring
   real(kind=8)   :: primsec(4,2), rho(npart), kappa(npart), temp(npart), u(npart), xyzh2(4,npart), vxyzu2(4,npart)
   real(kind=8), dimension(:), allocatable :: tau
   integer :: i,j,k,ierr,iu1,iu2,iu3,iu4, npart2
   integer :: start, finish, method, analyses, minOrder, maxOrder
   real :: totalTime, timeTau
   logical :: SPH, calcInwards

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
         kappa(i)=calc_kappa_bowen(temp(i))
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

   do i=1,npart2
      rho(i) = rhoh(xyzh2(4,i), particlemass)
   enddo

   call read_array_from_file(123,dumpfile,'x',primsec(1,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'y',primsec(2,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'z',primsec(3,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'h',primsec(4,:),ierr, 2)

   call set_linklist(npart2,npart2,xyzh2,vxyzu)

   print *,'What do you want to an analyses?'
   print *, '(1) Analysis'
   print *, '(2) Method'
   print *, '(3) Preloaded settings'
   read *,analyses
   ! analyses=3
   
   if (analyses==1) then
      print *,'Which analysis would you like to run?'
      print *, '(1) Analysis Outwards'
      print *, '(2) Analysis Adaptive'
      print *, '(3) Scaling'
      read *,method
      if (method == 1) then
         SPH = .true.
         calcInwards = .false.
         print *,'At which order would you like to start?'
         read *,minOrder
         print *,'At which order would you like to stop?'
         read *,maxOrder
      else if (method == 2) then
         SPH = .true.
         calcInwards = .false.
         print *,'At which order would you like to start?'
         read *,minOrder
         print *,'At which order would you like to stop?'
         read *,maxOrder
      else if (method == 3) then
         
      endif
   else if (analyses==2) then
      print *,'Which algorithm would you like to run?'
      print *, '(1) Inwards'
      print *, '(2) Outwards'
      print *, '(3) Adaptive'
      read *,method
      if (method == 1) then
         print *,'Do you want to use SPH neighbours? (T/F)'
         read*,SPH
      else if (method == 2) then
         print *,'What order do you want to run?'
         read*,j
         write(jstring,'(i0)') j
      else if (method == 3) then
         print *,'What order do you want to run? (integer below 7)'
         read*,j
         write(jstring,'(i0)') j
         print *,'What refinement level do you want to run? (integer below 7)'
         read*,k
         write(kstring,'(i0)') k
      endif
   endif

   if ((analyses==1 .and. calcInwards) .or. (analyses==2 .and. method==1)) then ! get neighbours
      if (SPH) then
         neighbourfile = 'neigh_'//TRIM(dumpfile)
         inquire(file=neighbourfile,exist = existneigh)
         if (existneigh.eqv..true.) then
            print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
            call read_neighbours(neighbourfile,npart2)
         else
            ! If there is no neighbour file, generate the list
            print*, 'No neighbour file found: generating'
            call system_clock(start)
            call generate_neighbour_lists(xyzh2,vxyzu2,npart2,dumpfile, .false.)
            call system_clock(finish)
            totalTime = (finish-start)/1000.
            print*,'Time = ',totalTime,' seconds.'
            call write_neighbours(neighbourfile, npart2)
            print*, 'Neighbour finding complete for file ', TRIM(dumpfile)
         endif
      else
         allocate(neighb(npart2+2,100))
         neighb = 0
         inquire(file='neighbors_tess.txt',exist = existneigh)
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
   endif

   if (analyses==1) then
      if (method == 1) then
         if (calcInwards) then
            print*,''
            print*, 'Start calculating optical depth inwards'
            call system_clock(start)
            call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, rho*kappa*1.496e+13, &
                                       2.37686663, tau, primsec(:,2),0.1)
            call system_clock(finish)
            timeTau = (finish-start)/1000.
            print*,'Time = ',timeTau,' seconds.'
            open(newunit=iu4, file='times_'//dumpfile//'.txt', status='replace', action='write')
               write(iu4, *) timeTau
            close(iu4)
            totalTime = timeTau
            open(newunit=iu2, file='taus_'//dumpfile//'_inwards.txt', status='replace', action='write')
            do i=1, size(tau)
               write(iu2, *) tau(i)
            enddo
            close(iu2)
            deallocate(neighb)
         else
            open(newunit=iu4, file='times_'//dumpfile//'.txt', status='replace', action='write')
               write(iu4, *) 0.
            close(iu4)
            totalTime=0
         endif
   
         do j = minOrder, maxOrder
            write(jstring,'(i0)') j
            print*,''
            print*, 'Start calculating optical depth outwards: ', trim(jstring)
            call system_clock(start)
            call get_all_tau_outwards(primsec(:,1), xyzh2, rho*kappa*1.496e+13, &
                                       2.37686663, j, tau, primsec(:,2),0.1)
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
      else if (method == 2) then
         if (calcInwards) then
            print*,''
            print*, 'Start calculating optical depth inwards'
            call system_clock(start)
            call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, rho*kappa*1.496e+13, &
                                       2.37686663, tau, primsec(:,2),0.1)
            call system_clock(finish)
            timeTau = (finish-start)/1000.
            print*,'Time = ',timeTau,' seconds.'
            open(newunit=iu4, file='times_'//dumpfile//'.txt', status='replace', action='write')
               write(iu4, *) timeTau
            close(iu4)
            totalTime = timeTau
            open(newunit=iu2, file='taus_'//dumpfile//'_inwards.txt', status='replace', action='write')
            do i=1, size(tau)
               write(iu2, *) tau(i)
            enddo
            close(iu2)
            deallocate(neighb)
         else
            open(newunit=iu4, file='times_opt_'//dumpfile//'.txt', status='replace', action='write')
               write(iu4, *) 0.
            close(iu4)
            totalTime=0
         endif
   
         do j = minOrder, maxOrder
            write(jstring,'(i0)') j
            do k = minOrder,maxOrder-j
               write(kstring,'(i0)') k
               print*,''
               print*, 'Start calculating optical depth outwards: minOrder = ', trim(jstring),', refineLevel = ', trim(kstring)
               call system_clock(start)
               call get_all_tau_adaptive(primsec(:,1), xyzh2, rho*kappa*1.496e+13, &
                                          2.37686663, j, k, tau, primsec(:,2),0.1)
               call system_clock(finish)
               timeTau = (finish-start)/1000.
               print*,'Time = ',timeTau,' seconds.'
               open(newunit=iu4, file='times_opt_'//dumpfile//'.txt',position='append', status='old', action='write')
               write(iu4, *) timeTau
               close(iu4)
               totalTime = totalTime + timeTau
               open(newunit=iu2, file='taus_'//dumpfile//'_opt_'//trim(jstring)// &
                     '_'//trim(kstring)//'.txt', status='replace', action='write')
               do i=1, size(tau)
                  write(iu2, *) tau(i)
               enddo
               close(iu2)
            enddo
         enddo
         print*,''
         print*,'Total time of the calculation = ',totalTime,' seconds.'
      else if (method == 3) then
         print*,'Start doing scaling analysis'
         open(newunit=iu4, file='times_'//dumpfile//'_scaling.txt', status='replace', action='write')
         close(iu4)
         do i=1, omp_get_max_threads()
            call omp_set_num_threads(i)
            call deallocate_linklist
            call allocate_linklist
            call set_linklist(npart2,npart2,xyzh2,vxyzu)
            call system_clock(start)
            call get_all_tau_outwards(primsec(:,1), xyzh2, rho*kappa*1.496e+13, &
                                       2.37686663, 5, tau, primsec(:,2),0.1)
            call system_clock(finish)
            timeTau = (finish-start)/1000.
            print*,'nthread = ',omp_get_max_threads(),': Time = ',timeTau,' seconds.'
            open(newunit=iu4, file='times_'//dumpfile//'_scaling.txt',position='append', status='old', action='write')
            write(iu4, *) omp_get_max_threads(), timeTau
            close(iu4)
         enddo
      endif
   else if (analyses==2) then
      if (method == 1) then
         print*,''
         print*, 'Start calculating optical depth inwards'
         call system_clock(start)
         call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, rho*kappa*1.496e+13, &
                                    2.37686663, tau, primsec(:,2),0.1)
         call system_clock(finish)
         timeTau = (finish-start)/1000.
         print*,'Time = ',timeTau,' seconds.'
         if (SPH) then
            open(newunit=iu2, file='taus_'//dumpfile//'_inwards.txt', status='replace', action='write')
         else
            open(newunit=iu2, file='taus_'//dumpfile//'_tess_inwards.txt', status='replace', action='write')
         endif
         do i=1, size(tau)
            write(iu2, *) tau(i)
         enddo
         close(iu2)
      else if (method == 2) then
         print*,''
         print*, 'Start calculating optical depth outwards: ', trim(jstring)
         call system_clock(start)
         call get_all_tau_outwards(primsec(:,1), xyzh2, rho*kappa*1.496e+13, 2.37686663, j, tau, primsec(:,2),0.1)
         call system_clock(finish)
         timeTau = (finish-start)/1000.
         print*,'Time = ',timeTau,' seconds.'
         open(newunit=iu2, file='taus_'//dumpfile//'_'//trim(jstring)//'.txt', status='replace', action='write')
         do i=1, size(tau)
            write(iu2, *) tau(i)
         enddo
         close(iu2)
      else if (method == 3) then
         print*,''
         print*, 'Start calculating optical depth adaptive: minOrder = ', trim(jstring),', refineLevel = ', trim(kstring)
         call system_clock(start)
         call get_all_tau_adaptive(primsec(:,1), xyzh2, rho*kappa*1.496e+13, &
                                    2.37686663, j, k, tau, primsec(:,2),0.1)
         call system_clock(finish)
         timeTau = (finish-start)/1000.
         print*,'Time = ',timeTau,' seconds.'
         totalTime = totalTime + timeTau
         open(newunit=iu2, file='taus_'//dumpfile//'_opt_'//trim(jstring)// &
               '_'//trim(kstring)//'.txt', status='replace', action='write')
         do i=1, size(tau)
            write(iu2, *) tau(i)
         enddo
         close(iu2)
      endif
   endif

   if (analyses == 3) then
      do i=1,npart
         if (norm2(xyzh2(1:3,i) - (/10.,10.,10./)) < 4.) then
            kappa(i) = 1e10
         endif
      enddo
      ! allocate(neighb(npart2+2,100))
      ! neighb = 0
      ! open(newunit=iu4, file='neighbors_tess.txt', status='old', action='read')
      ! do i=1, npart2+2
      !       read(iu4,*) neighb(i,:)
      ! enddo
      ! close(iu4)
      print*,''
      print*, 'Start calculating optical depth outwards'
      call system_clock(start)
      call get_all_tau_outwards(primsec(:,1), xyzh2, rho*kappa*1.496e+13, &
                                 2.37686663, 7, tau, primsec(:,2),0.1)
      call system_clock(finish)
      timeTau = (finish-start)/1000.
      print*,'Time = ',timeTau,' seconds.'
      totalTime = totalTime + timeTau
      open(newunit=iu2, file='taus_'//dumpfile//'_raypolation_7.txt', status='replace', action='write')
      do i=1, size(tau)
         write(iu2, *) tau(i)
      enddo
      close(iu2)
   endif

   if (.false.) then! Write out the location and density of all the points
      open(newunit=iu1, file='points_'//dumpfile//'.txt', status='replace', action='write')
      do i=1, npart2+2
         write(iu1, *) xyzh2(1:3,i)
      enddo
      close(iu1)

      open(newunit=iu3, file='rho_'//dumpfile//'.txt', status='replace', action='write')
      do i=1,npart2
         rho(i) = rhoh(xyzh2(4,i), particlemass)
         write(iu3, *) rho(i)
      enddo
      close(iu3)
   endif

end subroutine do_analysis
end module analysis