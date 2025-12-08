!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program testbin
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: testbin [no arguments]
!
! :Dependencies: datafiles, prompting, testbinary
!
 use testbinary, only:test_binary
 use prompting,  only:prompt
 use datafiles,  only:find_phantom_datafile
 use fileutils,  only:get_ncolumns,read_column_labels,find_column,load_data_file
 use setorbit,   only:orbit_t,set_orbit_elements
 use orbits,     only:orbit_is_parabolic
 implicit none
 integer :: j,ierr,ierr1,itex,iunit
 integer :: nargs,ncolumns,nheader,nlabels
 integer :: icol1,icol2,icol3,im1,im2,ichisq,iu_output
 real :: a,e,inc,o,w,f,m1,m2,chisq,rp
 character(len=120) :: filename
 character(len=32), allocatable :: labels(:) ! column labels
 real, allocatable :: array(:,:)
 type(orbit_t) :: orbit

 ! get filename from command line
 nargs = command_argument_count()
 if (nargs > 0) then
    call get_command_argument(1,filename)
 else
    filename = find_phantom_datafile('orbits.txt','orbits')
 endif

 itex = 0
 open(newunit=iunit,file=filename,status='old',iostat=ierr)
 if (ierr /= 0) then
    print*,' ERROR: could not open file ',trim(filename)
    stop
 endif
 call get_ncolumns(iunit,ncolumns,nheader)

 if (ncolumns == 8) then
    ! read orbits from orbits.txt file or similar
    if (ierr==0) then
       open(newunit=itex,file='orbits.tex',status='replace',iostat=ierr1)
       print "(a)",' writing to orbits.tex'
    endif
    j = 0
    do while(ierr==0)
       j = j + 1
       read(iunit,*,iostat=ierr) a,e,inc,o,w,f,m1,m2
       if (ierr==0) call test_binary(m1,m2,a,e,inc,o,w,f,j,itex)
    enddo
 else
    ! process orbits from optimisation runs, convert output from
    ! the Orbit Reconstructor^TM to orbital elements
    allocate(labels(ncolumns))
    icol2 = 0; icol3 = 0; im1 = 0; im2 = 0; ichisq = 0
    call read_column_labels(iunit,nheader,ncolumns,nlabels,labels)
    icol1 = find_column(labels,'binary_dvx',verbose=.true.)
    if (icol1 > 0) then
       orbit%input_type = 2
       orbit%obs%dx(:) = [' 346.','-242.','   0.']  
       orbit%obs%dv(:) = [' 0.040','-0.067',' 0.000']
       orbit%flyby%d = '1200.0'
       m1 = 1.0; m2 = 0.8
       chisq = 0.
       icol2 = find_column(labels,'binary_dvy',verbose=.true.)
       icol3 = find_column(labels,'binary_dz',verbose=.true.)
       im1 = find_column(labels,'m1',verbose=.true.)
       im2 = find_column(labels,'m2',verbose=.true.)
       ichisq = find_column(labels,'chisq',verbose=.true.)
    endif
    close(iunit)
    call load_data_file(filename,array,nheader)
    print*,' loaded ',size(array,1),' rows of data'
    print*,' first row of data: ',array(1,:)
    if (allocated(array)) then
       open(newunit=iu_output,file='orbits.csv',status='new',action='write')
       write(iu_output,"(a)") 'iteration,rp,e,i,O,w,f,m1,m2,chisq'
       do j=1,size(array,1)
          print*,' processing row ',j,' dvx=',array(j,icol1),' dvy=',array(j,icol2)
          if (im1 > 0) m1 = array(j,im1)
          if (im2 > 0) m2 = array(j,im2)
          if (icol1 > 0) write(orbit%obs%dv(1),"(g12.4)") array(j,icol1)
          if (icol2 > 0) write(orbit%obs%dv(2),"(g12.4)") array(j,icol2)
          if (icol3 > 0) write(orbit%obs%dx(3),"(g12.4)") array(j,icol3)
          if (ichisq > 0) chisq = array(j,ichisq)
          call set_orbit_elements(orbit,m1,m2,verbose=.true.)
          if (orbit_is_parabolic(orbit%e)) then
             rp = orbit%a
          else
             rp = orbit%a*(1.-orbit%e)
          endif
          if (rp > 10000.) then
             print*,' rp = ',rp,' e = ',orbit%e,' a = ',orbit%a
             read*
          endif
          write(iu_output,"(i0,',',8(g0,','),g0)") j,rp,orbit%e,orbit%i,orbit%O,orbit%w,orbit%f,m1,m2,chisq
       enddo
       close(iu_output)
    endif
   ! print*,' got ',nlabels,' labels: ',labels(1:nlabels)
 endif

 close(iunit)
 if (itex > 0) close(itex)

end program testbin
