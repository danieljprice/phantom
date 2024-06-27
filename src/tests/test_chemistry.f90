module testchemistry
!
! Unit tests of the equilibrium chemistry module
!
! :References: Bermudez-Bustamante et al. (2024)
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: 
!
 implicit none

contains
 
subroutine test_chemistry(ntests,npass)
 use chemistry, only:network,ncols,write_time_file,columns,nMolecules,iS,iSi,iSiO2,iSiS,iSiO
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 integer, intent(inout) :: ntests,npass
 integer, parameter :: ncheck = 5, ngrid = 100
 real :: T,rho_cgs,mu,gamma,abundance(ncols)
 real :: abundance_table(ncheck,ngrid),temp_grid(ngrid)
 real :: Tmin,Tmax,dT
 integer, parameter :: icols_to_check(ncheck) = [nMolecules-1+iS,&
                                                 nMolecules-1+iSi,&
                                                 iSiO2,&
                                                 iSiS,&
                                                 iSiO]
 integer :: num,i,j,k,nfailed(1)
 real, parameter :: tol = 1.e-2

 if (id==master) write(*,"(/,a,/)") '--> TESTING CHEMISTRY MODULE'

 rho_cgs = 2e-14
 num = 0
 Tmin = 600.
 Tmax = 2000.
 dT = (Tmax - Tmin)/(ngrid-1)

 do i=1,ngrid
    T = Tmin + (i-1)*dT
    temp_grid(i) = T
    call network(T,rho_cgs,mu,gamma,abundance,pressure_cgs=1.e-4)
    !print*,i,T,mu,gamma,abundance(1:3)

    ! Write the values to the file
    call write_time_file('abundances', columns, T, abundance, ncols, num)
    abundance_table(1:ncheck,i) = abundance(icols_to_check)
    num = num +1
 enddo

 !print*,'# T,S,Si,SiO2,SiS,SiO'
 !do i=1,ngrid
 !   print*,temp_grid(i),abundance_table(:,i)
 !enddo
 
 do j=1,ncheck
    k = icols_to_check(j)
   print*,'CHECKING '//trim(columns(k)),k
    call check_answers_compared_to_gail(ngrid,temp_grid,abundance_table(j,:),columns(k))
 enddo
 !  call checkval(ncols,abundance,ref_abundance,tol,nfailed(1),'abundances match Gail/Sedlmyer')
 !  call update_test_scores(ntests,nfailed(1:1),npass)
 !  enddo

 if (id==master) write(*,"(/,a)") '<-- CHEMISTRY TEST COMPLETE'

end subroutine test_chemistry

subroutine check_answers_compared_to_gail(ngrid,temp_grid,abundance,label)
 use table_utils, only:yinterp
 integer, intent(in)  :: ngrid
 real,    intent(in)  :: temp_grid(ngrid)
 real,    intent(out) :: abundance(ngrid)
 character(len=*), intent(in) :: label
 integer :: npts,i
 real, allocatable :: abundance_values(:),temp_values(:)
 real :: ref_abundance

 call read_gail_reference_abundances(npts,temp_values,abundance_values,label)

 do i=1,ngrid
    !ref_abundance = yinterp(abundance(:),temp_grid(:),temp_values(i))
    write(1,*) temp_grid(i),abundance(i)
    !write(2,*) temp_values(i),abundance_values(i)
    !print*,' HERE T= ',temp_values(i),'OURS=',ref_abundance,'GAIL=',abundance_values(i)
    !read*
 enddo
 stop
 !ref_abundance = 1. ! abundance

end subroutine check_answers_compared_to_gail

subroutine read_gail_reference_abundances(npts,temp_values,abundance_values,label)
 integer, intent(out) :: npts
 real, allocatable, intent(out) :: abundance_values(:),temp_values(:)
 character(len=*), intent(in) :: label
 integer :: i,iu,ierr

 npts = 100
 allocate(abundance_values(npts),temp_values(npts))
 open(newunit=iu,file=trim(label(2:))//'_abundance_vs_temperature.csv',status='old',action='read')
 do i=1,npts
    read(iu,*,iostat=ierr) temp_values(i),abundance_values(i)
   ! print*,i,' GOT ',temp_values(i),abundance_values(i)
 enddo
 close(iu)

end subroutine read_gail_reference_abundances

end module testchemistry
