module testcons2prim
 implicit none
contains

subroutine test_cons2prim_i(x,v,dens,u,p,ntests,npass)
 use cons2prim,    only:conservative_to_primitive,primitive_to_conservative
 use testutils,    only:checkval,checkvalbuf
 use testmetric,   only:test_metric_i
 use metric_tools, only:get_grpacki
 real, intent(in) :: x(1:3),v(1:3),dens,u,p
 integer, intent(inout) :: ntests,npass
 real :: grpacki(0:3,0:3,5)
 real :: rho,pmom(1:3),en
 real :: v_out(1:3),dens_out,u_out,p_out
 real, parameter :: tol = 4.e-12
 integer :: nerrors, ierr, j
 integer :: ncheck,dummy
 real :: errmax

 ntests = ntests + 1
 nerrors = 0

 ! Used for initial guess in conservative2primitive
 v_out    = v
 dens_out = dens
 u_out    = u
 p_out    = p

 call primitive_to_conservative(x,v,dens,u,P,rho,pmom,en)
 call get_grpacki(x,grpacki)
 call conservative_to_primitive(x,grpacki,v_out,dens_out,u_out,p_out,rho,pmom,en,ierr)

 ! call checkval(ierr,0,0,n_error,'ierr = 0 for convergence')
 call checkvalbuf(ierr,0,0,'[F]: ierr (convergence)',nerrors,ncheck)
 ! nerrors = nerrors + n_error

 ! call checkval(3,v_out,v,tol,n_error,'v_out = v')
 do j=1,3
    call checkvalbuf(v_out(j),v(j),tol,'[F]: v_out',nerrors,ncheck,errmax)
    ! nerrors = nerrors + n_error
 enddo

 ! call checkval(dens_out,dens,tol,n_error,'dens_out = dens')
 call checkvalbuf(dens_out,dens,tol,'[F]: dens_out',nerrors,ncheck,errmax)
 ! nerrors = nerrors + n_error

 ! call checkval(u_out,u,tol,n_error,'u_out = u')
 call checkvalbuf(u_out,u,tol,'[F]: u_out',nerrors,ncheck,errmax)
 ! nerrors = nerrors + n_error

 ! call checkval(p_out,p,tol,n_error,'p_out = p')
 call checkvalbuf(p_out,p,tol,'[F]: p_out',nerrors,ncheck,errmax)
 ! nerrors = nerrors + n_error

 if (nerrors/=0) then
    print*,'-- cons2prim test failed with'
    print*,'  - IN:'
    print*,'     x    =',x
    print*,'     v    =',v
    print*,'     dens =',dens
    print*,'     u    =',u
    print*,'     p    =',p
    print*,'  - OUT:'
    print*,'     v    =',v_out
    print*,'     dens =',dens_out
    print*,'     u    =',u_out
    print*,'     p    =',p_out
    print*,''
 else
    npass = npass + 1
 endif
end subroutine test_cons2prim_i

end module testcons2prim
