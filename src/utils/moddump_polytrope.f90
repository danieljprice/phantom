module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use prompting, only:prompt
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 real    :: vr,rhat(3),r,vmax,vmin,vmean,vi

!
!-- Find min, max, and mean velocties in star/polytrope
!
 vmin  = huge(vmin)
 vmax  = 0.
 vmean = 0.
 do i=1,npart
    vi    = sqrt(dot_product(vxyzu(1:3,i),vxyzu(1:3,i)))
    vmin  = min(vmin,vi)
    vmax  = max(vmax,vi)
    vmean = vmean + vi
 enddo
 vmean = vmean/npart

 print*,'--- Velocities in polytrope: ---'
 print*,'vmin  = ',vmin
 print*,'vmax  = ',vmax
 print*,'vmean = ',vmean
 print*,'--------------------------------'

!
!-- Set velocties in star/polytrope to desried value
!
 vr = -vmean
 call prompt('Enter the velocity you want the star to begin oscillating with',vr)
 do i=1,npart
    r     = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    rhat  = xyzh(1:3,i)/r
    vxyzu(1:3,i) = vr*rhat
 enddo

 return
end subroutine modify_dump

end module moddump
