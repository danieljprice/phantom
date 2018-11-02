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
 real    :: amp

 amp = 1.e-4

 call prompt('Enter the velocity amplitude you want the star to begin oscillating with',amp)

 do i=1,npart
    vxyzu(1:3,i) = amp*xyzh(1:3,i)
 enddo

 return
end subroutine modify_dump

end module moddump
