module damping
 implicit none

 public  :: calc_damp,apply_damp
 public  :: write_options_damping,read_options_damping

 private

 real, public :: tdyn_s      = 0.
 real, public :: damp        = 0.

contains


!-----------------------------------------------------------------------
!+
!  calculates damping factor
!+
!-----------------------------------------------------------------------
subroutine calc_damp(time, damp_fac, idamp)
 use units, only:utime
 real, intent(out)   :: damp_fac
 real, intent(in)    :: time
 integer, intent(in) :: idamp
 real                :: tau1, tau2, tdyn_star

 if (idamp == 0) then
    damp_fac = 0.
 else if (idamp == 1) then
    damp_fac = damp
 else if (idamp == 2) then
    tdyn_star = tdyn_s / utime 
    tau1 = tdyn_star * 0.1
    tau2 = tdyn_star
    if (time > 5. * tdyn_star) then
       damp_fac = 0.
    elseif (time > 2. * tdyn_star) then
       damp_fac = (tau1 * (tau2 / tau1)**((time - 2. * tdyn_star) / (3. * tdyn_star)))**(-1)
    else
       damp_fac = tau1**(-1)
    endif
 endif

end subroutine calc_damp

subroutine apply_damp(i, fextx, fexty, fextz, vxyzu, damp_fac)
 real, intent(inout) :: fextx, fexty, fextz
 real, intent(in)    :: vxyzu(:,:)
 real, intent(in)    :: damp_fac
 integer, intent(in) :: i

 fextx = fextx - damp_fac*vxyzu(1,i)
 fexty = fexty - damp_fac*vxyzu(2,i)
 fextz = fextz - damp_fac*vxyzu(3,i)

end subroutine apply_damp

!-----------------------------------------------------------------------
!+
!  applies damping to external force
!+
!-----------------------------------------------------------------------
subroutine write_options_damping(iunit, idamp)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit, idamp

 write(iunit,"(/,a)") '# options controlling damping'
 call write_inopt(idamp,'idamp','artificial damping of velocities (0=off, 1=constant, 2=star)',iunit)
 if (idamp > 0) then
    select case(idamp)
    case(1)
       call write_inopt(damp,'damp','artificial damping of velocities (if on, v=0 initially)',iunit)
    case(2)
       call write_inopt(tdyn_s,'tdyn_s','dynamical timescale of star in seconds - damping is dependent on it',iunit)
    end select
 endif

end subroutine write_options_damping



!-----------------------------------------------------------------------
!+
!  reads damping options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_damping(name,valstring,imatch,igotall,ierr,idamp)
 use io,            only:fatal
 use eos_helmholtz, only:eos_helmholtz_set_relaxflag
 character(len=*), intent(in)    :: name,valstring
 logical,          intent(out)   :: imatch,igotall
 integer,          intent(out)   :: ierr
 integer,          intent(inout) :: idamp
 integer,          save          :: ngot  = 0
 character(len=30), parameter    :: label = 'read_options_damp'
 integer :: tmp

 imatch  = .true.
 select case(trim(name))
 case('idamp')
    read(valstring,*,iostat=ierr) idamp
    ngot = ngot + 1
    if (idamp < 0) call fatal(label,'damping form choice out of range')
 case('damp')
    read(valstring,*,iostat=ierr) damp
    if (damp <= 0.)  call fatal(label,'damp <= 0: damping must be positive')
    ngot = ngot + 1
 case('tdyn_s')
    read(valstring,*,iostat=ierr) tdyn_s
    if (tdyn_s < 0.)  call fatal(label,'tdyn_s < 0: this makes no sense')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 !--make sure we have got all compulsory options (otherwise, rewrite input file)
 if (idamp > 0) then
    igotall = (ngot >= 2)
 else
    igotall = (ngot >= 1)
 endif

end subroutine read_options_damping

end module damping
