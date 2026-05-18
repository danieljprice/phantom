
! file: check_masunaga_vs_maxvals.f90
! Compares SPHNG maxvals.out columns (7=density [g/cm^3], 8=T [K])
! against Masunaga 2000 (digitized): (log10 rho, log10 T).
!
! For each row of maxvals.out, compute:
!   x_sim  = log10(density)
!   T_sim  = column 8
! Interpolate y_M(x) = log10(T_M) from Masunaga at x=x_sim (linear)
!   T_M = 10.0**y_M
! Check |T_sim - T_M|/T_sim <= tolerance
!
! Outputs a per-point result and a summary at the end.

program check_masunaga_vs_maxvals
  implicit none
  integer,parameter :: dp = kind(1.0d0)
  character(len=*), parameter :: f_masunaga   = "Masunaga_2000_collapse.csv"
  character(len=*), parameter :: f_maxvals   = "maxvals.out"   ! or "maxvals.out.txt" if that's your filename
  integer,parameter :: MAXN = 100000
  integer,parameter :: MAXMASU = 183
  ! Masunaga arrays (log10 rho, log10 T)
  real(dp) :: xM(MAXMASU), yM(MAXMASU)
  integer :: nM

  ! From maxvals: we will read only columns we need
  real(dp) :: cols(20)  ! there are up to 20 columns per header
  real(dp) :: rho,Tsim,xsim,yinterp,logTM,logTsim,relerr,tolerance
  integer  :: iostat,i,nOK,nFail,nTot
  logical  :: ok
  integer  :: u_csv, u_max

  tolerance = 13. ! per cent
  call read_masunaga(f_masunaga, xM, yM, nM)

  if (nM < 2) then
     write(*,*) "ERROR: Not enough points in ", trim(f_masunaga)
     stop 1
  end if
  write(*,'(A,I0)') "Masunaga points read: ", nM

  open(newunit=u_max, file=f_maxvals, status='old', action='read', iostat=iostat)
  if (iostat /= 0) then
     ! try alternative filename (as uploaded)
     open(newunit=u_max, file=f_maxvals//'.txt', status='old', action='read', iostat=iostat)
     if (iostat /= 0) then
        write(*,*) "ERROR: cannot open maxvals file: ", trim(f_maxvals), " or ", trim(f_maxvals)//".txt"
        stop 2
     end if
  end if

  nOK=0; nFail=0; nTot=0
  write(*,'(A,F4.1,A)') " idx     log10rho_sim     T_sim[K]    T_Masunaga[K]    rel_err(%)   within ", tolerance, "%"
  write(*,'(A)') "--------------------------------------------------------------------------"

  do
     read(u_max,*, iostat=iostat) cols
     if (iostat /= 0) then
        if (iostat < 0) exit   ! EOF
        cycle                   ! probably a short line; skip
     end if
     ! Skip comment lines (start with '#') by checking if col count failed;
     ! Our list-directed read above will fail on comment lines and continue.

     ! According to header: [07]=density [g/cm^3], [08]=temperature [K]
     rho  = cols(7)
     Tsim = cols(8)

     ! guard against non-physical entries
     if (rho <= 1.0e-17_dp .or. Tsim <= 0.0_dp) cycle

     xsim = log10(rho)

     ok = interp_linear(xM, yM, nM, xsim, yinterp)
     if (.not. ok) cycle   ! x outside tabulated range; skip or handle as you prefer

     !compare log10(temperature)
     logTM = yinterp
     logTsim = log10(cols(8))
     relerr = abs(logTsim - logTM) / logTsim * 100.0_dp

     nTot = nTot + 1
     if (relerr <= tolerance) then
        nOK = nOK + 1
        write(*,'(I5,1X,F14.6,1X,ES12.5,1X,ES14.5,1X,F10.4,1X,A)') &
             nTot, xsim, Tsim, 10**logTM, relerr, "YES"
     else
        nFail = nFail + 1
        write(*,'(I5,1X,F14.6,1X,ES12.5,1X,ES14.5,1X,F10.4,1X,A)') &
             nTot, xsim, Tsim, 10**logTM, relerr, "NO"
     end if
  end do
  close(u_max)

  write(*,*)
  write(*,'(A,I0)') "Checked points (within Masunaga x-range): ", nTot
  if (nFail == 0) then
	write(*,'(A)') "PASS"
	write(*,'(A,F4.1,A)') "All points within ", tolerance, "%"
  else
    write(*,'(A,F4.1,A)') "FAIL"
    write(*,'(A,F4.1,A,I0)') "Within ", tolerance," % : ", nOK
    write(*,'(A,F4.1,A,I0)') "Outside ", tolerance," % : ", nFail
    stop 1
  endif
contains

  subroutine read_masunaga(fname, x, y, n)
    character(len=*), intent(in) :: fname
    real(dp),intent(out) :: x(MAXMASU), y(MAXMASU)
    integer, intent(out) :: n
    integer :: u, ios, i
    character(len=80) :: line
    
    n = 0
    print *, "opening file", fname
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       print *, 'check_masunaga_vs_maxvals','file not found:',fname
       stop 1
    endif

    do
       read(u,'(A)', iostat=ios) line
       if (ios /= 0) exit
       call strip_spaces(line)
       if (len_trim(line) == 0) cycle
       if (line(1:1) == "#" .or. line(1:1) == "!") cycle
       n = n + 1
       ! Replace comma with space (robust list-directed read)
       call replace_char(line, ',', ' ')
       read(line, *, iostat=ios) x(n), y(n)
       if (ios /= 0) cycle
    end do
    close(u)

    if (n < 1) then
       print *, "check_masunaga_vs_maxvals","ERROR: no valid rows in ",trim(fname)
       stop
    end if

  end subroutine read_masunaga

  logical function interp_linear(x, y, n, xq, yq) result(ok)
    ! Piecewise-linear interpolation on (x,y). Assumes x is monotonic (increasing).
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in)  :: n
    real(dp), intent(in) :: xq
    real(dp), intent(out):: yq
    integer :: i, lo, hi
    ok = .false.
    if (n < 2) return

    ! ensure increasing order; if not, we can handle decreasing
    if (x(1) <= x(n)) then
       if (xq < x(1) .or. xq > x(n)) return
       lo = 1; hi = n
       ! binary search
       do
          if (hi - lo <= 1) exit
          i = (hi + lo)/2
          if (xq >= x(i)) then
             lo = i
          else
             hi = i
          end if
       end do
       yq = y(lo) + (y(hi)-y(lo)) * ((xq - x(lo)) / (x(hi)-x(lo)))
       ok = .true.
    else
       ! decreasing x
       if (xq > x(1) .or. xq < x(n)) return
       lo = 1; hi = n
       do
          if (hi - lo <= 1) exit
          i = (hi + lo)/2
          if (xq <= x(i)) then
             lo = i
          else
             hi = i
          end if
       end do
       yq = y(lo) + (y(hi)-y(lo)) * ((xq - x(lo)) / (x(hi)-x(lo)))
       ok = .true.
    end if
  end function interp_linear

  subroutine strip_spaces(s)
    character(len=40), intent(inout) :: s
    integer :: i
    ! trim only; keep inner spaces
    s = trim(adjustl(s))
  end subroutine strip_spaces

  subroutine replace_char(s, a, b)
    character(len=*), intent(inout) :: s
    character(len=1), intent(in)    :: a, b
    integer :: i
    do i = 1, len_trim(s)
       if (s(i:i) == a) s(i:i) = b
    end do
  end subroutine replace_char

end program check_masunaga_vs_maxvals
