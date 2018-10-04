!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testutils
!
!  DESCRIPTION:
!  This routine contains utility functions for use in
!  the testsuite modules
!
!  Requires mpi utility routines to print per-thread results
!
!  Also uses isdead routine from part module to identify
!  dead particles
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, mpiutils, part
!+
!--------------------------------------------------------------------------
module testutils
 use mpiutils, only:reduce_mpi
 use io,       only:id,master
 use part,     only:isdead
 implicit none
 public :: checkval,checkvalf,checkvalbuf,checkvalbuf_start,checkvalbuf_end

 private

 interface checkval
  module procedure checkvalconst,checkvalconstr4,checkvalconsti1
  module procedure checkval1_r4,checkval1_r8,checkval1_int,checkval1_int8,checkval1_logical
  module procedure checkval_r8arr,checkval_r4arr
 end interface checkval

 interface checkvalf
  module procedure checkvalfuncr8,checkvalfuncr4
 end interface

 interface checkvalbuf
  module procedure checkvalbuf_int,checkvalbuf_logical,checkvalbuf_real
 end interface checkvalbuf

 interface checkvalbuf_end
  module procedure checkvalbuf_end_int,checkvalbuf_end_real
 end interface checkvalbuf_end

 interface printerr
  module procedure printerr_real,printerr_int,printerr_logical
 end interface

 interface printresult
  module procedure printresult_real,printresult_int
 end interface

 real, parameter :: smallval = 1.e-6

contains

!----------------------------------------------------------------
!+
!  checks a constant
!+
!----------------------------------------------------------------
subroutine checkvalconst(n,x,val,tol,ndiff,label,checkmask)
 integer,          intent(in)  :: n
 real(kind=8),     intent(in)  :: x(:)
 real,             intent(in)  :: val,tol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 logical, optional,intent(in)  :: checkmask(:)
 integer      :: i
 real         :: erri,errmax

 call print_testinfo(trim(label))

 ndiff = 0
 errmax = 0.
 do i=1,n
    if (present(checkmask)) then
       if (.not. checkmask(i)) cycle
    endif
    if (.not.isdead(i)) then
       erri = abs(x(i)-val)
       if (abs(val) > epsilon(val)) erri = erri/abs(val)
       errmax = max(errmax,erri)

       if (erri > tol .or. erri /= erri) then
          ndiff = ndiff + 1
          if (ndiff==1) write(*,*)
          if (ndiff < 10) call printerr(label,real(x(i)),val,erri,tol,i)
       endif
    endif
 enddo

 call printresult(n,ndiff,errmax,tol)

 return
end subroutine checkvalconst

!----------------------------------------------------------------
!+
!  checks a constant (real4)
!+
!----------------------------------------------------------------
subroutine checkvalconstr4(n,x,val,tol,ndiff,label,checkmask)
 integer,          intent(in)  :: n
 real(kind=4),     intent(in)  :: x(:)
 real,             intent(in)  :: val,tol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 logical, optional,intent(in)  :: checkmask(:)
 integer :: i
 real    :: erri,errmax

 call print_testinfo(trim(label))

 ndiff = 0
 errmax = 0.
 do i=1,n
    if (present(checkmask)) then
       if (.not. checkmask(i)) cycle
    endif
    if (.not.isdead(i)) then
       erri = abs(x(i)-val)
       if (abs(val) > epsilon(val)) erri = erri/abs(val)
       errmax = max(errmax,erri)

       if (erri > tol .or. erri /= erri) then
          ndiff = ndiff + 1
          if (ndiff==1) write(*,*)
          if (ndiff < 10) call printerr(label,real(x(i)),val,erri,tol,i)
       endif
    endif
 enddo

 call printresult(n,ndiff,errmax,tol)

 return
end subroutine checkvalconstr4

!----------------------------------------------------------------
!+
!  checks an integer*1 array against a constant
!+
!----------------------------------------------------------------
subroutine checkvalconsti1(n,ix,ival,itol,ndiff,label,checkmask)
 integer,          intent(in)  :: n
 integer(kind=1),  intent(in)  :: ix(:)
 integer,          intent(in)  :: ival,itol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 logical, optional,intent(in)  :: checkmask(:)
 integer :: i
 integer :: erri,errmax

 call print_testinfo(trim(label))

 ndiff = 0
 errmax = 0
 do i=1,n
    if (present(checkmask)) then
       if (.not. checkmask(i)) cycle
    endif
    if (.not.isdead(i)) then
       erri = abs(ix(i)-ival)
       errmax = max(errmax,erri)

       if (erri > itol .or. erri /= erri) then
          ndiff = ndiff + 1
          if (ndiff==1) write(*,*)
          if (ndiff < 10) call printerr(label,int(ix(i)),ival,erri,itol)
       endif
    endif
 enddo

 call printresult(n,ndiff,errmax,itol)

 return
end subroutine checkvalconsti1

!----------------------------------------------------------------
!+
!  checks an array of values against a functional form
!+
!----------------------------------------------------------------
subroutine checkvalfuncr8(n,xyzhi,x,func,tol,ndiff,label,checkmask)
 integer,          intent(in)  :: n
 real,             intent(in)  :: xyzhi(:,:)
 real(kind=8),     intent(in)  :: x(:)
 real, external                          :: func
 real,             intent(in)  :: tol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 logical, optional,intent(in)  :: checkmask(:)
 integer :: i
 real(kind=8) :: erri,val,errmax
 real :: errmaxr

 call print_testinfo(trim(label))

 ndiff = 0
 errmax = 0.
 do i=1,n
    if (present(checkmask)) then
       if (.not. checkmask(i)) cycle
    endif
    if (.not.isdead(i)) then
       val = func(xyzhi(:,i))
       erri = abs(x(i)-val)
       if (abs(val) > smallval .and. erri > tol) erri = erri/abs(val)
!       if (abs(val) > tol) erri = erri/val

       if (erri > tol .or. erri /= erri) then
          ndiff = ndiff + 1
          if (ndiff==1) write(*,*)
          if (ndiff < 10 .or. erri > 2.*errmax) then
             call printerr(label,real(x(i)),real(val),real(erri),tol,i)
          endif
       endif
       errmax = max(errmax,erri)
    endif
 enddo

 errmaxr = real(errmax)
 call printresult(n,ndiff,errmaxr,real(tol))

 return
end subroutine checkvalfuncr8

!----------------------------------------------------------------
!+
!  as above but for a real*4 array
!+
!----------------------------------------------------------------
subroutine checkvalfuncr4(n,xyzhi,x,func,tol,ndiff,label,checkmask)
 integer,          intent(in)  :: n
 real,             intent(in)  :: xyzhi(:,:)
 real(kind=4),     intent(in)  :: x(:)
 real, external                         :: func
 real,             intent(in)  :: tol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 logical, optional,intent(in)  :: checkmask(:)
 integer :: i
 real    :: erri,val,errmax

 call print_testinfo(trim(label))

 ndiff = 0
 errmax = 0.
 do i=1,n
    if (present(checkmask)) then
       if (.not. checkmask(i)) cycle
    endif
    if (.not.isdead(i)) then
       val = func(xyzhi(:,i))
       erri = abs(x(i)-val)
       if (abs(val) > smallval .and. erri > tol) erri = erri/abs(val)
!       if (abs(val) > tol) erri = erri/val

       if (erri > tol .or. erri /= erri) then
          ndiff = ndiff + 1
          if (ndiff==1) write(*,*)
          if (ndiff < 10 .or. erri > 2.*errmax) then
             call printerr(label,real(x(i)),val,erri,tol,i)
          endif
       endif
       errmax = max(errmax,erri)
    endif
 enddo

 call printresult(n,ndiff,errmax,real(tol))

 return
end subroutine checkvalfuncr4

!----------------------------------------------------------------
!+
!  checks a single, scalar value
!+
!----------------------------------------------------------------
subroutine checkval1_r4(xi,val,tol,ndiff,label)
 real(kind=4),     intent(in)  :: xi
 real(kind=4),     intent(in)  :: val,tol
 character(len=*), intent(in)  :: label
 integer,          intent(out) :: ndiff
 real(kind=4) :: erri

 ndiff = 0
 call print_testinfo(trim(label))

 erri = abs(xi-val)
 if (abs(val) > smallval) erri = erri/abs(val)

 if (erri > tol .or. erri /= erri) then
    write(*,"(4(a,es10.3),a)") 'FAILED [got ',xi,' should be  ',val,', err =',erri,', tol =',tol,']'
    ndiff = 1
 else
    if (id==master)  write(*,"(a,2(es10.3,a))") 'OK     [max err =',erri,', tol =',tol,']'
    !call printresult(0,0,real(erri),real(tol))
 endif

 return
end subroutine checkval1_r4

!----------------------------------------------------------------
!+
!  checks a single, scalar value
!+
!----------------------------------------------------------------
subroutine checkval1_r8(xi,val,tol,ndiff,label)
 real(kind=8),     intent(in)  :: xi
 real(kind=8),     intent(in)  :: val,tol
 character(len=*), intent(in)  :: label
 integer,          intent(out) :: ndiff
 real(kind=8) :: erri

 ndiff = 0
 call print_testinfo(trim(label))

 erri = abs(xi-val)
 if (abs(val) > smallval) erri = erri/abs(val)

 if (erri > tol .or. erri /= erri) then
    write(*,"(4(a,es10.3),a)") 'FAILED [got ',xi,' should be  ',val,', err =',erri,', tol =',tol,']'
    ndiff = 1
 else
    if (id==master) write(*,"(a,2(es10.3,a))") 'OK     [max err =',erri,', tol =',tol,']'
    !call printresult(0,0,real(erri),real(tol))
 endif

 return
end subroutine checkval1_r8

!----------------------------------------------------------------
!+
!  checks a single, integer value
!+
!----------------------------------------------------------------
subroutine checkval1_logical(ix,ival,ndiff,label)
 logical,          intent(in)  :: ix
 logical,          intent(in)  :: ival
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label

 ndiff = 0
 call print_testinfo(trim(label))

 if (ix.eqv.ival) then
    if (id==master) write(*,"(a,l1,a,l1,a)") 'OK     [got ',ix,' should be ',ival,']'
 else
    write(*,"(a,l1,a,l1,a)") 'FAILED [got ',ix,' should be ',ival,']'
    ndiff = 1
 endif

 return
end subroutine checkval1_logical

!----------------------------------------------------------------
!+
!  checks a single, integer value
!+
!----------------------------------------------------------------
subroutine checkval1_int(ix,ival,itol,ndiff,label)
 integer,          intent(in)  :: ix
 integer,          intent(in)  :: ival,itol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 integer :: erri !,ncheck

 ndiff = 0
 call print_testinfo(trim(label))

 erri = abs(ix-ival)
 !if (erri > itol) ndiff = 1
 !ncheck = 1
 !call printresult(ncheck,ndiff,erri,itol)
 !return

 if (erri > itol) then
    write(*,"(a,i11,a,i12,a)") 'FAILED [got',ix,' should be',ival,']'
    ndiff = 1
 elseif (id==master) then
    if (itol==0) then
       write(*,"(a,i11,a,i12,a)") 'OK     [got',ix,' should be',ival,']'
    else
       write(*,"(a,i11,a,i12,a,i5,a)") 'OK     [got',ix,' should be',ival,', tol = ',itol,']'
    endif
 endif

 return
end subroutine checkval1_int

!----------------------------------------------------------------
!+
!  checks a single, integer*8 value
!+
!----------------------------------------------------------------
subroutine checkval1_int8(ix,ival,itol,ndiff,label)
 integer(kind=8),  intent(in)  :: ix
 integer(kind=8),  intent(in)  :: ival
 integer,          intent(in)  :: itol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 integer(kind=8) :: erri

 ndiff = 0
 call print_testinfo(trim(label))

 erri = abs(ix-ival)
 if (erri > itol) then
    write(*,"(a,i11,a,i12,a)") 'FAILED [got',ix,' should be',ival,']'
    ndiff = 1
 elseif (id==master) then
    if (itol==0) then
       write(*,"(a,i11,a,i12,a)") 'OK     [is ',ix,' should be',ival,']'
    else
       write(*,"(a,i11,a,i12,a,i5,a)") 'OK     [is ',ix,' should be',ival,', tol = ',itol,']'
    endif
 endif

 return
end subroutine checkval1_int8

!----------------------------------------------------------------
!+
!  checks an array of values against an array of expected answers
!+
!----------------------------------------------------------------
subroutine checkval_r8arr(n,x,xexact,tol,ndiff,label,checkmask)
 integer,          intent(in)  :: n
 real(kind=8),     intent(in)  :: x(:),xexact(:)
 real,             intent(in)  :: tol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 logical, optional,intent(in)  :: checkmask(:)
 integer :: i,nval
 real(kind=8) :: erri,val,errmax,valmax,errl2
 real :: errmaxr

 call print_testinfo(trim(label))

 ndiff = 0
 errmax = 0.
 errl2  = 0.
 valmax = 0.
 nval = 0
 do i=1,n
    if (present(checkmask)) then
       if (.not. checkmask(i)) cycle
    endif
    if (.not.isdead(i)) then
       val = xexact(i)
       erri = abs(x(i)-val)
       errl2 = errl2 + erri*erri
       valmax = max(val,valmax)
       if (abs(val) > smallval .and. erri > tol) erri = erri/abs(val)
!       if (abs(val) > tol) erri = erri/val

       if (erri > tol .or. erri /= erri) then
          ndiff = ndiff + 1
          if (ndiff==1) write(*,*)
          if (ndiff < 10 .or. erri > 2.*errmax) then
             call printerr(label,real(x(i)),real(val),real(erri),tol,i)
          endif
       endif
       nval = nval + 1
       errmax = max(errmax,erri)
    endif
 enddo

 errmaxr = real(errmax)
 call printresult(n,ndiff,errmaxr,real(tol),real(errl2),real(valmax),nval)

 return
end subroutine checkval_r8arr

!----------------------------------------------------------------
!+
!  checks an array of real*4 values against an array of expected answers
!+
!----------------------------------------------------------------
subroutine checkval_r4arr(n,x,xexact,tol,ndiff,label,checkmask)
 integer,          intent(in)  :: n
 real(kind=4),     intent(in)  :: x(:),xexact(:)
 real,             intent(in)  :: tol
 integer,          intent(out) :: ndiff
 character(len=*), intent(in)  :: label
 logical, optional,intent(in)  :: checkmask(:)
 integer :: i
 real(kind=4) :: erri,val,errmax
 real :: errmaxr

 call print_testinfo(trim(label))

 ndiff = 0
 errmax = 0.
 do i=1,n
    if (present(checkmask)) then
       if (.not. checkmask(i)) cycle
    endif
    if (.not.isdead(i)) then
       val = xexact(i)
       erri = abs(x(i)-val)
       if (abs(val) > smallval .and. erri > tol) erri = erri/abs(val)
!       if (abs(val) > tol) erri = erri/val

       if (erri > tol .or. erri /= erri) then
          ndiff = ndiff + 1
          if (ndiff==1) write(*,*)
          if (ndiff < 10 .or. erri > 2.*errmax) then
             call printerr(label,real(x(i)),real(val),real(erri),tol,i)
          endif
       endif
       errmax = max(errmax,erri)
    endif
 enddo

 errmaxr = errmax
 call printresult(n,ndiff,errmaxr,real(tol))

 return
end subroutine checkval_r4arr

!----------------------------------------------------------------
!+
!  start a buffered error check
!+
!----------------------------------------------------------------
subroutine checkvalbuf_start(label)
 character(len=*), intent(in) :: label

 call print_testinfo(trim(label))
 write(*,"(a)")

 return
end subroutine checkvalbuf_start

!----------------------------------------------------------------
!+
!  checks a single, integer value
!  (buffered: reports on errors only and ndiff is a running total)
!+
!----------------------------------------------------------------
subroutine checkvalbuf_int(ix,ival,itol,label,ndiff,ncheck,ierrmax)
 integer,          intent(in)    :: ix
 integer,          intent(in)    :: ival,itol
 character(len=*), intent(in)    :: label
 integer,          intent(inout) :: ndiff,ncheck
 integer,          intent(inout), optional :: ierrmax
 integer :: erri

 erri = abs(ix-ival)
 ncheck = ncheck + 1
 if (erri > itol) then
    ndiff = ndiff + 1
    if (ndiff < 10 .or. erri > 2*ierrmax) call printerr(label,ix,ival,erri,itol)
 endif
 if (present(ierrmax)) ierrmax = max(ierrmax,erri)

 return
end subroutine checkvalbuf_int

!----------------------------------------------------------------
!+
!  checks a single, real value
!  (buffered: reports on errors only and ndiff is a running total)
!+
!----------------------------------------------------------------
subroutine checkvalbuf_real(xi,val,tol,label,ndiff,ncheck,errmax)
 real,             intent(in)    :: xi
 real,             intent(in)    :: val,tol
 character(len=*), intent(in)    :: label
 integer,          intent(inout) :: ndiff,ncheck
 real,             intent(inout) :: errmax
 real :: erri

 erri = abs(xi-val)
 if (abs(val) > smallval .and. erri > tol) erri = erri/abs(val)

 ncheck = ncheck + 1
 if (erri > tol .or. erri /= erri) then
    ndiff = ndiff + 1
    if (ndiff < 10 .or. erri > 2.*errmax) call printerr(label,xi,val,erri,tol)
 endif
 errmax = max(errmax,erri)

 return
end subroutine checkvalbuf_real

!----------------------------------------------------------------
!+
!  checks a logical value
!  (buffered: reports on errors only and ndiff is a running total)
!+
!----------------------------------------------------------------
subroutine checkvalbuf_logical(lx,lval,label,ndiff,ncheck)
 logical,          intent(in)    :: lx,lval
 character(len=*), intent(in)    :: label
 integer,          intent(inout) :: ndiff,ncheck

 ncheck = ncheck + 1
 if (lval.neqv.lx) then
    ndiff = ndiff + 1
    if (ndiff < 10) call printerr(label,lx,lval)
 endif

 return
end subroutine checkvalbuf_logical

!----------------------------------------------------------------
!+
!  end a buffered error check (int)
!+
!----------------------------------------------------------------
subroutine checkvalbuf_end_int(label,n,ndiff,ierrmax,itol,ntot)
 character(len=*), intent(in) :: label
 integer,          intent(in) :: n,itol
 integer,          intent(inout) :: ndiff,ierrmax
 integer,          intent(in), optional :: ntot

 call print_testinfo(trim(label))
 if (present(ntot)) then
    call printresult(n,ndiff,ierrmax,itol,ntot)
 else
    call printresult(n,ndiff,ierrmax,itol)
 endif

 return
end subroutine checkvalbuf_end_int

!----------------------------------------------------------------
!+
!  end a buffered error check (real)
!+
!----------------------------------------------------------------
subroutine checkvalbuf_end_real(label,n,ndiff,errmax,tol)
 character(len=*), intent(in) :: label
 integer,          intent(in)    :: n
 integer,          intent(inout) :: ndiff
 real,             intent(inout) :: errmax
 real,             intent(in)    :: tol

 call print_testinfo(trim(label))
 call printresult(n,ndiff,errmax,tol)

 return
end subroutine checkvalbuf_end_real

!----------------------------------------------------------------
!+
!  formatting for printing errors in test results
!+
!----------------------------------------------------------------
subroutine printerr_real(label,x,val,erri,tol,i)
 character(len=*), intent(in) :: label
 real,             intent(in) :: x, val, erri, tol
 integer,          intent(in), optional :: i

 if (abs(val) > smallval) then
    if (present(i)) then
       write(*,"(1x,4(a,es10.3),a,i10,a)") &
            trim(label)//' = ',x,' should be ',val,' ratio =',x/val,' err =',erri,' (',i,')'
    else
       write(*,"(1x,5(a,es10.3),a)") &
            trim(label)//' = ',x,' should be ',val,' ratio =',x/val,' err =',erri,' (tol =',tol,')'
    endif
 else
    if (present(i)) then
       write(*,"(1x,3(a,es10.3),a,i10,a)") &
            trim(label)//' = ',x,' should be ',val,' err =',erri,' (',i,')'
    else
       write(*,"(1x,4(a,es10.3),a)") &
            trim(label)//' = ',x,' should be ',val,' err =',erri,' (tol =',tol,')'
    endif
 endif

 return
end subroutine printerr_real

!----------------------------------------------------------------
!+
!  formatting for printing errors in test results
!+
!----------------------------------------------------------------
subroutine printerr_int(label,ix,ival,erri,itol)
 character(len=*), intent(in) :: label
 integer,          intent(in) :: ix, ival, erri, itol

 if (itol > 0) then
    write(*,"(1x,4(a,i10),a)") &
      trim(label)//' = ',ix,' should be ',ival,' err =',erri,' (tol =',itol,')'
 else
    write(*,"(1x,3(a,i10),a)") &
      trim(label)//' = ',ix,' should be ',ival,' err =',erri
 endif

 return
end subroutine printerr_int

!----------------------------------------------------------------
!+
!  formatting for printing errors in test results
!+
!----------------------------------------------------------------
subroutine printerr_logical(label,lx,lval)
 character(len=*), intent(in) :: label
 logical,          intent(in) :: lx, lval

 write(*,"(1x,2(a,l1))") 'ERROR! '//trim(label)//' is ',lx,' should be ',lval

 return
end subroutine printerr_logical

!----------------------------------------------------------------
!+
!  formatting for initial test information
!+
!----------------------------------------------------------------
subroutine print_testinfo(string)
 character(len=*), intent(in) :: string
 character(len=20) :: fmtstring
 integer           :: ndots,istart

 if (id==master) then
    istart = 20
    ndots  = -1
    do while(ndots < 2 .and. istart < 60)
       istart = istart + 10
       ndots = istart - len_trim(string)
    enddo
    if (ndots < 2) ndots = 3

    write(fmtstring,"('(1x,a,',i2,'(''.''))')") ndots
    write(*,fmtstring,ADVANCE='NO') 'checking '//trim(string)
 endif

end subroutine print_testinfo

!----------------------------------------------------------------
!+
!  formatting for printing test results
!+
!----------------------------------------------------------------
subroutine printresult_real(npi,ndiff,errmax,tol,errl2i,valmaxi,nvali)
 integer, intent(in)    :: npi
 integer, intent(inout) :: ndiff
 real,    intent(inout) :: errmax
 real,    intent(in)    :: tol
 real,    intent(in), optional :: errl2i,valmaxi
 integer, intent(in), optional :: nvali
 integer(kind=8) :: np,nval
 real            :: valmax,errl2

 np     = reduce_mpi('+',npi)
 ndiff  = int(reduce_mpi('+',ndiff))
 errmax = reduce_mpi('max',errmax)

 if (present(errl2i)) then
    errl2 = reduce_mpi('+',errl2i)
    if (present(valmaxi) .and. present(nvali)) then
       valmax = reduce_mpi('max',valmaxi)
       nval   = reduce_mpi('+',nvali)
       if (nval > 0 .and. valmax > 0.) then
          errl2 = sqrt(errl2/(real(nval)*valmax*valmax))
       endif
    else
       errl2 = sqrt(errl2)
    endif
 endif

 if (id==master) then
    if (ndiff==0) then
       if (present(errl2i)) then
          write(*,"(a,3(es10.3,a))") 'OK     [max err =',errmax,', L2 err = ',errl2,' tol =',tol,']'
       else
          write(*,"(a,2(es10.3,a))") 'OK     [max err =',errmax,', tol =',tol,']'
       endif
    elseif (ndiff > 0) then
       if (present(errl2i)) then
          write(*,"(1x,2(a,i10),3(a,es10.3),a)") &
             'FAILED [on ',ndiff,' of ',np,' values, max err =',errmax,', L2 err = ',errl2,' tol =',tol,']'
       else
          write(*,"(1x,2(a,i10),2(a,es10.3),a)") 'FAILED [on ',ndiff,' of ',np,' values, max err =',errmax,', tol =',tol,']'
       endif
    else ! used for single values
       write(*,"(1x,a,es10.3,a,es10.3,a)") 'FAILED [max err =',errmax,', tol =',tol,']'
    endif
 endif

 return
end subroutine printresult_real

!----------------------------------------------------------------
!+
!  formatting for printing test results
!+
!----------------------------------------------------------------
subroutine printresult_int(nchecki,ndiff,ierrmax,itol,ntot)
 integer, intent(in)    :: nchecki
 integer, intent(inout) :: ndiff
 integer, intent(inout) :: ierrmax
 integer, intent(in)    :: itol
 integer, intent(in), optional :: ntot
 integer(kind=8) :: ncheck

 ncheck  = reduce_mpi('+',nchecki)
 ndiff   = int(reduce_mpi('+',ndiff))
 ierrmax = int(reduce_mpi('max',ierrmax))

 if (id==master) then
    if (ndiff==0) then
       if (ierrmax > 0) then
          write(*,"(a,i5,a,i2,a)") 'OK     [max err =',ierrmax,', tol =',itol,']'
       elseif (ncheck > 0) then
          if (present(ntot)) then
             if (ntot < 1e6 .and. ncheck < 1e6) then
                write(*,"(2(a,i5),a)")  'OK     [checked ',ncheck,' of ',ntot,' values]'
             else
                write(*,"(2(a,i10),a)") 'OK     [checked ',ncheck,' of ',ntot,' values]'
             endif
          else
             if (ncheck < 1e6) then
                write(*,"(a,i5,a)") 'OK     [checked ',ncheck,' values]'
             else
                write(*,"(a,i10,a)") 'OK     [checked ',ncheck,' values]'
             endif
          endif
       else
          write(*,"(a)") 'OK'
       endif
    elseif (ndiff > 0) then
       write(*,"(2(a,i10),a,i10,a)") 'FAILED [on ',ndiff,' of ',ncheck,' values, max err =',ierrmax,']'
    else ! this is used for single values
       write(*,"(1x,a,i5,a,i2,a)") 'FAILED [max err =',ierrmax,', tol =',itol,']'
    endif
 endif

 return
end subroutine printresult_int

end module testutils
