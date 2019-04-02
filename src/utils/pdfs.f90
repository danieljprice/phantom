!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: pdfs
!
!  DESCRIPTION:
!  module for probability distribution function calculation
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io
!+
!--------------------------------------------------------------------------
module pdfs
 implicit none
 public :: get_pdf,write_pdf

contains

subroutine get_pdf(xbin,pdf,nbins,npart,xpart,xmin,xmax)
 use io, only:fatal
 integer, intent(in)  :: npart,nbins
 real,    intent(in)  :: xpart(:)
 real,    intent(in)  :: xmin,xmax
 real,    intent(out) :: xbin(nbins),pdf(nbins)
 !integer, dimension(nbins) :: npartbin
 integer :: ibin,i
 real :: dx,totprob,fi,fprev,xi,xprev
 logical, parameter :: logx = .true.
 logical, parameter :: lnx  = .false.

 if (logx .and. lnx) call fatal('get_pdf','Cannot have both logx and lnx true. Pick one...')

 print*,'xmin,max = ',xmin,xmax
!
!--set up the bins in the quantity
!
 if (logx) then
    dx = (log10(xmax) - log10(xmin))/real(nbins-1)
    do ibin=1,nbins
       xbin(ibin) = xmin*10**((ibin-0.5)*dx)
    enddo
 elseif (lnx) then
    dx = (log(xmax) - log(xmin))/real(nbins-1)
    do ibin=1,nbins
       xbin(ibin) = xmin*exp((ibin-0.5)*dx)
    enddo
 else
    dx = (xmax - xmin)/(nbins - 1)
    do ibin=1,nbins
       xbin(ibin) = xmin + (ibin-0.5)*dx
    enddo
 endif

!
!--now calculate probability of finding a particle at each x
!
 pdf(:) = 0.
 do i=1,npart
    if (logx) then
       ibin = int((log10(xpart(i)) - log10(xmin))/dx) + 1
    elseif (lnx) then
       ibin = int((log(xpart(i)) - log(xmin))/dx) + 1
    else
       ibin = int((xpart(i) - xmin)/dx) + 1
    endif
    if (ibin < 1) ibin = 1
    if (ibin > nbins) ibin = nbins
    pdf(ibin) = pdf(ibin) + 1.
 enddo
!
!--get total area under pdf by trapezoidal rule
!
 totprob = 0.
 fprev = pdf(1)
 xprev = xbin(1)
 do ibin=2,nbins
    fi = pdf(ibin)
    xi = xbin(ibin)
    if (logx) then
       ! \int pdf dx = ln(10)*\int x*pdf d(log_10 x)
       totprob = totprob + 0.5*log(10.)*dx*xbin(ibin)*(fi + fprev)
    elseif (lnx) then
       ! \int pdf dx = \int x*pdf d(ln x)
       !   totprob = totprob + 0.5*dx*xbin(ibin)*(fi + fprev)
       totprob = totprob + 0.5*(fi + fprev)*(log(xi)-log(xprev))
    else
       totprob = totprob + 0.5*dx*(fi + fprev)
    endif
    fprev = fi
    xprev = xi
 enddo
!
!--normalise pdf so total area is unity
!
 print*,'normalisation factor = ',totprob,npart*dx
 pdf(1:nbins) = pdf(1:nbins)/totprob

end subroutine get_pdf

subroutine write_pdf(iunit,basename,variable,nbins,xbin,pb,time)
 integer,          intent(in) :: iunit
 character(len=*), intent(in) :: basename,variable
 integer,          intent(in) :: nbins
 real,             intent(in) :: xbin(nbins),pb(nbins)
 real,             intent(in), optional :: time
 integer :: i

 print "(a)",' writing to '//trim(basename)//'_'//trim(variable)//'.pdf'
 open(unit=iunit,file=trim(basename)//'_'//trim(variable)//'.pdf', &
      form='formatted',status='replace')
 if (present(time)) write(iunit,*) time
 do i=1,nbins
    write(iunit,*) xbin(i),pb(i)
 enddo
 close(iunit)

end subroutine write_pdf

end module pdfs
