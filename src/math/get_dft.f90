!!
!> DFT/IDFT Subroutine
!!
subroutine getDFT(pi,im,j,n,pij,ipij,kk,ikk)

  implicit none

  real*8,     intent(in)                  :: pi
  complex*16, intent(in)                  :: im
  integer,    intent(in)                  :: j,n
  integer,    dimension(n), intent(inout) :: kk
  complex*16, dimension(n), intent(inout) :: pij,ipij
  complex*16,               intent(out)   :: ikk
  !
  integer i,c
  complex*16 :: w
  !
  w = dcmplx(2.0d0*pi/dble(n))
  c    = 0
  if (mod(n,2).eq.0) then
     do i = -(n/2-1),n/2
        c = c + 1
        kk(c) = i
     enddo
  else
     do i = -(n-1)/2,(n-1)/2
        c = c + 1
        kk(c) = i
     enddo
  endif
!>
!> Store DFT in complex pij and IDFT in complex ipij
!>

!> N Even 
  if (mod(n,2).eq.0) then
!> Make IDFT 
     if (j.eq.n) then
        do i=1,n
           ipij(i)= cos(pi*(i-1))
        enddo
     else
        do i=1,n
           ipij(i)= cexp(  im*w*real(kk(j))*real(i-1) )
        enddo
     endif
!> DFT
     do i = 1,n-1
        pij(i) = cexp( -im*w*real(kk(i))*real(j-1) ) / real(n)
     enddo
     pij(n) = cos(pi*(j-1))/n
!> N Odd
  else
     do i=1,n
        pij(i) = cexp( -im*w*real(kk(i))*real(j-1) ) / real(n)
        ipij(i)= cexp(  im*w*real(kk(j))*real(i-1) )
     enddo
  endif

  ikk = dcmplx(real(kk(j)))
  !> Zero out derivative for oddball frequency (N/2) for N even if j = n
  if (j.eq.n.and.mod(n,2).eq.0) then
     ikk = dcmplx(0.)
  endif

  return
end subroutine getDFT
