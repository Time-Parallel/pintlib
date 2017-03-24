!!
!> Real-valued dft
!!
subroutine dftReal(pi,j,n,dft_j,idft_j,kk,ikk)
  !
  implicit none
  !
  real*8,  intent(in)                  :: pi
  integer, intent(in)                  :: j,n
  integer, dimension(n), intent(inout) :: kk
  real*8,  dimension(n), intent(inout) :: dft_j,idft_j
  real*8,                intent(out)   :: ikk
  !
  integer :: i,c
  real*8  :: w
  !
  w  = 2.0d0*pi/dble(n)
  c  = 1
  kk = 0
  if (mod(n,2).eq.0) then
     do i = 1,n/2-1
        kk(c+1) = i
        kk(c+2) = i
        c = c + 2
     enddo
  else
      do i = 1,(n-1)/2
         kk(c+1) = i
         kk(c+2) = i
         c = c +2
     enddo
  endif
  write(*,*) "kk: ", kk

!>
!> Store DFT in real pij and IDFT in real ipij
!>

!> N Even 
  if (mod(n,2).eq.0) then
     dft_j = 0.0d0
     idft_j = 0.0d0
  else
     ! IDFT
     if (mod(j,2).eq.0) then
        do i = 0,n-1
           idft_j(i+1) = sin(w*real(i)*kk(j))
        enddo
     else
        do i = 0,n-1
           idft_j(i+1) = cos(w*real(i)*kk(j)) 
        enddo
     endif
     ! DFT
     dft_j(1) = 0.5
     do i = 0,(n-1)/2-1
        dft_j(2*i+1+1) = sin(w*real(j-1)*kk(2*i+1+1))
        dft_j(2*i+2+1) = cos(w*real(j-1)*kk(2*i+2+1))
     enddo
  endif
  dft_j = (2.0d0/real(n))*dft_j
  !
  ikk = real(kk(j))
  !> Zero out derivative for oddball frequency (N/2) for N even if j = n
  if (j.eq.n.and.mod(n,2).eq.0) then
     ikk = 0.0d0
  endif
  !
  return
end subroutine dftReal
