!>
!> compute first order time-spectral derivative
!>
module temp_arrays_uns
 integer,    save        :: not_allocated_af=1
 integer,    allocatable :: kk(:)
 complex*16, allocatable :: pij(:),ipij(:)
 complex*16, allocatable :: qtmp1(:),qtmp2(:)
 complex*16, allocatable :: rtmp1(:),rtmp2(:),s(:)
end module temp_arrays_uns

subroutine afts_uns(q,rhs,vol,rank,n,h,w0,npts)
!
use temp_arrays_uns
!
implicit none
include 'mpif.h'

integer, intent(in)   :: npts
real*8, intent(in)    :: q(npts)
real*8, intent(in)    :: vol,w0
real*8, intent(inout) :: rhs(npts)
integer, intent(in)   :: rank
integer, intent(in)   :: n
real*8, intent(in)    :: h
!
integer   :: i,j,k,l,c,iq,n1
real*8    :: pi,w
integer   :: ierr
complex*16,parameter :: im=(0.0d0,1.0d0)
logical,parameter::test_case=.FALSE.
complex*16 :: num,den,dcmplxh,dcmplxw0,uhat,duhat,utmp1,utmp2,ikk
!
pi   = acos(-1.)
w    = dcmplx(2.0d0*pi/dble(n))
!npts = 5*jmax*kmax*lmax
!
if (not_allocated_af) then
   allocate(rtmp1(npts))
   allocate(rtmp2(npts))
   allocate(qtmp1(npts))
   allocate(qtmp2(npts))
   allocate(s(npts))
   allocate(pij(n))
   allocate(ipij(n))
   allocate(kk(n))
!
   j=rank+1
   
   call getDFT(pi,im,j,n,pij,ipij,kk,ikk)
   
   not_allocated_af=0
endif

! Apply Parallel DFT to q
do i=1,n
   qtmp1=pij(i)*dcmplx(vol)*dcmplx(q)
   call mpi_reduce(qtmp1,qtmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
enddo
! Apply Parallel DFT to rhs
do i=1,n
   rtmp1=pij(i)*dcmplx(rhs)
   call mpi_reduce(rtmp1,rtmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
enddo

! Compute implict update in transformed space via scalar division
ikk = dcmplx(real(kk(rank+1)))
n1 = n-1
!> Zero out derivative for oddball frequency (N/2) for N even if j = n
if (rank.eq.n1.and.mod(n,2).eq.0) then
   ikk = dcmplx(0.)
endif
dcmplxh  = dcmplx(h)
dcmplxw0 = dcmplx(w0)
do j=1,npts
      num = ( rtmp2(j) - im*dcmplxw0*ikk*qtmp2(j))
      den =  dcmplx(1.0d0) + dcmplxh*im*dcmplxw0*ikk
      s(j) = num/den
enddo

! Apply Parallel IDFT to update and place in rhs to be used for spatial implicit update
do i=1,n
   rtmp1=ipij(i)*s
   call mpi_reduce(rtmp1,rtmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
enddo

do j=1,npts
   rhs(j) = dble(rtmp2(j))
enddo
return
end subroutine afts_uns
