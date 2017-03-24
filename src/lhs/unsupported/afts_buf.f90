!>
!> compute first order time-spectral derivative
!>
module temp_arrays
 integer,    save        :: not_allocated_af=1
 integer,    allocatable :: kk(:)
 complex*16, allocatable :: pij(:),ipij(:)
 complex*16, allocatable :: qctmp1(:),qctmp2(:)
 complex*16, allocatable :: rctmp1(:),rctmp2(:),s(:)
end module temp_arrays
! 
subroutine afts_buffer(q,rhs,vol,rank,n,h,w0,tcomp,tcomm,jmax,kmax,lmax,ndof)
!
use temp_arrays
!
implicit none
include 'mpif.h'

integer, intent(in)   :: jmax,kmax,lmax,ndof
real*8, intent(in)    :: q(ndof)
real*8, intent(in)    :: vol,w0
real*8, intent(inout) :: rhs(ndof)
integer, intent(in)   :: rank
integer, intent(in)   :: n
real*8, intent(in)    :: h
real*8, intent(inout) :: tcomp,tcomm
!
integer   :: i,j,k,l,c,iq,n1,buffer_size,ncomm,nsent,ii
real*8    :: pi,w,t_start,t_end
integer   :: npts,ierr
complex*16,parameter :: im=(0.0d0,1.0d0)
logical,parameter::test_case=.FALSE.
complex*16 :: num,den,dcmplxh,dcmplxw0,uhat,duhat,utmp1,utmp2,ikk
!
!ncomm = floor(ndof/buffer_size)
!nlast = ndof - ncomm*buffer_size

pi   = acos(-1.)
w    = dcmplx(2.0d0*pi/dble(n))
npts = 5*jmax*kmax*lmax
!
if (not_allocated_af) then
   allocate(rctmp1(ndof))
   allocate(rctmp2(ndof))
   allocate(qctmp1(ndof))
   allocate(qctmp2(ndof))
   allocate(s(ndof))
   allocate(pij(n))
   allocate(ipij(n))
   allocate(kk(n))
!
   j=rank+1
   
   call getDFT(pi,im,j,n,pij,ipij,kk,ikk)
   
   not_allocated_af=0
endif

! Apply Parallel DFT to q and rhs
do i=1,n
   call cpu_time(t_start)
   qctmp1=pij(i)*dcmplx(vol)*dcmplx(q)
   rctmp1=pij(i)*dcmplx(rhs)
   call cpu_time(t_end)
   tcomp = tcomp + t_end - t_start

   call cpu_time(t_start)
   
   !do ii = 0,ncomm-1
   !   snd_buff(1:ncomm) = qctmp1(ii*buffer_size+1:(ii+1)*buffer_size)
   !   call mpi_reduce(snd_buff,rcv_buff,buffer_size,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
   !   qctmp2(ii*buffer_size+1:(ii+1)*buffer_size) = rcv_buff(1:ncomm)
   !enddo

   !snd_buff(1:nlast) = qctmp1((ncomm-1)*buffer_size+1:ndof)
   !call mpi_reduce(snd_buff,rcv_buff,nlast,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
   !qctmp2((ncomm-1)*buffer_size+1:ndof) = rcv_buff(1:nlast)
   call mpi_reduce(qctmp1,qctmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
   call mpi_reduce(rctmp1,rctmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
   call cpu_time(t_end)
   tcomm = tcomm + t_end - t_start
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
call cpu_time(t_start)
do i=1,ndof
   num = ( rctmp2(i) - im*dcmplxw0*ikk*qctmp2(i))
   den =  dcmplx(1.0d0) + dcmplxh*im*dcmplxw0*ikk
   s(i) = num/den
enddo
call cpu_time(t_end)
tcomp = tcomp + t_end - t_start
!s =  (rctmp2 - im*dcmplxw0*ikk*qctmp2)/(dcmplx(1.0d0) + dcmplxh*im*dcmplxw0*ikk)

! Apply Parallel IDFT to update and place in rhs to be used for spatial implicit update
do i=1,n
   call cpu_time(t_start)
   rctmp1=ipij(i)*s
   call cpu_time(t_end)
   tcomp = tcomp + t_end - t_start

   call cpu_time(t_start)
   call mpi_reduce(rctmp1,rctmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
   call cpu_time(t_end)
   tcomm = tcomm + t_end - t_start
enddo
call cpu_time(t_start)
do i=1,ndof
   rhs(i) = dble(rctmp2(i))
enddo
call cpu_time(t_end)
tcomp = tcomp + t_end - t_start
return
end subroutine afts_buffer
