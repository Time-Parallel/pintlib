!>
!> compute first order time-spectral derivative
!>
module temp_arrays
 integer,    save        :: not_allocated_af=1
 integer,    allocatable :: kk(:)
 complex*16, allocatable :: pij(:),ipij(:)
 complex*16, allocatable :: qctmp1(:,:,:,:),qctmp2(:,:,:,:)
 complex*16, allocatable :: rctmp1(:,:,:,:),rctmp2(:,:,:,:),s(:,:,:,:)
 complex*16, allocatable :: dft_snd(:,:,:,:,:),dft_rcv(:,:,:,:,:)
end module temp_arrays
! 
subroutine afts(q,rhs,vol,rank,n,h,w0,tcomp,tcomm,jmax,kmax,lmax,timeComm,nf)
!
use temp_arrays
!
implicit none
include 'mpif.h'

integer, intent(in)    :: jmax,kmax,lmax
real*8,  intent(in)    :: q(5,jmax,kmax,lmax)
real*8,  intent(in)    :: vol,w0
real*8,  intent(inout) :: rhs(5,jmax,kmax,lmax)
integer, intent(in)    :: rank
integer, intent(in)    :: n
real*8,  intent(in)    :: h
real*8,  intent(inout) :: tcomp,tcomm
integer, intent(in)    :: timeComm,nf
!
integer               :: i,j,k,l,c,iq,n1,js,je,ks,ke,ls,le
real*8                :: pi,w,t_start,t_end
integer               :: npts,ierr
complex*16, parameter :: one=(1.0d0,0.0d0),im=(0.0d0,1.0d0)
logical,    parameter :: test_case=.FALSE.
complex*16            :: num,den,iden,coef,dcmplxh,dcmplxw0,uhat,duhat,utmp1,utmp2,ikk
!
js = nf+1
je = jmax-nf
ks = nf+1
ke = kmax-nf
ls = nf+1
le = lmax-nf
!
pi   = acos(-1.)
w    = dcmplx(2.0d0*pi/dble(n))
npts = 5*jmax*kmax*lmax
!
if (not_allocated_af.eq.1) then
   allocate(rctmp1(5,jmax,kmax,lmax))
   allocate(rctmp2(5,jmax,kmax,lmax))
   allocate(qctmp1(5,jmax,kmax,lmax))
   allocate(qctmp2(5,jmax,kmax,lmax))
   allocate(s(5,jmax,kmax,lmax))
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
   call mpi_reduce(qctmp1,qctmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,timeComm,ierr)
   call mpi_reduce(rctmp1,rctmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,timeComm,ierr)
   call cpu_time(t_end)
   tcomm = tcomm + t_end - t_start
enddo

call cpu_time(t_start)
! Compute implict update in transformed space via scalar division
ikk = dcmplx(real(kk(rank+1)))
n1  = n-1

!> Zero out derivative for oddball frequency (N/2) for N even if j = n
if (rank.eq.n1.and.mod(n,2).eq.0) then
   ikk = dcmplx(0.)
endif
dcmplxh  = dcmplx(h)
dcmplxw0 = dcmplx(w0)
den = dcmplx(1.0d0) + dcmplxh*im*dcmplxw0*ikk
iden = one/den
coef = im*dcmplxw0*ikk

! Perform update
s = iden*(rctmp2-coef*qctmp2)

call cpu_time(t_end)
tcomp = tcomp + t_end - t_start

! Apply Parallel IDFT to update and place in rhs to be used for spatial implicit update
do i=1,n
   call cpu_time(t_start)
   rctmp1=ipij(i)*s
   call cpu_time(t_end)
   tcomp = tcomp + t_end - t_start

   call cpu_time(t_start)
   call mpi_reduce(rctmp1,rctmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,timeComm,ierr)
   call cpu_time(t_end)
   tcomm = tcomm + t_end - t_start
enddo

! Store in residual array
call cpu_time(t_start)
do l=ls,le
   do k=ks,ke
      do j=js,je
         rhs(:,j,k,l) = dble(rctmp2(:,j,k,l))
      enddo
   enddo
enddo
call cpu_time(t_end)
tcomp=tcomp+t_end-t_start
return
end subroutine afts

subroutine afts_double(q,rhs,vol,rank,n,h,w0,tcomp,tcomm,jmax,kmax,lmax,timeComm)
!
use temp_arrays
!
implicit none
include 'mpif.h'

integer, intent(in)   :: jmax,kmax,lmax
real*8, intent(in)    :: q(5,jmax,kmax,lmax)
real*8, intent(in)    :: vol,w0
real*8, intent(inout) :: rhs(5,jmax,kmax,lmax)
integer, intent(in)   :: rank
integer, intent(in)   :: n
real*8, intent(in)    :: h
real*8, intent(inout) :: tcomp,tcomm
integer, intent(in) :: timeComm
!
integer   :: i,j,k,l,c,iq,n1,js,je,ks,ke,ls,le
real*8    :: pi,w,t_start,t_end
integer   :: npts,ierr
complex*16,parameter :: im=(0.0d0,1.0d0)
logical,parameter::test_case=.FALSE.
complex*16 :: num,den,dcmplxh,dcmplxw0,uhat,duhat,utmp1,utmp2,ikk,pvol
!
js=2
je=jmax-1
ks=2
ke=kmax-1
ls=2
le=lmax-1
!
pi   = acos(-1.)
w    = dcmplx(2.0d0*pi/dble(n))
npts = 5*jmax*kmax*lmax
!
if (not_allocated_af.eq.1) then
   allocate(dft_snd(5,jmax,kmax,lmax,2))
   allocate(dft_rcv(5,jmax,kmax,lmax,2))
   allocate(qctmp1(5,jmax,kmax,lmax))
   allocate(qctmp2(5,jmax,kmax,lmax))
   allocate(rctmp2(5,jmax,kmax,lmax))
   allocate(s(5,jmax,kmax,lmax))
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
   dft_snd(:,:,:,:,1)=pij(i)*dcmplx(vol)*dcmplx(q)
   dft_snd(:,:,:,:,2)=pij(i)*dcmplx(rhs)
   call cpu_time(t_end)
   tcomp = tcomp + t_end - t_start

   call cpu_time(t_start)
   call mpi_reduce(dft_snd,dft_rcv,2*npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,timeComm,ierr)
   qctmp2 = dft_rcv(:,:,:,:,1)
   rctmp2 = dft_rcv(:,:,:,:,2)
   !call mpi_reduce(rctmp1,rctmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,mpi_comm_world,ierr)
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
do l=ls,le
   do k=ks,ke
      do j=js,je
         do iq=1,5
            num = ( rctmp2(iq,j,k,l) - im*dcmplxw0*ikk*qctmp2(iq,j,k,l))
            den =  dcmplx(1.0d0) + dcmplxh*im*dcmplxw0*ikk
            s(iq,j,k,l) = num/den
         enddo
      enddo
   enddo
enddo
call cpu_time(t_end)
tcomp = tcomp + t_end - t_start
!s =  (rctmp2 - im*dcmplxw0*ikk*qctmp2)/(dcmplx(1.0d0) + dcmplxh*im*dcmplxw0*ikk)

! Apply Parallel IDFT to update and place in rhs to be used for spatial implicit update
do i=1,n
   call cpu_time(t_start)
   qctmp1=ipij(i)*s
   call cpu_time(t_end)
   tcomp = tcomp + t_end - t_start

   call cpu_time(t_start)
   call mpi_reduce(qctmp1,qctmp2,npts,MPI_DOUBLE_COMPLEX,MPI_SUM,i-1,timeComm,ierr)
   call cpu_time(t_end)
   tcomm = tcomm + t_end - t_start
enddo

do l=ls,le
   do k=ks,ke
      do j=js,je
            rhs(:,j,k,l) = dble(qctmp2(:,j,k,l))
      enddo
   enddo
enddo
return
end subroutine afts_double
