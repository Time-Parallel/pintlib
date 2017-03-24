!>
!> Update Temporal ADI component using the Real-valued DFT
!>
module temp_arrays_real
 integer, save         :: not_allocated_af=1
 integer, allocatable  :: kk(:),sk(:)
 real*8,  allocatable  :: dft_j(:),idft_j(:)
 real*8,  allocatable  :: sbuf(:,:),rbuf(:,:)
end module temp_arrays_real
! 
subroutine updateTS_real(nq,nvar,r,vol,rank,n,h,w0,tcomp,tcomm,jmax,kmax,lmax,timeComm)
  !
  use temp_arrays_real
  !
  implicit none
  include 'mpif.h'
  
  integer, intent(in)    :: jmax,kmax,lmax
  real*8,  intent(in)    :: vol,w0
  real*8,  intent(inout) :: rhs(nq*jmax*kmax*lmax)
  integer, intent(in)    :: rank
  integer, intent(in)    :: n,nq,nvar
  real*8,  intent(in)    :: h
  real*8,  intent(inout) :: tcomp,tcomm
  integer, intent(in)    :: timeComm
  !
  integer   :: i,ii,j,k,l,c,iq,n1,js,je,ks,ke,ls,le
  real*8    :: pi,w,t_start,t_end
  integer   :: npts,ierr
  logical,parameter::test_case=.FALSE.
  real*8 :: num,den,ikk,iden,wk
  integer :: nbuf1,nbuf2,ndof
  real, parameter :: zero=0.0d0,one=1.0d0,two=2.0d0
  !
  ndof = nq*jmax*kmax*lmax
  !
  nbuf2 = floor(ndof/2) !< length of subset 2
  nbuf1 = ndof-nbuf2    !< length of subset 1 (nbuf1 >= nbuf2)
  !
  pi   = acos(-1.)
  w    = two*pi/dble(n)
  npts = nq*jmax*kmax*lmax
  !
  js = 2
  je = jmax-1
  ks = 2
  ke = kmax-1
  ls = 2
  le = lmax-1
  !
  if (not_allocated_af.eq.1) then
     allocate(dft_j(n))
     allocate(idft_j(n))
     allocate(kk(n),sk(n))
     !
     allocate(sbuf(nbuf1,2),rbuf(nbuf1,2))
     !
     sbuf = zero
     rbuf = zero
     !
     j=rank+1
     !
     rhs=0.0d0
     !
     call dftReal(pi,j,n,dft_j,idft_j,kk,ikk)
     !
     do i = 1,n
        sk(i) = kk(i)*(-1)**(i)
     enddo
     !
!     write(*,*) "dft(",j,"):", dft_j
!     write(*,*) "idft(",j,"):", idft_j
!     write(*,*) "sk: ", sk
     !
     !do i = 1,n
     !   qtarget(i) = i - 1
     !enddo
     !
     !rtarget=0
     !if (mod(n,2).eq.0) then
     !   i = 0
     !   do ii = 2,n/2-1
           !                                                                                                                                                                                                                                
     !      rtarget(i + 2) = i + 2
     !      rtarget(i + 3) = i + 1
           !                                                                                                                                                                                                                                
     !      i = i + 2
     !   enddo
        !
     !else
     !   i = 0
     !   do ii = 2,(n-1)/2+1
           !
     !      rtarget(i + 2) = i + 2
     !      rtarget(i + 3) = i + 1
           !
     !      i = i + 2
     !   enddo
     !endif
     !
     not_allocated_af=0
  endif
  !if (j.eq.1) then
  !   write(*,*) "qtarget: ",qtarget
  !   write(*,*) "rtarget: ",rtarget
  !endif
  
  
  ! Apply Parallel DFT to q and rhs
  do i=1,n
     call cpu_time(t_start)
     qhat_ij=dft_j(i)*vol*q
     rhat_ij=dft_j(i)*rhs
     call cpu_time(t_end)
     tcomp = tcomp + t_end - t_start
     !
     call cpu_time(t_start)
     ! Send qhat to qtarget(i)
     call mpi_reduce(qhat_ij,qhat,npts,MPI_DOUBLE,MPI_SUM,qtarget(i),timeComm,ierr)
     ! Send rhat to rtarget(i)
     call mpi_reduce(rhat_ij,rhat,npts,MPI_DOUBLE,MPI_SUM,rtarget(i),timeComm,ierr)
     call cpu_time(t_end)
     tcomm = tcomm + t_end - t_start
  enddo
  
  ! Compute implict update in transformed space via scalar division
  ikk = real(sk(rank+1))
  n1 = n-1
  !
  !> Zero out derivative for oddball frequency (N/2) for N even if j = n
  if (rank.eq.n1.and.mod(n,2).eq.0) then
     ikk = 0.0d0
  endif
  wk = real(w*ikk)
  den = (1.0d0 + h*wk)
  iden = 1.0d0/den
  call cpu_time(t_start)
  do l=ls,le
     do k=ks,ke
        do j=js,je
           do iq=1,nvar
              s(iq,j,k,l) = ( rhat(iq,j,k,l) - wk*qhat(iq,j,k,l) )*iden
           enddo
        enddo
     enddo
  enddo
  call cpu_time(t_end)
  tcomp = tcomp + t_end - t_start
  !
  ! Apply Parallel IDFT to update and place in rhs to be used for spatial implicit update
  !
  do i=1,n
     call cpu_time(t_start)
     rhat=idft_j(i)*s
     call cpu_time(t_end)
     tcomp = tcomp + t_end - t_start
     !
     call cpu_time(t_start)
     call mpi_reduce(rhat,rhat_ij,npts,MPI_DOUBLE,MPI_SUM,qtarget(i),timeComm,ierr)
     call cpu_time(t_end)
     tcomm = tcomm + t_end - t_start
  enddo
  !
  do l=ls,le
     do k=ks,ke
        do j=js,je
           do iq=1,nvar
              rhs(iq,j,k,l) = rhat_ij(iq,j,k,l)
           enddo
        enddo
     enddo
  enddo
  !
  return
end subroutine updateTS_real
