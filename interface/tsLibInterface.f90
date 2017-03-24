module libTSinterface
  use temporalCommunication, only : initSpaceTimeComm,getSpaceTimeRanks,&
                                    getRankInfo
  !
  include 'mpif.h'
  !
  integer :: numprocs,numprocs_spatial,numprocs_temporal,numprocs_active
  integer :: ninstances
  integer :: myid,myid_spatial,myid_temporal
  integer :: cartComm
  integer :: timeComm
  integer :: topologySpaceTime
  integer :: arithmetic
  !
contains
  !
  !> Initialize Data
  !
  subroutine libTS_init_data(n,tst,arithmetic_type)
    integer, intent(in)  :: n,tst
    integer, intent(in) :: arithmetic_type
    !
    ninstances        = n
    topologySpaceTime = tst
    arithmetic        = arithmetic_type
    !
    call initSpaceTimeComm(ninstances,topologySpaceTime)
    !
    call getRankInfo(numprocs_spatial,numprocs_temporal,myid_spatial,myid_temporal,cartComm,timeComm)
    !
  end subroutine libTS_init_data
  !
  !> Perform Time-Spectral Approximate Factorization Update
  !
  subroutine afts_update(q,rhs,vol,rank,ninstances,h,freq,tcomp_ts,tcomm_ts,nf,jmax,kmax,lmax)
    implicit none
    !
    integer, intent(in)    :: jmax,kmax,lmax
    real*8,  intent(in)    :: q(5,jmax,kmax,lmax)
    real*8,  intent(in)    :: vol,freq
    real*8,  intent(inout) :: rhs(5,jmax,kmax,lmax)
    integer, intent(in)    :: rank
    integer, intent(in)    :: ninstances
    real*8,  intent(in)    :: h
    real*8,  intent(inout) :: tcomp_ts,tcomm_ts
    integer, intent(in)    :: nf
    ! 
    !integer :: nq,nvar,ndof
    !
    !nq = 5
    !nvar=nq! read this in later but harcoded for now
    !ndof = nq*jmax*kmax*lmax
!    if (arithmetic.eq.0) then
!       call updateTS_real(nq,nvar,q,rhs,vol,myid_temporal,ninstances,h,freq,tcomp_ts,tcomm_ts,jmax,kmax,lmax,timeComm)
!    else
    call afts(q,rhs,vol,myid_temporal,ninstances,h,freq,tcomp_ts,tcomm_ts,jmax,kmax,lmax,timeComm,nf)
    !call afts_buffer(q,rhs,vol,rank,ninstances,h,freq,tcomp_ts,tcomm_ts,jmax,kmax,lmax,ndof)
    !call afts_uns(q,rhs,vol,rank,ninstances,h,freq,npts)
    return
  end subroutine afts_update
  !!
  !> clean up all the arrays
  !> and exit gracefully
  !!
  subroutine ts_cleanup
    implicit none
    !
    integer :: ierr
    
    !call mpi_finalize(ierr)
    !
  end subroutine ts_cleanup
  !
end module libTSinterface
