module temporalCommunication
  !
  use mpi
  implicit none
  !
  integer :: cartComm          !< spatial communicator
  integer :: timeComm          !< temporal communicator
  integer :: myid              !< myid in a global (space and time) sense
  integer :: myid_spatial      !< myid in local spatial communicator
  integer :: myid_temporal     !< myid in local temporal communicator
  integer :: numprocs          !< total number of processors
  integer :: numprocs_spatial  !< number of spatial procs for each time instance
  integer :: numprocs_temporal !< number of procs in spatial dimension
  integer :: numprocs_active   !< number of active processors (may be less than numprocs)
  integer :: globalGroup       !< MPI group containing all ranks
  integer :: cartElem          !< number of elements in spatial communicator
  integer :: timeElem          !< number of ranks in temporal communicator
  integer :: ierr
  !
  integer :: n1,n2
  integer :: myid_1,myid_2
  integer :: elem1,elem2
  integer :: group1,group2
  integer :: comm1,comm2
  integer :: topologySpaceTime
  integer, allocatable :: ranks1(:),ranks2(:)
  !
contains
  !
  subroutine initSpaceTimeComm(ninstances,topologySpaceTime)
    !
    implicit none
    !
    integer, intent(in) :: ninstances,topologySpaceTime
    !
    !> Determine Global Information (group,id,number of procs)
    !
    call mpi_comm_group(mpi_comm_world,globalGroup,ierr)
    call mpi_comm_rank (mpi_comm_world,myid,       ierr)
    call mpi_comm_size (mpi_comm_world,numprocs,   ierr)
    !
    !> Determine number of procs in spatial and temporal dimensions
    !> and total number of active procs 
    !
    numprocs_temporal = ninstances
    numprocs_spatial  = floor(real(numprocs)/real(numprocs_temporal))
    numprocs_active   = numprocs_spatial*numprocs_temporal
    if (topologySpaceTime.eq.0) then
       n1 = numprocs_spatial
       n2 = numprocs_temporal
    elseif (topologySpaceTime.eq.1) then
       n1 = numprocs_temporal
       n2 = numprocs_spatial
    endif
    !
    !> Generate MPI Groups for Space & Time
    !
    allocate(ranks1(n1),ranks2(n2))
    !
    call getSpaceTimeRanks(myid,n1,n2,numprocs_active,ranks1,ranks2,elem1,elem2)
    call mpi_group_incl(globalGroup,elem1,ranks1(1:elem1),group1,ierr)
    call mpi_group_incl(globalGroup,elem2,ranks2(1:elem2),group2,ierr)
    !
    !> Generate MPI Communicators for Space & Time
    !
    call mpi_comm_create(mpi_comm_world,group1,comm1,ierr)
    call mpi_comm_rank(comm1,myid_1,ierr)
    call mpi_comm_size(comm1,n1,ierr) 
    
    call mpi_comm_create(mpi_comm_world,group2,comm2,ierr)
    call mpi_comm_rank(comm2,myid_2,ierr)
    call mpi_comm_size(comm2,n2,ierr)
    !                                                                                                               
    !> Assign spatial and temporal communicators                                                                    
    !                                                                                                               
    if (topologySpaceTime.eq.0) then
       myid_spatial  = myid_1
       cartComm      = comm1
       cartElem      = elem1
       
       myid_temporal = myid_2
       timeComm      = comm2
       timeElem      = elem2
    else
       myid_spatial = myid_2
       cartComm     = comm2
       cartElem     = elem2
       
       myid_temporal = myid_1
       timeComm      = comm1
       timeElem      = elem1
    endif
    !    
    !if (myid.eq.0) then
    !   write(*,*) "cartcomm for myid0: ",cartComm
    !   write(*,*) "myid_temporal for myid=0",myid_temporal
    !endif
    
    deallocate(ranks1,ranks2)
    !
    return
  end subroutine initSpaceTimeComm
  !
  subroutine getSpaceTimeRanks(myid,n1,n2,numprocs_active,ranks1,ranks2,elem1,elem2)
    !                                                                                                                                                                           
    implicit none
    !                                                                                                                                                                          
    integer,                               intent(in)  :: myid
    integer,                               intent(in)  :: n1
    integer,                               intent(in)  :: n2
    integer,                               intent(in)  :: numprocs_active
    integer, dimension(n1),                intent(out) :: ranks1
    integer, dimension(n2),                intent(out) :: ranks2
    integer,                               intent(out) :: elem1,elem2
    !                                                                                                                                                                           
    integer :: i,j,fac1,fac2
    if (myid.lt.numprocs_active) then
       !                                                                                                            
       !> Compute Node                                                                                              
       !                                                                                                            
       elem1 = n1
       fac1 = floor(real(myid)/real(n1))
       do i = 0,n1-1
          ranks1(i+1) = i + fac1*n1
       enddo
       
       elem2 = n2
       fac2 = mod(myid,n1)
       do i = 0,n2-1
          ranks2(i+1) = i*n1 + fac2
       enddo
    else
       !                                                                                                            
       !> Non-compute Node                                                                                          
       !                                                                                                            
       elem1 = 1
       ranks1(1) = myid
       
       elem2 = 1
       ranks2(1)  = myid
    endif
    !                                                                                                               
    return
  end subroutine getSpaceTimeRanks
  !
  subroutine getRankInfo(numprocs_spatial1,numprocs_temporal1,myid_spatial1,myid_temporal1,cartComm1,timeComm1)
    !
    implicit none
    !
    integer, intent(out)  :: numprocs_spatial1
    integer, intent(out)  :: numprocs_temporal1
    integer, intent(out)  :: myid_spatial1
    integer, intent(out)  :: myid_temporal1
    integer, intent(out) :: cartComm1
    integer, intent(out) :: timeComm1
    !
    numprocs_spatial1  = numprocs_spatial
    numprocs_temporal1 = numprocs_temporal
    myid_spatial1      = myid_spatial
    myid_temporal1     = myid_temporal
    cartComm1          = cartComm
    timeComm1          = timeComm
    !
    return
  end subroutine getRankInfo
!
end module temporalCommunication
