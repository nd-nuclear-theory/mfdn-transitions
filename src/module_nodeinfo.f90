
module nodeinfo
  use MPI
  use omp_lib
  implicit none
  ! MPI variables
  integer :: icomm, myrank, nprocs, root, ierr
  ! OMP number of threads
  integer :: numthreads
  ! run control parameter
  logical :: checkpoint = .false.
  ! user-defined OMP enabled MPI_reduce operations
  integer :: My_MPI_Sum
  !
contains
  !
  subroutine Setup_Parallel
    integer :: multithreaded
    !$omp single
    ! Initialize MPI
    root = 0
    icomm = MPI_COMM_WORLD
    call MPI_Init_Thread(MPI_THREAD_FUNNELED, multithreaded, ierr)
    call MPI_Comm_Size(icomm, nprocs, ierr)
    call MPI_Comm_Rank(icomm, myrank, ierr)
    !
    ! Set/check number of threads
#ifdef _OPENMP
    numthreads = omp_get_max_threads()
#else
    numthreads = 1
#endif
    !$omp end single
    return
  end subroutine Setup_Parallel

  subroutine SplitProcsSeq(ncolprocs, nrowprocs, colstart, rowstart, nblocks)
    integer, intent(in) :: ncolprocs, nrowprocs
    integer, intent(out) :: colstart, rowstart, nblocks
    ! local variables
    integer :: ntotblocks, nextra, startblock
    !
    ! Determine assignments of work for each rank
    !
    ! Sequentially assign submatrices to ranks.  Each rank will have
    ! a start position and wrap around to the beginning of the next row
    ! when it reaches the end of the row.
    !
    ! Best performance anticipated for
    !   nprocs =< colprocs  with  ncolprocs / nprocs  integer
    !   nprocs >= colprocs  with  nprocs / ncolprocs  integer
    !
    ntotblocks = ncolprocs*nrowprocs
    nblocks = ntotblocks/nprocs
    nextra = mod(ntotblocks, nprocs)
    !
    if (myrank .lt. nextra) then
      nblocks = nblocks + 1
      startblock = myrank*nblocks
    else
      startblock = myrank*nblocks + nextra
    end if
    !
    colstart = startblock/nrowprocs
    rowstart = mod(startblock, nrowprocs)
    !
    return
  end subroutine SplitProcsSeq

  subroutine SplitProcs(ncolprocs, nrowprocs, colstart, rowstart, nblocks)
    integer, intent(in) :: ncolprocs, nrowprocs
    integer, intent(out) :: colstart, rowstart, nblocks
    ! limits to be distributed to procs
    integer :: i, j, k, indx, itmp, ntotblocks
    integer, dimension(:), allocatable :: colstart_assn, rowstart_assn
    integer, dimension(:), allocatable :: nblocks_assn
    !
    ! Determine assignments of work for each rank
    !
    !   Round-robin assign tiles to ranks. Each rank will have a start position
    !   and wrap around to the beginning of the next row when it reaches the
    !   end of the row.
    !
    if (myrank .eq. root) then
      allocate (colstart_assn(nprocs))
      allocate (rowstart_assn(nprocs))
      allocate (nblocks_assn(nprocs))
      ntotblocks = ncolprocs*nrowprocs
      ! if (ntotblocks .gt. nprocs) then
      !   cancelall(300)
      ! endif
      ! round-robin assignment of tiles
      do i = 0, ntotblocks - 1
        indx = mod(i, nprocs) + 1
        nblocks_assn(indx) = nblocks_assn(indx) + 1
      end do
      ! find starting tile locations
      itmp = 0
      do indx = 1, nprocs
        colstart_assn(indx) = itmp/nrowprocs
        rowstart_assn(indx) = mod(itmp, nrowprocs)
        itmp = itmp + nblocks_assn(indx)
      end do
    end if
    ! Scatter assignments to ranks
    call MPI_Scatter(colstart_assn, 1, MPI_INTEGER, colstart, 1, MPI_INTEGER, &
                     root, icomm, ierr)
    call MPI_Scatter(rowstart_assn, 1, MPI_INTEGER, rowstart, 1, MPI_INTEGER, &
                     root, icomm, ierr)
    call MPI_Scatter(nblocks_assn, 1, MPI_INTEGER, nblocks, 1, MPI_INTEGER, &
                     root, icomm, ierr)
    ! deallocate temporary arrays on root
    if (myrank .eq. root) then
      deallocate (colstart_assn)
      deallocate (rowstart_assn)
      deallocate (nblocks_assn)
    end if
    return
  end subroutine SplitProcs

  subroutine MFDn_Sum_OMP(inarray, outarray, length, datatype)
    implicit none
    real(kind=4), dimension(length), intent(in) :: inarray
    real(kind=4), dimension(length), intent(inout) :: outarray
    integer, intent(in) :: length, datatype
    integer :: i
    !$omp parallel do default(shared) private(i)
    do i = 1, length
      outarray(i) = outarray(i) + inarray(i)
    end do
    !$omp end parallel do
    return
  end subroutine MFDn_Sum_OMP
  !
  ! subroutine MFDn_Sum_saxpy(inarray, outarray, length, datatype)
  !   implicit none
  !   real(kind=4), dimension(length), intent(in) :: inarray
  !   real(kind=4), dimension(length), intent(inout) :: outarray
  !   integer, intent(in) :: length, datatype
  !   !
  !   call saxpy(length, 1.0, inarray, 1, outarray, 1)
  !   !
  !   return
  ! end subroutine MFDn_Sum_saxpy
  !
end module nodeinfo
