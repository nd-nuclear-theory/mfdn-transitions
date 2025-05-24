!
module SPbasis
  implicit none
  private
  public classoffsetgen
  public setupSPbasis
  public nclasses, classoffset, nparticles, norbt, norb_p, norb_n, nspstates
  public orb_sp, mj2_sp, j2_sp, next_sp_bin
  public n_orb, l_orb, j2_orb, pr_orb, wt_orb
  public OrbitalList_t, PartitioningInfo_t
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! nparticles doesn't really belong here...?
  integer, parameter :: nclasses=2
  integer, dimension(3) :: classoffset
  integer :: nparticles, norbt, norb_p, norb_n, nspstates
  ! Single-Particle state arrays
  integer, allocatable, dimension(:) :: next_sp_bin
  integer, allocatable, dimension(:) :: orb_sp, mj2_sp, j2_sp
  ! ORBital arrays
  integer, allocatable, dimension(:) :: n_orb, l_orb ! optional
  integer, allocatable, dimension(:) :: j2_orb, pr_orb
  real(4), allocatable, dimension(:) :: wt_orb
  !
  type :: OrbitalList_t
     integer :: norb_p, norb_n
     integer, allocatable, dimension(:) :: n_orb, l_orb ! optional
     integer, allocatable, dimension(:) :: j2_orb, pr_orb
     real(4), allocatable, dimension(:) :: wt_orb
  end type OrbitalList_t
  !
  type :: PartitioningInfo_t
     integer :: size_p, size_n
     integer, allocatable, dimension(:) :: partitions_p, partitions_n
  end type PartitioningInfo_t
  !
contains
  !
  ! computes class offsets
  pure function classoffsetgen(num_orb_p, num_orb_n, j2_arr) result(clsoffset)
    implicit none
    integer, intent(in) :: num_orb_p, num_orb_n
    integer, allocatable, dimension(:), intent(in) :: j2_arr
    integer, dimension(3) :: clsoffset
    integer :: i, nstates

    clsoffset(1) = 0
    nstates = 0
    do i = 1, num_orb_p
       nstates = nstates + j2_arr(i) + 1
    enddo
    clsoffset(2) = nstates
    do i = num_orb_p + 1, num_orb_p + num_orb_n
       nstates = nstates + j2_arr(i) + 1
    enddo
    clsoffset(3) = nstates
    return
  end function classoffsetgen
  !
  ! generates (n,l,j,m) Single-Particle arrays with m from -j to + j
  subroutine SPstatesgen
    ! local variables
    integer :: i, j2, m2, currentstate
    !
    classoffset = classoffsetgen(norb_p, norb_n, j2_orb)
    nspstates = classoffset(3)
    !
    allocate(j2_sp(nspstates))
    allocate(mj2_sp(nspstates))
    allocate(orb_sp(nspstates))
    allocate(next_sp_bin(nspstates))
    !
    currentstate = 0
    do i = 1, norbt
       j2 = j2_orb(i)
       do m2 = -j2, j2, 2
          currentstate = currentstate + 1
          j2_sp(currentstate) = j2
          mj2_sp(currentstate) = m2
          orb_sp(currentstate) = i
       enddo
    enddo
    !
    return
  end subroutine SPstatesgen
  !
  ! read partitioning file, and perform basic consistency check
  subroutine init_next_sp_bin_array(part_info)
    type(PartitioningInfo_t), intent(in) :: part_info
    ! local variables
    integer, dimension(:), allocatable :: partition
    integer :: partitionsize, i, k, next
    !
    partitionsize = part_info%size_p + part_info%size_n + 1
    allocate(partition(partitionsize))
    do i = 1, part_info%size_p
      partition(i) = part_info%partitions_p(i)
    end do
    do i = 1, part_info%size_n
      partition(i+part_info%size_p) = part_info%partitions_n(i) + classoffset(2)
    end do
    partition(partitionsize) = nspstates+1
    !
    ! minimal check partitioning for consistency
    k = 1
    do i = 1, partitionsize-1
       if (partition(i+1) .le. partition(i)) then
          print *, "ERROR in input partitioning"
          print *, partition(i), partition(i+1)
          print *, "should be monotonic"
          call cancelall(211)
       endif
    enddo
    if (partition(partitionsize) .ne. nspstates+1) then
       print *, "ERROR in input partitioning"
       print *, partition(partitionsize)
       print *, "should be number of SP states is", nspstates
       call cancelall(215)
    endif
    !
    next = 2
    next_sp_bin(1) = partition(next)
    do i = 2, nspstates
       next_sp_bin(i) = next_sp_bin(i-1)
       if (i .eq. next_sp_bin(i)) then
          next = next + 1
          next_sp_bin(i) = partition(next)
       endif
    enddo
    !
    deallocate(partition)
    return
  end subroutine init_next_sp_bin_array
  !
  ! read Single-Particle basis and corresponding partitioning
  subroutine setupSPbasis(orbital_list, part_info)
    type(OrbitalList_t), intent(in) :: orbital_list
    type(PartitioningInfo_t), intent(in) :: part_info
    ! local variables
    integer :: j

    norb_p = orbital_list%norb_p
    norb_n = orbital_list%norb_n
    norbt = norb_p + norb_n
    !
    allocate(j2_orb(norbt), pr_orb(norbt))
    allocate(n_orb(norbt), l_orb(norbt))
    allocate(wt_orb(norbt))
    do j = 1, norbt
      n_orb(j)  = orbital_list%n_orb(j)
      l_orb(j)  = orbital_list%l_orb(j)
      j2_orb(j)  = orbital_list%j2_orb(j)
      pr_orb(j)  = orbital_list%pr_orb(j)
      wt_orb(j)  = orbital_list%wt_orb(j)
    enddo
    !
    call SPstatesgen
    call init_next_sp_bin_array(part_info)
    !
    return
  end subroutine setupSPbasis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module SPbasis
