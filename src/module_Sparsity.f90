!
module Sparsity
  implicit none
  private
  public generateSparsity
  public tilediff, tileind, ntiles
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  integer(kind=4) :: ntiles
  integer(kind=4), dimension(:), allocatable :: tilediff, tileind
  !
contains
  !
  subroutine generateSparsity(nparticles, oprank, nspstates, n2offset, &
                              n1col, n2col, ncolgroupids, colgroupidlist, &
                              n1row, n2row, nrowgroupids, rowgroupidlist, &
                              coltileptr)
    !
    integer, intent(in) :: nparticles, oprank, nspstates, n2offset
    integer, intent(in) :: n1col, n2col, n1row, n2row, ncolgroupids, nrowgroupids
    integer(kind=2), dimension(nparticles, ncolgroupids), intent(in) :: colgroupidlist
    integer(kind=2), dimension(nparticles, nrowgroupids), intent(in) :: rowgroupidlist
    integer, dimension(ncolgroupids + 1), intent(out) :: coltileptr
    !
    ! local variables
    integer, dimension(0:oprank + 1) :: ndiffcounter
    integer :: statesize
    !
    integer(kind=8), dimension(ncolgroupids) :: colabbr1, colabbr2
    integer(kind=8), dimension(nrowgroupids) :: rowabbr1, rowabbr2
    !
    statesize = nparticles + oprank
    !
    call abbrstates(n1col, n2col, 0, n2offset, ncolgroupids, colgroupidlist, &
                    ncolgroupids, colabbr1, colabbr2)
    !
    call abbrstates(n1row, n2row, 0, n2offset, nrowgroupids, rowgroupidlist, &
                    nrowgroupids, rowabbr1, rowabbr2)
    !
    call countnonzerotilesoffdiag(statesize, nparticles, oprank, nspstates, &
                                  ncolgroupids, colgroupidlist, colabbr1, colabbr2, &
                                  nrowgroupids, rowgroupidlist, rowabbr1, rowabbr2, &
                                  coltileptr, ntiles)
    !
    ! Determine location and size of nonzero tiles
    allocate (tilediff(ntiles + 1))
    tilediff = 0
    allocate (tileind(ntiles + 1))
    tileind = 0
    ndiffcounter = 0
    !
    call matrixnonzerotilesoffdiag(statesize, nparticles, oprank, nspstates, &
                                   ncolgroupids, colgroupidlist, colabbr1, colabbr2, &
                                   nrowgroupids, rowgroupidlist, rowabbr1, rowabbr2, &
                                   coltileptr, ntiles, tileind, tilediff, ndiffcounter)
    !
    return
  end subroutine generateSparsity
  !
  !
  ! count number of nonzero tiles and set coltileptr
  subroutine countnonzerotilesoffdiag(statesize, nparticles, oprank, nspstates, &
                                      ncolgroupids, colgroupidlist, colabbr1, colabbr2, &
                                      nrowgroupids, rowgroupidlist, rowabbr1, rowabbr2, &
                                      coltileptr, ntiles)
    !
    integer, intent(in) :: statesize, nparticles, oprank, nspstates
    integer, intent(in) :: ncolgroupids, nrowgroupids
    integer(kind=2), dimension(nparticles, ncolgroupids), intent(in) :: colgroupidlist
    integer(kind=2), dimension(nparticles, nrowgroupids), intent(in) :: rowgroupidlist
    !
    integer(kind=8), dimension(ncolgroupids), intent(in) :: colabbr1, colabbr2
    integer(kind=8), dimension(nrowgroupids), intent(in) :: rowabbr1, rowabbr2
    !
    integer, intent(out) :: ntiles
    integer, dimension(ncolgroupids + 1), intent(out) :: coltileptr
    !
    ! local variables
    integer :: i, j, jmin, jj, nnz, nnonzero
    integer, dimension(ncolgroupids) :: colnonzero
    integer, dimension(statesize) :: colstate, rowstate
    integer(kind=8) :: ii1, ii2, ii3, ii4, XOR1, XOR2
    integer :: ndiffs
    !integer(kind=8), dimension(:), allocatable :: XOR1, XOR2
    !integer, dimension(:), allocatable :: rowdiflist, ndiflist
    !
    !$omp parallel default(shared) &
    !$omp          private(i, j, jmin, jj, nnz, nnonzero, &
    !$omp                  colstate, rowstate, ndiffs,    &
    !$omp                  ii1, ii2, ii3, ii4, XOR1, XOR2)
    !
    do i = nparticles + 1, statesize
      colstate(i) = nspstates + i
      rowstate(i) = nspstates + i
    end do
    !
    !$omp do schedule(dynamic)
    do i = 1, ncolgroupids
      nnz = 0
      nnonzero = 0
      colstate(1:nparticles) = colgroupidlist(1:nparticles, i)
      ii1 = colabbr1(i)
      ii2 = colabbr2(i)
      ! loop over row states
      do j = 1, nrowgroupids
        ii3 = rowabbr1(j)
        XOR1 = ieor(ii1, ii3)
        ii4 = rowabbr2(j)
        XOR2 = ieor(ii2, ii4)
        ndiffs = popcnt(XOR1) + popcnt(XOR2)
        if (ndiffs .le. 2*oprank) then
          ! detailed comparison
          rowstate(1:nparticles) = rowgroupidlist(1:nparticles, j)
          call groupIDnonzero(statesize, colstate, rowstate, ndiffs)
          ! add to counter
          if (ndiffs .le. oprank) nnonzero = nnonzero + 1
        end if
        !
      end do
      ! done
      colnonzero(i) = nnonzero
    end do
    !$omp end do
    !
    !$omp end parallel
    !
    coltileptr(1) = 1
    do i = 1, ncolgroupids
      coltileptr(i + 1) = coltileptr(i) + colnonzero(i)
    end do
    ! should be equal to nnonzero
    ntiles = coltileptr(ncolgroupids + 1) - 1
    !
    return
  end subroutine countnonzerotilesoffdiag
  !
  !
  ! construct number of nonzero tiles
  subroutine matrixnonzerotilesoffdiag(statesize, nparticles, oprank, nspstates, &
                                       ncolgroupids, colgroupidlist, colabbr1, colabbr2, &
                                       nrowgroupids, rowgroupidlist, rowabbr1, rowabbr2, &
                                       coltileptr, ntiles, tileind, tilediff, diffcounter)
    !
    integer, intent(in) :: statesize, nparticles, oprank, nspstates
    integer, intent(in) :: ncolgroupids, nrowgroupids
    integer(kind=2), dimension(nparticles, ncolgroupids), intent(in) :: colgroupidlist
    integer(kind=2), dimension(nparticles, nrowgroupids), intent(in) :: rowgroupidlist
    !
    integer(kind=8), dimension(ncolgroupids), intent(in) :: colabbr1, colabbr2
    integer(kind=8), dimension(nrowgroupids), intent(in) :: rowabbr1, rowabbr2
    !
    integer, intent(in) :: ntiles
    integer, dimension(ncolgroupids + 1), intent(in) :: coltileptr
    integer, dimension(ntiles + 1), intent(out) :: tileind, tilediff
    integer, dimension(0:oprank), intent(out) :: diffcounter
    !
    ! local variables
    integer :: i, j, jmin, jj, nnz, nnonzero, ndiffs
    integer, dimension(ncolgroupids) :: colnonzero
    integer, dimension(statesize) :: colstate, rowstate
    integer(kind=8) :: ii1, ii2, ii3, ii4, XOR1, XOR2
    !
    diffcounter(0:oprank) = 0
    !
    !$omp parallel default(shared) &
    !$omp          private(i, j, jmin, jj, nnz, nnonzero,  &
    !$omp                  colstate, rowstate, ndiffs,     &
    !$omp                  ii1, ii2, ii3, ii4, XOR1, XOR2) &
    !$omp          reduction(+:diffcounter)
    !
    do i = nparticles + 1, statesize
      colstate(i) = nspstates + i
      rowstate(i) = nspstates + i
    end do
    !
    !$omp do schedule(dynamic)
    do i = 1, ncolgroupids
      nnz = 0
      nnonzero = coltileptr(i)
      colstate(1:nparticles) = colgroupidlist(1:nparticles, i)
      ii1 = colabbr1(i)
      ii2 = colabbr2(i)
      ! loop over row states
      do j = 1, nrowgroupids
        ii3 = rowabbr1(j)
        XOR1 = ieor(ii1, ii3)
        ii4 = rowabbr2(j)
        XOR2 = ieor(ii2, ii4)
        ndiffs = popcnt(XOR1) + popcnt(XOR2)
        if (ndiffs .le. 2*oprank) then
          ! detailed comparison
          rowstate(1:nparticles) = rowgroupidlist(1:nparticles, j)
          call groupIDnonzero(statesize, colstate, rowstate, ndiffs)
          if (ndiffs .le. oprank) then
            tileind(nnonzero) = j
            tilediff(nnonzero) = ndiffs
            diffcounter(ndiffs) = diffcounter(ndiffs) + 1
            nnonzero = nnonzero + 1
          end if
        end if
        !
      end do
      ! done
      colnonzero(i) = nnonzero
    end do
    !$omp end do
    !
    !$omp end parallel
    !
    return
  end subroutine matrixnonzerotilesoffdiag
  !
  !
  ! might also be used elsewhere -- move to a different module or keep as subroutine?
  subroutine groupIDnonzero(statesize, colstate, rowstate, ndiffs)
    !
    integer, intent(in) :: statesize
    integer, dimension(statesize), intent(in) :: colstate, rowstate
    integer, intent(out) :: ndiffs
    !
    ! local variables
    integer :: i, irow, icol, nrowdiffs, ncoldiffs
    !
    irow = 1
    icol = 1
    nrowdiffs = 0
    ncoldiffs = 0
    !
    do i = 1, statesize
      if (rowstate(irow) .eq. colstate(icol)) then
        irow = irow + 1
        icol = icol + 1
      elseif (rowstate(irow) .gt. colstate(icol)) then
        icol = icol + 1
        ncoldiffs = ncoldiffs + 1
      else
        irow = irow + 1
        nrowdiffs = nrowdiffs + 1
      end if
    end do
    !
    ndiffs = max(ncoldiffs, nrowdiffs)
    !
    return
  end subroutine groupIDnonzero
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module Sparsity
