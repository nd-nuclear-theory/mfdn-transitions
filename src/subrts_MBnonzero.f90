
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MBnonzeroLocation(statesize, rowstate, colstate, rowdifloc, coldifloc, ndiffs)
  implicit none
  !
  ! Compare rowstate and colstate and return the locations wich differ
  !
  integer, intent(in) :: statesize
  integer(kind=2), dimension(statesize), intent(in) :: rowstate
  integer(kind=2), dimension(statesize), intent(in) :: colstate
  integer, dimension(statesize), intent(inout) :: rowdifloc
  integer, dimension(statesize), intent(inout) :: coldifloc
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
        ncoldiffs = ncoldiffs + 1
        coldifloc(ncoldiffs) = icol
        icol = icol + 1
     else
        nrowdiffs = nrowdiffs + 1
        rowdifloc(nrowdiffs) = irow
        irow = irow + 1
     endif
  enddo
  !
  ndiffs = max(ncoldiffs, nrowdiffs)
  !
  return
end subroutine MBnonzeroLocation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
