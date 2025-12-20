!
!      subroutine abbrstates
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine abbrstates(n1, n2, pntr1, pntr2, nstates, mbstates, nabbr, abbr1, abbr2)
  implicit none
  !
  integer, intent(in) :: n1, n2, pntr1, pntr2, nstates, nabbr
  integer(kind=2), dimension(n1 + n2, nstates), intent(in) :: mbstates
  integer(kind=8), dimension(nabbr), intent(out) :: abbr1, abbr2
  !
  ! local variables
  integer :: i, j, itest
  integer(kind=8), parameter :: allones = -1
  !
  !     set up abbreviated MBstates:
  !     bitrepresentation for the first 64 proton SPstates (n1)
  !     and first 64 neutron SPstates (n2) in two integer(8) arrays
  !
  do j = 1, nstates
    abbr1(j) = 0
    do i = 1, n1
      itest = mbstates(i, j) - pntr1
      if (itest .lt. 64) then
        abbr1(j) = IBSet(abbr1(j), itest)
      end if
    end do
  end do
  if (nabbr .gt. nstates) abbr1(nstates + 1:nabbr) = allones
  !
  do j = 1, nstates
    abbr2(j) = 0
    do i = n1 + 1, n1 + n2
      itest = mbstates(i, j) - pntr2
      if (itest .lt. 64) then
        abbr2(j) = IBSet(abbr2(j), itest)
      end if
    end do
  end do
  if (nabbr .gt. nstates) abbr2(nstates + 1:nabbr) = allones
  !
  return
end subroutine abbrstates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
