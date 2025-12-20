!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function myphase(i)
  implicit none
  integer, intent(in) :: i
  !
  !      myphase = (-1)**i
  !      myphase = 1 - 2 * mod(i,2)
  !      if (mod(i,2).eq.0) then
  !         myphase = 1
  !      else
  !         myphase = -1
  !      endif
  !
  myphase = 1.0 - 2.0*IAND(i, 1)
  !
  return
end function myphase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical function pairwiseless(a, b, c, d)
  implicit none
  integer, intent(in) :: a, b, c, d
  !
#ifdef DeltaTz
  pairwiseless = .true.
  return
#else
  if ((a .lt. c) .or. ((a .eq. c) .and. (b .le. d))) then
    pairwiseless = .true.
  else
    pairwiseless = .false.
  end if
#endif
  !
  return
end function pairwiseless

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer(kind=4) function intfour(neight, ierror)
  implicit none
  integer(kind=8), intent(in) :: neight
  integer(kind=4), intent(in) :: ierror
  !
  !     converts integer*8 to integer*4
  !     aborts MPI and stops if neight > 2^31
  !
  intfour = neight
  if (intfour .ne. neight) then
    write (0, *) ' Fatal error in program MFDn'
    write (0, *) ' integer*8 too big in function intfour', neight
    write (0, *) ' Error code', ierror
    call cancelall(ierror)
    stop
  end if
  !
  return
end function intfour

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cancelall(ierror)
  use nodeinfo
  use MPI
  implicit none
  integer(4), intent(in) :: ierror
  !
  !     writes error (ideally unique) code and cancels all processes
  !
  print *, ' Fatal error in program MFDn Transitions'
  print *, ' MFDn Transitions error code ', ierror
  print *, ' MPI rank', myrank
  !
  print *, ' Aborting MPI'
  call MPI_Abort(MPI_COMM_WORLD, ierror, ierr)
  !
  stop
  !
end subroutine cancelall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
