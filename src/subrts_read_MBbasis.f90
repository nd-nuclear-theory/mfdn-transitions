
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readMBgroupID_metadata( &
  fh, nclasses, nparticles, par, totMj, nsubvec, &
  indxNm, ngroupids, numstates, nblks)
  implicit none
  !
  integer, intent(in) :: fh, nclasses, nparticles, par, totMj, nsubvec
  integer, intent(out) :: indxNm, ngroupids, numstates, nblks
  !
  ! local variables
  integer, dimension(16) :: itemp
  integer :: versionnumber
  !
  versionnumber = 15003
  !
  read (fh) itemp(1:16)
  if (itemp(1) /= versionnumber) then
    print *, 'Wrong version MBgroup_info'
    call cancelall(121)
  end if
  if (itemp(2) /= nclasses) then
    print *, 'Wrong number of classes', itemp(2), nclasses
    call cancelall(122)
  end if
  if (itemp(3) /= nparticles) then
    print *, 'Wrong number of particles', itemp(3), nparticles
    call cancelall(123)
  end if
  if (itemp(5) /= par) then
    print *, 'Wrong parity', itemp(5), par
    call cancelall(124)
  end if
  if (itemp(6) /= totMj) then
    print *, 'Wrong total Mj', itemp(6), totMj
    call cancelall(125)
  end if
  indxNm = itemp(7)
  !
  if (itemp(8) /= nsubvec) then
    print *, 'Wrong number of subvectors', itemp(8), nsubvec
    call cancelall(126)
  end if
  !
  ngroupids = itemp(9)
  numstates = itemp(10)
  !
  nblks = itemp(12)
  !
  return

end subroutine readMBgroupID_metadata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
