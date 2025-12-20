!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     contains
!       subroutine ApplyTBop
!       subroutine ApplyTBopTileTD
!       subroutine ApplyTBopTileBU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ApplyTBop(maxnMstates, &
                     Z_col, N_col, ncolgroupids, colgroupidlist, &
                     TwoMj_col, colMstateptr, numcolstates, &
                     Z_row, N_row, nrowgroupids, rowgroupidlist, &
                     TwoMj_row, rowMstateptr, numrowstates, &
                     coltileptr, ntiles, rowind, tilediff, namps, colamp, rowamp)
  !
  use SPbasis, only: nparticles, nspstates, classoffset
  use SPbasis, only: orb_sp, mj2_sp, next_sp_bin, j2_orb, pr_orb
  use MjStates, only: MjStatesGen
#ifdef DeltaTz
  use TBME_Tz12, only: Jop, nTBMEs_max, Set_TBME_array
#else
  use TBME_Tz0, only: Jop, nTBMEs_max, Set_TBME_array
#endif
  implicit none
  integer, intent(in) :: maxnMstates, ntiles, namps
  integer, intent(in) :: Z_col, N_col, TwoMj_col, ncolgroupids, numcolstates
  integer, intent(in) :: Z_row, N_row, TwoMj_row, nrowgroupids, numrowstates
  integer(kind=2), dimension(nparticles, ncolgroupids), intent(in) :: colgroupidlist
  integer(kind=2), dimension(nparticles, nrowgroupids), intent(in) :: rowgroupidlist
  integer, dimension(ncolgroupids + 1), intent(in) :: colMstateptr
  integer, dimension(nrowgroupids + 1), intent(in) :: rowMstateptr
  integer, dimension(ncolgroupids + 1), intent(in) :: coltileptr
  integer, dimension(ntiles + 1), intent(in) :: rowind, tilediff
  real(kind=4), dimension(namps, numcolstates), intent(in) :: colamp
  real(kind=4), dimension(namps, numrowstates), intent(inout) :: rowamp
  !
  ! local variables
  logical, external :: pairwiseless
  logical :: abIDN, cdIDN, DIAG
  integer :: DeltaMj, TwoMj, offset1, offset2, statesize, nTBMEs, reorder
  integer :: i, j, k, k1, k2, ndiffs, iprev, icur, n_grp_SPstates
  integer :: ii, jj, kk
  integer :: icol, nabbc, colnumstates, jrow, nabbr, rownumstates
  integer :: orba, orbb, orbc, orbd, j2a, j2b, j2c, j2d, parab, parcd
  integer(kind=2), dimension(nparticles, maxnMstates) :: colMBstates, rowMBstates
  integer(kind=2), dimension(nparticles + 2) :: colstate, rowstate, tmpstate
  integer(kind=2), dimension(2) :: rowdiffs, coldiffs
  integer, dimension(nparticles + 2) :: rowdifloc, coldifloc
  real(kind=4), dimension(1, nTBMEs_max) :: TBMEarray
  !
  integer(kind=8), dimension(maxnMstates) :: colabbr1, colabbr2, rowabbr1, rowabbr2
  integer(kind=8), dimension(maxnMstates) :: colMBbitrep
  real(kind=4), dimension(namps, maxnMstates) :: tempamp
  integer, dimension(nparticles) :: colgrp_SPoffset
  !
  TwoMj = TwoMj_row - TwoMj_col
  DeltaMj = TwoMj/2
  offset1 = classoffset(1)
  offset2 = classoffset(2)
  statesize = nparticles + 2
  !
  !$omp parallel default(shared)                                   &
  !$omp          private(i, j, k, k1, k2,                          &
  !$omp                  ii, jj, kk,                               &
  !$omp                  icol, nabbc, colabbr1, colabbr2,          &
  !$omp                  colnumstates, colMBstates, colstate,      &
  !$omp                  colMBbitrep, colgrp_SPoffset,             &
  !$omp                  iprev, icur, n_grp_SPstates,              &
  !$omp                  jrow, nabbr, rowabbr1, rowabbr2,          &
  !$omp                  rownumstates, rowMBstates, rowstate,      &
  !$omp                  ndiffs, rowdifloc, coldifloc, tmpstate,   &
  !$omp                  orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d,     &
  !$omp                  parab, parcd, coldiffs, rowdiffs,         &
  !$omp                  reorder, abIDN, cdIDN, DIAG, nTBMEs, TBMEarray) &
  !$omp          private(tempamp)
  !
  colstate(nparticles + 1) = nspstates + 1
  colstate(nparticles + 2) = nspstates + 2
  rowstate(nparticles + 1) = nspstates + 1
  rowstate(nparticles + 2) = nspstates + 2
  tmpstate(nparticles + 1) = nspstates + 1
  tmpstate(nparticles + 2) = nspstates + 2
  !
  !$omp do schedule(dynamic)
  do i = 1, ncolgroupids
    icol = colMstateptr(i)
    colstate(1:nparticles) = colgroupidlist(:, i)
    !
    ! generate group of mbstates in M-scheme for colindex i
    colnumstates = maxnMstates
    call mjstatesgen(nparticles, nspstates, mj2_sp, TwoMj_col, &
                     next_sp_bin, colstate, colnumstates, colMBstates)
    !
    nabbc = colnumstates
    call abbrstates(Z_col, N_col, offset1, offset2, &
                    colnumstates, colMBstates, nabbc, colabbr1, colabbr2)
    !
    ! 'active' S.P. states in this colgroup
    iprev = colstate(1)
    colgrp_SPoffset(1) = iprev - 1
    n_grp_SPstates = next_sp_bin(iprev) - iprev
    do k = 2, nparticles
      icur = colstate(k)
      if (next_sp_bin(icur) .eq. next_sp_bin(iprev)) then
        ! icur and iprev are in the same SPbin (same partition)
        colgrp_SPoffset(k) = colgrp_SPoffset(k - 1)
      else
        ! icur is first (?) SPstate in new SPbin
        colgrp_SPoffset(k) = icur - n_grp_SPstates - 1
        n_grp_SPstates = n_grp_SPstates + next_sp_bin(icur) - icur
        iprev = icur
      end if
    end do
    !
    ! bit representation of colMBstates in terms of 'active' S.P. states
    do j = 1, colnumstates
      colMBbitrep(j) = 0
      do k = 1, nparticles
        icur = colMBstates(k, j) - colgrp_SPoffset(k)
        colMBbitrep(j) = IBset(colMBbitrep(j), icur)
      end do
    end do
    !
    ! loop over all interacting rowgroups in this colunm group
    k1 = coltileptr(i)
    k2 = coltileptr(i + 1) - 1
    do k = k1, k2
      if (tilediff(k) .gt. 2) cycle ! hard-wired for 2-body observables
      !
      j = rowind(k)
      jrow = rowMstateptr(j)
      rowstate(1:nparticles) = rowgroupidlist(:, j)
      !
      ! generate group of mbstates in M-scheme for rowindex
      rownumstates = maxnMstates
      call mjstatesgen(nparticles, nspstates, mj2_sp, TwoMj_row, &
                       next_sp_bin, rowstate, rownumstates, rowMBstates)
      tempamp(1:namps, 1:rownumstates) = 0.
      !
      if (tilediff(k) .lt. 2) then
        !
        nabbr = rownumstates
        call abbrstates(Z_row, N_row, offset1, offset2, &
                        rownumstates, rowMBstates, nabbr, rowabbr1, rowabbr2)
        !
        call ApplyTBopTileTD(statesize, nparticles, &
                             rowstate, rownumstates, rowMBstates, nabbr, rowabbr1, rowabbr2, &
                             tmpstate, colnumstates, colMBstates, nabbc, colabbr1, colabbr2, &
                             namps, &
                             colamp(1:namps, icol:icol + colnumstates - 1), &
                             tempamp(1:namps, 1:rownumstates))
        !
      elseif (tilediff(k) .eq. 2) then
        !
        ! compare complete row and column groupID
        call MBnonzeroLocation(statesize, rowstate, colstate, &
                               rowdifloc, coldifloc, ndiffs)
        !
        rowdiffs(1) = rowstate(rowdifloc(1))
        rowdiffs(2) = rowstate(rowdifloc(2))
        orba = orb_sp(rowdiffs(1))
        orbb = orb_sp(rowdiffs(2))
        parab = pr_orb(orba)*pr_orb(orbb)
        j2a = j2_orb(orba)
        j2b = j2_orb(orbb)
        !
        coldiffs(1) = colstate(coldifloc(1))
        coldiffs(2) = colstate(coldifloc(2))
        orbc = orb_sp(coldiffs(1))
        orbd = orb_sp(coldiffs(2))
        parcd = pr_orb(orbc)*pr_orb(orbd)
        j2c = j2_orb(orbc)
        j2d = j2_orb(orbd)
        !
        ! apply triangle inequality
        if (j2a + j2b + 2*Jop .lt. abs(j2c - j2d)) cycle
        if (j2c + j2d + 2*Jop .lt. abs(j2a - j2b)) cycle
        !
        ! lex order a, b, c, d, if necessary
        ! and collect relevant input matrix elements
        nTBMEs = nTBMEs_max
        if (pairwiseless(orba, orbb, orbc, orbd)) then
          call Set_TBME_array(parab, parcd, orba, orbb, orbc, orbd, &
                              1, nTBMEs, TBMEarray)
          reorder = 0
        else
          call Set_TBME_array(parcd, parab, orbc, orbd, orba, orbb, &
                              1, nTBMEs, TBMEarray)
          reorder = (-1)**DeltaMj
        end if
        if (nTBMEs .eq. 0) cycle
        !
        abIDN = orba .eq. orbb
        cdIDN = orbc .eq. orbd
        DIAG = (orba .eq. orbc) .and. (orbb .eq. orbd)
        !
        call ApplyTBopTileBU(j2a, j2b, j2c, j2d, TwoMj, &
                             rownumstates, rowMBstates, rowdiffs, rowdifloc, &
                             colnumstates, colMBbitrep, coldiffs, coldifloc, &
                             colgrp_SPoffset, abIDN, cdIDN, DIAG, reorder, &
                             nTBMEs, TBMEarray, namps, &
                             colamp(1:namps, icol:icol + colnumstates - 1), &
                             tempamp(1:namps, 1:rownumstates))
        !
      end if
      !omp critical
      !omp simd collapse(2)
      do ii = 1, namps
        do jj = 1, rownumstates
          rowamp(ii, jrow + jj - 1) = rowamp(ii, jrow + jj - 1) + tempamp(ii, jj)
        end do
      end do
      !omp end critical
      !
    end do
  end do
  !$omp end do
  !
  !$omp end parallel
  !
  return
end subroutine ApplyTBop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ApplyTBopTileTD(statesize, npart, &
                           rowstate, nrows, rowMBstates, nabbr, rowabbr1, rowabbr2, &
                           colstate, ncols, colMBstates, nabbc, colabbr1, colabbr2, &
                           namps, colamp, rowamp)
  use TBops, only: MBtwobodyObs
  implicit none
  !
  integer, intent(in) :: statesize, npart, nrows, ncols, nabbr, nabbc, namps
  integer(kind=2), dimension(statesize) :: colstate, rowstate
  integer(kind=2), dimension(npart, nrows), intent(in) :: rowMBstates
  integer(kind=2), dimension(npart, ncols), intent(in) :: colMBstates
  integer(kind=8), dimension(nabbr), intent(in) :: rowabbr1, rowabbr2
  integer(kind=8), dimension(nabbc), intent(in) :: colabbr1, colabbr2
  real(kind=4), dimension(namps, ncols), intent(in) :: colamp
  real(kind=4), dimension(namps, nrows), intent(inout) :: rowamp
  !
  ! local variables
  integer :: i, j, k, ndiffs, c, r, t
  integer(kind=8) :: ii1, ii2, ii3, ii4, XOR1, XOR2
  integer, dimension(statesize) :: rowdifloc, coldifloc
  real(kind=8), dimension(1) :: xTBops
  !
  ! loop over row states in this tile
  do j = 1, nrows
    rowstate(1:npart) = rowMBstates(1:npart, j)
    ii3 = rowabbr1(j)
    ii4 = rowabbr2(j)
    do i = 1, ncols           ! loop over column states
      ! camp(1:ncamps) = colamp(1:ncamps, i)
      ii1 = colabbr1(i)
      XOR1 = ieor(ii1, ii3)
      ii2 = colabbr2(i)
      XOR2 = ieor(ii2, ii4)
      ndiffs = popcnt(XOR1) + popcnt(XOR2)
      if (ndiffs .le. 4) then
        ! detailed comparison
        colstate(1:npart) = colMBstates(1:npart, i)
        call MBnonzeroLocation(statesize, rowstate, colstate, &
                               rowdifloc, coldifloc, ndiffs)
        !
        if (ndiffs .le. 2) then
          call MBtwobodyObs(rowstate, colstate, ndiffs, &
                            rowdifloc(1:2), coldifloc(1:2), 1, xTBops)
          !$omp simd
          do k = 1, namps
            rowamp(k, j) = rowamp(k, j) + colamp(k, i)*xTBops(1)
          end do
        end if
        !
      end if
    end do
    !
  end do
  !
  return
end subroutine ApplyTBopTileTD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ApplyTBopTileBU(j2a, j2b, j2c, j2d, TwoMj, &
                           nrows, rowMBstates, rowdiffs, rowdifloc, &
                           ncols, colMBbitrep, coldiffs, coldifloc, &
                           colgrp_SPoffset, abIDN, cdIDN, DIAG, reorder, &
                           nTBMEs, TBMEarray, namps, colamp, rowamp)
  !
  use SPbasis, only: nparticles, next_sp_bin, mj2_sp, orb_sp, pr_orb
  use TBops, only: TBops_Eval, TBops_EvalDiag
  use Wigner3J, only: Jmax_2body
  implicit none
  !
  logical, intent(in) :: abIDN, cdIDN, DIAG
  integer, intent(in) :: nrows, ncols, namps, reorder, nTBMEs
  integer(kind=2), dimension(nparticles, nrows), intent(in) :: rowMBstates
  integer(kind=2), dimension(2), intent(in) :: rowdiffs, coldiffs
  integer, dimension(2), intent(in) :: rowdifloc, coldifloc
  integer, intent(in) :: j2a, j2b, j2c, j2d, TwoMj
  integer(kind=8), dimension(ncols), intent(in) :: colMBbitrep
  integer, dimension(nparticles), intent(in) :: colgrp_SPoffset
  real(kind=4), dimension(namps, ncols), intent(in) :: colamp
  real(kind=4), dimension(1, nTBMEs), intent(in) :: TBMEarray
  real(kind=4), dimension(namps, nrows), intent(inout) :: rowamp
  !
  !     work arrays
  integer, dimension(nparticles) :: MBstate, rowgrp_SPoffset
  integer(kind=8), dimension(0:Jmax_2body) :: MBbitrep
  integer(kind=8), dimension(0:Jmax_2body) :: locbits1, locbits2
  integer, dimension(0:Jmax_2body) :: colflg, colindx
  !
  !     local variables
  integer(kind=2), dimension(2) :: tmprowdiffs, tmpcoldiffs
  integer, dimension(3) :: tmprowdifloc, tmpcoldifloc
  integer, dimension(2) :: diffsmjmin, diffsmjmax
  integer :: rowmjmin, rowmjmax, colmjmin, colmjmax
  integer :: flag, prev, mm, cj, c, r, t
  integer :: i, j, k, kk, icolloc, irowloc, difloc
  integer :: col1mjmin, col2mjmin, col1mjmax, col2mjmax
  integer :: mjrow, mjrow1, mjrow2, mjcol1, mjcol2
  integer :: colTBstate1, colTBstate2, numcolTBstates
  integer :: TBgrp_SPoffset1, TBgrp_SPoffset2, itest
  integer(kind=8) :: rowbitrep
  real(kind=8), dimension(1) :: xops
  real(kind=8) :: xamps
  real :: phasefac, fac
  real, external :: myphase
  real(kind=8), parameter :: sqr2 = sqrt(2.d0)
  !
  call set_extended_diffs_2NF(rowdiffs, rowdifloc, &
                              tmprowdiffs, tmprowdifloc, diffsmjmin, diffsmjmax)
  rowmjmin = diffsmjmin(1) + diffsmjmin(2)
  rowmjmax = diffsmjmax(1) + diffsmjmax(2)
  tmprowdifloc(3) = 0   ! serves as end-cap
  !
  call set_extended_diffs_2NF(coldiffs, coldifloc, &
                              tmpcoldiffs, tmpcoldifloc, diffsmjmin, diffsmjmax)
  col1mjmin = diffsmjmin(1)
  col2mjmin = diffsmjmin(2)
  colmjmin = col1mjmin + col2mjmin
  !
  col1mjmax = diffsmjmax(1)
  col2mjmax = diffsmjmax(2)
  colmjmax = col1mjmax + col2mjmax
  !
  if (colmjmin + twomj .gt. rowmjmax) return
  if (rowmjmin .gt. colmjmax + twomj) return
  !
  tmpcoldifloc(3) = 0   ! serves as end-cap
  !
  TBgrp_SPoffset1 = colgrp_SPoffset(coldifloc(1))
  TBgrp_SPoffset2 = colgrp_SPoffset(coldifloc(2))
  !
  ! flag to control additional passes through outer loop as needed
  ! for cases where there are several particles in the same rowgroupID
  flag = 1
  do while (flag .eq. 1)
    !
    kk = 1
    irowloc = 1
    icolloc = 1
    do k = 1, nparticles
      ! tmprowdifloc(oprank+1) = 0 serves as 'end-cap'
      if (k .eq. tmprowdifloc(irowloc)) then
        rowgrp_SPoffset(k) = 0
        irowloc = irowloc + 1
      else
        kk = k + icolloc - irowloc
        ! tmpcoldifloc(oprank+1) = 0 serves as 'end-cap'
        do while (kk .eq. tmpcoldifloc(icolloc))
          icolloc = icolloc + 1
          kk = kk + 1
        end do
        rowgrp_SPoffset(k) = colgrp_SPoffset(kk)
      end if
    end do
    !
    phasefac = myphase(tmprowdifloc(1) + tmprowdifloc(2))
    !
    do i = 1, nrows
      do k = 1, nparticles
        MBstate(k) = rowMBstates(k, i)
      end do
      !
      mjrow1 = mj2_sp(MBstate(tmprowdifloc(1)))
      mjrow2 = mj2_sp(MBstate(tmprowdifloc(2)))
      mjrow = mjrow1 + mjrow2
      !
      if (mjrow .lt. colmjmin + twomj) cycle
      if (mjrow .gt. colmjmax + twomj) cycle
      !
      rowbitrep = 0
      do k = 1, nparticles
        if ((k .eq. tmprowdifloc(1)) .or. (k .eq. tmprowdifloc(2))) cycle
        itest = MBstate(k) - rowgrp_SPoffset(k)
        rowbitrep = IBset(rowbitrep, itest)
      end do
      !
      ! set first colTBstate
      colTBstate1 = tmpcoldiffs(1)
      colTBstate2 = next_sp_bin(tmpcoldiffs(2)) - 1
      col1mjmin = diffsmjmin(1)
      col2mjmax = diffsmjmax(2)
      ! calculate number of colTBstates
      mm = mjrow - twomj - col1mjmin - col2mjmax
      if (mm .gt. 0) then ! increase col1mjmin
        col1mjmin = col1mjmin + mm
        colTBstate1 = colTBstate1 + mm/2
      elseif (mm .lt. 0) then
        ! decrease col2mjmax by mm (note: mm is negative)
        col2mjmax = col2mjmax + mm
        colTBstate2 = colTBstate2 + mm/2
      end if
      numcolTBstates = min(col1mjmax - col1mjmin, col2mjmax - col2mjmin)
      numcolTBstates = numcolTBstates/2
      if ((colTBstate1 + numcolTBstates) .ge. (colTBstate2 - numcolTBstates)) then
        numcolTBstates = colTBstate2 - (colTBstate1 + colTBstate2)/2 - 1
      end if
      !
      ! set MBbitrep(j)
      do j = 0, numcolTBstates
        itest = colTBstate1 + j - TBgrp_SPoffset1
        if (Btest(rowbitrep, itest)) then
          colflg(j) = 0
          MBbitrep(j) = 0
        else
          colflg(j) = 1
          locbits1(j) = IBits(rowbitrep, 1, itest)
          MBbitrep(j) = IBset(rowbitrep, itest)
        end if
        !
        itest = colTBstate2 - j - TBgrp_SPoffset2
        if (Btest(rowbitrep, itest)) then
          colflg(j) = 0
        else
          locbits2(j) = IBits(rowbitrep, 1, itest)
          MBbitrep(j) = IBset(MBbitrep(j), itest)
        end if
      end do
      !
      ! locate MBbitrep(j)
      do j = 0, numcolTBstates
        if (colflg(j) .eq. 1) then
          do k = 1, ncols
            if (MBbitrep(j) .eq. colMBbitrep(k)) then
              colindx(j) = k
              ! exit ! exit statement is not OpenMP standard-compliant
            end if
          end do
        end if
      end do
      !
      do j = 0, numcolTBstates
        if (colflg(j) .eq. 1) then
          mjcol1 = col1mjmin + 2*j
          mjcol2 = col2mjmax - 2*j
          !
          ! TO BE DONE?
          ! move part/most of TBops_Eval outside of main loop,
          ! e.g. determination jmin, jmax, jstep for (ab),
          ! symmetry factor, and myphase((j2a-j2b+j2c-j2d)/2);
          ! retrieval of Wigner sub-arrays
          !
          if (DIAG) then
            call TBops_EvalDiag(abIDN, cdIDN, &
                                j2a, j2b, j2c, j2d, mjrow1, mjrow2, mjcol1, mjcol2, &
                                nTBMEs, TBMEarray, 1, xops)
          else
            if (reorder .eq. 0) then
              call TBops_Eval(abIDN, cdIDN, &
                              j2a, j2b, j2c, j2d, mjrow1, mjrow2, mjcol1, mjcol2, &
                              nTBMEs, TBMEarray, 1, xops)
            else
              call TBops_Eval(cdIDN, abIDN, &
                              j2c, j2d, j2a, j2b, mjcol1, mjcol2, mjrow1, mjrow2, &
                              nTBMEs, TBMEarray, 1, xops)
              xops(1) = reorder*xops(1)
            end if
          end if
          !
          difloc = popcnt(locbits1(j)) + popcnt(locbits2(j)) + 1
          fac = myphase(difloc)
          !
          cj = colindx(j)
          !$omp simd
          do k = 1, namps
            rowamp(k, i) = rowamp(k, i) + fac*phasefac*colamp(k, cj)*xops(1)
          end do
        end if
      end do
      !
    end do
    !
    ! Additional passes through outer loop as needed for cases
    ! where there are several particles in the same rowgroupID
    flag = 0
    !
    if (tmprowdiffs(2) .lt. rowdiffs(2)) then
      tmprowdiffs(2) = tmprowdiffs(2) + 1
      tmprowdifloc(2) = tmprowdifloc(2) + 1
      flag = 1
    end if
    !
    if (flag .eq. 0) then
      if (tmprowdiffs(1) .lt. rowdiffs(1)) then
        tmprowdiffs(1) = tmprowdiffs(1) + 1
        tmprowdifloc(1) = tmprowdifloc(1) + 1
        !
        prev = tmprowdiffs(1) + 1
        do while ((tmprowdiffs(2) .gt. prev) .and. &
                  (next_sp_bin(tmprowdiffs(2)) .eq. next_sp_bin(tmprowdiffs(2) - 1)))
          tmprowdiffs(2) = tmprowdiffs(2) - 1
          tmprowdifloc(2) = tmprowdifloc(2) - 1
        end do
        flag = 1
      end if
    end if
    !
  end do
  !
  return
end subroutine ApplyTBopTileBU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
