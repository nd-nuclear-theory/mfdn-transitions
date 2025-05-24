!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     contains
!       subroutine EvalOBDME
!       subroutine ComputeOBDMETileDiag
!       subroutine ComputeOBDMETileBU
!       subroutine set_extended_diffs_OB
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalOBDME(maxnMstates,                         &
     Z_col, N_col, ncolgroupids, colgroupidlist,          &
     TwoMj_col, colMstateptr, numcolstates, ncamps, camp, &
     Z_row, N_row, nrowgroupids, rowgroupidlist,          &
     TwoMj_row, rowMstateptr, numrowstates, nramps, ramp, &
     coltileptr, ntiles, rowind, tilediff,                &
     mnK, mxK, nobdme, obdmeK_ptr, obdme)
  !
  use SPbasis, only: norbt, nparticles, nspstates
  use SPbasis, only: orb_sp, mj2_sp, next_sp_bin, j2_orb
  use MjStates, only: MjStatesGen
  use Wigner3J, only: Retrieve_Wigner3J_Array
  implicit none
  integer, intent(in) :: maxnMstates, ntiles, mnK, mxK, nobdme
  integer, intent(in) :: Z_col, N_col, TwoMj_col, ncolgroupids, numcolstates, ncamps
  integer, intent(in) :: Z_row, N_row, TwoMj_row, nrowgroupids, numrowstates, nramps
  integer(kind=2), dimension(nparticles, ncolgroupids), intent(in) :: colgroupidlist
  integer(kind=2), dimension(nparticles, nrowgroupids), intent(in) :: rowgroupidlist
  integer, dimension(ncolgroupids+1), intent(in) :: colMstateptr
  integer, dimension(nrowgroupids+1), intent(in) :: rowMstateptr
  integer, dimension(ncolgroupids+1), intent(in) :: coltileptr
  integer, dimension(ntiles+1), intent(in) :: rowind, tilediff
  !
  real(kind=4), dimension(ncamps, numcolstates), intent(in) :: camp
  real(kind=4), dimension(nramps, numrowstates), intent(in) :: ramp
  real(kind=8), dimension(nobdme, nramps, ncamps), intent(out) :: obdme
  integer, dimension(norbt, norbt), intent(in) :: obdmeK_ptr
  !
  ! local variables
  integer :: statesize, deltamj
  integer(kind=2), dimension(nparticles, maxnMstates) :: colMBstates, rowMBstates
  integer(kind=2), dimension(nparticles+1) :: colstate, rowstate, tmpstate
  integer, dimension(nparticles+1) :: rowdifloc, coldifloc
  integer, dimension(nparticles) :: twjk, indx, colgrp_SPoffset
  integer :: i, j, ii, jj, ik, k, k1, k2, ndiffs, r, c
  integer :: iprev, icur, itest, n_grp_SPstates
  integer :: icol, jrow, korb, mj2, mntwj, mxtwj
  integer(kind=2) :: aa, cc
  integer :: orba, j2a, orbc, j2c, indxca
  integer :: colnumstates, rownumstates
  !
  real(kind=8), dimension(0:mxK) :: wigarray
  real(kind=8) :: phase, fac, xamps
  real, external :: myphase
  !
  integer(kind=8), dimension(maxnMstates) :: colMBbitrep, rowMBbitrep
  !
  obdme(1:nobdme, 1:nramps, 1:ncamps) = 0.d0
  statesize = nparticles + 1
  deltamj = TwoMj_row - TwoMj_col
  !
  !$omp parallel default(shared)                                               &
  !$omp          private(i, j, k, k1, k2, ik, ii, jj, r, c,                    &
  !$omp                  ndiffs, colMBbitrep, rowMBbitrep,                     &
  !$omp                  icol, colnumstates, colMBstates, colstate, coldifloc, &
  !$omp                  jrow, rownumstates, rowMBstates, rowstate, rowdifloc, &
  !$omp                  iprev, icur, itest, n_grp_SPstates, colgrp_SPoffset,  &
  !$omp                  aa, orba, j2a, cc, orbc, j2c, indxca,                 &
  !$omp                  korb, mj2, mntwj, mxtwj, twjk, indx,                  &
  !$omp                  wigarray, phase, fac, xamps, tmpstate)                &
  !$omp          reduction(+: obdme)
  !
  colstate(nparticles+1) = nspstates + 1
  rowstate(nparticles+1) = nspstates + 1
  tmpstate(nparticles+1) = nspstates + 1
  !
  !$omp do schedule(dynamic)
  do i = 1, ncolgroupids
     icol = colMstateptr(i)
     colstate(1:nparticles) = colgroupidlist(1:nparticles, i)
     !
     ! generate group of mbstates in M-scheme for colindex i
     colnumstates = maxnMstates
     call mjstatesgen(nparticles, nspstates, mj2_sp, TwoMj_col, &
          next_sp_bin, colstate, colnumstates, colMBstates)
     !
     ! Bit representation of colMBstates in terms of 'active' S.P. states
     ! Determine 'active' S.P. states
     iprev = colstate(1)
     colgrp_SPoffset(1) = iprev - 1
     n_grp_SPstates = next_sp_bin(iprev) - iprev
     do k = 2, nparticles
        icur = colstate(k)
        if ( next_sp_bin(icur) .eq. next_sp_bin(iprev) ) then
           ! icur and iprev are in the same SPbin (same partition)
           colgrp_SPoffset(k) = colgrp_SPoffset(k-1)
        else
           ! icur is first (?) SPstate in new SPbin
           colgrp_SPoffset(k) = icur - n_grp_SPstates - 1
           n_grp_SPstates = n_grp_SPstates+next_sp_bin(icur)-icur
           iprev = icur
        endif
     enddo
     ! check on n_grp_SPstates
     if (n_grp_SPstates > 64) then
        print*, n_grp_SPstates
        call cancelall(350)
     endif
     !
     ! Bit representation of colMBstates
     !$omp simd
     do j = 1, colnumstates
        colMBbitrep(j) = 0
        do k = 1, nparticles
           itest = colMBstates(k, j) - colgrp_SPoffset(k)
           colMBbitrep(j) = IBset(colMBbitrep(j), itest)
        enddo
     enddo
     !
     ! loop over all interacting rowgroups in this column group
     k1 = coltileptr(i)
     k2 = coltileptr(i+1) - 1
     do k = k1, k2
        if (tilediff(k) .eq. 0) then
           ! can only happen if Z_row=Z_col, N_row=N_col, Delta_Par=1
           if (TwoMj_row .eq. TwoMj_col) then
              ! true diagonal, sum over particles
              do jj = 1, nparticles
                 korb = orb_sp(colstate(jj))
                 twjk(jj) = j2_orb(korb)
                 indx(jj) = obdmeK_ptr(korb, korb)
              enddo
              !
              do ii = 0, colnumstates-1 ! columns
                 tmpstate(1:nparticles) = colMBstates(1:nparticles, ii+1)
                 do jj = 1, nparticles
                    mj2 = mj2_sp(tmpstate(jj))
                    j2a = twjk(jj)
                    mxtwj = min(j2a, mxK)
                    call Retrieve_Wigner3J_Array(j2a, j2a, mj2, -mj2,    &
                         0, mxtwj, wigarray(0:mxtwj) )
                    phase = myphase((j2a - mj2)/2)
                    ! phase = myphase((j2a + mj2)/2)
                    !$omp simd collapse(3)
                    do ik = mnK, mxtwj
                       do c = 1, ncamps
                          do r = 1, nramps
                             j = indx(jj) - mnK + ik
                             fac = wigarray(ik) * phase
                             xamps = camp(c, icol+ii) * ramp(r, icol+ii)
                             obdme(j,r,c) = obdme(j,r,c) + fac * xamps
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              !
           else
              j = rowind(k)
              jrow = rowMstateptr(j)
              rowstate(1:nparticles) = rowgroupidlist(:,j)
              !
              ! generate group of mbstates in M-scheme for rowindex
              rownumstates = maxnMstates
              call mjstatesgen(nparticles, nspstates, mj2_sp, TwoMj_row, &
                   next_sp_bin, rowstate, rownumstates, rowMBstates)
              !
              ! Bit representation of rowMBstates
              do jj = 1, rownumstates
                 rowMBbitrep(jj) = 0
                 do ik = 1, nparticles
                    itest = rowMBstates(ik, jj) - colgrp_SPoffset(ik)
                    rowMBbitrep(jj) = IBset(rowMBbitrep(jj), itest)
                 enddo
              enddo
              !
              call ComputeOBDMETileDiag(statesize, obdmeK_ptr,       &
                   rowstate, rownumstates, rowMBstates, rowMBbitrep, &
                   tmpstate, colnumstates, colMBstates, colMBbitrep, &
                   nramps, ramp(1:nramps, jrow:jrow+rownumstates-1), &
                   ncamps, camp(1:ncamps, icol:icol+colnumstates-1), &
                   mnK, mxK, nobdme, obdme)
           endif
           !
        elseif (tilediff(k) .eq. 1) then
           j = rowind(k)
           jrow = rowMstateptr(j)
           rowstate(1:nparticles) = rowgroupidlist(1:nparticles, j)
           !
           call MBnonzeroLocation(statesize, rowstate, colstate, &
                rowdifloc, coldifloc, ndiffs)
           !
           aa = rowstate(rowdifloc(1))
           cc = colstate(coldifloc(1))
           orba = orb_sp(aa)
           orbc = orb_sp(cc)
           j2a  = j2_orb(orba)
           j2c  = j2_orb(orbc)
           !
           mntwj = max(abs(j2a-j2c)/2, mnK)
           if (mntwj .gt. mxK) cycle
           !
           mxtwj = min((j2a+j2c)/2, mxK)
           if (mxtwj .lt. mnK) cycle
           !
           indxca = obdmeK_ptr(orbc, orba) - mntwj
           !
           ! generate group of mbstates in M-scheme for rowindex
           rownumstates = maxnMstates
           call mjstatesgen(nparticles, nspstates, mj2_sp, TwoMj_row, &
                next_sp_bin, rowstate, rownumstates, rowMBstates)
           !
           call ComputeOBDMETileBU(j2a, j2c, mntwj, mxtwj, deltamj,           &
                rownumstates, rowMBstates, aa, rowdifloc(1),                  &
                colnumstates, colMBbitrep, cc, coldifloc(1), colgrp_SPoffset, &
                nramps, ramp(1:nramps, jrow:jrow+rownumstates-1),             &
                ncamps, camp(1:ncamps, icol:icol+colnumstates-1),             &
                mnK, mxK, nobdme, obdme, indxca)
        endif
        !
     enddo
  enddo
  !$omp end do
  !
  !$omp end parallel
  !
  return
end subroutine EvalOBDME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ComputeOBDMETileDiag(statesize, obdmeK_ptr,   &
     rowstate, nrows, rowMBstates, rowMBbitrep,          &
     colstate, ncols, colMBstates, colMBbitrep,          &
     nramps, rowamp, ncamps, colamp, mnK, mxK, nobdme, obdme)
  !
  ! Given a tile of ncols X nrows
  ! evaluate OneBodyDensityMatrixElements < ramp | O | camp>
  !
  use SPbasis, only: norbt, nparticles, mj2_sp, orb_sp, j2_orb
  use Wigner3J, only: Retrieve_Wigner3J_Array
  implicit none
  !
  integer, intent(in) :: statesize, nrows, ncols, mnK, mxK, nobdme, nramps, ncamps
  integer, dimension(norbt, norbt), intent(in) :: obdmeK_ptr
  integer(kind=2), dimension(statesize) :: colstate, rowstate
  integer(kind=2), dimension(nparticles,nrows), intent(in) :: rowMBstates
  integer(kind=2), dimension(nparticles,ncols), intent(in) :: colMBstates
  integer(kind=8), dimension(nrows), intent(in) :: rowMBbitrep
  integer(kind=8), dimension(ncols), intent(in) :: colMBbitrep
  real(kind=4), dimension(nramps, nrows), intent(in) :: rowamp
  real(kind=4), dimension(ncamps, ncols), intent(in) :: colamp
  real(kind=8), dimension(nobdme, nramps, ncamps), intent(inout) :: obdme
  !     local variables
  integer :: i, j, k, kk, ndiffs, mxtwj, indx, r, c
  integer(kind=8) :: iicol, jjrow, XOR
  integer, dimension(statesize) :: rowdifloc, coldifloc
  integer :: aa, orba, j2a, m2a, cc, orbc, j2c, m2c
  real(kind=8), dimension(0:mxK) :: wigarray
  real(kind=8) :: phase, fac, ramp(1:nramps), xamps
  real, external :: myphase
  !
  ! loop over row states in this tile
  do j = 1, nrows
     rowstate(1:nparticles) = rowMBstates(1:nparticles, j)
     ramp(1:nramps) = rowamp(1:nramps, j)
     jjrow = rowMBbitrep(j)
     do i = 1, ncols           ! loop over columns
        ! camp(1:ncamps) = colamp(1:ncamps,i)
        iicol = colMBbitrep(i)
        XOR = ieor(iicol, jjrow)
        ndiffs = popcnt(XOR)
        if (ndiffs .eq. 2) then
           ! detailed comparison
           colstate(1:nparticles) = colMBstates(1:nparticles, i)
           call MBnonzeroLocation(statesize, rowstate, colstate,  &
                rowdifloc, coldifloc, ndiffs)
           !
           if (ndiffs .eq. 1) then ! should be superfluous
              phase = myphase(rowdifloc(1) - coldifloc(1))
              !
              cc = colstate(coldifloc(1))
              aa = rowstate(rowdifloc(1))
              !
              orbc = orb_sp(cc)
              orba = orb_sp(aa)
              j2c = j2_orb(orbc)
              j2a = j2_orb(orba)
              !
              mxtwj = min((j2a+j2c)/2, mxK)
              if (mxtwj .lt. mnK) cycle
              !
              indx = obdmeK_ptr(orbc, orba) - mnK
              !
              m2c = mj2_sp(cc)
              m2a = mj2_sp(aa)
              !
              call Retrieve_Wigner3J_Array(j2c, j2a, m2c, -m2a,   &
                   mnK, mxtwj, wigarray(mnK:mxtwj) )
              ! phase = myphase( (j2c - m2c)/2 ) !
              phase = phase * myphase( (j2a - m2a)/2 )
              !! phase = phase * myphase( (j2c-j2a-m2c+m2a)/2)
              ! phase = myphase( (j2a + m2a)/2 ) !
              !
              !$omp simd collapse(3)
              do k = mnK, mxtwj
                 do c = 1, ncamps
                    do r = 1, nramps
                       kk = indx + k
                       fac = phase * wigarray(k)
                       xamps = ramp(r)* colamp(c,i)
                       obdme(kk,r,c) = obdme(kk,r,c) + fac * xamps
                    enddo
                 enddo
              enddo
           end if
           !
        end if
     end do
  end do
  !
  return
  !
end subroutine ComputeOBDMETileDiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ComputeOBDMETileBU(j2row, j2col, mntwj, mxtwj, deltamj, &
     nrows, rowMBstates, rowdiff, rowdifloc,                       &
     ncols, colMBbitrep, coldiff, coldifloc, colgrp_SPoffset,      &
     nramps, rowamp, ncamps, colamp, mnK, mxK, nobdme, obdme, indx)
  !
  ! Given a tile of ncols X nrows
  ! evaluate OneBodyDensityMatrixElements < ramp | O | camp>
  !
  use SPbasis, only: nparticles, nspstates, mj2_sp, next_sp_bin
  use Wigner3J, only: Retrieve_Wigner3J_Array
  implicit none
  !
  integer, intent(in) :: j2row, j2col, mntwj, mxtwj, deltamj, indx
  integer, intent(in) :: nrows, nramps, ncols, ncamps, mnK, mxK, nobdme
  integer(kind=2), dimension(nparticles,nrows), intent(in) :: rowMBstates
  integer(kind=2), intent(in) :: rowdiff, coldiff
  integer, intent(in) :: rowdifloc, coldifloc
  integer, dimension(nparticles), intent(in) :: colgrp_SPoffset
  integer(kind=8), dimension(ncols), intent(in) :: colMBbitrep
  real(kind=4), dimension(nramps, nrows), intent(in) :: rowamp
  real(kind=4), dimension(ncamps, ncols), intent(in) :: colamp
  real(kind=8), dimension(nobdme, nramps, ncamps), intent(inout) :: obdme
  !     local variables
  integer(kind=2) :: tmprowdiff, tmpcoldiff, colOBstate
  integer :: i, j, k, kk, r, c
  integer :: tmprowdifloc, tmpcoldifloc, difloc
  integer :: rowmjmin, rowmjmax, colmjmin, colmjmax
  integer :: OBgrp_SPoffset, flag, irowloc, icolloc
  integer :: mjrow, mjcol, itest, colindx
  integer(kind=8) :: rowbitrep, MBbitrep, locbits
  integer, dimension(nparticles) :: rowgrp_SPoffset
  integer(kind=2), dimension(nparticles) :: MBstate
  !
  real(kind=8), dimension(0:mxK) :: wigarray
  real(kind=8) :: phase, fac, xamps
  real, external :: myphase
  !
  call set_extended_diffs_OB(nspstates, mj2_sp, next_sp_bin,                 &
       rowdiff, rowdifloc, tmprowdiff, tmprowdifloc, rowmjmin, rowmjmax)
  !
  call set_extended_diffs_OB(nspstates, mj2_sp, next_sp_bin,                 &
       coldiff, coldifloc, tmpcoldiff, tmpcoldifloc, colmjmin, colmjmax)
  !
  if (colmjmin + deltamj .gt. rowmjmax) return
  if (rowmjmin .gt. colmjmax + deltamj) return
  !
  OBgrp_SPoffset = colgrp_SPoffset(coldifloc)
  !
  ! flag to control additional passes through outer loop as needed
  ! for cases where there are several particles in the same rowgroupID
  flag = 1
  !
  do while (flag .eq. 1)
     !
     kk = 1
     irowloc = 1
     icolloc = 1
     do k = 1, nparticles
        if (k .eq. tmprowdifloc) then
           rowgrp_SPoffset(k) = 0
           irowloc = irowloc+1
        else
           kk = k+icolloc-irowloc
           if (kk .eq. tmpcoldifloc) then
              icolloc = icolloc+1
              kk = kk + 1
           endif
           rowgrp_SPoffset(k) = colgrp_SPoffset(kk)
        endif
     enddo
     !
     do i = 1, nrows
        do k = 1, nparticles
           MBstate(k) = rowMBstates(k, i)
        enddo
        !
        mjrow = mj2_sp(MBstate(tmprowdifloc))
        !
        if (mjrow .lt. colmjmin + deltamj) cycle
        if (mjrow .gt. colmjmax + deltamj) cycle
        !
        rowbitrep = 0
        do k = 1, nparticles
           if (k .eq. tmprowdifloc) cycle
           itest = MBstate(k) - rowgrp_SPoffset(k)
           rowbitrep = IBset(rowbitrep, itest)
        enddo
        !
        ! set colOBstate
        mjcol = mjrow - deltamj
        colOBstate = tmpcoldiff + (mjcol - colmjmin)/2
        !
        ! set MBbitrep
        itest = colOBstate - OBgrp_SPoffset
        if ( Btest(rowbitrep, itest) ) then
           cycle  ! state already occupied
        else
           MBbitrep = IBset(rowbitrep, itest)
           locbits = IBits(rowbitrep, 1, itest)
           difloc = popcnt(locbits) + 1
           phase = myphase(tmprowdifloc - difloc)
        endif
        phase = phase * myphase( (j2row - mjrow)/2 ) ! quantum numbers of row state
        !
        ! locate MBbitrep
        do k = 1, ncols
           if (MBbitrep .eq. colMBbitrep(k)) then
              colindx = k
              ! exit ! exit statement is not OpenMP standard-compliant
           endif
        enddo
        !
        call Retrieve_Wigner3J_Array(j2col, j2row, mjcol,   &
             -mjrow, mntwj, mxtwj, wigarray(mntwj:mxtwj) )
        !
        !$omp simd collapse(3)
        do j = mntwj, mxtwj
           do c = 1, ncamps
              do r = 1, nramps
                 kk = indx + j
                 fac = phase * wigarray(j)
                 xamps = colamp(c, colindx) * rowamp(r,i)
                 obdme(kk,r,c) = obdme(kk,r,c) + fac * xamps
              enddo
           enddo
        enddo
        !
     enddo
     !
     ! Additional passes through outer loop as needed for cases
     ! where there are several particles in the same rowgroupID
     flag = 0
     !
     if (flag .eq. 0) then
        if ( tmprowdiff .lt. rowdiff ) then
           tmprowdiff = tmprowdiff + 1
           tmprowdifloc = tmprowdifloc + 1
           !
           flag = 1
        endif
     endif
     !
  enddo
  !
  return
  !
end subroutine ComputeOBDMETileBU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_extended_diffs_OB(nspstates, mj, nextbin,     &
     diffs, difloc, extdiffs, extdifloc, mjmin, mjmax)
  implicit none
  !
  integer, intent(in) :: nspstates
  integer, dimension(nspstates), intent(in) :: mj, nextbin
  integer(kind=2), intent(in) :: diffs
  integer, intent(in) :: difloc
  integer(kind=2), intent(out) :: extdiffs
  integer, intent(out) :: extdifloc
  integer, intent(out) :: mjmin, mjmax
  !
  integer :: state1, iloc1
  !
  ! if necessary, extend diffs to include
  ! the first SPstate/location of a (partitioned) orbital
  iloc1 = difloc
  state1 = diffs
  if (nextbin(state1) .eq. nextbin(1)) then
     state1 = 1
     iloc1 = 1
  else
     do while (nextbin(state1) .eq. nextbin(state1-1) )
        state1 = state1 - 1
        iloc1 = iloc1 - 1
     enddo
  endif
  mjmin = mj(state1)
  mjmax = mj(nextbin(state1)-1)

  extdiffs = state1
  extdifloc = iloc1
  !
  return
  !
end subroutine set_extended_diffs_OB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
