!
!     contains
!
!       subroutine CountOBDMEs
!       subroutine SetOBDMEindices
!       subroutine ReduceOBDME
!
!       subroutine CountOBDMEs_Tz
!       subroutine SetOBDMEindices_Tz
!       subroutine ReduceOBDME_Tz
!
!       subroutine WriteOBDME
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     obdmeKvals(nobdme)
!     ordered according to index =
!     = obdmeK_ptr(orba, orbc) : obdmeK_ptr(orba, orbc) + mxtwj-mintwj
!
!     robdmes(nobdme)
!      ordered according to
!       K = 0:  orba = 1:norb_p, orbc = 1:norb_p
!       K = 1:  orba = 1:norb_p, orbc = 1:norb_p
!       K = 2:  orba = 1:norb_p, orbc = 1:norb_p
!       etc.
!       K = 0:  orba = 1:norb_n, orbc = 1:norb_n
!       K = 1:  orba = 1:norb_n, orbc = 1:norb_n
!       K = 2:  orba = 1:norb_n, orbc = 1:norb_n
!       etc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CountOBDMEs(mnK, mxK, Delta_Par, nobdme_p, nobdme_n, &
     obdmeOrb1offset, obdmeOrb2offset)
  use SPbasis, only: norb_p, norb_n, j2_orb, pr_orb
  implicit none
  integer, intent(in) :: mnK, mxK, Delta_Par
  integer, intent(out) :: nobdme_p, nobdme_n
  integer, dimension(0:mxK+1), intent(out) :: obdmeOrb1offset, obdmeOrb2offset
  !
  integer :: ia, ipa, j2a, ic, mnj, mxj, k, nobdme
  integer, dimension(0:mxK) :: NobdmeOrb
  !
  nobdme = 0
  NobdmeOrb(0:mxK) = 0
  !$omp parallel do default(shared)                        &
  !$omp             private(ia, ipa, j2a, ic, mnj, mxj, k) &
  !$omp             reduction(+: nobdme, NobdmeOrb)
  do ia = 1, norb_p        ! loop over a^tilde
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do ic = 1, norb_p     ! loop over a^dagger
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        if (mxj .lt. mnK) cycle
        mxj = min(mxj, mxK)
        nobdme = nobdme + 1 + mxj - mnj
        do k = mnj, mxj
           NobdmeOrb(k) = NobdmeOrb(k) + 1
        enddo
     enddo
  enddo
  !$omp end parallel do
  nobdme_p = nobdme
  obdmeOrb1offset(0) = 0
  do k = 1, mxK+1
     obdmeOrb1offset(k) = obdmeOrb1offset(k-1) + NobdmeOrb(k-1)
  enddo
  !
  nobdme = 0
  NobdmeOrb(0:mxK) = 0
  !$omp parallel do default(shared)                        &
  !$omp             private(ia, ipa, j2a, ic, mnj, mxj, k) &
  !$omp             reduction(+: nobdme, NobdmeOrb)
  do ia = norb_p + 1, norb_p + norb_n
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do ic = norb_p + 1, norb_p + norb_n
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        if (mxj .lt. mnK) cycle
        mxj = min(mxj, mxK)
        nobdme = nobdme + 1 + mxj - mnj
        do k = mnj, mxj
           NobdmeOrb(k) = NobdmeOrb(k) + 1
        enddo
     enddo
  enddo
  !$omp end parallel do
  nobdme_n = nobdme
  obdmeOrb2offset(0) = nobdme_p
  do k = 1, mxK+1
     obdmeOrb2offset(k) = obdmeOrb2offset(k-1) + NobdmeOrb(k-1)
  enddo
  !
  return
end subroutine CountOBDMEs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetOBDMEindices(mnK, mxK, Delta_Par, obdmeOrb1offset, obdmeOrb2offset, &
     nobdme, obdmeK_ptr, obdmeOrbKbraket, obdmeOrbKindex)
  use SPbasis, only: norbt, norb_p, norb_n, j2_orb, pr_orb
  implicit none
  integer, intent(in) :: mnK, mxK, Delta_Par, nobdme
  integer, dimension(0:mxK+1), intent(in) :: obdmeOrb1offset, obdmeOrb2offset
  integer, dimension(norbt, norbt), intent(out) :: obdmeK_ptr
  integer, dimension(2, nobdme), intent(out) :: obdmeOrbKbraket
  integer, dimension(nobdme), intent(out) :: obdmeOrbKindex
  ! local variables
  integer :: ia, ipa, j2a, ic, mnj, mxj, k, n
  integer, dimension(0:mxK) :: NobdmeOrb
  !
  n = 1
  !
  NobdmeOrb(0:mxK) = obdmeOrb1offset(0:mxK)
  do ia = 1, norb_p        ! loop over a^tilde
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do ic = 1, norb_p     ! loop over a^dagger
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        if (mxj .lt. mnK) cycle
        mxj = min(mxj, mxK)
        obdmeK_ptr(ic, ia) = n
        do k = mnj, mxj
           NobdmeOrb(k) = NobdmeOrb(k) + 1
           obdmeOrbKbraket(1, NobdmeOrb(k)) = ia
           obdmeOrbKbraket(2, NobdmeOrb(k)) = ic
           obdmeOrbKindex(n) = NobdmeOrb(k)
           n = n + 1
        enddo
     enddo
  enddo
  !
  NobdmeOrb(0:mxK) = obdmeOrb2offset(0:mxK)
  do ia = norb_p + 1, norb_p + norb_n
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do ic = norb_p + 1, norb_p + norb_n
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        if (mxj .lt. mnK) cycle
        mxj = min(mxj, mxK)
        obdmeK_ptr(ic, ia) = n
        do k = mnj, mxj
           NobdmeOrb(k) = NobdmeOrb(k) + 1
           obdmeOrbKbraket(1, NobdmeOrb(k)) = ia
           obdmeOrbKbraket(2, NobdmeOrb(k)) = ic
           obdmeOrbKindex(n) = NobdmeOrb(k)
           n = n + 1
        enddo
     enddo
  enddo
  !
  return
end subroutine SetOBDMEindices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReduceOBDME(TwoJ_bra, TwoJ_ket, TwoMj_bra, TwoMj_ket, &
     mnK, mxK, Delta_Par, nobdme, obdme, obdmeK_ptr, obdmeOrbKindex, &
     reducefac, robdme)
  use SPbasis, only: norbt, norb_p, norb_n, pr_orb, j2_orb
  use Wigner3J, only: Wig3J
  implicit none
  integer, intent(in) :: TwoJ_bra, TwoJ_ket, TwoMj_bra, TwoMj_ket, mnK, mxK, Delta_Par, nobdme
  real(kind=8), dimension(nobdme), intent(in) :: obdme
  real(kind=8), dimension(nobdme), intent(out) :: robdme
  real(kind=8), dimension(0:mxK), intent(out) :: reducefac
  integer, dimension(norbt, norbt), intent(in) :: obdmeK_ptr
  integer, dimension(nobdme), intent(in) :: obdmeOrbKindex
  !
  ! local variables
  integer :: k, twok, mnj, mxj, indx, ii, ri, d2MK
  integer :: ia, j2a, ipa, ic
  real, external :: myphase
  !
  ! Wigner-Eckart reduction factor
  d2MK = TwoMj_bra - TwoMj_ket
  do k = mnK, mxK
     twok = 2 * k
     reducefac(k) = Wig3J(TwoJ_bra, twok, TwoJ_ket, &
                        -TwoMj_bra, d2MK, TwoMj_ket) * myphase( (TwoJ_bra-TwoMj_bra)/2 )

     if (reducefac(k) .eq. 0.d0) then
        print*, 'WARNING: Wigner-Eckart theorem gives 0/0 for K is', K
     else
        ! reducefac(k) = 1.d0/(rt4pi * reducefac(k))
        ! reducefac(k) = sqrt(2.d0*k + 1.d0) / (rt4pi * reducefac(k))
        reducefac(k) = sqrt(2.d0*k + 1.d0) / reducefac(k)
     endif
  enddo
  !
  ! reduced pp OBDMEs
  robdme(1:nobdme) = 0.d0
  !$omp parallel default(shared)                                      &
  !$omp          private(ia, ipa, j2a, ic, mnj, mxj, indx, k, ii, ri) &
  !$omp          reduction(+: robdme)
  !$omp do
  do ia = 1, norb_p         ! loop over a^tilde
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do ic = 1, norb_p      ! loop over a^dagger
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        mxj = min(mxj, mxK)
        indx = obdmeK_ptr(ic, ia)
        do k = mnj, mxj
           ii = indx - mnj + k
           ri = obdmeOrbKindex(ii)
           robdme(ri) = obdme(ii) * reducefac(k)
        enddo
     enddo
  enddo
  !$omp end do
  !
  ! reduced nn OBDMEs
  !$omp do
  do ia = norb_p + 1, norb_p + norb_n
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do ic = norb_p + 1, norb_p + norb_n
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        mxj = min(mxj, mxK)
        indx = obdmeK_ptr(ic, ia)
        do k = mnj, mxj
           ii = indx - mnj + k
           ri = obdmeOrbKindex(ii)
           robdme(ri) = obdme(ii) * reducefac(k)
        enddo
     enddo
  enddo
  !$omp end do
  !$omp end parallel
  !
  return
end subroutine ReduceOBDME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CountOBDMEs_Tz(mnK, mxK, Delta_Par,                         &
     offset_a, norb_a, offset_c, norb_c, nobdme, obdmeOrb_offset)
  use SPbasis, only: j2_orb, pr_orb
  implicit none
  integer, intent(in) :: mnK, mxK, Delta_Par
  integer, intent(in) :: offset_a, norb_a, offset_c, norb_c
  integer, intent(out) :: nobdme
  integer, dimension(0:mxK+1), intent(out) :: obdmeOrb_offset
  !
  integer :: a, ia, ipa, j2a, c, ic, mnj, mxj, k
  integer, dimension(0:mxK) :: NobdmeOrb
  !
  nobdme = 0
  NobdmeOrb(0:mxK) = 0
  !$omp parallel do default(shared)                              &
  !$omp             private(a, ia, ipa, j2a, c, ic, mnj, mxj, k) &
  !$omp             reduction(+: nobdme, NobdmeOrb)
  do a = 1, norb_a        ! loop over a^tilde
     ia = offset_a + a
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do c = 1, norb_c     ! loop over a^dagger
        ic = offset_c + c
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        if (mxj .lt. mnK) cycle
        mxj = min(mxj, mxK)
        nobdme = nobdme + 1 + mxj - mnj
        do k = mnj, mxj
           NobdmeOrb(k) = NobdmeOrb(k) + 1
        enddo
     enddo
  enddo
  !$omp end parallel do
  obdmeOrb_offset(0) = 0
  do k = 1, mxK+1
     obdmeOrb_offset(k) = obdmeOrb_offset(k-1) + NobdmeOrb(k-1)
  enddo
  !
  return
end subroutine CountOBDMEs_Tz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetOBDMEindices_Tz(mnK, mxK, Delta_Par,                   &
     offset_a, norb_a, offset_c, norb_c, obdmeOrb_offset,            &
     nobdme, obdmeK_ptr, obdmeOrbKbraket, obdmeOrbKindex)
  use SPbasis, only: norbt, j2_orb, pr_orb
  implicit none
  integer, intent(in) :: mnK, mxK, Delta_Par, nobdme
  integer, intent(in) :: offset_a, norb_a, offset_c, norb_c
  integer, dimension(0:mxK+1), intent(in) :: obdmeOrb_offset
  integer, dimension(norbt, norbt), intent(out) :: obdmeK_ptr
  integer, dimension(2, nobdme), intent(out) :: obdmeOrbKbraket
  integer, dimension(nobdme), intent(out) :: obdmeOrbKindex
  ! local variables
  integer :: a, ia, ipa, j2a, c, ic, mnj, mxj, k, n
  integer, dimension(0:mxK) :: NobdmeOrb
  !
  n = 1
  !
  NobdmeOrb(0:mxK) = obdmeOrb_offset(0:mxK)
  do a = 1, norb_a        ! loop over a^tilde
     ia = offset_a + a
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do c = 1, norb_c     ! loop over a^dagger
        ic = offset_c + c
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        if (mxj .lt. mnK) cycle
        mxj = min(mxj, mxK)
        obdmeK_ptr(ic, ia) = n
        do k = mnj, mxj
           NobdmeOrb(k) = NobdmeOrb(k) + 1
           obdmeOrbKbraket(1, NobdmeOrb(k)) = ia ! or a ?
           obdmeOrbKbraket(2, NobdmeOrb(k)) = ic ! or c ?
           obdmeOrbKindex(n) = NobdmeOrb(k)
           n = n + 1
        enddo
     enddo
  enddo
  !
  return
end subroutine SetOBDMEindices_Tz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReduceOBDME_Tz(TwoJ_bra, TwoJ_ket, TwoMj_bra, TwoMj_ket, &
     mnK, mxK, Delta_Par, offset_a, norb_a, offset_c, norb_c,       &
     nobdme, obdme, obdmeK_ptr, obdmeOrbKindex, reducefac, robdme)
  use SPbasis, only: norbt, pr_orb, j2_orb
  use Wigner3J, only: Wig3J
  implicit none
  integer, intent(in) :: TwoJ_bra, TwoJ_ket, TwoMj_bra, TwoMj_ket, nobdme
  integer, intent(in) :: mnK, mxK, Delta_Par, offset_a, norb_a, offset_c, norb_c
  real(kind=8), dimension(nobdme), intent(in) :: obdme
  real(kind=8), dimension(nobdme), intent(out) :: robdme
  real(kind=8), dimension(0:mxK), intent(out) :: reducefac
  integer, dimension(norbt, norbt), intent(in) :: obdmeK_ptr
  integer, dimension(nobdme), intent(in) :: obdmeOrbKindex
  !
  ! local variables
  integer :: k, twok, mnj, mxj, indx, ii, ri, d2MK
  integer :: a, ia, j2a, ipa, c, ic
  real, external :: myphase
  !
  print*, 'in ReduceOBDME_Tz'
  ! Wigner-Eckart reduction factor
  d2MK = TwoMj_bra - TwoMj_ket
  do k = mnK, mxK
     twok = 2 * k
     reducefac(k) = Wig3J(TwoJ_bra, twok, TwoJ_ket, &
                        -TwoMj_bra, d2MK, TwoMj_ket) * myphase( (TwoJ_bra-TwoMj_bra)/2 )

     if (reducefac(k) .eq. 0.d0) then
        print*, 'WARNING: Wigner-Eckart theorem gives 0/0 for K is', K
     else
        ! reducefac(k) = 1.d0/(rt4pi * reducefac(k))
        ! reducefac(k) = sqrt(2.d0*k + 1.d0) / (rt4pi * reducefac(k))
        reducefac(k) = sqrt(2.d0*k + 1.d0) / reducefac(k)
     endif
  enddo
  !
  robdme(1:nobdme) = 0.d0
  do a = 1, norb_a         ! loop over a^tilde
     ia = offset_a + a
     ipa = pr_orb(ia)
     j2a = j2_orb(ia)
     do c = 1, norb_c      ! loop over a^dagger
        ic = offset_c + c
        if (pr_orb(ic) .ne. ipa*Delta_Par) cycle
        mnj = abs(j2a - j2_orb(ic)) / 2
        if (mnj .gt. mxK) cycle
        mnj = max(mnj, mnK)
        mxj = (j2a + j2_orb(ic)) / 2
        mxj = min(mxj, mxK)
        indx = obdmeK_ptr(ic, ia)
        do k = mnj, mxj
           ii = indx - mnj + k
           ri = obdmeOrbKindex(ii)
           robdme(ri) = obdme(ii) * reducefac(k)
        enddo
     enddo
  enddo
  !
  return
end subroutine ReduceOBDME_Tz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteOBDME(seq_bra, seq_ket, mnK, mxK, nobdme, nobdme_p, nobdme_n, &
     obdmeOrb1offset, obdmeOrb2offset, obdmeOrbKbraket, reducefac, vals, nket)
  use SPbasis, only: norb_p, norb_n, n_orb, l_orb, j2_orb, wt_orb
  use Wavefunctions_m, only: wf_info_bra, wf_info_ket
  implicit none
  integer, intent(in) :: seq_bra, seq_ket, nket
  integer, intent(in) :: mnK, mxK, nobdme, nobdme_p, nobdme_n
  integer, dimension(0:mxK+1), intent(in) :: obdmeOrb1offset, obdmeOrb2offset
  integer, dimension(2, nobdme), intent(in) :: obdmeOrbKbraket
  real(kind=8), dimension(nobdme), intent(in) :: vals
  real(kind=8), dimension(0:mxK), intent(in) :: reducefac
  !
  ! local variables
  integer, parameter :: fh=110
  integer :: ia, ic, k, j, i, versionno, tz
  integer :: nobdme_p_adj, nobdme_n_adj
  real, external :: myphase
  character(LEN=256) :: obdmefile
  !
  ! example file name: mfdn.robdme.2J03.g1.n01.2J01.g1.n01
  write(obdmefile, '(A,A,I0.2,A,I0.1,A,I0.2,A,I0.2,A,I0.1,A,I0.2)') 'transitions.robdme', &
    '.2J', wf_info_bra%TwoJ(seq_bra), '.g', (1-wf_info_bra%parity)/2, '.n', wf_info_bra%Jseq(seq_bra), &
    '.2J', wf_info_ket%TwoJ(seq_ket), '.g', (1-wf_info_ket%parity)/2, '.n', wf_info_ket%Jseq(seq_ket)
  open(unit=fh, file=TRIM(obdmefile), status='replace', action='write')
  versionno = 1520
  ! hack-y temporary fix to write correct number of OBDMEs
  ! TODO: don't even store/calculate OBDMEs which are necessarily zero
  nobdme_p_adj = nobdme_p
  nobdme_n_adj = nobdme_n
  do k = mnK, mxK
     if (reducefac(k) .eq. 0d0) then
        nobdme_p_adj = nobdme_p_adj - (obdmeOrb1offset(k+1) - obdmeOrb1offset(k))
        nobdme_n_adj = nobdme_n_adj - (obdmeOrb2offset(k+1) - obdmeOrb2offset(k))
     endif
  enddo
  ! store source quantum numbers for independent life of OBDME file
1 format(2i3,x,i4,2i4,sp,i4,ss,x,i4,2f12.5,'  ! bra Z N seq 2J 2Mj par n T En')
2 format(2i3,x,i4,2i4,sp,i4,ss,x,i4,2f12.5,'  ! ket Z N seq 2J 2Mj par n T En')
  !
11 format('#  LABELS for Single-Particle Orbitals')
12 format('#  ia, na, la, 2ja, 2tz, wt')
16 format(i8, 2x, 4i4, 3x, f12.6)
  !
21 format('#  LABELS for Reduced Multipole expansion of OBDMEs')
22 format('#  ia, ib,   K,  ROBDME')
26 format(2i8, 4x, i4, 3x, e16.8)
  !
  write(fh,*) versionno,     '  ! version number'
  write(fh,1) wf_info_bra%num_protons, wf_info_bra%num_neutrons, seq_bra, &
              wf_info_bra%TwoJ(seq_bra), wf_info_bra%TwoMj, wf_info_bra%parity, &
              wf_info_bra%Jseq(seq_bra), wf_info_bra%isospin(seq_bra), wf_info_bra%eigval(seq_bra)
  write(fh,2) wf_info_ket%num_protons, wf_info_ket%num_neutrons, seq_ket, &
              wf_info_ket%TwoJ(seq_ket), wf_info_ket%TwoMj, wf_info_ket%parity, &
              wf_info_ket%Jseq(seq_ket), wf_info_ket%isospin(seq_ket), wf_info_ket%eigval(seq_ket)
  write(fh,*) mnK,           '  ! min. K in multipole expansion'
  write(fh,*) mxK,           '  ! max. K in multipole expansion'
  write(fh,*) norb_p, norb_n,   '  ! number of p and n sp orbitals'
  write(fh,*) nobdme_p_adj,nobdme_n_adj,'  ! number of p and n OBDMEs'
  !
  write(fh, 11)
  write(fh, 12)
  tz = 1
  do i = 1, norb_p
     write(fh, 16) i, n_orb(i), l_orb(i), j2_orb(i), tz, wt_orb(i)
  enddo
  tz = -1
  do i = norb_p+1, norb_p+norb_n
     write(fh, 16) i, n_orb(i), l_orb(i), j2_orb(i), tz, wt_orb(i)
  enddo
  !
  write(fh, 21)
  write(fh, 22)
  i = 0
  do k = mnK, mxK
     if (reducefac(k) .eq. 0.d0) cycle
     do j = obdmeOrb1offset(k)+1, obdmeOrb1offset(k+1)
        ia = obdmeOrbKbraket(1, j)
        ic = obdmeOrbKbraket(2, j)
        write(fh, 26) ia, ic, k, vals(j)
     enddo
  enddo
  !
  write(fh, *)
  !
  do k = mnK, mxK
     if (reducefac(k) .eq. 0.d0) cycle
     do j = obdmeOrb2offset(k)+1, obdmeOrb2offset(k+1)
        ia = obdmeOrbKbraket(1, j)
        ic = obdmeOrbKbraket(2, j)
        write(fh, 26) ia, ic, k, vals(j)
     enddo
  enddo
  !
  close(unit=fh, status='keep')
  !
  return
end subroutine WriteOBDME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
