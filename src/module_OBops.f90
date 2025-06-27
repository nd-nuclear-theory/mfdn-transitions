
module OBobservables
  use Wigner3J, only: Wig3J
  implicit none
  private
  public CalcOBobs, CalcFGTobs
  !
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine CalcOBobs(Par, mxK, nobdme, robdme, obdmeOrb_offset, &
                       obdmeOrbKbraket, totR20, totM1L, totM1S, totE2Q)
    !
    use SPbasis, only: n_orb, l_orb, j2_orb
    integer, intent(in) :: Par, mxK, nobdme
    real(kind=8), dimension(nobdme), intent(in) :: robdme
    integer, dimension(0:mxK + 1), intent(in) :: obdmeOrb_offset
    integer, dimension(2, nobdme), intent(in) :: obdmeOrbKbraket
    real(kind=8), intent(out) :: totR20, totM1L, totM1S, totE2Q
    !
    ! local variables
    integer :: ia, na, la, j2a, ic, nc, lc, j2c, i
    real(kind=8) :: xR20op, xM1Lop, xM1Sop, xE2Qop
    !
    totR20 = 0.d0
    totM1L = 0.d0
    totM1S = 0.d0
    totE2Q = 0.d0
    !
    !$omp parallel default(shared)                               &
    !$omp          private(ia, na, la, j2a, ic, nc, lc, j2c, i,  &
    !$omp                  xR20op, xM1Lop, xM1Sop, xE2Qop)       &
    !$omp          reduction(+: totR20, totM1L, totM1S, totE2Q)
    !
    if (Par .eq. 1) then
      !$omp do
      do i = obdmeOrb_offset(0) + 1, obdmeOrb_offset(1)
        ia = obdmeOrbKbraket(1, i)
        na = n_orb(ia)
        la = l_orb(ia)
        j2a = j2_orb(ia)
        !
        ic = obdmeOrbKbraket(2, i)
        nc = n_orb(ic)
        lc = l_orb(ic)
        j2c = j2_orb(ic)
        !
        xR20op = R2OBoperator(na, la, j2a, nc, lc, j2c)
        totR20 = totR20 + robdme(i)*xR20op
      end do
      !$omp end do nowait
    end if
    !
    if (mxK .ge. 1) then
      !$omp do
      do i = obdmeOrb_offset(1) + 1, obdmeOrb_offset(2)
        ia = obdmeOrbKbraket(1, i)
        na = n_orb(ia)
        la = l_orb(ia)
        j2a = j2_orb(ia)
        !
        ic = obdmeOrbKbraket(2, i)
        nc = n_orb(ic)
        lc = l_orb(ic)
        j2c = j2_orb(ic)
        !
        if (Par .eq. 1) then
          call M1operators(na, la, j2a, nc, lc, j2c, xM1Lop, xM1Sop)
          totM1L = totM1L + robdme(i)*xM1Lop
          totM1S = totM1S + robdme(i)*xM1Sop
        else
          xE2Qop = E1operator(na, la, j2a, nc, lc, j2c)
          totE2Q = totE2Q + robdme(i)*xE2Qop
        end if
      end do
      !$omp end do nowait
    end if
    !
    if (mxK .ge. 2) then
      !$omp do
      do i = obdmeOrb_offset(2) + 1, obdmeOrb_offset(3)
        ia = obdmeOrbKbraket(1, i)
        na = n_orb(ia)
        la = l_orb(ia)
        j2a = j2_orb(ia)
        !
        ic = obdmeOrbKbraket(2, i)
        nc = n_orb(ic)
        lc = l_orb(ic)
        j2c = j2_orb(ic)
        !
        if (Par .eq. 1) then
          xE2Qop = E2operator(na, la, j2a, nc, lc, j2c)
          totE2Q = totE2Q + robdme(i)*xE2Qop
        else
          call M2operators(na, la, j2a, nc, lc, j2c, xM1Lop, xM1Sop)
          totM1L = totM1L + robdme(i)*xM1Lop
          totM1S = totM1S + robdme(i)*xM1Sop
        end if
      end do
      !$omp end do
    end if
    !$omp end parallel
    !
    return
  end subroutine CalcOBobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine CalcFGTobs(Par, mxK, nobdme, robdme, obdmeOrb_offset, &
                        obdmeOrbKbraket, totF0, totGT)
    !
    use SPbasis, only: n_orb, l_orb, j2_orb
    integer, intent(in) :: Par, mxK, nobdme
    real(kind=8), dimension(nobdme), intent(in) :: robdme
    integer, dimension(0:mxK + 1), intent(in) :: obdmeOrb_offset
    integer, dimension(2, nobdme), intent(in) :: obdmeOrbKbraket
    real(kind=8), intent(out) :: totF0, totGT
    !
    ! local variables
    integer :: ia, na, la, j2a, ic, nc, lc, j2c, i
    real(kind=8) :: xF0op, xGTop
    !
    totF0 = 0.d0
    totGT = 0.d0
    !
    if (Par .eq. -1) return
    !
    !$omp parallel default(shared)                                             &
    !$omp          private(ia, na, la, j2a, ic, nc, lc, j2c, i,  xF0op, xGTop) &
    !$omp          reduction(+: totF0, totGT)
    !
    !$omp do
    do i = obdmeOrb_offset(0) + 1, obdmeOrb_offset(1)
      ia = obdmeOrbKbraket(1, i)
      na = n_orb(ia)
      la = l_orb(ia)
      j2a = j2_orb(ia)
      !
      ic = obdmeOrbKbraket(2, i)
      nc = n_orb(ic)
      lc = l_orb(ic)
      j2c = j2_orb(ic)
      !
      xF0op = F0operator(na, la, j2a, nc, lc, j2c)
      totF0 = totF0 + robdme(i)*xF0op
    end do
    !$omp end do nowait
    !
    !$omp do
    do i = obdmeOrb_offset(1) + 1, obdmeOrb_offset(2)
      ia = obdmeOrbKbraket(1, i)
      na = n_orb(ia)
      la = l_orb(ia)
      j2a = j2_orb(ia)
      !
      ic = obdmeOrbKbraket(2, i)
      nc = n_orb(ic)
      lc = l_orb(ic)
      j2c = j2_orb(ic)
      !
      xGTop = GToperator(na, la, j2a, nc, lc, j2c)
      totGT = totGT + robdme(i)*xGTop
      !
    end do
    !$omp end do
    !
    !$omp end parallel
    !
    return
  end subroutine CalcFGTobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=8) function R2OBoperator(nna, lla, j2a, nnc, llc, j2c)
    integer, intent(in) :: nna, lla, j2a, nnc, llc, j2c
    !
    real(kind=8) :: r2me
    integer :: nn
    !
    if ((lla .eq. llc) .and. (j2a .eq. j2c)) then
      if (nna .eq. nnc) then
        r2me = (2*nna + lla + 1.5d0)
      elseif (abs(nna - nnc) .eq. 1) then
        nn = max(nna, nnc)
        r2me = -sqrt(nn*(nn + lla + 0.5d0))
      else
        r2me = 0.d0
      end if
      R2OBoperator = r2me*sqrt(j2a + 1.d0)
    else
      R2OBoperator = 0.d0
    end if
    !
    return
  end function R2OBoperator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=8) function E2operator(nna, lla, j2a, nnc, llc, j2c)
    integer, intent(in) :: nna, lla, j2a, nnc, llc, j2c
    !
    integer :: nn
    real(kind=8) :: r2me
    real, external :: myphase
    !
    if ((abs(j2a - j2c) .le. 4) .and. ((j2a + j2c) .ge. 4)) then
      !
      r2me = 0.d0
      !
      if (lla .eq. llc) then
        if (nna .eq. nnc) then
          r2me = (2*nna + lla + 1.5d0)
        elseif (abs(nna - nnc) .eq. 1) then
          nn = max(nna, nnc)
          r2me = -sqrt(nn*(nn + lla + 0.5d0))
        end if
        !
      elseif ((lla - llc) .eq. 2) then
        if (nna .eq. nnc) then
          r2me = sqrt((nnc + lla + 0.5d0)*(nnc + lla - 0.5d0))
        elseif ((nna - nnc) .eq. -1) then
          r2me = -2*sqrt(nnc*(nnc + lla - 0.5d0))
        elseif ((nna - nnc) .eq. -2) then
          r2me = sqrt(nnc*(nnc - 1.0d0))
        end if
        !
      elseif ((lla - llc) .eq. -2) then
        if (nna .eq. nnc) then
          r2me = sqrt((nna + llc + 0.5d0)*(nna + llc - 0.5d0))
        elseif ((nna - nnc) .eq. 1) then
          r2me = -2*sqrt(nna*(nna + llc - 0.5d0))
        elseif ((nna - nnc) .eq. 2) then
          r2me = sqrt(nna*(nna - 1.0d0))
        end if
      end if
      !
      if (r2me .ne. 0.d0) then
        E2operator = r2me*sqrt((j2a + 1.d0)*(j2c + 1.d0)) &
                     *Wig3J(j2a, 4, j2c, 1, 0, -1)*myphase((j2a + 1)/2)
      else
        E2operator = 0.d0
      end if
      !
    else
      E2operator = 0.d0
    end if
    !
    return
  end function E2operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine M1operators(nna, lla, j2a, nnc, llc, j2c, M1Lop, M1Sop)
    integer, intent(in) :: nna, lla, j2a, nnc, llc, j2c
    real(kind=8), intent(out) :: M1Lop, M1Sop
    !
    real(kind=8) :: factor, kappa
    real, external :: myphase
    !
    if ((nna .eq. nnc) .and. (lla .eq. llc)) then
      factor = sqrt((j2a + 1.d0)*(j2c + 1.d0)) &
               *Wig3J(j2a, 2, j2c, 1, 0, -1)*myphase((j2a + 1)/2)
      kappa = 0.5d0*(j2a + 1.d0)*myphase(lla + (j2a + 1)/2) &
              + 0.5d0*(j2c + 1.d0)*myphase(llc + (j2c + 1)/2)
      M1Lop = (1.d0 + 0.5d0*kappa)*(1.d0 - kappa)*factor
      M1Sop = -0.5d0*(1.d0 - kappa)*factor
    else
      M1Lop = 0.d0
      M1Sop = 0.d0
    end if
    !
    return
  end subroutine M1operators

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=8) function E1operator(nna, lla, j2a, nnc, llc, j2c)
    integer, intent(in) :: nna, lla, j2a, nnc, llc, j2c
    !
    real(kind=8) :: r1me
    real, external :: myphase
    !
    ! see Suhonen (6.40) -- (6.45)
    if ((abs(j2a - j2c) .le. 2) .and. ((j2a + j2c) .ge. 2)) then
      !
      r1me = 0.d0
      !
      if (lla .eq. llc + 1) then
        if (nna .eq. nnc) then
          r1me = sqrt(lla + nna + 0.5d0)
        elseif (nna .eq. nnc - 1) then
          r1me = -sqrt(nna + 1.0d0)
        end if
      elseif (llc .eq. lla + 1) then
        if (nnc .eq. nna) then
          r1me = sqrt(llc + nnc + 0.5d0)
        elseif (nnc .eq. nna - 1) then
          r1me = -sqrt(nnc + 1.0d0)
        end if
      end if
      !
      if (r1me .ne. 0.d0) then
        E1operator = r1me*sqrt((j2a + 1.d0)*(j2c + 1.d0)) &
                     *Wig3J(j2a, 2, j2c, 1, 0, -1)*myphase((j2a + 1)/2)
      else
        E1operator = 0.d0
      end if
      !
    else
      E1operator = 0.d0
    end if
    !
    return
  end function E1operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine M2operators(nna, lla, j2a, nnc, llc, j2c, M2Lop, M2Sop)
    integer, intent(in) :: nna, lla, j2a, nnc, llc, j2c
    real(kind=8), intent(out) :: M2Lop, M2Sop
    !
    real(kind=8) :: r2me, factor, kappa
    real, external :: myphase
    !
    ! see Suhonen (6.40) -- (6.45)
    if ((abs(j2a - j2c) .le. 4) .and. ((j2a + j2c) .ge. 4)) then
      !
      r2me = 0.d0
      if (lla .eq. llc + 1) then
        if (nna .eq. nnc) then
          r2me = sqrt(lla + nna + 0.5d0)
        elseif (nna .eq. nnc - 1) then
          r2me = -sqrt(nna + 1.0d0)
        end if
      elseif (llc .eq. lla + 1) then
        if (nnc .eq. nna) then
          r2me = sqrt(llc + nnc + 0.5d0)
        elseif (nnc .eq. nna - 1) then
          r2me = -sqrt(nnc + 1.0d0)
        end if
      end if
      !
      if (r2me .ne. 0.d0) then
        factor = r2me*sqrt((j2a + 1.d0)*(j2c + 1.d0)) &
                 *Wig3J(j2a, 4, j2c, 1, 0, -1)*myphase((j2a + 1)/2)
        kappa = 0.5d0*(j2a + 1.d0)*myphase(lla + (j2a + 1)/2) &
                + 0.5d0*(j2c + 1.d0)*myphase(llc + (j2c + 1)/2)
        M2Lop = (1.d0 + kappa/3.d0)*(2.d0 - kappa)*factor
        M2Sop = -0.5d0*(2.d0 - kappa)*factor
      else
        M2Lop = 0.d0
        M2Sop = 0.d0
      end if
      !
    else
      M2Lop = 0.d0
      M2Sop = 0.d0
    end if
    !
    return
  end subroutine M2operators

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=8) function F0operator(nna, lla, j2a, nnc, llc, j2c)
    integer, intent(in) :: nna, lla, j2a, nnc, llc, j2c
    !
    real, external :: myphase
    !
    if ((nna .eq. nnc) .and. (lla .eq. llc) .and. (j2a .eq. j2c)) then
      F0operator = sqrt(j2a + 1.d0)
    else
      F0operator = 0.d0
    end if
    !
    return
  end function F0operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=8) function GToperator(nna, lla, j2a, nnc, llc, j2c)
    integer, intent(in) :: nna, lla, j2a, nnc, llc, j2c
    !
    real(kind=8) :: factor, kappa
    real, external :: myphase
    !
    if ((nna .eq. nnc) .and. (lla .eq. llc)) then
      factor = sqrt((j2a + 1.d0)*(j2c + 1.d0)/3.d0) &
               *Wig3J(j2a, 2, j2c, 1, 0, -1)*myphase((j2a + 1)/2)
      kappa = 0.5d0*(j2a + 1.d0)*myphase(lla + (j2a + 1)/2) &
              + 0.5d0*(j2c + 1.d0)*myphase(llc + (j2c + 1)/2)
      ! M1Lop = (1.d0 + 0.5d0*kappa) * (1.d0 - kappa) * factor
      GToperator = -1.0d0*(1.d0 - kappa)*factor
    else
      ! M1Lop = 0.d0
      GToperator = 0.d0
    end if
    !
    return
  end function GToperator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module OBobservables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
