
module TBops
  use SPbasis, only: nparticles, norb_p, norb_n, mj2_sp, orb_sp, pr_orb, j2_orb
  use TBME, only: TBMEfull_pp, TBMEfull_nn, TBMEfull_pn, &
       ntbme_pp, ntbme_nn, ntbme_pn, &
       J2max_pp, J2max_nn, J2max_pn, &
       ntpsJ_pp, ntpsJ_nn, ntpsJ_pn, ntbmeJ_pp, ntbmeJ_nn, ntbmeJ_pn, &
       tpsJindx_pp, tpsJindx_nn, tpsJindx_pn
  implicit none
  private
  public MBtwobodyObs, MBtwobodyObsDiag
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
contains
  !
  subroutine MBtwobodyObs(rowstate, colstate, ndiffs, &
       rowdifloc, coldifloc, nTBops, xTBops)
    !
    ! SPECIFICALLY for SCALAR 2-BODY OPERATORS preserving Tz
    !
    integer, intent(in) :: nTBops, ndiffs
    integer(kind=2), dimension(nparticles), intent(in) :: rowstate, colstate
    integer, dimension(2), intent(in) :: rowdifloc, coldifloc
    real(kind=8), dimension(nTBops), intent(out) :: xTBops
    !
    ! local variables
    integer :: aa, bb, cc, dd, kk, i
    integer :: m2a, orba, pra, j2a
    integer :: m2b, orbb, prb, j2b
    integer :: m2c, orbc, prc, j2c
    integer :: m2d, orbd, prd, j2d
    integer :: m2k, orbk, prk, j2k
    real(kind=8), dimension(nTBops) :: xtmp
    real :: phase
    real, external :: myphase
    !
    if (ndiffs .eq. 0) then
       !
       call MBtwobodyObsDiag(rowstate, nTBops, xTBops)
       !
    elseif (ndiffs .eq. 1) then
       !
       aa = rowstate(rowdifloc(1))
       m2a  = mj2_sp(aa)
       orba = orb_sp(aa)
       pra  = pr_orb(orba)
       j2a  = j2_orb(orba)
       !
       cc = colstate(coldifloc(1))
       m2c  = mj2_sp(cc)
       orbc = orb_sp(cc)
       prc  = pr_orb(orbc)
       j2c  = j2_orb(orbc)
       !
       ! sum over one spectator
       ! exclude SP state based on column state
       ! (or equivalently, can be done based on row state)
       xTBops(1:nTBops) = 0.d0
       do i = 1, coldifloc(1) - 1
          kk = colstate(i)
          m2k  = mj2_sp(kk)
          orbk = orb_sp(kk)
          prk  = pr_orb(orbk)
          j2k  = j2_orb(orbk)
          !
          if (aa .le. kk) then
             call TBMEops(nTBops, xtmp,                     &
                  orba, pra, j2a, m2a, orbk, prk, j2k, m2k, &
                  orbk, prk, j2k, m2k, orbc, prc, j2c, m2c )
             xTBops(1:nTBops) = xTBops(1:nTBops) - xtmp(1:nTBops)
          else
             call TBMEops(nTBops, xtmp,                     &
                  orbk, prk, j2k, m2k, orba, pra, j2a, m2a, &
                  orbk, prk, j2k, m2k, orbc, prc, j2c, m2c )
             xTBops(1:nTBops) = xTBops(1:nTBops) + xtmp(1:nTBops)
          endif
       enddo
       do i = coldifloc(1) + 1, nparticles
          kk = colstate(i)
          m2k  = mj2_sp(kk)
          orbk = orb_sp(kk)
          prk  = pr_orb(orbk)
          j2k  = j2_orb(orbk)
          !
          if (aa .le. kk) then
             call TBMEops(nTBops, xtmp,                     &
                  orba, pra, j2a, m2a, orbk, prk, j2k, m2k, &
                  orbc, prc, j2c, m2c, orbk, prk, j2k, m2k )
             xTBops(1:nTBops) = xTBops(1:nTBops) + xtmp(1:nTBops)
          else
             call TBMEops(nTBops, xtmp,                     &
                  orbk, prk, j2k, m2k, orba, pra, j2a, m2a, &
                  orbc, prc, j2c, m2c, orbk, prk, j2k, m2k )
             xTBops(1:nTBops) = xTBops(1:nTBops) - xtmp(1:nTBops)
          endif
       enddo
       phase = myphase(rowdifloc(1) - coldifloc(1))
       xTBops(1:nTBops) = phase * xTBops(1:nTBops)
       !
    else
       !
       aa = rowstate(rowdifloc(1))
       bb = rowstate(rowdifloc(2))
       cc = colstate(coldifloc(1))
       dd = colstate(coldifloc(2))
       !
       orba = orb_sp(aa)
       orbb = orb_sp(bb)
       orbc = orb_sp(cc)
       orbd = orb_sp(dd)
       !
       pra  = pr_orb(orba)
       prb  = pr_orb(orbb)
       prc  = pr_orb(orbc)
       prd  = pr_orb(orbd)
       !
       j2a  = j2_orb(orba)
       j2b  = j2_orb(orbb)
       j2c  = j2_orb(orbc)
       j2d  = j2_orb(orbd)
       !
       m2a  = mj2_sp(aa)
       m2b  = mj2_sp(bb)
       m2c  = mj2_sp(cc)
       m2d  = mj2_sp(dd)
       !     
       call TBMEops(nTBops, xTBops, &
            orba, pra, j2a, m2a, orbb, prb, j2b, m2b, &
            orbc, prc, j2c, m2c, orbd, prd, j2d, m2d )
       !     
       phase = myphase(rowdifloc(1) + rowdifloc(2) &
            - coldifloc(1) - coldifloc(2) )
       xTBops(1:nTBops) = phase * xTBops(1:nTBops)
       !
    endif
    !
    return
  end subroutine MBtwobodyObs
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine MBtwobodyObsDiag(state, nTBops, xTBops)
    !
    ! SPECIFICALLY for SCALAR 2-BODY OPERATORS preserving Tz
    !
    integer, intent(in) :: nTBops
    integer(kind=2), dimension(nparticles), intent(in) :: state
    real(kind=8), dimension(nTBops), intent(out) :: xTBops
    !
    ! local variables
    integer :: aa, bb
    integer :: i, j
    integer :: m2a, orba, pra, j2a
    integer :: m2b, orbb, prb, j2b
    real(kind=8), dimension(nTBops) :: xtmp
    !
    xTBops(1:nTBops) = 0.d0
    !
    do i = 1, nparticles-1
       aa = state(i)
       m2a  = mj2_sp(aa)
       orba = orb_sp(aa)
       pra  = pr_orb(orba)
       j2a  = j2_orb(orba)
       do j = i+1, nparticles
          bb = state(j)
          m2b  = mj2_sp(bb)
          orbb = orb_sp(bb)
          prb  = pr_orb(orbb)
          j2b  = j2_orb(orbb)
          call TBMEops(nTBops, xtmp, &
               orba, pra, j2a, m2a, orbb, prb, j2b, m2b, &
               orba, pra, j2a, m2a, orbb, prb, j2b, m2b )
          xTBops(1:nTBops) = xTBops(1:nTBops) + xtmp(1:nTBops)
       enddo
    enddo
    !
    return
  end subroutine MBtwobodyObsDiag
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine TBMEops(nops, xops, &
       orba, pra, j2a, m2a, orbb, prb, j2b, m2b, &
       orbc, prc, j2c, m2c, orbd, prd, j2d, m2d )
    !
    !     single-particle indiced aa < bb ; cc < dd
    !
    integer, intent(in) :: nops
    real(kind=8), dimension(nops), intent(out) :: xops
    integer, intent(in) :: orba, orbb, orbc, orbd
    integer, intent(in) :: pra, j2a, m2a, prb, j2b, m2b
    integer, intent(in) :: prc, j2c, m2c, prd, j2d, m2d
    !
    ! local variables
    integer :: mmjj, jmin, jmax, tz2a, tz2b, pt
    real, external :: myphase
    real :: fac
    !
    mmjj = m2a + m2b
    ! Test on J (scalar operators only for now)
    jmin = max(abs(j2a-j2b), abs(j2c-j2d), abs(mmjj)) / 2
    jmax = min( (j2a+j2b), (j2c+j2d)) / 2
    if (jmin .gt. jmax) then
       xops(1:nops) = 0.d0
       return
    endif
    !
    pt = (1-pra*prb)/2
    if (orbb .le. norb_p) then
       call TBops_IDN(orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d,     &
            m2a,m2b,m2c,m2d,  pt, jmin, jmax, J2max_pp, norb_p, &
            ntbme_pp, tpsJindx_pp, ntpsJ_pp, ntbmeJ_pp,         &
            TBMEfull_pp, nops, xops)
       tz2a = 1
       tz2b = 1
    elseif (orba .gt. norb_p) then
       call TBops_IDN(orba-norb_p, orbb-norb_p,                 &
            orbc-norb_p, orbd-norb_p, j2a,j2b,j2c,j2d,          &
            m2a,m2b,m2c,m2d,  pt, jmin, jmax, J2max_nn, norb_n, &
            ntbme_nn, tpsJindx_nn, ntpsJ_nn, ntbmeJ_nn,         &
            TBMEfull_nn, nops, xops)                 
       tz2a = -1
       tz2b = -1
    else
       call TBops_DIS(orba, orbb-norb_p,                        &
            orbc, orbd-norb_p, j2a,j2b,j2c,j2d,                 &
            m2a,m2b,m2c,m2d,  pt, jmin, jmax, J2max_pn, norb_p,norb_n, &
            ntbme_pn, tpsJindx_pn, ntpsJ_pn, ntbmeJ_pn,         &
            TBMEfull_pn, nops, xops)
       tz2a = 1
       tz2b = -1
    endif
    !
    ! additional overall phase factors
    ! (due to conventions for Clebsch-Gordan vs. Wigner3J )
    fac = myphase((j2a-j2b+j2c-j2d)/2)
    xops(1:nops) = fac * xops(1:nops)
    !
    return
  end subroutine TBMEops
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine TBops_IDN(orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d, &
       m2a,m2b,m2c,m2d, pt, jmin, jmax, Jtot2max, norb, &
       ntbme, tpsJindx, ntpsJ, ntbmeJ, TBMEfull, nops, xops)
    use Wigner3J, only: Retrieve_Wigner3J_Array
    !
    integer, intent(in) :: orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d, &
         m2a,m2b,m2c,m2d, pt, jmin, jmax, Jtot2max, norb, ntbme, nops
    integer, dimension(0:Jtot2max, norb*(norb+1)/2) :: tpsJindx
    integer, dimension(0:Jtot2max+1, 0:1) :: ntpsJ, ntbmeJ
    real(kind=4), dimension(nops, ntbme), intent(in) :: TBMEfull
    real(kind=8), dimension(nops), intent(out) :: xops
    ! local variables
    integer :: jt, nj12, nj34, indx
    real(kind=8) :: factor
    real(kind=8), dimension(Jtot2max+1) :: wigab, wigcd
    real(kind=8), parameter :: sqr2 = sqrt(2.d0)
    logical, external :: pairwiseless
    !
    call Retrieve_Wigner3j_array(j2a,j2b,m2a,m2b,jmin,jmax, wigab)
    call Retrieve_Wigner3j_array(j2c,j2d,m2c,m2d,jmin,jmax, wigcd)
    !
    xops(1:nops) = 0.0d0
    do jt = jmin, jmax
       if (pairwiseless(orba, orbb, orbc, orbd)) then
          nj12 = tpsJindx(jt, (orba+orbb*(orbb-1)/2))
          nj34 = tpsJindx(jt, (orbc+orbd*(orbd-1)/2))
       else
          nj12 = tpsJindx(jt, (orbc+orbd*(orbd-1)/2))
          nj34 = tpsJindx(jt, (orba+orbb*(orbb-1)/2))
       endif
       if (nj12*nj34 .eq. 0) cycle
       !
       indx = ntbmeJ(jt,pt) + ntpsJ(jt,pt) * (nj12-1) + nj34 - nj12 * (nj12-1) / 2
       !
       factor = wigab(jt-jmin+1) * wigcd(jt-jmin+1) * (2.0*jt+1.0)
       xops(1:nops) = xops(1:nops) + factor * TBMEfull(1:nops, indx)
    enddo
    ! Needed for identical particles
    if (orba.eq.orbb) then
       xops(1:nops) = xops(1:nops) * sqr2
    endif
    if (orbc.eq.orbd) then
       xops(1:nops) = xops(1:nops) * sqr2
    endif
    !
    return
  end subroutine TBops_IDN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine TBops_DIS(orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d, &
       m2a,m2b,m2c,m2d, pt, jmin, jmax, Jtot2max, norb1, norb2, &
       ntbme, tpsJindx, ntpsJ, ntbmeJ, TBMEfull, nops, xops)
    use Wigner3J, only: Retrieve_Wigner3J_Array
    !
    integer, intent(in) :: orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d, &
         m2a,m2b,m2c,m2d, pt, jmin, jmax, Jtot2max, norb1, norb2, ntbme, nops
    integer, dimension(0:Jtot2max, norb1, norb2) :: tpsJindx
    integer, dimension(0:Jtot2max+1, 0:1) :: ntpsJ, ntbmeJ
    real(kind=4), dimension(nops, ntbme), intent(in) :: TBMEfull
    real(kind=8), dimension(nops), intent(out) :: xops
    ! local variables
    integer :: jt, nj12, nj34, indx
    real(kind=8) :: factor
    real(kind=8), dimension(Jtot2max+1) :: wigab, wigcd
    logical, external :: pairwiseless
    !
    call Retrieve_Wigner3j_array(j2a,j2b,m2a,m2b,jmin,jmax, wigab)
    call Retrieve_Wigner3j_array(j2c,j2d,m2c,m2d,jmin,jmax, wigcd)
    !
    xops(1:nops) = 0.0d0
    do jt = jmin, jmax
       if (pairwiseless(orba, orbb, orbc, orbd)) then
          nj12 = tpsJindx(jt, orba, orbb)
          nj34 = tpsJindx(jt, orbc, orbd)
       else
          nj12 = tpsJindx(jt, orbc, orbd)
          nj34 = tpsJindx(jt, orba, orbb)
       endif
       if (nj12*nj34 .eq. 0) cycle
       !
       indx = ntbmeJ(jt,pt) + ntpsJ(jt,pt) * (nj12-1) + nj34 - nj12 * (nj12-1) / 2
       !
       factor = wigab(jt-jmin+1) * wigcd(jt-jmin+1) * (2.0*jt+1.0)
       xops(1:nops) = xops(1:nops) + factor * TBMEfull(1:nops, indx)
    enddo
    !
    return
  end subroutine TBops_DIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module TBops
