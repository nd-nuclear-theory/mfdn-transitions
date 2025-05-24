
module TBops
  use SPbasis, only: nparticles, norb_p, norb_n, mj2_sp, orb_sp, j2_orb, pr_orb
#ifdef DeltaTz
  use TBME_Tz12, only: Jop, Parop
  use TBME_Tz12, only: nTBMEs_max, Set_TBME_array
#else
  use TBME_Tz0, only: Jop, Parop
  use TBME_Tz0, only: nTBMEs_max, Set_TBME_array
#endif
  implicit none
  private
  public MBtwobodyObs, TBops_Eval, TBops_EvalDiag
  !
  ! subroutine MBtwobodyObs
  ! subroutine MBtwobodyObsDiag
  ! subroutine TBMEops
  ! subroutine TBops_Eval
  ! subroutine TBops_EvalDiag
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
contains
  !
  subroutine MBtwobodyObs(rowstate, colstate, ndiffs, &
       rowdifloc, coldifloc, nTBops, xTBops)
    !
    integer, intent(in) :: nTBops, ndiffs
    integer(kind=2), dimension(nparticles), intent(in) :: rowstate, colstate
    integer, dimension(2), intent(in) :: rowdifloc, coldifloc
    real(kind=8), dimension(nTBops), intent(out) :: xTBops
    !
    ! local variables
    integer :: aa, bb, cc, dd, kk, i
    integer :: m2a, orba, j2a, pra
    integer :: m2b, orbb, j2b, prb
    integer :: m2c, orbc, j2c, prc
    integer :: m2d, orbd, j2d, prd
    integer :: m2k, orbk, j2k, prk
    real(kind=8), dimension(nTBops) :: xtmp
    real :: phase
    real, external :: myphase
    !
    if (ndiffs .eq. 0) then
       !
       if (Parop .eq. 1) then
          call MBtwobodyObsDiag(rowstate, nTBops, xTBops)
       else
          xTBops(1:nTBops) = 0.d0
          return
       endif
       !
    elseif (ndiffs .eq. 1) then
       !
       aa = rowstate(rowdifloc(1))
       m2a  = mj2_sp(aa)
       orba = orb_sp(aa)
       j2a  = j2_orb(orba)
       pra  = pr_orb(orba)
       !
       cc = colstate(coldifloc(1))
       m2c  = mj2_sp(cc)
       orbc = orb_sp(cc)
       j2c  = j2_orb(orbc)
       prc  = pr_orb(orbc)
       !
       if (pra*prc .ne. Parop) then
          xTBops = 0.d0
          return
       endif
       !
       ! sum over one spectator
       ! exclude SP state based on column state
       ! (or equivalently, can be done based on row state)
       xTBops(1:nTBops) = 0.d0
       do i = 1, coldifloc(1) - 1
          kk = colstate(i)
          m2k  = mj2_sp(kk)
          orbk = orb_sp(kk)
          j2k  = j2_orb(orbk)
          prk  = pr_orb(orbk)
          !
          if (aa .le. kk) then
             call TBMEops(nTBops, xtmp, pra*prk, prk*prc, &
                  orba, j2a, m2a, orbk, j2k, m2k,         &
                  orbk, j2k, m2k, orbc, j2c, m2c)
             xTBops(1:nTBops) = xTBops(1:nTBops) - xtmp(1:nTBops)
          else
             call TBMEops(nTBops, xtmp, prk*pra, prk*prc, &
                  orbk, j2k, m2k, orba, j2a, m2a,         &
                  orbk, j2k, m2k, orbc, j2c, m2c)
             xTBops(1:nTBops) = xTBops(1:nTBops) + xtmp(1:nTBops)
          endif
       enddo
       do i = coldifloc(1) + 1, nparticles
          kk = colstate(i)
          m2k  = mj2_sp(kk)
          orbk = orb_sp(kk)
          j2k  = j2_orb(orbk)
          prk  = pr_orb(orbk)
          !
          if (aa .le. kk) then
             call TBMEops(nTBops, xtmp, pra*prk, prc*prk, &
                  orba, j2a, m2a, orbk, j2k, m2k,         &
                  orbc, j2c, m2c, orbk, j2k, m2k)
             xTBops(1:nTBops) = xTBops(1:nTBops) + xtmp(1:nTBops)
          else
             call TBMEops(nTBops, xtmp, prk*pra, prc*prk, &
                  orbk, j2k, m2k, orba, j2a, m2a,         &
                  orbc, j2c, m2c, orbk, j2k, m2k)
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
       if (pra*prb*prc*prd .ne. Parop) then
          xTBops = 0.d0
          return
       endif
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
       call TBMEops(nTBops, xTBops, pra*prb, prc*prd, &
            orba, j2a, m2a, orbb, j2b, m2b,       &
            orbc, j2c, m2c, orbd, j2d, m2d)
       !
       phase = myphase(rowdifloc(1) + rowdifloc(2) &
            - coldifloc(1) - coldifloc(2) )
       xTBops(1:nTBops) = phase * xTBops(1:nTBops)
       !
    endif
    !
    return
  end subroutine MBtwobodyObs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine MBtwobodyObsDiag(state, nTBops, xTBops)
    !
    integer, intent(in) :: nTBops
    integer(kind=2), dimension(nparticles), intent(in) :: state
    real(kind=8), dimension(nTBops), intent(out) :: xTBops
    !
    ! local variables
    logical :: abIDN
    integer :: i, j, aa, bb, parab, nTBMEs
    integer :: m2a, orba, j2a, pra
    integer :: m2b, orbb, j2b, prb
    real(kind=4), dimension(nTBops, nTBMEs_max) :: TBMEarray
    real(kind=8), dimension(nTBops) :: xtmp
    !
    xTBops(1:nTBops) = 0.d0
    !
    do i = 1, nparticles-1
       aa = state(i)
       m2a  = mj2_sp(aa)
       orba = orb_sp(aa)
       j2a  = j2_orb(orba)
       pra  = pr_orb(orba)
       do j = i+1, nparticles
          bb = state(j)
          m2b  = mj2_sp(bb)
          orbb = orb_sp(bb)
          j2b  = j2_orb(orbb)
          prb  = pr_orb(orbb)
          parab = pra * prb
          !
          nTBMEs = nTBMEs_max
          call Set_TBME_array(parab, parab, orba, orbb, orba, orbb,   &
               nTBops, nTBMEs, TBMEarray)
          !
          abIDN = orba .eq. orbb
          call TBops_EvalDiag(abIDN, abIDN,       &
               j2a,j2b,j2a,j2b, m2a,m2b,m2a,m2b,  &
               nTBMEs, TBMEarray, nTBops, xtmp)
          !
          xTBops(1:nTBops) = xTBops(1:nTBops) + xtmp(1:nTBops)
       enddo
    enddo
    !
    return
  end subroutine MBtwobodyObsDiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine TBMEops(nTBops, xops, parab, parcd, &
       orba, j2a, m2a, orbb, j2b, m2b, orbc, j2c, m2c, orbd, j2d, m2d)
    !
    !     single-orbital indices aa =< bb and cc =< dd
    !
    real(kind=8), dimension(nTBops), intent(out) :: xops
    integer, intent(in) :: nTBops, parab, parcd, &
         orba, j2a, m2a, orbb, j2b, m2b, orbc, j2c, m2c, orbd, j2d, m2d
    !
    ! local variables
    logical :: abIDN, cdIDN
    integer :: nTBMEs, twomj
    real(kind=4), dimension(nTBops, nTBMEs_max) :: TBMEarray
    real, external :: myphase
    logical, external :: pairwiseless
    !
    nTBMEs = nTBMEs_max
    !
    abIDN = orba .eq. orbb
    cdIDN = orbc .eq. orbd
    !
    !if (Tzop .ne. 0) then
    !   ! always ordered (pp2pn, pn2nn, or pp2nn)
    !   ! DIAG = .false.       ! no diagonal tiles
    !   call Set_TBME_array(parab, parcd, orba, orbb, orbc, orbd,    &
    !        nTBops, nTBMEs, TBMEarray)
    !   call TBops_EvalDiag(abIDN, cdIDN,                            &
    !        j2a,j2b,j2c,j2d, m2a,m2b,m2c,m2d, twomj, m2ab, m2cd,    &
    !        nTBMEs, TBMEarray, nTBops, xops)
    !   return
    !endif
    !
    if ((orba.eq.orbc) .and. (orbb.eq.orbd)) then
       call Set_TBME_array(parab, parcd, orba, orbb, orbc, orbd,    &
            nTBops, nTBMEs, TBMEarray)
       call TBops_EvalDiag(abIDN, cdIDN,       &
            j2a,j2b,j2c,j2d, m2a,m2b,m2c,m2d,  &
            nTBMEs, TBMEarray, nTBops, xops)
    elseif (pairwiseless(orba,orbb,orbc,orbd)) then
       call Set_TBME_array(parab, parcd, orba, orbb, orbc, orbd,    &
            nTBops, nTBMEs, TBMEarray)
       call TBops_Eval(abIDN, cdIDN,           &
            j2a,j2b,j2c,j2d, m2a,m2b,m2c,m2d,  &
            nTBMEs, TBMEarray, nTBops, xops)
    else
       call Set_TBME_array(parcd, parab, orbc, orbd, orba, orbb,    &
            nTBops, nTBMEs, TBMEarray)
       call TBops_Eval(cdIDN, abIDN,           &
            j2c,j2d,j2a,j2b, m2c,m2d,m2a,m2b,  &
            nTBMEs, TBMEarray, nTBops, xops)
       twomj = m2a + m2b - m2c - m2d
       xops(1:nTBops) = myphase(twomj/2) * xops(1:nTBops)
    endif
    !
    return
  end subroutine TBMEops

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine TBops_Eval(abIDN, cdIDN,          &
       j2a,j2b,j2c,j2d, m2a,m2b,m2c,m2d,       &
       nTBMEs, TBMEarray, nTBops, xops )
    use Wigner3J, only: Retrieve_Wigner3J_Array, Wig3J_Array_JJK, Jmax_2body
    use Wigner3J, only: Wig3J
    !
    !     single-orbital indices aa =< bb and cc =< dd ; (aa, bb) =< (cc, dd)
    !
    logical, intent(in) :: abIDN, cdIDN
    integer, intent(in) :: j2a,j2b,j2c,j2d, m2a, m2b, m2c, m2d, nTBMEs, nTBops
    real(kind=4), dimension(nTBops, nTBMEs), intent(in) :: TBMEarray
    real(kind=8), dimension(nTBops), intent(out) :: xops
    ! local variables
    real(kind=8), parameter :: sqr2 = sqrt(2.d0)
    integer :: jabmin, jabmax, jabstep, jcdmin, jcdmax, jcdstep,    &
          iab, icd, icd_start, icd_stop, indx, deltamj, mab, mcd, idj
    real(kind=8) :: facab, faccd, facJop
    real(kind=8), dimension(Jmax_2body+1) :: wigab, wigcd
    real, external :: myphase
    !
    xops(1:nTBops) = 0.0d0
    !
    mab = (m2a + m2b)/2
    mcd = (m2c + m2d)/2
    deltamj = mab - mcd
    !
    jabmin = abs(j2a-j2b)/2
    jabmax = (j2a+j2b)/2
    jabstep = 1
    if (abIDN) then
       jabmin = jabmin + mod(jabmin,2)
       jabmax = jabmax - mod(jabmax,2)
       jabstep = 2
    endif
    call Retrieve_Wigner3j_array(j2a,j2b,m2a,m2b,jabmin,jabmax, wigab)
    !
    jcdmin = abs(j2c-j2d)/2
    jcdmax = (j2c+j2d)/2
    jcdstep = 1
    if (cdIDN) then
       jcdmin = jcdmin + mod(jcdmin,2)
       jcdmax = jcdmax - mod(jcdmax,2)
       jcdstep = 2
    endif
    call Retrieve_Wigner3j_array(j2c,j2d,m2c,m2d,jcdmin,jcdmax, wigcd)
    !
    indx = 0
    do iab = jabmin, jabmax, jabstep
       facab = wigab(iab-jabmin+1) * sqrt(2.d0*iab + 1.d0) * myphase(iab)
       !
       icd_start = max(abs(iab-Jop), jcdmin)
       icd_stop  = min(iab+Jop, jcdmax)
       if (cdIDN) then
          icd_start = icd_start + mod(icd_start,2)
          icd_stop = icd_stop - mod(icd_stop,2)
       endif
       !
       do icd = icd_start, icd_stop, jcdstep
          indx = indx + 1
          !if (indx .gt. nTBMEs) then
          !   print*, 'error', indx, nTBMEs, jabstep, jcdstep
          !   print*, 'error', j2a,m2a,j2b,m2b, iab
          !   print*, 'error', j2c,m2c,j2d,m2d, icd
          !   call cancelall(999)
          !endif
          faccd = wigcd(icd-jcdmin+1) * sqrt(2.d0*icd + 1.d0)
          !
          idj = icd - iab
          if (deltamj .ge. 0) then
             facJop = Wig3J_Array_JJK(idj, iab, mab)
          else
             facJop = Wig3J_Array_JJK(idj, iab, -mab) *    &
                  myphase(iab + Jop + icd)
          endif
          !
          ! vectorizes
          xops(1:nTBops) = xops(1:nTBops) + &
               facab * faccd * facJop * TBMEarray(1:nTBops, indx)
       end do
    end do
    !
    if (abIDN) then
       xops(1:nTBops) = xops(1:nTBops) * sqr2
    endif
    if (cdIDN) then
       xops(1:nTBops) = xops(1:nTBops) * sqr2
    endif
    !
    xops(1:nTBops) = myphase((j2a-j2b+j2c-j2d)/2 - mcd) * xops(1:nTBops)
    return
    !
  end subroutine TBops_Eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine TBops_EvalDiag(abIDN, cdIDN,      &
       j2a,j2b,j2c,j2d, m2a,m2b,m2c,m2d,       &
       nTBMEs, TBMEarray, nTBops, xops )
    use Wigner3J, only: Retrieve_Wigner3J_Array, Wig3J_Array_JJK, Jmax_2body
    use Wigner3J, only: Wig3J
    !
    !     single-orbital indices aa =< bb and cc =< dd ; (aa, bb) =< (cc, dd)
    !
    logical, intent(in) :: abIDN, cdIDN
    integer, intent(in) :: j2a,j2b,j2c,j2d, m2a, m2b, m2c, m2d, nTBMEs, nTBops
    real(kind=4), dimension(nTBops, nTBMEs), intent(in) :: TBMEarray
    real(kind=8), dimension(nTBops), intent(out) :: xops
    ! local variables
    real(kind=8), parameter :: sqr2 = sqrt(2.d0)
    integer :: jabmin, jabmax, jabstep, jcdmin, jcdmax, jcdstep,    &
          iab, icd, icd_start, icd_stop, indx, deltamj, mab, mcd, idj
    real(kind=8) :: facab, faccd, facJop, facabp
    real(kind=8), dimension(Jmax_2body+1) :: wigab, wigcd
    real, external :: myphase
    !
    xops(1:nTBops) = 0.0d0
    !
    mab = (m2a + m2b)/2
    mcd = (m2c + m2d)/2
    deltamj = mab - mcd
    !
    jabmin = abs(j2a-j2b)/2
    jabmax = (j2a+j2b)/2
    jabstep = 1
    if (abIDN) then
       jabmin = jabmin + mod(jabmin,2)
       jabmax = jabmax - mod(jabmax,2)
       jabstep = 2
    endif
    call Retrieve_Wigner3j_array(j2a,j2b,m2a,m2b,jabmin,jabmax, wigab)
    !
    jcdmin = abs(j2c-j2d)/2
    jcdmax = (j2c+j2d)/2
    jcdstep = 1
    if (cdIDN) then
       jcdmin = jcdmin + mod(jcdmin,2)
       jcdmax = jcdmax - mod(jcdmax,2)
       jcdstep = 2
    endif
    call Retrieve_Wigner3j_array(j2c,j2d,m2c,m2d,jcdmin,jcdmax, wigcd)
    !
    indx = 0
    do iab = jabmin, jabmax, jabstep
       facab = wigab(iab-jabmin+1) * sqrt(2.d0*iab + 1.d0) * myphase(iab)
       !
       icd_start = max(abs(iab-Jop), iab)  ! specific for DIAG
       icd_stop  = min(iab+Jop, jcdmax)
       if (cdIDN) then
          icd_start = icd_start + mod(icd_start,2)
          icd_stop = icd_stop - mod(icd_stop,2)
       endif
       !
       do icd = icd_start, icd_stop, jcdstep
          indx = indx + 1
          faccd = wigcd(icd-jcdmin+1) * sqrt(2.d0*icd + 1.d0)
          !
          idj = icd - iab
          if (deltamj .ge. 0) then
             facJop = Wig3J_Array_JJK(idj, iab, mab)
          else
             facJop = Wig3J_Array_JJK(idj, iab, -mab) *    &
                  myphase(iab + Jop + icd)
          endif
          !
          ! vectorizes
          xops(1:nTBops) = xops(1:nTBops) + &
               facab * faccd * facJop * TBMEarray(1:nTBops, indx)
          ! handle diag tile specially, considering the lower triangle
          if (icd .gt. iab) then
            ! note: this phase is (-1)^iab rather than (-1)^icd because
            ! conjugating the matrix element gives a phase (-1)^(icd-iab)
            facabp = wigab(icd-jabmin+1) * sqrt(2.d0*icd+1.d0)
            faccd = wigcd(iab-jcdmin+1) * sqrt(2.d0*iab+1.d0) * myphase(iab)
            if (deltamj .ge. 0) then
               facJop = Wig3J_Array_JJK(-idj, icd, mab)
            else
               facJop = Wig3J_Array_JJK(-idj, icd, -mab) *     &
                    myphase(iab+Jop+icd)
            endif
            xops(1:nTBops) = xops(1:nTBops) + &
            facabp * faccd * facJop * TBMEarray(1:nTBops, indx)
          endif
       end do
    end do
    !
    if (abIDN) then
       xops(1:nTBops) = xops(1:nTBops) * sqr2
    endif
    if (cdIDN) then
       xops(1:nTBops) = xops(1:nTBops) * sqr2
    endif
    !
    ! note: myphase((j2a-j2b+j2c-j2d)/2)=1 on diag because j2a==j2c and j2b==j2d
    xops(1:nTBops) = myphase(-mcd) * xops(1:nTBops)
    return
    !
  end subroutine TBops_EvalDiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module TBops
