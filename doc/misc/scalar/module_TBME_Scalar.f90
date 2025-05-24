
module TBME
  implicit none
  private
  public read_TBME_ascii, read_TBME_bin
  public Init_H2full, Finalize_H2full
  public ntbme_pp, ntbme_nn, ntbme_pn
  public H2full_pp, H2full_nn, H2full_pn
  public TBMEfull_pp, TBMEfull_nn, TBMEfull_pn
  public J2max_pp, J2max_nn, J2max_pn, WT2mx_pp, WT2mx_nn, WT2mx_pn
  !public ntps_pp, ntps_nn, ntps_pn
  public ntpsJ_pp, ntpsJ_nn, ntpsJ_pn, ntbmeJ_pp, ntbmeJ_nn, ntbmeJ_pn
  public tpsJindx_pp, tpsJindx_nn, tpsJindx_pn
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  integer :: J2max_pp, J2max_nn, J2max_pn
  real(4) :: WT2mx_pp, WT2mx_nn, WT2mx_pn
  integer :: ntbme_pp, ntps_pp, ntbme_nn, ntps_nn, ntbme_pn, ntps_pn
  integer, dimension(:,:), allocatable :: ntbmeJ_pp, ntbmeJ_nn, ntbmeJ_pn
  integer, dimension(:,:), allocatable :: ntpsJ_pp, ntpsJ_nn, ntpsJ_pn
  integer, dimension(:,:), allocatable :: tpsJindx_pp, tpsJindx_nn
  integer, dimension(:,:,:), allocatable :: tpsJindx_pn
  integer(kind=2), dimension(:,:), allocatable :: TPJstates_pp, TPJstates_nn, TPJstates_pn
  real(4), parameter :: zero=0.0
  real(4), dimension(:), allocatable ::  H2full_pp, H2full_nn, H2full_pn
  real(4), dimension(:,:), allocatable :: TBMEfull_pp, TBMEfull_nn, TBMEfull_pn
  !
contains
  !
  subroutine read_TBME_ascii(hamfile) 
    character(LEN=128), intent(in) :: hamfile
    !
    integer :: versionnumber, fh
    character(LEN=128) :: comments
    character(LEN=1) :: key
    !
    fh = 17
    open(unit=fh, file=TRIM(hamfile)//'.dat', status='old', action='read')
    ! allow for comments
    key = '#'
    do while (key .eq. '#')
       read(unit=fh, fmt='(a)', advance='NO') key
       backspace(fh)         
       if (key .eq. '#') then
          read(fh, *) comments
       endif
    enddo
    ! header info
    read(fh, *) versionnumber
    if (versionnumber .eq. 15000) then
       call read_TBME_ascii15000(fh)
    elseif (versionnumber .eq. 15099) then
       call read_TBME_ascii15099(fh)
    else 
       call cancelall(340)
    endif
    !
    ! done reading
    close(unit=fh, status='keep')
    !
    return
  end subroutine read_TBME_ascii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine read_TBME_ascii15000(fh)
    use SPbasis, only: norbt, norb_p, norb_n, n_orb, l_orb, j2_orb, pr_orb, wt_orb
    integer, intent(in) :: fh
    !
    integer :: i, j, k
    integer :: a, na, la, j2a, pra, cls, tmplabels
    integer :: ia,ib,ic, id
    integer :: jt, pt, indx
    real :: wta, x, y, z, phase
    real(4) :: matel
    !
    if (allocated(j2_orb)) then
       read(fh, *) i, j       ! norb has to be the same as input norb
       if (i .ne. norb_p) call cancelall(341)
       if (j .ne. norb_n) call cancelall(342)
       read(fh, *) tmplabels
       if (tmplabels .eq. 7) then
          do j = 1, norbt
             read(fh, *) a, na, la, j2a, pra, cls, wta
          enddo
       elseif (tmplabels .eq. 5) then
          do j = 1, norbt
             read(fh, *) a, j2a, pra, cls, wta
          enddo
       else
          call cancelall(349)
       endif
    else
       read(fh, *) norb_p, norb_n
       norbt = norb_p + norb_n
       read(fh,*) tmplabels
       allocate(j2_orb(norbt), pr_orb(norbt))
       allocate(wt_orb(norbt))
       if (tmplabels .eq. 7) then
          allocate(n_orb(norbt), l_orb(norbt))
          do j = 1, norbt
             read(fh, *) a, na, la, j2a, pra, cls, wta
             if (cls .eq. 1) then
                n_orb(a)  = na 
                l_orb(a)  = la
                j2_orb(a) = j2a
                pr_orb(a) = pra 
                wt_orb(a) = wta
             elseif (cls .eq. 2) then
                ! label a runs from 1 to norb_n for IO of neutron SP info,
                ! but from norb_p+1 to norb_p+norb_n in MFDn
                a = a + norb_p
                n_orb(a)  = na 
                l_orb(a)  = la
                j2_orb(a) = j2a
                pr_orb(a) = pra 
                wt_orb(a) = wta
             else
                call cancelall(347)
             endif
          enddo
       elseif (tmplabels .eq. 5) then
          do j = 1, norbt
             read(fh, *) a, j2a, pra, cls, wta
             if (cls .eq. 1) then
                j2_orb(a) = j2a
                pr_orb(a) = pra 
                wt_orb(a) = wta
             elseif (cls .eq. 2) then
                a = a + norb_p
                j2_orb(a) = j2a
                pr_orb(a) = pra 
                wt_orb(a) = wta
             else
                call cancelall(345)
             endif
          enddo
       else
          call cancelall(349)
       endif
    endif
    !
    read(fh, *) i, j, k
    i = i/2
    j = j/2
    k = k/2
    read(fh, *) x, y, z
    !
    if (.not. allocated(ntbmeJ_pp)) then
       J2max_pp = i
       J2max_nn = j
       J2max_pn = k
       !
       WT2mx_pp = x
       WT2mx_nn = y
       WT2mx_pn = z
       !
       call Init_H2full
       !
    elseif ( (i.gt.J2max_pp).or.(j.gt.J2max_nn).or.(k.gt.J2max_pn) .or. &
         ((x-WT2mx_pp).gt.1d-4) .or. ((y-WT2mx_nn).gt.1d-4) .or. ((z-WT2mx_pn).gt.1d-4) ) then
       !
       print*, 'Inconsistent interaction files'
       print*, i, j, k
       print*, J2max_pp, J2max_nn, J2max_pn
       print*, x, y, z
       print*, WT2mx_pp, WT2mx_nn, WT2mx_pn
       call cancelall(351)
       !
    endif
    !
    read(fh, *) i, j, k
    if (i .ne. ntbme_pp) call cancelall(353)
    if (j .ne. ntbme_nn) call cancelall(354)
    if (k .ne. ntbme_pn) call cancelall(355)
    !
    ! actual matrix elements
    do i = 1, ntbme_pp
       read(fh, *) ia, ib, ic, id, jt, cls, matel
       jt = jt / 2
       pt = pr_orb(ia) * pr_orb(ib)
       pt = (1-pt) / 2        ! convert +1 to 0 and -1 to 1
       call retrieveTBMEindex_IDN(ia, ib, ic, id, jt, pt, J2max_pp, 0, norb_p, &
            j2_orb, ntbmeJ_pp, ntpsJ_pp, tpsJindx_pp, phase, indx)
       H2full_pp(indx) = phase * matel
    enddo
    !
    if (norb_n .eq. 0) then
       ! done reading
       close(unit=fh, status='keep')
       return
    endif
    !
    do i = 1, ntbme_nn
       read(fh, *) ia, ib, ic, id, jt, cls, matel
       jt = jt / 2
       pt = pr_orb(ia+norb_p) * pr_orb(ib+norb_p)
       pt = (1-pt) / 2        ! convert +1 to 0 and -1 to 1
       call retrieveTBMEindex_IDN(ia, ib, ic, id, jt, pt, J2max_nn, norb_p, norb_n, &
            j2_orb, ntbmeJ_nn, ntpsJ_nn, tpsJindx_nn, phase, indx)
       H2full_nn(indx) = phase * matel
    enddo
    !
    do i = 1, ntbme_pn
       read(fh, *) ia, ib, ic, id, jt, cls, matel
       jt = jt / 2
       pt = pr_orb(ia) * pr_orb(ib+norb_p)
       pt = (1-pt) / 2        ! convert +1 to 0 and -1 to 1
       call retrieveTBMEindex_DIS(ia, ib, ic, id, jt, pt, J2max_pn, norb_p, norb_n, &
            ntbmeJ_pn, ntpsJ_pn, tpsJindx_pn, indx)
       H2full_pn(indx) = matel
    enddo
    !
    return
  end subroutine read_TBME_ascii15000
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine read_TBME_ascii15099(fh)
    use SPbasis, only: norbt, norb_p, norb_n, n_orb, l_orb, j2_orb, pr_orb, wt_orb
    integer, intent(in) :: fh
    !  
    integer :: i, j, k
    integer :: a, na, la, j2a, cls
    integer :: ia, ib, ic, id, jtab, jtcd
    integer :: jt, pt, indx
    integer :: TwoJ_op, Parity_op, Tz_op
    real :: WT1mx_p, WT1mx_n
    real :: wta, x, y, z, phase
    real(4) :: matel
    !
    if (allocated(j2_orb)) then
       read(fh, *) i, j       ! norb has to be the same as input norb
       if (i .ne. norb_p) call cancelall(341)
       if (j .ne. norb_n) call cancelall(342)
       do j = 1, norbt
          read(fh, *) a, na, la, j2a, cls, wta
       enddo
    else
       read(fh, *) norb_p, norb_n
       norbt = norb_p + norb_n
       allocate(j2_orb(norbt), pr_orb(norbt))
       allocate(n_orb(norbt), l_orb(norbt))
       allocate(wt_orb(norbt))
       do j = 1, norbt
          read(fh, *) a, na, la, j2a, cls, wta
          if (cls .eq. 1) then
             n_orb(a)  = na 
             l_orb(a)  = la
             j2_orb(a) = j2a
             pr_orb(a) = (-1)**la
             wt_orb(a) = wta
          elseif (cls .eq. 2) then
             ! label a runs from 1 to norb_n for IO of neutron SP info,
             ! but from norb_p+1 to norb_p+norb_n in MFDn
             a = a + norb_p
             n_orb(a)  = na 
             l_orb(a)  = la
             j2_orb(a) = j2a
             pr_orb(a) = (-1)**la
             wt_orb(a) = wta
          else
             call cancelall(347)
          endif
       enddo
    endif
    !
    read(fh, *) TwoJ_op, Parity_op, Tz_op ! currently only scalar pos. parity Tz=0 implemented
    read(fh, *) WT1mx_p, WT1mx_n
    read(fh, *) x, y, z
    read(fh, *) i, j, k
    i = i/2
    j = j/2
    k = k/2
    !
    if (.not. allocated(ntbmeJ_pp)) then
       J2max_pp = i
       J2max_nn = j
       J2max_pn = k
       !
       WT2mx_pp = x
       WT2mx_nn = y
       WT2mx_pn = z
       !
       call Init_H2full
       !
    elseif ( (i.gt.J2max_pp).or.(j.gt.J2max_nn).or.(k.gt.J2max_pn) .or. &
         ((x-WT2mx_pp).gt.1d-4) .or. ((y-WT2mx_nn).gt.1d-4) .or. ((z-WT2mx_pn).gt.1d-4) ) then
       !
       print*, 'Inconsistent interaction files'
       print*, i, j, k
       print*, J2max_pp, J2max_nn, J2max_pn
       print*, x, y, z
       print*, WT2mx_pp, WT2mx_nn, WT2mx_pn
       call cancelall(351)
       !
    endif
    !
    read(fh, *) i, j, k
    if (i .ne. ntbme_pp) call cancelall(353)
    if (j .ne. ntbme_nn) call cancelall(354)
    if (k .ne. ntbme_pn) call cancelall(355)
    !
    ! actual matrix elements
    do i = 1, ntbme_pp
       read(fh, *) ia, ib, ic, id, jtab, jtcd, cls, matel
       jt = jtab / 2
       pt = pr_orb(ia) * pr_orb(ib)
       pt = (1-pt) / 2        ! convert +1 to 0 and -1 to 1
       call retrieveTBMEindex_IDN(ia, ib, ic, id, jt, pt, J2max_pp, 0, norb_p, &
            j2_orb, ntbmeJ_pp, ntpsJ_pp, tpsJindx_pp, phase, indx)
       H2full_pp(indx) = phase * matel
    enddo
    !
    if (norb_n .eq. 0) then
       ! done reading
       close(unit=fh, status='keep')
       return
    endif
    !
    do i = 1, ntbme_nn
       read(fh, *) ia, ib, ic, id, jtab, jtcd, cls, matel
       jt = jtab / 2
       pt = pr_orb(ia+norb_p) * pr_orb(ib+norb_p)
       pt = (1-pt) / 2        ! convert +1 to 0 and -1 to 1
       call retrieveTBMEindex_IDN(ia, ib, ic, id, jt, pt, J2max_nn, norb_p, norb_n, &
            j2_orb, ntbmeJ_nn, ntpsJ_nn, tpsJindx_nn, phase, indx)
       H2full_nn(indx) = phase * matel
    enddo
    !
    do i = 1, ntbme_pn
       read(fh, *) ia, ib, ic, id, jtab, jtcd, cls, matel
       jt = jtab / 2
       pt = pr_orb(ia) * pr_orb(ib+norb_p)
       pt = (1-pt) / 2        ! convert +1 to 0 and -1 to 1
       call retrieveTBMEindex_DIS(ia, ib, ic, id, jt, pt, J2max_pn, norb_p, norb_n, &
            ntbmeJ_pn, ntpsJ_pn, tpsJindx_pn, indx)
       H2full_pn(indx) = matel
    enddo
    !
    return
  end subroutine read_TBME_ascii15099

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_TBME_bin(hamfile) 
    character(LEN=128), intent(in) :: hamfile
    !
    integer :: versionnumber, fh
    !
    fh = 17
    open(unit=fh, file=TRIM(hamfile)//'.bin', status='old', action='read', form='unformatted')
    ! header info
    read(fh) versionnumber
    if (versionnumber .eq. 15000) then
       call read_TBME_bin15000(fh)
    elseif (versionnumber .eq. 15099) then
       call read_TBME_bin15099(fh)
    else 
       call cancelall(370)
    endif
    !
    ! done reading
    close(unit=fh, status='keep')      
    !
    return
  end subroutine read_TBME_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_TBME_bin15000(fh) 
    use SPbasis, only: norbt, norb_p, norb_n, n_orb, l_orb, j2_orb, pr_orb, wt_orb
    integer, intent(in) :: fh
    !     
    integer :: itemp(3), idummy(norbt), i, nlabels
    real(kind=4) :: xtemp(3)
    !
    print*, ' reading 15000 binary (obsolete)'
    !
    if (allocated(j2_orb)) then
       read(fh) itemp(1:3)    ! norb has to be the same as input norb
       if (itemp(1) .ne. norb_p) call cancelall(371)
       if (itemp(2) .ne. norb_n) call cancelall(372)
       do i = 1, itemp(3)
          read(fh) idummy(1:norb_p)
       enddo
       if (norb_n .gt. 0) then
          do i = 1, itemp(3)
             read(fh) idummy(1:norb_n)
          enddo
       endif
    else
       read(fh) norb_p, norb_n, nlabels
       norbt = norb_p + norb_n
       if (nlabels .eq. 5) then
          allocate(n_orb(norbt), l_orb(norbt))
          read(fh) n_orb(1:norb_p)
          read(fh) l_orb(1:norb_p)
       endif
       allocate(j2_orb(norbt), pr_orb(norbt))
       allocate(wt_orb(norbt))
       read(fh) j2_orb(1:norb_p)
       read(fh) pr_orb(1:norb_p)
       read(fh) wt_orb(1:norb_p)
       if (norb_n .gt. 0) then
          if (nlabels .eq. 5) then
             read(fh) n_orb(norb_p+1:norbt)
             read(fh) l_orb(norb_p+1:norbt)
          endif
          read(fh) j2_orb(norb_p+1:norbt)
          read(fh) pr_orb(norb_p+1:norbt)
          read(fh) wt_orb(norb_p+1:norbt)
       endif
    endif
    !
    read(fh) itemp(1:3)
    itemp = itemp / 2
    read(fh) xtemp(1:3)
    !
    if ( (itemp(1).lt.J2max_pp) .or. (itemp(2).lt.J2max_nn) .or. (itemp(3).lt.J2max_pn) &
         .or. ((xtemp(1)-WT2mx_pp).lt.-1d-4) .or. ((xtemp(2)-WT2mx_nn).lt.-1d-4)        &
         .or. ((xtemp(3)-WT2mx_pn).lt.-1d-4) ) then
       print*, 'input file too small'
       print*, itemp(1), itemp(2), itemp(3)
       print*, J2max_pp, J2max_nn, J2max_pn
       print*, xtemp(1), xtemp(2), xtemp(3)
       print*, WT2mx_pp, WT2mx_nn, WT2mx_pn
       call cancelall(380)
    endif
    !
    if (.not. allocated(ntbmeJ_pp)) then
       J2max_pp = itemp(1)
       J2max_nn = itemp(2)
       J2max_pn = itemp(3)
       !
       WT2mx_pp = xtemp(1)
       WT2mx_nn = xtemp(2)
       WT2mx_pn = xtemp(3)
       !
       call Init_H2full
       !
    elseif ( (itemp(1).gt.J2max_pp) .or. (itemp(2).gt.J2max_nn) .or.    &
         (itemp(3).gt.J2max_pn) .or. ((xtemp(1)-WT2mx_pp).gt.1d-4) .or. &
         ((xtemp(2)-WT2mx_nn).gt.1d-4) .or. ((xtemp(3)-WT2mx_pn).gt.1d-4) ) then
       !
       print*, 'Inconsistent interaction files'
       print*, itemp(1), itemp(2), itemp(3)
       print*, J2max_pp, J2max_nn, J2max_pn
       print*, xtemp(1), xtemp(2), xtemp(3)
       print*, WT2mx_pp, WT2mx_nn, WT2mx_pn
       call cancelall(381)
       !
    endif
    !
    read(fh) itemp(1:3)
    if (itemp(1) .ne. ntbme_pp) call cancelall(383)     
    if (itemp(2) .ne. ntbme_nn) call cancelall(384)
    if (itemp(3) .ne. ntbme_pn) call cancelall(385)
    ! actual matrix elements
    read(fh) H2full_pp(1:ntbme_pp)
    if (norb_n .gt. 0) then
       read(fh) H2full_nn(1:ntbme_nn)
       read(fh) H2full_pn(1:ntbme_pn)
    endif
    !
    return
  end subroutine read_TBME_bin15000

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_TBME_bin15099(fh) 
    use SPbasis, only: norbt, norb_p, norb_n, n_orb, l_orb, j2_orb, pr_orb, wt_orb
    integer, intent(in) :: fh
    !  
    integer :: norb_1, norb_2, i, j, k, TwoJ_op, Parity_op, Tz_op
    integer, dimension(norbt) :: dummy_n, dummy_l, dummy_j
    real(4) :: x, y, z, WT1mx_p, WT1mx_n
    real(4), dimension(norbt) :: wtdummy
    !  
    print*, ' reading 15099 binary'
    !
    read(fh) norb_1, norb_2
    if (allocated(j2_orb)) then
       if (norb_1 .ne. norb_p) call cancelall(371)
       if (norb_2 .ne. norb_n) call cancelall(372)
       read(fh) dummy_n(1:norb_p)
       read(fh) dummy_l(1:norb_p)
       read(fh) dummy_j(1:norb_p)
       read(fh) wtdummy(1:norb_p)
       if (norb_n .gt. 0) then
          read(fh) dummy_n(norb_p+1:norbt)
          read(fh) dummy_l(norb_p+1:norbt)
          read(fh) dummy_j(norb_p+1:norbt)
          read(fh) wtdummy(norb_p+1:norbt)
       endif
       do i = 1, norbt
          if (dummy_n(i) .ne. n_orb(i)) call cancelall(374)
       enddo
       do i = 1, norbt
          if (dummy_l(i) .ne. l_orb(i)) call cancelall(374)
       enddo
       do i = 1, norbt
          if (dummy_j(i) .ne. j2_orb(i)) call cancelall(374)
       enddo
       do i = 1, norbt
          if (wtdummy(i) .ne. wt_orb(i)) call cancelall(374)
       enddo
    else
       norbt = norb_p + norb_n
       allocate(n_orb(norbt), l_orb(norbt))
       allocate(j2_orb(norbt), pr_orb(norbt))
       allocate(wt_orb(norbt))
       read(fh) n_orb(1:norb_p)
       read(fh) l_orb(1:norb_p)
       read(fh) j2_orb(1:norb_p)
       read(fh) wt_orb(1:norb_p)
       if (norb_n .gt. 0) then
          read(fh) n_orb(norb_p+1:norbt)
          read(fh) l_orb(norb_p+1:norbt)
          read(fh) j2_orb(norb_p+1:norbt)
          read(fh) wt_orb(norb_p+1:norbt)
       endif
       do i = 1, norbt
          pr_orb(i) = (-1)**l_orb(i)
       enddo
    endif
    !
    read(fh) TwoJ_op, Parity_op, Tz_op ! currently on scalar pos. parity Tz=0 implemented
    read(fh) WT1mx_p, WT1mx_n
    read(fh) x, y, z
    read(fh) i, j, k
    i = i/2
    j = j/2
    k = k/2
    !
    if ( (i.lt.J2max_pp).or.(j.lt.J2max_nn).or.(k.lt.J2max_pn)    &
         .or. ((WT2mx_pp-x).gt.1d-4) .or. ((WT2mx_nn-y).gt.1d-4)  &
         .or. ((WT2mx_pn-z).gt.1d-4) ) then
       print*, 'input file too small'
       print*, 'read', x, y, z
       print*, 'expected', WT2mx_pp, WT2mx_nn, WT2mx_pn
       print*, 'read', i, j, k
       print*, 'expected', J2max_pp, J2max_nn, J2max_pn
       call cancelall(380)
    endif
    !
    if (.not. allocated(ntbmeJ_pp)) then
       J2max_pp = i
       J2max_nn = j
       J2max_pn = k
       !
       WT2mx_pp = x
       WT2mx_nn = y
       WT2mx_pn = z
       !
       call Init_H2full
       !
    elseif ( (i.gt.J2max_pp).or.(j.gt.J2max_nn).or.(k.gt.J2max_pn) &
         .or. ((x-WT2mx_pp).gt.1d-4) .or. ((y-WT2mx_nn).gt.1d-4)   &
         .or. ((z-WT2mx_pn).gt.1d-4) ) then
       !
       print*, 'Inconsistent interaction files'
       print*, i, j, k
       print*, J2max_pp, J2max_nn, J2max_pn
       print*, x, y, z
       print*, WT2mx_pp, WT2mx_nn, WT2mx_pn
       call cancelall(381)
       !
    endif
    !
    read(fh) i, j, k
    if (i .ne. ntbme_pp) call cancelall(383)     
    if (j .ne. ntbme_nn) call cancelall(384)
    if (k .ne. ntbme_pn) call cancelall(385)
    ! actual matrix elements
    read(fh) H2full_pp(1:ntbme_pp)
    if (norb_n .gt. 0) then
       read(fh) H2full_nn(1:ntbme_nn)
       read(fh) H2full_pn(1:ntbme_pn)
    endif
    !
    return
  end subroutine read_TBME_bin15099

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Init_H2full
    use SPbasis, only: norb_p, norb_n, j2_orb, pr_orb, wt_orb
    !
    allocate(ntbmeJ_pp(0:J2max_pp+1, 0:1))
    allocate(ntpsJ_pp(0:J2max_pp+1, 0:1))
    allocate(tpsJindx_pp(0:J2max_pp, norb_p*(norb_p+1)/2) )
    !
    call setTBMEindex_IDN(J2max_pp, WT2mx_pp, 0, norb_p, wt_orb, j2_orb, pr_orb, &
         ntbmeJ_pp, ntpsJ_pp, tpsJindx_pp, ntps_pp, ntbme_pp)
    allocate(TPJstates_pp(2, ntps_pp))
    call setTPJstates_IDN(J2max_pp, WT2mx_pp, 0, norb_p, wt_orb, j2_orb, pr_orb, &
         ntps_pp, TPJstates_pp)
    !
    allocate(ntbmeJ_nn(0:J2max_nn+1, 0:1))
    allocate(ntpsJ_nn(0:J2max_nn+1, 0:1))
    allocate(tpsJindx_nn(0:J2max_nn, norb_n*(norb_n+1)/2) )
    !
    call setTBMEindex_IDN(J2max_nn, WT2mx_nn, norb_p, norb_n, wt_orb, j2_orb, pr_orb, &
         ntbmeJ_nn, ntpsJ_nn, tpsJindx_nn, ntps_nn, ntbme_nn)
    !
    allocate(TPJstates_nn(2, ntps_nn))
    call setTPJstates_IDN(J2max_nn, WT2mx_nn, norb_p, norb_n, wt_orb, j2_orb, pr_orb, &
         ntps_nn, TPJstates_nn)
    !
    allocate(ntbmeJ_pn(0:J2max_pn+1, 0:1))
    allocate(ntpsJ_pn(0:J2max_pn+1, 0:1))
    allocate(tpsJindx_pn(0:J2max_pn, norb_p, norb_n))
    !
    call setTBMEindex_DIS(J2max_pn, WT2mx_pn, norb_p, norb_n, wt_orb, j2_orb, pr_orb, &
         ntbmeJ_pn, ntpsJ_pn, tpsJindx_pn, ntps_pn, ntbme_pn)
    !
    allocate(TPJstates_pn(2, ntps_pn))
    call setTPJstates_DIS(J2max_pn, WT2mx_pn, norb_p, norb_n, wt_orb, j2_orb, pr_orb, &
         ntps_pn, TPJstates_pn)
    !
    allocate(H2full_pp(ntbme_pp))
    H2full_pp(1:ntbme_pp) = zero
    allocate(H2full_nn(ntbme_nn))
    H2full_nn(1:ntbme_nn) = zero
    allocate(H2full_pn(ntbme_pn))
    H2full_pn(ntbme_pn) = zero
    !  
    return
  end subroutine Init_H2full

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine Finalize_H2full
    ! de-allocate H2full arrays
    deallocate(ntbmeJ_pp)
    deallocate(ntpsJ_pp)
    deallocate(tpsJindx_pp)
    deallocate(TPJstates_pp)
    if (allocated(H2full_pp)) deallocate(H2full_pp)
    if (allocated(TBMEfull_pp)) deallocate(TBMEfull_pp)
    !
    deallocate(ntbmeJ_nn)
    deallocate(ntpsJ_nn)
    deallocate(tpsJindx_nn)
    deallocate(TPJstates_nn)
    if (allocated(H2full_nn)) deallocate(H2full_nn)
    if (allocated(TBMEfull_nn)) deallocate(TBMEfull_nn)
    !
    deallocate(ntbmeJ_pn)
    deallocate(ntpsJ_pn)
    deallocate(tpsJindx_pn)
    deallocate(TPJstates_pn)
    if (allocated(H2full_pn)) deallocate(H2full_pn)
    if (allocated(TBMEfull_pn)) deallocate(TBMEfull_pn)
    !
    return
  end subroutine Finalize_H2full

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine setTBMEindex_IDN(jmx, wtmx, offset, norbt, wt_orb, j2_orb, pr_orb, &
       ntbmeJ, ntpsJ, tpsJindex, ntps, ngs)
    !
    integer, intent(in) :: jmx, offset, norbt
    real(4), intent(in) :: wtmx
    real(4), dimension(offset+norbt), intent(in) :: wt_orb
    integer, dimension(offset+norbt), intent(in) :: j2_orb, pr_orb
    integer, dimension(0:jmx+1,0:1), intent(out) :: ntbmeJ, ntpsJ
    integer, dimension(0:jmx, norbt*(norbt+1)/2), intent(out) :: tpsJindex
    integer, intent(out) :: ntps, ngs
    ! local variables
    integer :: j, p, njp, ia, j2a, ipa, ib, j2b, ipb, mnj, mxj, ka, kb, par
    real(4) :: wta, wtb
    !
    tpsJindex(0:jmx, 1:norbt*(norbt+1)/2) = 0
    ntps = 0
    ngs = 0
    do j = 0, jmx
       do p = 0, 1
          ! par = (-1)**p
          njp = 0             ! Initialize njp to zero
          do ia = offset+1, offset+norbt
             wta = wt_orb(ia)
             j2a = j2_orb(ia)
             ipa = pr_orb(ia)
             ka = ia - offset
             do ib = ia, offset+norbt
                if ((ib.eq.ia).and.mod(j,2).eq.1) cycle
                wtb = wt_orb(ib)
                if (wta+wtb .gt. wtmx) cycle
                j2b = j2_orb(ib)
                ipb = pr_orb(ib)
                ! if (par*ipa*ipb .eq. -1) cycle ! parity constraint
                if (abs(ipa-ipb)/2 .ne. p) cycle ! parity constraint
                !
                mnj = abs(j2a - j2b) / 2
                mxj = ( j2a + j2b ) / 2
                !
                if (j.lt.mnj.or.j.gt.mxj) cycle
                !
                njp = njp + 1
                kb = ib - offset
                ! Use the 'wrong order' for indices, but better for retrieval:
                tpsJindex(j,(ka+kb*(kb-1)/2)) = njp
             enddo
          enddo
          ! Offset for TBMEs 
          ntbmeJ(j,p) = ngs
          ! Number of TwoParticleStates coupled to J
          ntpsJ(j,p) = njp
          ! Cumulative number of TwoParticleStates and TBMEs
          ntps = ntps + njp
          ngs = ngs + njp*(njp+1)/2
       enddo
    enddo
    ntpsJ(jmx+1,0:1) = 0
    ntbmeJ(jmx+1,0:1) = ngs
    ! print*, ' Total number of two-body states :', ntps
    ! print*, ' Total number of two-body mat.els:', ngs
    !
    return
  end subroutine setTBMEindex_IDN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine setTPJstates_IDN(jmx, wtmx, offset, norbt, &
       wt_orb, j2_orb, pr_orb, ntps, TPJstates)
    !
    integer, intent(in) :: jmx, offset, norbt, ntps
    real(4), intent(in) :: wtmx
    real(4), dimension(offset+norbt), intent(in) :: wt_orb
    integer, dimension(offset+norbt), intent(in) :: j2_orb, pr_orb
    integer(kind=2), dimension(2,ntps), intent(out) :: TPJstates
    ! local variables
    integer :: j, p, nj, ia, j2a, ipa, ib, j2b, ipb, mnj, mxj, par
    real(4) :: wta, wtb
    !
    TPJstates(1:2, 1:ntps) = 0
    nj = 0
    do j = 0, jmx
       do p = 0, 1
          ! par = (-1)**p
          do ia = offset+1, offset+norbt
             wta = wt_orb(ia)
             j2a = j2_orb(ia)
             ipa = pr_orb(ia)
             do ib = ia, offset+norbt
                if ((ia.eq.ib).and.mod(j,2).eq.1) cycle
                wtb = wt_orb(ib)
                if (wta+wtb .gt. wtmx) cycle
                j2b = j2_orb(ib)
                ipb = pr_orb(ib)
                ! if (par*ipa*ipb .eq. -1) cycle ! parity constraint
                if (abs(ipa-ipb)/2 .ne. p) cycle ! parity constraint
                !
                mnj = abs(j2a - j2b) / 2
                mxj = ( j2a + j2b ) / 2
                !
                if (j.lt.mnj.or.j.gt.mxj) cycle
                !
                nj = nj + 1                  
                TPJstates(1, nj) = ia - offset
                TPJstates(2, nj) = ib - offset
                ! TPJstates(3, nj) = j
             enddo
          enddo
       enddo
    enddo
    !
    return
  end subroutine setTPJstates_IDN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine setTBMEindex_DIS(jmx, wtmx, norb1, norb2, &
       wt_orb, j2_orb, pr_orb, ntbmeJ, ntpsJ, tpsJindex, ntps, ngs)
    !
    integer, intent(in) :: jmx, norb1, norb2
    real(4), intent(in) :: wtmx
    real(4), dimension(norb1+norb2), intent(in) :: wt_orb
    integer, dimension(norb1+norb2), intent(in) :: j2_orb, pr_orb
    integer, dimension(0:jmx+1,0:1), intent(out) :: ntbmeJ, ntpsJ
    integer, dimension(0:jmx, norb1, norb2), intent(out) :: tpsJindex
    integer, intent(out) :: ntps, ngs
    ! local variables
    integer :: j, p, njp, ia, j2a, ipa, ib, j2b, ipb, mnj, mxj, kb, par
    real(4) :: wta, wtb
    !
    tpsJindex(0:jmx, 1:norb1, 1:norb2) = 0
    ntps = 0
    ngs = 0
    do j = 0, jmx
       do p = 0, 1
          ! par = (-1)**p
          njp = 0             ! Initialize njp to zero
          do ia = 1, norb1
             wta = wt_orb(ia)
             j2a = j2_orb(ia)
             ipa = pr_orb(ia)
             do kb = 1, norb2
                ib = norb1 + kb
                wtb = wt_orb(ib)
                if (wta+wtb .gt. wtmx) cycle
                j2b = j2_orb(ib)
                ipb = pr_orb(ib)
                ! if (par*ipa*ipb .eq. -1) cycle ! parity constraint
                if (abs(ipa-ipb)/2 .ne. p) cycle ! parity constraint
                !
                mnj = abs(j2a - j2b) / 2
                mxj = ( j2a + j2b ) / 2
                !
                if (j.lt.mnj.or.j.gt.mxj) cycle
                njp = njp + 1
                ! Use the 'wrong order' for indices, but better for retrieval:
                tpsJindex(j,ia,kb) = njp
             enddo
          enddo
          ! Offset for TBMEs 
          ntbmeJ(j,p) = ngs
          ! Number of TwoParticleStates coupled to J
          ntpsJ(j,p) = njp
          ! Cumulative number of ParticleStates and TBMEs
          ntps = ntps + njp
          ngs = ngs + njp*(njp+1)/2
       enddo
    enddo
    ntpsJ(jmx+1,0:1) = 0
    ntbmeJ(jmx+1,0:1) = ngs
    ! print*, ' Total number of two-body states :', ntps
    ! print*, ' Total number of two-body mat.els:', ngs
    !
    return
  end subroutine setTBMEindex_DIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine setTPJstates_DIS(jmx, wtmx, norb1, norb2, &
       wt_orb, j2_orb, pr_orb, ntps, TPJstates)
    integer, intent(in) :: jmx, norb1, norb2, ntps
    real(4), intent(in) :: wtmx
    real(4), dimension(norb1+norb2), intent(in) :: wt_orb
    integer, dimension(norb1+norb2), intent(in) :: j2_orb, pr_orb
    integer(kind=2), dimension(2,ntps), intent(out) :: TPJstates
    ! local variables
    integer :: j, p, nj, ia, j2a, ipa, ib, j2b, ipb, mnj, mxj, par
    real(4) :: wta, wtb
    !
    TPJstates(1:2, 1:ntps) = 0
    nj = 0
    do j = 0, jmx
       do p = 0, 1
          ! par = (-1)**p
          do ia = 1, norb1
             wta = wt_orb(ia)
             j2a = j2_orb(ia)
             ipa = pr_orb(ia)
             do ib = norb1+1, norb1+norb2
                wtb = wt_orb(ib)
                if (wta+wtb .gt. wtmx) cycle
                j2b = j2_orb(ib)
                ipb = pr_orb(ib)
                ! if (par*ipa*ipb .eq. -1) cycle ! parity constraint
                if (abs(ipa-ipb)/2 .ne. p) cycle ! parity constraint
                !
                mnj = abs(j2a - j2b) / 2
                mxj = ( j2a + j2b ) / 2
                !
                if (j.lt.mnj.or.j.gt.mxj) cycle
                !
                nj = nj + 1                  
                TPJstates(1, nj) = ia
                TPJstates(2, nj) = ib - norb1
                ! TPJstates(3, nj) = j
             enddo
          enddo
       enddo
    enddo
    !
    return
  end subroutine setTPJstates_DIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine retrieveTBMEindex_IDN(iaa, ibb, icc, idd, jt, pt, jmx, &
       offset, nobt, j2_orb, ntbmeJ, ntpsJ, tpsJindx, phase, indx)
    !
    integer, intent(in) :: jt, pt, jmx, offset, nobt
    integer, intent(in) :: iaa, ibb, icc, idd
    integer, dimension(nobt+offset), intent(in) :: j2_orb
    integer, dimension(0:jmx+1,0:1), intent(in) :: ntbmeJ, ntpsJ
    integer, dimension(0:jmx,nobt*(nobt+1)/2), intent(in) :: tpsJindx
    real, intent(out) :: phase
    integer, intent(out) :: indx
    ! local variables
    integer :: ia, ib, ic, id, itmp
    integer :: iph, njab, njcd
    real, external :: myphase
    !
    ! order ia < ib
    if (iaa .gt. ibb) then
       ia = ibb
       ib = iaa
       iph = jt + (j2_orb(ia+offset)+j2_orb(ib+offset))/2 + 1
    else
       ia = iaa
       ib = ibb
       iph = 0
    endif
    ! order ic < id
    if (icc .gt. idd) then
       ic = idd
       id = icc
       iph = iph + jt + (j2_orb(ic+offset)+j2_orb(id+offset))/2 + 1
    else
       ic = icc
       id = idd
    endif
    phase = myphase(iph)
    ! order (ia,ib) < (ic,id)
    if ((ia.gt.ic).or.((ia.eq.ic).and.(ib.gt.id))) then
       itmp = ia
       ia = ic
       ic = itmp
       itmp = ib
       ib = id
       id = itmp
    endif
    !
    ! actual index calculation
    njab = tpsJindx(jt, (ia+ib*(ib-1)/2))
    njcd = tpsJindx(jt, (ic+id*(id-1)/2))
    !
    indx = ntbmeJ(jt,pt) &     ! offset for number of TBMEs
         + ntpsJ(jt,pt) * (njab-1) + njcd - njab * (njab-1) / 2
    !
    return
  end subroutine retrieveTBMEindex_IDN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine retrieveTBMEindex_DIS(iaa, ibb, icc, idd, jt, pt, jmx, &
       nobt1, nobt2, ntbmeJ, ntpsJ, tpsJindx, indx)
    !
    integer, intent(in) :: jt, pt, jmx, nobt1, nobt2
    integer, intent(in) :: iaa, ibb, icc, idd
    integer, dimension(0:jmx+1,0:1), intent(in) :: ntbmeJ, ntpsJ
    integer, dimension(0:jmx,nobt1,nobt2), intent(in) :: tpsJindx
    integer, intent(out) :: indx
    ! local variables
    integer :: ia, ib, ic, id, njab, njcd
    !
    ! order (ia,ib) < (ic,id)
    if ((iaa.gt.icc).or.((iaa.eq.icc).and.(ibb.gt.idd))) then
       ia = icc
       ic = iaa
       ib = idd
       id = ibb
    else
       ia = iaa
       ib = ibb
       ic = icc
       id = idd
    endif
    !
    ! actual index calculation
    njab = tpsJindx(jt, ia, ib)
    njcd = tpsJindx(jt, ic, id)
    !
    indx = ntbmeJ(jt,pt) &     ! offset for number of TBMEs
         + ntpsJ(jt,pt) * (njab-1) + njcd - njab * (njab-1) / 2
    !
    return
  end subroutine retrieveTBMEindex_DIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine tbme_IDN(orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d, &
       m2a,m2b,m2c,m2d, pt, jmin, jmax, Jtot2max, norb, &
       ntbme, tpsJindx, ntpsJ, ntbmeJ, TBMEfull, tbme)
    use Wigner3J, only: Retrieve_Wigner3J_Array
    !
    integer, intent(in) :: orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d, &
         m2a,m2b,m2c,m2d, pt, jmin, jmax, Jtot2max, norb, ntbme
    integer, dimension(0:Jtot2max, norb*(norb+1)/2) :: tpsJindx
    integer, dimension(0:Jtot2max+1, 0:1) :: ntpsJ, ntbmeJ
    real(kind=4), dimension(ntbme), intent(in) :: TBMEfull
    real(kind=4), intent(out) :: tbme
    ! local variables
    integer :: jt, njab, njcd, indx
    real(kind=8) :: factor
    real(kind=8), dimension(Jtot2max+1) :: wigab, wigcd
    real(kind=8), parameter :: sqr2 = sqrt(2.d0)
    !
    call Retrieve_Wigner3j_array(j2a,j2b,m2a,m2b,jmin,jmax, wigab)
    call Retrieve_Wigner3j_array(j2c,j2d,m2c,m2d,jmin,jmax, wigcd)
    !
    tbme = 0.0
    do jt = jmin, jmax
       njab = tpsJindx(jt, (orba+orbb*(orbb-1)/2))
       njcd = tpsJindx(jt, (orbc+orbd*(orbd-1)/2))
       if (njab*njcd .eq. 0) cycle
       !
       indx = ntbmeJ(jt,pt) &
            + ntpsJ(jt,pt) * (njab-1) + njcd - njab * (njab-1) / 2
       !
       factor = wigab(jt-jmin+1) * wigcd(jt-jmin+1) * (2.0*jt+1.0)
       tbme = tbme + factor * TBMEfull(indx)
    enddo
    ! Needed for identical particles
    if (orba.eq.orbb) then
       tbme = tbme * sqr2
    endif
    if (orbc.eq.orbd) then
       tbme = tbme * sqr2
    endif
    !
    return
  end subroutine tbme_IDN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine tbme_DIS(orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d, &
       m2a,m2b,m2c,m2d, pt, jmin, jmax, Jtot2max, norb1, norb2, &
       ntbme, tpsJindx, ntpsJ, ntbmeJ, TBMEfull, tbme)
    use Wigner3J, only: Retrieve_Wigner3J_Array
    !
    integer, intent(in) :: orba,orbb,orbc,orbd, j2a,j2b,j2c,j2d, &
         m2a,m2b,m2c,m2d, pt, jmin, jmax, Jtot2max, norb1, norb2, ntbme
    integer, dimension(0:Jtot2max, norb1, norb2) :: tpsJindx
    integer, dimension(0:Jtot2max+1, 0:1) :: ntpsJ, ntbmeJ
    real(kind=4), dimension(ntbme), intent(in) :: TBMEfull
    real(kind=4), intent(out) :: tbme
    ! local variables
    integer :: jt, njab, njcd, indx
    real(kind=8) :: factor
    real(kind=8), dimension(Jtot2max+1) :: wigab, wigcd
    !
    call Retrieve_Wigner3j_array(j2a,j2b,m2a,m2b,jmin,jmax, wigab)
    call Retrieve_Wigner3j_array(j2c,j2d,m2c,m2d,jmin,jmax, wigcd)
    !
    tbme = 0.0
    do jt = jmin, jmax
       njab = tpsJindx(jt, orba, orbb)
       njcd = tpsJindx(jt, orbc, orbd)
       if (njab*njcd .eq. 0) cycle
       !
       indx = ntbmeJ(jt,pt) &
            + ntpsJ(jt,pt) * (njab-1) + njcd - njab * (njab-1) / 2
       !
       factor = wigab(jt-jmin+1) * wigcd(jt-jmin+1) * (2.0*jt+1.0)
       tbme = tbme + factor * TBMEfull(indx)
    enddo
    !
    return
  end subroutine tbme_DIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module TBME
