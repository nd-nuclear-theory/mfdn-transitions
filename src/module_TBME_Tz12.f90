
module TBME_Tz12
  use SPbasis, only: norb_p, norb_n, norbt, n_orb, l_orb, j2_orb, pr_orb, wt_orb
  implicit none
  private
  public readTBMEascii, readTBMEbin
  public RetrieveTBMEoffset, RetrieveTBMEindx, Set_TBME_array
  public Jop, Parop, Tzop, nTBMEs_max
  public norb2_pp, nTPpos_pp, nTPneg_pp, nTPorb_pp, negoffset_pp, TPorb_pp
  public norb2_nn, nTPpos_nn, nTPneg_nn, nTPorb_nn, negoffset_nn, TPorb_nn
  public norb2_pn, nTPpos_pn, nTPneg_pn, nTPorb_pn, negoffset_pn, TPorb_pn
  public nTBMEindx_pp, nTBMEs_pp, numTBME_pp, TBME_offset_pp, TBME_pp
  public nTBMEindx_nn, nTBMEs_nn, numTBME_nn, TBME_offset_nn, TBME_nn
  public nTBMEindx_pn, nTBMEs_pn, numTBME_pn, TBME_offset_pn, TBME_pn
  public TBMEarray_pp, TBMEarray_pn, TBMEarray_nn
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  logical :: IDN_pp=.true., IDN_pn=.false., IDN_nn=.true.
  integer :: Jop=-1, Parop=1, gop=0, Tzop=0, nTBMEs_max=0
  integer :: norb2_pp, nTPpos_pp, nTPneg_pp, nTPorb_pp, negoffset_pp
  integer :: norb2_nn, nTPpos_nn, nTPneg_nn, nTPorb_nn, negoffset_nn
  integer :: norb2_pn, nTPpos_pn, nTPneg_pn, nTPorb_pn, negoffset_pn
  real(kind=4) :: WT2max_pp=-1, WT2max_pn=-1, WT2max_nn=-1
  integer :: nTBMEindx_pp, nTBMEs_pp
  integer :: nTBMEindx_nn, nTBMEs_nn
  integer :: nTBMEindx_pn, nTBMEs_pn
  integer, allocatable, dimension(:) :: TPorb_pp, numTBME_pp, TBME_offset_pp
  integer, allocatable, dimension(:) :: TPorb_nn, numTBME_nn, TBME_offset_nn
  integer, allocatable, dimension(:) :: TPorb_pn, numTBME_pn, TBME_offset_pn
  real(kind=4), allocatable, dimension(:) :: TBME_pp, TBME_nn, TBME_pn
  real(kind=4), allocatable, dimension(:,:) :: TBMEarray_pp, TBMEarray_pn, TBMEarray_nn
  !
contains
  !
  subroutine readTBMEascii(TBMEinputfile)
    character(LEN=*), intent(in) :: TBMEinputfile
    ! local variables
    integer :: fh=42, versionnumber
    !
    open(unit=fh, file=TRIM(TBMEinputfile)//'.dat', status='old', action='read')
    read(fh, *) versionnumber
    !
    if (versionnumber .eq. 15200) then
       call readTBMEascii_15200(fh)
    else
       print*, versionnumber
       call cancelall(100)
    endif
    !
    close(unit=fh, status='keep')
    !
    return
  end subroutine readTBMEascii

   subroutine readTBMEbin(TBMEinputfile)
      character(LEN=*), intent(in) :: TBMEinputfile
      ! local variables
      integer :: fh=42, versionnumber
      !
      ! print*, 'Binary input for Tz changing operators not yet implemented'
      ! call cancelall(112)
      !
      open(unit=fh, file=TRIM(TBMEinputfile)//'.bin', &
        status='old', action='read', form='unformatted')
      read(fh) versionnumber
      !
      if (versionnumber .eq. 15200) then
        call readTBMEbin_15200(fh)
      else
        print*, versionnumber
        call cancelall(112)
      endif
      !
      close(unit=fh, status='keep')
      !
      return
   end subroutine readTBMEbin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine readTBMEascii_15200(fh)
    integer, intent(in) :: fh
    ! local variables
    real :: wta
    real, dimension(-1:1) :: WT2max
    integer, dimension(-1:1) :: nTBMEs
    integer :: versionnumber
    integer :: i, a, na, la, tz, indx, offset_indx
    integer :: orba, orbb, orbc, orbd, j2a, j2b, j2c, j2d, noa, nob, noc, nod
    integer :: orbab, orbcd, jab, jcd, parab, parcd, TPorbab, TPorbcd
    integer :: Jop_f, gop_f, Tzop_f
    real(kind=4) :: matel
    logical, external :: pairwiseless
    real, external :: myphase
    !
    if (.not.allocated(j2_orb)) then
       read(fh, *) norb_p, norb_n
       norbt = norb_p + norb_n
       !
       allocate(j2_orb(norbt), pr_orb(norbt))
       allocate(n_orb(norbt), l_orb(norbt))
       allocate(wt_orb(norbt))
       do i = 1, norbt
          read(fh, *) a, na, la, j2a, tz, wta
          n_orb(a)  = na
          l_orb(a)  = la
          j2_orb(a) = j2a
          pr_orb(a) = (-1)**la
          wt_orb(a) = wta
       enddo
       !
    else
       read(fh, *) orba, orbb
       if (orba .ne. norb_p) then
          print*, 'incorrect norb_p', orba, norb_p
          call cancelall(001)
       endif
       if (orbb .ne. norb_n) then
          print*, 'incorrect norb_n', orbb, norb_n
          call cancelall(002)
       endif
       do i = 1, orba+orbb
          read(fh, *) a, na, la, j2a, tz, wta
       enddo
       !
    endif
    !
    read(fh, *) Jop_f, gop_f, Tzop_f
    if (Jop == -1) then
       Jop = Jop_f
       gop = gop_f
       Tzop = Tzop_f
       Parop = (-1)**gop
    else
       if ((Jop_f/=Jop).or.(gop_f/=gop).or.(Tzop_f/=Tzop)) then
          print*, 'incompatible operator quantum numbers'
          print*, Jop, gop, Tzop
          print*, Jop_f, gop_f, Tzop_f
          call cancelall(005)
       endif
    endif
    !
    read(fh, *) WT2max(1), WT2max(0), WT2max(-1)  ! two-body truncation weights
    read(fh, *) nTBMEs(1), nTBMEs(0), nTBMEs(-1)  ! used as consistency check
    !
    if (.not.allocated(TPorb_pp)) then
       ! Set up TwoBody indexing arrays
       WT2max_pp = WT2max(1)
       WT2max_pn = WT2max(0)
       WT2max_nn = WT2max(-1)
       nTBMEs_pp = nTBMEs(1)
       nTBMEs_pn = nTBMEs(0)
       nTBMEs_nn = nTBMEs(-1)
       call setup_TBME_indexing
    else
      if ((abs(WT2max(1)-WT2max_pp) .gt. 1e-6) &
            .or. (abs(WT2max(0)-WT2max_pn) .gt. 1e-6) &
            .or. (abs(WT2max(-1)-WT2max_nn) .gt. 1e-6)) &
      then
         print*, 'tbme file wrong size'
         print*, WT2max(1), WT2max(0), WT2max(-1)
         print*, WT2max_pp, WT2max_pn, WT2max_nn
         call cancelall(006)
      endif
      if ((nTBMEs(1) /= nTBMEs_pp) .or. (nTBMEs(0) /= nTBMEs_pn) &
            .or. (nTBMEs(-1) /= nTBMEs_nn)) &
      then
         print*, 'tbme file wrong size'
         print*, nTBMEs(1), nTBMEs(0), nTBMEs(-1)
         print*, nTBMEs_pp, nTBMEs_pn, nTBMEs_nn
         call cancelall(006)
      endif
    endif
    !
    ! For all operators (Tz=0, Tz=1, Tz=2)
    !  < pp || Tz=0 || pp > matrix elements -- see module_TBME_Tz0
    !  < pp || Tz=1 || pn > matrix elements
    !  < pp || Tz=2 || nn > matrix elements
    do i = 1, nTBMEs_pp
       read(fh, *) orba, orbb, orbc, orbd, jab, jcd, matel
       !
       parab = pr_orb(orba)*pr_orb(orbb)
       parcd = pr_orb(orbc)*pr_orb(orbd)
       j2a = j2_orb(orba)
       j2b = j2_orb(orbb)
       j2c = j2_orb(orbc)
       j2d = j2_orb(orbd)
       jab = jab/2
       jcd = jcd/2
       !
       orbab = norb2_pp - (norb_p-orba+1)*(norb_p-orba+2)/2 + (orbb-orba) + 1
       TPorbab = TPorb_pp(orbab)
       !
       if (Tzop .eq. 1) then
          nod = orbd - norb_p
          orbcd = (orbc-1)*norb_n + nod
          TPorbcd = TPorb_pn(orbcd)
          call retrieveTBMEoffset(parab, parcd,              &
               TPorbab, nTPpos_pp, nTPneg_pp,                &
               TPorbcd, nTPpos_pn, nTPneg_pn,                &
               negoffset_pp, offset_indx)
          call retrieveTBMEindx(orba, orbb, orbc, orbd,      &
               j2a, j2b, j2c, j2d, jab, jcd, Jop, indx)
          !
       elseif (Tzop .eq. 2) then
          noc = orbc - norb_p
          nod = orbd - norb_p
          orbcd = norb2_nn - (norb_n-noc+1)*(norb_n-noc+2)/2 + (nod-noc) + 1
          TPorbcd = TPorb_nn(orbcd)
          call retrieveTBMEoffset(parab, parcd,              &
               TPorbab, nTPpos_pp, nTPneg_pp,                &
               TPorbcd, nTPpos_nn, nTPneg_nn,                &
               negoffset_pp, offset_indx)
          call retrieveTBMEindx(orba, orbb, orbc, orbd,      &
               j2a, j2b, j2c, j2d, jab, jcd, Jop, indx)
          !
       else
          print*, 'invalid Tz', Tzop
          call cancelall(100)
          !
       end if
       !
       indx = indx + TBME_offset_pp(offset_indx)
       TBME_pp(indx) = matel * sqrt(2.d0*jab + 1.d0)
       !
    end do
    !
    ! no pn matrix elements for Tz0 = 2
    if (Tzop .eq. 2) return
    !
    ! Only for Tz=0 and Tz=1 operators
    !  < pn || Tz=0 || pn > matrix elements -- see module_TBME_Tz0
    !  < pn || Tz=1 || nn > matrix elements
    ! nTBMEs_pn should be 0 for Tz=2
    do i = 1, nTBMEs_pn
       read(fh, *) orba, orbb, orbc, orbd, jab, jcd, matel
       !
       parab = pr_orb(orba)*pr_orb(orbb)
       parcd = pr_orb(orbc)*pr_orb(orbd)
       j2a = j2_orb(orba)
       j2b = j2_orb(orbb)
       j2c = j2_orb(orbc)
       j2d = j2_orb(orbd)
       jab = jab/2
       jcd = jcd/2
       !
       nob = orbb - norb_p
       orbab = (orba-1)*norb_n + nob
       TPorbab = TPorb_pn(orbab)
       !
       if (Tzop .eq. 1) then
          noc = orbc - norb_p
          nod = orbd - norb_p
          orbcd = norb2_nn - (norb_n-noc+1)*(norb_n-noc+2)/2 + (nod-noc) + 1
          TPorbcd = TPorb_nn(orbcd)
          call retrieveTBMEoffset(parab, parcd,              &
               TPorbab, nTPpos_pn, nTPneg_pn,                &
               TPorbcd, nTPpos_nn, nTPneg_nn,                &
               negoffset_pn, offset_indx)
          call retrieveTBMEindx(orba, orbb, orbc, orbd,      &
               j2a, j2b, j2c, j2d, jab, jcd, Jop, indx)
          !
       else
          print*, 'invalid Tz', Tzop
          call cancelall(100)
          !
       end if
       !
       indx = indx + TBME_offset_pn(offset_indx)
       TBME_pn(indx) = matel * sqrt(2.d0*jab + 1.d0)
       !
    end do
    !
    ! Only for  Tz=0 operators
    !  < nn || Tz=0 || nn > matrix elements -- see module_TBME_Tz0
    !
    return
  end subroutine readTBMEascii_15200

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine readTBMEbin_15200(fh)
   use SPbasis, only: OrbitalList_t
   implicit none
   integer, intent(in) :: fh
   ! local variables
   integer, dimension(-1:1) :: J2max
   real(kind=4), dimension(-1:1) :: WT2max
   real(kind=4) :: wtab, wtcd
   type(OrbitalList_t) :: orblist
   integer, allocatable, dimension(:) :: tz2_orb
   integer :: versionnumber, i, indx, offset_indx
   integer :: Jab, Jcd, gab, gcd, parab, parcd
   integer :: orblist_norbt, orbc_min, orbc_max
   integer :: orba, orbb, orbc, orbd, j2a, j2b, j2c, j2d, noa, nob, noc, nod
   integer :: orbab, orbcd, TPorbab, TPorbcd
   integer :: Jop_f, gop_f, Tzop_f
   integer(kind=8) :: nTBMEs_pp_i8, nTBMEs_pn_i8, nTBMEs_nn_i8
   integer :: tmp_index
   real(kind=4), allocatable, dimension(:) :: tmp_mat
   character(len=MAX_PATHLENGTH) :: filename
   logical, external :: pairwiseless
   real, external :: myphase
   !
   ! reopen file with stream access, and throw away version header
   !
   inquire(unit=fh, name=filename)
   close(unit=fh, status='keep')
   open(unit=fh, file=TRIM(filename), &
        status='old', action='read', access='stream', form='unformatted')
   read(fh) i, versionnumber, i
   !
   ! check that we haven't arrived here by mistake
   !
   if (versionnumber /= 15200) then
      call cancelall(101)
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read orbital info
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! get number of orbitals
   read(fh) orblist%norb_p, orblist%norb_n
   orblist_norbt = orblist%norb_p + orblist%norb_n

   ! allocate arrays for orbital quantum numbers
   allocate(orblist%j2_orb(orblist_norbt), orblist%pr_orb(orblist_norbt))
   allocate(orblist%n_orb(orblist_norbt), orblist%l_orb(orblist_norbt))
   allocate(tz2_orb(orblist_norbt), orblist%wt_orb(orblist_norbt))

   ! get orbital quantum numbers
   read(fh) orblist%n_orb(1:orblist_norbt)
   read(fh) orblist%l_orb(1:orblist_norbt)
   read(fh) orblist%j2_orb(1:orblist_norbt)
   read(fh) tz2_orb(1:orblist_norbt)
   read(fh) orblist%wt_orb(1:orblist_norbt)

   ! assign orbital parities from orbital am l
   do i = 1, orblist_norbt
      orblist%pr_orb(i) = (-1)**orblist%l_orb(i)
   enddo

   ! check that orbitals are sorted by tz and match p and n counts
   do i = 1, orblist%norb_p
      if (tz2_orb(i) /= 1) then
         print*, 'not a proton orbital', i,                         &
            orblist%n_orb(i), orblist%l_orb(i),             &
            orblist%j2_orb(i), tz2_orb(i), orblist%wt_orb(i)
         call cancelall(104)
      endif
   enddo
   do i = orblist%norb_p+1, orblist_norbt
      if (tz2_orb(i) /= -1) then
         print*, 'not a neutron orbital', i,                        &
            orblist%n_orb(i), orblist%l_orb(i),             &
            orblist%j2_orb(i), tz2_orb(i), orblist%wt_orb(i)
         call cancelall(105)
      endif
   enddo

   ! check that orbitals are a subset of SPbasis orbitals
   if(allocated(j2_orb)) then
      if (orblist%norb_p < norb_p) then
         print*, 'insufficient norb_p', orblist%norb_p, norb_p
         call cancelall(001)
      endif
      if (orblist%norb_n < norb_n) then
         print*, 'insufficient norb_n', orblist%norb_n, norb_n
         call cancelall(002)
      endif
      do i = 1, norb_p
         if ((orblist%n_orb(i) /= n_orb(i)) .or.         &
             (orblist%l_orb(i) /= l_orb(i)) .or.         &
             (orblist%j2_orb(i) /= j2_orb(i)) .or.       &
             (orblist%pr_orb(i) /= pr_orb(i)) .or.       &
             (abs(orblist%wt_orb(i) - wt_orb(i)) > 1e-4) &
         ) then
            print*, 'mismatched proton orbitals'
            print*, i, orblist%n_orb(i), orblist%l_orb(i),  &
               orblist%j2_orb(i), orblist%pr_orb(i),        &
               orblist%wt_orb(i)
            print*, i, n_orb(i), l_orb(i), j2_orb(i), pr_orb(i), wt_orb(i)
            call cancelall(003)
         endif
      enddo
      offset_indx = orblist%norb_p
      do i = 1, norb_n
         if ((orblist%n_orb(i+offset_indx) /= n_orb(i+norb_p)) .or.   &
             (orblist%l_orb(i+offset_indx) /= l_orb(i+norb_p)) .or.   &
             (orblist%j2_orb(i+offset_indx) /= j2_orb(i+norb_p)) .or. &
             (orblist%pr_orb(i+offset_indx) /= pr_orb(i+norb_p)) .or. &
             (abs(orblist%wt_orb(i+offset_indx)-wt_orb(i+norb_p)) > 1e-4) &
         ) then
            print*, 'mismatched neutron orbitals'
            print*, i+offset_indx,                &
               orblist%n_orb(i+offset_indx),  &
               orblist%l_orb(i+offset_indx),  &
               orblist%j2_orb(i+offset_indx), &
               orblist%pr_orb(i+offset_indx), &
               orblist%wt_orb(i+offset_indx)
            print*, i+norb_p, n_orb(i+norb_p), l_orb(i+norb_p), &
               j2_orb(i+norb_p), pr_orb(i+norb_p), wt_orb(i+norb_p)
            call cancelall(004)
         endif
      enddo
      !
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read operator info
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! operator quantum numbers
   read(fh) Jop_f, gop_f, Tzop_f
   if (Jop == -1) then
      Jop = Jop_f
      gop = gop_f
      Tzop = Tzop_f
      Parop = (-1)**gop
   else
      if ((Jop_f/=Jop).or.(gop_f/=gop).or.(Tzop_f/=Tzop)) then
         print*, 'incompatible operator quantum numbers'
         print*, Jop, gop, Tzop
         print*, Jop_f, gop_f, Tzop_f
         call cancelall(005)
      endif
   endif

   if (Tzop .eq. 0) then
      print*, 'invalid Tzop', Tzop
      call cancelall(100)
   endif

   ! two-body truncation weights
   !   WT2max is indexed by Tz, i.e. 1->pp, 0->pn, -1->nn
   read(fh) WT2max(1:-1:-1)

   if (WT2max_pp == -1) then
      WT2max_pp = WT2max(1)
      WT2max_pn = WT2max(0)
      WT2max_nn = WT2max(-1)
   else
      if ((WT2max(1) < WT2max_pp) .or. (WT2max(0) < WT2max_pn) &
            .or. (WT2max(-1) < WT2max_nn)) &
      then
         print*, 'tbme file too small'
         print*, WT2max(1), WT2max(0), WT2max(-1)
         print*, WT2max_pp, WT2max_pn, WT2max_nn
         call cancelall(006)
      endif
   endif

   ! number of input TBMEs
   read(fh) nTBMEs_pp_i8, nTBMEs_pn_i8, nTBMEs_nn_i8  ! used as consistency check

   !
   if (.not.allocated(TPorb_pp)) then
      ! indexing only supports up to 2**31 - 1 TBMEs
      if (nTBMEs_pp_i8 >= 2_8**31) then
         print*, 'too many pp TBMEs', nTBMEs_pp_i8
         call cancelall(007)
      else
         nTBMEs_pp = nTBMEs_pp_i8
      endif
      if (nTBMEs_pn_i8 >= 2_8**31) then
         print*, 'too many pn TBMEs', nTBMEs_pn_i8
         call cancelall(008)
      else
         nTBMEs_pn = nTBMEs_pn_i8
      endif
      if (nTBMEs_nn_i8 >= 2_8**31) then
         print*, 'too many nn TBMEs', nTBMEs_nn_i8
         call cancelall(009)
      else
         nTBMEs_nn = nTBMEs_nn_i8
      endif

      ! Set up TwoBody indexing arrays
      call setup_TBME_indexing
   else
      if (nTBMEs_pp_i8 < TBME_offset_pp(nTBMEindx_pp+1)) then
         print*, 'input pp tbme too small'
         print*, nTBMEs_pp, TBME_offset_pp(nTBMEindx_pp+1)
         call cancelall(011)
      endif
      !
      if ((Tzop .eq. 1) .and. (nTBMEs_pn_i8 < TBME_offset_pn(nTBMEindx_pn+1))) then
         print*, 'input pn tbme too small'
         print*, nTBMEs_pn, TBME_offset_pn(nTBMEindx_pn+1)
         call cancelall(012)
      endif
      !
   endif
   !

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read TBMEs <ab;Jab||O||cd;Jcd>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   J2max(1) = maxval(orblist%j2_orb(1:orblist%norb_p))
   J2max(1) = J2max(1) - mod(J2max(1),2)  ! only even J for like particles

   J2max(0) = (maxval(orblist%j2_orb(1:orblist%norb_p))       &
              + maxval(orblist%j2_orb(orblist%norb_p+1:orblist_norbt)))/2

   J2max(-1) = maxval(orblist%j2_orb(orblist%norb_p+1:orblist_norbt))
   J2max(-1) = J2max(-1) - mod(J2max(-1),2)  ! only even J for like particles

   ! pp
   if (Tzop .eq. 1) then
      ! orbc is proton
      orbc_min = 1
      orbc_max = orblist%norb_p
   elseif (Tzop .eq. 2) then
      ! orbc is neutron
      orbc_min = orblist%norb_p+1
      orbc_max = orblist_norbt
   else
      print*, 'invalid Tz', Tzop
      call cancelall(100)
   endif
   ! read matrix elements into temporary array
   allocate(tmp_mat(nTBMEs_pp_i8))
   read(fh) tmp_mat(1:nTBMEs_pp_i8)
   tmp_index = 0
   do Jab = 0, J2max(1)
      do gab = 0, 1
         parab = (-1)**gab
         ! parcd is uniquely constrained by operator quantum number
         parcd = parab*Parop
         ! loop over ket J-subspaces
         do Jcd = abs(Jab-Jop), min((Jab+Jop),J2max(1-Tzop))
            ! loop over states in bra subspace
            do orba = 1, orblist%norb_p
               do orbb = orba, orblist%norb_p
                  if (orblist%pr_orb(orba)*orblist%pr_orb(orbb) /= parab) cycle ! parity
                  if ((orba == orbb) .and. (mod(Jab,2) /= 0)) cycle    ! antisymmetry
                  j2a = orblist%j2_orb(orba)
                  j2b = orblist%j2_orb(orbb)
                  if ((abs(j2a-j2b)/2 > Jab).or.(Jab > (j2a+j2b)/2)) cycle ! triangularity
                  wtab = orblist%wt_orb(orba)+orblist%wt_orb(orbb)
                  if (wtab > WT2max(1)) cycle ! WTmax

                  ! loop over states in ket subspace
                  do orbc = orbc_min, orbc_max
                     do orbd = max(orbc,orblist%norb_p+1), orblist_norbt
                        if (orblist%pr_orb(orbc)*orblist%pr_orb(orbd) /= parcd) cycle ! parity
                        if ((orbc == orbd) .and. (mod(Jcd,2) /= 0)) cycle    ! antisymmetry
                        j2c = orblist%j2_orb(orbc)
                        j2d = orblist%j2_orb(orbd)
                        if ((abs(j2c-j2d)/2 > Jcd).or.(Jcd > (j2c+j2d)/2)) cycle ! triangularity
                        wtcd = orblist%wt_orb(orbc)+orblist%wt_orb(orbd)
                        if (wtcd > WT2max(1-Tzop)) cycle ! WTmax

                        ! increment counter
                        tmp_index = tmp_index + 1

                        ! skip matrix element if larger than target indexing
                        if ((orba .gt. norb_p) .or. (orbb .gt. norb_p)) cycle
                        if ((Tzop .eq. 1) .and. (orbc > norb_p)) cycle
                        if ((Tzop .eq. 2) .and. (orbc-orblist%norb_p > norb_n)) cycle
                        if (orbd-orblist%norb_p > norb_n) cycle
                        if (wtab > WT2max_pp) cycle
                        if ((Tzop .eq. 1) .and. (wtcd > WT2max_pn)) cycle
                        if ((Tzop .eq. 2) .and. (wtcd > WT2max_nn)) cycle

                        ! get target indexing
                        orbab = norb2_pp - (norb_p-orba+1)*(norb_p-orba+2)/2 + (orbb-orba) + 1
                        TPorbab = TPorb_pp(orbab)
                        if (Tzop .eq. 1) then
                           nod = orbd - orblist%norb_p  ! orbd is relative to file's orbitals
                           orbcd = (orbc-1)*norb_n + nod
                           TPorbcd = TPorb_pn(orbcd)
                           call retrieveTBMEoffset(parab, parcd,              &
                                TPorbab, nTPpos_pp, nTPneg_pp,                &
                                TPorbcd, nTPpos_pn, nTPneg_pn,                &
                                negoffset_pp, offset_indx)
                           call retrieveTBMEindx(orba, orbb, orbc, orbd,      &
                                j2a, j2b, j2c, j2d, jab, jcd, Jop, indx)
                           !
                        elseif (Tzop .eq. 2) then
                           noc = orbc - orblist%norb_p  ! orbc is relative to file's orbitals
                           nod = orbd - orblist%norb_p  ! orbd is relative to file's orbitals
                           orbcd = norb2_nn - (norb_n-noc+1)*(norb_n-noc+2)/2 + (nod-noc) + 1
                           TPorbcd = TPorb_nn(orbcd)
                           call retrieveTBMEoffset(parab, parcd,              &
                                TPorbab, nTPpos_pp, nTPneg_pp,                &
                                TPorbcd, nTPpos_nn, nTPneg_nn,                &
                                negoffset_pp, offset_indx)
                           call retrieveTBMEindx(orba, orbb, orbc, orbd,      &
                                j2a, j2b, j2c, j2d, jab, jcd, Jop, indx)
                           !
                        end if

                        indx = indx + TBME_offset_pp(offset_indx)

                        ! convert to Edmonds convention and store
                        TBME_pp(indx) = tmp_mat(tmp_index) * sqrt(2.d0*Jab + 1.d0)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   ! deallocate temporary array
   deallocate(tmp_mat)

   ! no pn matrix elements if Tzop is 2
   if (Tzop .eq. 2) return

   ! pn
   ! read matrix elements into temporary array
   allocate(tmp_mat(nTBMEs_pn_i8))
   read(fh) tmp_mat(1:nTBMEs_pn_i8)
   tmp_index = 0
   do Jab = 0, J2max(0)
      do gab = 0, 1
         parab = (-1)**gab
         ! parcd is uniquely constrained by operator quantum number
         parcd = parab*Parop
         ! loop over ket J-subspaces
         do Jcd = abs(Jab-Jop), min((Jab+Jop),J2max(-1))
            ! loop over states in bra subspace
            do orba = 1, orblist%norb_p
               do orbb = orblist%norb_n+1, orblist_norbt
                  if (orblist%pr_orb(orba)*orblist%pr_orb(orbb) /= parab) cycle ! parity
                  if ((orba == orbb) .and. (mod(Jab,2) /= 0)) cycle    ! antisymmetry
                  j2a = orblist%j2_orb(orba)
                  j2b = orblist%j2_orb(orbb)
                  if ((abs(j2a-j2b)/2 > Jab).or.(Jab > (j2a+j2b)/2)) cycle ! triangularity
                  wtab = orblist%wt_orb(orba)+orblist%wt_orb(orbb)
                  if (wtab > WT2max(1)) cycle ! WTmax

                  ! loop over states in ket subspace
                  do orbc = orblist%norb_n+1, orblist_norbt
                     do orbd = orbc, orblist_norbt
                        if (orblist%pr_orb(orbc)*orblist%pr_orb(orbd) /= parcd) cycle ! parity
                        if ((orbc == orbd) .and. (mod(Jcd,2) /= 0)) cycle    ! antisymmetry
                        j2c = orblist%j2_orb(orbc)
                        j2d = orblist%j2_orb(orbd)
                        if ((abs(j2c-j2d)/2 > Jcd).or.(Jcd > (j2c+j2d)/2)) cycle ! triangularity
                        wtcd = orblist%wt_orb(orbc)+orblist%wt_orb(orbd)
                        if (wtcd > WT2max(1-Tzop)) cycle ! WTmax

                        ! increment counter
                        tmp_index = tmp_index + 1

                        ! skip matrix element if larger than target indexing
                        if ((orba > norb_p) .or. (orbb-orblist%norb_p > norb_n)) cycle
                        if ((orbc-orblist%norb_p > norb_n) .or. (orbd-orblist%norb_p > norb_n)) cycle
                        if ((wtab > WT2max_pn) .or. (wtcd > WT2max_nn)) cycle

                        ! get target indexing
                        nob = orbb - orblist%norb_p  ! orbb is relative to file's orbitals
                        orbab = (orba-1)*norb_n + nob
                        TPorbab = TPorb_pn(orbab)
                        noc = orbc - orblist%norb_p  ! orbc is relative to file's orbitals
                        nod = orbd - orblist%norb_p  ! orbd is relative to file's orbitals
                        orbcd = norb2_nn - (norb_n-noc+1)*(norb_n-noc+2)/2 + (nod-noc) + 1
                        TPorbcd = TPorb_nn(orbcd)
                        call retrieveTBMEoffset(parab, parcd,              &
                              TPorbab, nTPpos_pn, nTPneg_pn,                &
                              TPorbcd, nTPpos_nn, nTPneg_nn,                &
                              negoffset_pn, offset_indx)
                        call retrieveTBMEindx(orba, orbb, orbc, orbd,      &
                              j2a, j2b, j2c, j2d, jab, jcd, Jop, indx)

                        indx = indx + TBME_offset_pn(offset_indx)

                        ! convert to Edmonds convention and store
                        TBME_pn(indx) = tmp_mat(tmp_index) * sqrt(2.d0*Jab + 1.d0)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   ! deallocate temporary array
   deallocate(tmp_mat)
   !
   return
end subroutine readTBMEbin_15200

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setup_TBME_indexing
   norb2_pp = norb_p * (norb_p+1)/2
   allocate(TPorb_pp(norb2_pp))
   call setTPOindx(0,      norb_p, 0,      norb_p, norb2_pp, WT2max_pp, &
        TPorb_pp, nTPpos_pp, nTPneg_pp, nTPorb_pp, IDN_pp)
   !
   norb2_pn = norb_p * norb_n
   allocate(TPorb_pn(norb2_pn))
   call setTPOindx(0,      norb_p, norb_p, norb_n, norb2_pn, WT2max_pn, &
        TPorb_pn, nTPpos_pn, nTPneg_pn, nTPorb_pn, IDN_pn)
   !
   norb2_nn = norb_n * (norb_n+1)/2
   allocate(TPorb_nn(norb2_nn))
   call setTPOindx(norb_p, norb_n, norb_p, norb_n, norb2_nn, WT2max_nn, &
        TPorb_nn, nTPpos_nn, nTPneg_nn, nTPorb_nn, IDN_nn)
   !
   if (Tzop .eq. 1) then
      ! convention/notation: _pp for < pp || Tz=1 || pn >
      ! convention/notation: _pn for < pn || Tz=1 || nn >
      if (Parop .eq. 1) then
         negoffset_pp = nTPpos_pp * nTPpos_pn
         nTBMEindx_pp = negoffset_pp + nTPneg_pp * nTPneg_pn
         negoffset_pn = nTPpos_pn * nTPpos_nn
         nTBMEindx_pn = negoffset_pn + nTPneg_pn * nTPneg_nn
      else
         negoffset_pp = nTPpos_pp * nTPneg_pn
         nTBMEindx_pp = negoffset_pp + nTPneg_pp * nTPpos_pn
         negoffset_pn = nTPpos_pn * nTPneg_nn
         nTBMEindx_pn = negoffset_pn + nTPneg_pn * nTPpos_nn
      endif
      !
      ! for < pp || Tz=1 || pn > matrix elements
      allocate(numTBME_pp(nTBMEindx_pp))
      allocate(TBME_offset_pp(nTBMEindx_pp+1))
      call setTBMEoffset(Jop, Parop, nTBMEindx_pp, negoffset_pp,        &
           0,      norb_p, 0,      norb_p, WT2max_pp, IDN_pp,           &
           0,      norb_p, norb_p, norb_n, WT2max_pn, IDN_pn,           &
           norb2_pp, TPorb_pp, nTPpos_pp, nTPneg_pp,                    &
           norb2_pn, TPorb_pn, nTPpos_pn, nTPneg_pn,                    &
           numTBME_pp, TBME_offset_pp)
      !
      if (nTBMEs_pp .lt. TBME_offset_pp(nTBMEindx_pp+1)) then
         print*, 'pp2pn error', nTBMEs_pp, TBME_offset_pp(nTBMEindx_pp+1)
         call cancelall(111)
      endif
      allocate(TBME_pp(nTBMEs_pp))
      !
      ! for < pn || Tz=1 || nn > matrix elements
      allocate(numTBME_pn(nTBMEindx_pn))
      allocate(TBME_offset_pn(nTBMEindx_pn+1))
      call setTBMEoffset(Jop, Parop, nTBMEindx_pn, negoffset_pn,        &
           0,      norb_p, norb_p, norb_n, WT2max_pn, IDN_pn,           &
           norb_p, norb_n, norb_p, norb_n, WT2max_nn, IDN_nn,           &
           norb2_pn, TPorb_pn, nTPpos_pn, nTPneg_pn,                    &
           norb2_nn, TPorb_nn, nTPpos_nn, nTPneg_nn,                    &
           numTBME_pn, TBME_offset_pn)
      !
      if (nTBMEs_pn .lt. TBME_offset_pn(nTBMEindx_pn+1)) then
         print*, 'pn2nn error', nTBMEs_pn, TBME_offset_pn(nTBMEindx_pn+1)
         call cancelall(112)
      endif
      allocate(TBME_pn(nTBMEs_pn))
      !
      ! no < nn || Tz=1 || ?? > matrix elements
      if (nTBMEs_nn .gt. 0) then
         print*, 'Tz = 1 error', nTBMEs_nn
         call cancelall(122)
      endif
      !
      !
   elseif (Tzop .eq. 2) then
      ! convention/notation: _pp for < pp || Tz=2 || nn >
      if (Parop .eq. 1) then
         negoffset_pp = nTPpos_pp * nTPpos_nn
         nTBMEindx_pp = negoffset_pp + nTPneg_pp * nTPneg_nn
      else
         negoffset_pp = nTPpos_pp * nTPneg_nn
         nTBMEindx_pp = negoffset_pp + nTPneg_pp * nTPpos_nn
      endif
      !
      ! for < pp || Tz=2 || nn > matrix elements
      allocate(numTBME_pp(nTBMEindx_pp))
      allocate(TBME_offset_pp(nTBMEindx_pp+1))
      call setTBMEoffset(Jop, Parop, nTBMEindx_pp, negoffset_pp,  &
           0,      norb_p, 0,      norb_p, WT2max_pp, IDN_pp,     &
           norb_p, norb_n, norb_p, norb_n, WT2max_nn, IDN_nn,     &
           norb2_pp, TPorb_pp, nTPpos_pp, nTPneg_pp,              &
           norb2_nn, TPorb_nn, nTPpos_nn, nTPneg_nn,              &
           numTBME_pp, TBME_offset_pp)
      !
      if (nTBMEs_pp .lt. TBME_offset_pp(nTBMEindx_pp+1)) then
         print*, 'pp2nn error', nTBMEs_pp, TBME_offset_pp(nTBMEindx_pp+1)
         call cancelall(211)
      endif
      allocate(TBME_pp(nTBMEs_pp))
      !
      ! no < pn || Tz=2 || ?? > matrix elements
      if (nTBMEs_pn .gt. 0) then
         print*, 'Tz = 2 error', nTBMEs_pn
         call cancelall(212)
      endif
      ! no < nn || Tz=2 || ?? > matrix elements
      if (nTBMEs_nn .gt. 0) then
         print*, 'Tz = 2 error', nTBMEs_nn
         call cancelall(222)
      endif
      !
   else
      print*, 'invalid Tz', Tzop
      call cancelall(100)
      !
   end if
end subroutine setup_TBME_indexing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setTPOindx(offseta, norba, offsetb, norbb, norb2, WT2max, &
       TPorb, nTPpos, nTPneg, nTPorb, IDN)
    logical, intent(in) :: IDN
    integer, intent(in) :: offseta, norba, offsetb, norbb, norb2
    real, intent(in) :: WT2max
    integer, dimension(norb2), intent(out) :: TPorb
    integer, intent(out) :: nTPpos, nTPneg, nTPorb
    ! local variables
    integer :: ipos, ineg, ia, orba, ib, ibstart, orbb, iab, parab
    real :: wtab
    !
    TPorb(1:norb2) = 0
    !
    ibstart = 1
    !
    ipos = 0
    ineg = 0
    do ia = 1, norba
       orba = offseta + ia
       if (IDN) ibstart = ia
       do ib = ibstart, norbb
          orbb = offsetb + ib
          !
          wtab = wt_orb(orba) + wt_orb(orbb)
          if (wtab .gt. WT2max) cycle
          !
          if (IDN) then
             iab = norb2 - (norba-ia+1)*(norba-ia+2)/2 + ib - ia + 1
          else
             iab = (ia-1) * norbb + ib
          endif
          !
          parab = pr_orb(orba) * pr_orb(orbb)
          if (parab .eq. 1) then
             ipos = ipos + 1
             TPorb(iab) = ipos
          else
             ineg = ineg + 1
             TPorb(iab) = ineg
          endif
       enddo
    enddo
    !
    nTPpos = ipos
    nTPneg = ineg
    nTPorb = ipos+ineg
    !
    return
  end subroutine setTPOindx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setTBMEoffset(Jop, Parop, nTBMEindx, negoffset,                &
       offseta, norba, offsetb, norbb, WT_ab_max, IDN_ab,                   &
       offsetc, norbc, offsetd, norbd, WT_cd_max, IDN_cd,                   &
       norb2_ab, TPorb_ab, nTPpos_ab, nTPneg_ab,                            &
       norb2_cd, TPorb_cd, nTPpos_cd, nTPneg_cd,                            &
       numTBME, TBME_offset)
    logical, intent(in) :: IDN_ab, IDN_cd
    integer, intent(in) :: Jop, Parop, nTBMEindx, negoffset
    integer, intent(in) :: offseta, norba, offsetb, norbb, norb2_ab
    integer, intent(in) :: offsetc, norbc, offsetd, norbd, norb2_cd
    real, intent(in) :: WT_ab_max, WT_cd_max
    integer, dimension(norb2_ab), intent(in) :: TPorb_ab
    integer, dimension(norb2_cd), intent(in) :: TPorb_cd
    integer, intent(in) :: nTPpos_ab, nTPneg_ab
    integer, intent(in) :: nTPpos_cd, nTPneg_cd
    integer, dimension(nTBMEindx), intent(out) :: numTBME
    integer, dimension(nTBMEindx+1), intent(out) :: TBME_offset
    ! local variables
    integer :: ipos, ineg, indx, i, iastart, ibstart, icstart, idstart
    integer :: ia, ib, id, ic, orba, orbb, orbc, orbd,  j2a, j2b, j2c, j2d
    integer :: iab, icd, jab, jabmin, jabmax, jcd, jcdmin, jcdmax, parab, parcd
    real :: wtab, wtcd
    !
    iastart = 1
    ibstart = 1
    icstart = 1
    idstart = 1
    !
    ipos = 0
    ineg = 0
    !
    do ia = iastart, norba
       orba = offseta + ia
       !
       if (IDN_ab) ibstart = ia
       do ib = ibstart, norbb
          orbb = offsetb + ib
          !
          wtab = wt_orb(orba) + wt_orb(orbb)
          if (wtab .gt. WT_ab_max) cycle
          !
          if (IDN_ab) then
             iab = norb2_ab - (norba-ia+1)*(norba-ia+2)/2 + ib - ia + 1
          else
             iab = (ia-1) * norbb + ib
          endif
          !
          j2a = j2_orb(orba)
          j2b = j2_orb(orbb)
          jabmin = abs(j2a-j2b)/2
          jabmax = (j2a+j2b)/2
          if (IDN_ab .and. (orba.eq.orbb)) jabmax = jabmax - 1
          !
          parab = pr_orb(orba) * pr_orb(orbb)
          !
          do ic = icstart, norbc
             orbc = offsetc + ic
             !
             if (IDN_cd) idstart = ic
             do id = idstart, norbd
                orbd = offsetd + id
                !
                wtcd = wt_orb(orbc) + wt_orb(orbd)
                if (wtcd .gt. WT_cd_max) cycle
                !
                if (IDN_cd) then
                   icd = norb2_cd - (norbc-ic+1)*(norbc-ic+2)/2 + id - ic + 1
                else
                   icd = (ic-1) * norbd + id
                endif
                !
                j2c = j2_orb(orbc)
                j2d = j2_orb(orbd)
                jcdmin = abs(j2c-j2d)/2
                jcdmax = (j2c+j2d)/2
                if (IDN_cd .and. (orbc.eq.orbd)) jcdmax = jcdmax - 1
                !
                parcd = pr_orb(orbc) * pr_orb(orbd)
                !
                indx = 0
                if (parab*parcd .eq. Parop) then
                   do jab = jabmin, jabmax
                      if (IDN_ab .and. (orba.eq.orbb)       &
                           .and. (mod(jab,2).eq.1)) cycle
                      do jcd = max(jab-Jop, jcdmin), min(jab+Jop, jcdmax)
                         if (IDN_cd .and. (orbc.eq.orbd)    &
                              .and. (mod(jcd,2).eq.1)) cycle
                         if (abs(jab-jcd) .gt. Jop) cycle
                         if ((jab+jcd) .lt. Jop) cycle
                         indx = indx + 1
                      enddo
                   enddo
                   if ((parab .eq. 1).and.(parcd .eq. 1)) then
                      i = (TPorb_ab(iab)-1) * nTPpos_cd + TPorb_cd(icd)
                   elseif ((parab .eq. -1).and.(parcd .eq. -1)) then
                      i = negoffset +    &
                           (TPorb_ab(iab)-1) * nTPneg_cd + TPorb_cd(icd)
                   elseif ((parab .eq. 1).and.(parcd .eq. -1)) then
                      i = (TPorb_ab(iab)-1) * nTPneg_cd + TPorb_cd(icd)
                   elseif ((parab .eq. -1).and.(parcd .eq. 1)) then
                      i = negoffset +    &
                           (TPorb_ab(iab)-1) * nTPpos_cd + TPorb_cd(icd)
                   endif
                   numTBME(i) = indx
                endif
                !
             enddo
          enddo
       enddo
    enddo
    !
    TBME_offset(1) = 0
    do i = 1, nTBMEindx
       TBME_offset(i+1) = TBME_offset(i) + numTBME(i)
       nTBMEs_max = max(nTBMEs_max, numTBME(i))
    enddo
    !
    return
  end subroutine setTBMEoffset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine retrieveTBMEoffset(parab, parcd,                              &
       TPorbab, nTPpos_ab, nTPneg_ab, TPorbcd, nTPpos_cd, nTPneg_cd,       &
       negoffset, indx)
    integer, intent(in) :: parab, parcd, negoffset
    integer, intent(in) :: TPorbab, nTPpos_ab, nTPneg_ab
    integer, intent(in) :: TPorbcd, nTPpos_cd, nTPneg_cd
    integer, intent(out) :: indx
    !
    if ((TPorbab .eq. 0) .or. (TPorbcd .eq. 0)) then
       indx = 0
       return
    endif
    if ( (parab.eq.1) .and. (parcd.eq.1) ) then
       indx = (TPorbab - 1) * nTPpos_cd + TPorbcd
    elseif ((parab.eq.-1).and. (parcd.eq.-1)) then
       indx = negoffset + (TPorbab - 1) * nTPneg_cd + TPorbcd
    elseif ((parab.eq.1) .and. (parcd.eq.-1)) then
       indx = (TPorbab - 1) * nTPneg_cd + TPorbcd
    elseif ((parab.eq.-1).and. (parcd.eq.1) ) then
       indx = negoffset + (TPorbab - 1) * nTPpos_cd + TPorbcd
    endif
    !
    return
  end subroutine retrieveTBMEoffset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine retrieveTBMEindx(orba, orbb, orbc, orbd,        &
       j2a, j2b, j2c, j2d, jab, jcd, Jop, indx)
    integer, intent(in) :: orba, orbb, orbc, orbd
    integer, intent(in) :: j2a, j2b, j2c, j2d, jab, jcd, Jop
    integer, intent(out) :: indx
    ! local variables
    logical :: IDN_ab, IDN_cd, DIAG
    integer :: jabmin, jabmax, jabstep, iab
    integer :: jcdmin, jcdmax, jcdstep, icd, icd_start, icd_stop
    !
    IDN_ab = orba .eq. orbb
    IDN_cd = orbc .eq. orbd
    DIAG  = (orba.eq.orbc) .and. (orbb.eq.orbd)
    !
    jabmin = abs(j2a-j2b)/2
    jabmax = (j2a+j2b)/2
    jabstep = 1
    if (IDN_ab) then
       jabmin = jabmin + mod(jabmin,2)
       jabmax = jabmax - mod(jabmax,2)
       jabstep = 2
    endif
    !
    jcdmin = abs(j2c-j2d)/2
    jcdmax = (j2c+j2d)/2
    jcdstep = 1
    if (IDN_cd) jcdstep = 2
    !
    indx = 0
    do iab = jabmin, jabmax, jabstep
       if (DIAG) then
          icd_start = max(abs(iab-Jop), iab)
       else
          icd_start = max(abs(iab-Jop), jcdmin)
       endif
       icd_stop  = min(iab+Jop, jcdmax)
       if (IDN_cd) then
          icd_start = icd_start + mod(icd_start,2)
          icd_stop = icd_stop - mod(icd_stop,2)
       endif
       do icd = icd_start, icd_stop, jcdstep
          !
          indx = indx + 1
          if ((iab.eq.jab).and.(icd.eq.jcd)) return
       end do
    end do
    !
    indx = -1
    return
  end subroutine retrieveTBMEindx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Set_TBME_array(parab, parcd,                  &
       orba, orbb, orbc, orbd, nops, nTBMEs, TBMEarray)
    !
    ! single-orbital indices aa =< bb and cc =< dd ; (aa, bb) =< (cc, dd)
    !
    real(kind=4), dimension(nops, nTBMEs), intent(out) :: TBMEarray
    integer, intent(in) :: nops, parab, parcd, orba, orbb, orbc, orbd
    integer, intent(inout) :: nTBMEs
    !
    ! local variables
    integer :: orbab, TPorbab, orbcd, TPorbcd, ii, offset
    logical, external :: pairwiseless
    !
    if (orbb .le. norb_p) then
       ! pp
       orbab = norb2_pp - (norb_p-orba+1)*(norb_p-orba+2)/2 + (orbb-orba) + 1
       TPorbab = TPorb_pp(orbab)
       !
       if (orbc .le. norb_p) then
          ! pp2pn
          orbcd = (orbc-1)*norb_n + (orbd - norb_p)
          TPorbcd = TPorb_pn(orbcd)
          !
          call retrieveTBMEoffset(parab, parcd,                              &
               TPorbab, nTPpos_pp, nTPneg_pp, TPorbcd, nTPpos_pn, nTPneg_pn, &
               negoffset_pp, ii)
       else
          ! pp2nn
          orbcd = norb2_nn - &
               (norb_n+norb_p-orbc+1)*(norb_n+norb_p-orbc+2)/2 + (orbd-orbc) + 1
          TPorbcd = TPorb_nn(orbcd)
          !
          call retrieveTBMEoffset(parab, parcd,                              &
               TPorbab, nTPpos_pp, nTPneg_pp, TPorbcd, nTPpos_nn, nTPneg_nn, &
               negoffset_pp, ii)
       endif
       if (ii .eq. 0) then
          offset = 0
          nTBMEs = 0
       else
          offset = TBME_offset_pp(ii)
          nTBMEs = TBME_offset_pp(ii+1) - offset
          TBMEarray(1:nops, 1:nTBMEs) = TBMEarray_pp(1:nops, offset+1:offset+nTBMEs)
       endif
       !
    else
       ! pn
       orbab = (orba-1)*norb_n + (orbb - norb_p)
       TPorbab = TPorb_pn(orbab)
       !
       ! pn2nn
       orbcd = norb2_nn - &
            (norb_n+norb_p-orbc+1)*(norb_n+norb_p-orbc+2)/2 + (orbd-orbc) + 1
       TPorbcd = TPorb_nn(orbcd)
       !
       call retrieveTBMEoffset(parab, parcd,                              &
            TPorbab, nTPpos_pn, nTPneg_pn, TPorbcd, nTPpos_nn, nTPneg_nn, &
            negoffset_pn, ii)
       !
       if (ii .eq. 0) then
          offset = 0
          nTBMEs = 0
       else
          offset = TBME_offset_pn(ii)
          nTBMEs = TBME_offset_pn(ii+1) - offset
          TBMEarray(1:nops, 1:nTBMEs) = TBMEarray_pn(1:nops, offset+1:offset+nTBMEs)
       endif
       !
    endif
    !
    return
  end subroutine Set_TBME_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module TBME_Tz12
