!
!     read multiple TBME files for 2-body observables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_TBME_multi(nTBME, TBME_opfile)
#ifdef DeltaTz
  use TBME_Tz12, only: nTBMEs_pp, nTBMEs_nn, nTBMEs_pn, TBME_pp, TBME_nn, TBME_pn
  use TBME_Tz12, only: TBMEarray_pp, TBMEarray_pn, TBMEarray_nn
  use TBME_Tz12, only: readTBMEascii, readTBMEbin
#else
  use TBME_Tz0, only: nTBMEs_pp, nTBMEs_nn, nTBMEs_pn, TBME_pp, TBME_nn, TBME_pn
  use TBME_Tz0, only: TBMEarray_pp, TBMEarray_pn, TBMEarray_nn
  use TBME_Tz0, only: readTBMEascii, readTBMEbin
#endif
  implicit none
  !
  integer, intent(in) :: nTBME
  character(LEN=*), dimension(nTBME), intent(in) :: TBME_opfile
  !
  !     local variables
  integer :: iTB
  character(LEN=:), allocatable :: TBMEfile
  logical :: binary
  !
  do iTB = 1, nTBME
     TBMEfile = TBME_opfile(iTB)
     inquire(file=TRIM(TBMEfile)//'.bin', exist=binary)
     !
     if (binary) then
        !     actual binary read statements
        call readTBMEbin(TBMEfile)
     else
        !     actual ascii read statements
        call readTBMEascii(TBMEfile)
     endif
     if (iTB .eq. 1) then
        if (nTBMEs_pp.gt.0) allocate(TBMEarray_pp(nTBME, nTBMEs_pp))
        if (nTBMEs_pn.gt.0) allocate(TBMEarray_pn(nTBME, nTBMEs_pn))
        if (nTBMEs_nn.gt.0) allocate(TBMEarray_nn(nTBME, nTBMEs_nn))
     endif
     !
     if (nTBMEs_pp.gt.0) TBMEarray_pp(iTB, 1:nTBMEs_pp) = TBME_pp(1:nTBMEs_pp)
     if (nTBMEs_pn.gt.0) TBMEarray_pn(iTB, 1:nTBMEs_pn) = TBME_pn(1:nTBMEs_pn)
     if (nTBMEs_nn.gt.0) TBMEarray_nn(iTB, 1:nTBMEs_nn) = TBME_nn(1:nTBMEs_nn)
  enddo
  !
  return
end subroutine read_TBME_multi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
