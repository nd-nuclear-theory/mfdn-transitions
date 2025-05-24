!
!     read multiple TBME files for 2-body observables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_TBME_multi(nTBME, TBME_opfile)
  use TBME, only: ntbme_pp, ntbme_nn, ntbme_pn, read_TBME_ascii, read_TBME_bin, &
       H2full_pp, H2full_nn, H2full_pn, TBMEfull_pp, TBMEfull_nn, TBMEfull_pn
  implicit none
  !
  integer, intent(in) :: nTBME
  character(LEN=128), dimension(nTBME), intent(in) :: TBME_opfile
  !
  !     local variables
  integer :: iTB, i
  character(LEN=14) :: outputname= 'mfdn_trans.res'
  character(LEN=128) :: hamfile
  logical :: binary
  !
  open(file=outputname, unit=10, position='append', action='write')
1 format('', A)
8 format('', A,' = ',i8)
  write(10, 1) '[Observables]'
  write(10, 8) '  numTBops   ', nTBME
  !     
  do iTB = 1, nTBME
     hamfile = TBME_opfile(iTB)
     inquire(file=TRIM(hamfile)//'.bin', exist=binary) 
     !
     if (binary) then
        write(10, 1) 'TBMEfile(', iTB, TRIM(hamfile)//'.bin'
        !     actual binary read statements
        call read_TBME_bin(hamfile)
     else
        write(10, 1) 'TBMEfile(', iTB, TRIM(hamfile)//'.dat'
        !     actual ascii read statements
        call read_TBME_ascii(hamfile) 
     endif
     if (iTB .eq. 1) then
        allocate(TBMEfull_pp(nTBME, ntbme_pp))
        allocate(TBMEfull_nn(nTBME, ntbme_nn))
        allocate(TBMEfull_pn(nTBME, ntbme_pn))
     endif
     !
     TBMEfull_pp(iTB, 1:ntbme_pp) = H2full_pp(1:ntbme_pp)
     TBMEfull_nn(iTB, 1:ntbme_nn) = H2full_nn(1:ntbme_nn)
     TBMEfull_pn(iTB, 1:ntbme_pn) = H2full_pn(1:ntbme_pn)
  enddo
  write(10,*)
  close(unit=10, status='keep')
  !
  return
end subroutine read_TBME_multi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
