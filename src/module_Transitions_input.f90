
module Transition_input
  implicit none
  !
  ! Basis parameters
  real(4) :: fmass = 938.92, hbomeg = 0.0
  !
  ! Number of transition operators
  integer :: numTBtrans = 0
  !
  ! Transition OBDMEs
  logical :: obdme = .true.
  integer :: max2K = 4
  !
  ! <bra| and |ket>  state
  integer :: TwoJ_bra
  integer :: n_bra = 1
  integer, dimension(MaxNumKets) :: TwoJ_ket = -1
  integer, dimension(MaxNumKets) :: n_ket = 0
  !
  character(LEN=MAX_PATHLENGTH) :: infofilename_bra = 'mfdn_smwf.info'
  character(LEN=MAX_PATHLENGTH) :: infofilename_ket = 'mfdn_smwf.info'
  character(LEN=MAX_PATHLENGTH) :: basisfilename_bra = 'mfdn_MBgroups'
  character(LEN=MAX_PATHLENGTH) :: basisfilename_ket = 'mfdn_MBgroups'
  character(LEN=MAX_PATHLENGTH) :: smwffilename_bra = 'mfdn_smwf'
  character(LEN=MAX_PATHLENGTH) :: smwffilename_ket = 'mfdn_smwf'
  !
  character(LEN=MAX_PATHLENGTH), dimension(MaxNumTBMEops) :: TBMEoperators
  !
end module Transition_input
