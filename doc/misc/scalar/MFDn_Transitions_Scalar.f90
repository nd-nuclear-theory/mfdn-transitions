!
! Calculates 
!
! (1) if (obdme) Reduced One-Body Transition Densities
!
! (2) if (numTBtrans > 0) numTBtrans Two-Body Transition matrix elements
!  (NOTE: initial implementation only for scalar Two-Body operators)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program MFDn_Transitions
  !use MPI
  use omp_lib
  use SPbasis, only: readSPbasis, nclasses, classoffset, nspstates, norbt, nparticles, j2_orb
  use Sparsity, only: generateSparsity, coltileptr, tilediff, rowind, ntiles
  use Wigner3J, only: Init_Wigner3J
  use TBME, only: Finalize_H2full
  use Transition_input
  ! use OBobservables, only: CalcOBobs
  implicit none
  !
  integer :: numthreads, root, myrank, nprocs
  !
  integer(kind=8) :: coldim, rowdim
  integer :: oprank, statesize
  integer :: TwoMj_col, colparity, ncolgroupids, numcolstates, ncolblks
  integer :: TwoMj_row, rowparity, nrowgroupids, numrowstates, nrowblks
  integer :: colrank, ncolprocs, rowrank, nrowprocs
  integer :: maxnMstates_col, maxnMstates_row, maxnMstates
  integer :: Delta_Tz, Delta_TwoMj, Delta_TwoJ, Delta_Par
  integer :: i, j, ifile, jfile, fh, indxNm_col, indxNm_row, Jtot1mx, Jtot2mx
  logical :: input
  !
  integer :: mnK, mxK, nobdme, nobdme_p, nobdme_n
  integer, dimension(:,:), allocatable :: obdmeK_ptr,obdmeOrbKbraket
  integer, dimension(:), allocatable :: obdmeOrb1offset, obdmeOrb2offset, obdmeOrbKindx
  real(kind=8), dimension(:), allocatable :: obdmeKvals, robdmes, reducefac
  real(kind=8) :: totR20, totM1L, totM1S, totE2Q
  !
  integer(kind=2), dimension(:,:), allocatable :: colgroupidlist
  integer(kind=2), dimension(:,:), allocatable :: rowgroupidlist
  integer, dimension(:), allocatable :: colMstateptr, colblkgroup, colblkoffset
  integer, dimension(:), allocatable :: rowMstateptr, rowblkgroup, rowblkoffset
  integer, dimension(:), allocatable :: ncolblksNm, nrowblksNm
  !
  real(kind=4), dimension(:), allocatable :: colamp, rowamp
  real(kind=8), dimension(:), allocatable :: transobs, transition_observables
  real(kind=8) :: colsum, rowsum
  ! File names
  character(LEN=15) :: outfilename = 'transitions.out'
  character(LEN=17) :: datfilename = 'transitions.input'
  character(LEN=131):: mbfile, smwffile
  !
  namelist /transition_data/ fmass, hbomeg, numTBtrans, obdme, max2K, &
       Zprotons_bra, Nneutrons_bra, Zprotons_ket, Nneutrons_ket, &
       brastate, ketstate, TwoJ_bra, TwoJ_ket, &
       colbasisfilename, rowbasisfilename, &
       colsmwffilename, rowsmwffilename, &
       TBMEoperators

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Initialize MPI
  !call MPI_Init_Thread(MPI_THREAD_FUNNELED, multithreaded, ierr)
  !call MPI_Comm_Size(MPI_COMM_WORLD, nprocs, ierr)
  !call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  root = 0
  nprocs = 1
  myrank = 0
  !     
  ! Set/check number of threads
  !numthreads = omp_get_max_threads()
  !
  ! Reading input file
  inquire(file=datfilename, exist=input)
  if (input) then
     open(unit=11, file=datfilename, status='old', action='read')
     read(11, nml=transition_data)
     close(unit=11, status='keep')
  else
     print*, 'No inputfile...'
     call cancelall(000)
  endif
  !
  nparticles = Zprotons_bra + Nneutrons_bra
  if (nparticles .ne. Zprotons_ket + Nneutrons_ket) then
     print*, 'Different number of particles for bra and ket state'
     print*, Zprotons_bra, Nneutrons_bra
     print*, Zprotons_ket, Nneutrons_ket
     call cancelall(002)
  endif
  Delta_Tz = Zprotons_bra - Zprotons_ket
  if (Delta_Tz .gt. 2) then
     print*, 'Delta Tz larger two'
     print*, Zprotons_bra,  Zprotons_ket
     print*, Nneutrons_bra, Nneutrons_ket
     call cancelall(004)
  elseif (Delta_Tz .lt. 0) then
     print*, 'Delta Tz smaller than zero'
     print*, Zprotons_bra,  Zprotons_ket
     print*, Nneutrons_bra, Nneutrons_ket
     call cancelall(005)
  endif
  !
  ! Use the same SP basis for row (<bra|, psibar) and column (|ket>, psi)
  call readSPbasis
  !
  ! Read first column and row basis files, to determine ncolprocs and nrowprocs
  fh = 17
  ifile = 1
  mbfile = TRIM(colbasisfilename)//achar(ifile/100+48)&
       //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
  open(unit=fh, file=TRIM(mbfile), status='old', action='read', form='unformatted')
  call readMBgroupID_metadata(fh, nclasses, nparticles, colparity, TwoMj_col, &
       indxNm_col, ncolgroupids, numcolstates, ncolblks, ncolprocs)
  close(fh, status='keep')
  allocate(ncolblksNm(indxNm_col))    ! not used?
  !
  mbfile = TRIM(rowbasisfilename)//achar(ifile/100+48)&
       //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
  open(unit=fh, file=TRIM(mbfile), status='old', action='read', form='unformatted')
  call readMBgroupID_metadata(fh, nclasses, nparticles, rowparity, TwoMj_row, &
       indxNm_row, nrowgroupids, numrowstates, nrowblks, nrowprocs)
  close(fh, status='keep')
  allocate(nrowblksNm(indxNm_row))    ! not used?
  !
  Delta_TwoMj = TwoMj_col - TwoMj_row
  Delta_TwoJ = abs(TwoJ_ket - TwoJ_bra)
  Delta_Par = colparity * rowparity
  !
  ! Set up OBDME indices and arrays
  if (obdme) then
     if (max2K .lt. Delta_TwoJ) then
        obdme = .false.
     elseif (Delta_Tz .gt. 1) then
        obdme = .false.
     elseif (Delta_Tz .eq. 1) then
        ! Delta_Tz OBDME < p | ( a^\dagger \tilde{a} )^K | n >
        ! to be implemented
     else
        ! <p| (a^\dagger \tilde{a})^K |p> and <n| (a^\dagger \tilde{a})^K |n> 
        mnK = abs(Delta_TwoJ) / 2
        mxK = max2K / 2
        allocate(obdmeOrb1offset(0:mxK+1))
        allocate(obdmeOrb2offset(0:mxK+1))
        if (Delta_Tz .eq. 0) then
           call CountOBDMEs(mnK, mxK, Delta_Par, nobdme_p, nobdme_n, &
                obdmeOrb1offset, obdmeOrb2offset)
           nobdme = nobdme_p + nobdme_n
        endif
        !
        allocate(obdmeK_ptr(norbt, norbt))
        allocate(obdmeOrbKbraket(2, nobdme))
        allocate(obdmeOrbKindx(nobdme))
        call SetOBDMEindices(mnK, mxK, Delta_Par, obdmeOrb1offset, obdmeOrb2offset, &
             nobdme, obdmeK_ptr, obdmeOrbKbraket, obdmeOrbKindx)
        !
        allocate(obdmeKvals(nobdme), robdmes(nobdme))
        obdmeKvals(1:nobdme) = 0.0d0
        !
     endif
  endif
  !
  ! Read Two-Body matrix elements
  if (numTBtrans .gt. 0) then
     !
     ! Read Two-Body Matrix Elements (currently limited to 16 operators)
     call read_TBME_multi(numTBtrans, TBMEoperators)
     !
     oprank = 2
     allocate(transobs(numTBtrans), transition_observables(numTBtrans))
     transition_observables(1:numTBtrans) = 0.0d0
  else
     oprank = 1
  endif
  statesize = nparticles + oprank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank .eq. root) then
     open(unit=7, file=outfilename, status='unknown', action='write')
     write(7,*)
     write(7,*) '  Output from MFDn Transitions in M-scheme basis'
     write(7,*)
     write(7,*) ' Total number of particles ', nparticles
     write(7,*) ' Number of protons  <bra| and |ket>', Zprotons_bra, Zprotons_ket
     write(7,*) ' Number of neutrons <bra| and |ket>', Nneutrons_bra, Nneutrons_ket
     write(7,*)
     write(7,*) ' Total spin  J  <bra| and |ket>', 0.5*TwoJ_bra,  0.5*TwoJ_ket
     write(7,*) ' Projection M_J <bra| and |ket>', 0.5*TwoMj_row, 0.5*TwoMj_col 
     write(7,*) ' Parity         <bra| and |ket>', rowparity, colparity
     write(7,*)
     if (obdme) then
        write(7,*) ' Calculation Reduced One-Body Density Matrix Elements'
        write(7,*) '   Lower limit on K', mnK
        write(7,*) '   Upper limit on K', mxK
        if (Delta_Par .ne. 1) write(7,*) '   Parity changing ', Delta_Par
        if (Delta_Tz  .eq. 0) then
           write(7,*)
           write(7,*) '   Number of proton OBDMEs ', nobdme_p
           write(7,*) '   Number of neutron OBDMEs', nobdme_n
        elseif (Delta_Tz  .eq. 1) then
           write(7,*) '   Isospin changing', Delta_Tz
           write(7,*)
           write(7,*) '   Number of <p|O|n> OBDMEs', nobdme
        endif
        write(7,*)
     endif
     close(unit=7, status='keep')
  endif

  Jtot1mx = (maxval(j2_orb)+1) / 2
  Jtot2mx = maxval(j2_orb)
  call Init_Wigner3j(Jtot1mx, Jtot2mx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Loop over row (< bra |, psibar) and column (| ket >, psi)
  coldim = 0
  colsum = 0.d0
  do colrank = 0, ncolprocs-1
     !
     ! Read Many-body basis for column ( | ket >, psi)
     ifile = colrank+1
     mbfile = TRIM(colbasisfilename)//achar(ifile/100+48)&
          //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
     open(unit=fh, file=TRIM(mbfile), status='old', action='read', form='unformatted')
     !
     call readMBgroupID_metadata(fh, nclasses, nparticles, colparity, TwoMj_col, &
          indxNm_col, ncolgroupids, numcolstates, ncolblks, ncolprocs)
     !
     read(fh) ncolblksNm(1:indxNm_col)    ! not used?
     !
     allocate(colMstateptr(ncolgroupids+1))
     read(fh) colMstateptr(1:ncolgroupids+1)
     !
     allocate(colgroupidlist(nparticles, ncolgroupids))
     read(fh) colgroupidlist(1:nparticles, 1:ncolgroupids)
     !
     allocate(colblkgroup(ncolblks+1))
     read(fh) colblkgroup(1:ncolblks+1)
     !
     allocate(colblkoffset(ncolblks+1))
     read(fh) colblkoffset(1:ncolblks+1)
     !
     close(fh, status='keep')
     !
     ! Set max number of M-scheme states in largest group
     maxnMstates_col = 0
     do i = 1, ncolgroupids
        maxnMstates_col = max(maxnMstates_col, colMstateptr(i+1)-colMstateptr(i))
     enddo
     !
     ! Read wavefunction | ket >
     allocate(colamp(numcolstates))
     smwffile = TRIM(colsmwffilename)//achar(ifile/100+48)&
          //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
     open(unit=fh, file=TRIM(smwffile), status='unknown', action='read', form='unformatted')
     do i = 1, ketstate
        read(fh) colamp(1:numcolstates)
     enddo
     close(fh, status='keep')
     ! Test on norm
     do j = 1, numcolstates
        colsum = colsum + colamp(j)**2
     enddo
     coldim = coldim + numcolstates
     !
     allocate(coltileptr(ncolgroupids+1))
     !
     rowdim = 0
     rowsum = 0.0
     do rowrank = 0, nrowprocs-1
        !
        ! Read Many-body basis for row ( < bra |, psibar)
        jfile = rowrank+1
        mbfile = TRIM(rowbasisfilename)//achar(jfile/100+48)&
             //achar(mod(jfile/10,10)+48)//achar(mod(jfile,10)+48)
        open(unit=fh, file=TRIM(mbfile), status='old', action='read', form='unformatted')
        call readMBgroupID_metadata(fh, nclasses, nparticles, rowparity, TwoMj_row, &
             indxNm_row, nrowgroupids, numrowstates, nrowblks, nrowprocs)
        !
        read(fh) nrowblksNm(1:indxNm_row)
        !
        allocate(rowMstateptr(nrowgroupids+1))
        read(fh) rowMstateptr(1:nrowgroupids+1)
        !
        allocate(rowgroupidlist(nparticles, nrowgroupids))
        read(fh) rowgroupidlist(1:nparticles, 1:nrowgroupids)
        !
        allocate(rowblkgroup(nrowblks+1))
        read(fh) rowblkgroup(1:nrowblks+1)
        !
        allocate(rowblkoffset(nrowblks+1))
        read(fh) rowblkoffset(1:nrowblks+1)
        !
        close(fh, status='keep')
        !
        ! Set max number of M-scheme states in largest group
        maxnMstates_row = 0
        do i = 1, nrowgroupids
           maxnMstates_row = max(maxnMstates_row, rowMstateptr(i+1)-rowMstateptr(i))
        enddo
        maxnMstates = max(maxnMstates_col, maxnMstates_row)
        !
        ! Read wavefunction < bra |
        allocate(rowamp(numrowstates))
        smwffile = TRIM(rowsmwffilename)//achar(jfile/100+48)&
             //achar(mod(jfile/10,10)+48)//achar(mod(jfile,10)+48)
        open(unit=fh, file=TRIM(smwffile), status='unknown', action='read', form='unformatted')
        do i = 1, brastate
           read(fh) rowamp(1:numrowstates)
        enddo
        close(fh, status='keep')
        ! Test on norm
        do j = 1, numrowstates
           rowsum = rowsum + rowamp(j)**2
        enddo
        rowdim = rowdim + numrowstates
        !
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! NOTE: difference with MFDn -- TwoMj could be different for < bra | and | ket >

        call generateSparsity(statesize, nparticles, oprank, nspstates, classoffset(2), &
             Zprotons_bra, Nneutrons_bra, Zprotons_ket, Nneutrons_ket, &
             nrowgroupids, rowgroupidlist, ncolgroupids, colgroupidlist )
        !
        if (obdme) then
           call EvalOBDME(Zprotons_ket, Nneutrons_ket, Zprotons_bra, Nneutrons_bra, maxnMstates, &
                TwoMj_col, ncolgroupids, colgroupidlist, colMstateptr, numcolstates, colamp, &
                TwoMj_row, nrowgroupids, rowgroupidlist, rowMstateptr, numrowstates, rowamp, &
                coltileptr, ntiles, rowind, tilediff, mnK, mxK, nobdme, obdmeK_ptr, robdmes)
           obdmeKvals(1:nobdme) = obdmeKvals(1:nobdme) + robdmes(1:nobdme)
        endif
        !
        if (numTBtrans .gt. 0) then
           call EvalTBobs(Zprotons_ket, Nneutrons_ket, Zprotons_bra, Nneutrons_bra, maxnMstates, &
                TwoMj_col, ncolgroupids, colgroupidlist, colMstateptr, numcolstates, colamp, &
                TwoMj_row, nrowgroupids, rowgroupidlist, rowMstateptr, numrowstates, rowamp, &
                coltileptr, ntiles, rowind, tilediff, numTBtrans, transobs )
           transition_observables(1:numTBtrans) = & 
                transition_observables(1:numTBtrans) + transobs(1:numTBtrans)
        endif
        !
        deallocate(tilediff, rowind) ! allocated in generateSparsity
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(rowMstateptr)
        deallocate(rowgroupidlist)
        deallocate(rowblkgroup)
        deallocate(rowblkoffset)
        deallocate(rowamp)
        
     enddo
     
     deallocate(coltileptr)

     deallocate(colMstateptr)
     deallocate(colgroupidlist)
     deallocate(colblkgroup)
     deallocate(colblkoffset)
     deallocate(colamp)

  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (obdme) then
     allocate(reducefac(0:mxK))
     call ReduceOBDME(TwoJ_bra, TwoJ_ket, TwoMj_row, TwoMj_col, mnK, mxK, &
          Delta_Par, nobdme, obdmeKvals, obdmeK_ptr, obdmeOrbKindx, reducefac, robdmes)
     call WriteOBDME(mnK, mxK, Delta_Par, nobdme, nobdme_p, nobdme_n, &
          obdmeOrb1offset, obdmeOrb2offset, obdmeOrbKbraket, reducefac, robdmes)
     deallocate(reducefac)
  endif

  if (numTBtrans .gt. 0) then
     do i = 1, numTBtrans
        print*, i, TRIM(TBMEoperators(i)), transition_observables(i)
     enddo
  endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  deallocate(nrowblksNm, ncolblksNm)
  
  if (obdme) then
     deallocate(obdmeOrb1offset)
     deallocate(obdmeOrb2offset)
     deallocate(obdmeK_ptr)
     deallocate(obdmeOrbKbraket)
     deallocate(obdmeOrbKindx)
     deallocate(obdmeKvals, robdmes)
  endif

  if (numTBtrans .gt. 0) then
     deallocate(transobs, transition_observables)
     call Finalize_H2full
  endif

  !call MPI_Finalize(ierr)
  stop
end program MFDn_Transitions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
