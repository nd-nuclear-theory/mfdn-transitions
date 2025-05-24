!
! Calculates
!
! (1) if (obdme) Reduced One-Body Transition Densities
!
! (2) if (numTBtrans > 0) numTBtrans Two-Body Transition matrix elements
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program MFDn_Transitions
  use MPI
  use nodeinfo
  use SPbasis, only: nclasses, classoffset, nspstates, &
       norbt, norb_p, norb_n, nparticles, j2_orb, classoffsetgen
  use Sparsity, only: generateSparsity, tilediff, tileind, ntiles
  use Wigner3J, only: Wig3J, Init_Wigner3J, Init_Wigner3J_Array_JJK
#ifdef DeltaTz
  use TBME_Tz12, only: Jop, Parop, Tzop
#else
  use TBME_Tz0, only: Jop, Parop, Tzop
#endif
  use OBobservables, only: CalcOBobs, CalcFGTobs
  use Transition_input
  use Wavefunctions_m
  implicit none
  !
  integer(kind=8) :: tprempi, tpostmpi, tcountrate
  real(8) :: tstart, tbefore, tbefore1, tbefore2, tafter,              &
             attest, ttestmin, ttestavg, ttestmax, &
             tketread, tbraread, tobdme, ttbo, ttbo_generateSparsity, ttbo_EvalTBobs,    &
             tketreadmax, tbrareadmax, tobdmemax, ttbomax
  integer :: nketread, nbraread
  integer :: subvec_ket_start, subvec_bra_start, num_submat, submat_ctr
  !
  integer(kind=8) :: dim_ket, dim_bra
  integer :: oprank, statesize, zero=0
  integer :: numgroupids_ket, numstates_ket, numblks_ket
  integer :: numgroupids_bra, numstates_bra, numblks_bra
  integer, dimension(3) :: classoffset_ket, classoffset_bra
  integer :: n_group_offset_ket, n_group_offset_bra
  integer :: subvec_ket, subvec_bra
  integer :: maxnMstates_ket, maxnMstates_bra, maxnMstates
  integer :: TwoJop, Delta_TwoMj, Delta_Par, Delta_Tz
  integer :: i, j, k, ifile, jfile, fh, indxNm_ket, indxNm_bra, Jtot1mx, Jtot2mx
  integer :: brastate, numbras=1, numkets, duplicate, numTBtransTot
  integer, dimension(MaxNumKets) :: ketstate
  integer, dimension(:), allocatable :: Sum_TwoJ, Delta_TwoJ
  logical :: input
  !
  integer :: mnK, mxK, nobdme, nobdme_p, nobdme_n, nobdmeTot
  integer, dimension(:), allocatable :: obdmeOrbKindx
  integer, dimension(:), allocatable :: obdmeOrb_p_offset, obdmeOrb_n_offset
  integer, dimension(:,:), allocatable :: obdmeK_ptr,obdmeOrbKbraket
  real(kind=8), dimension(:), allocatable :: obdmeKvals, robdmes, reducefac
  !
  integer(kind=2), dimension(:,:), allocatable :: groupidlist_ket
  integer(kind=2), dimension(:,:), allocatable :: groupidlist_bra
  integer, dimension(:), allocatable :: Mstateptr_ket, Mstateptr_bra
!   integer, dimension(:), allocatable :: blkgroup_ket, blkoffset_ket
!   integer, dimension(:), allocatable :: blkgroup_bra, blkoffset_bra
  ! integer, dimension(:), allocatable :: numblksNm_ket, numblksNm_bra
  integer, dimension(:), allocatable :: tileptr_ket, tileptr_bra
  !
  real(kind=4), dimension(:), allocatable :: amp_bra
  real(kind=4), dimension(:,:), allocatable :: amp_ket
  real(kind=8) :: brasum, tmp, factEn, factor
  real(kind=8), dimension(:), allocatable :: transobs, transition_observables, TBMEredfac
  real(kind=8), dimension(:), allocatable :: ketsum,            &
       totR20_p, totM1L_p, totM1S_p, totE1Q_p, totE2Q_p, totF0, &
       totR20_n, totM1L_n, totM1S_n, totE1Q_n, totE2Q_n, totGT
  real(kind=8) :: rt4pi = 3.54490775113405
  !
  ! File names
  character(LEN=15) :: resfilename = 'transitions.res'
  character(LEN=15) :: outfilename = 'transitions.out'
  character(LEN=17) :: datfilename = 'transitions.input'
  character(:), allocatable :: mbfile, smwffile
  real, external :: myphase
  !
  namelist /transition_data/ fmass, hbomeg, numTBtrans, obdme, max2K, &
       TwoJ_bra, TwoJ_ket, n_bra, n_ket,     &
       infofilename_bra, infofilename_ket,   &
       basisfilename_bra, basisfilename_ket, &
       smwffilename_bra, smwffilename_ket,   &
       TBMEoperators

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_clock(tprempi)
  !
  call Setup_Parallel
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
  ! Get info about input wave functions
  call InitializeWavefunctions(infofilename_bra, infofilename_ket)
  !
  nparticles = wf_info_bra%num_protons + wf_info_bra%num_neutrons
  if (nparticles /= wf_info_ket%num_protons + wf_info_ket%num_neutrons) then
     print*, 'Different number of particles for bra and ket state'
     print*, wf_info_bra%num_protons, wf_info_bra%num_neutrons
     print*, wf_info_ket%num_protons, wf_info_ket%num_neutrons
     call cancelall(002)
  endif
  Delta_Tz = wf_info_bra%num_protons - wf_info_ket%num_protons
  if (abs(Delta_Tz) .gt. 2) then
     print*, 'Delta Tz larger two'
     print*, wf_info_bra%num_protons,  wf_info_ket%num_protons
     print*, wf_info_bra%num_neutrons, wf_info_ket%num_neutrons
     call cancelall(004)
  endif
  !
  ! determine desired state indices
  brastate = 0
  do i = 1, wf_info_bra%num_states
     if ((wf_info_bra%TwoJ(i) == TwoJ_bra) .and. (wf_info_bra%Jseq(i) == n_bra)) then
        brastate = i
        exit
     end if
  end do
  if (brastate == 0) then
    print*, 'bra state not found with J,n equals', TwoJ_bra, n_bra
    call cancelall(005)
  end if
  numkets = 0
  ketstate(1:MaxNumKets) = 0
  do k = 1, MaxNumKets
     if (n_ket(k) == 0) then
        if (myrank .eq. root) print*, 'WARNING: skipped'
        if (myrank .eq. root) print*, ' ket = ', k, 'because n_ket is zero'
        cycle  ! skip default-initialized n_ket
     end if
     if (TwoJ_ket(k) < 0) then
        if (myrank .eq. root) print*, 'WARNING: skipped'
        if (myrank .eq. root) print*, ' ket = ', k, 'because TwoJ_ket less than zero'
        cycle  ! skip unphysical (and/or default-initialized) TwoJ_ket
     end if
     do i = 1, wf_info_ket%num_states
        if ((wf_info_ket%TwoJ(i) == TwoJ_ket(k)) .and.         &
             (wf_info_ket%Jseq(i) == n_ket(k))) then
           numkets = numkets + 1
           ketstate(numkets) = i
           ! numkets should always be less than or equal to k
           if (numkets < k) then
              TwoJ_ket(numkets) = TwoJ_ket(k)   ! reset TwoJ_ket if necessary
              n_ket(numkets) = n_ket(k)         ! reset n_ket if necessary
              if (myrank .eq. root) print*, 'resetting TwoJ_ket and n_ket because of skipped state(s)'
           endif
           !
           exit
        end if
     end do
     if (myrank .eq. root) print*, ketstate(k), TwoJ_ket(k), n_ket(k), ketstate(k)
     if (ketstate(k) == 0) then
      if (myrank .eq. root) print*, 'ket state not found with J,n equals', TwoJ_ket(k), n_ket(k)
     endif
  end do
  if (numkets == 0) then
     if (myrank .eq. root) print*, 'none of the requested ket states found'
     call cancelall(006)
  end if
  !
  ! Identify and remove duplicate kets at end of ket-array
  duplicate = 0
  duploop: do k = numkets, 1, -1
     do i = k-1, 1, -1
        if (ketstate(i) == ketstate(k)) then
           duplicate = numkets - k + 1
           cycle duploop
        end if
     end do
     if (duplicate > 0) then
        numkets = numkets - duplicate
        if (myrank .eq. root) print*, 'WARNING: removed', duplicate, 'duplicate ket states'
        exit duploop
     end if
  end do duploop
  !
  ! Identify (but not remove) additional duplicate kets
  do k = 1, numkets
     do i = k+1, numkets
        if (ketstate(i) == ketstate(k)) then
           if (myrank .eq. root) print*, 'WARNING: duplicate ket states', k, i
        end if
     end do
  end do
  !
  if (myrank .eq. root) then
     open(unit=7, file=outfilename, status='unknown', action='write')
     write(7,*)
     write(7,*) '  Output from MFDn Transitions in M-scheme basis'
     write(7,*)
     write(7,*) ' Total number of particles ', nparticles
     write(7,*) ' Number of protons  <bra| and |ket>', wf_info_bra%num_protons, wf_info_ket%num_protons
     write(7,*) ' Number of neutrons <bra| and |ket>', wf_info_bra%num_neutrons, wf_info_ket%num_neutrons
     write(7,*)
     write(7,*) ' Parity         <bra| and |ket>', wf_info_bra%parity, wf_info_ket%parity
     write(7,*) ' Projection M_J <bra| and |ket>', 0.5*wf_info_bra%TwoMj, 0.5*wf_info_ket%TwoMj
     write(7,*)
     write(7,*) ' Number of |bra> states ', numbras
     write(7,*) ' Total spin J, n(J), wfn(n) |bra>', 0.5*TwoJ_bra, n_bra, brastate
     write(7,*) ' Number of |ket> states ', numkets
     do k = 1, numkets
        write(7,*) ' Total spin J, n(J), wfn(n) |ket>', 0.5*TwoJ_ket(k), n_ket(k), ketstate(k)
     end do
     write(7,*)
     close(unit=7, status='keep')
   endif
  !
  ! start timing
  tstart = MPI_WTIME()
  !
  ! Read first bra and ket basis files, for consistency checks
  fh = 17
  ifile = 1
  mbfile = TRIM(basisfilename_ket)//achar(ifile/100+48)&
      //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
  open(unit=fh, file=TRIM(mbfile), status='old', action='read', form='unformatted')
  call readMBgroupID_metadata(                      &
      fh, nclasses, nparticles, wf_info_ket%parity, &
      wf_info_ket%TwoMj, wf_info_ket%num_subvec,    &
      indxNm_ket, numgroupids_ket, numstates_ket, numblks_ket)
  close(fh, status='keep')
  !
  mbfile = TRIM(basisfilename_bra)//achar(ifile/100+48)&
       //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
  open(unit=fh, file=TRIM(mbfile), status='old', action='read', form='unformatted')
  call readMBgroupID_metadata(                      &
      fh, nclasses, nparticles, wf_info_bra%parity, &
      wf_info_bra%TwoMj, wf_info_bra%num_subvec,    &
      indxNm_bra, numgroupids_bra, numstates_bra, numblks_bra)
  close(fh, status='keep')
  !
  allocate(Sum_TwoJ(numkets), Delta_TwoJ(numkets))
  do k = 1, numkets
     Sum_TwoJ(k) = TwoJ_ket(k) + TwoJ_bra
     Delta_TwoJ(k) = abs(TwoJ_ket(k) - TwoJ_bra)
  end do
  Delta_TwoMj = wf_info_bra%TwoMj - wf_info_ket%TwoMj
  Delta_Par = wf_info_bra%parity * wf_info_ket%parity
  !
  ! Set up OBDME indices and arrays
  if (max2K .lt. minval(Delta_TwoJ) ) obdme = .false.
  if (abs(Delta_Tz) .gt. 1) obdme = .false.
  !
  if (obdme) then
     !
     mnK = minval(Delta_TwoJ) / 2
     mxK = max2K / 2
     !
     allocate(obdmeK_ptr(norbt, norbt))
     obdmeK_ptr(1:norbt, 1:norbt) = 0
     allocate(obdmeOrb_p_offset(0:mxK+1))
     obdmeOrb_p_offset(0:mxK+1) = 0
     allocate(obdmeOrb_n_offset(0:mxK+1))
     obdmeOrb_n_offset(0:mxK+1) = 0
     !
     if (Delta_Tz .eq. 0) then
        ! <p| (a^\dagger \tilde{a})^K |p> and <n| (a^\dagger \tilde{a})^K |n>
        call CountOBDMEs(mnK, mxK, Delta_Par, nobdme_p, nobdme_n, &
             obdmeOrb_p_offset, obdmeOrb_n_offset)
        nobdme = nobdme_p + nobdme_n
        !
        allocate(obdmeOrbKbraket(2, nobdme))
        allocate(obdmeOrbKindx(nobdme))
        call SetOBDMEindices(mnK, mxK, Delta_Par,    &
             obdmeOrb_p_offset, obdmeOrb_n_offset,   &
             nobdme, obdmeK_ptr, obdmeOrbKbraket, obdmeOrbKindx)
        !
     elseif (Delta_Tz .eq. 1) then
        ! Delta_Tz OBDME < p | ( a^\dagger \tilde{a} )^K | n >
        call CountOBDMEs_Tz(mnK, mxK, Delta_Par,                      &
             zero,   norb_p, norb_p, norb_n, nobdme, obdmeOrb_p_offset)
        nobdme_p = nobdme
        nobdme_n = 0
        !
        allocate(obdmeOrbKbraket(2, nobdme))
        allocate(obdmeOrbKindx(nobdme))
        call SetOBDMEindices_Tz(mnK, mxK, Delta_Par,             &
             zero,   norb_p, norb_p, norb_n, obdmeOrb_p_offset,  &
             nobdme, obdmeK_ptr, obdmeOrbKbraket, obdmeOrbKindx)
        !
     elseif (Delta_Tz .eq. -1) then
        ! Delta_Tz OBDME < n | ( a^\dagger \tilde{a} )^K | p >
        !
        call CountOBDMEs_Tz(mnK, mxK, Delta_Par,                      &
             norb_p, norb_n, zero,   norb_p, nobdme, obdmeOrb_n_offset)
        nobdme_p = 0
        nobdme_n = nobdme
        !
        allocate(obdmeOrbKbraket(2, nobdme))
        allocate(obdmeOrbKindx(nobdme))
        call SetOBDMEindices_Tz(mnK, mxK, Delta_Par,             &
             norb_p, norb_n, zero,   norb_p, obdmeOrb_n_offset,  &
             nobdme, obdmeK_ptr, obdmeOrbKbraket, obdmeOrbKindx)
        !
     endif
     !
     nobdmeTot = nobdme * numbras * numkets
     allocate(obdmeKvals(nobdmeTot), robdmes(nobdmeTot))
     obdmeKvals(1:nobdmeTot) = 0.0d0
     !
  endif
  !
  ! Read Two-Body matrix elements
  if (numTBtrans .gt. 0) then
     !
     ! Read Two-Body Matrix Elements (max. number set in preprocessor)
     call read_TBME_multi(numTBtrans, TBMEoperators)
     !
     ! Consistency checks
     TwoJop = 2 * Jop
     if (myrank .eq. root) then
       if (TwoJop .lt. minval(Delta_TwoJ)) then
          print*, 'TwoJop too small for min(Delta_TwoJ)', TwoJop, minval(Delta_TwoJ)
          call cancelall(200)
       endif
       if (TwoJop .gt. maxval(Sum_TwoJ)) then
          print*, 'TwoJop too large for max(Sum_TwoJ)', TwoJop, maxval(Sum_TwoJ)
          call cancelall(202)
       endif
       if (TwoJop .lt. abs(Delta_TwoMj)) then
          print*, 'TwoJop too small for Delta_TwoMj', TwoJop, Delta_TwoMj
          call cancelall(205)
       endif
       if (Parop .ne. Delta_Par) then
          print*, 'inconsistent Parop', Parop, Delta_Par
          call cancelall(210)
       endif
       if (Tzop .ne. abs(Delta_Tz)) then
          print*, 'inconsistent Tzop', Tzop, Delta_Tz
          call cancelall(220)
       endif
     endif
     !
     oprank = 2
     numTBtransTot = numTBtrans * numbras * numkets
     allocate(transobs(numTBtransTot))
     allocate(transition_observables(numTBtransTot))
     transition_observables(1:numTBtransTot) = 0.0d0
     !
  else
     oprank = 1
  endif
  statesize = nparticles + oprank
  !
  deallocate(Sum_TwoJ, Delta_TwoJ)

  ! calculate possible groupid offsets
  classoffset_bra = classoffsetgen(wf_info_bra%orbital_list%norb_p,            &
                                   wf_info_bra%orbital_list%norb_n,            &
                                   wf_info_bra%orbital_list%j2_orb)
     !
  classoffset_ket = classoffsetgen(wf_info_ket%orbital_list%norb_p,            &
                                   wf_info_ket%orbital_list%norb_n,            &
                                   wf_info_ket%orbital_list%j2_orb)
  n_group_offset_bra = 0
  n_group_offset_ket = 0
  if (classoffset_bra(2) > classoffset_ket(2)) then
     n_group_offset_ket = classoffset_bra(2) - classoffset_ket(2)
  end if
  if (classoffset_bra(2) < classoffset_ket(2)) then
     n_group_offset_bra = classoffset_ket(2) - classoffset_bra(2)
  end if
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank .eq. root) then
     open(unit=7, file=outfilename, status='old', action='write', position='append')
     if (obdme) then
        write(7,*) ' Calculation Reduced One-Body Density Matrix Elements'
        write(7,*) '   Lower limit on K', mnK
        write(7,*) '   Upper limit on K', mxK
        if (Delta_Par .ne. 1) write(7,*) '   Parity changing ', Delta_Par
        if (Delta_Tz  .eq. 0) then
           write(7,*)
           write(7,*) '   Number of proton OBDMEs ', nobdme_p
           write(7,*) '   Number of neutron OBDMEs', nobdme_n
        elseif (abs(Delta_Tz) .eq. 1) then
           write(7,*) '   Isospin changing', Delta_Tz
           write(7,*)
           write(7,*) '   Number of <p|O|n> OBDMEs', nobdme
        endif
        write(7,*)
     endif
     if (numTBtrans .gt. 0) then
        write(7,*) ' Calculation Matrix Elements for Two-Body Operators'
        write(7,*)
        write(7,*) ' Spin of Operator', Jop
        write(7,*) '  M_J of Operator', Delta_TwoMj/2
        write(7,*) '  Parity of Oprtr', Parop
        write(7,*)
     endif
     !
     close(unit=7, status='keep')
     !
     open(unit=10, file=resfilename, status='unknown', action='write')
1100 format('', A)
1101 format('', A,' = ', A)
1108 format('', A,' = ',i8)
1118 format('', A,' = ',i18)
1128 format('', A,' = ',2i8)
1110 format('', A,' = ',f10.1)
     write(10, *)
     write(10, 1100) '[MFDn Transitions]'
     write(10, 1101) 'Version', '0'
     write(10, 1101) 'Revision', VCS_REVISION
     write(10, 1101) 'Platform'
     write(10, 1101) 'Username'
     write(10, 1108) 'MPIranks ', nprocs
     write(10, 1108) 'OMPthreads', numthreads
     write(10, *)
     write(10, 1100) '[PARAMETERS]'
     write(10, *)
     write(10, 1100) '[Bra basis]'
     write(10, 1101) 'basisfilename', TRIM(basisfilename_bra)
     write(10, 1108) 'num_protons  ', wf_info_bra%num_protons
     write(10, 1108) 'num_neutrons ', wf_info_bra%num_neutrons
     write(10, 1110) 'Mj           ', 0.5*wf_info_bra%TwoMj
     write(10, 1108) 'parity       ', wf_info_bra%parity
     write(10, *)
     write(10, 1100) '[Ket basis]'
     write(10, 1101) 'basisfilename', TRIM(basisfilename_ket)
     write(10, 1108) 'num_protons  ', wf_info_ket%num_protons
     write(10, 1108) 'num_neutrons ', wf_info_ket%num_neutrons
     write(10, 1110) 'Mj           ', 0.5*wf_info_ket%TwoMj
     write(10, 1108) 'parity       ', wf_info_ket%parity
     write(10, *)
     write(10, 1100) '[State information]'
     write(10, 1101) 'smwffilename_bra', TRIM(smwffilename_bra)
     write(10, 1108) 'num_bras', numbras
     write(10, 1108) 'brastate', brastate
     write(10, 1110) 'J_bra   ', 0.5*TwoJ_bra
     write(10, 1108) 'n_bra   ', n_bra
     write(10, 1101) 'smwffilename_ket', TRIM(smwffilename_ket)
     write(10, 1108) 'num_kets', numkets
     do k = 1, numkets
        write(10, 1108) 'ketstate', ketstate(k)
        write(10, 1110) 'J_ket   ', 0.5*TwoJ_ket(k)
        write(10, 1108) 'n_ket   ', n_ket(k)
     end do
     write(10, *)
     close(unit=10, status='keep')
  endif

  Jtot1mx = (maxval(j2_orb)+1) / 2
  Jtot2mx = maxval(j2_orb)
  call Init_Wigner3j(Jtot1mx, Jtot2mx)
  call Init_Wigner3J_Array_JJK(Jop, abs(Delta_TwoMj/2))

  if (myrank .eq. root) then
      print*, 'initial setup time', MPI_WTIME()-tstart
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generate and distribute work assignments
  ! call SplitProcs(wf_info_ket%num_subvec, wf_info_bra%num_subvec, &
  !                 subvec_ket_start, subvec_bra_start, num_submat)
  !
  call SplitProcsSeq(wf_info_ket%num_subvec, wf_info_bra%num_subvec, &
       subvec_ket_start, subvec_bra_start, num_submat)
  !
  if (myrank .eq. root) then
     open(unit=7, file=outfilename, status='old', action='write', position='append')
     !
     write(7, *) ' Number of MPI ranks           ', nprocs
     write(7, *) ' Number of submatrices per rank', num_submat
     if ( nprocs * num_submat .eq. wf_info_ket%num_subvec * wf_info_bra%num_subvec) then
        write(7, *) ' Balanced number of MPI ranks'
     else
        write(7, *) ' Not quite balanced number of MPI ranks'
     endif
     write(7, *)
     write(7, *) ' Number of OMP threads per rank', numthreads
     write(7, *)
     !
     close(unit=7, status='keep')
  endif
  !
  ! Loop over assigned row (< bra |, psibar) and column (| ket >, psi) ranks
  dim_ket = 0
  allocate(ketsum(numkets))
  ketsum(1:numkets) = 0.d0
  submat_ctr = 0
  tketread = 0.d0
  tbraread = 0.d0
  tobdme = 0.d0
  ttbo = 0.d0
  ttbo_generateSparsity = 0.d0
  ttbo_EvalTBobs = 0.d0
  tketreadmax = 0.d0
  tbrareadmax = 0.d0
  tobdmemax = 0.d0
  ttbomax = 0.d0
  nketread = 0
  nbraread = 0
  ketloop: do subvec_ket = subvec_ket_start, wf_info_ket%num_subvec-1
     ! start timing for col read
     tbefore = MPI_WTIME()
     !
     ! Read Many-body basis for column ( | ket >, psi)
     ifile = subvec_ket+1
     mbfile = TRIM(basisfilename_ket)//achar(ifile/100+48)&
          //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
     open(unit=fh, file=TRIM(mbfile), status='old', action='read', form='unformatted')
     !
     call readMBgroupID_metadata(fh, nclasses, nparticles,              &
        wf_info_ket%parity, wf_info_ket%TwoMj, wf_info_ket%num_subvec, &
        indxNm_ket, numgroupids_ket, numstates_ket, numblks_ket)
     !
     read(fh) ! numblksNm_ket(1:indxNm_ket)    ! not used
     !
     allocate(Mstateptr_ket(numgroupids_ket+1))
     read(fh) Mstateptr_ket(1:numgroupids_ket+1)
     !
     allocate(groupidlist_ket(nparticles, numgroupids_ket))
     read(fh) groupidlist_ket(1:nparticles, 1:numgroupids_ket)
     !
     ! allocate(blkgroup_ket(numblks_ket+1))
     read(fh) ! blkgroup_ket(1:numblks_ket+1)  ! not used
     !
     ! allocate(blkoffset_ket(numblks_ket+1))
     read(fh) ! blkoffset_ket(1:numblks_ket+1)  ! not used
     !
     close(fh, status='keep')
     !
     ! apply offset to neutrons in ket groupids
     if (n_group_offset_ket > 0) then
      do i = 1, numgroupids_ket
         do j = wf_info_ket%num_protons+1, nparticles
            groupidlist_ket(j, i) = groupidlist_ket(j, i) + n_group_offset_ket
         end do
      end do
   end if
   !
   ! Set max number of M-scheme states in largest group
     maxnMstates_ket = 0
     do i = 1, numgroupids_ket
        maxnMstates_ket = max(maxnMstates_ket, Mstateptr_ket(i+1)-Mstateptr_ket(i))
     enddo
     !
     ! Read wavefunction | ket >
     allocate(amp_ket(numkets, numstates_ket))
     smwffile = TRIM(smwffilename_ket)//achar(ifile/100+48)&
          //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
     open(unit=fh, file=TRIM(smwffile), status='unknown', action='read', form='unformatted')
     k = 0
     ketread: do i = 1, wf_info_ket%num_states
        do j = 1, numkets
           if (ketstate(j) == i) then
              read(fh) amp_ket(j, 1:numstates_ket)
              k = k + 1
              if (k < numkets) then
                 cycle ketread
              else
                 exit ketread
              endif
           endif
        end do
        read(fh)
     enddo ketread
     close(fh, status='keep')
     ! Test on norm
     !$omp parallel do default(shared)          &
     !$omp             private(j, k)            &
     !$omp             reduction(+: ketsum)
     do j = 1, numstates_ket
        do k = 1, numkets
           ketsum(k) = ketsum(k) + amp_ket(k,j)**2
        enddo
     enddo
     !$omp end parallel do
     dim_ket = dim_ket + numstates_ket
     !
     allocate(tileptr_ket(numgroupids_ket+1))
     tafter = MPI_WTIME()
     tketread = tketread + tafter-tbefore
     tketreadmax = max(tketreadmax, tafter-tbefore)
     nketread = nketread + 1
     !
     dim_bra = 0
     brasum = 0.0d0
     braloop: do subvec_bra = subvec_bra_start, wf_info_bra%num_subvec-1
        ! start timing for row read
        tbefore = MPI_WTIME()
        !
        !!print*, "rank", myrank, "starting block", subvec_ket, subvec_bra
        !
        ! Read Many-body basis for row ( < bra |, psibar)
        jfile = subvec_bra+1
        mbfile = TRIM(basisfilename_bra)//achar(jfile/100+48)&
             //achar(mod(jfile/10,10)+48)//achar(mod(jfile,10)+48)
        open(unit=fh, file=TRIM(mbfile), status='old', action='read', form='unformatted')
        call readMBgroupID_metadata(fh, nclasses, nparticles,               &
            wf_info_bra%parity, wf_info_bra%TwoMj, wf_info_bra%num_subvec, &
            indxNm_bra, numgroupids_bra, numstates_bra, numblks_bra)
        !
        read(fh) ! numblksNm_bra(1:indxNm_bra)  ! not used
        !
        allocate(Mstateptr_bra(numgroupids_bra+1))
        read(fh) Mstateptr_bra(1:numgroupids_bra+1)
        !
        allocate(groupidlist_bra(nparticles, numgroupids_bra))
        read(fh) groupidlist_bra(1:nparticles, 1:numgroupids_bra)
        !
        ! allocate(blkgroup_bra(numblks_bra+1))
        read(fh) ! blkgroup_bra(1:numblks_bra+1)  ! not used
        !
        ! allocate(blkoffset_bra(numblks_bra+1))
        read(fh) ! blkoffset_bra(1:numblks_bra+1)  ! not used
        !
        close(fh, status='keep')
        !
        ! apply offset to neutrons in bra groupids
        if (n_group_offset_bra > 0) then
           do i = 1, numgroupids_bra
              do j = wf_info_bra%num_protons+1, nparticles
                 groupidlist_bra(j, i) = groupidlist_bra(j, i) + n_group_offset_bra
              end do
           end do
        end if
        !
        ! Set max number of M-scheme states in largest group
        maxnMstates_bra = 0
        do i = 1, numgroupids_bra
           maxnMstates_bra = max(maxnMstates_bra, Mstateptr_bra(i+1)-Mstateptr_bra(i))
        enddo
        maxnMstates = max(maxnMstates_ket, maxnMstates_bra)
        !
        ! Read wavefunction < bra |
        allocate(amp_bra(numstates_bra))
        smwffile = TRIM(smwffilename_bra)//achar(jfile/100+48)&
             //achar(mod(jfile/10,10)+48)//achar(mod(jfile,10)+48)
        open(unit=fh, file=TRIM(smwffile), status='unknown', action='read', form='unformatted')
        do i = 1, brastate-1
           read(fh)
        enddo
        read(fh) amp_bra(1:numstates_bra)
        close(fh, status='keep')
        ! Test on norm
        do j = 1, numstates_bra
           brasum = brasum + amp_bra(j)**2
        enddo
        dim_bra = dim_bra + numstates_bra
        !
        allocate(tileptr_bra(numgroupids_bra+1))
        !
        tafter = MPI_WTIME()
        tbraread = tbraread + tafter-tbefore
        tbrareadmax = max(tbrareadmax, tafter-tbefore)
        nbraread = nbraread + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! NOTE: difference with MFDn -- TwoMj could be different for < bra | and | ket >

        if (obdme) then
           tbefore = MPI_WTIME()
           call generateSparsity(nparticles, 1, nspstates, classoffset(2), &
                wf_info_ket%num_protons, wf_info_ket%num_neutrons, numgroupids_ket, groupidlist_ket, &
                wf_info_bra%num_protons, wf_info_bra%num_neutrons, numgroupids_bra, groupidlist_bra, &
                tileptr_ket)                  ! also sets tileind and tilediff
           !
           call EvalOBDME(maxnMstates,                                     &
                wf_info_ket%num_protons, wf_info_ket%num_neutrons, numgroupids_ket, groupidlist_ket, &
                wf_info_ket%TwoMj, Mstateptr_ket, numstates_ket, numkets, amp_ket,                   &
                wf_info_bra%num_protons, wf_info_bra%num_neutrons, numgroupids_bra, groupidlist_bra, &
                wf_info_bra%TwoMj, Mstateptr_bra, numstates_bra, numbras, amp_bra,                   &
                tileptr_ket, ntiles, tileind, tilediff,                     &
                mnK, mxK, nobdme, obdmeK_ptr, robdmes)
           obdmeKvals(1:nobdmeTot) = obdmeKvals(1:nobdmeTot) + robdmes(1:nobdmeTot)
           !
           deallocate(tilediff, tileind)     ! allocated in generateSparsity
           tobdme = tobdme + MPI_WTIME()-tbefore
        endif
        !
        if (numTBtrans .gt. 0) then
           tbefore = MPI_WTIME()
           if (Delta_Tz .lt. 0) then
              call generateSparsity(nparticles, oprank, nspstates, classoffset(2), &
                   wf_info_bra%num_protons, wf_info_bra%num_neutrons, numgroupids_bra, groupidlist_bra, &
                   wf_info_ket%num_protons, wf_info_ket%num_neutrons, numgroupids_ket, groupidlist_ket, &
                   tileptr_bra)                  ! also sets tileind and tilediff
              !
              ttbo_generateSparsity = ttbo_generateSparsity + MPI_WTIME()-tbefore
              tbefore1 = MPI_WTIME()
              call EvalTBobs(maxnMstates,                                     &
                   wf_info_bra%num_protons, wf_info_bra%num_neutrons, numgroupids_bra, groupidlist_bra, &
                   wf_info_bra%TwoMj, Mstateptr_bra, numstates_bra, numbras, amp_bra,                   &
                   wf_info_ket%num_protons, wf_info_ket%num_neutrons, numgroupids_ket, groupidlist_ket, &
                   wf_info_ket%TwoMj, Mstateptr_ket, numstates_ket, numkets, amp_ket,                   &
                   tileptr_bra, ntiles, tileind, tilediff, numTBtrans, transobs)
              ttbo_EvalTBobs = ttbo_EvalTBobs + MPI_WTIME()-tbefore1
              ! apply conjugation phase from exchanging bra and ket
              transobs(1:numTBtransTot) = myphase(Delta_TwoMj/2 + Delta_Tz) * transobs(1:numTBtransTot)
           else
              call generateSparsity(nparticles, oprank, nspstates, classoffset(2), &
                   wf_info_ket%num_protons, wf_info_ket%num_neutrons, numgroupids_ket, groupidlist_ket, &
                   wf_info_bra%num_protons, wf_info_bra%num_neutrons, numgroupids_bra, groupidlist_bra, &
                   tileptr_ket)                  ! also sets tileind and tilediff
              !
              ttbo_generateSparsity = ttbo_generateSparsity + MPI_WTIME()-tbefore
              tbefore1 = MPI_WTIME()
              call EvalTBobs(maxnMstates,                                     &
                   wf_info_ket%num_protons, wf_info_ket%num_neutrons, numgroupids_ket, groupidlist_ket, &
                   wf_info_ket%TwoMj, Mstateptr_ket, numstates_ket, numkets, amp_ket,                   &
                   wf_info_bra%num_protons, wf_info_bra%num_neutrons, numgroupids_bra, groupidlist_bra, &
                   wf_info_bra%TwoMj, Mstateptr_bra, numstates_bra, numbras, amp_bra,                   &
                   tileptr_ket, ntiles, tileind, tilediff, numTBtrans, transobs)
              ttbo_EvalTBobs = ttbo_EvalTBobs + MPI_WTIME()-tbefore1
           endif
           transition_observables(1:numTBtransTot) =                          &
                transition_observables(1:numTBtransTot) + transobs(1:numTBtransTot)
           !
           deallocate(tilediff, tileind)       ! allocated in generateSparsity
           ttbo = ttbo + MPI_WTIME()-tbefore
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate(Mstateptr_bra)
        deallocate(groupidlist_bra)
        ! deallocate(blkgroup_bra)
        ! deallocate(blkoffset_bra)
        deallocate(amp_bra)
        deallocate(tileptr_bra)

        submat_ctr = submat_ctr + 1
        !!print*, "rank", myrank, "finished block", subvec_ket, subvec_bra
        !!print*, myrank, subvec_ket, subvec_bra, transition_observables(1)
        if (submat_ctr .ge. num_submat) exit braloop

     enddo braloop
     subvec_bra_start = 0

     deallocate(Mstateptr_ket)
     deallocate(groupidlist_ket)
     ! deallocate(blkgroup_ket)
     ! deallocate(blkoffset_ket)
     deallocate(amp_ket)
     deallocate(tileptr_ket)

     if (submat_ctr .ge. num_submat) exit ketloop

  enddo ketloop

  if (myrank .eq. root) then
     tbefore = MPI_WTIME()
     call MPI_Reduce(MPI_IN_PLACE, tketread, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
     call MPI_Reduce(MPI_IN_PLACE, tbraread, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)

     call MPI_Reduce(MPI_IN_PLACE, tketreadmax, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, root, icomm, ierr)
     call MPI_Reduce(MPI_IN_PLACE, tbrareadmax, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, root, icomm, ierr)

     call MPI_Reduce(MPI_IN_PLACE, nketread, 1, &
          MPI_INTEGER, MPI_SUM, root, icomm, ierr)
     call MPI_Reduce(MPI_IN_PLACE, nbraread, 1, &
          MPI_INTEGER, MPI_SUM, root, icomm, ierr)

     call MPI_Reduce(MPI_IN_PLACE, brasum, 1,       &
          MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
     call MPI_Reduce(MPI_IN_PLACE, ketsum, numkets, &
          MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)

     if (obdme) then
        call MPI_Reduce(MPI_IN_PLACE, obdmeKvals, nobdmeTot,      &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
        call MPI_Reduce(MPI_IN_PLACE, tobdme, 1, &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
     endif
     if (numTBtrans .gt. 0) then
        call MPI_Reduce(MPI_IN_PLACE, transition_observables, numTBtransTot, &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
        call MPI_Reduce(MPI_IN_PLACE, ttbo, 1, &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
        call MPI_Reduce(MPI_IN_PLACE, ttbo_generateSparsity, 1, &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
        call MPI_Reduce(MPI_IN_PLACE, ttbo_EvalTBobs, 1, &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
     endif
     print*, 'time collecting', MPI_WTIME()-tbefore
  else
     call MPI_Reduce(tketread, 0, 1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
     call MPI_Reduce(tbraread, 0, 1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)

     call MPI_Reduce(tketreadmax, 0, 1,        &
          MPI_DOUBLE_PRECISION, MPI_MAX, root, icomm, ierr)
     call MPI_Reduce(tbrareadmax, 0, 1,        &
          MPI_DOUBLE_PRECISION, MPI_MAX, root, icomm, ierr)

     call MPI_Reduce(nketread, 0, 1,        &
          MPI_INTEGER, MPI_SUM, root, icomm, ierr)
     call MPI_Reduce(nbraread, 0, 1,        &
          MPI_INTEGER, MPI_SUM, root, icomm, ierr)

     call MPI_Reduce(brasum, 0, 1,          &
          MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
     call MPI_Reduce(ketsum, 0, numkets,    &
          MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)

     if (obdme) then
        call MPI_Reduce(obdmeKvals, 0, nobdmeTot,                 &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
        call MPI_Reduce(tobdme, 0, 1,      &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
     endif

     if (numTBtrans .gt. 0) then
        call MPI_Reduce(transition_observables, 0, numTBtransTot, &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
        call MPI_Reduce(ttbo, 0, 1,      &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
        call MPI_Reduce(ttbo_generateSparsity, 0, 1,      &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
        call MPI_Reduce(ttbo_EvalTBobs, 0, 1,      &
             MPI_DOUBLE_PRECISION, MPI_SUM, root, icomm, ierr)
     endif
  endif

  if (myrank .eq. root) then
     print*, '    num ket reads', nketread
     print*, 'avg ket read time', tketread/nketread
     print*, 'max ket read time', tketreadmax
     print*, '    num bra reads', nbraread
     print*, 'avg bra read time', tbraread/nbraread
     print*, 'max bra read time', tbrareadmax
     print*, '   avg obdme time', tobdme/(wf_info_ket%num_subvec*wf_info_bra%num_subvec)
     print*, '     avg tbo time', ttbo/(wf_info_ket%num_subvec*wf_info_bra%num_subvec)
     print*, '       avg tbo_generateSparsity time', ttbo_generateSparsity/(wf_info_ket%num_subvec*wf_info_bra%num_subvec)
     print*, '       avg tbo_EvalTBobs time', ttbo_EvalTBobs/(wf_info_ket%num_subvec*wf_info_bra%num_subvec)
     print*, '   num smwf pairs', wf_info_ket%num_subvec*wf_info_bra%num_subvec
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank .eq. root) then
     tbefore = MPI_WTIME()
     open(unit=7, file=outfilename, status='old', action='write', position='append')
     !
     if (abs(brasum-1.d0) > 1d-6) then
        write(7, *) 'WARNING: bra norm incorrect', brasum
        print*, 'WARNING: bra norm incorrect', brasum
     endif
     do k = 1, numkets
        if (abs(ketsum(k)-1.d0) > 1d-6) then
           write(7, *) 'WARNING: ket norm incorrect', k, ketsum(k)
           print*, 'WARNING: ket norm incorrect', k, ketsum(k)
        endif
     end do
     !
102  format('',f8.1, i8, 2(2x,f10.4) )
104  format('',f8.1, i8, 4(2x,f10.4) )
201  format(f8.1, i8, e16.6, MaxNumKets(x,e15.5) )
     !
     if (obdme) then
        allocate(reducefac(0:mxK))
        if (Delta_Tz .eq. 0) then
           allocate(totR20_p(numkets), totM1L_p(numkets), totM1S_p(numkets), totE2Q_p(numkets))
           allocate(totR20_n(numkets), totM1L_n(numkets), totM1S_n(numkets), totE2Q_n(numkets))
        else
           allocate(totF0(numkets), totGT(numkets))
        endif
        !
        do k = 1, numkets
           j = (k-1)*nobdme + 1
           if (Delta_Tz .eq. 0) then
              call ReduceOBDME(TwoJ_bra, TwoJ_ket(k), wf_info_bra%TwoMj, wf_info_ket%TwoMj,  &
                   mnK, mxK, Delta_Par, nobdme, obdmeKvals(j), obdmeK_ptr,   &
                   obdmeOrbKindx, reducefac, robdmes)
              !
              ! Calculate static OB observables
              call CalcOBobs(Delta_Par, mxK, nobdme, robdmes, obdmeOrb_p_offset,&
                   obdmeOrbKbraket, totR20_p(k), totM1L_p(k), totM1S_p(k), totE2Q_p(k))
              !
              call CalcOBobs(Delta_Par, mxK, nobdme, robdmes, obdmeOrb_n_offset,&
                   obdmeOrbKbraket, totR20_n(k), totM1L_n(k), totM1S_n(k), totE2Q_n(k))
              !
           elseif (Delta_Tz .eq. 1) then
              call ReduceOBDME_Tz(TwoJ_bra, TwoJ_ket(k), wf_info_bra%TwoMj, wf_info_ket%TwoMj, &
                   mnK, mxK, Delta_Par, zero,   norb_p, norb_p, norb_n, nobdme, &
                   obdmeKvals(j), obdmeK_ptr, obdmeOrbKindx, reducefac, robdmes)
              !
              call CalcFGTobs(Delta_Par, mxK, nobdme, robdmes, obdmeOrb_p_offset,&
                   obdmeOrbKbraket, totF0(k), totGT(k))
              !
           elseif (Delta_Tz .eq. -1) then
              call ReduceOBDME_Tz(TwoJ_bra, TwoJ_ket(k), wf_info_bra%TwoMj, wf_info_ket%TwoMj, &
                   mnK, mxK, Delta_Par, norb_p, norb_n, zero,   norb_p, nobdme, &
                   obdmeKvals(j), obdmeK_ptr, obdmeOrbKindx, reducefac, robdmes)
              !
              call CalcFGTobs(Delta_Par, mxK, nobdme, robdmes, obdmeOrb_n_offset,&
                   obdmeOrbKbraket, totF0(k), totGT(k))
              ! apply conjugation phase from exchanging bra and ket
              totF0(k) = myphase((TwoJ_bra - TwoJ_ket(k))/2 + Delta_Tz) * totF0(k)
              totGT(k) = myphase((TwoJ_bra - TwoJ_ket(k))/2 + Delta_Tz) * totGT(k)
              !
           endif
           !
           call WriteOBDME(brastate, ketstate(k), mnK, mxK, nobdme, nobdme_p, nobdme_n, &
                obdmeOrb_p_offset, obdmeOrb_n_offset, obdmeOrbKbraket,      &
                reducefac, robdmes, k)
        end do
        !
        write(7,*) 'Reduced OBDMs written to file for'
        do k = mnK, mxK
           if (reducefac(k) .eq. 0.d0) then
              write(7,*) 'WARNING: skipped'
              write(7,*) '  K = ', k, 'because Wig3j is zero'
              print*, 'WARNING: skipped'
              print*, '  K = ', k, 'because Wig3j is zero'
           else
              write(7,*) '  K = ', k
           endif
        enddo
        write(7,*)
        deallocate(reducefac)
        !
        if (hbomeg .gt. 0) then
           write(7,*) 'OBDM observables'
           write(7, *) ' J <bra|   ', 0.5*TwoJ_bra
           write(7, *) ' n <bra|   ', n_bra
        write(7, *)
        write(7, *) ' |ket> J   n(J)    <bra || OBop || ket> '
        if (Delta_Tz .eq. 0) then
           factEn = (197.327d0)**2/(fmass*hbomeg)        ! scale factor for E2
           if (Delta_par .eq. -1) factEn = sqrt(factEn)  ! scale factor for E1 and M2
           !
           ! see Suhonen 6.23 and 6.24 for the factor of sqrt(4\pi)
           factor = 1.d0/rt4pi
           !
           ! for M1 moments (see Suhonen 6.52)
           ! factor = sqrt(TwoJ_bra / (3.d0*(TwoJ_bra+2.0d0)*(TwoJ_bra+1.0d0)))
           totM1L_p(1:numkets) = totM1L_p(1:numkets) * factor
           totM1S_p(1:numkets) = totM1S_p(1:numkets) * factor
           totM1L_n(1:numkets) = totM1L_n(1:numkets) * factor
           totM1S_n(1:numkets) = totM1S_n(1:numkets) * factor
           if (Delta_par .eq. -1) then
              write(7,*) 'Reduced M2 transition matrix elements (orbital p, n; spin p, n)'
              do k = 1, numkets
                 write(7, 104) 0.5*TwoJ_ket(k), n_ket(k),   &
                      factEn * totM1L_p(k), factEn * totM1L_n(k), factEn * totM1S_p(k), factEn * totM1S_n(k)
              end do
           else
              write(7,*) 'Reduced M1 transition matrix elements (orbital p, n; spin p, n)'
              do k = 1, numkets
                 write(7, 104) 0.5*TwoJ_ket(k), n_ket(k),   &
                      totM1L_p(k), totM1L_n(k), totM1S_p(k), totM1S_n(k)
              end do
           endif
           write(7,*)
           !
           ! for E2 moments (see Suhonen 6.53)
           ! factor = 2.d0 * sqrt(TwoJ_bra*(TwoJ_bra-1.0d0) /                   &
           !      (5.d0*(TwoJ_bra+2.0d0) * (TwoJ_bra+1.0d0)*(TwoJ_bra+3.0d0)) )
           totE2Q_p(1:numkets) = totE2Q_p(1:numkets) * factor * factEn
           totE2Q_n(1:numkets) = totE2Q_n(1:numkets) * factor * factEn
           if (Delta_par .eq. -1) then
              write(7,*) 'Reduced E1 transition matrix elements (proton, neutron)'
           else
              write(7,*) 'Reduced E2 transition matrix elements (proton, neutron)'
           endif
           do k = 1, numkets
              write(7, 102) 0.5*TwoJ_ket(k), n_ket(k), totE2Q_p(k), totE2Q_n(k)
           end do
           write(7,*)
           !
        elseif (abs(Delta_Tz) .eq. 1) then
           !
           ! totF0(1:numkets) = totF0(1:numkets) * factor
           ! totGT(1:numkets) = totGT(1:numkets) * factor
           write(7,*) 'Reduced Fermi and GT transition matrix elements (?)'
           do k = 1, numkets
              write(7, 102) 0.5*TwoJ_ket(k), n_ket(k), totF0(k), totGT(k)
           end do
              write(7,*)
              !
           endif
        endif
        !
        if (Delta_Tz .eq. 0) then
           deallocate(totR20_p, totM1L_p, totM1S_p, totE2Q_p)
           deallocate(totR20_n, totM1L_n, totM1S_n, totE2Q_n)
        else
           deallocate(totF0, totGT)
        endif
        !
     endif
     !
     if (numTBtrans .gt. 0) then
        !
        write(7, *) ' J operator', Jop
        do i = 1, numTBtrans
           write(7, 1101) ' TBMEfilename', TRIM(TBMEoperators(i))
        end do
        write(7, *)
        !
        write(7, *) ' J <bra|   ', 0.5*TwoJ_bra
        write(7, *) ' n(J) <bra|', n_bra
        write(7, *)
        write(7, *) ' |ket> J    n(J)  W3J(Jbra,Jop,Jket)  <bra | TBMEfile | ket> '
        allocate(TBMEredfac(numkets))
        do k = 1, numkets
           j = (k-1)*numTBtrans
           !
           tmp = Wig3J(TwoJ_bra, TwoJop, TwoJ_ket(k), -wf_info_bra%TwoMj, Delta_TwoMj, wf_info_ket%TwoMj)
           write(7, 201) 0.5*TwoJ_ket(k), n_ket(k), tmp, transition_observables(j+1:j+numTBtrans)
           !
           if (abs(tmp).gt.1d-12) then
              TBMEredfac(k) = 1.d0 / (tmp * myphase((TwoJ_bra-wf_info_bra%TwoMj)/2))
           else
              TBMEredfac(k) = 0.d0
              write(7,*) 'WARNING: Wigner3J(J_bra, J_op, J_ket) is zero for 2J_ket', TwoJ_ket(k), n_ket(k)
              print*, 'WARNING: Wigner3J is zero', tmp
              print*, 'WARNING: for (2J_bra, 2J_op, 2J_ket)', TwoJ_bra, TwoJop, TwoJ_ket(k)
              print*, 'WARNING: and ( -2Mj, 2Mj_op,   2Mj )', -wf_info_bra%TwoMj, Delta_TwoMj, wf_info_ket%TwoMj
              print*, 'WARNING: Matrix elements of Two-Body operator should be zero', transition_observables(j+1:j+numTBtrans)
           endif
           !
           ! do i = 1, numTBtrans
           !   j = (k-1)*numTBtrans + i
           !   write(7, 22) TRIM(TBMEoperators(i)), transition_observables(j), &
           !   TBMEredfac * transition_observables(j)
           !   write(10, 22) TRIM(TBMEoperators(i)), transition_observables(j), &
           !   TBMEredfac * transition_observables(j)
           ! enddo
           !
        end do
        write(7,*)
        close(unit=7, status='keep')
        !
        open(unit=10, file=resfilename, status='old', action='write', position='append')
        write(10, 1100) '[RESULTS]'
        write(10, *)
2200    format('# ',A3,' ',A3,' ',A3,'  ',A)
2201    format('  ',i3,' ',i3,' ',i3,'  ',A)
2210    format('# ',A4,  ' ',A3,' ',A3,'  ',A4,  ' ',A3,' ',A3,'  ',A15)
2211    format('  ',f4.1,' ',i3,' ',i3,'  ',f4.1,' ',i3,' ',i3,'  ',e15.8)
        do i = 1, numTBtrans
           write(10, 1100) '[Two-body observable]'
           write(10, 2200) 'J0', 'g0', 'Tz0', 'name'
           write(10, 2201) Jop, (1-parop)/2, Tzop, TRIM(TBMEoperators(i))
           write(10, *)
           write(10, 2210) 'Jf', 'gf', 'nf', 'Ji', 'gi', 'ni', 'rme'
           do k = 1, numkets
              j = (k-1)*numTBtrans+i
              if (TBMEredfac(k) /= 0.0) then
                 write(10, 2211) 0.5*TwoJ_bra, (1-wf_info_bra%parity)/2, n_bra,       &
                                 0.5*TwoJ_ket(k), (1-wf_info_ket%parity)/2, n_ket(k), &
                                 TBMEredfac(k) * transition_observables(j)
              endif
           end do
           write(10,*)
         end do
        close(unit=10, status='keep')
        deallocate(TBMEredfac)
     endif
     !
     open(unit=7, file=outfilename, status='old', action='write', position='append')
     print*, 'time writing output', MPI_WTIME()-tbefore
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! deallocate(numblksNm_bra, numblksNm_ket)

  if (obdme) then
     deallocate(obdmeK_ptr)
     deallocate(obdmeOrbKbraket)
     deallocate(obdmeOrbKindx)
     deallocate(obdmeKvals, robdmes)
     if (allocated(obdmeOrb_p_offset)) deallocate(obdmeOrb_p_offset)
     if (allocated(obdmeOrb_n_offset)) deallocate(obdmeOrb_n_offset)
  endif

  if (numTBtrans .gt. 0) then
     deallocate(transobs, transition_observables)
  endif

  if (myrank .eq. root) print*, 'total time', MPI_WTIME()-tstart
  call MPI_Finalize(ierr)
  !
  call system_clock(tpostmpi, tcountrate)
  if (myrank .eq. root) print*, 'total time with MPI', dble(tpostmpi-tprempi)/dble(tcountrate)
  !
  stop
  !
end program MFDn_Transitions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
