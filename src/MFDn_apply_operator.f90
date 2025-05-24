!
! Calculates
!
! (1) if (obdme) Reduced One-Body Transition Densities
!
! (2) if (numTBtrans > 0) numTBtrans Two-Body Transition matrix elements
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program MFDn_Apply
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
  real(8) :: tstart, tbefore, tafter,              &
             attest, ttestmin, ttestavg, ttestmax, &
             tketread, tbraread, tobdme, ttbo,     &
             tketreadmax, tbrareadmax, tobdmemax, ttbomax
  integer :: nketread, nbraread
  integer :: subvec_ket_start, subvec_bra_start, num_submat, submat_ctr
  !
  integer, dimension(MaxNumKets) :: TwoJ_out=-1
  integer, dimension(MaxNumKets) :: n_out=0
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
  logical :: input
  logical :: normalize=.false.
  !
  integer(kind=2), dimension(:,:), allocatable :: groupidlist_ket
  integer(kind=2), dimension(:,:), allocatable :: groupidlist_bra
  integer, dimension(:), allocatable :: Mstateptr_ket, Mstateptr_bra
!   integer, dimension(:), allocatable :: blkgroup_ket, blkoffset_ket
!   integer, dimension(:), allocatable :: blkgroup_bra, blkoffset_bra
  ! integer, dimension(:), allocatable :: numblksNm_ket, numblksNm_bra
  integer, dimension(:), allocatable :: tileptr_ket, tileptr_bra
  !
  real(kind=4), dimension(:,:), allocatable :: amp_bra
  real(kind=4), dimension(:,:), allocatable :: amp_ket
  real(kind=8) :: tmp, factEn, factor
  real(kind=8), dimension(:), allocatable :: ketsum, brasum, normfac, &
       totR20_p, totM1L_p, totM1S_p, totE1Q_p, totE2Q_p, totF0,       &
       totR20_n, totM1L_n, totM1S_n, totE1Q_n, totE2Q_n, totGT
  real(kind=8) :: rt4pi = 3.54490775113405
  !
  ! File names
  character(LEN=15) :: resfilename = 'apply.res'
  character(LEN=15) :: outfilename = 'apply.out'
  character(LEN=17) :: datfilename = 'apply.input'
  character(:), allocatable :: mbfile, smwffile
  real, external :: myphase
  !
  namelist /transition_data/ fmass, hbomeg,         &
       TwoJ_ket, n_ket, TwoJ_out, n_out, normalize, &
       infofilename_bra, infofilename_ket,          &
       basisfilename_bra, basisfilename_ket,        &
       smwffilename_bra, smwffilename_ket,          &
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
  numTBtrans = 1
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
  numbras = numkets
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
  Delta_TwoMj = wf_info_bra%TwoMj - wf_info_ket%TwoMj
  Delta_Par = wf_info_bra%parity * wf_info_ket%parity
  !
  ! Read Two-Body matrix elements
  if (numTBtrans .ne. 1) then
     call cancelall(199)
  endif
  !
  ! Read Two-Body Matrix Elements (max. number set in preprocessor)
  call read_TBME_multi(numTBtrans, TBMEoperators)
  !
  ! Consistency checks
  TwoJop = 2 * Jop
  if (myrank .eq. root) then
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
    if (TwoJop /= abs(Delta_TwoMj)) then
       print*, 'WARNING TwoJop not equal to Delta_TwoMj; final state may not have good J'
    endif
    do k = 1, numkets
       if (TwoJ_out(k) < 0) then
          print*, 'invalid TwoJ_out', k, TwoJ_out(k)
          call cancelall(201)
       endif
       if (TwoJ_out(k) > (TwoJ_ket(k)+TwoJop)) then
          print*, 'TwoJ_out greater than TwoJ_ket + TwoJop', &
                  k, TwoJ_out(k), TwoJ_ket(k), TwoJop
          call cancelall(200)  ! same error code as MFDn_Transitions
       endif
       if (TwoJ_out(k) < abs(TwoJ_ket(k)-TwoJop)) then
          print*, 'TwoJ_out less than abs(TwoJ_ket - TwoJop)', &
                  k, TwoJ_out(k), TwoJ_ket(k), TwoJop
          call cancelall(202)  ! same error code as MFDn_Transitions
       endif
       if ((TwoJop > 0) .and. (TwoJ_ket(k) < wf_info_ket%TwoMj)) then
             print*, 'WARNING TwoJop>0 and TwoJ not equal to TwoMj;', &
                     ' final state may not have good J',              &
                     k, TwoJ_ket(k), wf_info_ket%TwoMj
       endif
    end do
  endif
  !
  oprank = 2
  !
  statesize = nparticles + oprank
  !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  ! write out new mfdn_smwf.info file
  if (myrank .eq. root) then
     print*, 'writing updated info to ', TRIM(infofilename_bra)
9901 format(i8, '    ! Version Number' )
9902 format(3i4,'    ! Z, N, 2Mj' )
     open(unit=99, file=TRIM(infofilename_bra), status='replace', action='write')
     write(99, 9901) 15200  ! versionno
     write(99, 9902) wf_info_bra%num_protons, wf_info_bra%num_neutrons, wf_info_bra%TwoMj
     write(99, *) TRIM(wf_info_bra%uniqueID)
9920 format(2i8, '    ! Number of proton and neutron orbitals')
9921 format(i4, 4x, 4i4, f12.6)
     write(99, *)
     write(99, 9920) wf_info_bra%orbital_list%norb_p, wf_info_bra%orbital_list%norb_n
     do k = 1, wf_info_bra%orbital_list%norb_p
        write(99, 9921) &
             k,                                      &
             wf_info_bra%orbital_list%n_orb(k),      &
             wf_info_bra%orbital_list%l_orb(k),      &
             wf_info_bra%orbital_list%j2_orb(k),     &
             +1,                                     & ! tz
             wf_info_bra%orbital_list%wt_orb(k)
     enddo
     do k = wf_info_bra%orbital_list%norb_p+1, wf_info_bra%orbital_list%norb_p+wf_info_bra%orbital_list%norb_n
        write(99, 9921) &
             k,                                      &
             wf_info_bra%orbital_list%n_orb(k),      &
             wf_info_bra%orbital_list%l_orb(k),      &
             wf_info_bra%orbital_list%j2_orb(k),     &
             -1,                                     & ! tz
             wf_info_bra%orbital_list%wt_orb(k)
     enddo
     write(99, *)
9930 format(2i8, '    ! Number of proton and neutron partitions')
     write(99, 9930) wf_info_bra%part_info%size_p, wf_info_bra%part_info%size_n
     write(99, *) wf_info_bra%part_info%partitions_p
     write(99, *) wf_info_bra%part_info%partitions_n
     write(99, *)
9911 format(i4,f12.4,i16, i8,'    ! Par, WTm, Dim, Npe')
     write(99, *)
     write(99, 9911) wf_info_bra%parity, &
                     wf_info_bra%WTmax, &
                     wf_info_bra%totalMdim, &
                     wf_info_bra%num_subvec
9992 format(i8,                                                                &
          '    ! n_states, followed by (i, 2J, nJ, T, -Eb, res)')
9999 format(3i8, f8.2, f16.4, e16.2)
     write(99, *)
     write(99, 9992) numbras
     do k = 1, numbras
        write(99, 9999) k, TwoJ_out(k), n_out(k), &
           wf_info_bra%isospin(k), wf_info_bra%eigval(k), wf_info_bra%residue(k)
     enddo
     write(99, *)
     close(unit=99, status='keep')
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank .eq. root) then
     open(unit=7, file=outfilename, status='unknown', action='write')
     write(7,*)
     write(7,*) '  Output from MFDn Apply in M-scheme basis'
     write(7,*)
     write(7,*) '  Revision ', VCS_REVISION
     write(7,*)
     write(7,*) ' Total number of particles ', nparticles
     write(7,*) ' Number of protons  <bra| and |ket>', wf_info_bra%num_protons, wf_info_ket%num_protons
     write(7,*) ' Number of neutrons <bra| and |ket>', wf_info_bra%num_neutrons, wf_info_ket%num_neutrons
     write(7,*)
     write(7,*) ' Parity         <bra| and |ket>', wf_info_bra%parity, wf_info_ket%parity
     write(7,*) ' Projection M_J <bra| and |ket>', 0.5*wf_info_bra%TwoMj, 0.5*wf_info_ket%TwoMj
     write(7,*)
     write(7,*) ' Number of output states ', numbras
     write(7,*) ' Total spin J, n(J), wfn(n) |bra>', 0.5*TwoJ_bra, n_bra, brastate
     write(7,*) ' Number of input states ', numkets
     do k = 1, numkets
        write(7,*) ' Total spin J, n(J), wfn(n) |ket>', 0.5*TwoJ_ket(k), n_ket(k), ketstate(k)
     end do
     write(7,*)
     if (numTBtrans .gt. 0) then
        write(7,*) ' Applying Two-Body Operator'
        write(7,*)
        write(7,*) ' Spin of Operator', Jop
        write(7,*) '  M_J of Operator', Delta_TwoMj/2
        write(7,*) '  Parity of Oprtr', Parop
        write(7,*)
     endif
     !
     close(unit=7, status='keep')
     !
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
  if (nprocs .ne. wf_info_bra%num_subvec) then
     print*, 'nprocs not equal to number of bra subvecs', nprocs, wf_info_bra%num_subvec
     call cancelall(567)
  endif
  !
  if (myrank .eq. root) then
     open(unit=7, file=outfilename, status='old', action='write', position='append')
     !
     write(7, *) ' Number of MPI ranks           ', nprocs
     write(7, *)
     write(7, *) ' Number of OMP threads per rank', numthreads
     write(7, *)
     !
     close(unit=7, status='keep')
  endif
  !
  allocate(ketsum(numkets), brasum(numbras))
  dim_ket = 0
  ketsum(1:numkets) = 0.d0
  dim_bra = 0
  brasum(1:numkets) = 0.d0
  submat_ctr = 0
  tketread = 0.d0
  tbraread = 0.d0
  tobdme = 0.d0
  ttbo = 0.d0
  tketreadmax = 0.d0
  tbrareadmax = 0.d0
  tobdmemax = 0.d0
  ttbomax = 0.d0
  nketread = 0
  nbraread = 0
  !
  ! start timing for row read
  tbefore = MPI_WTIME()
  !
  !!print*, "rank", myrank, "starting block", subvec_ket, subvec_bra
  !
  subvec_bra = myrank
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
!   allocate(blkgroup_bra(numblks_bra+1))
  read(fh) ! blkgroup_bra(1:numblks_bra+1)  ! not used
  !
!   allocate(blkoffset_bra(numblks_bra+1))
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
  !
  ! Initialize wavefunction < bra |
  allocate(amp_bra(numkets, numstates_bra))
  amp_bra(1:numkets, 1:numstates_bra) = 0.d0
  dim_bra = dim_bra + numstates_bra
  !
  allocate(tileptr_bra(numgroupids_bra+1))
  !
  tafter = MPI_WTIME()
  tbraread = tbraread + tafter-tbefore
  tbrareadmax = max(tbrareadmax, tafter-tbefore)
  nbraread = nbraread + 1
  ! Loop over column (| ket >, psi) ranks
  ketloop: do subvec_ket = 0, wf_info_ket%num_subvec-1
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
     maxnMstates = max(maxnMstates_ket, maxnMstates_bra)
     !
     ! Read wavefunction | ket >
     allocate(amp_ket(numkets, numstates_ket))
     smwffile = TRIM(smwffilename_ket)//achar(ifile/100+48)&
          //achar(mod(ifile/10,10)+48)//achar(mod(ifile,10)+48)
     print*, 'reading ', smwffile
     open(unit=fh, file=TRIM(smwffile), status='unknown', action='read', form='unformatted')
     k = 0
     ketread: do i = 1, wf_info_ket%num_states
        do j = 1, numkets
           if (ketstate(j) == i) then
              read(fh) amp_ket(j, 1:numstates_ket)
              print*, j, numstates_ket
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! NOTE: difference with MFDn -- TwoMj could be different for < bra | and | ket >
     tbefore = MPI_WTIME()
     call generateSparsity(nparticles, oprank, nspstates, classoffset(2), &
          wf_info_ket%num_protons, wf_info_ket%num_neutrons, numgroupids_ket, groupidlist_ket, &
          wf_info_bra%num_protons, wf_info_bra%num_neutrons, numgroupids_bra, groupidlist_bra, &
          tileptr_ket)                  ! also sets tileind and tilediff
     !
     call ApplyTBop(maxnMstates,                                     &
          wf_info_ket%num_protons, wf_info_ket%num_neutrons, numgroupids_ket, groupidlist_ket, &
          wf_info_ket%TwoMj, Mstateptr_ket, numstates_ket,                   &
          wf_info_bra%num_protons, wf_info_bra%num_neutrons, numgroupids_bra, groupidlist_bra, &
          wf_info_bra%TwoMj, Mstateptr_bra, numstates_bra,                   &
          tileptr_ket, ntiles, tileind, tilediff, numkets, amp_ket, amp_bra)
     !
     deallocate(tilediff, tileind)       ! allocated in generateSparsity
     ttbo = ttbo + MPI_WTIME()-tbefore

     deallocate(groupidlist_ket)
     ! deallocate(blkgroup_ket)
     ! deallocate(blkoffset_ket)
     deallocate(amp_ket)
     deallocate(tileptr_ket)
     deallocate(Mstateptr_ket)
  enddo ketloop

  ! check norm on output wave functions
  !$omp parallel do simd collapse(2) default(none)          &
  !$omp             shared(amp_bra, numbras, numstates_bra) &
  !$omp             reduction(+: brasum)
  do j = 1, numstates_bra
     do k = 1, numbras
        brasum(k) = brasum(k) + amp_bra(k,j)**2
     enddo
  enddo
  !$omp end parallel do simd
  if (myrank .eq. root) then
     tbefore = MPI_WTIME()
     call MPI_AllReduce(MPI_IN_PLACE, brasum, numbras, &
          MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
     call MPI_AllReduce(MPI_IN_PLACE, ketsum, numkets, &
          MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
     print*, 'time collecting', MPI_WTIME()-tbefore
  else
     call MPI_AllReduce(MPI_IN_PLACE, brasum, numbras, &
          MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
     call MPI_AllReduce(MPI_IN_PLACE, ketsum, numkets, &
          MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
  endif


  if (myrank .eq. root) then
     tbefore = MPI_WTIME()
     open(unit=7, file=outfilename, status='old', action='write', position='append')
     !
     do k = 1, numkets
        if (abs(ketsum(k)-1.d0) > 1d-6) then
           write(7, *) 'WARNING: input norm incorrect', k, ketsum(k)
           print*, 'WARNING: input norm incorrect', k, ketsum(k)
        endif
     end do
  endif

  if (normalize) then
     allocate(normfac(numbras))
     if (myrank .eq. root) print*, 'normalizing output vectors'
     normfac(1:numbras) = 1./sqrt(brasum(1:numbras))
     if (myrank .eq. root) print*, normfac
     !$omp parallel do simd collapse(2) default(none)                        &
     !$omp                  shared(amp_bra, normfac, numbras, numstates_bra)
     do j = 1, numstates_bra
        do k = 1, numbras
           amp_bra(k, j) = amp_bra(k, j) * normfac(k)
        end do
     end do
     deallocate(normfac)
  else if (myrank .eq. root) then
     do k = 1, numbras
        if (abs(brasum(k)-1.d0) > 1d-6) then
           write(7, *) 'WARNING: output norm incorrect', k, brasum(k)
           print*, 'WARNING: output norm incorrect', k, brasum(k)
        endif
     end do
  endif

  smwffile = TRIM(smwffilename_bra)//achar(jfile/100+48)&
       //achar(mod(jfile/10,10)+48)//achar(mod(jfile,10)+48)
  print*, 'writing ', smwffile
  open(unit=fh, file=TRIM(smwffile), status='unknown', action='write', form='unformatted')
  do j = 1, numbras
     print*, j, numstates_bra
     write(fh) amp_bra(j, 1:numstates_bra)
  end do
  close(fh, status='keep')
!   ! Test on norm
!   do j = 1, numstates_bra
!      brasum = brasum + amp_bra(j)**2
!   enddo

  deallocate(Mstateptr_bra)
  deallocate(groupidlist_bra)
!   deallocate(blkgroup_bra)
!   deallocate(blkoffset_bra)
  deallocate(amp_bra)
  deallocate(tileptr_bra)

  if (myrank .eq. root) print*, 'total time', MPI_WTIME()-tstart
  call MPI_Finalize(ierr)
  !
  call system_clock(tpostmpi, tcountrate)
  if (myrank .eq. root) print*, 'total time with MPI', dble(tpostmpi-tprempi)/dble(tcountrate)
  !
  stop
  !
end program MFDn_Apply

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
