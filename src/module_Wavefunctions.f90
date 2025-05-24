module Wavefunctions_m
  use SPbasis, only: OrbitalList_t, PartitioningInfo_t, setupSPbasis
  implicit none
  private
  public InitializeWavefunctions
  public WaveFunctionInfo_t
  public wf_info_bra, wf_info_ket

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! diagonalization parameters
  type :: WaveFunctionInfo_t
     ! basis parameters
     character(len=129) :: uniqueID
     type(OrbitalList_t) :: orbital_list
     type(PartitioningInfo_t) :: part_info
     integer :: num_protons, num_neutrons, TwoMj, parity
     real(kind=4) :: WTmax
     integer(kind=8) :: totalMdim
     integer :: num_subvec
     ! state parameters
     integer :: num_states
     integer, allocatable, dimension(:) :: TwoJ
     integer, allocatable, dimension(:) :: Jseq
     real(kind=4), allocatable, dimension(:) :: isospin
     real(kind=4), allocatable, dimension(:) :: eigval
     real(kind=4), allocatable, dimension(:) :: residue
   contains
     procedure :: InitFromFile => WaveFunctionInfo_InitFromFile
  end type WaveFunctionInfo_t
  !
  ! wave function parameters
  type(WaveFunctionInfo_t), protected :: wf_info_bra, wf_info_ket

contains

    subroutine InitializeWavefunctions(info_filename_bra, info_filename_ket)
        use nodeinfo
        character(len=*), intent(in) :: info_filename_bra, info_filename_ket
        ! local temporaries
        type(OrbitalList_t) :: orb_list_bra, orb_list_ket, orbital_list
        type(PartitioningInfo_t) :: part_info_bra, part_info_ket, part_info
        integer :: j, k, offset_n, norb_tot
        logical :: b

        call wf_info_bra%InitFromFile(info_filename_bra)
        orb_list_bra = wf_info_bra%orbital_list
        part_info_bra = wf_info_bra%part_info
        call wf_info_ket%InitFromFile(info_filename_ket)
        orb_list_ket = wf_info_ket%orbital_list
        part_info_ket = wf_info_ket%part_info

        if (myrank .eq. root) then
            ! compare orbitals and partitioning, check that they are equal
            if (orb_list_bra%norb_p /= orb_list_ket%norb_p) then
                print *, 'WARNING unequal number of proton orbitals', &
                        orb_list_bra%norb_p, orb_list_ket%norb_p
            end if
            if (orb_list_bra%norb_n /= orb_list_ket%norb_n) then
                print *, 'WARNING unequal number of neutron orbitals', &
                        orb_list_bra%norb_n, orb_list_ket%norb_n
            end if
            do j = 1, min(orb_list_bra%norb_p, orb_list_ket%norb_p)
                b = (orb_list_bra%n_orb(j) == orb_list_ket%n_orb(j))
                b = b .and. (orb_list_bra%l_orb(j) == orb_list_ket%l_orb(j))
                b = b .and. (orb_list_bra%j2_orb(j) == orb_list_ket%j2_orb(j))
                b = b .and. (orb_list_bra%pr_orb(j) == orb_list_ket%pr_orb(j))
                b = b .and. (orb_list_bra%wt_orb(j) - orb_list_ket%wt_orb(j)) <= 1.d-6
                if (.not. b) then
                    print *, 'ERROR non-matching proton orbitals'
                    print *, 'bra:',                &
                            orb_list_bra%n_orb(j),  &
                            orb_list_bra%l_orb(j),  &
                            orb_list_bra%j2_orb(j), &
                            orb_list_bra%pr_orb(j), &
                            orb_list_bra%wt_orb(j)
                    print *, 'ket:',                &
                            orb_list_ket%n_orb(j),  &
                            orb_list_ket%l_orb(j),  &
                            orb_list_ket%j2_orb(j), &
                            orb_list_ket%pr_orb(j), &
                            orb_list_ket%wt_orb(j)
                    call cancelall(233)
                end if
            end do
            do j = 1, min(orb_list_bra%norb_n, orb_list_ket%norb_p)
                b =         (orb_list_bra%n_orb(j+orb_list_bra%norb_p)             &
                                == orb_list_ket%n_orb(j+orb_list_ket%norb_p))
                b = b .and. (orb_list_bra%l_orb(j+orb_list_bra%norb_p)             &
                                == orb_list_ket%l_orb(j+orb_list_ket%norb_p))
                b = b .and. (orb_list_bra%j2_orb(j+orb_list_bra%norb_p)            &
                                == orb_list_ket%j2_orb(j+orb_list_ket%norb_p))
                b = b .and. (orb_list_bra%pr_orb(j+orb_list_bra%norb_p)            &
                                == orb_list_ket%pr_orb(j+orb_list_ket%norb_p))
                b = b .and. (orb_list_bra%wt_orb(j+orb_list_bra%norb_p)            &
                                - orb_list_ket%wt_orb(j+orb_list_ket%norb_p)) <= 1.d-6
                if (.not. b) then
                    print *, 'ERROR non-matching neutron orbitals'
                    print *, 'bra:',                                      &
                            orb_list_bra%n_orb(j+orb_list_bra%norb_p),    &
                            orb_list_bra%l_orb(j+orb_list_bra%norb_p),    &
                            orb_list_bra%j2_orb(j+orb_list_bra%norb_p),   &
                            orb_list_bra%pr_orb(j+orb_list_bra%norb_p),   &
                            orb_list_bra%wt_orb(j+orb_list_bra%norb_p)
                    print *, 'ket:',                                      &
                            orb_list_ket%n_orb(j+orb_list_ket%norb_p),    &
                            orb_list_ket%l_orb(j+orb_list_ket%norb_p),    &
                            orb_list_ket%j2_orb(j+orb_list_ket%norb_p),   &
                            orb_list_ket%pr_orb(j+orb_list_ket%norb_p),   &
                            orb_list_ket%wt_orb(j+orb_list_ket%norb_p)
                    call cancelall(233)
                end if
            end do
        endif

        orbital_list%norb_p = max(orb_list_bra%norb_p, orb_list_ket%norb_p)
        orbital_list%norb_n = max(orb_list_bra%norb_n, orb_list_ket%norb_n)
        norb_tot = orbital_list%norb_p + orbital_list%norb_n
        allocate(                           &
            orbital_list%j2_orb(norb_tot),  &
            orbital_list%pr_orb(norb_tot),  &
            orbital_list%n_orb(norb_tot),   &
            orbital_list%l_orb(norb_tot),   &
            orbital_list%wt_orb(norb_tot)   &
            )

        ! cobble together an orbital list
        !
        !   first we copy the larger of the two proton lists
        offset_n = orb_list_bra%norb_p - orb_list_ket%norb_p
        if (orb_list_bra%norb_p >= orb_list_ket%norb_p) then
            ! if bra's proton orbital list is bigger, copy it
            do j = 1, orbital_list%norb_p
                orbital_list%n_orb(j)  = orb_list_bra%n_orb(j)
                orbital_list%l_orb(j)  = orb_list_bra%l_orb(j)
                orbital_list%j2_orb(j) = orb_list_bra%j2_orb(j)
                orbital_list%pr_orb(j) = orb_list_bra%pr_orb(j)
                orbital_list%wt_orb(j) = orb_list_bra%wt_orb(j)
            end do
        else ! orb_list_bra%norb_p < orb_list_ket%norb_p
            ! if the ket's proton orbital list is bigger, copy it
            do j = 1, orbital_list%norb_p
                orbital_list%n_orb(j)  = orb_list_ket%n_orb(j)
                orbital_list%l_orb(j)  = orb_list_ket%l_orb(j)
                orbital_list%j2_orb(j) = orb_list_ket%j2_orb(j)
                orbital_list%pr_orb(j) = orb_list_ket%pr_orb(j)
                orbital_list%wt_orb(j) = orb_list_ket%wt_orb(j)
            end do
        endif
        !   next we copy the larger of the two neutron lists
        if (orb_list_bra%norb_n >= orb_list_ket%norb_n) then
            ! if bra's neutron orbital list is bigger, copy it
            do j = 1, orbital_list%norb_n
                orbital_list%n_orb(orbital_list%norb_p+j)         &
                    = orb_list_bra%n_orb(orb_list_bra%norb_p+j)
                orbital_list%l_orb(orbital_list%norb_p+j)         &
                    = orb_list_bra%l_orb(orb_list_bra%norb_p+j)
                orbital_list%j2_orb(orbital_list%norb_p+j)        &
                    = orb_list_bra%j2_orb(orb_list_bra%norb_p+j)
                orbital_list%pr_orb(orbital_list%norb_p+j)        &
                    = orb_list_bra%pr_orb(orb_list_bra%norb_p+j)
                orbital_list%wt_orb(orbital_list%norb_p+j)        &
                    = orb_list_bra%wt_orb(orb_list_bra%norb_p+j)
            end do
        else ! orb_list_bra%norb_n < orb_list_ket%norb_n
            ! if the ket's neutron orbital list is bigger, copy it
            do j = 1, orbital_list%norb_n
                orbital_list%n_orb(orbital_list%norb_p+j)         &
                    = orb_list_ket%n_orb(orb_list_ket%norb_p+j)
                orbital_list%l_orb(orbital_list%norb_p+j)         &
                    = orb_list_ket%l_orb(orb_list_ket%norb_p+j)
                orbital_list%j2_orb(orbital_list%norb_p+j)        &
                    = orb_list_ket%j2_orb(orb_list_ket%norb_p+j)
                orbital_list%pr_orb(orbital_list%norb_p+j)        &
                    = orb_list_ket%pr_orb(orb_list_ket%norb_p+j)
                orbital_list%wt_orb(orbital_list%norb_p+j)        &
                    = orb_list_ket%wt_orb(orb_list_ket%norb_p+j)
            end do
        endif

        if (myrank .eq. root) then
            if (part_info_bra%size_p /= part_info_ket%size_p) then
                print *, 'WARNING unequal number of proton partitions', &
                    part_info_bra%size_p, part_info_ket%size_p
                ! call cancelall(235)
            end if
            if (part_info_bra%size_n /= part_info_ket%size_n) then
                print *, 'WARNING unequal number of neutron partitions', &
                    part_info_bra%size_n, part_info_ket%size_n
                ! call cancelall(236)
            end if
            do j = 1, min(part_info_bra%size_p, part_info_ket%size_p)
                if (part_info_bra%partitions_p(j) /= part_info_ket%partitions_p(j)) then
                    print *, 'ERROR non-matching proton partitions', j
                    print *, 'bra:', part_info_bra%partitions_p(j)
                    print *, 'ket:', part_info_ket%partitions_p(j)
                    call cancelall(237)
                end if
            end do
            do j = 1, min(part_info_bra%size_n, part_info_ket%size_n)
                if (part_info_bra%partitions_n(j) /= part_info_ket%partitions_n(j)) then
                    print *, 'ERROR non-matching neutron partitions', j
                    print *, 'bra:', part_info_bra%partitions_n(j)
                    print *, 'ket:', part_info_ket%partitions_n(j)
                    call cancelall(238)
                end if
            end do
        endif

        ! copy the larger partitioning
        part_info%size_p = max(part_info_bra%size_p, part_info_ket%size_p)
        allocate(part_info%partitions_p(part_info%size_p))
        part_info%size_n = max(part_info_bra%size_n, part_info_ket%size_n)
        allocate(part_info%partitions_n(part_info%size_n))

        if (part_info_bra%size_p >= part_info_ket%size_p) then
            do j = 1, part_info_bra%size_p
                part_info%partitions_p(j) = part_info_bra%partitions_p(j)
            end do
        else
            do j = 1, part_info_ket%size_p
                part_info%partitions_p(j) = part_info_ket%partitions_p(j)
            end do
        end if

        if (part_info_bra%size_n >= part_info_ket%size_n) then
            do j = 1, part_info_bra%size_n
                part_info%partitions_n(j) = part_info_bra%partitions_n(j)
            end do
        else
            do j = 1, part_info_ket%size_n
                part_info%partitions_n(j) = part_info_ket%partitions_n(j)
            end do
        end if


        ! initialize single-particle basis and partitioning
        call setupSPbasis(orbital_list, part_info)

        return
    end subroutine InitializeWavefunctions

    subroutine WaveFunctionInfo_InitFromFile(self, info_filename)
        class(WaveFunctionInfo_t), intent(inout) :: self
        character(len=*), intent(in) :: info_filename

        ! local temporaries
        integer :: fh, j, k, version_number
        integer :: a, na, la, j2a, tz2a, norb_tot
        real(kind=4) :: wta

        fh=51
        open(unit=fh, file=info_filename, status='old', action='read')

        read(fh, *) version_number
        if (version_number .eq. 15200) then
            read(fh, *) self%num_protons, self%num_neutrons, self%TwoMj
            read(fh, *) self%uniqueID

            ! read orbital definitions
            read(fh, *) self%orbital_list%norb_p, self%orbital_list%norb_n
            norb_tot = self%orbital_list%norb_p + self%orbital_list%norb_n
            allocate( &
                self%orbital_list%j2_orb(norb_tot),                            &
                self%orbital_list%pr_orb(norb_tot),                            &
                self%orbital_list%n_orb(norb_tot),                             &
                self%orbital_list%l_orb(norb_tot),                             &
                self%orbital_list%wt_orb(norb_tot)                             &
                )
            do j = 1, norb_tot
                read(fh, *) a, na, la, j2a, tz2a, wta
                if ((tz2a == +1) .and. &
                    ( .not. (j <= self%orbital_list%norb_p))) then
                    print*, 'invalid orbital', a, na, la, j2a, tz2a, wta
                    call cancelall(241)
                else if ((tz2a == -1) .and. &
                         ( .not. (j > self%orbital_list%norb_p .and. j <= norb_tot))) then
                    print*, 'invalid orbital', a, na, la, j2a, tz2a, wta
                    call cancelall(242)
                endif
                self%orbital_list%n_orb(a)  = na
                self%orbital_list%l_orb(a)  = la
                self%orbital_list%j2_orb(a) = j2a
                self%orbital_list%pr_orb(a) = (-1)**la
                self%orbital_list%wt_orb(a) = wta
            enddo

            ! read partitioning info
            read(fh, *) self%part_info%size_p, self%part_info%size_n
            allocate(self%part_info%partitions_p(self%part_info%size_p))
            allocate(self%part_info%partitions_n(self%part_info%size_n))
            read(fh, *) self%part_info%partitions_p(1:self%part_info%size_p)
            ! minimal check partitioning for consistency
            do j = 1, self%part_info%size_p-1
               if (self%part_info%partitions_p(j+1) <= self%part_info%partitions_p(j)) then
                  print *, "ERROR in input partitioning"
                  print *, self%part_info%partitions_p(j), self%part_info%partitions_p(j+1)
                  print *, "should be monotonic"
                  call cancelall(243)
               endif
            enddo
            read(fh, *) self%part_info%partitions_n(1:self%part_info%size_n)
            do j = 1, self%part_info%size_n-1
                if (self%part_info%partitions_n(j+1) <= self%part_info%partitions_n(j)) then
                   print *, "ERROR in input partitioning"
                   print *, self%part_info%partitions_n(j), self%part_info%partitions_n(j+1)
                   print *, "should be monotonic"
                   call cancelall(243)
                endif
             enddo

            ! read run info
            read(fh, *) self%parity, self%WTmax, self%totalMdim, self%num_subvec
            read(fh, *) self%num_states
            allocate(self%TwoJ(1:self%num_states))
            allocate(self%Jseq(1:self%num_states))
            allocate(self%eigval(1:self%num_states))
            allocate(self%isospin(1:self%num_states))
            allocate(self%residue(1:self%num_states))
            do j = 1, self%num_states
                read(fh, *) k, self%TwoJ(j), self%Jseq(j), self%isospin(j), &
                        self%eigval(j), self%residue(j)
                if (j /= k) then
                    print *, 'ERROR reading eigenstates', k, self%TwoJ(j), &
                        self%Jseq(j), self%isospin(j), self%eigval(j), self%residue(j)
                    call cancelall(244)
                end if
            end do

        else
            print*, 'unknown wave function info file version', version_number
            call cancelall(251)
        end if
        close(fh)

        return
    end subroutine WaveFunctionInfo_InitFromFile
end module Wavefunctions_m
