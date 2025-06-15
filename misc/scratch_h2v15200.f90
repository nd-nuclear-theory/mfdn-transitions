      use SPbasis, only: OrbitalList_t
      type(OrbitalList_t) :: orblist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! loop over bra J-subspaces
      J2max = 2*maxval(orblist%j2_orb)
      do Tzab = 1, -1, -1
        ! Tzcd is uniquely constrained by operator quantum number
        Tzcd = Tzab - Tzop
        if (Tzcd < -1) cycle
        do Jab = 0, J2max
           do gab = 0, 1
              parab = (-1)**gab
              ! parcd is uniquely constrained by operator quantum number
              parcd = parab*Parop
              ! loop over ket J-subspaces
              do Jcd = Jab, (Jab+Jop)
                 ! loop over states in bra subspace
                 do orba = 1, norbt
                    do orbb = orba, norbt
                       if ((tz2_orb(orba)+tz2_orb(orbb))/2 /= Tzab) cycle ! Tz
                       if (orblist%pr_orb(orba)*orblist%pr_orb(orbb) /= parab) cycle ! parity
                       if ((orba == orbb) .and. (mod(Jab,2) /= 0)) cycle    ! antisymmetry
                       j2a = orblist%j2_orb(orba)
                       j2b = orblist%j2_orb(orbb)
                       if ((abs(j2a-j2b)/2 > Jab).or.(Jab > (j2a+j2b)/2)) cycle ! triangularity
                       wtab = orblist%wt_orb(orba)+orblist%wt_orb(orbb)
                       if (wtab > WT2max(Tzab)) cycle ! WTmax

                       ! loop over states in ket subspace
                       do orbc = 1, norbt
                          do orbd = orbc, norbt
                             if (.not.pairwiseless(orba, orbb,orbc, orbd)) cycle
                             if ((tz2_orb(orbc)+tz2_orb(orbd))/2 /= Tzcd) cycle ! Tz
                             if (orblist%pr_orb(orbc)*orblist%pr_orb(orbd) /= parcd) cycle ! parity
                             if ((orbc == orbd) .and. (mod(Jcd,2) /= 0)) cycle    ! antisymmetry
                             j2c = orblist%j2_orb(orbc)
                             j2d = orblist%j2_orb(orbd)
                             if ((abs(j2c-j2d)/2 > Jcd).or.(Jcd > (j2c+j2d)/2)) cycle ! triangularity
                             wtcd = orblist%wt_orb(orbc)+orblist%wt_orb(orbd)
                             if (wtcd > WT2max(Tzcd)) cycle ! WTmax

                             ! read and place matrix element
                             read(fh) matel

                             ! throw away matrix element if larger than target indexing
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! loop over bra J-subspaces
      J2max = 2*maxval(orblist%j2_orb)
      Tzab = 1
      Tzcd = Tzab - Tzop
      excluded orbitals: (orba > norb_p).and.(orba<orblist%norb_p).or.(orba-orblist%norb_p > norb_n).and.(orba > orblist%norb_p)
      do Jab = 0, J2max
         do gab = 0, 1
            parab = (-1)**gab
            ! parcd is uniquely constrained by operator quantum number
            parcd = parab*Parop
            ! loop over ket J-subspaces
            do Jcd = Jab, (Jab+Jop)
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
                     !   depending on Tzop, orbc can be either p or n
                     !   however, orbd must always be n in this module
                     do orbc = 1, norbt
                        do orbd = max(orbc,orblist%norb_p+1), norbt
                           if (.not.pairwiseless(orba, orbb, orbc, orbd)) cycle
                           if ((tz2_orb(orbc)+tz2_orb(orbd))/2 /= Tzcd) cycle ! Tz
                           if (orblist%pr_orb(orbc)*orblist%pr_orb(orbd) /= parcd) cycle ! parity
                           if ((orbc == orbd) .and. (mod(Jcd,2) /= 0)) cycle    ! antisymmetry
                           j2c = orblist%j2_orb(orbc)
                           j2d = orblist%j2_orb(orbd)
                           if ((abs(j2c-j2d)/2 > Jcd).or.(Jcd > (j2c+j2d)/2)) cycle ! triangularity
                           wtcd = orblist%wt_orb(orbc)+orblist%wt_orb(orbd)
                           if (wtcd > WT2max(Tzcd)) cycle ! WTmax

                           ! read and place matrix element
                           read(fh) matel
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
