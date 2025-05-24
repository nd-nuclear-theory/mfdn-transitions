
module MjStates
  implicit none
  private
  public MjStatesCount, MjStatesGen
  !
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine MjStatesCount(nparticles, nspstates, mj2_sp, twomj, nextbin, groupID, numstates)
    !
    ! Given a groupID, this counts the number of many-body states
    ! belonging to that groupID that satisfy the given mj constraint
    !
    integer, intent(in) :: nparticles, nspstates, twomj
    integer, dimension(nspstates), intent(in) :: mj2_sp, nextbin
    integer(kind=2), dimension(nparticles), intent(in) :: groupID
    integer, intent(out) :: numstates
    ! local variables
    integer, dimension(nparticles) :: mbstate, mbgroup
    integer :: flag
    !
    numstates = 0
    mbgroup(1:nparticles) = groupID(1:nparticles)
    mbstate(1:nparticles) = mbgroup(1:nparticles)
    !      
    call setlastmj(nparticles, nspstates, mj2_sp, twomj, nextbin, mbstate, flag)
    if (flag .eq. 0) numstates = numstates + 1
    !
    flag = 0
    do while (flag .eq. 0)
       flag = -1
       do while (flag .eq. -1)
          call incrementmj(nparticles, nspstates, mj2_sp, twomj, nextbin, mbgroup, mbstate, flag)
       enddo
       if (flag .eq. 0) numstates = numstates + 1
    enddo
    !
    return
  end subroutine MjStatesCount
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine MjStatesGen(nparticles, nspstates, mj2_sp, twomj, nextbin, groupID, numstates, mbstatelist)
    !
    ! Given a groupID, this generates all many-body states
    ! belonging to that groupID that satisfy the given mj constraint
    !
    integer, intent(in) :: nparticles, nspstates, twomj
    integer, intent(inout) :: numstates
    integer, dimension(nspstates), intent(in) :: mj2_sp, nextbin
    integer(kind=2), dimension(nparticles), intent(in) :: groupID
    integer(kind=2), dimension(nparticles, numstates), intent(out) :: mbstatelist
    ! local variables
    integer, dimension(nparticles) :: mbstate, mbgroup
    integer :: currentstate, flag
    !
    currentstate = 0
    mbgroup(1:nparticles) = groupID(1:nparticles)
    mbstate(1:nparticles) = mbgroup(1:nparticles)
    !      
    call setlastmj(nparticles, nspstates, mj2_sp, twomj, nextbin, mbstate, flag)
    !
    if (flag .eq. 0) then
       currentstate = currentstate + 1
       mbstatelist(1:nparticles, currentstate) = mbstate(1:nparticles)
    endif
    !
    flag = 0
    do while (flag .eq. 0)
       flag = -1
       do while (flag .eq. -1)
          call incrementmj(nparticles, nspstates, mj2_sp, twomj, nextbin, mbgroup, mbstate, flag)
       enddo
       if (flag .eq. 0) then
          currentstate = currentstate + 1
          mbstatelist(1:nparticles, currentstate) = mbstate(1:nparticles)
       endif
    enddo
    !
    numstates = currentstate
    return
  end subroutine MjStatesGen
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setlastmj(nparticles, nspstates, mj2_sp, twomj, nextbin, mbstate, flag)
    !
    ! Attempts to change the last particle in the many-body state to
    ! satisfy the mj constraint within the same group of many-body states.  
    ! On exit, 
    !   flag = 1  : have pulled too much mj from last state
    !   flag = 0  : success!
    !   flag = -1 : need to pull more mj from last state
    !
    ! NOTE: this code relies on complete set (m_j) values in each orbital
    !   and order from low to high:   
    !   m_j = -j, -j+1, ....., j-1, j 
    !
    integer, intent(in) :: nparticles, nspstates, twomj
    integer, dimension(nspstates), intent(in) :: mj2_sp, nextbin
    integer, dimension(nparticles), intent(inout) :: mbstate
    integer, intent(out) :: flag
    ! local variables
    integer :: i, deltamj, ilast
    ! 
   deltamj = twomj
    do i = 1, nparticles
       deltamj = deltamj - mj2_sp(mbstate(i))
    enddo
    !
    if (deltamj .lt. 0) then
       flag = 1
       return
    endif
    !
    ilast = mbstate(nparticles) + deltamj / 2
    if (ilast .lt. nextbin(mbstate(nparticles))) then  
       mbstate(nparticles) = ilast
       flag = 0
    else
       flag = -1
    endif
    !
    return
  end subroutine setlastmj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine incrementmj(nparticles, nspstates, mj2_sp, twomj, nextbin, mbgroup, mbstate, flag)
    !
    ! Given a many-body state that satisfies the mj constraint,
    ! this replaces it with the lexicographically next state
    ! within the same mbgroup that also satisfies the mj constraint.
    !
    integer, intent(in) :: nparticles, nspstates, twomj
    integer, dimension(nspstates), intent(in) :: mj2_sp, nextbin 
    integer, dimension(nparticles), intent(in) :: mbgroup
    integer, dimension(nparticles), intent(inout) :: mbstate
    integer, intent(inout) :: flag
    ! local variables
    integer :: i, j
    !
    do i = nparticles - 1, 1, -1
       if (mbstate(i) .lt. (nextbin(mbgroup(i))-1)) then
          mbstate(i) = mbstate(i) + 1
          do j = i + 1, nparticles
             if ((mbstate(j-1)+1) .ge. nextbin(mbgroup(j))) then
                flag = 2
                exit
             else
                mbstate(j) = mbstate(j-1) + 1
                mbstate(j) = max(mbstate(j), mbgroup(j))
             endif
          enddo
          if (flag .eq. 2) then
             flag = 0
             cycle
          endif
          call setlastmj(nparticles, nspstates, mj2_sp, twomj, nextbin, mbstate, flag)
          if (flag .eq. 1) then
             cycle
          else
             return
          endif
       endif
    enddo
    !
    flag = 1
    return      
  end subroutine incrementmj
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module MjStates
