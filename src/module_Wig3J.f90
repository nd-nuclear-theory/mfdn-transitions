!
module Wigner3J
  implicit none
  private
  public Wig3J, Jmax_2body
  public Init_Wigner3J
  public Retrieve_Wigner3J, Retrieve_Wigner3J_Array
  public Init_Wigner3J_Array_JJK, Wig3J_Array_JJK
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  integer :: Jmax_2body, Nj2
  real(kind=8), dimension(:), allocatable :: Wigner3Jarray_hhi
  real(kind=8), dimension(:, :, :), allocatable :: Wig3J_Array_JJK
  !
contains
  !
  ! initialization of Wigner arrays, specific for MFDn
  subroutine Init_Wigner3J(Jtot1mx, Jtot2mx)
    !
    ! Jtot1mx = (J2totmax(1)+1)/2  ! highest J+1/2 value for 1-body system
    ! Jtot2mx =  J2totmax(2) / 2   ! highest J value for 2-body system
    !
    integer, intent(in) :: Jtot1mx, Jtot2mx
    ! local variables
    integer :: j1, j2, j3, m1, m2
    integer :: index, Nj1m1, Nj2m2, isize
    !
    Jmax_2body = Jtot2mx
    Nj2 = Jtot1mx*(Jtot1mx + 1)
    !
    isize = Nj2*Nj2*(Jmax_2body + 1)
    !
    ! allocate array at appropriate size
    allocate (wigner3jarray_hhi(0:isize))
    wigner3jarray_hhi = 0.d0
    !
    ! size for Jtot1mx (= nshell) = 20, Jmax_2body = 20
    ! memory (double precision): 30 MB
    ! (in principle, further reduction by factor of 2 possible
    ! but that complicates the retrieval process significantly)
    !
    Nj1m1 = 0
    do j1 = 1, Jtot1mx
      do m1 = -j1 + 1, j1
        Nj2m2 = 0
        do j2 = 1, Jtot1mx
          do m2 = -j2 + 1, j2
            do j3 = 0, Jmax_2body
              index = j3 + (Jmax_2body + 1)* &
                      (Nj2m2 + (m2 + j2 - 1) + Nj2*(Nj1m1 + (m1 + j1 - 1)))
              wigner3jarray_hhi(index) = &
                wig3j(2*j1 - 1, 2*j2 - 1, 2*j3, 2*m1 - 1, 2*m2 - 1, -(2*m1 + 2*m2 - 2))
              !
              ! print*,index, wigner3jarray_hhi(index)
              !
            end do
          end do
          Nj2m2 = Nj2m2 + 2*j2
        end do
      end do
      Nj1m1 = Nj1m1 + 2*j1
    end do
    !
    return
  end subroutine Init_Wigner3J
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! retrieve single Wigner 3J element from array
  subroutine Retrieve_Wigner3J(jj1, jj2, jj3, mm1, mm2, mm3, wignervalue)
    !
    ! j1, j2, m1, m2 -- twice half-integer j and m values
    ! j3             -- twice integer j value
    ! output: wigner 3-j value (double precision)
    !
    integer, intent(in) :: jj1, jj2, jj3, mm1, mm2, mm3
    real(kind=8), intent(out) :: wignervalue
    ! local variables
    integer :: Nj1m1, Nj2m2, index
    integer :: j1, j2, j3, m1, m2
    !
    j1 = jj1 + 1
    j1 = j1/2
    j2 = jj2 + 1
    j2 = j2/2
    j3 = jj3/2
    m1 = mm1 + 1
    m1 = m1/2
    m2 = mm2 + 1
    m2 = m2/2
    !
    Nj1m1 = j1*(j1 - 1)
    Nj2m2 = j2*(j2 - 1)
    index = j3 + (Jmax_2body + 1)*(Nj2m2 + (m2 + j2 - 1) + Nj2*(Nj1m1 + (m1 + j1 - 1)))
    wignervalue = wigner3jarray_hhi(index)
    ! wignervalue = wig3j(jj1, jj2, jj3, mm1, mm2, mm3)
    !
    return
  end subroutine Retrieve_Wigner3J
  !
  !
  ! retrieve set of neighboring Wigner 3J elements from array
  subroutine Retrieve_Wigner3J_Array(jj1, jj2, mm1, mm2, j3min, j3max, wignerarray)
    !
    ! j1, j2, m1, m2 -- twice half-integer j and m values
    ! j3min, j3max   -- range of j3, twice integer j values
    ! output: array of length j3min-j3max+1 with Wigner 3J elements
    !
    implicit none
    integer, intent(in) :: jj1, jj2, mm1, mm2, j3min, j3max
    real(kind=8), intent(out) :: wignerarray(1 + j3max - j3min)
    !
    ! local variables
    integer :: Nj1m1, Nj2m2, index
    integer :: j1, j2, m1, m2
    !
    j1 = jj1 + 1
    j1 = j1/2
    j2 = jj2 + 1
    j2 = j2/2
    m1 = mm1 + 1
    m1 = m1/2
    m2 = mm2 + 1
    m2 = m2/2
    !
    Nj1m1 = j1*(j1 - 1)
    Nj2m2 = j2*(j2 - 1)
    index = (Jmax_2body + 1)*(Nj2m2 + (m2 + j2 - 1) + Nj2*(Nj1m1 + (m1 + j1 - 1)))
    wignerarray(1:1 + j3max - j3min) = wigner3jarray_hhi(j3min + index:j3max + index)
    !
    return
  end subroutine Retrieve_Wigner3J_Array
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  ! initialization of Wigner arrays, specific for Post-Processor
  subroutine Init_Wigner3J_Array_JJK(Jop, Mjop)
    !
    integer, intent(in) :: Jop, Mjop
    ! local variables
    integer :: j1, j2, m1, m2, idj
    !
    ! allocate array (over-allocated, but convenient for now)
    allocate (Wig3J_Array_JJK(-Jop:Jop, 0:Jmax_2body, -Jmax_2body:Jmax_2body))
    !
    Wig3J_Array_JJK = 0.d0
    !
    do j1 = 0, Jmax_2body
      do m1 = -j1, j1
        do idj = -Jop, Jop
          j2 = j1 + idj
          m2 = m1 - Mjop
          !
          if (j2 .lt. 0) cycle
          if (abs(m2) .gt. j2) cycle
          !
          Wig3J_Array_JJK(idj, j1, m1) = &
            wig3j(2*j1, 2*Jop, 2*j2, -2*m1, 2*Mjop, 2*m2)
          !
          ! print*, j1, m1, j2, m2, Wig3J_Array_JJK(idj, j1, m1)
        end do
      end do
    end do
    !
    return
  end subroutine Init_Wigner3J_Array_JJK
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! evaluation of Wig3J symbol

  real(kind=8) function Wig3J(j1, j2, j3, m1, m2, m3)
    !
    ! Input: integer j1, j2, j3, m1, m2, m3, twice the physical parameters
    ! Output: real wig3j, the value of the Wigner 3J coefficient
    !
    use iso_c_binding
    implicit none
    interface
      real(c_double) function gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc) bind(c)
        ! https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_coupling_3j
        use iso_c_binding
        integer(c_int), value :: two_ja, two_jb, two_jc, two_ma, two_mb, two_mc
      end function gsl_sf_coupling_3j
    end interface
    integer(c_int), intent(in) :: j1, j2, j3, m1, m2, m3

    ! short circuit nonphysical m values
    if (.not. ((abs(m1) .le. j1) .and. (abs(m2) .le. j2) .and. (abs(m3) .le. j3))) then
       !! write(0,*) "nonphysical 3j: ", j1, j2, j3, m1, m2, m3
      wig3j = 0
      return
    end if

    wig3j = gsl_sf_coupling_3j(j1, j2, j3, m1, m2, m3)
    return
  end function Wig3J
  !
  !
end module Wigner3J

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
