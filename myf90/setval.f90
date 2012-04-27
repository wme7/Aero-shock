module setval

  !make the global data available
  use dg

  implicit none

contains
  subroutine setv(a,b)
    real :: a, b

    ! set arguments to global values
    a = a_g; b = b_g
  end subroutine setv
end module setval
