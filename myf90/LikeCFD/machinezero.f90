!*******************************************************************************
!* This program finds the machine zero.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*******************************************************************************
 program machine_zero
 implicit none
 integer, parameter :: sp = kind(1.0)
 integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))

 real(dp) :: temp, one, two ! <-- Replace dp by sp to run with single precision.
 integer  :: i

!Assign constants:
!(These will be single precision if `one' and `two' are declared so.)
   one = 1.0_dp
   two = 2.0_dp

!Find the machine zero.
    i = 0
  do
    i = i + 1
    temp = one + two**(-i)
    if (temp == one) exit  ! This means that 2**(-i) is considered as zero.
  end do

  write(*,*) "Machine zero(single) = ", two**(-i)
  write(*,*) "which is equal to 2^i: ", " i =", i

 end program machine_zero
!*******************************************************************************



