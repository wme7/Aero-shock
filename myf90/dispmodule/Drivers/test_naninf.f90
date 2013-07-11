PROGRAM TEST_NANINF

  ! PROGRAM TO TEST DISPLAY OF NOT-A-NUMBER-S AND INFINITIES BY DISPMODULE
  ! (Appropriate for compilers that support IEEE-arithmetic as described in the Fortran 2003 standard)

  USE TEST_NANINF_MOD
  ! USE DISP_R16MOD ! uncomment this line if testing of disp_r16mod (quad precision) is required
  !
  ! Copyright (c) 2008, Kristján Jónasson, Dept. of Computer Science, University of
  ! Iceland (jonasson@hi.is). This software is free. For details see the file README.

  implicit none
  real(srk) :: zero = 0._srk, big = 1.e20, nan, inf, minf
  nan = zero/zero
  inf = exp(big)
  minf = -inf

  call test_ni(nan, inf, minf)

END PROGRAM TEST_NANINF

