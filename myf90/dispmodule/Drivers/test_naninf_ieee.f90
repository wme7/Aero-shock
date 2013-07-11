
PROGRAM TEST_NANINF

  ! PROGRAM TO TEST DISPLAY OF NOT-A-NUMBER-S AND INFINITIES BY DISPMODULE
  ! (Appropriate for compilers that support IEEE-arithmetic as described in the Fortran 2003 standard)

  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE TEST_NANINF_MOD
  ! USE DISP_R16MOD ! uncomment this line if testing of disp_r16mod (quad precision) is required
  !
  ! Copyright (c) 2008, Kristján Jónasson, Dept. of Computer Science, University of
  ! Iceland (jonasson@hi.is). This software is free. For details see the file README.

  implicit none
  real(srk) nan, inf, minf
  
  nan = ieee_value(0._srk, ieee_quiet_nan)
  inf = ieee_value(0._srk, ieee_positive_inf)
  minf = ieee_value(0._srk, ieee_negative_inf)

  call test_ni(nan, inf, minf)

END PROGRAM TEST_NANINF

