PROGRAM DISPDEMO
  !
  ! Example driver for DISPMODULE. To test this program issue from a command window:
  !
  !    $ make demo compiler=xxx platform=yyy
  !    $ dispdemo
  !
  ! where $ is the shell prompt, xxx is the compiler name (e.g. gfortran) and platform is
  ! either windows or unix. The default compiler name is f95 and the default platform is
  ! unix. In some cases it may be necessary to edit the makefile (for example to change
  ! the object file extension from .o to .obj). On some platforms, use gmake instead of make.
  !
  ! The program should display the output shown in comments at the bottom of this file.
  ! The examples given here are from the paper submitted to TOMS in April 2008 with two
  ! minor changes: Advance = 'double' has been set here, and the unit for the fourth example
  ! is -3 (for the screen) instead of 8.
  !
  USE DISPMODULE
  ! USE DISP_R16MOD
  IMPLICIT NONE
  INTEGER, PARAMETER :: RK = SELECTED_REAL_KIND(6), N = 3
  REAL(RK) :: A(N,N), B(N,N), X
  INTEGER I, J, K(5)
  CALL DISP_SET(ADVANCE = 'DOUBLE')
  FORALL(I=1:N, J=1:N)
    A(I,J) = EXP(REAL(I+J-1, RK))
    B(I,J) = EXP(REAL(I**J, RK))
  END FORALL
  CALL DISP('A = ', A)
  CALL DISP(B)
  CALL DISP(A(1:2,:),'F0.5')
  CALL DISP('MATRIX', A, STYLE='UNDERLINE & NUMBER', UNIT=-3, DIGMAX=4)
  K = (/-3,0,12,14,0/)
  CALL DISP('K', K, STYLE='PAD', ORIENT='ROW', SEP=' ', ZEROAS='.')
  X = 1.5
  CALL DISP('The square of '//TOSTRING(X)//' is '//TOSTRING(X*X))
  CALL DISP_SET(MATSEP = ' | ')
  CALL DISP((/11,12,13/), ADVANCE='NO')
  CALL DISP((/.TRUE., .FALSE., .TRUE./), ADVANCE='NO')
  CALL DISP((/'A','B','C'/))
END PROGRAM DISPDEMO
