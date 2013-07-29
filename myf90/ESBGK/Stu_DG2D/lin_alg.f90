module lin_alg
  implicit none


contains

  subroutine diag_matrix_invert(A,A_inv,n)
    ! of course this is dumb...but until I get something put together for
    ! non-diagonal matrices, this will have to do.
    integer, intent(in) :: n
    real(kind=4),dimension(n,n),intent(in) :: A
    real(kind=4),dimension(n,n),intent(inout) :: A_inv

    integer :: i
    real(kind=4) :: a_tmp

    A_inv = 0.

    do i = 1,n
       a_tmp = A(i,i)
       A_inv(i,i) = 1./a_tmp
    end do

  end subroutine diag_matrix_invert

  subroutine matrix_invert(A,A_inv,n)
    integer, intent(in) :: n
    real(kind=4),dimension(n,n),intent(inout) :: A
    real(kind=4),dimension(n,n),intent(inout) :: A_inv
    
    ! local variables
    integer, dimension(n) :: indx

    call MIGS(A,n,A_inv,indx)



  end subroutine matrix_invert

  ! Updated 10/24/2001.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                       !
  ! Please Note:                                                          !
  !                                                                       !
  ! (1) This computer program is written by Tao Pang in conjunction with  !
  !     his book, "An Introduction to Computational Physics," published   !
  !     by Cambridge University Press in 1997.                            !
  !                                                                       !
  ! (2) No warranties, express or implied, are made for this program.     !
  !                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE MIGS (A,N,X,INDX)
    !
    ! Subroutine to invert matrix A(N,N) with the inverse stored
    ! in X(N,N) in the output.  Copyright (c) Tao Pang 2001.
    !
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: N
    INTEGER :: I,J,K
    INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
    REAL(kind=4), INTENT (INOUT), DIMENSION (N,N):: A
    REAL(kind=4), INTENT (OUT), DIMENSION (N,N):: X
    REAL(kind=4), DIMENSION (N,N) :: B
    !
    DO I = 1, N
       DO J = 1, N
          B(I,J) = 0.0
       END DO
    END DO
    DO I = 1, N
       B(I,I) = 1.0
    END DO
    !
    CALL ELGS (A,N,INDX)
    !
    DO I = 1, N-1
       DO J = I+1, N
          DO K = 1, N
             B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
          END DO
       END DO
    END DO
    !
    DO I = 1, N
       X(N,I) = B(INDX(N),I)/A(INDX(N),N)
       DO J = N-1, 1, -1
          X(J,I) = B(INDX(J),I)
          DO K = J+1, N
             X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
          END DO
          X(J,I) =  X(J,I)/A(INDX(J),J)
       END DO
    END DO
  END SUBROUTINE MIGS
  !
  SUBROUTINE ELGS (A,N,INDX)
    !
    ! Subroutine to perform the partial-pivoting Gaussian elimination.
    ! A(N,N) is the original matrix in the input and transformed matrix
    ! plus the pivoting element ratios below the diagonal in the output.
    ! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
    !
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: N
    INTEGER :: I,J,K,ITMP
    INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
    REAL(kind=4) :: C1,PI,PI1,PJ
    REAL(kind=4), INTENT (INOUT), DIMENSION (N,N) :: A
    REAL(kind=4), DIMENSION (N) :: C
    !
    ! Initialize the index
    !
    DO I = 1, N
       INDX(I) = I
    END DO
    !
    ! Find the rescaling factors, one from each row
    !
    DO I = 1, N
       C1= 0.0
       DO J = 1, N
          !C1 = AMAX1(C1,ABS(A(I,J)))
          C1 = max(C1,ABS(A(I,J)))
       END DO
       C(I) = C1
    END DO
    !
    ! Search the pivoting (largest) element from each column
    !
    DO J = 1, N-1
       PI1 = 0.0
       DO I = J, N
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
             PI1 = PI
             K   = I
          ENDIF
       END DO
       !
       ! Interchange the rows via INDX(N) to record pivoting order
       !
       ITMP    = INDX(J)
       INDX(J) = INDX(K)
       INDX(K) = ITMP
       DO I = J+1, N
          PJ  = A(INDX(I),J)/A(INDX(J),J)
          !
          ! Record pivoting ratios below the diagonal
          !
          A(INDX(I),J) = PJ
          !
          ! Modify other elements accordingly
          !
          DO K = J+1, N
             A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
          END DO
       END DO
    END DO
    !
  END SUBROUTINE ELGS

end module lin_alg
