subroutine diff_x1(N_max,N,A,DA,Domain_Num)
  use Legendre     ! Legendre.f90
  use Metric_Var   ! Metric_Var.f90
  implicit none
  ! Declare subroutine arguments
    integer:: N_max(1:2)
    integer:: N(1:2)
    integer:: Domain_Num
    real(kind=8)::  A(0:N_max(1),0:N_max(2))
    real(kind=8):: DA(0:N_max(1),0:N_max(2))
  !
  ! Declare local arguments
    integer:: i,j
    integer:: N1,N2
  !
    N1=N(1); N2=N(2);

    DA(0:N1,0:N2) = dxi1_dx1(0:N1,0:N2,Domain_Num) * &
         matmul(Diff_xi1(0:N1,0:N1,N1), A(0:N1,0:N2)) &
                  + dxi2_dx1(0:N1,0:N2,Domain_Num) * &
         matmul(A(0:N1,0:N2), Diff_xi2(0:N2,0:N2,N2)) 

  return

end subroutine diff_x1

!========================================================================

subroutine diff_x2(N_max,N,A,DA,Domain_Num)
  use Legendre     ! Legendre.f90
  use Metric_Var   ! Metric_Var.f90
  implicit none

  ! Declare subroutine arguments
    integer:: N_max(1:2)
    integer:: N(1:2)
    integer:: Domain_Num
    real(kind=8)::  A(0:N_max(1),0:N_max(2))
    real(kind=8):: DA(0:N_max(1),0:N_max(2))

  ! Declare local arguments
    integer:: i,j
    integer:: N1,N2

  ! subroutine begin
  !
    N1=N(1); N2=N(2);

    DA(0:N1,0:N2) = dxi1_dx2(0:N1,0:N2,Domain_Num) * &
         matmul(Diff_xi1(0:N1,0:N1,N1), A(0:N1,0:N2)) &
                  + dxi2_dx2(0:N1,0:N2,Domain_Num) * &
         matmul(A(0:N1,0:N2), Diff_xi2(0:N2,0:N2,N2)) 

  return

end subroutine diff_x2

!========================================================================
subroutine differentiate_xi1(N_max,N,A,DA)
  use Legendre     ! Legendre.f90
  implicit none
  ! Declare subroutine arguments
    integer:: N_max(1:2)
    integer:: N(1:2)
    real(kind=8)::  A(0:N_max(1),0:N_max(2))
    real(kind=8):: DA(0:N_max(1),0:N_max(2))
  !
  ! Declare local arguments
    integer:: i
    integer:: N1,N2
  !
  ! subroutine begin
  !
    N1=N(1); N2=N(2);

    DA(0:N1,0:N2)=matmul(Diff_xi1(0:N1,0:N1,N1),A(0:N1,0:N2))

  return

end subroutine differentiate_xi1
!====================================================================!

subroutine differentiate_xi2(N_max,N,A,DA)
  use Legendre     ! Legendre.f90
  implicit none
  ! Declare subroutine arguments
    integer:: N_max(1:2)
    integer:: N(1:2)
    real(kind=8)::  A(0:N_max(1),0:N_max(2))
    real(kind=8):: DA(0:N_max(1),0:N_max(2))
  !
  ! Declare local arguments
    integer:: i
    integer:: N1,N2
  !
  ! subroutine begin
  !
    N1=N(1); N2=N(2);
    DA(0:N1,0:N2)=matmul(A(0:N1,0:N2),Diff_xi2(0:N2,0:N2,N2))

  return

end subroutine differentiate_xi2

!====================================================================!
