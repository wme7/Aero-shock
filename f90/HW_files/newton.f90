program newton
  !
  ! Solves F(x)=0 by Newtons's Method
  !

  implicit none

  integer :: its = 0       ! iteration counter
  integer :: maxits = 20   ! maximun number of iterations
  integer :: converged = 0 ! convergence flag
  real :: eps = 1.0e-6     ! maximun error
  real :: x = 2.0          ! starting to guess

  ! introduce a new form of the do loop
  do while (converged == 0 .and. its < maxits)
     x = x - f(x) / df(x)
     write(*,*) x, f(x)
     its = its + 1
     if (abs(f(x)) <= eps) converged = 1
  end do
  if (converged == 1) then
     write (*,*) 'Newton converged'
  else
     write (*,*) 'Newton did not converged'
  end if

Contains
  function f(x)
    real :: f, x
    f = x**3 + x - 3.0
  end function f

  function df(x)
    !first derivative of f(x)
    real :: df, x
    df = 3*x**2 + 1
  end function df

end program newton
